###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2018'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'

import os
import gzip
import hashlib
import logging
import multiprocessing as mp
import shutil
import subprocess
import tempfile
from typing import Dict, List, NamedTuple, Optional, Sequence, Tuple, Union

import numpy as np

from tqdm import tqdm

from gtdb_migration_tk.biolib_lite.common import check_file_exists, remove_extension, make_sure_path_exists
from gtdb_migration_tk.biolib_lite.external.execute import check_on_path
from gtdb_migration_tk.biolib_lite.seq_tk import genome_size


# Below this many bases Prodigal cannot train on the genome itself, and its
# precalculated parameters are used instead.
MIN_BASES_TO_TRAIN = 100000

# The tables tried when nothing says which one a genome uses.
CANDIDATE_TABLES = (4, 11)

# Table 4 is taken only where it codes appreciably more of the genome.
DENSITY_MARGIN = 0.05
DENSITY_FLOOR = 0.7

# Files are read and hashed a block at a time.
CHUNK = 1 << 20


class ProdigalTask(NamedTuple):
    """Everything the gene calling of ONE genome needs.

    The pool pickles the callable and its argument once per task, so a worker is
    given its own genome's table rather than the mapping for the whole release.
    The output paths are given rather than composed here: what the files are
    called is the caller's convention, not this wrapper's.
    """

    genome_id: str
    genome_file: str
    translation_table: Optional[int]
    aa_gene_file: str
    nt_gene_file: str
    gff_file: str
    checksum_file: Optional[str] = None
    tmp_root: Optional[str] = None
    called_genes: bool = False
    meta: bool = False
    closed_ends: bool = False


class ConsumerData(NamedTuple):
    """What the caller learns about one genome's called genes."""

    aa_gene_file: str
    nt_gene_file: str
    gff_file: str
    best_translation_table: int
    coding_density_4: float
    coding_density_11: float
    checksum: Optional[str] = None


def run_prodigal(cmd: List[str], genome_file: str, gff_file: str, tmp_dir: str) -> None:
    """Run Prodigal over one genome, feeding it the genome on stdin.

    Prodigal reads stdin when given no -i, so a gzipped genome needs no
    uncompressed copy on disk. A non-zero exit raises, carrying what Prodigal said.

    Parameters
    ----------
    cmd : list of str
        Prodigal and its arguments, without -i.
    genome_file : str
        Genome to feed it, optionally gzipped.
    gff_file : str
        File its GFF output is written to.
    tmp_dir : str
        Directory its stderr is collected in.

    @return: None
    """

    opener = gzip.open if genome_file.endswith('.gz') else open
    stderr_file = os.path.join(tmp_dir, 'prodigal.stderr')

    # stdout and stderr both go to files, so neither can fill a pipe and deadlock
    # the process while this one is still writing the genome into its stdin
    with open(gff_file, 'wb') as gff_out, open(stderr_file, 'wb') as err_out:
        proc = subprocess.Popen(cmd, stdin=subprocess.PIPE,
                                stdout=gff_out, stderr=err_out)
        try:
            with opener(genome_file, 'rb') as genome:
                shutil.copyfileobj(genome, proc.stdin)
        finally:
            proc.stdin.close()
            proc.wait()

    if proc.returncode != 0:
        with open(stderr_file) as handle:
            complaint = ' '.join(handle.read().split())
        raise RuntimeError('prodigal failed on {} with exit code {}: {}'.format(
            genome_file, proc.returncode, complaint or '(no output on stderr)'))


def scratch_dir(tmp_root: Optional[str], results: Sequence[Optional[str]]) -> str:
    """Make a private directory for the intermediate files of one genome.

    Everything in it is deleted when the genome is done, so it must not be able to
    contain a result. The results are now written where they belong -- a genome's
    own prodigal directory -- and a scratch directory that turned out to be an
    ancestor of one would take the release's gene calls with it. mkdtemp() makes a
    fresh, uniquely named directory, so this cannot happen; it is checked rather
    than assumed, because the cost of being wrong is not recoverable.

    Parameters
    ----------
    tmp_root : str
        Directory to create the scratch directory in, or None for the default.
    results : sequence of str
        The files this genome will produce, which must lie outside it.

    @return: path of the scratch directory, which the caller must remove.
    """

    scratch = tempfile.mkdtemp(dir=tmp_root)
    root = os.path.realpath(scratch) + os.sep

    for result in results:
        if result and os.path.realpath(result).startswith(root):
            shutil.rmtree(scratch, ignore_errors=True)
            raise RuntimeError(
                'refusing to run: {} would be written inside the scratch '
                'directory {}, which is deleted when the genome is done'.format(
                    result, scratch))

    return scratch


def call_genes(task: ProdigalTask) -> Tuple[str, str, str, str, int, float, float, Optional[str]]:
    """Call the genes of one genome and report what they were called under.

    A module-level function rather than a method, so that what crosses to the
    worker is one genome's task and not whatever an instance happens to hold.

    The results are compressed straight to the paths the task names, and the
    checksum of the protein file is taken from the bytes already passing through:
    writing them to a scratch directory for the caller to move and re-read costs a
    copy of every file and a decompression of every protein file, the second of
    them in the caller and so in one process however many are calling genes.

    Parameters
    ----------
    task : ProdigalTask
        The genome, its table or None, and where the results go.

    @return: (genome_id, aa file, nt file, gff file, table, density 4, density 11,
             checksum), the densities -1 where they were not measured.
    """

    best_translation_table = -1
    table_coding_density = {4: -1, 11: -1}

    for path in (task.aa_gene_file, task.nt_gene_file, task.gff_file):
        make_sure_path_exists(os.path.dirname(path))

    if task.called_genes:
        shutil.copyfile(os.path.abspath(task.genome_file), task.aa_gene_file)
        return (task.genome_id, task.aa_gene_file, task.nt_gene_file, task.gff_file,
                best_translation_table, table_coding_density[4],
                table_coding_density[11], None)

    scratch = scratch_dir(task.tmp_root, (task.aa_gene_file, task.nt_gene_file,
                                          task.gff_file, task.checksum_file))
    try:
        total_bases = genome_size(task.genome_file)

        translation_tables = ([task.translation_table] if task.translation_table
                              else list(CANDIDATE_TABLES))

        # the density exists to choose between tables; with one table there is
        # nothing to choose, so the GFF is not parsed and no mask is built
        measure_density = len(translation_tables) > 1

        for translation_table in translation_tables:
            table_dir = os.path.join(scratch, str(translation_table))
            os.makedirs(table_dir)

            aa_gene_file_tmp = os.path.join(table_dir, 'genes.faa')
            nt_gene_file_tmp = os.path.join(table_dir, 'genes.fna')
            gff_file_tmp = os.path.join(table_dir, 'genes.gff')

            # too small to train on the genome itself, so use Prodigal's own parameters
            proc_str = 'meta' if (total_bases < MIN_BASES_TO_TRAIN or task.meta) else 'single'

            cmd = ['prodigal', '-m', '-p', proc_str, '-q', '-f', 'gff',
                   '-g', str(translation_table),
                   '-a', aa_gene_file_tmp, '-d', nt_gene_file_tmp]
            if task.closed_ends:
                cmd.append('-c')

            run_prodigal(cmd, task.genome_file, gff_file_tmp, table_dir)

            if measure_density:
                parser = ProdigalGeneFeatureParser(gff_file_tmp)
                coding_bases = sum(parser.coding_bases(seq_id) for seq_id in parser.genes)
                table_coding_density[translation_table] = float(coding_bases) / total_bases

        if measure_density:
            best_translation_table = 11
            if (table_coding_density[4] - table_coding_density[11] > DENSITY_MARGIN
                    and table_coding_density[4] > DENSITY_FLOOR):
                best_translation_table = 4
        else:
            best_translation_table = task.translation_table

        best_dir = os.path.join(scratch, str(best_translation_table))
        checksum = None
        for produced, final in (('genes.faa', task.aa_gene_file),
                                ('genes.fna', task.nt_gene_file),
                                ('genes.gff', task.gff_file)):
            digest = compress_to(os.path.join(best_dir, produced), final)
            if final == task.aa_gene_file:
                checksum = digest

        if task.checksum_file and checksum:
            with open(task.checksum_file, 'w') as handle:
                handle.write(checksum)
    finally:
        shutil.rmtree(scratch, ignore_errors=True)

    return (task.genome_id, task.aa_gene_file, task.nt_gene_file, task.gff_file,
            best_translation_table, table_coding_density[4],
            table_coding_density[11], checksum)


def compress_to(source: str, destination: str) -> str:
    """Gzip one file to its destination, hashing what goes in.

    The digest is of the UNCOMPRESSED bytes, which is what vouches for the genes
    rather than for the gzip container, and it costs nothing here: the bytes are
    being read anyway.

    SHA-1, because that is what biolib_lite.checksum.sha256_rb() computes despite
    its name, and that function is what reads these digests back to decide a
    genome can be skipped. Every .sha256 file of every past release holds a SHA-1;
    computing a real SHA-256 here would match none of them and would have the next
    run call the genes of the entire release again.

    Parameters
    ----------
    source : str
        Uncompressed file Prodigal produced.
    destination : str
        Path of the gzipped copy.

    The gzip is built beside the destination and moved onto it, rather than
    written into it. os.replace() is atomic within a filesystem, so a reader --
    the next run deciding whether this genome still needs calling, or the command
    that takes the release -- sees either the whole file or the one that was there
    before, and never a half-written one.

    That matters because these files are the release itself rather than scratch,
    and two writers can meet on one genome: prodigal divides a release into
    batches of disjoint genomes, so it takes two machines holding the SAME batch,
    which --reclaim does deliberately and an expired lease does to a machine that
    has stalled rather than stopped. Writing into the destination directly left
    them interleaving into one file and leaving a corrupt gzip behind.

    It does not make two writers safe, and is not meant to: the proteins of one
    and the .sha256 of the other still pair up wrongly. What it buys is that the
    mismatch is DETECTED -- the digest will not match the file, so the genome is
    called again by the next run -- where a corrupt gzip was simply carried into
    the release.

    @return: digest of the uncompressed bytes, in the form sha256_rb() returns.
    """

    digest = hashlib.sha1()

    # beside the destination, so the move is a rename within one filesystem and
    # not a copy; the leading dot keeps it out of anything globbing the directory
    handle, staged = tempfile.mkstemp(
        dir=os.path.dirname(destination),
        prefix='.{}.'.format(os.path.basename(destination)))
    os.close(handle)

    try:
        with open(source, 'rb') as f_in, gzip.open(staged, 'wb') as f_out:
            while True:
                block = f_in.read(CHUNK)
                if not block:
                    break
                digest.update(block)
                f_out.write(block)

        os.replace(staged, destination)
    except BaseException:
        # including KeyboardInterrupt: a half-written file left beside the
        # genome would be swept up by nothing, this directory not being scratch
        try:
            os.unlink(staged)
        except OSError:
            pass
        raise

    return digest.hexdigest()


class Prodigal(object):
    """Wrapper for running Prodigal in parallel."""

    def __init__(self, cpus: int, verbose: bool = True) -> None:
        """Initialization.

        Prodigal is checked for here rather than when the first genome is reached.

        Parameters
        ----------
        cpus : int
            Number of genomes to process at once.
        verbose : bool
            Report progress.

        @return: None
        """

        self.logger = logging.getLogger('timestamp')

        check_on_path('prodigal')

        self.cpus = cpus
        self.verbose = verbose

    def run(self, tasks: Sequence[ProdigalTask]) -> Dict[str, ConsumerData]:
        """Call genes for a set of genomes, one task per genome.

        The caller says where each genome's results go and what they are called;
        this only spreads the work and collects what came back.

        Parameters
        ----------
        tasks : sequence of ProdigalTask
            One per genome, each naming its own output paths.

        @return: genome ID to the summary statistics of its called genes.
        """

        file_type = 'scaffolds' if (tasks and tasks[0].meta) else 'genomes'

        if self.verbose:
            self.logger.info('Identifying genes within %s: ' % file_type)

        # imap_unordered, so that a genome Prodigal failed on raises here rather
        # than leaving the run to report success with a genome missing. The results
        # are collected in the parent, as the vendored Parallel class collected them
        summary_stats = {}
        with mp.Pool(processes=self.cpus) as pool:
            results = pool.imap_unordered(call_genes, tasks)
            for produced in tqdm(results, total=len(tasks),
                                 ncols=100,
                                 unit=file_type.rstrip('s'),
                                 disable=not self.verbose):
                genome_id, *stats = produced
                summary_stats[genome_id] = ConsumerData(*stats)

        return summary_stats


class ProdigalGeneFeatureParser():
    """Parses prodigal gene feature files (GFF) output."""

    def __init__(self, filename: str) -> None:
        """Initialization.

        Parameters
        ----------
        filename : str
            GFF file to parse.

        @return: None
        """
        check_file_exists(filename)

        self.genes: Dict[str, List[List[int]]] = {}
        self.last_coding_base: Dict[str, int] = {}

        self.__parseGFF(filename)

        self.coding_base_masks: Dict[str, np.ndarray] = {}
        for seq_id in self.genes:
            self.coding_base_masks[seq_id] = self.__build_coding_base_mask(seq_id)

    def __parseGFF(self, filename: str) -> None:
        """Read the genes of each contig as a list of [start, end] intervals.

        Intervals rather than gene IDs of this module's own making: the counter
        those IDs came from was reset only when a contig was first met, so a GFF
        returning to an earlier contig overwrote that contig's own genes. Nothing
        read the IDs.

        Parameters
        ----------
        filename : str
            GFF file to parse.

        @return: None
        """
        with open(filename) as handle:
            for line in handle:
                if line.startswith('#'):
                    continue

                line_split = line.split('\t')
                seq_id = line_split[0]
                if seq_id not in self.genes:
                    self.genes[seq_id] = []
                    self.last_coding_base[seq_id] = 0

                start = int(line_split[3])
                end = int(line_split[4])

                self.genes[seq_id].append([start, end])
                self.last_coding_base[seq_id] = max(self.last_coding_base[seq_id], end)

    def __build_coding_base_mask(self, seq_id: str) -> np.ndarray:
        """Mark which bases of a contig are coding.

        A mask rather than a sum of gene lengths, so that overlapping genes are
        counted once.

        Parameters
        ----------
        seq_id : str
            Unique id of sequence.

        @return: one entry per base, True where the base is within a gene.
        """

        # bool, not the float64 np.zeros() gives by default: one byte a base
        # rather than eight, and every worker holds a genome's worth.
        # last_coding_base + 1 because a GFF counts bases from 1: without the
        # extra entry the final base of the rightmost gene on every contig fell
        # off the end of the array and was never counted as coding
        coding_base_mask = np.zeros(self.last_coding_base[seq_id] + 1, dtype=bool)
        for pos in self.genes[seq_id]:
            coding_base_mask[pos[0]:pos[1] + 1] = True

        return coding_base_mask

    def coding_bases(self, seq_id: str, start: int = 0, end: Optional[int] = None) -> int:
        """Number of coding bases in a contig between [start, end).

        Parameters
        ----------
        seq_id : str
            Unique id of sequence.
        start : int
            Start calculation at this position in sequence.
        end : int
            End calculation just before this position; None for the last gene's end.

        @return: number of coding bases, 0 for a contig with no genes.
        """

        # check if sequence has any genes
        if seq_id not in self.genes:
            return 0

        # set end to just past the last coding base if not specified, the range
        # being half open and the bases counted from 1
        if end is None:
            end = self.last_coding_base[seq_id] + 1

        return int(np.count_nonzero(self.coding_base_masks[seq_id][start:end]))
