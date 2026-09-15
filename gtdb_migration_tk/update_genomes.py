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

"""
update_genomes.py -- update the GTDB genome directories to match the NCBI mirror.

Each GTDB release is built by comparing the genomes held by the previous release
with the genomes NCBI currently offers on its FTP site. Genomes NCBI no longer
offers are removed, genomes new to NCBI are copied across, and genomes held by
both are checked for a changed genomic FASTA since the previous release. The
bookkeeping of which genomes fall into which of those three groups is
GenomeManager's; the copying, comparing and reporting that follows is FTPTools'.

GenomeManager decides nothing beyond what the genome directory files of the
mirror and of the previous release already say. The mirror is built from the
selection, and the selection is where a genome is accepted or passed over, so
every genome the mirror holds is wanted. RefSeq and GenBank are handled in
separate runs of the same code, told apart by the accession prefix (GCF or GCA),
because the update of each is reported separately; nothing else differs.
"""

import os
import shutil
import logging
import multiprocessing as mp
from contextlib import ExitStack
from multiprocessing.queues import Queue
from typing import Dict, List, TextIO

from tqdm import tqdm

from gtdb_migration_tk import config
from gtdb_migration_tk.ncbi_utils import MD5_LINE_RE
from gtdb_migration_tk.utils.common import count_lines


# Outcomes a genome held by both the previous release and NCBI can have, as
# written to the report. Constants rather than literals so the tests, and any
# reader of the report, name the same strings this module writes. A dry run
# reports the same outcomes as a real one; it only refrains from copying.
STATUS_FASTA_UNCHANGED = 'genomic FASTA file unchanged'
STATUS_FASTA_CHANGED = 'genomic FASTA file changed'

# NCBI's manifest of the files it serves for a genome, one "<md5>  ./<name>"
# line per file. Both the mirror and the previous release carry a copy.
MD5_MANIFEST = 'md5checksums.txt'

# Suffix NCBI appends to the assembly name for the genome assembly itself.
GENOMIC_FASTA_EXT = '_genomic.fna.gz'


class GenomeManager:
    """Update the GTDB copy of one NCBI database (RefSeq or GenBank) from the mirror.

    Every genome held by the FTP mirror is of interest, the mirror being a copy
    of the genomes selected for the release, so genome selection is simply a
    matter of reading its genome directory file. Genomes are then added,
    removed, or compared relative to the previous GTDB release. Genomes are
    tracked as dictionaries mapping an accession to its genome directory.

    One instance handles one database, named by the accession prefix it is
    given: only genomes with that prefix are read from either genome directory
    file, and the reports carry the prefix in their names so the RefSeq and
    GenBank runs of a release sit side by side in one output directory.
    """

    def __init__(self,
                 accession_prefix: str,
                 new_genome_dir: str,
                 dry_run: bool = False,
                 cpus: int = 1) -> None:
        """Record which database is handled and where the release is written.

        Parameters
        ----------
        accession_prefix : str
            Accession prefix of the database of interest, REFSEQ_PREFIX or
            GENBANK_PREFIX.
        new_genome_dir : str
            Output directory for the new release, where reports are written.
        dry_run : bool
            Report the changes that would be made without modifying any files.
        cpus : int
            Number of processes used when comparing genomes.
        """

        self.accession_prefix = accession_prefix
        self.new_genome_dir = new_genome_dir
        self.dry_run = dry_run
        self.cpus = cpus
        self.logger = logging.getLogger('timestamp')

    def report_file(self) -> str:
        """Report recording the fate of every genome of this database.

        @return: path of the report, named for the accession prefix.
        """

        return os.path.join(self.new_genome_dir,
                            'report_{}.log'.format(self.accession_prefix.lower()))

    def review_file(self) -> str:
        """Report recording genomes of this database needing manual attention.

        @return: path of the report, named for the accession prefix.
        """

        return os.path.join(self.new_genome_dir,
                            '{}_to_review.log'.format(self.accession_prefix.lower()))

    def load_genome_dirs(self, genome_dirs_file: str) -> Dict[str, str]:
        """Read the genomes of this database from a genome directory file.

        The file describes a whole release, RefSeq and GenBank together; only
        the genomes with this manager's accession prefix are kept.

        Parameters
        ----------
        genome_dirs_file : str
            Genome directory file (accession, path) of the mirror or of a release.

        @return: dict of accession to genome directory.
        """

        self.logger.info('Reading {} genomes from {}:'.format(
            self.accession_prefix, genome_dirs_file))

        genome_paths = {}
        with open(genome_dirs_file, 'r') as f:
            for line in tqdm(f, total=count_lines(genome_dirs_file)):
                gid, path, *_ = line.split('\t')
                if gid.startswith(self.accession_prefix):
                    genome_paths[gid] = path.strip()

        self.logger.info(' - identified {:,} genomes'.format(len(genome_paths)))

        return genome_paths

    def generate_genomes_to_remove(self,
                                   new_genomes: Dict[str, str],
                                   old_genomes: Dict[str, str]) -> Dict[str, str]:
        """Identify genomes present in the previous release, but no longer on the NCBI FTP site.

        Parameters
        ----------
        new_genomes : dict
            Accession to genome directory for genomes currently on the FTP site.
        old_genomes : dict
            Accession to genome directory for genomes in the previous release.

        @return: dict of accession to genome directory for genomes to remove.
        """

        removed_genomes = {gid: path for gid, path in old_genomes.items()
                           if gid not in new_genomes}
        self.logger.info('Identified {:,} genomes to remove.'.format(len(removed_genomes)))

        return removed_genomes

    def generate_genomes_to_add(self,
                                new_genomes: Dict[str, str],
                                old_genomes: Dict[str, str]) -> Dict[str, str]:
        """Identify genomes new to the NCBI FTP site since the previous release.

        Parameters
        ----------
        new_genomes : dict
            Accession to genome directory for genomes currently on the FTP site.
        old_genomes : dict
            Accession to genome directory for genomes in the previous release.

        @return: dict of accession to genome directory for genomes to add.
        """

        added_genomes = {gid: path for gid, path in new_genomes.items()
                         if gid not in old_genomes}
        self.logger.info('Identified {:,} genomes to add.'.format(len(added_genomes)))

        return added_genomes

    def generate_genomes_to_compare(self,
                                    new_genomes: Dict[str, str],
                                    old_genomes: Dict[str, str]) -> List[str]:
        """Identify genomes common to the previous release and the NCBI FTP site.

        These genomes are candidates for an update as their files may have
        changed since the previous release.

        Parameters
        ----------
        new_genomes : dict
            Accession to genome directory for genomes currently on the FTP site.
        old_genomes : dict
            Accession to genome directory for genomes in the previous release.

        @return: list of accessions to compare between the two releases.
        """

        shared_genomes = list(old_genomes.keys() & new_genomes.keys())
        self.logger.info('Identified {:,} genomes to compare.'.format(len(shared_genomes)))

        return shared_genomes

    def run_comparison(self,
                       ftp_dir: str,
                       ftp_genome_dirs: str,
                       old_genome_dirs: str) -> None:
        """Update the GTDB genome directories of this database to match the mirror.

        The genomes of the mirror are compared to those of the previous GTDB
        release. Genomes no longer at NCBI are recorded as removed, new genomes
        are copied across, and genomes common to both are checked for a changed
        genomic FASTA and carried over with their derived data when it is not.

        Parameters
        ----------
        ftp_dir : str
            Root of the NCBI FTP mirror, replaced by the new release directory
            to place each genome.
        ftp_genome_dirs : str
            Genome directory file (accession, path) for the mirror.
        old_genome_dirs : str
            Genome directory file (accession, path) for the previous release.
        """

        self.logger.info('Updating {} genomes.'.format(self.accession_prefix))

        # reports are opened for the duration of the update so they are closed,
        # and their contents kept, if the update fails part way through
        with ExitStack() as reports:
            report = reports.enter_context(open(self.report_file(), 'w', 1))
            genomes_to_review = reports.enter_context(open(self.review_file(), 'w', 1))

            old_genomes = self.load_genome_dirs(old_genome_dirs)
            new_genomes = self.load_genome_dirs(ftp_genome_dirs)

            ftptools = FTPTools(report, genomes_to_review, self.dry_run)

            removed_genomes = self.generate_genomes_to_remove(new_genomes, old_genomes)
            ftptools.remove_genomes(removed_genomes)

            added_genomes = self.generate_genomes_to_add(new_genomes, old_genomes)
            ftptools.add_genomes(added_genomes, ftp_dir, self.new_genome_dir)

            shared_genomes = self.generate_genomes_to_compare(new_genomes, old_genomes)
            ftptools.compare_genomes(shared_genomes, old_genomes, new_genomes,
                                     ftp_dir, self.new_genome_dir, self.cpus)


class FTPTools():
    """Carry out the genome directory changes required by a GTDB release.

    The decision about which genomes belong in a release is made by the managers
    in ncbi_ftp_manager.py; this class performs the resulting work. Genomes new to
    NCBI are copied into the new release, genomes that have gone are recorded,
    and for genomes held by both the previous release and NCBI the mirror's copy
    is taken and, if the genomic FASTA is unchanged, the derived data of the
    previous release (config.GTDB_DERIVED_DIRS_TO_COPY) is carried across with
    it, as that data took the most effort to produce and is still valid.

    Every genome handled is described in the report file. A genome held by both
    sides has one of the outcomes named by the STATUS_* constants above, or
    'to_curate;<exception>' if it could not be compared at all.

    Set dry_run to describe the changes in the reports without touching any file.
    """

    def __init__(self,
                 report: TextIO,
                 genomes_to_review: TextIO,
                 dry_run: bool) -> None:
        """Set up the file suffix vocabulary and the reports to be written.

        Parameters
        ----------
        report : file
            Open file recording the fate of every genome in the release.
        genomes_to_review : file
            Open file recording genomes needing manual attention.
        dry_run : bool
            Report the changes that would be made without modifying any files.
        """

        self.report = report
        self.genomes_to_review = genomes_to_review
        self.dry_run = dry_run

        # this is being maintained for backwards compatibility with previous GTDB releases, but it is not 
        # required starting with GTDB r237 as we now only sync the exact files required by GTDB
        self.ignore_extensions_compressed = ["*_assembly_structure","*_cds_from_genomic.fna.gz","*_genomic_gaps.txt.gz",
                                "*_genomic.gtf.gz","*_rna_from_genomic.fna.gz","*_translated_cds.faa.gz",
                                "*_protein.faa.gz","*_feature_count.txt.gz","*_feature_table.txt.gz",
                                "*_protein.gpff.gz"]

        self.ignore_extensions_uncompressed = ["*_cds_from_genomic.fna","*_genomic_gaps.txt",
                                "*_genomic.gtf","*_rna_from_genomic.fna","*_translated_cds.faa",
                                "*_protein.faa","*_feature_count.txt","*_feature_table.txt",
                                "*_protein.gpff"]

        self.ignore_extensions = self.ignore_extensions_compressed + self.ignore_extensions_uncompressed

    def add_genomes(self,
                    added_genomes: Dict[str, str],
                    ftp_dir: str,
                    new_directory: str) -> None:
        """Copy genomes new to NCBI into the new release.

        These genomes are on the FTP site but were not in the previous release,
        so there is nothing to compare against and the NCBI directory is taken
        whole, less the files GTDB does not keep.

        Parameters
        ----------
        added_genomes : dict
            Accession to genome directory for the genomes to add.
        ftp_dir : str
            Base directory of the FTP mirror, replaced to form the target path.
        new_directory : str
            Base directory of the new release.
        """

        for gid, path_record in tqdm(added_genomes.items(), desc='Adding new genomes', ncols=100):
            target_dir = os.path.join(new_directory, os.path.relpath(path_record, ftp_dir))
            self.report.write("{0}\tnew\n".format(gid))
            if not self.dry_run:
                shutil.copytree(path_record, target_dir, symlinks=True,
                                ignore=shutil.ignore_patterns(*self.ignore_extensions))

    def remove_genomes(self, removed_genomes: Dict[str, str]) -> None:
        """Record the genomes NCBI no longer offers.

        These genomes are in the previous release but have gone from the FTP
        site. Nothing is deleted here: the new release is assembled in a fresh
        directory, so a genome is dropped by not being copied into it, and this
        report is the record of which genomes that applies to.

        Parameters
        ----------
        removed_genomes : dict
            Accession to genome directory for the genomes to drop.
        """

        for gid in removed_genomes:
            self.report.write("{0}\tremoved\n".format(gid))

    def compare_genomes(self,
                        shared_genomes: List[str],
                        old_genome_dirs: Dict[str, str],
                        new_genome_dirs: Dict[str, str],
                        ftp_directory: str,
                        new_directory: str,
                        threads: int) -> None:
        """Compare genomes held by both the previous release and the NCBI FTP site.

        Each genome is examined by a worker process which decides whether the
        derived data (e.g. Prodigal results) from the previous GTDB release should
        be retained. All NCBI data files are always copied from the NCBI FTP mirror 
        as these are the latest versions of these files.

        Parameters
        ----------
        shared_genomes : list
            Accessions held by both the previous release and the FTP site.
        old_genome_dirs : dict
            Accession to genome directory for the previous release.
        new_genome_dirs : dict
            Accession to genome directory for the FTP mirror.
        ftp_directory : str
            Base directory of the FTP mirror, replaced to form the target path.
        new_directory : str
            Base directory of the new release.
        threads : int
            Number of worker processes.
        """

        # populate worker queue with data to process
        worker_queue = mp.Queue()
        writer_queue = mp.Queue()

        for gca_record in shared_genomes:
            gtdb_dir = old_genome_dirs.get(gca_record)
            ftp_dir = new_genome_dirs.get(gca_record)
            target_dir = os.path.join(new_directory, os.path.relpath(ftp_dir, ftp_directory))

            worker_queue.put((gtdb_dir, ftp_dir, target_dir, gca_record))

        for _ in range(threads):
            worker_queue.put((None, None, None, None))

        # bound before the try, so the handler cannot fail with NameError when
        # creating the processes is itself what raised
        worker_proc, write_proc = [], None

        try:
            worker_proc = [mp.Process(target=self.__worker_thread, args=(
                worker_queue, writer_queue)) for _ in range(threads)]
            write_proc = mp.Process(target=self.__listener, args=(len(shared_genomes), writer_queue))
            write_proc.start()

            for p in worker_proc:
                p.start()

            for p in worker_proc:
                p.join()

            writer_queue.put(None)
            write_proc.join()
        except Exception:
            for p in worker_proc:
                p.terminate()

            if write_proc is not None:
                write_proc.terminate()

            # genomes are left uncompared, which must not read as a completed run
            raise

        # a worker that died took its share of the genomes with it, and those
        # genomes are simply absent from the report rather than flagged
        failed = [p.exitcode for p in worker_proc if p.exitcode != 0]
        if failed:
            raise RuntimeError(
                '{} of {} comparison processes failed (exit codes: {}); '
                'the report is incomplete'.format(
                    len(failed), len(worker_proc), ', '.join(str(c) for c in failed)))

    def __worker_thread(self, queue_in: Queue, queue_out: Queue) -> None:
        """Compare one genome at a time until the queue is exhausted.

        Parameters
        ----------
        queue_in : multiprocessing.Queue
            Genomes to compare, as (gtdb_dir, ftp_dir, target_dir, accession);
            an accession of None ends the worker.
        queue_out : multiprocessing.Queue
            Report rows produced by the comparisons.
        """

        while True:
            gtdb_dir, ftp_dir, target_dir, gca_record = queue_in.get(
                block=True, timeout=None)
            if gca_record is None:
                break

            # a genome that cannot be compared is reported for curation; letting
            # it kill the worker silently drops every genome queued behind it
            try:
                status_gca = self.compare_genome_directories(gtdb_dir, ftp_dir, target_dir, gca_record)
            except Exception as e:
                status_gca = "{0}\tto_curate;{1}\n".format(gca_record, type(e).__name__)

            queue_out.put(status_gca)

    def __listener(self, num_genomes: int, writer_queue: Queue) -> None:
        """Write the outcome of every comparison to the report.

        The report is written by this process alone, so rows from the workers
        cannot interleave.

        Parameters
        ----------
        num_genomes : int
            Number of genomes being compared, used to size the progress bar.
        writer_queue : multiprocessing.Queue
            Report rows from the workers, terminated by None.
        """

        pbar = tqdm(total=num_genomes, desc='Comparing shared genomes', ncols=100)
        for item in iter(writer_queue.get, None):
            self.report.write(item)
            pbar.update()

    def compare_genome_directories(self,
                                   prev_gtdb_dir: str,
                                   ftp_dir: str,
                                   target_dir: str,
                                   genome_record: str) -> str:
        """Build the new release's copy of a genome held by both GTDB and NCBI.

        The mirror's directory is copied whole, md5checksums.txt included: what
        NCBI serves for a genome is what the new release carries, whether or not
        anything changed. The question is only whether the derived data of the
        previous release can come with it, and the genomic FASTA decides that.
        Its MD5 is read from the md5checksums.txt of each directory rather than
        computed, so the file is never opened and nothing has to be decompressed;
        NCBI published both sums, and the sync verified the mirror against its
        copy. If the two agree the sequence is unchanged, so gene calls, rRNA
        classifications and tRNA scans made on it still hold and the directories
        named by config.GTDB_DERIVED_DIRS_TO_COPY are copied across from the
        previous release. If they differ the derived data is left behind, to be
        regenerated from the new sequence.

        A derived directory the previous release lacks is not an error, as a
        genome may not have had every step run on it; it is noted in the review
        report so that the step can be run this time.

        Under dry_run the comparison is made in full and reported, so the report
        of a dry run is the report the real run would write; only the copying is
        withheld, from the mirror and from the previous release alike.

        Parameters
        ----------
        prev_gtdb_dir : str
            Genome directory in the previous GTDB release.
        ftp_dir : str
            Genome directory on the NCBI FTP mirror.
        target_dir : str
            Genome directory to create for the new release.
        genome_record : str
            Accession of the genome.

        @return: report row of accession and outcome, one of the STATUS_* constants.
        """

        ftp_md5 = self.genomic_fasta_md5(ftp_dir)
        prev_md5 = self.genomic_fasta_md5(prev_gtdb_dir)

        if not self.dry_run:
            # a rerun must not inherit derived data from a run made before the
            # FASTA changed, so an existing target is replaced rather than added to
            if os.path.exists(target_dir):
                shutil.rmtree(target_dir)
            shutil.copytree(ftp_dir, target_dir, symlinks=True)

        if ftp_md5 != prev_md5:
            return '{}\t{}\n'.format(genome_record, STATUS_FASTA_CHANGED)

        for derived in config.GTDB_DERIVED_DIRS_TO_COPY:
            source = os.path.join(prev_gtdb_dir, derived)
            if not os.path.isdir(source):
                # noted on a dry run too: it is part of what the run would report
                self.genomes_to_review.write(
                    '{}\tno {} directory in the previous release: {}\n'.format(
                        genome_record, derived, prev_gtdb_dir))
                continue
            if not self.dry_run:
                # symlinks kept as symlinks: prodigal/ holds version-free links to the
                # marker results beside them, and following them would copy each
                # Pfam and TIGRFAM file twice
                shutil.copytree(source, os.path.join(target_dir, derived), symlinks=True)

        return '{}\t{}\n'.format(genome_record, STATUS_FASTA_UNCHANGED)

    def genomic_fasta_md5(self, genome_dir: str) -> str:
        """Read the MD5 NCBI publishes for the genomic FASTA of a genome.

        NCBI names the file for the assembly, <accession>_<asm_name>_genomic.fna.gz,
        and names the genome directory <accession>_<asm_name>, so the entry wanted
        is known exactly and is looked up by name; this is how the rest of the
        toolkit finds the file too. The manifest is read with the same line
        pattern ncbi_genome_sync uses to mirror and verify it.

        Parameters
        ----------
        genome_dir : str
            Genome directory holding an md5checksums.txt.

        @return: hex MD5 of the genomic FASTA, as recorded in the manifest.
        """

        assembly = os.path.basename(os.path.normpath(genome_dir))
        wanted = assembly + GENOMIC_FASTA_EXT
        manifest = os.path.join(genome_dir, MD5_MANIFEST)

        with open(manifest) as handle:
            for line in handle:
                match = MD5_LINE_RE.match(line.strip())
                if not match:
                    continue
                name = match.group(2).strip()
                if name.startswith('./'):
                    name = name[2:]
                if name == wanted:
                    return match.group(1)

        raise ValueError('{} has no entry for {}'.format(manifest, wanted))
