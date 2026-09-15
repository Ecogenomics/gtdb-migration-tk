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


import os
import gzip
import multiprocessing as mp
import shutil
import tarfile
import urllib.request
from multiprocessing.queues import Queue
from typing import Dict, List, TextIO, Tuple

from tqdm import tqdm

from gtdb_migration_tk import config
from gtdb_migration_tk.ncbi_genome_sync import MD5_LINE_RE


# Table of the genomes selected for a new GTDB release, written by the
# select_genomes command into its output directory. It is gzipped, as the files
# of a release are once they are part of GTDB, and because the table runs to one
# row per genome in NCBI: a few million lines of highly repetitive accessions and
# FTP paths, which compress to roughly a tenth of their size.
SELECTED_GENOMES_FILE = 'gtdb_selected_genomes.tsv.gz'

# One row of that table: accession, ftp_path, version_status,
# excluded_from_refseq, gbrs_paired_asm, notes.
SelectedRow = Tuple[str, str, str, str, str, str]

# Header of the table. It is '#'-prefixed so the file is read by the same
# readers as an NCBI assembly summary file, which skip comment lines, while
# still naming its columns for anyone opening it.
#
# The first four columns are exactly the four ncbi_genome_sync reads, in the
# order it writes its own .fail and .bad files, so this table can be handed
# straight to it as the list of genomes to mirror. It needs assembly_accession
# and ftp_path, and renders assembly_status.txt from version_status and
# excluded_from_refseq; the last two columns it simply ignores, as every reader
# of these tables locates columns by name.
#
# The notes column carries what would otherwise be a per genome warning in the
# log. It exists because NCBI's summary files routinely pair a GenBank assembly
# with a RefSeq assembly they do not themselves list: too common to report one
# line at a time, and too material to drop, since it is the reason a genome that
# appears to be covered by RefSeq was taken from GenBank instead.
SELECTED_GENOMES_HEADER = ('#assembly_accession\tftp_path\tversion_status'
                           '\texcluded_from_refseq\tgbrs_paired_asm\tnotes')


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


# --------------------------------------------------------------- NCBI metadata sync

# Where NCBI publishes the data a GTDB release is built from.
NCBI_FTP = 'https://ftp.ncbi.nlm.nih.gov'
TAXDUMP_URL = NCBI_FTP + '/pub/taxonomy/taxdump.tar.gz'

# The two NCBI databases every group is taken from.
NCBI_DATABASES = ('refseq', 'genbank')

# The two groups of organisms GTDB builds from, and the directories NCBI serves
# each one's assembly summaries from (genomes/<database>/<domain>). They are
# kept apart because they are handled differently at every step after this one:
# prokaryotes are the release, fungi are selected, assessed (busco) and
# curated by a procedure of their own. One run of ncbi_metadata_sync downloads
# one group, so a fungal run cannot quietly rewrite the prokaryotic taxonomy a
# release has already been built on, and the two can be refreshed on their own
# schedules.
GROUP_PROK = 'PROK'
GROUP_FUNGI = 'FUNGI'
NCBI_GROUPS = (GROUP_PROK, GROUP_FUNGI)

NCBI_PROK_DOMAINS = ('archaea', 'bacteria')
NCBI_FUNGI_DOMAINS = ('fungi',)

NCBI_GROUP_DOMAINS = {GROUP_PROK: NCBI_PROK_DOMAINS,
                      GROUP_FUNGI: NCBI_FUNGI_DOMAINS}

# Whether the standardised taxonomy of a group keeps NCBI's subranks. The
# prokaryotic taxonomy is the 7 ranks GTDB curates, and a subphylum or subclass
# in it would be a rank GTDB has no name for. Fungal classification leans on
# those intermediate ranks, so they are kept, and the taxonomy runs to the 13
# ranks standardize_taxonomy() writes when told to.
GROUP_KEEP_SUBRANKS = {GROUP_PROK: False,
                       GROUP_FUNGI: True}

# What a group's taxonomy files are named for: ncbi_r237_prok_*.tsv against
# ncbi_r237_fungi_*.tsv, so both groups can be downloaded into the one release
# directory without either overwriting the other.
GROUP_FILE_TAG = {GROUP_PROK: 'prok',
                  GROUP_FUNGI: 'fungi'}

# Read a request in 1 MiB blocks: the assembly summary of GenBank bacteria alone
# is well over a gigabyte, so nothing may be held in memory whole.
DOWNLOAD_BLOCK = 1024 * 1024

# A download that stalls outright must fail rather than hold the release up
# overnight; NCBI answers in well under this even when busy.
DOWNLOAD_TIMEOUT = 300


def assembly_summary_downloads(group: str) -> List[Tuple[str, str, str, str]]:
    """The assembly summary files of one group, and the names to save them under.

    NCBI calls every one of these files assembly_summary.txt, distinguishing them
    only by the directory they sit in, so downloading them into one directory
    means putting the database and domain back into the name. The names built
    here are the ones GTDB has always used, and are the names select_genomes
    reads a file's database from, so the two must agree. They carry the domain
    and not the group, so a fungal file is assembly_summary_fungi_refseq.txt.gz:
    what a file holds is the domain, and the group is only which of them are
    downloaded together.

    The database and domain are returned alongside, as the taxonomy step keys the
    files it was given by them. The names end in .gz because the files are
    compressed as they are downloaded; every reader of an assembly summary goes
    through ncbi_utils.open_summary(), which takes either form.

    Parameters
    ----------
    group : str
        Group to download, GROUP_PROK or GROUP_FUNGI.

    @return: list of (database, domain, url, file name), RefSeq before GenBank.
    """

    return [(database, domain,
             '{}/genomes/{}/{}/assembly_summary.txt'.format(NCBI_FTP, database, domain),
             'assembly_summary_{}_{}.txt.gz'.format(domain, database))
            for database in NCBI_DATABASES
            for domain in NCBI_GROUP_DOMAINS[group]]


def file_checksum(file_path: str, checksum) -> str:
    """Feed a file to a hash object a block at a time.

    Parameters
    ----------
    file_path : str
        File to checksum.
    checksum : hashlib hash
        Hash object to update.

    @return: hex digest of the file.
    """

    try:
        with open(file_path, 'rb') as file_reader:
            for block in iter(lambda: file_reader.read(DOWNLOAD_BLOCK), b''):
                checksum.update(block)
    except OSError as e:
        raise OSError('cannot read {}'.format(file_path)) from e

    return checksum.hexdigest()


def download_file(url: str,
                  output_file: str,
                  quiet: bool = False,
                  compress: bool = False) -> int:
    """Download a URL to a file, leaving nothing behind if it fails.

    The bytes go to a neighbouring .partial file which is renamed into place only
    once the transfer has finished, so an interrupted download cannot leave a
    truncated file that every later step will read as complete. The rename is
    atomic within a directory, which is why the temporary file is not in /tmp.

    With compress set the file is gzipped as it arrives rather than afterwards,
    so the uncompressed form never has to exist on disk -- which for the GenBank
    bacteria summary is a gigabyte and a half that would be written only to be
    read back and thrown away.

    Parameters
    ----------
    url : str
        URL to download.
    output_file : str
        File to write.
    quiet : bool
        Suppress the progress bar.
    compress : bool
        Gzip the file as it is written.

    @return: number of bytes received, before any compression.
    """

    partial = output_file + '.partial'
    written = 0
    open_output = gzip.open if compress else open

    try:
        with urllib.request.urlopen(url, timeout=DOWNLOAD_TIMEOUT) as response:
            # absent on a chunked response, in which case the bar shows a rate
            # and a byte count but no percentage
            total = int(response.headers.get('Content-Length') or 0)
            with open_output(partial, 'wb') as handle, tqdm(
                    total=total or None, unit='B', unit_scale=True, unit_divisor=1024,
                    desc=os.path.basename(output_file), disable=quiet) as progress:
                for block in iter(lambda: response.read(DOWNLOAD_BLOCK), b''):
                    handle.write(block)
                    written += len(block)
                    progress.update(len(block))

        if total and written != total:
            raise OSError('expected {:,} bytes but received {:,}'.format(total, written))
    except Exception:
        if os.path.exists(partial):
            os.remove(partial)
        raise

    os.replace(partial, output_file)

    return written


def extract_tarball(tarball: str, output_dir: str) -> None:
    """Extract a gzipped tarball into a directory.

    Parameters
    ----------
    tarball : str
        Gzipped tar archive to extract.
    output_dir : str
        Directory to extract into; created if it does not exist.
    """

    os.makedirs(output_dir, exist_ok=True)

    with tarfile.open(tarball, 'r:gz') as archive:
        try:
            # refuses members that would write outside output_dir; the argument
            # is only available from Python 3.11.4, and is the default from 3.14
            archive.extractall(path=output_dir, filter='data')
        except TypeError:
            archive.extractall(path=output_dir)


def write_selected_genomes(selected: List[SelectedRow], output_file: str) -> None:
    """Write the table of genomes selected for a new GTDB release.

    Rows are sorted by accession rather than left in the order the summary files
    were read, so that the tables of two releases can be compared directly to
    see what the new release gained and lost.

    Parameters
    ----------
    selected : list
        Selected genomes as (accession, ftp_path, version_status,
        excluded_from_refseq, gbrs_paired_asm, notes) rows.
    output_file : str
        Gzipped table to write.
    """

    with gzip.open(output_file, 'wt') as table:
        table.write(SELECTED_GENOMES_HEADER + '\n')
        for row in sorted(selected):
            table.write('\t'.join(row) + '\n')


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
