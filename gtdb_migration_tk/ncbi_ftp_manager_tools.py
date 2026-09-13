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
import glob
import gzip
import hashlib
import multiprocessing as mp
import shutil
import tarfile
import tempfile
import urllib.request
from multiprocessing.queues import Queue
from pathlib import Path
from typing import Dict, List, TextIO, Tuple

from tqdm import tqdm

from gtdb_migration_tk import config


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


# --------------------------------------------------------------- NCBI metadata sync

# Where NCBI publishes the data a GTDB release is built from.
NCBI_FTP = 'https://ftp.ncbi.nlm.nih.gov'
TAXDUMP_URL = NCBI_FTP + '/pub/taxonomy/taxdump.tar.gz'

# The two NCBI databases and the two domains GTDB takes from them. Fungal
# genomes are downloaded by a separate procedure and are not included here.
NCBI_DATABASES = ('refseq', 'genbank')
NCBI_DOMAINS = ('archaea', 'bacteria')

# Read a request in 1 MiB blocks: the assembly summary of GenBank bacteria alone
# is well over a gigabyte, so nothing may be held in memory whole.
DOWNLOAD_BLOCK = 1024 * 1024

# A download that stalls outright must fail rather than hold the release up
# overnight; NCBI answers in well under this even when busy.
DOWNLOAD_TIMEOUT = 300


def assembly_summary_downloads() -> List[Tuple[str, str, str, str]]:
    """The assembly summary files to download, and the names to save them under.

    NCBI calls every one of these files assembly_summary.txt, distinguishing them
    only by the directory they sit in, so downloading the four into one directory
    means putting the database and domain back into the name. The names built
    here are the ones GTDB has always used, and are the names select_genomes
    reads a file's database from, so the two must agree.

    The database and domain are returned alongside, as the taxonomy step wants
    these four files individually rather than as a list. The names end in .gz
    because the files are compressed as they are downloaded; every reader of an
    assembly summary goes through ncbi_utils.open_summary(), which takes either
    form.

    @return: list of (database, domain, url, file name), RefSeq before GenBank.
    """

    return [(database, domain,
             '{}/genomes/{}/{}/assembly_summary.txt'.format(NCBI_FTP, database, domain),
             'assembly_summary_{}_{}.txt.gz'.format(domain, database))
            for database in NCBI_DATABASES
            for domain in NCBI_DOMAINS]


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
    and genomes held by both the previous release and NCBI are compared file by
    file so that only those whose sequence data has actually changed are taken
    from NCBI again. Everything else is carried over from the previous release,
    which preserves the derived files (Prodigal, Pfam, TIGRFAM) that took the
    most effort to produce.

    Every genome handled is described in the report file, and the outcome of a
    comparison is one or more of:

        incomplete    the two directories do not hold the same FASTA files
        modified      sequence data has changed, so NCBI's copy is taken
        unmodified    sequence data is unchanged, so the GTDB copy is carried over
        new_metadata  a non-sequence file has changed and was refreshed
        new_hashes    the FTP directory gained a hashes file
        old_folder_dir  neither directory has a hashes file
        to_curate     the directories do not match any expected arrangement

    Set dry_run to describe the changes in the reports without touching any file.
    """

    def __init__(self,
                 report: TextIO,
                 genomes_to_review: TextIO,
                 genome_domain_dict: Dict[str, str],
                 dry_run: bool) -> None:
        """Set up the file suffix vocabulary and the reports to be written.

        Parameters
        ----------
        report : file
            Open file recording the fate of every genome in the release.
        genomes_to_review : file
            Open file recording genomes needing manual attention.
        genome_domain_dict : dict
            Accession to domain, used to label rows of the report.
        dry_run : bool
            Report the changes that would be made without modifying any files.
        """

        self.genomic_ext = ("_genomic.fna",)
        self.from_genomic_ext= ("_cds_from_genomic.fna","_rna_from_genomic.fna")

        self.extensions = ("_genomic.gbff", "_genomic.gff", "_wgsmaster.gbff")
        self.reports = ("_assembly_report.txt", "_assembly_stats.txt", "_hashes.txt")

        self.hmmer_exts_to_gzip = config.HMMER_EXTS_TO_GZIP
        self.prodigal_exts_to_gzip = ("_protein.faa", "_protein.fna", "_protein.gff")
        self.exts_to_gzip = self.genomic_ext + self.extensions + self.hmmer_exts_to_gzip + self.prodigal_exts_to_gzip
        self.all_but_fasta = self.extensions + self.reports

        self.report = report
        self.genomes_to_review = genomes_to_review
        self.genome_domain_dict = genome_domain_dict
        self.dry_run = dry_run

        self.ignore_extensions = ["*_assembly_structure","*_cds_from_genomic.fna.gz","*_genomic_gaps.txt.gz",
                                "*_genomic.gtf.gz","*_rna_from_genomic.fna.gz","*_translated_cds.faa.gz",
                                "*_protein.faa.gz","*_feature_count.txt.gz","*_feature_table.txt.gz",
                                "*_protein.gpff.gz"]

        self.ignore_extensions_not_archived = ["*_cds_from_genomic.fna","*_genomic_gaps.txt",
                                "*_genomic.gtf","*_rna_from_genomic.fna","*_translated_cds.faa",
                                "*_protein.faa","*_feature_count.txt","*_feature_table.txt",
                                "*_protein.gpff"]

        # We can ignore to copy folders that are no longer in use for the current GTDB release
        self.deprecated_folders=['rna_silva','pfam_27','rna_silva_132','lsu_5S','rna_silva_138',
                                 'ssu_gg','ssu_gg_2013_08','ssu_silva_199_gg_taxonomy','prokka','*ltp_132_deprecated.tar.gz']

    def rreplace(self, s: str, old: str, new: str, occurrence: int) -> str:
        """Replace the last occurrences of a substring, rather than the first.

        str.replace() works from the left, which is wrong for file names: some
        assembly directories carry a .gz in the middle of their name, so only
        the final .gz denotes compression and may be stripped.

        Parameters
        ----------
        s : str
            String to modify.
        old : str
            Substring to replace.
        new : str
            Replacement substring.
        occurrence : int
            Number of occurrences to replace, counting from the right.

        @return: string with the trailing occurrences replaced.
        """
        li = s.rsplit(old, occurrence)
        return new.join(li)

    def compare_genomes(self,
                        intersect_list: List[str],
                        old_dict: Dict[str, str],
                        new_dict: Dict[str, str],
                        ftp_directory: str,
                        new_directory: str,
                        threads: int) -> None:
        """Compare genomes held by both the previous release and the NCBI FTP site.

        Each genome is examined by a worker process, which decides whether the
        NCBI copy or the previous GTDB copy should form the new release, and the
        outcome is written to the report by a single listener process.

        Parameters
        ----------
        intersect_list : list
            Accessions held by both the previous release and the FTP site.
        old_dict : dict
            Accession to genome directory for the previous release.
        new_dict : dict
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

        for gca_record in intersect_list:
            gtdb_dir = old_dict.get(gca_record)
            ftp_dir = new_dict.get(gca_record)
            target_dir = os.path.join(
                new_directory, os.path.relpath(ftp_dir, ftp_directory))

            worker_queue.put((gtdb_dir, ftp_dir, target_dir,
                             gca_record))

        for _ in range(threads):
            worker_queue.put((None, None, None, None))

        # bound before the try, so the handler cannot fail with NameError when
        # creating the processes is itself what raised
        worker_proc, write_proc = [], None

        try:
            worker_proc = [mp.Process(target=self.__worker_thread, args=(
                worker_queue, writer_queue)) for _ in range(threads)]
            write_proc = mp.Process(target=self.__listener, args=(len(intersect_list), writer_queue))
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
                status_gca = self.compare_genome_directories(
                    gtdb_dir, ftp_dir, target_dir, gca_record)
            except Exception as e:
                status_gca = "{0}\t{1}\tto_curate;{2}\n".format(
                    self.genome_domain_dict.get(gca_record, 'UNDEFINED').upper(),
                    gca_record, type(e).__name__)

            queue_out.put(status_gca)

    def __listener(self, numgenometoprocess: int, writer_queue: Queue) -> None:
        """Write the outcome of every comparison to the report.

        The report is written by this process alone, so rows from the workers
        cannot interleave.

        Parameters
        ----------
        numgenometoprocess : int
            Number of genomes being compared, used to size the progress bar.
        writer_queue : multiprocessing.Queue
            Report rows from the workers, terminated by None.
        """

        pbar = tqdm(total=numgenometoprocess)
        for item in iter(writer_queue.get, None):
            self.report.write(item)
            pbar.update()

    def add_genomes(self,
                    added_dict: Dict[str, str],
                    ftp_dir: str,
                    new_directory: str,
                    genome_domain_dict: Dict[str, str]) -> None:
        """Copy genomes new to NCBI into the new release.

        These genomes are on the FTP site but were not in the previous release,
        so there is nothing to compare against and the NCBI directory is taken
        whole, less the files GTDB does not keep.

        :TODO: Check if the new genome is a new version of an existing genome. in
        that case we overwrite the previous one and keep the same database id.
        This will cause a conflict with the remove_genomes function.

        Parameters
        ----------
        added_dict : dict
            Accession to genome directory for the genomes to add.
        ftp_dir : str
            Base directory of the FTP mirror, replaced to form the target path.
        new_directory : str
            Base directory of the new release.
        genome_domain_dict : dict
            Accession to domain, used to label rows of the report.
        """

        for gcf_record,path_record in tqdm(added_dict.items(), desc='Adding new genomes',ncols=100):

            target_dir = os.path.join(
                new_directory, os.path.relpath(path_record, ftp_dir))
            self.report.write("{0}\t{1}\tnew\n".format(
                genome_domain_dict.get(gcf_record, 'UNDEFINED').upper(), gcf_record))
            if not self.dry_run:
                shutil.copytree(path_record, target_dir, ignore=shutil.ignore_patterns(*self.ignore_extensions))

    def remove_genomes(self, removed_dict: Dict[str, str]) -> None:
        """Record the genomes NCBI no longer offers.

        These genomes are in the previous release but have gone from the FTP
        site. Nothing is deleted here: the new release is assembled in a fresh
        directory, so a genome is dropped by not being copied into it, and this
        report is the record of which genomes that applies to.

        Parameters
        ----------
        removed_dict : dict
            Accession to genome directory for the genomes to drop.
        """

        for gca_record in removed_dict:
            self.report.write("UNDEFINED\t{0}\tremoved\n".format(gca_record))

    def compare_genome_directories(self,
                          gtdb_dir: str,
                          ftp_dir: str,
                          target_dir: str,
                          genome_record: str) -> str:
        """Decide whether a genome is taken from NCBI or carried over from GTDB.

        Both directories are copied to temporary space and decompressed, so that
        genomes compare equal whether or not they happen to be gzipped, and the
        files of each are then checksummed.

        The FASTA files decide the outcome. If they differ between the two, the
        NCBI directory becomes the new release and the genome is marked modified,
        as its derived files must be regenerated. If they match, the previous
        GTDB directory is carried over instead, keeping the derived files
        (Prodigal, Pfam, TIGRFAM) that were expensive to produce, and only those
        non-sequence files whose checksums have changed are refreshed from NCBI.
        Differing sets of FASTA files mean the genome is incomplete and is
        recorded for review.

        Under dry_run nothing is copied or compared, and the genome is reported
        as an undefined comparison.

        Parameters
        ----------
        gtdb_dir : str
            Genome directory in the previous GTDB release.
        ftp_dir : str
            Genome directory on the NCBI FTP mirror.
        target_dir : str
            Genome directory to create for the new release.
        genome_record : str
            Accession of the genome.

        @return: report row of domain, accession, and the outcomes observed.
        """

        pathftpmd5 = os.path.join(ftp_dir, "md5checksums.txt")
        pathgtdbmd5 = os.path.join(gtdb_dir, "md5checksums.txt")
        target_pathnewmd5 = os.path.join(target_dir, "md5checksums.txt")
        status = []
        if not self.dry_run:
            # 12: a context manager, so the two genome copies are removed
            # even if the comparison below raises
            with tempfile.TemporaryDirectory() as tmp_ftp_dir:
                tmp_ftp_target = os.path.join(tmp_ftp_dir,'ftp', os.path.basename(target_dir))
                tmp_existing_target = os.path.join(tmp_ftp_dir, 'previous_release', os.path.basename(target_dir))

                shutil.copytree(ftp_dir, tmp_ftp_target, symlinks=True,
                                ignore=shutil.ignore_patterns("*_assembly_structure",'prokka','*_deprecated.tar.gz'))
                shutil.copytree(gtdb_dir, tmp_existing_target, symlinks=True,
                                ignore=shutil.ignore_patterns("*_assembly_structure",'prokka','*_deprecated.tar.gz'))
                for tmp_target in [tmp_ftp_target, tmp_existing_target]:
                    for compressed_file in glob.glob(tmp_target + "/*.gz"):
                        if os.path.isdir(compressed_file) is False:
                            try:
                                with gzip.open(compressed_file, 'rb') as f_in:
                                    with open(self.rreplace(compressed_file, ".gz", "", 1), 'wb') as f_out:
                                        shutil.copyfileobj(f_in, f_out)
                            except OSError as e:
                                raise OSError('failed to decompress {}'.format(
                                    compressed_file)) from e

                            os.remove(compressed_file)

                ftpdict, ftpdict_fasta = self.checksum_genome_dir(tmp_ftp_target)
                gtdbdict, gtdbdict_fasta = self.checksum_genome_dir(tmp_existing_target)

                # if the genomic.fna.gz or the protein.faa.gz are missing, we set this
                # record as incomplete
                if len(list(set(ftpdict_fasta.keys()).symmetric_difference(set(gtdbdict_fasta.keys())))) > 0:
                    self.genomes_to_review.write(
                        "ftp_dir:{}\nftpdict_fasta.keys():{}\n"
                        "gtdb_dir:{}\ngtdbdict_fasta.keys():{}\n\n".format(
                            ftp_dir, sorted(ftpdict_fasta), gtdb_dir, sorted(gtdbdict_fasta)))
                    status.append("incomplete")
                    shutil.copytree(ftp_dir, target_dir, symlinks=True, dirs_exist_ok=True,
                                    ignore=shutil.ignore_patterns(*self.ignore_extensions))
                else:
                    ftp_folder = False
                    # check if genomic.fna.gz and protein.faa.gz are similar between
                    # previous GTDB and ftp
                    for key, value in ftpdict_fasta.items():
                        if value != gtdbdict_fasta.get(key):
                            ftp_folder = True

                    # if one of the 2 files is different than the previous version , we
                    # use the ftp record over the previous GTDB one , we then need to
                    # re run the metadata generation
                    if ftp_folder:
                        if os.path.exists(target_dir):
                            shutil.rmtree(target_dir)
                        shutil.copytree(
                            ftp_dir, target_dir, symlinks=True,
                            ignore=shutil.ignore_patterns(*self.ignore_extensions))
                        for name in glob.glob(os.path.join(target_dir, '*')):
                            if name.endswith(self.exts_to_gzip):
                                with open(name, 'rb') as f_in, gzip.open(name+'.gz','wb') as f_out:
                                    f_out.writelines(f_in)
                                os.remove(name)
                        status.append("modified")

                    else:
                        # The 2 main fasta files haven't changed so we can copy the old
                        # GTDB folder over
                        extensions_to_ignore = self.ignore_extensions + self.ignore_extensions_not_archived+self.deprecated_folders
                        if os.path.exists(target_dir):
                            shutil.rmtree(target_dir)
                        shutil.copytree(
                            gtdb_dir, target_dir, symlinks=True,
                            ignore=shutil.ignore_patterns(*extensions_to_ignore))
                        # little hack here, there is 2 _protein.faa files in each folder one from NCBI, one generated by
                        # prodigal,we want to copy the one from prodigal
                        shutil.copyfile(
                            os.path.join(gtdb_dir,'prodigal',genome_record+'_protein.faa.gz'),
                            os.path.join(target_dir, 'prodigal', genome_record+'_protein.faa.gz'))

                        for path in Path(target_dir).rglob('*'):
                            name = str(path)
                            if name.endswith(self.exts_to_gzip):
                                with open(name, 'rb') as f_in, gzip.open(name+'.gz','wb') as f_out:
                                    f_out.writelines(f_in)
                                try:
                                    os.remove(name)
                                except OSError as e:
                                    print("Failed with:", e.strerror)
                                    print("Error code:", e.errno)
                            if os.path.islink(name):
                                os.unlink(name)

                        """Process each data item in parallel."""
                        # create symlink in prodigal_folder
                        new_hit_link = os.path.join(target_dir,'prodigal', genome_record + config.TIGRFAM_SYMLINK_EXT)
                        new_tophit_link = os.path.join(target_dir,'prodigal', genome_record + config.TIGRFAM_TOPHIT_SYMLINK_EXT)
                        new_out_link = os.path.join(target_dir,'prodigal', genome_record + config.TIGRFAM_OUT_SYMLINK_EXT)

                        # Symlink needs to be relative to avoid pointing to previous version of Tigrfam when we copy folder
                        output_hit_file_relative = os.path.join('.', config.TIGRFAM_MARKER_DIR, genome_record + config.TIGRFAM_EXT)
                        tigrfam_tophit_file_relative = os.path.join('.', config.TIGRFAM_MARKER_DIR, genome_record + config.TIGRFAM_TOPHIT_EXT)
                        tigrfam_out_file_relative = os.path.join('.', config.TIGRFAM_MARKER_DIR, genome_record + config.TIGRFAM_OUT_EXT)

                        os.symlink(output_hit_file_relative, new_hit_link)
                        os.symlink(tigrfam_tophit_file_relative, new_tophit_link)
                        os.symlink(tigrfam_out_file_relative, new_out_link)

                        """Process each data item in parallel."""

                        # create symlink in prodigal_folder
                        new_hit_link = os.path.join(target_dir,'prodigal', genome_record + config.PFAM_SYMLINK_EXT)
                        new_tophit_link = os.path.join(target_dir,'prodigal', genome_record + config.PFAM_TOPHIT_SYMLINK_EXT)

                        # Symlink needs to be relative to avoid pointing to previous version of Pfam when we copy folder
                        output_hit_file_relative = os.path.join('.', config.PFAM_MARKER_DIR, genome_record + config.PFAM_EXT)
                        pfam_tophit_file_relative = os.path.join('.', config.PFAM_MARKER_DIR, genome_record + config.PFAM_TOPHIT_EXT)

                        os.symlink(output_hit_file_relative, new_hit_link)
                        os.symlink(pfam_tophit_file_relative, new_tophit_link)

                        status.append("unmodified")

                        # We check if all other file of this folder are the same.
                        checksum_changed = False

                        for key, value in ftpdict.items():
                            if value != gtdbdict.get(key):
                                checksum_changed = True
                                shutil.copy2(
                                    os.path.join(tmp_ftp_target, key), os.path.join(target_dir, key))
                                if key.endswith(self.exts_to_gzip):
                                    with open(os.path.join(target_dir, key), 'rb') as f_in, gzip.open(os.path.join(target_dir, key+'.gz'),'wb') as f_out:
                                        f_out.writelines(f_in)
                                    os.remove(os.path.join(target_dir, key))
                                status.append("new_metadata")

                        # we copy the new checksum
                        if checksum_changed:
                            try:
                                shutil.copy2(pathftpmd5, target_pathnewmd5)
                            except IOError:
                                os.chmod(target_pathnewmd5, 0o664)
                                shutil.copy2(pathftpmd5, target_pathnewmd5)
                        # Only reached for an unmodified genome, and deliberately so.
                        # The target directory here was carried over from the previous
                        # release, so its copies of these NCBI files may be years old and
                        # must be checked against the mirror. The modified and incomplete
                        # branches copy the mirror wholesale, so their copies are current
                        # by construction and there is nothing to compare. This is also
                        # why _hashes.txt (NCBI's annotation_hashes.txt) is asked about
                        # only here: old_folder_dir means the previous release predates
                        # the file, and new_hashes that NCBI has since added one -- both
                        # questions about the previous release, not about the mirror.
                        for report in self.reports:
                            target_files = glob.glob(
                                os.path.join(target_dir, "*" + report))
                            ftp_files = glob.glob(os.path.join(ftp_dir, "*" + report))
                            if len(target_files) == 1 and len(ftp_files) == 1:
                                status = self.compare_md5(
                                    ftp_files[0], target_files[0], status)
                            elif len(target_files) == 0 and len(ftp_files) == 0 and report == '_hashes.txt':
                                status.append("old_folder_dir")
                            elif len(target_files) == 0 and len(ftp_files) == 1 and report == '_hashes.txt':
                                shutil.copy2(ftp_files[0], target_dir)
                                status.append("new_hashes")
                            else:
                                print("########")
                                print(target_dir)
                                print(target_files)
                                print(ftp_dir)
                                print(ftp_files)
                                print(f"IT SHOULDN'T HAPPEN ({target_dir},{ftp_dir}) ")
                                print("########")
                                status.append("to_curate")

                # 3: sorted so two runs over the same data produce identical reports
                status_record = "{0}\t{1}\t{2}\n".format(
                    self.genome_domain_dict.get(genome_record, 'UNDEFINED').upper(),
                    genome_record, ';'.join(sorted(set(status))))
                return status_record
        
        return "{0}\t{1}\t{2}\n".format(
            self.genome_domain_dict.get(genome_record, 'UNDEFINED').upper(),
            genome_record, 'undefined comparison')

    def compare_md5(self,
                    ftp_file: str,
                    target_file: str,
                    status: List[str]) -> List[str]:
        """Refresh a report file from NCBI if it has changed.

        Used for the small per-assembly reports NCBI ships alongside a genome.
        The file is replaced in place when the two copies differ, and the genome
        is then marked as carrying new metadata.

        Parameters
        ----------
        ftp_file : str
            File on the NCBI FTP mirror.
        target_file : str
            Corresponding file in the new release, overwritten if it differs.
        status : list
            Outcomes observed for this genome so far.

        @return: the status list, with new_metadata appended if the file changed.
        """

        if self.md5_calculator(ftp_file) != self.md5_calculator(target_file):
            try:
                shutil.copy2(ftp_file, target_file)
            except IOError:
                os.chmod(target_file, 0o664)
                shutil.copy2(ftp_file, target_file)
            status.append("new_metadata")
        return status

    def checksum_genome_dir(self, pathtodir: str) -> Tuple[Dict[str, str], Dict[str, str]]:
        """Checksum the files of a genome directory, FASTA apart from the rest.

        The FASTA files are kept separate because they alone decide whether a
        genome has really changed; everything else is metadata that can be
        refreshed on its own. Derived FASTA files (_cds_from_genomic,
        _rna_from_genomic) are not genome assemblies and are excluded.

        Checksums are computed here rather than read from md5checksums.txt, as
        the files have been decompressed and no longer match the sums NCBI
        published for them.

        Parameters
        ----------
        pathtodir : str
            Genome directory to checksum.

        @return: (metadata file to checksum, FASTA file to checksum) dicts.
        """

        out_dict, out_dict_fasta= {}, {}

        for name in glob.glob(os.path.join(pathtodir, '*')):
            if name.endswith(self.genomic_ext) and not name.endswith(self.from_genomic_ext):
                out_dict_fasta[os.path.basename(
                    name)] = self.sha256_calculator(name)
                os.chmod(name, 0o664)
            elif name.endswith(self.all_but_fasta):
                out_dict[os.path.basename(name)] = self.sha256_calculator(name)
                os.chmod(name, 0o664)
        return (out_dict, out_dict_fasta)

    def md5_calculator(self, file_path: str) -> str:
        """Compute the MD5 checksum of a file.

        Used for the small report files NCBI ships alongside a genome. The file
        is read a block at a time so that its size never dictates memory use.

        Parameters
        ----------
        file_path : str
            File to checksum.

        @return: hex digest of the file.
        """

        return self._checksum(file_path, hashlib.md5())

    def sha256_calculator(self, file_path: str) -> str:
        """Compute the SHA-256 checksum of a file.

        The file is read a line at a time so that genome assemblies too large to
        hold in memory can be checksummed.

        Parameters
        ----------
        file_path : str
            File to checksum.

        @return: hex digest of the file.
        """

        return self._checksum(file_path, hashlib.sha256())

    def _checksum(self, file_path: str, checksum) -> str:
        """Feed a file to a hash object a block at a time.

        Parameters
        ----------
        file_path : str
            File to checksum.
        checksum : hashlib hash
            Hash object to update.

        @return: hex digest of the file.
        """

        return file_checksum(file_path, checksum)
