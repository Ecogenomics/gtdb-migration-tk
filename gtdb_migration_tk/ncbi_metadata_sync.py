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
ncbi_metadata_sync.py -- download the NCBI metadata a GTDB release is built from.

NCBIMetadataSync comes before anything else in a release: it downloads the NCBI
taxonomy and the assembly summary files a release is selected from, and decides
nothing. The genomes are then chosen by select_genomes.py from the summary files
alone, and mirrored by ncbi_genome_sync.py from that selection.
"""

import os
import sys
import gzip
import datetime
import hashlib
import logging
import tarfile
import urllib.request
from typing import Dict, List, Tuple

from tqdm import tqdm

from gtdb_migration_tk.ncbi_tax_manager import TaxonomyNCBI


# Subdirectories of the root NCBI directory of a release. Everything derived from
# the NCBI taxonomy lives under TAXONOMY_DIR, the standardised 7 rank form
# included, so the taxonomy of a release is one directory to copy or archive.
TAXONOMY_DIR = 'taxonomy'
STANDARDISED_TAXONOMY_DIR = 'standardised_taxonomy'

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


class NCBIMetadataSync:
    """Download the NCBI metadata one group of a GTDB release is built from.

    A release starts from two things NCBI publishes and GTDB only reads: the
    taxonomy, and the assembly summary files describing every assembly NCBI
    holds. This fetches both into one directory, which is then the input to
    select_genomes and, later, to the commands that attach NCBI taxonomy to the
    genomes chosen.

    One run covers one group -- PROK (archaea and bacteria) or FUNGI -- because
    everything downstream of here treats the two differently. Both groups are
    downloaded the same way and both get a standardised taxonomy; what the group
    decides is which of NCBI's directories are fetched, and whether that taxonomy
    keeps NCBI's subranks. The prokaryotic taxonomy is the 7 ranks GTDB curates;
    fungal classification leans on the intermediate ranks, so the fungal taxonomy
    keeps them.

    This is the procedure the GTDB wiki gives as "Download latest NCBI taxonomy",
    the first step of "Download the latest RefSeq and GenBank assembly data", and
    "Generate 7-rank NCBI taxonomy", run in one go. The output directory is the
    root NCBI directory of a release, and is laid out as the wiki leaves it:

        <output_dir>/assembly_summary_<domain>_<database>.txt.gz
        <output_dir>/taxonomy/taxdump_<date>/
        <output_dir>/taxonomy/standardised_taxonomy/ncbi_r<release>_<group>_*.tsv

    Both groups can therefore be run into the one release directory: their
    summary files are named for their domains, and their taxonomy files for their
    group, so neither run overwrites the other's output. A group run after
    another on the same day reuses the taxdump already extracted rather than
    downloading and verifying the same 80 MB archive twice; the two groups are
    then placed against the same taxonomy, which is what a release wants.

    Three things differ from doing it by hand. NCBI names every assembly summary
    file assembly_summary.txt, so the database and domain are put back into the
    name as they are saved, and they are gzipped on the way in rather than kept
    as the 1.8 GB of text NCBI serves. The taxonomy dump is checked against the
    MD5 NCBI publishes beside it, since a truncated taxdump is not obviously
    broken until a release has been built on it, and the archive is then
    discarded: it has been verified and unpacked, and nothing reads it again.

    Nothing here decides anything: unlike the other managers in this module it
    only fetches, and hands what it fetched to the NCBI taxonomy parser.
    """

    def __init__(self, output_dir: str, group: str) -> None:
        """Record the group to download, and where it is to be written.

        Parameters
        ----------
        output_dir : str
            Output directory for the downloaded NCBI metadata.
        group : str
            Group to download, PROK or FUNGI.

        @raise ValueError: if the group is not one this toolkit knows.
        """

        if group not in NCBI_GROUPS:
            raise ValueError('unknown group {}; expected one of {}'.format(
                group, ', '.join(NCBI_GROUPS)))

        self.output_dir = output_dir
        self.group = group
        self.logger = logging.getLogger('timestamp')

    def file_tag(self) -> str:
        """What this group's taxonomy files are named for.

        @return: the group's file name tag, e.g. 'prok'.
        """

        return GROUP_FILE_TAG[self.group]

    def taxonomy_dir(self) -> str:
        """Directory holding everything derived from the NCBI taxonomy.

        @return: the taxonomy subdirectory of the output directory.
        """

        return os.path.join(self.output_dir, TAXONOMY_DIR)

    def _download(self, url: str, output_file: str, compress: bool = False) -> int:
        """Download one file into the output directory, failing the run if it cannot.

        Parameters
        ----------
        url : str
            URL to download.
        output_file : str
            File to write.
        compress : bool
            Gzip the file as it is written.

        @return: number of bytes received, before any compression.
        """

        if os.path.exists(output_file):
            self.logger.warning('Replacing existing {}.'.format(os.path.basename(output_file)))

        self.logger.info('Downloading {}'.format(url))
        try:
            written = download_file(url, output_file, compress=compress)
        except Exception as exc:
            # a half-downloaded release is worse than none: stop at the first
            # failure rather than leaving the operator to notice a missing file
            self.logger.error('Failed to download {}: {}'.format(url, exc))
            sys.exit()

        if compress:
            self.logger.info('Wrote {} ({:,} bytes, {:,} compressed).'.format(
                os.path.basename(output_file), written, os.path.getsize(output_file)))
        else:
            self.logger.info('Wrote {} ({:,} bytes).'.format(
                os.path.basename(output_file), written))

        return written

    def download_taxonomy(self, taxonomy_dir: str, date_stamp: str) -> str:
        """Download and extract the NCBI taxonomy, unless today's is already here.

        The archive is verified against NCBI's published MD5, unpacked, and then
        removed: only the extracted directory is ever read again.

        A taxdump already extracted under today's date is reused. Both groups
        need the taxonomy, and running one after the other would otherwise fetch
        and verify the same 80 MB archive twice and, worse, place the two groups
        against two different downloads of the taxonomy. A run on a later date
        names a different directory and so still gets a fresh dump.

        Parameters
        ----------
        taxonomy_dir : str
            Directory to download the taxonomy into.
        date_stamp : str
            Date the download was made, as YYYYMMDD.

        @return: directory the taxonomy was extracted into.
        """

        os.makedirs(taxonomy_dir, exist_ok=True)

        extracted = os.path.join(taxonomy_dir, 'taxdump_{}'.format(date_stamp))
        if all(os.path.exists(os.path.join(extracted, dmp))
               for dmp in ('names.dmp', 'nodes.dmp')):
            self.logger.info('Reusing the NCBI taxonomy already extracted to {}.'.format(
                extracted))
            return extracted

        tarball = os.path.join(taxonomy_dir, 'taxdump_{}.tar.gz'.format(date_stamp))
        self._download(TAXDUMP_URL, tarball)
        self._verify_taxonomy(tarball)

        taxdump_dir = extracted
        self.logger.info('Extracting {} to {}'.format(
            os.path.basename(tarball), os.path.basename(taxdump_dir)))
        try:
            extract_tarball(tarball, taxdump_dir)
        except (OSError, tarfile.TarError) as exc:
            self.logger.error('Failed to extract {}: {}'.format(tarball, exc))
            sys.exit()

        self.logger.info('Taxonomy holds {:,} files, including {}.'.format(
            len(os.listdir(taxdump_dir)),
            ', '.join(name for name in ('names.dmp', 'nodes.dmp')
                      if os.path.exists(os.path.join(taxdump_dir, name)))))

        # The archive has served its purpose: its MD5 has been checked against
        # NCBI and its contents are on disk. Keeping it would hold 80 MB of a
        # release for a check that can no longer fail.
        for spent in (tarball, tarball + '.md5'):
            os.remove(spent)
            self.logger.info('Removed {}.'.format(os.path.basename(spent)))

        return taxdump_dir

    def _verify_taxonomy(self, tarball: str) -> None:
        """Check the taxonomy archive against the MD5 NCBI publishes beside it.

        A truncated or corrupted taxdump reads as a valid, smaller taxonomy, and
        the damage only shows up as genomes that mysteriously lost their NCBI
        lineage several commands later.

        Parameters
        ----------
        tarball : str
            Downloaded taxonomy archive.
        """

        md5_file = tarball + '.md5'
        self._download(TAXDUMP_URL + '.md5', md5_file)

        with open(md5_file) as handle:
            # NCBI writes "<md5>  taxdump.tar.gz"
            published = handle.read().split()[0].strip().lower()

        observed = file_checksum(tarball, hashlib.md5())
        if observed != published:
            self.logger.error(
                'MD5 of {} is {}, but NCBI publishes {}; the download is corrupt.'.format(
                    os.path.basename(tarball), observed, published))
            sys.exit()

        self.logger.info('MD5 of {} matches the one published by NCBI.'.format(
            os.path.basename(tarball)))

    def download_assembly_summaries(self) -> Dict[Tuple[str, str], str]:
        """Download the assembly summary file of each database and domain of this group.

        These land in the root of the output directory, which is where
        select_genomes and the taxonomy step below both expect to find them.
        Their names carry the domain rather than the group, so running the other
        group into the same directory adds files beside these rather than
        replacing them.

        @return: dict of (database, domain) to the file downloaded.
        """

        downloaded = {}
        for database, domain, url, name in assembly_summary_downloads(self.group):
            output_file = os.path.join(self.output_dir, name)
            self._download(url, output_file, compress=True)
            downloaded[(database, domain)] = output_file

        return downloaded

    def generate_standardised_taxonomy(self,
                                       taxdump_dir: str,
                                       summaries: Dict[Tuple[str, str], str],
                                       release_number: int) -> str:
        """Produce the standardised NCBI taxonomy of the genomes just downloaded.

        This is the NCBI taxonomy parser, run over the files just downloaded
        rather than over files named by hand, and over this group's files alone.
        Whether the taxonomy keeps NCBI's subranks is the group's: the 7 ranks
        GTDB curates for prokaryotes, the 13 that carry NCBI's intermediate
        ranks for fungi.

        The output prefix is a path, so the files land in
        taxonomy/standardised_taxonomy/ without the working directory being
        changed, and it names the group, so the other group's run into the same
        directory neither overwrites these files nor their filter report.

        Parameters
        ----------
        taxdump_dir : str
            Directory holding the extracted nodes.dmp and names.dmp.
        summaries : dict
            (database, domain) to assembly summary file, as downloaded for this
            group.
        release_number : int
            GTDB release number, which names the output files.

        @return: directory the standardised taxonomy was written to.
        """

        keep_subranks = GROUP_KEEP_SUBRANKS[self.group]

        output_dir = os.path.join(self.taxonomy_dir(), STANDARDISED_TAXONOMY_DIR)
        os.makedirs(output_dir, exist_ok=True)
        output_prefix = os.path.join(output_dir, 'ncbi_r{}_{}'.format(
            release_number, self.file_tag()))

        self.logger.info('Generating {} rank NCBI taxonomy as {}_*.tsv'.format(
            13 if keep_subranks else 7, output_prefix))
        try:
            TaxonomyNCBI().parse_ncbi_taxonomy(
                taxdump_dir,
                [summaries[key] for key in sorted(summaries)],
                keep_subranks,
                output_prefix,
                os.path.join(output_dir, '{}_failed_filters.tsv'.format(self.file_tag())))
        except Exception as exc:
            self.logger.error('Failed to parse the NCBI taxonomy: {}'.format(exc))
            sys.exit()

        self.logger.info('Wrote {:,} file(s) to {}'.format(
            len(os.listdir(output_dir)), output_dir))

        return output_dir

    def run(self, release_number: int) -> None:
        """Download the NCBI metadata of one group and standardise its taxonomy.

        Parameters
        ----------
        release_number : int
            GTDB release number, which names the standardised taxonomy files.
        """

        date_stamp = datetime.date.today().strftime('%Y%m%d')
        self.logger.info('Downloading {} NCBI metadata for release {} to {}'.format(
            self.group, release_number, self.output_dir))

        taxdump_dir = self.download_taxonomy(self.taxonomy_dir(), date_stamp)
        summaries = self.download_assembly_summaries()

        self.logger.info('Downloaded {:,} assembly summary file(s): {}'.format(
            len(summaries), ', '.join(sorted(os.path.basename(f)
                                             for f in summaries.values()))))

        self.generate_standardised_taxonomy(taxdump_dir, summaries, release_number)
