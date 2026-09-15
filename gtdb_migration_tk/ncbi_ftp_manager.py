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
ncbi_ftp_manager.py -- select the NCBI genomes belonging in a GTDB release, and update
the GTDB genome directories to match.

Each GTDB release is built by comparing the genomes held by the previous release
with the genomes NCBI currently offers on its FTP site. Genomes NCBI no longer
offers are removed, genomes new to NCBI are copied across, and genomes held by
both are checked for a changed genomic FASTA since the previous release. The
copying and comparing is done by FTPTools in ncbi_ftp_manager_tools.py; what
lives here is the bookkeeping of which genomes fall into which of those three
groups, done by GenomeManager.

GenomeManager decides nothing beyond what the genome directory files of the
mirror and of the previous release already say. The mirror is built from the
selection, and the selection is where a genome is accepted or passed over, so
every genome the mirror holds is wanted. RefSeq and GenBank are handled in
separate runs of the same code, told apart by the accession prefix (GCF or GCA),
because the update of each is reported separately; nothing else differs.

MetadataSyncManager comes before any of that: it downloads the NCBI taxonomy and
the assembly summary files a release is selected from, and decides nothing.

SelectedGenomesManager answers only the first half of that question, and answers it
from the summary files alone: no previous release, no FTP mirror, no genome
directories. It is what the select_genomes command runs, so the genomes wanted in
a release can be settled and reviewed before any of them is copied.

Which genomes NCBI considers latest, which belong to large multi-isolate projects,
and which are paired with a RefSeq assembly are all read from the NCBI assembly
summary files by ncbi_utils.py, which ncbi_genome_sync.py reads them with too.
"""

import os
import sys
import datetime
import hashlib
import logging
import tarfile
from contextlib import ExitStack
from typing import Dict, List, Set, Tuple

from tqdm import tqdm

from gtdb_migration_tk.biolib_lite.common import canonical_gid
from gtdb_migration_tk.ncbi_ftp_manager_tools import (
    FTPTools, SELECTED_GENOMES_FILE, TAXDUMP_URL, SelectedRow,
    assembly_summary_downloads, download_file, extract_tarball, file_checksum,
    write_selected_genomes)
from gtdb_migration_tk.ncbi_tax_manager import TaxonomyNCBI
from gtdb_migration_tk.ncbi_utils import count_summary_rows, read_assembly_summary
from gtdb_migration_tk.utils.common import count_lines


# Accession prefixes of the two NCBI databases.
REFSEQ_PREFIX = 'GCF'
GENBANK_PREFIX = 'GCA'

# Suffixes of the assembly summary file names of the two NCBI databases. NCBI
# names each file for the database it describes, and the per-domain copies GTDB
# works from keep that suffix: assembly_summary_archaea_refseq.txt,
# assembly_summary_bacteria_genbank.txt. Reading the suffix is what lets each
# file be parsed once, RefSeq files before GenBank files, rather than every file
# being parsed once per database.
REFSEQ_SUFFIX = '_refseq.txt'
GENBANK_SUFFIX = '_genbank.txt'

# Values of excluded_from_refseq marking an assembly as one of the thousands of
# near-identical isolates NCBI sequences for outbreak and pathogen surveillance.
#
# NCBI renamed this annotation: summary files written before the change carry
# 'derived from surveillance project', and those written after it carry 'large
# multi-isolate project', for the same genomes. Both are matched, so a selection
# made from an archived set of summary files excludes what a selection made from
# a current set excludes. The column holds several such annotations separated by
# semicolons, hence the substring test.
MULTI_ISOLATE_TAGS = ('large multi-isolate project', 'surveillance')


# Placeholder for the notes column of the selected genome table, matching the
# convention NCBI uses for an empty field in the summary files themselves.
NO_NOTE = 'na'

# Subdirectories of the root NCBI directory of a release. Everything derived from
# the NCBI taxonomy lives under TAXONOMY_DIR, the standardised 7 rank form
# included, so the taxonomy of a release is one directory to copy or archive.
TAXONOMY_DIR = 'taxonomy'
STANDARDISED_TAXONOMY_DIR = 'standardised_taxonomy'


def has_ftp_path(ftp_path: str) -> bool:
    """Report whether NCBI serves a directory for an assembly.

    NCBI writes 'na' in ftp_path for an assembly it lists but does not serve, and the
    column can be empty in an older file. Such a genome cannot be mirrored, so it is not
    selected: the table select_genomes writes is the list ncbi_genome_sync fetches from,
    and a row with nothing to fetch would be reported as skipped by every run of it
    forever.

    This is deliberately the same test ncbi_genome_sync.read_assembly_summary() applies
    when deciding a row has "no usable ftp_path". The two must agree, or the selection
    would promise genomes the sync then refuses.

    Parameters
    ----------
    ftp_path : str
        Value of the ftp_path column of an assembly summary file.

    @return: True if the assembly has a directory at NCBI.
    """

    return bool(ftp_path) and ftp_path.lower() != 'na'


def is_multi_isolate(excluded_from_refseq: str) -> bool:
    """Report whether an assembly belongs to a large multi-isolate project.

    Parameters
    ----------
    excluded_from_refseq : str
        Value of the excluded_from_refseq column of an assembly summary file.

    @return: True if the assembly is annotated as part of such a project.
    """

    return any(tag in excluded_from_refseq for tag in MULTI_ISOLATE_TAGS)


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


class SelectedGenomesManager:
    """Select the NCBI genomes belonging in a GTDB release, from the summary files alone.

    This is the one decision in a release about which genomes are wanted, and
    it is made without reference to a previous release or to the FTP mirror: the
    NCBI assembly summary files are the only input, and the answer is a table of
    the genomes wanted. Nothing is copied, compared, or deleted, so the selection
    can be produced and reviewed before any part of a release is built; the sync
    then mirrors exactly this table, and GenomeManager updates the release from
    that mirror without deciding anything further.

    Every RefSeq assembly is wanted. A GenBank assembly is wanted only where
    RefSeq does not already cover the genome and the assembly is not part of a
    large multi-isolate project, those projects being the thousands of
    near-identical isolates NCBI excludes from RefSeq and which GTDB does not
    want either.

    A genome NCBI lists but does not serve -- ftp_path 'na' -- is not selected
    from either database: there is nothing to mirror, and the table this writes
    is what ncbi_genome_sync fetches from. That cuts both ways, and the second
    way is the point: RefSeq covers a genome only if the RefSeq assembly was
    itself selected, so a GenBank assembly whose RefSeq counterpart has no
    ftp_path is NOT covered and is selected in its place. The alternative would
    drop the genome from the release entirely -- NCBI holds it, under its
    GenBank accession, and nothing else would carry it.

    Two assumptions the rule rests on are verified rather than trusted, since a
    summary file that breaks either yields a plausible looking but wrong
    selection: that NCBI flags no RefSeq assembly as a large multi-isolate
    project, and that every assembly listed is the latest version of itself. A
    genome breaking either is reported and left out of the selection; neither is
    treated as fatal, so a single odd row cannot cost a whole run.

    Which database a summary file describes is read from its name, so each file
    is parsed once: the RefSeq files first, since a GenBank assembly cannot be
    judged until the genomes RefSeq covers are known, then the GenBank files.

    The table it writes carries the four columns ncbi_genome_sync reads, so the
    selection is directly the list of genomes to mirror: nothing has to be joined
    back against the NCBI summary files to sync what was chosen.

A third discrepancy is recorded but not acted on. NCBI regularly pairs a
    GenBank assembly with a RefSeq assembly that its own summary files do not
    list, which is an error in those files rather than anything about the genome.
    The genome is selected, since nothing else covers it, and the pairing is
    written to the notes column of its row: reporting each one in the log would
    bury the failures that do need attention. A GenBank assembly taken in place
    of a RefSeq assembly that WAS listed but not selected carries a note too, and
    a different one, naming the reason the RefSeq assembly was passed over --
    those two look identical in the table otherwise, and only one of them is
    NCBI's mistake.
    """

    def __init__(self, output_dir: str) -> None:
        """Record where the selection and its log are to be written.

        Parameters
        ----------
        output_dir : str
            Output directory for the selected genome table.
        """

        self.output_dir = output_dir
        self.logger = logging.getLogger('timestamp')

    def _read_summary(self, assembly_summary: str):
        """Read the fields of an assembly summary file needed to select genomes.

        Parameters
        ----------
        assembly_summary : str
            NCBI assembly summary file.

        @return: iterator over (accession, version_status, gbrs_paired_asm,
            excluded_from_refseq, ftp_path), one tuple per genome.
        """

        records = read_assembly_summary(assembly_summary,
                                        'assembly_accession',
                                        'version_status',
                                        'gbrs_paired_asm',
                                        'excluded_from_refseq',
                                        'ftp_path')

        return tqdm(records,
                    total=count_summary_rows(assembly_summary),
                    desc='Selecting genomes from {}'.format(
                        os.path.basename(assembly_summary)))

    def _is_latest(self, accession: str, version_status: str) -> bool:
        """Report whether NCBI considers an assembly the latest version of itself.

        A superseded assembly in the input is a sign the summary files describe
        something other than the genomes of the new release, so it is reported
        rather than quietly passed over.

        Parameters
        ----------
        accession : str
            Assembly accession of the genome.
        version_status : str
            Value of the version_status column.

        @return: True if the assembly is flagged as the latest version.
        """

        if version_status == 'latest':
            return True

        self.logger.warning(
            '{} has version_status "{}" rather than "latest" and was not selected.'.format(
                accession, version_status or 'na'))

        return False

    def group_by_database(self, new_list_genomes: List[str]) -> Tuple[List[str], List[str]]:
        """Split the assembly summary files into the RefSeq files and the GenBank files.

        A file whose name says neither is fatal rather than ignored: skipping it
        would leave every genome it lists out of the release, which looks exactly
        like a successful run of a smaller release.

        Parameters
        ----------
        new_list_genomes : list of str
            NCBI assembly summary files describing the new release.

        @return: tuple of the RefSeq files and the GenBank files.
        """

        refseq_files, genbank_files = [], []

        for assembly_summary in new_list_genomes:
            # the suffix is what names the database, and sits behind the .gz on a
            # compressed file; ncbi_metadata_sync writes these gzipped, older
            # releases hold them plain, and both are read the same way
            name = os.path.basename(assembly_summary)
            if name.endswith('.gz'):
                name = name[:-len('.gz')]

            if name.endswith(REFSEQ_SUFFIX):
                refseq_files.append(assembly_summary)
            elif name.endswith(GENBANK_SUFFIX):
                genbank_files.append(assembly_summary)
            else:
                self.logger.error(
                    'Cannot tell which NCBI database {} describes: the name of an assembly '
                    'summary file must end with {} or {}, optionally gzipped.'.format(
                        os.path.basename(assembly_summary), REFSEQ_SUFFIX, GENBANK_SUFFIX))
                sys.exit()

        # either list being empty is legal, but it is far more often a mistake in
        # the file list than a release genuinely drawn from one database
        for database, files in (('RefSeq', refseq_files), ('GenBank', genbank_files)):
            if not files:
                self.logger.warning(
                    'No {} assembly summary files were given; the selection will '
                    'contain no {} genomes.'.format(database, database))

        return refseq_files, genbank_files

    def select_refseq(self,
                      refseq_files: List[str]
                      ) -> Tuple[List[SelectedRow], Set[str], Dict[str, str]]:
        """Select the RefSeq genomes of the new release.

        Every RefSeq assembly is wanted, so this is a matter of reading the
        summary files. NCBI is not expected to flag any RefSeq assembly as a
        large multi-isolate project; one that is flagged is reported and
        dropped, as it would otherwise enter the release by the very route the
        GenBank rule exists to close.

        This runs before select_genbank, which needs the genomes RefSeq covers.

        Parameters
        ----------
        refseq_files : list of str
            NCBI RefSeq assembly summary files describing the new release.

        @return: tuple of the selected rows (accession, ftp_path, version_status,
            excluded_from_refseq, gbrs_paired_asm, notes), the canonical accessions
            they cover, and the canonical accession of every RefSeq genome that was
            read but NOT selected, mapped to why. select_genbank needs both: the
            first says which genomes need no GenBank copy, the second explains, in
            the notes column, why a genome that appears to have a RefSeq assembly
            was taken from GenBank anyway.
        """

        selected = []
        covered = set()
        dropped = {}
        read = 0
        misfiled = 0
        no_ftp_path = 0

        for assembly_summary in refseq_files:
            for accession, version_status, paired_asm, excluded, ftp_path in self._read_summary(
                    assembly_summary):

                if not accession.startswith(REFSEQ_PREFIX):
                    # the file is named for RefSeq but this row is not a RefSeq
                    # assembly, so the file list says something the table does not
                    misfiled += 1
                    continue

                read += 1

                gid = canonical_gid(accession)

                if not self._is_latest(accession, version_status):
                    dropped[gid] = 'version_status is "{}"'.format(version_status or 'na')
                    continue

                if is_multi_isolate(excluded):
                    self.logger.warning(
                        '[Unexpected] RefSeq genome {} is annotated as "{}" and was not '
                        'selected.'.format(accession, excluded))
                    dropped[gid] = 'annotated as "{}"'.format(excluded)
                    continue

                if not has_ftp_path(ftp_path):
                    # NCBI lists it but serves no directory, so there is nothing to
                    # mirror. Recorded in `dropped` so that its GenBank counterpart,
                    # which may well have a directory, is selected in its place.
                    no_ftp_path += 1
                    dropped[gid] = 'NCBI serves no directory for it'
                    continue

                selected.append((accession, ftp_path, version_status, excluded,
                                 paired_asm, NO_NOTE))
                covered.add(gid)

        self.logger.info('RefSeq: selected {:,} of {:,} genomes.'.format(len(selected), read))
        if no_ftp_path:
            self.logger.warning('RefSeq: {:,} genomes ignored, NCBI serves no directory for '
                                'them (ftp_path is na); where a GenBank assembly exists it is '
                                'selected instead.'.format(no_ftp_path))
        if misfiled:
            self.logger.warning('RefSeq: {:,} rows in the RefSeq assembly summary files are '
                                'not RefSeq assemblies and were ignored.'.format(misfiled))

        return selected, covered, dropped

    def select_genbank(self,
                       genbank_files: List[str],
                       covered: Set[str],
                       dropped: Dict[str, str]) -> List[SelectedRow]:
        """Select the GenBank genomes needed to supplement the RefSeq genomes.

        A GenBank assembly is wanted when RefSeq does not already cover the
        genome and the assembly is not part of a large multi-isolate project.
        Coverage is decided by canonical accession against the RefSeq genomes
        actually selected above, not by the gbrs_paired_asm field, so a GenBank
        assembly paired with a RefSeq assembly that is absent from the new
        release is still selected -- including when that assembly was read and
        passed over because NCBI serves no directory for it, which is the case
        this rule exists to catch. Where the two disagree NCBI describes a RefSeq
        assembly its own summary files do not list; that is common enough to be
        recorded in the notes column of the affected row rather than reported
        genome by genome.

        Parameters
        ----------
        genbank_files : list of str
            NCBI GenBank assembly summary files describing the new release.
        covered : set
            Canonical accessions of the selected RefSeq genomes.
        dropped : dict
            Canonical accession to the reason its RefSeq assembly was not selected.

        @return: list of selected rows (accession, ftp_path, version_status,
            excluded_from_refseq, gbrs_paired_asm, notes).
        """

        selected = []
        read = 0
        multi_isolate = 0
        unpaired = 0
        misfiled = 0
        no_ftp_path = 0
        rescued = 0

        for assembly_summary in genbank_files:
            for accession, version_status, paired_asm, excluded, ftp_path in self._read_summary(
                    assembly_summary):

                if not accession.startswith(GENBANK_PREFIX):
                    misfiled += 1
                    continue

                read += 1

                gid = canonical_gid(accession)

                if not self._is_latest(accession, version_status):
                    continue

                if gid in covered:
                    continue

                if not has_ftp_path(ftp_path):
                    # nothing to mirror. Unlike the RefSeq case this rescues
                    # nothing: RefSeq was already asked first and did not cover it
                    no_ftp_path += 1
                    continue

                if is_multi_isolate(excluded):
                    multi_isolate += 1
                    continue

                note = NO_NOTE
                if paired_asm.startswith(REFSEQ_PREFIX):
                    # The row appears to have a RefSeq counterpart, yet here it is
                    # being taken from GenBank; the note says why, and the two
                    # reasons are not the same kind of thing.
                    reason = dropped.get(gid)
                    if reason is None:
                        # NCBI names a RefSeq assembly that its own summary files
                        # do not list: an error in those files, nothing about the
                        # genome
                        unpaired += 1
                        note = ('paired RefSeq assembly {} absent from the assembly '
                                'summary files'.format(paired_asm))
                    else:
                        # the RefSeq assembly was listed and passed over, so the
                        # GenBank copy is the only one this release can carry
                        rescued += 1
                        note = ('paired RefSeq assembly {} was not selected: {}'
                                .format(paired_asm, reason))

                selected.append((accession, ftp_path, version_status, excluded,
                                 paired_asm, note))

        self.logger.info('GenBank: selected {:,} of {:,} genomes.'.format(len(selected), read))
        self.logger.info('GenBank: {:,} genomes skipped as part of a large multi-isolate '
                         'project.'.format(multi_isolate))
        if no_ftp_path:
            self.logger.warning('GenBank: {:,} genomes ignored, NCBI serves no directory for '
                                'them (ftp_path is na).'.format(no_ftp_path))
        if rescued:
            self.logger.info('GenBank: {:,} genomes selected in place of a RefSeq assembly '
                             'that was read but not selected; the reason is in the notes '
                             'column of each.'.format(rescued))
        if misfiled:
            self.logger.warning('GenBank: {:,} rows in the GenBank assembly summary files are '
                                'not GenBank assemblies and were ignored.'.format(misfiled))
        if unpaired:
            # common enough that it is reported as a count rather than per genome:
            # the affected rows carry a note, and the fault is NCBI's
            self.logger.info('GenBank: {:,} selected genomes name a paired RefSeq assembly '
                             'absent from the assembly summary files; each is noted in the '
                             'table.'.format(unpaired))

        return selected

    def run(self, new_list_genomes: List[str]) -> None:
        """Write the table of genomes selected for the new release.

        Parameters
        ----------
        new_list_genomes : list of str
            NCBI assembly summary files describing the new release.
        """

        refseq_files, genbank_files = self.group_by_database(new_list_genomes)
        self.logger.info('Reading {:,} RefSeq and {:,} GenBank assembly summary file(s).'.format(
            len(refseq_files), len(genbank_files)))

        refseq, covered, dropped = self.select_refseq(refseq_files)
        genbank = self.select_genbank(genbank_files, covered, dropped)

        output_file = os.path.join(self.output_dir, SELECTED_GENOMES_FILE)
        write_selected_genomes(refseq + genbank, output_file)

        self.logger.info('Selected {:,} genomes for the new release: {}'.format(
            len(refseq) + len(genbank), output_file))


class MetadataSyncManager:
    """Download the NCBI metadata a GTDB release is built from.

    A release starts from two things NCBI publishes and GTDB only reads: the
    taxonomy, and the assembly summary files describing every assembly NCBI
    holds. This fetches both into one directory, which is then the input to
    select_genomes and, later, to the commands that attach NCBI taxonomy to the
    genomes chosen.

    This is the procedure the GTDB wiki gives as "Download latest NCBI taxonomy",
    the first step of "Download the latest RefSeq and GenBank assembly data", and
    "Generate 7-rank NCBI taxonomy", run in one go. The output directory is the
    root NCBI directory of a release, and is laid out as the wiki leaves it:

        <output_dir>/assembly_summary_<domain>_<database>.txt.gz
        <output_dir>/taxonomy/taxdump_<date>/
        <output_dir>/taxonomy/standardised_taxonomy/ncbi_r<release>_*.tsv

    Three things differ from doing it by hand. NCBI names every assembly summary
    file assembly_summary.txt, so the database and domain are put back into the
    name as they are saved, and they are gzipped on the way in rather than kept
    as the 1.8 GB of text NCBI serves. The taxonomy dump is checked against the
    MD5 NCBI publishes beside it, since a truncated taxdump is not obviously
    broken until a release has been built on it, and the archive is then
    discarded: it has been verified and unpacked, and nothing reads it again.

    The fungal assembly summaries are downloaded alongside the prokaryotic ones,
    so that a release holds the table NCBI was serving when it was built. Fungal
    genomes are otherwise a separate procedure: they are not given to
    select_genomes, and the 7 rank taxonomy below is built from archaea and
    bacteria only.

    Nothing here decides anything: unlike the other managers in this module it
    only fetches, and hands what it fetched to the NCBI taxonomy parser.
    """

    def __init__(self, output_dir: str) -> None:
        """Record where the metadata and its log are to be written.

        Parameters
        ----------
        output_dir : str
            Output directory for the downloaded NCBI metadata.
        """

        self.output_dir = output_dir
        self.logger = logging.getLogger('timestamp')

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
        """Download and extract the NCBI taxonomy.

        The archive is verified against NCBI's published MD5, unpacked, and then
        removed: only the extracted directory is ever read again.

        Parameters
        ----------
        taxonomy_dir : str
            Directory to download the taxonomy into.
        date_stamp : str
            Date the download was made, as YYYYMMDD.

        @return: directory the taxonomy was extracted into.
        """

        os.makedirs(taxonomy_dir, exist_ok=True)

        tarball = os.path.join(taxonomy_dir, 'taxdump_{}.tar.gz'.format(date_stamp))
        self._download(TAXDUMP_URL, tarball)
        self._verify_taxonomy(tarball)

        taxdump_dir = os.path.join(taxonomy_dir, 'taxdump_{}'.format(date_stamp))
        self.logger.info('Extracting {} to {}'.format(
            os.path.basename(tarball), os.path.basename(taxdump_dir)))
        try:
            extract_tarball(tarball, taxdump_dir)
        except (OSError, tarfile.TarError) as exc:
            self.logger.error('Failed to extract {}: {}'.format(tarball, exc))
            sys.exit()

        self.logger.info('Taxonomy holds {} files, including {}.'.format(
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
        """Download the assembly summary file of each database and domain.

        These land in the root of the output directory, which is where
        select_genomes and the taxonomy step below both expect to find them. The
        fungal summaries land there too, and are simply not among the keys
        either of those steps asks for.

        @return: dict of (database, domain) to the file downloaded.
        """

        downloaded = {}
        for database, domain, url, name in assembly_summary_downloads():
            output_file = os.path.join(self.output_dir, name)
            self._download(url, output_file, compress=True)
            downloaded[(database, domain)] = output_file

        return downloaded

    def generate_standardised_taxonomy(self,
                                       taxdump_dir: str,
                                       summaries: Dict[Tuple[str, str], str],
                                       release_number: int) -> str:
        """Produce the 7 rank NCBI taxonomy of the genomes NCBI holds.

        This is the NCBI taxonomy parser, run over the files just downloaded
        rather than over files named by hand. The parser takes the four
        prokaryotic summaries by name, so the fungal ones alongside them are
        left out of the taxonomy rather than filtered out of it. The output
        prefix is a path, so the files land in taxonomy/standardised_taxonomy/
        without the working directory being changed.

        Parameters
        ----------
        taxdump_dir : str
            Directory holding the extracted nodes.dmp and names.dmp.
        summaries : dict
            (database, domain) to assembly summary file; the four prokaryotic
            entries are used.
        release_number : int
            GTDB release number, which names the output files.

        @return: directory the standardised taxonomy was written to.
        """

        output_dir = os.path.join(self.taxonomy_dir(), STANDARDISED_TAXONOMY_DIR)
        os.makedirs(output_dir, exist_ok=True)
        output_prefix = os.path.join(output_dir, 'ncbi_r{}'.format(release_number))

        self.logger.info('Generating 7 rank NCBI taxonomy as {}_*.tsv'.format(output_prefix))
        try:
            TaxonomyNCBI().parse_ncbi_taxonomy(
                taxdump_dir,
                summaries[('refseq', 'archaea')],
                summaries[('refseq', 'bacteria')],
                summaries[('genbank', 'archaea')],
                summaries[('genbank', 'bacteria')],
                False,                           # subranks are not kept, as the wiki has it
                output_prefix)
        except Exception as exc:
            self.logger.error('Failed to parse the NCBI taxonomy: {}'.format(exc))
            sys.exit()

        self.logger.info('Wrote {} file(s) to {}'.format(
            len(os.listdir(output_dir)), output_dir))

        return output_dir

    def run(self, release_number: int) -> None:
        """Download the NCBI metadata of a release and standardise its taxonomy.

        Parameters
        ----------
        release_number : int
            GTDB release number, which names the standardised taxonomy files.
        """

        date_stamp = datetime.date.today().strftime('%Y%m%d')
        self.logger.info('Downloading NCBI metadata for release {} to {}'.format(
            release_number, self.output_dir))

        taxdump_dir = self.download_taxonomy(self.taxonomy_dir(), date_stamp)
        summaries = self.download_assembly_summaries()

        self.logger.info('Downloaded {} assembly summary file(s): {}'.format(
            len(summaries), ', '.join(sorted(os.path.basename(f)
                                             for f in summaries.values()))))

        self.generate_standardised_taxonomy(taxdump_dir, summaries, release_number)
