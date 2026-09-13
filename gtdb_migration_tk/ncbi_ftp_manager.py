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
both are checked for files that have changed since the previous release. The
copying, deleting, and checksumming is done by FTPTools in ncbi_ftp_manager_tools.py;
what lives here is the decision about which genomes are wanted in the first place.

RefSeq and GenBank are handled separately, by RefSeqManager and GenBankManager,
because that decision differs between them. Every RefSeq assembly flagged by NCBI
as the latest version is wanted. GenBank is only consulted where RefSeq falls
short, so a GenBank assembly is wanted when it has no RefSeq counterpart, or when
that counterpart is missing from the FTP site or holds no genome assembly. Each
GenBank decision is recorded in gca_selection.log, as the selection is the part of
a release that is hardest to reconstruct after the fact.

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
from gtdb_migration_tk.ncbi_utils import genome_assembly_file, read_assembly_summary
from gtdb_migration_tk.utils.common import count_lines


# Domain labels held in the genome domain map passed to FTPTools.
ARCHAEA = 'Archaea'
BACTERIA = 'Bacteria'

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


def is_multi_isolate(excluded_from_refseq: str) -> bool:
    """Report whether an assembly belongs to a large multi-isolate project.

    Parameters
    ----------
    excluded_from_refseq : str
        Value of the excluded_from_refseq column of an assembly summary file.

    @return: True if the assembly is annotated as part of such a project.
    """

    return any(tag in excluded_from_refseq for tag in MULTI_ISOLATE_TAGS)


class GenericManager:
    """Base class for comparing a GTDB release to the genomes held at NCBI.

    Subclasses (RefSeqManager, GenBankManager) are responsible for deciding
    which accessions at NCBI are of interest for the new release. This class
    provides the operations common to both: reading the genomes of the previous
    release and of the FTP site, and determining which genomes must be removed,
    added, or checked for updates relative to that release. Genomes are tracked
    as dictionaries mapping an accession to its genome directory.
    """

    def __init__(self):
        """Initialise the domain map shared with the FTP tools."""

        self.genome_domain_dict = {}
        self.logger = logging.getLogger('timestamp')

    def load_previous_records(self, old_genome_dirs: str) -> Dict[str, str]:
        """Read the genome directory file of the previous GTDB release.

        Parameters
        ----------
        old_genome_dirs : str
            Genome directory file (accession, path) for the previous release.

        @return: dict of accession to genome directory for the previous release.
        """

        with open(old_genome_dirs, 'r') as old_file:
            old_genomes = {}
            for old_line in old_file:
                accession, path, *_ = old_line.split('\t')
                old_genomes[accession] = path.strip()

        self.logger.info('Previous release: {:,} genomes.'.format(len(old_genomes)))

        return old_genomes

    def load_ftp_records(self,
                         ftp_genome_dirs: str,
                         accession_prefix: str,
                         accessions: Set[str]) -> Dict[str, str]:
        """Read the genome directory file of the NCBI FTP site.

        Only genomes selected for the new release are retained, so this is the
        intersection of the genomes held on the FTP site with the accessions
        identified from the NCBI assembly summary files.

        Parameters
        ----------
        ftp_genome_dirs : str
            Genome directory file (accession, path) for the FTP mirror.
        accession_prefix : str
            Accession prefix of the database of interest, i.e. GCF or GCA.
        accessions : set
            Accessions selected for the new release.

        @return: dict of accession to genome directory for genomes on the FTP site.
        """

        self.logger.info('Reading genomes on the FTP site.')

        new_genomes = {}
        with open(ftp_genome_dirs, 'r') as new_genome_dirs_file:
            for new_line in tqdm(new_genome_dirs_file, total=count_lines(ftp_genome_dirs)):
                gid, path, *_ = new_line.split('\t')
                if gid.startswith(accession_prefix) and gid in accessions:
                    new_genomes[gid] = path.strip()

        self.logger.info('FTP site: {:,} genomes.'.format(len(new_genomes)))

        return new_genomes

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

        self.logger.info('Remove Genome Step')
        removed_dict = {gid: old_genomes[gid]
                        for gid in old_genomes.keys() - new_genomes.keys()}
        self.logger.info('{0:,} genomes to remove'.format(len(removed_dict)))

        return removed_dict

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

        self.logger.info('Add Genome Step')
        added_dict = {gid: new_genomes[gid]
                      for gid in new_genomes.keys() - old_genomes.keys()}
        self.logger.info('{0:,} genomes to add'.format(len(added_dict)))

        return added_dict

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

        self.logger.info('Update Genome Step')
        intersect_list = list(old_genomes.keys() & new_genomes.keys())
        self.logger.info('{:,} genomes to compare'.format(len(intersect_list)))

        return intersect_list


class RefSeqManager(GenericManager):
    """Update the GTDB copy of RefSeq (GCF) genomes from the NCBI FTP site.

    All RefSeq assemblies flagged as the latest version are of interest, so
    genome selection is simply a matter of reading the NCBI assembly summary
    files. Genomes are then added, removed, or refreshed relative to the
    previous GTDB release.
    """

    def __init__(self,
                 new_refseq_dir: str,
                 dry_run: bool = False,
                 cpus: int = 1) -> None:
        """Record where the new release and its reports are to be written.

        Parameters
        ----------
        new_refseq_dir : str
            Output directory for the new release, where reports are written.
        dry_run : bool
            Report the changes that would be made without modifying any files.
        cpus : int
            Number of processes used when comparing genomes.
        """

        super().__init__()
        self.new_refseq_dir = new_refseq_dir
        self.dry_run = dry_run
        self.cpus = cpus

    def parse_assembly_summary(self, assembly_summary: str) -> List[str]:
        """Identify the latest assembly version of each RefSeq genome.

        Parameters
        ----------
        assembly_summary : str
            NCBI assembly summary file for a single domain.

        @return: list of accessions flagged by NCBI as the latest version.
        """

        records = read_assembly_summary(assembly_summary,
                                        'assembly_accession',
                                        'version_status')

        return [accession for accession, version_status in records
                if version_status == 'latest']

    def run_comparison(self,
                       ftp_refseq: str,
                       new_refseq: str,
                       ftp_genome_dirs: str,
                       old_genome_dirs: str,
                       archaea_assembly_summary: str,
                       bacteria_assembly_summary: str) -> None:
        """Update the GTDB genome directories to reflect the current RefSeq holdings.

        Genomes on the FTP site that are flagged as the latest version are
        compared to the previous GTDB release. Genomes no longer at NCBI are
        removed, new genomes are copied across, and genomes common to both are
        checked for changed files. Only directories containing
        latest_assembly_versions, and not _assembly_structure, are of interest.

        Parameters
        ----------
        ftp_refseq : str
            Local mirror of the RefSeq portion of the NCBI FTP site.
        new_refseq : str
            Output directory for the new release.
        ftp_genome_dirs : str
            Genome directory file (accession, path) for the FTP mirror.
        old_genome_dirs : str
            Genome directory file (accession, path) for the previous release.
        archaea_assembly_summary : str
            NCBI assembly summary file for archaeal genomes.
        bacteria_assembly_summary : str
            NCBI assembly summary file for bacterial genomes.
        """

        # reports are opened for the duration of the update so they are closed,
        # and their contents kept, if the update fails part way through
        with ExitStack() as reports:
            self.report_gcf = reports.enter_context(
                open(os.path.join(self.new_refseq_dir, 'report_gcf.log'), 'w', 1))
            self.genomes_to_review = reports.enter_context(
                open(os.path.join(self.new_refseq_dir, 'gid_to_review.log'), 'w', 1))

            old_genomes = self.load_previous_records(old_genome_dirs)

            # all genomes flagged by NCBI as the latest assembly version are of interest
            for assembly_summary, domain in ((archaea_assembly_summary, ARCHAEA),
                                             (bacteria_assembly_summary, BACTERIA)):
                for accession in self.parse_assembly_summary(assembly_summary):
                    self.genome_domain_dict[accession] = domain
            self.logger.info('NCBI: {:,} latest assemblies.'.format(len(self.genome_domain_dict)))

            new_genomes = self.load_ftp_records(ftp_genome_dirs,
                                                'GCF',
                                                set(self.genome_domain_dict))

            ftptools = FTPTools(self.report_gcf,
                                self.genomes_to_review,
                                self.genome_domain_dict,
                                self.dry_run)

            # delete genomes from the Database
            removed_dict = self.generate_genomes_to_remove(new_genomes, old_genomes)
            ftptools.remove_genomes(removed_dict)

            # new genomes in FTP
            added_dict = self.generate_genomes_to_add(new_genomes, old_genomes)
            ftptools.add_genomes(added_dict, ftp_refseq, new_refseq, self.genome_domain_dict)

            intersect_list = self.generate_genomes_to_compare(new_genomes, old_genomes)
            ftptools.compare_genomes(intersect_list, old_genomes, new_genomes,
                                     ftp_refseq, new_refseq, self.cpus)


class GenBankManager(GenericManager):
    """Update the GTDB copy of GenBank (GCA) genomes from the NCBI FTP site.

    RefSeq is preferred over GenBank, so only a subset of GenBank assemblies
    are of interest: those with no RefSeq counterpart, and those whose RefSeq
    counterpart is missing or incomplete on the FTP site. Selection decisions
    are recorded in the gca_selection.log report.
    """

    def __init__(self,
                 new_genbank_dir: str,
                 dry_run: bool = False,
                 cpus: int = 1) -> None:
        """Record where the new release and its reports are to be written.

        Parameters
        ----------
        new_genbank_dir : str
            Output directory for the new release, where reports are written.
        dry_run : bool
            Report the changes that would be made without modifying any files.
        cpus : int
            Number of processes used when comparing genomes.
        """

        super().__init__()
        self.new_genbank_dir = new_genbank_dir
        self.dry_run = dry_run
        self.cpus = cpus

    def select_genbank_genomes(self,
                               gbk_arc_assembly: str,
                               gbk_bac_assembly: str,
                               new_refseq_genome_dirs: str) -> List[str]:
        """Identify GenBank genomes required to supplement the RefSeq genomes.

        A GenBank assembly is retained when it is flagged by NCBI as the latest
        version, is not a surveillance genome, and either has no RefSeq
        counterpart, or has a counterpart that is absent from the FTP site. Each 
        decision is written to gca_selection.log.

        Parameters
        ----------
        gbk_arc_assembly : str
            NCBI GenBank assembly summary file for archaeal genomes.
        gbk_bac_assembly : str
            NCBI GenBank assembly summary file for bacterial genomes.
        new_refseq_genome_dirs : str
            Genome directory file for the RefSeq genomes in the new release.

        @return: list of GenBank accessions to include in the new release.
        """

        selected_gca = []
        refseq_dirs = self._populate_genomes_dict(new_refseq_genome_dirs)
        self.logger.info('Indexed {:,} RefSeq genome directories.'.format(len(refseq_dirs)))

        for domain, assembly_file in ((ARCHAEA, gbk_arc_assembly),
                                      (BACTERIA, gbk_bac_assembly)):
            records = read_assembly_summary(assembly_file,
                                            'assembly_accession',
                                            'version_status',
                                            'gbrs_paired_asm',
                                            'excluded_from_refseq')

            for gca_accession, version_status, paired_asm, excluded_from_refseq in tqdm(
                    records,
                    total=count_lines(assembly_file),
                    desc='Selecting {} genomes'.format(domain.lower())):

                if 'surveillance' in excluded_from_refseq:
                    continue
                elif version_status == 'latest':
                    paired_gcf = canonical_gid(paired_asm)

                    if paired_asm.startswith('GCF'):
                        if paired_gcf in refseq_dirs:
                            # if the RefSeq directory does not contain a genome assembly, we copy the GenBank directory instead
                            if not os.path.exists(genome_assembly_file(refseq_dirs[paired_gcf])):
                                self.select_gca.write(
                                    '[Unexpected] {0} associated with {1}: {1} missed files in NCBI FTP directory\n'.format(gca_accession, paired_gcf))
                                selected_gca.append(gca_accession)
                        else:
                            # if the RefSeq directory is not present, we copy the GenBank directory instead
                            self.select_gca.write(
                                '[Unexpected] {0} associated with {1}: {1} not present in NCBI FTP directory\n'.format(
                                    gca_accession, paired_gcf))
                            selected_gca.append(gca_accession)
                    else:
                        if canonical_gid(gca_accession) in refseq_dirs:
                            # THIS SEEMS LIKE A LOGICAL ERROR SINCE WE DON'T REMOVE GENOMES FROM THE NCBI DIRECTORY!
                            self.select_gca.write(
                                '[Unexpected - Logical Error?] {0} skipped because {1} in RefSeq (although {0} has no paired assembly)\n'.format(
                                    gca_accession, paired_gcf))
                            continue
                        else:
                            # this is the expected case where the GenBank assembly does not have a paired RefSeq assembly
                            # so we must use the GenBank assembly
                            selected_gca.append(gca_accession)

                self.genome_domain_dict[gca_accession] = domain

        self.logger.info('Selected {:,} GenBank genomes.'.format(len(selected_gca)))

        return selected_gca

    def run_comparison(self,
                       ftp_genbank: str,
                       new_genbank: str,
                       ftp_genbank_genome_dirs: str,
                       old_genbank_genome_dirs: str,
                       new_refseq_genome_dirs: str,
                       gbk_arc_assembly: str,
                       gbk_bac_assembly: str) -> None:
        """Update the GTDB genome directories to reflect the current GenBank holdings.

        GenBank genomes not already covered by RefSeq are compared to the
        previous GTDB release. Genomes no longer at NCBI are removed, new
        genomes are copied across, and genomes common to both are checked for
        changed files. Only directories containing latest_assembly_versions,
        and not _assembly_structure, are of interest.

        Parameters
        ----------
        ftp_genbank : str
            Local mirror of the GenBank portion of the NCBI FTP site.
        new_genbank : str
            Output directory for the new release.
        ftp_genbank_genome_dirs : str
            Genome directory file (accession, path) for the FTP mirror.
        old_genbank_genome_dirs : str
            Genome directory file (accession, path) for the previous release.
        new_refseq_genome_dirs : str
            Genome directory file for the RefSeq genomes in the new release.
        gbk_arc_assembly : str
            NCBI GenBank assembly summary file for archaeal genomes.
        gbk_bac_assembly : str
            NCBI GenBank assembly summary file for bacterial genomes.
        """

        # reports are opened for the duration of the update so they are closed,
        # and their contents kept, if the update fails part way through
        with ExitStack() as reports:
            self.report = reports.enter_context(
                open(os.path.join(self.new_genbank_dir, 'extra_gbk_report_gcf.log'), 'w', 1))
            self.genomes_to_review = reports.enter_context(
                open(os.path.join(self.new_genbank_dir, 'gcaid_to_review.log'), 'w', 1))
            self.select_gca = reports.enter_context(
                open(os.path.join(self.new_genbank_dir, 'gca_selection.log'), 'w', 1))

            old_genomes = self.load_previous_records(old_genbank_genome_dirs)

            selected_gca = self.select_genbank_genomes(gbk_arc_assembly,
                                                       gbk_bac_assembly,
                                                       new_refseq_genome_dirs)

            new_genomes = self.load_ftp_records(ftp_genbank_genome_dirs,
                                                'GCA',
                                                set(selected_gca))

            ftptools = FTPTools(self.report,
                                self.genomes_to_review,
                                self.genome_domain_dict,
                                self.dry_run)

            # delete genomes from the Database
            removed_dict = self.generate_genomes_to_remove(new_genomes, old_genomes)
            ftptools.remove_genomes(removed_dict)

            # new genomes in FTP
            added_dict = self.generate_genomes_to_add(new_genomes, old_genomes)
            ftptools.add_genomes(added_dict, ftp_genbank, new_genbank, self.genome_domain_dict)

            intersect_list = self.generate_genomes_to_compare(new_genomes, old_genomes)
            ftptools.compare_genomes(intersect_list, old_genomes, new_genomes,
                                     ftp_genbank, new_genbank, self.cpus)

    def _populate_genomes_dict(self, genome_dirs_file: str) -> Dict[str, str]:
        """Index genome directories by their canonical accession.

        NCBI gives paired GenBank and RefSeq assemblies the same 9 digit number,
        so reducing an accession to its canonical form (GCF_005435135.1 and
        GCA_005435135.1 both become G005435135) is what allows a GenBank genome
        to be matched to its RefSeq counterpart, and vice versa.

        Parameters
        ----------
        genome_dirs_file : str
            Genome directory file (accession, path).

        @return: dict of canonical accession to genome directory.
        """

        genome_dirs = {}
        with open(genome_dirs_file, 'r') as list_dirs:
            for line in list_dirs:
                accession, path, *_ = line.split('\t')
                genome_dirs[canonical_gid(accession)] = path.rstrip()

        return genome_dirs


class SelectedGenomesManager:
    """Select the NCBI genomes belonging in a GTDB release, from the summary files alone.

    This answers the same question as RefSeqManager and GenBankManager, but
    without reference to a previous release or to the FTP mirror: the NCBI
    assembly summary files are the only input, and the answer is a table of the
    genomes wanted. Nothing is copied, compared, or deleted, so the selection can
    be produced and reviewed before any part of a release is built.

    The rule is the one GenBankManager applies, reduced to what a summary file
    alone can support. Every RefSeq assembly is wanted. A GenBank assembly is
    wanted only where RefSeq does not already cover the genome and the assembly
    is not part of a large multi-isolate project, those projects being the
    thousands of near-identical isolates NCBI excludes from RefSeq and which
    GTDB does not want either. The checks GenBankManager makes against the
    contents of a RefSeq genome directory have no counterpart here, as no
    directory is consulted.

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
    bury the failures that do need attention.
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
                    total=count_lines(assembly_summary),
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
                      refseq_files: List[str]) -> Tuple[List[SelectedRow], Set[str]]:
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
            excluded_from_refseq, gbrs_paired_asm, notes) and the canonical
            accessions they cover.
        """

        selected = []
        covered = set()
        read = 0
        misfiled = 0

        for assembly_summary in refseq_files:
            for accession, version_status, paired_asm, excluded, ftp_path in self._read_summary(
                    assembly_summary):

                if not accession.startswith(REFSEQ_PREFIX):
                    # the file is named for RefSeq but this row is not a RefSeq
                    # assembly, so the file list says something the table does not
                    misfiled += 1
                    continue

                read += 1

                if not self._is_latest(accession, version_status):
                    continue

                if is_multi_isolate(excluded):
                    self.logger.warning(
                        '[Unexpected] RefSeq genome {} is annotated as "{}" and was not '
                        'selected.'.format(accession, excluded))
                    continue

                selected.append((accession, ftp_path, version_status, excluded,
                                 paired_asm, NO_NOTE))
                covered.add(canonical_gid(accession))

        self.logger.info('RefSeq: selected {:,} of {:,} genomes.'.format(len(selected), read))
        if misfiled:
            self.logger.warning('RefSeq: {:,} rows in the RefSeq assembly summary files are '
                                'not RefSeq assemblies and were ignored.'.format(misfiled))

        return selected, covered

    def select_genbank(self,
                       genbank_files: List[str],
                       covered: Set[str]) -> List[SelectedRow]:
        """Select the GenBank genomes needed to supplement the RefSeq genomes.

        A GenBank assembly is wanted when RefSeq does not already cover the
        genome and the assembly is not part of a large multi-isolate project.
        Coverage is decided by canonical accession against the RefSeq genomes
        actually selected above, not by the gbrs_paired_asm field, so a GenBank
        assembly paired with a RefSeq assembly that is absent from the new
        release is still selected. Where the two disagree NCBI describes a RefSeq
        assembly its own summary files do not list; that is common enough to be
        recorded in the notes column of the affected row rather than reported
        genome by genome.

        Parameters
        ----------
        genbank_files : list of str
            NCBI GenBank assembly summary files describing the new release.
        covered : set
            Canonical accessions of the selected RefSeq genomes.

        @return: list of selected rows (accession, ftp_path, version_status,
            excluded_from_refseq, gbrs_paired_asm, notes).
        """

        selected = []
        read = 0
        multi_isolate = 0
        unpaired = 0
        misfiled = 0

        for assembly_summary in genbank_files:
            for accession, version_status, paired_asm, excluded, ftp_path in self._read_summary(
                    assembly_summary):

                if not accession.startswith(GENBANK_PREFIX):
                    misfiled += 1
                    continue

                read += 1

                if not self._is_latest(accession, version_status):
                    continue

                if canonical_gid(accession) in covered:
                    continue

                if is_multi_isolate(excluded):
                    multi_isolate += 1
                    continue

                note = NO_NOTE
                if paired_asm.startswith(REFSEQ_PREFIX):
                    # NCBI names a RefSeq assembly that its own summary files do
                    # not list. The GenBank genome is selected regardless, as
                    # nothing else covers the genome, but the row says why it is
                    # here despite appearing to have a RefSeq counterpart.
                    unpaired += 1
                    note = ('paired RefSeq assembly {} absent from the assembly '
                            'summary files'.format(paired_asm))

                selected.append((accession, ftp_path, version_status, excluded,
                                 paired_asm, note))

        self.logger.info('GenBank: selected {:,} of {:,} genomes.'.format(len(selected), read))
        self.logger.info('GenBank: {:,} genomes skipped as part of a large multi-isolate '
                         'project.'.format(multi_isolate))
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

        refseq, covered = self.select_refseq(refseq_files)
        genbank = self.select_genbank(genbank_files, covered)

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

    Three things differ from doing it by hand. NCBI names all four assembly
    summary files assembly_summary.txt, so the database and domain are put back
    into the name as they are saved, and they are gzipped on the way in rather
    than kept as the 1.8 GB of text NCBI serves. The taxonomy dump is checked
    against the MD5 NCBI publishes beside it, since a truncated taxdump is not
    obviously broken until a release has been built on it, and the archive is
    then discarded: it has been verified and unpacked, and nothing reads it
    again. Fungal genomes are a separate procedure and are not downloaded here.

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
        select_genomes and the taxonomy step below both expect to find them.

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
        rather than over files named by hand. The output prefix is a path, so the
        files land in taxonomy/standardised_taxonomy/ without the working
        directory being changed.

        Parameters
        ----------
        taxdump_dir : str
            Directory holding the extracted nodes.dmp and names.dmp.
        summaries : dict
            (database, domain) to assembly summary file.
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
