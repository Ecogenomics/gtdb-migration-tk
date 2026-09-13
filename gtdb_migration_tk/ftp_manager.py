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
ftp_manager.py -- select the NCBI genomes belonging in a GTDB release, and update
the GTDB genome directories to match.

Each GTDB release is built by comparing the genomes held by the previous release
with the genomes NCBI currently offers on its FTP site. Genomes NCBI no longer
offers are removed, genomes new to NCBI are copied across, and genomes held by
both are checked for files that have changed since the previous release. The
copying, deleting, and checksumming is done by FTPTools in ftp_manager_tools.py;
what lives here is the decision about which genomes are wanted in the first place.

RefSeq and GenBank are handled separately, by RefSeqManager and GenBankManager,
because that decision differs between them. Every RefSeq assembly flagged by NCBI
as the latest version is wanted. GenBank is only consulted where RefSeq falls
short, so a GenBank assembly is wanted when it has no RefSeq counterpart, or when
that counterpart is missing from the FTP site or holds no genome assembly. Each
GenBank decision is recorded in gca_selection.log, as the selection is the part of
a release that is hardest to reconstruct after the fact.

Which genomes NCBI considers latest, which are surveillance genomes, and which are
paired with a RefSeq assembly are all read from the NCBI assembly summary files by
ncbi_utils.py, which ncbi_sync.py reads them with too.
"""

import os
import logging
from contextlib import ExitStack
from typing import Dict, List, Set

from tqdm import tqdm

from gtdb_migration_tk.biolib_lite.common import canonical_gid
from gtdb_migration_tk.ftp_manager_tools import FTPTools
from gtdb_migration_tk.ncbi_utils import genome_assembly_file, read_assembly_summary
from gtdb_migration_tk.utils.common import count_lines


# Domain labels held in the genome domain map passed to FTPTools.
ARCHAEA = 'Archaea'
BACTERIA = 'Bacteria'


class GenericDatabaseManager:
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

        self.logger.info('Previous release: {} genomes.'.format(len(old_genomes)))

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

        self.logger.info('FTP site: {} genomes.'.format(len(new_genomes)))

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
        self.logger.info('{0} genomes to remove'.format(len(removed_dict)))

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
        self.logger.info('{0} genomes to add'.format(len(added_dict)))

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
        self.logger.info('{} genomes to compare'.format(len(intersect_list)))

        return intersect_list


class RefSeqManager(GenericDatabaseManager):
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
            self.logger.info('NCBI: {} latest assemblies.'.format(len(self.genome_domain_dict)))

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


class GenBankManager(GenericDatabaseManager):
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
        self.logger.info('Indexed {} RefSeq genome directories.'.format(len(refseq_dirs)))

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

        self.logger.info('Selected {} GenBank genomes.'.format(len(selected_gca)))

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
