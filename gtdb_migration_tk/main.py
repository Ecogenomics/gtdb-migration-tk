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

import sys
import logging

from biolib.common import check_file_exists, make_sure_path_exists

from gtdb_migration_tk.lpsn import LPSN
from gtdb_migration_tk.bacdive import BacDive
from gtdb_migration_tk.propagate_taxonomy import Propagate
from gtdb_migration_tk.strains import Strains
from gtdb_migration_tk.ncbi_strain_summary import NCBIStrainParser
from gtdb_migration_tk.tools import Tools
from gtdb_migration_tk.genome_manager import DirectoryManager
from gtdb_migration_tk.ftp_manager import RefSeqManager, GenBankManager
from gtdb_migration_tk.prodigal_manager import ProdigalManager
from gtdb_migration_tk.marker_manager import MarkerManager
from gtdb_migration_tk.metadata_manager import MetadataManager, MetadataTable
from gtdb_migration_tk.metadata_ncbi_manager import NCBIMeta, NCBIMetaDir
from gtdb_migration_tk.rna_manager import RnaManager
from gtdb_migration_tk.rna_ltp_manager import RnaLTPManager
from gtdb_migration_tk.trnascan_manager import tRNAScan
from gtdb_migration_tk.checkm_manager import CheckMManager
from gtdb_migration_tk.ncbi_tax_manager import TaxonomyNCBI
from gtdb_migration_tk.database_manager import DatabaseManager
from gtdb_migration_tk.curation_lists import CurationLists
from gtdb_migration_tk.checkm_database_manager import CheckMDatabaseManager
from gtdb_migration_tk.metadata_database_manager import MetadataDatabaseManager, NCBITaxDatabaseManager
from gtdb_migration_tk.ncbi_genome_category import GenomeType


class OptionsParser():
    def __init__(self):
        """Initialization"""

        self.logger = logging.getLogger('timestamp')

    def full_lpsn_wf(self, options):
        """Full workflow to parse LPSN."""
        make_sure_path_exists(options.output_dir)
        p = LPSN(options.output_dir)
        p.full_lpsn_wf()

    def pull_html(self, options):
        """Pull all genus.html files."""
        make_sure_path_exists(options.output_dir)
        p = LPSN(options.output_dir)
        p.download_family_lpsn_html()
        p.download_genus_lpsn_html()
        p.download_species_lpsn_html()
        p.download_subspecies_lpsn_html()

    def parse_html(self, options):
        """Parse all html files."""
        make_sure_path_exists(options.output_dir)
        p = LPSN(options.output_dir)
        p.parse_html(options.input_dir)

    def download_strains(self, options):
        make_sure_path_exists(options.output_dir)
        p = BacDive(options.output_dir, options.username, options.pwd)
        p.download_strains()

    def generate_date_table(self, options):
        p = Strains()
        p.generate_date_table(options.lpsn_scraped_species_info,
                              options.lpsn_gss_file,
                              options.output_file)

    def generate_ncbi_strains_summary(self, options):
        p = NCBIStrainParser(options.assembly_summary_bacteria_genbank, options.assembly_summary_archaea_genbank,
                             options.assembly_summary_bacteria_refseq, options.assembly_summary_archaea_refseq)
        p.generate_ncbi_strains_summary(
            options.genome_dir_file, options.out_dir)

    def generate_type_table(self, options):
        p = Strains(options.output_dir, options.cpus)
        p.generate_type_strain_table(options.metadata_file,
                                     options.ncbi_names,
                                     options.ncbi_nodes,
                                     options.lpsn_dir,
                                     options.dsmz_dir,
                                     options.year_table,
                                     options.source_strain)

    def compare_metadata(self, options):
        p = Tools()
        p.compare_metadata(options.previous_metadata_file,
                           options.new_metadata_file,
                           options.only_ncbi,
                           options.use_formatted_id)

    def compare_markers(self, options):
        check_file_exists(options.first_domain_report)
        check_file_exists(options.second_domain_report)
        p = Tools()
        p.compare_markers(options.first_domain_report,
                           options.second_domain_report,
                          options.output_file,
                           options.only_ncbi,
                           options.use_formatted_id)

    def compare_selected_data(self, options):
        p = Tools()
        p.compare_selected_data(options.previous_metadata_file,
                                options.new_metadata_file,
                                options.field_of_interest,
                                options.output_file, options.only_ncbi)

    def clean_ftp(self, options):
        p = DirectoryManager()
        p.clean_ftp(options.new_list_genomes,
                    options.ftp_genome_dir_file, options.ftp_genome_dir,
                    options.report_dir,
                    options.taxonomy_file)

    def parse_genome_directory(self, options):
        p = DirectoryManager()
        p.run(options.genome_dir, options.output_file)

    def update_refseq_from_ftp_files(self, options):
        p = RefSeqManager(options.output_dir, options.dry_run, options.cpus)
        p.runComparison(
            options.ftp_refseq, options.output_dir, options.ftp_genome_dirs, options.old_genome_dirs,
            options.arc_assembly_summary, options.bac_assembly_summary)

    def update_genbank_from_ftp_files(self, options):
        p = GenBankManager(options.output_dir, options.dry_run, options.cpus)
        p.runComparison(
            options.ftp_genbank, options.output_dir, options.ftp_genbank_genome_dirs,
            options.old_genbank_genome_dirs, options.new_refseq_genome_dirs,
            options.arc_assembly_summary, options.bac_assembly_summary)

    def run_prodigal(self, options):
        p = ProdigalManager(options.tmp_dir, options.cpus)
        p.run(options.gtdb_genome_path_file, options.all_genomes)

    def run_prodigal_check(self, options):
        p = ProdigalManager()
        p.run_prodigal_check(options.gtdb_genome_path_file)

    def run_hmmsearch(self, options):
        p = MarkerManager(options.tmp_dir, options.cpus)
        p.run_hmmsearch(options.gtdb_genome_path_file,
                        options.report, options.db)
        self.logger.info('Done.')

    def run_tophit(self, options):
        p = MarkerManager('/tmp', options.cpus)
        p.run_tophit(options.gtdb_genome_path_file, options.db)

    def generate_metadata(self, options):
        p = MetadataManager(options.cpus)
        p.generate_metadata(options.gtdb_genome_path_file)

    def create_metadata_tables(self, options):
        p = MetadataTable(options.silva_version)
        p.create_metadata_tables(
            options.gtdb_genome_path_file, options.output_dir)

    def parse_assemblies(self, options):
        p = NCBIMeta()
        p.parse_assemblies(options.rb,
                           options.ra,
                           options.gb,
                           options.ga, options.genome_list, options.output_file)

    def generate_rna_silva(self, options):
        p = RnaManager(options.cpus, options.version,
                       options.rnapath, options.rna_gene)
        p.generate_rna_silva(options.gtdb_genome_path_file)

    def update_silva(self, options):
        p = RnaManager(1, None, None, None)
        p.update_silva(options.ssu_ref, options.lsu_ref, options.out_dir)

    def generate_rna_ltp(self, options):
        p = RnaLTPManager(options.cpus, options.ltp_version,
                          options.ssu_version, options.rnapath)
        p.generate_rna_ltp(options.gtdb_genome_path_file, options.all_genomes)

    def generate_checkm_data(self, options):
        p = CheckMManager(options.cpus)
        p.run_checkm(options.gtdb_genome_path_file, options.report,
                     options.output_dir, options.all_genomes)

    def join_checkm_files(self, options):
        p = CheckMManager()
        p.join_checkm_files_releases(options.checkm_files, options.output_file)

    def update_db(self, options):
        p = DatabaseManager(options.user, options.hostname,
                            options.db, options.password,
                            options.ftp_download_date,
                            options.repository,
                            options.output_dir, options.cpus)
        p.runUpdate(options.checkm_profile_new_genomes,
                    options.genome_dirs_file, options.ftp_download_date)

    def update_checkm_db(self, options):
        p = CheckMDatabaseManager(options.hostname, options.user, options.password, options.db)
        p.add_checkm_to_db(options.checkm_profiles, options.check_qa, options.metadata)

    def update_metadata_db(self, options):
        p = MetadataDatabaseManager(options.hostname, options.user, options.password, options.db)
        p.process_metadata_files(options.genome_list, do_not_null_field=options.do_not_null_field,
                                 table_folder=options.input_folder, table_file=options.metadata_table,
                                 table_file_desc=options.metadata_table_desc)

    def update_ncbitax_db(self, options):
        p = NCBITaxDatabaseManager(options.hostname, options.user, options.password, options.db)
        p.update_ncbitax_db(options.organism_name_file, options.filtered, options.unfiltered, options.genome_list,
                            options.do_not_null_field)

    def propagate_gtdb_taxonomy(self, options):
        p = Propagate()
        p.propagate_taxonomy(options.gtdb_metadata_prev, options.gtdb_metadata_cur, options.taxonomy_file,
                             options.rep_file)

    def update_propagated_tax(self, options):
        p = Propagate(options.hostname, options.user, options.password, options.db)
        p.add_propagated_taxonomy(options.taxonomy_file, options.metadata_file, options.genome_list,
                                  options.truncate_taxonomy, options.rep_id_file)

    def set_gtdb_domain(self, options):
        p = Propagate(options.hostname, options.user, options.password, options.db)
        p.set_gtdb_domain()

    def ncbi_genome_category(self, options):
        p = GenomeType()
        p.run(options.genbank_assembly_summary,
              options.refseq_assembly_summary,
              options.genome_file,
              options.output_file)

    def generate_trnascan_data(self, options):
        p = tRNAScan(options.gbk_arc_assembly_file, options.gbk_bac_assembly_file,
                     options.rfq_arc_assembly_file, options.rfq_bac_assembly_file, options.cpus)
        p.run(options.gtdb_genome_path_file)

    def parse_ncbi_dir(self, options):
        p = NCBIMetaDir()
        p.parse_ncbi_dir(options.gtdb_genome_path_file, options.output_file)

    def parse_ncbi_taxonomy(self, options):
        p = TaxonomyNCBI()
        p.parse_ncbi_taxonomy(options.taxonomy_dir,
                              options.ra, options.rb, options.ga, options.gb,
                              options.output_prefix)

    def curation_lists(self, options):
        check_file_exists(options.gtdb_init_taxonomy)
        check_file_exists(options.gtdb_sp_clusters)
        check_file_exists(options.gtdb_prev_sp_clusters)
        check_file_exists(options.gtdb_decorate_table)
        make_sure_path_exists(options.output_dir)

        p = CurationLists(options.domain, options.output_dir)
        p.run(options.gtdb_init_taxonomy,
              options.gtdb_sp_clusters,
              options.gtdb_prev_sp_clusters,
              options.gtdb_decorate_table)



    def check_unique_strains(self, options):
        check_file_exists(options.node)
        check_file_exists(options.name)
        check_file_exists(options.metadata_file)
        p = Tools()
        p.parse_ncbi_names_and_nodes(options.name, options.node, options.metadata_file, options.output_file)

    def parse_options(self, options):
        """Parse user options and call the correct pipeline(s)"""
        if options.subparser_name == 'list_genomes':
            self.parse_genome_directory(options)
        elif options.subparser_name == 'prodigal':
            self.run_prodigal(options)
        elif options.subparser_name == 'clean_ftp':
            self.clean_ftp(options)
        elif options.subparser_name == 'prodigal_check':
            self.run_prodigal_check(options)
        elif options.subparser_name == 'hmmsearch':
            self.run_hmmsearch(options)
        elif options.subparser_name == 'top_hit':
            self.run_tophit(options)
        elif options.subparser_name == 'metadata':
            self.generate_metadata(options)
        elif options.subparser_name == 'create_tables':
            self.create_metadata_tables(options)
        elif options.subparser_name == 'parse_assemblies':
            self.parse_assemblies(options)
        elif options.subparser_name == "parse_ncbi_dir":
            self.parse_ncbi_dir(options)
        elif options.subparser_name == 'parse_ncbi_taxonomy':
            self.parse_ncbi_taxonomy(options)
        elif options.subparser_name == 'rna_silva':
            self.generate_rna_silva(options)
        elif options.subparser_name == 'update_silva':
            self.update_silva(options)
        elif options.subparser_name == 'rna_ltp':
            self.generate_rna_ltp(options)
        elif options.subparser_name == 'trnascan':
            self.generate_trnascan_data(options)
        elif options.subparser_name == 'join_checkm':
            self.join_checkm_files(options)
        elif options.subparser_name == 'checkm':
            self.generate_checkm_data(options)
        elif options.subparser_name == 'update_checkm_db':
            self.update_checkm_db(options)
        elif options.subparser_name == 'update_db':
            self.update_db(options)
        elif options.subparser_name == 'propagate_gtdb_taxonomy':
            self.propagate_gtdb_taxonomy(options)
        elif options.subparser_name == 'update_propagated_tax':
            self.update_propagated_tax(options)
        elif options.subparser_name == 'set_gtdb_domain':
            self.set_gtdb_domain(options)
        elif options.subparser_name == 'ncbi_genome_category':
            self.ncbi_genome_category(options)
        elif options.subparser_name == 'update_metadata_db':
            self.update_metadata_db(options)
        elif options.subparser_name == 'update_ncbitax_db':
            self.update_ncbitax_db(options)
        elif options.subparser_name == 'update_refseq':
            self.update_refseq_from_ftp_files(options)
        elif options.subparser_name == 'update_genbank':
            self.update_genbank_from_ftp_files(options)
        elif options.subparser_name == 'lpsn':
            if options.lpsn_subparser_name == 'lpsn_wf':
                self.full_lpsn_wf(options)
            elif options.lpsn_subparser_name == 'parse_html':
                self.parse_html(options)
            elif options.lpsn_subparser_name == 'pull_html':
                self.pull_html(options)
            else:
                self.logger.error('Unknown command: ' +
                                  options.lpsn_subparser_name + '\n')
        elif options.subparser_name == 'bacdive':
            if options.bacdive_subparser_name == 'download_strains':
                self.download_strains(options)
        elif options.subparser_name == 'ncbi_strains':
            self.generate_ncbi_strains_summary(options)
        elif options.subparser_name == 'strains':
            if options.strains_subparser_name == 'date_table':
                self.generate_date_table(options)
            if options.strains_subparser_name == 'type_table':
                self.generate_type_table(options)
        elif options.subparser_name == 'overview':
            self.compare_metadata(options)
        elif options.subparser_name == 'compare_markers':
            self.compare_markers(options)
        elif options.subparser_name == 'compare_field':
            self.compare_selected_data(options)
        elif options.subparser_name == 'curation_lists':
            self.curation_lists(options)
        elif options.subparser_name == 'check_unique_strains':
            self.check_unique_strains(options)
        else:
            self.logger.error('Unknown command: ' +
                              options.subparser_name + '\n')
            sys.exit()

        return 0
