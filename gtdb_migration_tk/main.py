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

from gtdb_migration_tk.biolib_lite.common import check_file_exists, make_sure_path_exists
from gtdb_migration_tk.busco_manager import BuscoManager
from gtdb_migration_tk.checkm_database_manager import CheckM2DatabaseManager, CheckMDatabaseManager
from gtdb_migration_tk.checkm_manager import CheckM, CheckM2, CheckMManager
from gtdb_migration_tk.curation_lists import CurationLists
from gtdb_migration_tk import config
from gtdb_migration_tk.database_manager import DatabaseManager
from gtdb_migration_tk.directory_manager import DirectoryManager
from gtdb_migration_tk.trans_table import GTranslate
from gtdb_migration_tk.lpsn import LPSN
from gtdb_migration_tk.marker_alignment_manager import AlignmentError, MarkerAlignmentManager
from gtdb_migration_tk.marker_manager import BadHmmDatabase, MarkerManager
from gtdb_migration_tk.metadata_database_manager import MetadataDatabaseManager, NCBITaxDatabaseManager, MetadataTableError, confirm_partial_load
from gtdb_migration_tk.metadata_manager import EmptyGenomeDirs, MetadataManager, MetadataTable
from gtdb_migration_tk.metadata_ncbi_manager import NCBIMeta, NCBIMetaDir
from gtdb_migration_tk.ncbi_genome_category import GenomeCategoryError, GenomeType
from gtdb_migration_tk.ncbi_utils import BadInput
from gtdb_migration_tk.ncbi_strain_summary import NCBIStrainParser
from gtdb_migration_tk.ncbi_genome_sync import NCBIGenomeSync
from gtdb_migration_tk.ncbi_metadata_sync import NCBIMetadataSync
from gtdb_migration_tk.ncbi_tax_manager import TaxonomyNCBI
from gtdb_migration_tk.prodigal_manager import ProdigalManager
from gtdb_migration_tk.propagate_taxonomy import Propagate, PropagationError
from gtdb_migration_tk.rna_manager_ltp import RnaManagerLTP
from gtdb_migration_tk.rna_manager_silva import RnaManagerSILVA
from gtdb_migration_tk.select_genomes import SelectGenomes
from gtdb_migration_tk.seqcode_manager import SeqCodeError, SeqCodeManager
from gtdb_migration_tk.strains import Strains
from gtdb_migration_tk.trnascan_manager import tRNAScan
from gtdb_migration_tk.update_genomes import UpdateGenomes
from gtdb_migration_tk.utils.common import database_keywords
from gtdb_migration_tk.utils.tools import Tools


class OptionsParser():
    def __init__(self):
        """Initialization"""

        self.logger = logging.getLogger('timestamp')

    def full_lpsn_wf(self, options):
        """Full workflow to parse LPSN."""
        make_sure_path_exists(options.output_dir)
        p = LPSN(False, options.output_dir)
        p.full_lpsn_wf()

    def add_lpsn_metadata(self, options):
        p = LPSN(False, None)
        p.add_lpsn_metadata(database_keywords(options), options.lpsn_metadata_file)

    def update_taxid_to_db(self, options):
        p = TaxonomyNCBI()
        p.update_taxid_to_db(database_keywords(options), options.input_file)

    def pull_html(self, options):
        """Pull all genus.html files."""
        make_sure_path_exists(options.output_dir)
        p = LPSN(options.skip_taxa_per_letter_dl, options.output_dir)
        for rk in ['phylum', 'class', 'order', 'family', 'genus', 'species']:
           p.download_rank_lpsn_html(rk)
        p.download_subspecies_lpsn_html()

    def parse_html(self, options):
        """Parse all html files."""
        make_sure_path_exists(options.output_dir)
        p = LPSN(False, options.output_dir)
        p.parse_html(options.input_dir, options.lpsn_gss_file)

    def generate_date_table(self, options):
        p = Strains()
        p.generate_date_table(options.lpsn_scraped_species_info,
                              options.lpsn_gss_file,
                              options.output_file)

    def generate_ncbi_strains_summary(self, options):
        for assembly_summary in options.new_list_genomes:
            check_file_exists(assembly_summary)
        make_sure_path_exists(options.output_dir)
        p = NCBIStrainParser(options.new_list_genomes, options.cpus)
        p.generate_ncbi_strains_summary(
            options.gtdb_genome_path_file, options.output_dir)

    def generate_type_table(self, options):
        for assembly_summary in options.new_list_genomes:
            check_file_exists(assembly_summary)
        check_file_exists(options.seqcode_table)
        p = Strains(options.output_dir, options.cpus)
        p.generate_type_strain_table(options.gtdb_genome_path_file,
                                     options.new_list_genomes,
                                     options.ncbi_names,
                                     options.ncbi_nodes,
                                     options.lpsn_gss_file,
                                     options.lpsn_dir,
                                     options.year_table,
                                     options.seqcode_table)

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

    def ncbi_metadata_sync(self, options):
        make_sure_path_exists(options.output_dir)
        p = NCBIMetadataSync(options.output_dir, options.group)
        p.run(options.release_number)

    def select_genomes(self, options):
        for assembly_summary in options.new_list_genomes:
            check_file_exists(assembly_summary)
        make_sure_path_exists(options.output_dir)
        p = SelectGenomes(options.output_dir,
                          min_genome_size=options.min_genome_size,
                          max_genome_size=options.max_genome_size)
        p.run(options.new_list_genomes)

    def parse_genome_directory(self, options):
        p = DirectoryManager()
        check_file_exists(options.gtdb_selected_genomes)
        p.generate_genome_dir_file(options.genome_dir,
                                   options.output_file,
                                   options.gtdb_selected_genomes,
                                   options.cpus)

    def update_genomes(self, options):
        check_file_exists(options.ftp_genome_dirs)

        if not options.fresh:
            # argparse cannot make one argument required by the absence of another,
            # so the previous release is asked for here; checked before the output
            # directory is made, a run that cannot start leaving nothing behind
            if not options.old_genome_dirs:
                self.logger.error(
                    '--old_genome_dirs_file is required unless --fresh is given.\n')
                sys.exit(1)
            check_file_exists(options.old_genome_dirs)

        make_sure_path_exists(options.output_dir)

        # RefSeq and GenBank in one pass, with one pair of reports for the release
        p = UpdateGenomes(options.output_dir, options.dry_run, options.cpus,
                          options.resume)

        if options.fresh:
            # the mirror is the release: nothing compared, nothing inherited
            p.run_fresh(options.ftp_dir, options.ftp_genome_dirs)
        else:
            p.run_comparison(options.ftp_dir, options.ftp_genome_dirs,
                             options.old_genome_dirs)

    def run_trans_table(self, options):
        p = GTranslate(options.cpus,
                       options.batch_size,
                       options.tmp_dir,
                       options.force,
                       options.keep_called_genes,
                       options.prefix,
                       options.custom_model_path,
                       options.reclaim,
                       options.lease * 60 * 60)
        check_file_exists(options.gtdb_genome_path_file)
        check_file_exists(options.taxonomy_file)
        make_sure_path_exists(options.output_dir)
        p.run(options.gtdb_genome_path_file, options.taxonomy_file, options.output_dir)

    def run_prodigal(self, options):
        p = ProdigalManager(options.tmp_dir,
                            options.cpus,
                            options.batch_size,
                            options.reclaim,
                            options.lease * 60 * 60)
        check_file_exists(options.gtdb_genome_path_file)
        check_file_exists(options.trans_table_file)
        if options.tt_override:
            check_file_exists(options.tt_override)
        make_sure_path_exists(options.output_dir)
        p.run(options.gtdb_genome_path_file,
              options.trans_table_file,
              options.output_dir,
              options.tt_override,
              options.all_genomes)

    def marker_dir_suffix(self, options):
        """The --dir_suffix given, or the one config.py declares for --db.

        Parameters
        ----------
        options : argparse.Namespace
            Options of hmmsearch or top_hit.

        @return: suffix naming the marker directory and files, e.g. 33.1_lite.
        """

        if options.dir_suffix:
            return options.dir_suffix

        dir_suffix = config.MARKER_DIR_SUFFIX[options.db]
        self.logger.info('Using the {} version declared in config.py: --dir_suffix {}'.format(
            options.db, dir_suffix))
        return dir_suffix

    def run_hmmsearch(self, options):
        p = MarkerManager(options.tmp_dir,
                          options.cpus,
                          options.batch_size,
                          options.reclaim,
                          options.lease * 60 * 60)
        check_file_exists(options.gtdb_genome_path_file)
        check_file_exists(options.report)
        make_sure_path_exists(options.output_dir)
        try:
            p.run_hmmsearch(options.gtdb_genome_path_file,
                            options.report, options.db, self.marker_dir_suffix(options),
                            options.hmm_db_path,
                            options.output_dir,
                            options.all_genomes)
        except BadHmmDatabase as exc:
            # a mistyped --hmm_db_path is the user's to fix and says so in one
            # line, rather than coming out as a traceback through the argparse
            # frames; nothing has been claimed or written by this point
            self.logger.error(str(exc))
            sys.exit(1)

    def align_marker_genes(self, options):
        check_file_exists(options.gtdb_genome_path_file)
        make_sure_path_exists(options.output_dir)
        p = MarkerAlignmentManager(database_keywords(options), options.cpus, options.batch_size,
                                   options.tmp_dir, options.reclaim, options.lease * 60 * 60)
        try:
            finished = p.run(options.marker_set_ids, options.all_genomes, options.gtdb_genome_path_file,
                             options.output_dir)
        except AlignmentError as exc:
            # the markers or the genome_dirs file to put right, said in one line;
            # nothing has been aligned
            self.logger.error(str(exc))
            sys.exit(1)
        if not finished:
            sys.exit(1)

    def run_tophit(self, options):
        p = MarkerManager('/tmp', options.cpus)
        p.run_tophit(options.gtdb_genome_path_file, options.db, self.marker_dir_suffix(options))

    def generate_metadata(self, options):
        p = MetadataManager(options.cpus)
        check_file_exists(options.gtdb_genome_path_file)
        make_sure_path_exists(options.output_dir)
        p.generate_metadata(options.gtdb_genome_path_file, options.output_dir)

    def create_metadata_tables(self, options):
        p = MetadataTable(options.silva_version)
        try:
            p.create_metadata_tables(
                options.gtdb_genome_path_file, options.output_dir, options.cpus)
        except EmptyGenomeDirs as exc:
            # the wrong file, or one an earlier step left empty: the user's to
            # fix, said in one line rather than as a traceback
            self.logger.error(str(exc))
            sys.exit(1)

    def parse_ncbi_assemblies(self, options):
        for assembly_summary in options.new_list_genomes:
            check_file_exists(assembly_summary)
        make_sure_path_exists(options.output_dir)
        p = NCBIMeta()
        p.parse_assemblies(options.new_list_genomes, options.output_dir)

    def generate_rna_silva(self, options):
        check_file_exists(options.gtdb_genome_path_file)
        check_file_exists(options.gtdb_domain_file)
        check_file_exists(options.taxonomy_file)
        make_sure_path_exists(options.output_dir)
        p = RnaManagerSILVA(options.rna_version,
                            options.rnapath,
                            options.rna_gene,
                            options.gtdb_domain_file,
                            options.taxonomy_file,
                            options.cpus,
                            options.tmp_dir,
                            options.batch_size,
                            options.reclaim,
                            options.lease * 60 * 60,
                            max_genome_size=options.max_genome_size)
        p.run(options.gtdb_genome_path_file,
              options.output_dir,
              options.all_genomes,
              options.remove)

    def update_silva(self, options):
        RnaManagerSILVA.update_silva(options.ssu_ref, options.lsu_ref, options.output_dir)

    def generate_rna_ltp(self, options):
        check_file_exists(options.gtdb_genome_path_file)
        make_sure_path_exists(options.output_dir)
        p = RnaManagerLTP(options.ltp_version,
                          options.ssu_version,
                          options.rnapath,
                          options.cpus,
                          options.tmp_dir,
                          options.batch_size,
                          options.reclaim,
                          options.lease * 60 * 60,
                          max_genome_size=options.max_genome_size)
        p.run(options.gtdb_genome_path_file,
              options.output_dir,
              options.all_genomes,
              options.remove)

    def generate_checkm_data(self, options, program=CheckM):
        p = program(options.cpus,
                    options.tmp_dir,
                    options.batch_size,
                    options.reclaim,
                    options.lease * 60 * 60,
                    max_genome_size=options.max_genome_size)
        p.run(options.gtdb_genome_path_file,
              options.report,
              options.output_dir,
              options.all_genomes)

    def generate_checkm2_data(self, options):
        self.generate_checkm_data(options, program=CheckM2)


    def join_checkm_files(self, options):
        p = CheckMManager()
        p.join_checkm_files_releases(options.checkm_files, options.output_file)

    def generate_busco_data(self, options):
        """Estimate quality of fungal genomes using BUSCO."""

        p = BuscoManager(options.cpus)
        p.run_busco(options.gtdb_genome_path_file, 
                    options.report,
                    options.output_dir, 
                    options.all_genomes)

    def update_db(self, options):
        p = DatabaseManager(database_keywords(options), options.cpus)
        p.run(options.gtdb_genome_path_file,
              options.report,
              options.ftp_download_date,
              options.output_dir,
              rehash_all=options.rehash_all,
              dry_run=options.dry_run)

    def update_checkm_db(self, options):
        make_sure_path_exists(options.output_dir)
        p = CheckMDatabaseManager(database_keywords(options))
        p.add_checkm_to_db(options.checkm_profile_file, options.check_qa, options.not_assessed,
                           options.output_dir)

    def update_checkm2_db(self, options):
        make_sure_path_exists(options.output_dir)
        p = CheckM2DatabaseManager(database_keywords(options))
        p.add_checkm2_to_db(options.checkm2_profile_file, options.not_assessed, options.output_dir)

    def update_metadata_db(self, options):
        if options.genome_list and not options.do_not_null_field:
            # asked before the database is reached, so no transaction waits on the answer
            confirm_partial_load(options.genome_list)
        p = MetadataDatabaseManager(database_keywords(options))
        try:
            p.process_metadata_files(options.genome_list, do_not_null_field=options.do_not_null_field,
                                     table_folder=options.input_folder, table_file=options.metadata_table,
                                     table_file_desc=options.metadata_table_desc)
        except MetadataTableError as exc:
            # a table or description to put right, said in one line; nothing was written
            self.logger.error(str(exc))
            sys.exit(1)
        self.logger.info('Update metadata Done.')

    def update_reps_db(self, options):
        p = MetadataDatabaseManager(database_keywords(options))
        p.update_reps(options.final_cluster_file)

    def update_ncbi_tax_db(self, options):
        for path in (options.organism_names, options.filtered, options.unfiltered):
            check_file_exists(path)
        if options.genome_list and not options.do_not_null_field:
            # asked before the database is reached, so no transaction waits on the answer
            confirm_partial_load(options.genome_list, 'update_ncbi_tax_db')
        make_sure_path_exists(options.output_dir)
        p = NCBITaxDatabaseManager(database_keywords(options))
        try:
            p.update_ncbi_tax_db(options.organism_names, options.filtered, options.unfiltered,
                                 options.genome_list, options.output_dir, options.do_not_null_field)
        except MetadataTableError as exc:
            # a file to put right, said in one line; nothing was written
            self.logger.error(str(exc))
            sys.exit(1)

    def add_surveillance_genomes(self, options):
        for assembly_summary in options.new_list_genomes:
            check_file_exists(assembly_summary)
        p = MetadataDatabaseManager(database_keywords(options))
        try:
            p.add_surveillance_genomes(options.new_list_genomes)
        except (MetadataTableError, BadInput) as exc:
            # the files to put right, said in one line; nothing was written
            self.logger.error(str(exc))
            sys.exit(1)

    def add_names_dmp(self, options):
        p = TaxonomyNCBI()
        p.populate_names_dmp_table(options.taxonomy_dir,
                                   options.ra, options.rb, options.ga, options.gb,
                                   options.output_file)

    def propagate_gtdb_taxonomy(self, options):
        for metadata_file in options.gtdb_metadata_prev:
            check_file_exists(metadata_file)
        make_sure_path_exists(options.output_dir)
        p = Propagate(database_keywords(options))
        try:
            p.propagate_taxonomy(options.gtdb_metadata_prev, options.output_dir)
        except PropagationError as exc:
            # the files or the database to put right, said in one line; nothing was written
            self.logger.error(str(exc))
            sys.exit(1)

    def propagate_curated_taxonomy(self, options):
        p = Propagate()
        p.propagate_taxonomy_from_reps_to_cluster(options.taxonomy_file, options.metadata, options.output_file)

    def add_taxonomy_to_database(self, options):
        p = Propagate(database_keywords(options))
        p.add_taxonomy_to_database(options.taxonomy_file, options.metadata, options.truncate_taxonomy)

    def update_propagated_tax(self, options):
        p = Propagate(database_keywords(options))
        p.add_propagated_taxonomy(options.input_dir)

    def set_gtdb_domain(self, options):
        make_sure_path_exists(options.output_dir)
        p = Propagate(database_keywords(options))
        p.set_gtdb_domain(options.output_dir)

    def parse_ncbi_genome_category(self, options):
        for assembly_summary in options.new_list_genomes:
            check_file_exists(assembly_summary)
        check_file_exists(options.gtdb_genome_path_file)
        make_sure_path_exists(options.output_dir)
        p = GenomeType(options.cpus)
        try:
            p.run(options.new_list_genomes, options.gtdb_genome_path_file, options.output_dir)
        except GenomeCategoryError as exc:
            # a value of the summaries for a person to place, said in one line;
            # nothing was written
            self.logger.error(str(exc))
            sys.exit(1)

    def generate_trnascan_data(self, options):
        check_file_exists(options.gtdb_genome_path_file)
        check_file_exists(options.gtdb_domain_file)
        check_file_exists(options.taxonomy_file)
        make_sure_path_exists(options.output_dir)
        p = tRNAScan(options.gtdb_domain_file,
                     options.taxonomy_file,
                     options.cpus,
                     options.tmp_dir,
                     options.batch_size,
                     options.reclaim,
                     options.lease * 60 * 60,
                     max_genome_size=options.max_genome_size)
        p.run(options.gtdb_genome_path_file,
              options.output_dir,
              options.all_genomes)

    def parse_ncbi_dir(self, options):
        check_file_exists(options.gtdb_genome_path_file)
        make_sure_path_exists(options.output_dir)
        p = NCBIMetaDir(options.cpus)
        p.parse_ncbi_dir(options.gtdb_genome_path_file, options.output_dir)

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
        check_file_exists(options.metadata)
        p = Tools()
        p.parse_ncbi_names_and_nodes(options.name, options.node, options.metadata, options.output_file)

    def compare_metadata_genome_dir(self, options):
        check_file_exists(options.metadata)
        check_file_exists(options.gtdb_genome_path_file)
        p = Tools()
        p.compare_metadata_genome_dir(options.metadata, options.gtdb_genome_path_file)

    def generate_ltp_db(self, options):
        check_file_exists(options.fasta)
        check_file_exists(options.csv)
        check_file_exists(options.compressed_fasta)
        p = Tools()
        p.generate_ltp_db(options.csv, options.compressed_fasta, options.fasta,options.output_dir, options.output_prefix)

    def download_seqcode_data(self, options):
        check_file_exists(options.gtdb_genome_path_file)
        for assembly_summary in options.new_list_genomes:
            check_file_exists(assembly_summary)
        make_sure_path_exists(options.output_dir)
        p = SeqCodeManager()
        try:
            p.run(options.gtdb_genome_path_file, options.output_dir, options.new_list_genomes,
                  options.species_cache)
        except (SeqCodeError, BadInput) as exc:
            # the Registry or NCBI did not answer, or a summary lacks a column,
            # said in one line; no table was written
            self.logger.error(str(exc))
            sys.exit(1)

    def check_db_population(self, options):
        p = Tools()
        p.check_db_population(options.metadata, options.id_last_genome, options.log)

    def ncbi_genome_sync(self, options):
        """Sync a local mirror of NCBI genomes against the selected genomes table.

        Returns the exit code rather than raising: callers distinguish 75 (locked, retry
        later), 130/143 (signalled) and 74 (I/O) from a plain failure.
        """
        p = NCBIGenomeSync(options)
        return p.run()

    def parse_options(self, options):
        """Parse user options and call the correct pipeline(s)"""
        if options.subparser_name == 'ncbi_genome_sync':
            return self.ncbi_genome_sync(options)
        elif options.subparser_name == 'list_genomes':
            self.parse_genome_directory(options)
        elif options.subparser_name == 'generate_ltp_db':
            self.generate_ltp_db(options)
        elif options.subparser_name == 'trans_table':
            self.run_trans_table(options)
        elif options.subparser_name == 'prodigal':
            self.run_prodigal(options)
        elif options.subparser_name == 'ncbi_metadata_sync':
            self.ncbi_metadata_sync(options)
        elif options.subparser_name == 'select_genomes':
            self.select_genomes(options)
        elif options.subparser_name == 'hmmsearch':
            self.run_hmmsearch(options)
        elif options.subparser_name == 'align_marker_genes':
            self.align_marker_genes(options)
        elif options.subparser_name == 'top_hit':
            self.run_tophit(options)
        elif options.subparser_name == 'genomic_metadata':
            self.generate_metadata(options)
        elif options.subparser_name == 'create_tables':
            self.create_metadata_tables(options)
        elif options.subparser_name == 'download_seqcode_data':
            self.download_seqcode_data(options)
        elif options.subparser_name == 'parse_ncbi_assemblies':
            self.parse_ncbi_assemblies(options)
        elif options.subparser_name == "parse_ncbi_dir":
            self.parse_ncbi_dir(options)
        elif options.subparser_name == 'update_taxid_to_db':
            self.update_taxid_to_db(options)
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
        elif options.subparser_name == 'checkm2':
            self.generate_checkm2_data(options)
        elif options.subparser_name == 'update_checkm_db':
            self.update_checkm_db(options)
        elif options.subparser_name == 'update_checkm2_db':
            self.update_checkm2_db(options)
        elif options.subparser_name == 'busco':
            self.generate_busco_data(options)
        elif options.subparser_name == 'add_surveillance_genomes':
            self.add_surveillance_genomes(options)
        elif options.subparser_name == 'add_names_dmp':
            self.add_names_dmp(options)
        elif options.subparser_name == 'check_db_population':
            self.check_db_population(options)
        elif options.subparser_name == 'update_db':
            self.update_db(options)
        elif options.subparser_name == 'propagate_gtdb_taxonomy':
            self.propagate_gtdb_taxonomy(options)
        elif options.subparser_name == 'propagate_curated_taxonomy':
            self.propagate_curated_taxonomy(options)
        elif options.subparser_name == 'add_taxonomy_to_database':
            self.add_taxonomy_to_database(options)
        elif options.subparser_name == 'update_propagated_tax':
            self.update_propagated_tax(options)
        elif options.subparser_name == 'set_gtdb_domain':
            self.set_gtdb_domain(options)
        elif options.subparser_name == 'parse_ncbi_genome_category':
            self.parse_ncbi_genome_category(options)
        elif options.subparser_name == 'update_metadata_db':
            self.update_metadata_db(options)
        elif options.subparser_name == 'update_reps_db':
            self.update_reps_db(options)
        elif options.subparser_name == 'update_ncbi_tax_db':
            self.update_ncbi_tax_db(options)
        elif options.subparser_name == 'update_genomes':
            self.update_genomes(options)
        elif options.subparser_name == 'lpsn':
            if options.lpsn_subparser_name == 'lpsn_wf':
                self.full_lpsn_wf(options)
            elif options.lpsn_subparser_name == 'add_metadata':
                self.add_lpsn_metadata(options)
            elif options.lpsn_subparser_name == 'parse_html':
                self.parse_html(options)
            elif options.lpsn_subparser_name == 'pull_html':
                self.pull_html(options)
            else:
                self.logger.error('Unknown command: ' +
                                  options.lpsn_subparser_name + '\n')
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
        elif options.subparser_name == 'compare_metadata_genome_dir':
            self.compare_metadata_genome_dir(options)
        else:
            self.logger.error('Unknown command: ' +
                              options.subparser_name + '\n')
            sys.exit()

        self.logger.info('Done.')

        return 0




