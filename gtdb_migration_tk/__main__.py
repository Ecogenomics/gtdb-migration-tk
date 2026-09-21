#!/usr/bin/env python3
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #r
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = "Pierre Chaumeil"
__copyright__ = "Copyright 2019"
__credits__ = ["Donovan Parks", "Pierre Chaumeil", "Aaron Mussig"]
__license__ = "GPL3"
__maintainer__ = "Pierre Chaumeil"
__email__ = "uqpchaum@uq.edu.au"
__status__ = "Development"

import os
import sys
import argparse
from contextlib import contextmanager
from datetime import datetime

from gtdb_migration_tk import __version__
from gtdb_migration_tk.biolib_lite.custom_help_formatter import CustomHelpFormatter
from gtdb_migration_tk.biolib_lite.logger import logger_setup
from gtdb_migration_tk.main import OptionsParser
from gtdb_migration_tk.ncbi_metadata_sync import NCBI_GROUPS
from gtdb_migration_tk.ncbi_genome_sync import add_sync_arguments


def print_help():
    """Help menu."""

    print('')
    print('                ...::: GTDB Migration Toolkit v' +
          __version__ + ' :::...''')
    print('''\

    NCBI data sync:
      ncbi_metadata_sync -> Sync NCBI metadata to local directory
      select_genomes     -> Select NCBI genomes which will comprise the new GTDB release
      ncbi_genome_sync   -> Sync NCBI data to local directory

    NCBI folder to GTDB folder:
      list_genomes   -> Produce file indicating the directory of each genome
      update_genomes -> Update RefSeq and GenBank genomes from the NCBI FTP mirror

    Call genes:
      call_genes_wf  -> Full call genes workflow (prodigal -> hmmsearch -> top_hit)
      trans_table    -> Predict translation table for each genome using gTranslate
      prodigal       -> Call genes using Prodigal
      hmmsearch      -> Search Tigrfam/Pfam markers genes and generate tophit files
      top_hit        -> Generate tophit files
      metadata       -> Generate metadata derived from nucleotide (e.g., GC) and protein (e.g., gene count) files.
      rna_silva      -> Identify, extract, and taxonomically classify 16S, 23S, and 5S rRNA genes in genomes against SILVA
      rna_ltp        -> Identify, extract, and taxonomically classify 16S rRNA genes against the LTP DB
      trnascan       -> Identifies tRNAs in genomes
      join_checkm    -> Join CheckM output files for different releases
      checkm         -> Estimates the quality of the new genomes
      busco          -> Estimate quality of new fungal genomes
      
    Access to Database:
     update_db          -> Update the GTDB database
     update_checkm_db   -> Import CheckM estimates
     update_metadata_db -> Update metadata in database
     update_reps_db     -> Update species cluter representatives in database

    Metadata:
      create_tables     -> Create tables with metadata for all genomes (currently only NCBI)
      parse_assemblies  -> Create tables with metadata for all NCBI genomes from assembly summaries
      parse_ncbi_dir    -> Create tables with metadata for all NCBI genomes from directories

    GTDB Taxonomy:
      propagate_gtdb_taxonomy -> Propagating GTDB taxonomy to new release
      update_propagate_tax    -> Push propagated taxonomy to new DB

    Information from online resources:
      lpsn         -> Process steps for LPSN
      bacdive      -> Process steps for BacDive [In Dev]
      strains      -> Set of tools to combined information from LPSN and DSMZ
      ncbi_strains -> Parse the assembly report file, the genomic.gbff file and the wgsmaster.gbff to find all strain ids

    Curation files
      curation_lists -> Lists and pseudo-trees for new representatives, polyphyletic taxa, rogue genomes, and genomes with modified NCBI names

    Miscellaneous commands:
      generate_ltp_db -> Generate LTP database

    Test suite for data validation:
      overview             -> Compare the Metadata file from the previous version with the new one
      compare_field        -> Compare a specific metadata field between to metadata files
      check_unique_strains -> Check if a genomes has to different strains from a same collection
      check_db_population  -> Check if the database is populated with the correct number of genomes
      

  Use: gtdb_migration_tk <command> -h for command specific help.

  Feature requests or bug reports can be sent to Donovan Parks (donovan.parks@gmail.com)
    or posted on GitHub (https://github.com/Ecogenomics/gtdb_migration_tk).
    
    ''')

    #get python filename
    print(f'run from {os.getcwd()}, {os.path.basename(__file__)} ')


@contextmanager
def subparser(parser, name, desc):
    yield parser.add_parser(name, conflict_handler='resolve', help=desc,
                            formatter_class=CustomHelpFormatter)


@contextmanager
def mutex_group(parser, required):
    group = parser.add_argument_group(f'mutually exclusive {"required" if required else "optional"} arguments')
    yield group.add_mutually_exclusive_group(required=required)


@contextmanager
def arg_group(parser, name):
    yield parser.add_argument_group(name)


def valid_date(s):
    try:
        return datetime.strptime(s, "%Y-%m-%d")
    except ValueError:
        msg = "Not a valid date: '{0}'.".format(s)
        raise argparse.ArgumentTypeError(msg)


def __all_genomes(group):
    group.add_argument('--all', dest='all_genomes', help="Re-run all genomes", action='store_true')




def __checkm_files(group, required):
    group.add_argument('-f',
                       '--checkm_files', help='Output CheckM files from different releases.', nargs="+",
                       required=required)


def __checkm_profile(group, required):
    group.add_argument('-c', '--checkm_profile', help='CheckM profile for new genomes', required=required)


def __checkm_qa(group, required):
    group.add_argument('-q', '--check_qa',
                       help='CheckM QA file for 100%% strain heterogeneity for all genomes of interest.',
                       required=required)


def __cpus(group, default=1):
    group.add_argument('-c', '--cpus', type=int ,default=default, help='Number of threads.')


def __database_setup(group, required):
    group.add_argument('--hostname', help='Hostname', required=required)
    group.add_argument('-u', '--user', help='PostgreSQL username', required=required)
    group.add_argument('-d', '--db', help='Database name.', required=required)
    group.add_argument('-p', '--password', help='Password for psql user', required=required)


def __do_not_null_field(group):
    group.add_argument('--do_not_null_field', help='Do not set fields to NULL for all genomes before updating values',
                       action='store_true')


def __domain(group, required):
    group.add_argument('--domain', required=required, help='Domain to append to output files', choices=['bac', 'ar'])


def __dry_run(group):
    group.add_argument('--dry_run', action='store_true',
                       help='Report what the run would do without copying any file. '
                            'Every genome is still compared -- the genomic FASTA MD5 of '
                            'the mirror against the previous release, and where those '
                            'differ the sequences themselves -- so the reports and the '
                            'counts of changed and unchanged genomes are the ones a real '
                            'run would produce. With --fresh there is nothing to compare, '
                            'and the reports describe the release that would be built.')


def __dsmz_directory(group, required):
    group.add_argument('--dsmz_dir',
                       help='Directory including the 3 DSMZ result files (dsmz_strains.tsv, dsmz_species.tsv, '
                            'and dsmz_genera.tsv ).',
                       required=required)


def __field_of_interest(group, required):
    group.add_argument('--field_of_interest',
                       help='Common field to compare between files.',
                       required=required)


def __filtered_taxonomy(group, required):
    group.add_argument('--filtered', help='Filtered taxonomy file.', required=required)


def __first_domain_report(group, required):
    group.add_argument('--first_domain_report', required=required,
                       help='File generated from GTDB power domain report from early release.')


def __ftp_download_date(group, required):
    group.add_argument('-f', '--ftp_download_date', type=valid_date,
                       help='Date when data was downloaded.format YYYY-MM-DD.', required=required)







def __gbk_arc_assembly_file(group, required):
    group.add_argument('--ga','--gbk_arc_assembly_file', required=required, help="Archaeal Assembly summary file from GenBank")


def __gbk_bac_assembly_file(group, required):
    group.add_argument('--gb','--gbk_bac_assembly_file', required=required,
                       help="Bacterial Assembly summary file from GenBank")


def __genbank_assembly_summary(group, required):
    group.add_argument('-g', '--genbank_assembly_summary', required=required,
                       help='File from NCBI indicating metadata for genome assemblies in GenBank '
                            '(ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt)')


def __genome_directory(group, required):
    group.add_argument('-d', '--genome_dir', required=required,
                       help='Base directory leading to NCBI archaeal and bacterial genome assemblies')


def __gtdb_selected_genomes(group, required):
    group.add_argument('-g', '--gtdb_selected_genomes', required=required,
                       help='Table of the genomes selected for the new release, written by '
                            'select_genomes (gtdb_selected_genomes.tsv.gz).')


def __genome_list(group, required):
    group.add_argument('--genome_list',
                       help='Only process genomes in this list ( we can use the metadata file exported from GTDB',
                       required=required, default=None)


def __gtdb_comparison_report(group, required):
    group.add_argument('--report', help='Report log indicating new, modified, unmodified, ..., genomes',
                       required=required)


def __gtdb_decorate_table(group, required):
    group.add_argument('--gtdb_decorate_table', required=required,
                       help='Decoration table produced by PhyloRank decorate')


def __gtdb_genome_path_file(group, required):
    group.add_argument('-g', '--gtdb_genome_path_file', help='Genome paths to GTDB genomes.', required=required)


def __gtdb_domain_file(group, required):
    group.add_argument('-d', '--gtdb_domain_file',
                        help='File indicating predicted domain for each GTDB genomes', required=required)


def __gtdb_init_taxonomy(group, required):
    group.add_argument('--gtdb_init_taxonomy', required=required, help='Initial taxonomy for latest release')


def __gtdb_metadata_current_release(group, required):
    group.add_argument('--gtdb_metadata_cur', required=required,
                       help='GTDB metadata for current NCBI release.')


def __gtdb_metadata_previous_release(group, required):
    group.add_argument('--gtdb_metadata_prev', required=required,
                       help='GTDB metadata for previous NCBI release.')


def __gtdb_prev_sp_clusters(group, required):
    group.add_argument('--gtdb_prev_sp_clusters', required=required, help='Species clusters for previous release')


def __gtdb_sp_clusters(group, required):
    group.add_argument('--gtdb_sp_clusters', required=required, help='Species clusters for latest release')


def __input_dir(group, required):
    group.add_argument('--input_dir', '--in_dir', help='Input directory.', required=required)

def __input_file(group, required):
    group.add_argument('--input_file', '--in_file', help='Input file.', required=required)


def __log_file(group, required):
    group.add_argument('-l', '--log', required=required, help='Log file.')

def __csv_file(grp, required):
    grp.add_argument('--csv', required=required, help='CSV file.')

def __fasta_file(grp, required):
    grp.add_argument('--fasta', required=required, help='FASTA file.')

def __compressed_fasta(grp, required):
    grp.add_argument('--compressed_fasta', required=required, help='Compressed fasta file.')


def __lpsn_directory(group, required):
    group.add_argument('--lpsn_dir',
                       help='Directory including the 3 LPSN result files (lpsn_genera.tsv, lpsn_species.tsv, '
                            'and lpsn_strains.tsv ).',
                       required=required)


def __lpsn_gss_file(group, required):
    group.add_argument('--lpsn_gss_file',
                       help="Table from lpsn.dsmz.de with nomenclature information (lpsn_gss_<date>.csv)",
                       required=required)

def __lpsn_metadata_file(group, required):
    group.add_argument('--lpsn_metadata_file',
                       help="TSV file generated from LPSN parse command.",
                       required=required)


def __lpsn_scraped_species_info(group, required):
    group.add_argument('--lpsn_scraped_species_info',
                       help="LPSN species file created by LPSN website parsing (lpsn_species.tsv "
                            "from 'gtdb_migration_tk lpsn parse_html').",
                       required=required)


def __lsu_ref(group, required):
    group.add_argument('--lsu_ref', help='SILVA_xxx_LSURef_tax_silva.fasta file downloaded from Silva.',
                       required=required)


def __ltp_version(group, required):
    group.add_argument('-v', '--ltp_version', help='LTP version to use.', required=required)


def __marker_db(group, required):
    group.add_argument('-d', '--db', required=required, help='Pfam or tigrfam.', choices=['pfam', 'tigrfam'])


def __metadata_file(group, required):
    group.add_argument('-m', '--metadata', help='Metadata file generated from "gtdb metadata export".',
                       required=required)

def __report_dir(group, required):
    group.add_argument('--report_dir', required=required,
                       help='Output directory to list reports.')


def __release_number(group, required):
    group.add_argument('-r', '--release_number', type=int, required=required,
                       help='GTDB release number, e.g. 237.')


def __ncbi_group(group, required):
    # 'group' here is the argparse group; the argument itself is the group of
    # organisms, PROK or FUNGI, which are downloaded and curated separately
    group.add_argument('--group', required=required, choices=NCBI_GROUPS,
                       help='Group of organisms to download (PROK or FUNGI).')


def __metadata_input_folder(group):
    group.add_argument('-i', '--input_folder', help='Directory folder with all the tables created by metadata.')


def __metadata_table(group):
    group.add_argument('--metadata_table', help='Specific metadata table.')


def __metadata_table_description(group):
    group.add_argument('--metadata_table_desc', help='Specific metadata table description file.')


def __name(group, required):
    group.add_argument('--name', required=required, help='Species clusters for latest release')


def __ncbi_names(group, required):
    group.add_argument('--ncbi_names',
                       help='NCBI names.dmp file',
                       required=required)


def __ncbi_nodes(group, required):
    group.add_argument('--ncbi_nodes',
                       help='NCBI nodes.dmp file',
                       required=required)



def __new_list_genomes(group, required):
    group.add_argument('-n','--new_list_genomes', required=required, nargs='+',
                       help='NCBI assembly summary files indicating the genomes present in the new release.')


def __new_metadata_file(group, required):
    group.add_argument('--new_metadata_file',
                       help='File indicating metadata of each genome in latest GTDB version.',
                       required=required)




def __node(group, required):
    group.add_argument('--node', required=required, help='Initial taxonomy for latest release')



def __old_genome_dirs_file(group, required):
    group.add_argument('--old_genome_dirs_file', dest="old_genome_dirs", required=required,
                       help='Genome directory file of the previous GTDB release (generated by list_genomes). '
                            'Required unless --fresh is given, a fresh release having no previous release to read.')


def __ftp_directory(group, required):
    group.add_argument('--ftp_directory', dest="ftp_dir", required=required,
                       help='Root directory of the NCBI FTP mirror (the --root given to ncbi_genome_sync).')


def __new_directory(group, required):
    group.add_argument('--new_directory', dest="output_dir", required=required,
                       help='Root directory of the new GTDB release. RefSeq genomes are written under refseq/ and GenBank genomes under genbank/, each keeping NCBI\'s nesting (genbank/GCA/047/639/395/GCA_047639395.1_ASM4763939v1), alongside the reports and the genome_dirs.tsv of the new release.')


def __ftp_genome_dirs_file(group, required):
    group.add_argument('--ftp_genome_dirs_file', dest="ftp_genome_dirs", required=required,
                       help='Genome directory file of the NCBI FTP mirror (generated by list_genomes).')


def __fresh(group):
    group.add_argument('--fresh', action='store_true',
                       help='Start the release from the NCBI genome data alone. Every genome '
                            'of the mirror is copied into the new release and reported as new; '
                            'nothing is compared to the previous GTDB release and no derived '
                            'data (Prodigal results, marker hits, rRNA and tRNA calls) is '
                            'carried across, so everything derived from the genomes is to be '
                            'regenerated. --old_genome_dirs_file is not read.')


def __resume(group):
    group.add_argument('--resume', action='store_true',
                       help='Continue a run that was interrupted, reading the '
                            'genome_dirs.tsv it left in --new_directory and handling '
                            'again only the genomes it does not name. A genome is '
                            'written to that file once its directory is there, so the '
                            'genomes the run failed on or never reached are the ones '
                            'missing from it. The reports of the interrupted run are '
                            'carried forward for the genomes it finished, so the '
                            'release ends up described once through. Combine with '
                            '--dry_run to report what is left without building it.')


def __only_ncbi(group):
    group.add_argument('--only_ncbi', help='Only process NCBI genomes.',
                       action='store_true')


def __organism_names(group, required):
    group.add_argument('-o', '--organism_names', help='NCBI Organism name file.', default=None, required=required)


def __output_dir(group, required):
    group.add_argument('-o', '--output_dir', '--out_dir', help='Output directory.', required=required)


def __output_file(group, required):
    group.add_argument('-o', '--output_file', help='Output file.', required=required)


def __output_prefix(group, required):
    group.add_argument('-p', '--output_prefix', required=required, help='Output prefix')


def __password(group, required):
    group.add_argument('-p', '--pwd', '--password', help='Password', required=required)


def __previous_metadata_file(group, required):
    group.add_argument('--previous_metadata_file',
                       help='File indicating metadata of each genome in previous GTDB version.',
                       required=required)


def __refseq_assembly_summary(group, required):
    group.add_argument('-g', '--refseq_assembly_summary', required=required,
                       help='File from NCBI indicating metadata for genome assemblies in RefSeq '
                            '(ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt)')


def __report_file(group, required):
    group.add_argument('-r', '--report', help='Report file generated by update steps.', required=required)


def __report_folder(group, required):
    group.add_argument('-r', '--report_folder', help='Path to report directory.', required=required)


def __repository(group, required):
    group.add_argument('-r', '--repository', help='NCBI repository.', choices=['refseq', 'genbank'], required=required)


def __representative_file(group, required):
    group.add_argument('--rep_file', required=required,
                       help='GTDB representatives file.')


def __rfq_arc_assembly_file(group, required):
    group.add_argument('--ra','--rfq_arc_assembly_file', required=required, help="Archaeal Assembly summary file from RefSeq")


def __rfq_bac_assembly_file(group, required):
    group.add_argument('--rb','--rfq_bac_assembly_file', required=required, help="Bacterial Assembly summary file from RefSeq")


def __rna_file_path(group):
    group.add_argument('-p', '--rnapath', help='Path to rna Silva file',
                       default='/srv/db/silva/')


def __rna_gene(group, required):
    group.add_argument('-r', '--rna_gene', required=required,
                       choices=['ssu', 'lsu_23S', 'lsu_5S'], help="RRNA gene to process")


def __rna_version(group, required):
    group.add_argument('-v', '--rna_version', help='RNA Silva version.', required=required)


def __second_domain_report(group, required):
    group.add_argument('--second_domain_report', required=required,
                       help='File generated from GTDB power domain report from latest release.')


def __trans_table_file(group, required):
    group.add_argument('--trans_table', dest='trans_table_file', required=required,
                       help='Translation table of each genome '
                            '(gtranslate.translation_table_summary.tsv.gz, written by '
                            'trans_table; gzipped or not is read either way).')


def __tt_override(group):
    group.add_argument('--tt_override',
                       help='Corrections to --trans_table: a TSV with genome_id and '
                            'translation_table columns, overriding the prediction.')


def __batch_size(group, default=10000):
    group.add_argument('--batch_size', type=int, default=default,
                       help='Number of genomes per batch, each batch being the unit of '
                            'restart and of sharing between machines.')


def __reclaim(group):
    group.add_argument('--reclaim', action='store_true',
                       help='Take over a batch another machine holds before its claim has '
                            'expired. A claim untouched for --lease hours, and one this host '
                            'left behind in a process that has gone, are taken without it.')


def __lease(group, default=2.0):
    group.add_argument('--lease', type=float, default=default,
                       help='Hours a claim on a batch survives without the machine holding '
                            'it saying so, after which another machine may take the batch.')


def __gtranslate_no_force(group):
    group.add_argument('--no_force', dest='force', action='store_false',
                       help='Stop a batch when a single genome fails. gTranslate is given '
                            '--force otherwise, a genome it cannot process being named in '
                            'no_prediction.tsv rather than costing the whole batch.')


def __keep_called_genes(group):
    group.add_argument('--keep_called_genes', action='store_true',
                       help='Keep the genes gTranslate called under the predicted translation table.')


def __gtranslate_prefix(group):
    group.add_argument('--prefix',
                       help='Prefix of the gTranslate output files (default: gtranslate).')


def __custom_model_path(group):
    group.add_argument('--custom_model_path',
                       help='Classifiers to predict with (default: the models named by GTRANSLATE_MODEL_PATH).')


def __silent(group):
    group.add_argument('--silent', help="Suppress output", action='store_true')

def __truncate_taxonomy(group):
        group.add_argument('--truncate_taxonomy', help="Truncate taxonomy in database to only gtdb_domain", action='store_true')

def __silva_version(group, required):
    group.add_argument('-v', '--silva_version', help='Silva version to use.', required=required)


def __skip_taxa_per_letter_dl(group):
    group.add_argument('--skip_taxa_per_letter_dl', action='store_true',
                       help='Skip downloading the set of taxa under each rank; '
                            'assumes this information has previously been downloaded.')


def __ssu_ref(group, required):
    group.add_argument('--ssu_ref', required=required,
                       help='SILVA_xxx_SSURef_NR99_tax_silva.fasta file downloaded from Silva.')


def __ssu_version(group, required):
    group.add_argument('-v', '--ssu_version', help='SSu version to use.', required=required)


def __surveillance_list(group, required):
    group.add_argument('--genome_list', required=required, help='Surveillance genomes.')


def __taxonomy_directory(group, required):
    group.add_argument('-t', '--taxonomy_dir', required=required,
                       help='Directory containing NCBI taxonomy files (dmp files)')


def __taxonomy_file(group, required):
    group.add_argument('-t','--taxonomy_file', required=required,
                       help='Standardised taxonomy file from NCBI.')


def __tmp_dir(group):
    group.add_argument('--tmp_dir', help='Directory for scratch files; no results are written here', default='/tmp')


def __truncate_taxonomy(group):
    group.add_argument('--truncate_taxonomy',
                       help='Truncate current taxonomy strings to just the domain before updating taxonomy',
                       action='store_true')


def __unfiltered_taxonomy(group, required):
    group.add_argument('--unfiltered', help='Unfiltered taxonomy file.', required=required)


def __use_formatted_id(group):
    group.add_argument('--use_formatted_id',
                       help='Use formatted id to compared GCA and GCF ids',
                       action='store_true')


def __username(group, required):
    group.add_argument('-u', '--username', help='Username', required=required)


def __year_table_file(group, required):
    group.add_argument('--year_table',
                       help='Date table generated by generate_date_table.py.',
                       required=required)

def __folder_suffix(group, required=False):
    group.add_argument('--folder_suffix', default=None, required=required,
                       help='Suffix of the marker directory and files written inside prodigal/, e.g. '
                            '33.1_lite gives pfam_33.1_lite/. Defaults to the version declared for '
                            '--db in config.py, with _lite appended.')

def __canonical_gid_table(group, required):
    group.add('--canonical_gids', help='tsv table listing both the canonical ids and genome_ids',required=required)

def __gtdb_sp_clusters_file(group, required):
    group.add('--gtdb_sp_clusters_file', help='gtdb_cluster_de_novo.tsv file generate when we recreated clusters',required=required)

def __hmm_db_path(group, required):
    group.add_argument('--hmm_db_path',
                       help="Path to hmm db folder. for pfam_33.1 : '/srv/db/pfam/33.1/' , "
                            "for pfam_33.1_lite : '/srv/db/gtdb/marker_genes/hmms_extended_pfam33.1_tigr15/',"
                            "for tigrfam_15.0 : '/srv/db/tigrfam/15.0/TIGRFAMs_15.0_HMM/tigrfam.hmm',"
                            "for tigrfam_15.0_lite : '/srv/db/gtdb/marker_genes/hmms_extended_pfam33.1_tigr15/tigrfam.hmm')",
                       required=required)

def __rerun(group):
    group.add_argument('--rerun', help='Rerun all genomes', action='store_true')
    pass


def __remove(group, db_name):
    group.add_argument('--remove',
                       help=f'Remove all previous results in {db_name} directory (use with caution!)',
                       action='store_true')


def __checkm_summary_refseq(grp, required):
    grp.add_argument('--checkm_summary_refseq', required=required, help='CheckM summary file for RefSeq genomes.')

def __checkm_summary_genbank(grp, required):
    grp.add_argument('--checkm_summary_genbank', required=required, help='CheckM summary file for GenBank genomes.')

def __checkm2_output_dir(grp, required):
    grp.add_argument('--checkm2_output_dir', required=required, help='CheckM v2 output directory with batches.')

def __report(grp, required):
    grp.add_argument('--report', required=required,
                     help='Report log indicating new, modified, unmodified, ..., genomes')



def __id_last_genome(grp, required):
    grp.add_argument('--id_last_genome', required=required, help='ID of the last genome in the previous release.',
                        type=int,default= 0)
    
def __final_cluster_file(grp, required):
    grp.add_argument('--final_cluster_file', required=required,
                     help="Clusters for named species")


def get_main_parser():
    # Setup the main, and sub parsers.
    main_parser = argparse.ArgumentParser(prog='gtdb_migration_tk', add_help=False, conflict_handler='resolve')
    sub_parsers = main_parser.add_subparsers(help="--", dest='subparser_name')

    # Download the NCBI metadata a release is built from
    with subparser(sub_parsers, 'ncbi_metadata_sync',
                   'Sync NCBI metadata to local directory.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __output_dir(grp, required=True)
            __release_number(grp, required=True)
            __ncbi_group(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            __silent(grp)

    # Select the genomes comprising the new release
    with subparser(sub_parsers, 'select_genomes',
                   'Select NCBI genomes which will comprise the new GTDB release.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __new_list_genomes(grp, required=True)
            __output_dir(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            __silent(grp)

    # Call genes with prodigal
    with subparser(sub_parsers, 'trans_table',
                   'Predict the translation table of each genome using gTranslate.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __gtdb_genome_path_file(grp, required=True)
            __taxonomy_file(grp, required=True)
            __output_dir(grp, required=True)
            __log_file(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            __silent(grp)
            __tmp_dir(grp)
            __cpus(grp)
            __batch_size(grp)
            __reclaim(grp)
            __lease(grp)
            __gtranslate_no_force(grp)
            __keep_called_genes(grp)
            __gtranslate_prefix(grp)
            __custom_model_path(grp)

    with subparser(sub_parsers, 'prodigal', 'Call genes using Prodigal.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __gtdb_genome_path_file(grp, required=True)
            __trans_table_file(grp, required=True)
            __output_dir(grp, required=True)
            __log_file(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            __silent(grp)
            __tt_override(grp)
            __tmp_dir(grp)
            __cpus(grp)
            __batch_size(grp)
            __reclaim(grp)
            __lease(grp)
            __all_genomes(grp)

    # Search Pfam Tigrfam markers
    with subparser(sub_parsers, 'top_hit', 'Generate TopHit file for Tigrfam or Pfam.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __gtdb_genome_path_file(grp, required=True)
            __log_file(grp, required=True)
            __marker_db(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            __folder_suffix(grp)
            __silent(grp)
            __cpus(grp)

    # Generate metadata for each genome
    with subparser(sub_parsers, 'metadata',
                   'Generate metadata derived from nucleotide (e.g., GC) '
                   'and protein (e.g., gene count) files.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __gtdb_genome_path_file(grp, required=True)
            __log_file(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            __silent(grp)
            __cpus(grp)

    with subparser(sub_parsers, 'update_silva',
                   'Update the Taxonomy files and Blast database based on the latest SILVA release') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __ssu_ref(grp, required=True)
            __lsu_ref(grp, required=True)
            __output_dir(grp, required=True)
            __log_file(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            __silent(grp)

    with subparser(sub_parsers,
                   'rna_silva',
                   'Identify, extracts and taxonomically classifies 16S '
                   '23S, and 5S rRNA genes in genomes against SILVA') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __gtdb_genome_path_file(grp, required=True)
            __gtdb_domain_file(grp, required=True)
            __rna_version(grp, required=True)
            __rna_gene(grp, required=True)
            __log_file(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            __rna_file_path(grp)
            __silent(grp)
            __rerun(grp)
            __remove(grp, 'SILVA')
            __cpus(grp)
            __all_genomes(grp)

    with subparser(sub_parsers, 'rna_ltp',
                   'Identify, extracts and taxonomically classifies 16S '
                   'rRNA genes in genomes against LTP') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __gtdb_genome_path_file(grp, required=True)
            __gtdb_domain_file(grp, required=True)
            __ltp_version(grp, required=True)
            __ssu_version(grp, required=True)
            __log_file(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            __rna_file_path(grp)
            __silent(grp)
            __cpus(grp)
            __all_genomes(grp)
            __remove(grp, 'LTP')

    with subparser(sub_parsers, 'trnascan',
                   'Identifies tRNAs in genomes') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __gtdb_genome_path_file(grp, required=True)
            __gbk_bac_assembly_file(grp, required=True)
            __gbk_arc_assembly_file(grp, required=True)
            __rfq_arc_assembly_file(grp, required=True)
            __rfq_bac_assembly_file(grp, required=True)
            __log_file(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            __cpus(grp)
            __silent(grp)
            __all_genomes(grp)

    with subparser(sub_parsers, 'join_checkm',
                   'Join checkm output file for different versions of GTDB') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __checkm_files(grp, required=True)
            __output_file(grp, required=True)
            __log_file(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            __cpus(grp)
            __silent(grp)

    with subparser(sub_parsers, 'checkm',
                   'Run CheckM on new and modified genomes') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __gtdb_genome_path_file(grp, required=True)
            __gtdb_comparison_report(grp, required=True)
            __output_dir(grp, required=True)
            __log_file(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            __all_genomes(grp)
            __cpus(grp)
            __silent(grp)

    with subparser(sub_parsers, 'busco', 'Estimate quality of new fungal genomes') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __gtdb_genome_path_file(grp, required=True)
            __report(grp, required=True)
            __output_dir(grp, required=True)
            __log_file(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            __all_genomes(grp)
            __cpus(grp)
            __silent(grp)



    with subparser(sub_parsers, 'prepare_checkm2',
                      'Prepare files to run CheckM2 for the new release') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __checkm_summary_genbank(grp, required=True)
            __checkm_summary_refseq(grp, required=True)
            __gtdb_genome_path_file(grp, required=True)
            __metadata_file(grp, required=True)
            __output_dir(grp, required=True)
            __log_file(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            __silent(grp)
            __cpus(grp)

    with subparser(sub_parsers, 'join_checkm2',
                        'Join CheckM2 output files for different batches') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __checkm2_output_dir(grp, required=True)
            __output_dir(grp, required=True)
            __log_file(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            __silent(grp)
            __cpus(grp)
    with subparser(sub_parsers, 'update_db',
                   'Update the Postgres database.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __database_setup(grp, required=True)
            __checkm_profile(grp, required=True)
            __ftp_download_date(grp, required=True)
            __gtdb_genome_path_file(grp, required=True)
            __log_file(grp, required=True)
            __repository(grp, required=True)
            __report_dir(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            __silent(grp)
            __cpus(grp)

    # Commands to Update CheckM value in DB
    with subparser(sub_parsers, 'update_checkm_db', 'Update the CheckM value in Postgres database.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __database_setup(grp, required=True)
            __checkm_profile(grp, required=True)
            __checkm_qa(grp, required=True)
            __metadata_file(grp, required=True)
            __log_file(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            __silent(grp)

    # Commands to Update CheckM value in DB
    with subparser(sub_parsers, 'update_metadata_db', 'Update metadata information in the database.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __database_setup(grp, required=True)
            __log_file(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            __metadata_input_folder(grp)
            __metadata_table(grp)
            __metadata_table_description(grp)
            __do_not_null_field(grp)
            __genome_list(grp, required=False)
            __silent(grp)

    with subparser(sub_parsers,'update_reps_db', 'Update the species cluster representatives in database.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __database_setup(grp, required=True)
            __final_cluster_file(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            __silent(grp)

    with subparser(sub_parsers, 'update_type_designation', 'Update type_designation columns when all Seqcode,NCBI and LPSN infos are in the db.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __database_setup(grp, required=True)
            __log_file(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            __silent(grp)


    with subparser(sub_parsers, 'update_ncbitax_db', 'Update Organism name in the database.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __database_setup(grp, required=True)
            __organism_names(grp, required=True)
            __filtered_taxonomy(grp, required=True)
            __unfiltered_taxonomy(grp, required=True)
            __log_file(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            __do_not_null_field(grp)
            __genome_list(grp, required=False)
            __silent(grp)

    # Search Pfam Tigrfam markers
    with subparser(sub_parsers, 'hmmsearch', 'Call Hmmsearch on new and modified genomes.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __gtdb_genome_path_file(grp, required=True)
            __marker_db(grp, required=True)
            __report_file(grp, required=True)
            __hmm_db_path(grp, required=True)
            __log_file(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            __folder_suffix(grp)
            __cpus(grp)
            __tmp_dir(grp)
            __silent(grp)

    # Create metadata tables
    with subparser(sub_parsers, 'create_tables', 'Create tables with metadata for all NCBI genomes.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __gtdb_genome_path_file(grp, required=True)
            __output_dir(grp, required=True)
            __silva_version(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            __silent(grp)

    # Create metadata tables from NCBI assemblies
    with subparser(sub_parsers, 'parse_assemblies',
                   'Parse NCBI assembly summary files to generate metadata.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __rfq_bac_assembly_file(grp, required=True)
            __rfq_arc_assembly_file(grp, required=True)
            __gbk_bac_assembly_file(grp, required=True)
            __gbk_arc_assembly_file(grp, required=True)
            __metadata_file(grp, required=True)
            __output_file(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            __silent(grp)

    # Parse GTDB directory to generate extra NCBI metadata
    with subparser(sub_parsers, 'parse_ncbi_dir', 'Parse GTDB directory to generate extra NCBI metadata.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __gtdb_genome_path_file(grp, required=True)
            __output_file(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            __silent(grp)
            __cpus(grp)

    # # Parse NCBI Taxonomy files
    with subparser(sub_parsers, 'list_genomes', 'Produce file indicating the directory of each genome.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __genome_directory(grp, required=True)
            __output_file(grp, required=True)
            __gtdb_selected_genomes(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            # the mirror is on NFS, so the walk is bound by round trip latency
            # rather than by CPU; a handful of threads saturates the client
            __cpus(grp, default=8)
            __silent(grp)

    with subparser(sub_parsers, 'ncbi_genome_sync', 'Sync NCBI data to local directory.') as parser:
        # ncbi_genome_sync owns its whole interface, -l/--log included: it places its
        # .fail/.bad/.rm outputs beside the log, so the flag is defined where it is used
        add_sync_arguments(parser)

    with subparser(sub_parsers, 'update_genomes',
                   'Update RefSeq and GenBank genomes from the NCBI FTP mirror.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __ftp_directory(grp, required=True)
            __new_directory(grp, required=True)
            __ftp_genome_dirs_file(grp, required=True)
            # not required with --fresh, which reads no previous release; enforced
            # in main.py, argparse having no way to say one argument requires another
            __old_genome_dirs_file(grp, required=False)
        with arg_group(parser, 'options arguments') as grp:
            __silent(grp)
            __fresh(grp)
            __resume(grp)
            __dry_run(grp)
            # neither half of the work is bound by the CPU: copying a genome is
            # round trips to the file server, and a FASTA is read only where NCBI
            # reissued one. Measured against a mirror on NFS, 200 genomes of the
            # size NCBI serves take 17.1s copied one at a time and 4.5s copied 16
            # at a time, and a comparison in which every FASTA had been reissued
            # 46.0s against 4.3s. Past 16 the copying is bound by the link and the
            # hashing turns back down -- 64 was slower than 32 -- so it does not
            # follow the core count
            __cpus(grp, default=16)

    # # Steps to propagate GTDB Taxonomy
    with subparser(sub_parsers, 'propagate_gtdb_taxonomy', 'Propagating GTDB taxonomy to new release.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __gtdb_metadata_previous_release(grp, required=True)
            __gtdb_metadata_current_release(grp, required=True)
            __taxonomy_file(grp, required=True)
            __representative_file(grp, required=True)
            __log_file(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            __silent(grp)

    with subparser(sub_parsers, 'propagate_curated_taxonomy', 'Propagate curated taxonomy from representative genomes to all genomes in clusters.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __taxonomy_file(grp, required=True)
            __metadata_file(grp,required=True)
            __output_file(grp,required=True)
        with arg_group(parser, 'options arguments') as grp:
            __silent(grp)

    with subparser(sub_parsers, 'add_taxonomy_to_database', 'Update the taxonomy in the database, this is usually '
                                                            'used after propagate_taxonomy_from_reps_to_cluster '
                                                            'function') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __database_setup(grp, required=True)
            __taxonomy_file(grp, required=True)
            __metadata_file(grp,required=True)
        with arg_group(parser, 'options arguments') as grp:
            __silent(grp)
            __truncate_taxonomy(grp)


    # # Update surveillance genome list
    with subparser(sub_parsers, 'add_surveillance_genomes', 'Add surveillance genome to a table in GTDB.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __database_setup(grp, required=True)
            __surveillance_list(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            __silent(grp)

    with subparser(sub_parsers, 'add_names_dmp', 'Parse dmp file to a table') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __gbk_arc_assembly_file(grp, required=True)
            __gbk_bac_assembly_file(grp, required=True)
            __rfq_arc_assembly_file(grp, required=True)
            __rfq_bac_assembly_file(grp, required=True)
            __output_file(grp, required=True)
            __taxonomy_directory(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            __silent(grp)

    with subparser(sub_parsers, 'update_taxid_to_db', 'add taxid for each rank of each genomes to generate link to ncbi') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __database_setup(grp, required=True)
            __input_file(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            __silent(grp)

    with subparser(sub_parsers, 'generate_ltp_db',
        'Generate blast database based on information from the LTP website.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __csv_file(grp, required=True)
            __compressed_fasta(grp, required=True)
            __fasta_file(grp, required=True)
            __output_prefix(grp, required=True)
            __output_dir(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            __silent(grp)




    with subparser(sub_parsers, 'update_propagated_tax',
                   'Push changed from propagated taxonomy to new database.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __database_setup(grp, required=True)
            __taxonomy_file(grp, required=True)
            __metadata_file(grp, required=True)
            __genome_list(grp, required=True)
            __representative_file(grp, required=True)
            __log_file(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            __truncate_taxonomy(grp)
            __silent(grp)

    with subparser(sub_parsers, 'set_gtdb_domain',
                   'Set missing GTDB domain information to reflect NCBI domain.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __database_setup(grp, required=True)
            __log_file(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            __silent(grp)

    with subparser(sub_parsers, 'ncbi_genome_category',
                   'Identify genomes marked by NCBI as being a MAG or SAG.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __genbank_assembly_summary(grp, required=True)
            __refseq_assembly_summary(grp, required=True)
            __gtdb_genome_path_file(grp, required=True)
            __output_file(grp, required=True)
            __log_file(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            __silent(grp)
            __cpus(grp)

    # Generate metadata table for genomes in Seqcode
    with subparser(sub_parsers, 'generate_seqcode_table', 'Generate metadata table for genomes in Seqcode.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __output_dir(grp, required=True)
            __gtdb_genome_path_file(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            __silent(grp)
            __cpus(grp)

    with subparser(sub_parsers, 'lpsn', 'Steps to update LPSN Metadata.') as lpsn_parser:
        lpsn_sub_parsers = lpsn_parser.add_subparsers(help="--", dest='lpsn_subparser_name')

        with subparser(lpsn_sub_parsers, 'lpsn_wf', 'Full Pipeline Pull HTML and Parse HTML.') as parser:
            with arg_group(parser, 'required named arguments') as grp:
                __output_dir(grp, required=True)
                __lpsn_gss_file(grp, required=True)
                __unfiltered_taxonomy(grp, required=True)
            with arg_group(parser, 'options arguments') as grp:
                __silent(grp)

        with subparser(lpsn_sub_parsers, 'pull_html', 'Get files from LPSN listing all genera and species.') as parser:
            with arg_group(parser, 'required named arguments') as grp:
                __output_dir(grp, required=True)
            with arg_group(parser, 'options arguments') as grp:
                __skip_taxa_per_letter_dl(grp)
                __silent(grp)

        with subparser(lpsn_sub_parsers, 'parse_html', 'Parse HTML files.') as parser:
            with arg_group(parser, 'required named arguments') as grp:
                __input_dir(grp, required=True)
                __output_dir(grp, required=True)
                __lpsn_gss_file(grp, required=True)
            with arg_group(parser, 'options arguments') as grp:
                __silent(grp)

        with subparser(lpsn_sub_parsers, 'add_metadata', 'Add a lot of LPSN metadata to Database as a separate table') as parser:
            with arg_group(parser, 'required named arguments') as grp:
                __database_setup(grp, required=True)
                __lpsn_metadata_file(grp, required=True)
            with arg_group(parser, 'options arguments') as grp:
                __silent(grp)

    with subparser(sub_parsers, 'bacdive', 'Steps to update DSMZ Metadata.') as dsmz_parser:
        dsmz_sub_parsers = dsmz_parser.add_subparsers(help="--", dest='bacdive_subparser_name')

        with subparser(dsmz_sub_parsers, 'download_strains',
                       'Produce metadata files describing type genera, species, '
                       'and strains according to DSMZ.') as parser:
            with arg_group(parser, 'required named arguments') as grp:
                __username(grp, required=True)
                __password(grp, required=True)
                __output_dir(grp, required=True)
            with arg_group(parser, 'options arguments') as grp:
                __silent(grp)

    # # NCBI Strain
    with subparser(sub_parsers, 'ncbi_strains', 'NCBI Strain Parser.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __gtdb_genome_path_file(grp, required=True)
            __gbk_bac_assembly_file(grp, required=True)
            __gbk_arc_assembly_file(grp, required=True)
            __rfq_bac_assembly_file(grp, required=True)
            __rfq_arc_assembly_file(grp, required=True)
            __log_file(grp, required=True)
            __output_dir(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            __silent(grp)
            __cpus(grp)

    # # Tools to combine Medata from different databases.
    with subparser(sub_parsers, 'strains', 'Tools to combine metadata from different databases.') as strains_parser:
        strains_sub_parsers = strains_parser.add_subparsers(help="--", dest='strains_subparser_name')

        with subparser(strains_sub_parsers, 'date_table',
                       'Generate table with LPSN year or priority for species and subspecies names.') as parser:
            with arg_group(parser, 'required named arguments') as grp:
                __lpsn_scraped_species_info(grp, required=True)
                __lpsn_gss_file(grp, required=True)
                __output_file(grp, required=True)
            with arg_group(parser, 'options arguments') as grp:
                __silent(grp)

        with subparser(strains_sub_parsers, 'type_table',
                       'Generate table indicating metadata from LPSN for each species name.') as parser:
            with arg_group(parser, 'required named arguments') as grp:
                __lpsn_gss_file(grp, required=True)
                __lpsn_directory(grp, required=True)
                __year_table_file(grp, required=True)
                __metadata_file(grp, required=True)
                __ncbi_names(grp, required=True)
                __ncbi_nodes(grp, required=True)
                __output_dir(grp, required=True)
            with arg_group(parser, 'options arguments') as grp:
                __silent(grp)
                __cpus(grp)

    # # Overview
    with subparser(sub_parsers, 'overview',
                   'Compare the Metadata file from the previous version with the new one.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __previous_metadata_file(grp, required=True)
            __new_metadata_file(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            __only_ncbi(grp)
            __use_formatted_id(grp)
            __silent(grp)

    # # Compare Field of interest between metadata files
    with subparser(sub_parsers, 'compare_field',
                   'Compare a specific metadata field between to metadata files.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __previous_metadata_file(grp, required=True)
            __new_metadata_file(grp, required=True)
            __field_of_interest(grp, required=True)
            __output_file(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            __only_ncbi(grp)
            __silent(grp)

    with subparser(sub_parsers, 'compare_markers', 'Compare marker frequencies between 2 releases.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __first_domain_report(grp, required=True)
            __second_domain_report(grp, required=True)
            __output_file(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            __only_ncbi(grp)
            __use_formatted_id(grp)
            __silent(grp)

    # Lists and pseudo-trees for new representatives, polyphyletic taxa, rogue genomes,
    # and genomes with modified NCBI names
    with subparser(sub_parsers, 'curation_lists',
                   'Lists and pseudo-trees for new representatives, polyphyletic taxa, '
                   'rogue genomes, and genomes with modified NCBI names.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __gtdb_init_taxonomy(grp, required=True)
            __gtdb_sp_clusters(grp, required=True)
            __gtdb_prev_sp_clusters(grp, required=True)
            __gtdb_decorate_table(grp, required=True)
            __domain(grp, required=True)
            __output_dir(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            __silent(grp)

    with subparser(sub_parsers, 'check_unique_strains',
                   'Lists and pseudo-trees for new representatives, polyphyletic taxa, '
                   'rogue genomes, and genomes with modified NCBI names.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __node(grp, required=True)
            __name(grp, required=True)
            __metadata_file(grp, required=True)
            __output_file(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            __silent(grp)

    with subparser(sub_parsers, 'check_db_population',
                      'Check if the database is populated correctly.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __metadata_file(grp, required=True)
            __log_file(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            __id_last_genome(grp, required=True)
            __silent(grp)


    with subparser(sub_parsers, 'compare_metadata_genome_dir',
                   'Compare list of genomes in metadata file and genome_dirs') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __metadata_file(grp, required=True)
            __gtdb_genome_path_file(grp, required=True)
        with arg_group(parser, 'options arguments') as grp:
            __silent(grp)

    return main_parser


def main():
    # -------------------------------------------------
    # get and check options
    software_name = 'gtdb_migration_tk'
    if len(sys.argv) == 1:
        print_help()
        sys.exit(0)
    elif sys.argv[1] in {'-v', '--v', '-version', '--version'}:
        print(f"{software_name}: version %s %s %s" % (__version__,
                                                       __copyright__,
                                                       __author__))
        sys.exit(0)
    elif sys.argv[1] in {'-h', '--h', '-help', '--help'}:
        print_help()
        sys.exit(0)
    else:
        args = get_main_parser().parse_args()

        silent = False
        if hasattr(args, 'silent'):
            silent = args.silent

        if args.subparser_name in ('select_genomes', 'ncbi_metadata_sync'):
            # these have no --log flag: the log is part of the output of the command,
            # recording what was downloaded, and which genomes the selection rejected
            # and why; logger_setup() creates the directory before opening the file
            args.log = os.path.join(args.output_dir, 'gtdb_migration_tk.log')

        try:
            # dirname('sync.log') is '' and logger_setup() treats a falsy directory as
            # "no log file", so a bare filename silently produced no log at all
            logger_setup(os.path.dirname(args.log) or '.',
                         os.path.basename(args.log),
                         'GTDB Migration Tk',
                         software_name,
                         __version__,
                         silent)
        except:
            logger_setup('.',
                         'gtdb_migration_tk.log',
                         'GTDB Migration Tk',
                         software_name,
                         __version__,
                         silent)

        # do what we came here to do
        rtn_code = 0
        try:
            parser = OptionsParser()
            if False:
                import cProfile

                cProfile.run('parser.parse_options(args)', 'prof')
            elif False:
                import pdb

                pdb.run(parser.parse_options(args))
            else:
                rtn_code = parser.parse_options(args)
        except SystemExit:
            print("\n  Controlled exit resulting from an unrecoverable error or warning.")
        except:
            print("\nUnexpected error:", sys.exc_info()[0])
            raise

        # commands that report a meaningful exit code (ncbi_genome_sync) must not have it
        # swallowed; everything else returns 0 and exits normally
        if rtn_code:
            sys.exit(rtn_code)



if __name__ == '__main__':
    main()
