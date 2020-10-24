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
import sys
import argparse
import tempfile
import logging
import glob

from collections import defaultdict

from gtdb_migration_tk.gtdb_lite.gtdb_importer import GTDBImporter
from gtdb_migration_tk.biolib_lite.taxonomy import Taxonomy

from gtdb_migration_tk.database_configuration import GenomeDatabaseConnectionFTPUpdate

class MetadataDatabaseManager(object):

    def __init__(self,hostname,user,password,db):
        """Initialization."""
        self.logger = logging.getLogger('timestamp')
        self.description_table = {'metadata_gene.tsv':['metadata_gene.desc.tsv'],
                                  'metadata_nt.tsv':['metadata_nt.desc.tsv'],
                                  'metadata_ssu_gg.tsv':['metadata_rna.table.desc.tsv'],
                                  'metadata_ssu_silva.tsv':['metadata_rna.table.desc.tsv','metadata_sequence.desc.tsv'],
                                  'metadata_lsu_silva_23s.tsv':['metadata_rna.table.desc.tsv','metadata_sequence.desc.tsv'],
                                  'metadata_lsu_5S.tsv':['metadata_rna.table.desc.tsv','metadata_sequence.desc.tsv'],
                                  'metadata_ssu_silva_count.tsv':['metadata_ssu_count.desc.tsv'],
                                  'metadata_lsu_silva_23s_count.tsv':['metadata_ssu_count.desc.tsv'],
                                  'metadata_lsu_5S_count.tsv':['metadata_ssu_count.desc.tsv'],
                                  'metadata_trna_count.tsv':['metadata_trna.desc.tsv'],
                                  'ncbi_assembly_summary.tsv':['metadata_ncbi_assembly_file.desc.tsv'],
                                  'strain_summary_file.tsv':['metadata_ncbi_assembly.desc.tsv','metadata_ncbi_assembly_file.desc.tsv'],
                                  'ncbi_assembly_metadata.tsv':['metadata_ncbi_assembly.desc.tsv']
                                  }

        self.password = password
        self.hostname = hostname
        self.user = user
        self.db = db

        self.temp_con = GenomeDatabaseConnectionFTPUpdate.GenomeDatabaseConnectionFTPUpdate(
            hostname, user, password, db)
        self.temp_con.MakePostgresConnection()
        self.temp_cur = self.temp_con.cursor()

    def process_metadata_files(self,genome_list_file,do_not_null_field=False,table_folder=None,table_file=None,table_file_desc=None):
        file_dir = os.path.dirname(os.path.realpath(__file__))
        desc_table_dir = os.path.join(file_dir, 'data_files', 'table_description')
        if table_folder is not None:
            list_tsv_files = glob.glob(os.path.join(table_folder, '*.tsv'))
            for tsv_file in list_tsv_files:
                if os.path.basename(tsv_file) not in self.description_table:
                    print(f'{os.path.basename(tsv_file)} is not a standard table')
                    sys.exit(-1)
            for tsv_file in list_tsv_files:
                for desc_file in self.description_table.get(os.path.basename(tsv_file)):
                    metadata_desc_file = os.path.join(desc_table_dir,desc_file)
                    self.update_metadata_db(tsv_file,metadata_desc_file,genome_list_file,do_not_null_field)
        elif table_file is not None:
            self.logger.info("Processing specific file")
            self.update_metadata_db(table_file, table_file_desc, genome_list_file, do_not_null_field)

    def update_metadata_db(self,metadata_file,metadata_desc_file,genome_list_file,do_not_null_field):
        # get fields in metadata file
        gtdbimporter = GTDBImporter(self.temp_cur)
        with open(metadata_file) as f:
            metadata_fields = f.readline().strip().split('\t')[1:]
        self.logger.info(
            'Metadata file contains {} fields.'.format(len(metadata_fields)))
        self.logger.info('Fields: %s' % ', '.join(metadata_fields))

        # get database table and data type of each metadata field
        metadata_type = {}
        metadata_table = {}
        with open(metadata_desc_file) as f:
            for line in f:
                line_split = line.strip('\n').split('\t')
                field = line_split[0]
                if field in metadata_fields:
                    metadata_type[field] = line_split[2]
                    metadata_table[field] = line_split[3]
        self.logger.info('Identified {} matching fields in metadata description file.'.format(
            len(metadata_table)))
        self.logger.info('Fields: %s' % ', '.join(metadata_table))

        # set fields to NULL if requested
        if not do_not_null_field:
            response = ''
            while response.lower() not in ['y', 'n']:
                response = input(
                    "Set fields to NULL for all genomes [y/n]: ")

            if response.lower() == 'y':

                for field in metadata_table:
                    q = ("UPDATE {} SET {} = NULL".format(
                        metadata_table[field], field))
                    print(q)
                    self.temp_cur.execute(q)
                self.temp_con.commit()

            elif response.lower() == 'n':
                pass
            else:
                self.logger.error('Unrecognized input.')
                sys.exit(-1)

        # get genomes to process
        genome_list = set()
        if genome_list_file:
            for line in open(genome_list_file):
                if len(line.split('\t')) >= len(line.split(',')):
                    genome_list.add(line.rstrip().split('\t')[0])
                else:
                    genome_list.add(line.rstrip().split(',')[0])

        # read metadata file
        metadata = defaultdict(lambda: defaultdict(str))
        with open(metadata_file) as f:
            fields = [x.strip() for x in f.readline().split('\t')]

            for line in f:
                line_split = line.rstrip('\n').split('\t')

                genome_id = line_split[0]
                # print line_split
                for i, value in enumerate(line_split[1:]):
                    metadata[fields[i + 1]][genome_id] = value

        # add each field to the database
        for field in metadata:
            data_to_commit = []

            if field not in metadata_type:
                continue

            data_type = metadata_type[field]
            table = metadata_table[field]

            records_to_update = 0
            for orig_genome_id, value in metadata[field].items():

                try:
                    if float(value) and data_type in ['INT', 'INTEGER']:
                        # assume specified data type is correct and that we may need
                        # to cast floats to integers
                        value = str(int(float(value)))
                except:
                    pass

                if value.strip():
                    genome_id = str(orig_genome_id)
                    if genome_id.startswith('GCA_'):
                        genome_id = 'GB_' + genome_id
                    elif genome_id.startswith('GCF_'):
                        genome_id = 'RS_' + genome_id

                    if (not genome_list
                            or genome_id in genome_list
                            or orig_genome_id in genome_list):
                        data_to_commit.append((genome_id, value))
                        records_to_update += 1

            self.logger.info('Updating {} for {} genomes.'.format(
                field, records_to_update))

            gtdbimporter.importMetadata(table, field, data_type, data_to_commit)
            self.temp_con.commit()


class NCBITaxDatabaseManager(object):
    """Add organism name to GTDB."""

    def __init__(self,hostname,user,password,db):
        """Initialization."""
        self.logger = logging.getLogger('timestamp')

        self.password = password
        self.hostname = hostname
        self.user = user
        self.db = db

        self.temp_con = GenomeDatabaseConnectionFTPUpdate.GenomeDatabaseConnectionFTPUpdate(
            hostname, user, password, db)
        self.temp_con.MakePostgresConnection()
        self.temp_cur = self.temp_con.cursor()


    # set fields to NULL if requested
    def set_field_to_null(self,metadata_table,field):
        response = ''
        while response.lower() not in ['y', 'n']:
            response = input(
                "Set fields to NULL for all genomes [y/n]: ")

        if response.lower() == 'y':
            q = ("UPDATE {} SET {} = NULL".format(
                metadata_table, field))
            print(q)
            self.temp_cur.execute(q)
            self.temp_con.commit()

        elif response.lower() == 'n':
            pass
        else:
            self.logger.error('Unrecognized input.')
            sys.exit(-1)

    def update_ncbitax_db(self, organism_name_file,filtered_file,unfiltered_file, genome_list_file,do_not_null_field=False):
        """Add organism name to database."""
        gtdbimporter = GTDBImporter(self.temp_cur)

        genome_list = set()
        data_to_commit = []
        if genome_list_file:
            for line in open(genome_list_file):
                if '\t' in line:
                    genome_list.add(line.rstrip().split('\t')[0])
                else:
                    genome_list.add(line.rstrip().split(',')[0])

        # add full taxonomy string to database
        records_to_update =0
        for line in open(organism_name_file):
            line_split = line.strip().split('\t')

            gid = line_split[0]
            org_name = line_split[1]
            if genome_list_file and gid not in genome_list:
                continue

            data_to_commit.append((gid, org_name))
            records_to_update += 1

        if not do_not_null_field:
            self.set_field_to_null('metadata_ncbi', 'ncbi_organism_name')
        self.logger.info('Updating {} for {} genomes.'.format(
            'ncbi_organism_name', records_to_update))
        gtdbimporter.importMetadata('metadata_ncbi', 'ncbi_organism_name', 'TEXT', data_to_commit)
        self.temp_con.commit()

        taxonomy = Taxonomy().read(filtered_file)
        data_filtered_to_commit = []
        records_to_update =0
        # add full taxonomy string to database
        for genome_id, taxa in taxonomy.items():
            if genome_id.startswith('GCA_'):
                genome_id = 'GB_' + genome_id
            elif genome_id.startswith('GCF_'):
                genome_id = 'RS_' + genome_id

            if genome_list_file and genome_id not in genome_list:
                continue
            taxa_str = ';'.join(taxa)
            data_filtered_to_commit.append((genome_id, taxa_str))
            records_to_update += 1
        if not do_not_null_field:
            self.set_field_to_null('metadata_taxonomy', 'ncbi_taxonomy')
        self.logger.info('Updating {} for {} genomes.'.format(
            'ncbi_taxonomy', records_to_update))
        gtdbimporter.importMetadata('metadata_taxonomy', 'ncbi_taxonomy', 'TEXT', data_filtered_to_commit)
        self.temp_con.commit()

        # read taxonomy file
        unfiltered_taxonomy = Taxonomy().read(unfiltered_file)
        data_unfiltered_to_commit = []

        # add full taxonomy string to database
        records_to_update =0
        for genome_id, taxa in unfiltered_taxonomy.items():
            if genome_id.startswith('GCA_'):
                genome_id = 'GB_' + genome_id
            elif genome_id.startswith('GCF_'):
                genome_id = 'RS_' + genome_id

            if genome_list_file and genome_id not in genome_list:
                continue

            taxa_str = ';'.join(taxa)
            data_unfiltered_to_commit.append((genome_id, taxa_str))
            records_to_update += 1
        if not do_not_null_field:
            self.set_field_to_null('metadata_taxonomy', 'ncbi_taxonomy_unfiltered')
        self.logger.info('Updating {} for {} genomes.'.format(
            'ncbi_taxonomy_unfiltered', records_to_update))
        gtdbimporter.importMetadata('metadata_taxonomy', 'ncbi_taxonomy_unfiltered', 'TEXT', data_unfiltered_to_commit)
        self.temp_con.commit()



