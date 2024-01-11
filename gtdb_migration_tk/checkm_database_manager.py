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

from gtdb_migration_tk.gtdb_lite.gtdb_importer import GTDBImporter
from gtdb_migration_tk.database_configuration import GenomeDatabaseConnectionFTPUpdate


class CheckMDatabaseManager(object):

    def __init__(self,hostname,user,password,db):
        """Initialization."""
        self.logger = logging.getLogger('timestamp')
        self.metadata = {'Completeness': ['checkm_completeness', 'FLOAT'],
                         'Contamination': ['checkm_contamination', 'FLOAT'],
                         'Strain heterogeneity': ['checkm_strain_heterogeneity', 'FLOAT'],
                         'Marker lineage': ['checkm_marker_lineage', 'TEXT'],
                         '# genomes': ['checkm_genome_count', 'INT'],
                         '# markers': ['checkm_marker_count', 'INT'],
                         '# marker sets': ['checkm_marker_set_count', 'INT']}

        self.password = password
        self.hostname = hostname
        self.user = user
        self.db = db

        self.temp_con = GenomeDatabaseConnectionFTPUpdate.GenomeDatabaseConnectionFTPUpdate(
            hostname, user, password, db)
        self.temp_con.MakePostgresConnection()
        self.temp_cur = self.temp_con.cursor()


    def add_checkm_to_db(self,checkm_profile_file, checkm_qa_sh100_file, genome_list_file):
        # get genomes to process
        genome_list = set()
        gtdbimporter = GTDBImporter(self.temp_cur)
        if genome_list_file:
            for line in open(genome_list_file):
                if line[0] == '#':
                    continue
                if len(line.split('\t')) >= len(line.split(',')):
                    genome_list.add(line.rstrip().split('\t')[0])
                else:
                    genome_list.add(line.rstrip().split(',')[0])

        # add CheckM profile fields
        for header, data in self.metadata.items():
            db_header, data_type = data

            num_genomes = 0
            data_to_commit = []
            with open(checkm_profile_file) as f:
                headers = f.readline().rstrip().split('\t')
                col_index = headers.index(header)

                for line in f:
                    line_split = line.rstrip().split('\t')
                    orig_genome_id = line_split[0]
                    genome_id = orig_genome_id[0:orig_genome_id.find('_', 4)]

                    if genome_id.startswith('GCA_'):
                        genome_id = 'GB_' + genome_id
                    elif genome_id.startswith('GCF_'):
                        genome_id = 'RS_' + genome_id

                    if genome_id not in genome_list:
                        print('Skipping genome: {}'.format(genome_id))
                        continue

                    data = line_split[col_index]
                    data_to_commit.append(((genome_id, data)))
                    num_genomes += 1

            if num_genomes == 0:
                print('No genomes identified.')
                sys.exit(-1)

            gtdbimporter.importMetadata('metadata_genes', db_header, data_type, data_to_commit)

        # add strain heterogeneity results at 100%
        data_to_commit = []
        with open(checkm_qa_sh100_file) as f:
            headers = f.readline().rstrip().split('\t')
            sh_index = headers.index('Strain heterogeneity')

            for line in f:
                line_split = line.rstrip().split('\t')
                genome_id = line_split[0]
                genome_id = genome_id[0:genome_id.find('_', 4)]

                if genome_id.startswith('GCA_'):
                    genome_id = 'GB_' + genome_id
                elif genome_id.startswith('GCF_'):
                    genome_id = 'RS_' + genome_id

                if genome_id not in genome_list:
                    continue

                sh = line_split[sh_index]
                if '949373205' in genome_id:
                    print('949373205: {}'.format(sh))
                data_to_commit.append(((genome_id, sh)))

        db_header = 'checkm_strain_heterogeneity_100'
        data_type = 'FLOAT'
        gtdbimporter.importMetadata('metadata_genes', db_header, data_type, data_to_commit)
        self.temp_con.commit()


