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
import psycopg2
from gtdb_migration_tk.database_configuration import GenomeDatabaseConnectionFTPUpdate


class GTDBImporter(object):

    def __init__(self,temp_cur):

        self.temp_cur = temp_cur

    def importMetadata(self, table=None, field=None, typemeta=None, data_list=None):
        '''
        Function importMetadata
        import one field of Metadata for a list of Genomes
        :param table: Table where the column is located
        :param field: Name of the Column
        :param typemeta: Data type of the column
        :param metafile: TSV file with the format (Genome_id \t Value)
        '''
        try:
            data_zip = list(zip(*data_list))
            genome_id = list(data_zip[0])
            meta_value = list(data_zip[1])
            for n, i in enumerate(genome_id):
                new_i = i.split("_", 1)[1]
                genome_id[n] = new_i
            query = "SELECT upsert('{0}','{1}','{2}',%s,%s)".format(
                table, field, typemeta)

            self.temp_cur.execute(query, (genome_id, meta_value))
        except Exception as e:
            print(e)