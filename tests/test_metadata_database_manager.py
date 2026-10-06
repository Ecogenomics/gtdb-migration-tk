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

"""update_metadata_db reads a metadata table gzipped or not: strains type_table
writes gtdb_type_strain_summary.tsv.gz. The database is stood in for by mocks,
the test being of what is read and handed to the importer."""

import gzip
import logging
import os
import shutil
import tempfile
import unittest
from unittest import mock

from gtdb_migration_tk import metadata_database_manager as M

SUMMARY = ('accession\tgtdb_type_designation_ncbi_taxa\tlpsn_priority_year\n'
           'RS_GCF_000000001.1\ttype strain of species\t1919\n'
           'GB_GCA_000000002.1\tnot type material\t\n')
DESCRIPTION = ('gtdb_type_designation_ncbi_taxa\tdesc\tTEXT\tmetadata_type_material\n'
               'lpsn_priority_year\tdesc\tINTEGER\tmetadata_type_material\n')


class LoadingATable(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='metadata_database_manager_test.')
        self.addCleanup(shutil.rmtree, self.dir, True)
        logging.getLogger('timestamp').addHandler(logging.NullHandler())
        self.description = os.path.join(self.dir, 'desc.tsv')
        with open(self.description, 'w') as handle:
            handle.write(DESCRIPTION)

    def load(self, metadata_file):
        manager = M.MetadataDatabaseManager.__new__(M.MetadataDatabaseManager)
        manager.logger = logging.getLogger('timestamp')
        manager.temp_cur, manager.temp_con = mock.Mock(), mock.Mock()
        with mock.patch.object(M, 'GTDBImporter') as importer:
            manager.update_metadata_db(metadata_file, self.description, None, True)
        return {call.args[1]: sorted(call.args[3])
                for call in importer.return_value.import_metadata_to_db.call_args_list}

    def test_a_gzipped_table_is_read_as_a_plain_one_is(self):
        plain = os.path.join(self.dir, 'gtdb_type_strain_summary.tsv')
        with open(plain, 'w') as handle:
            handle.write(SUMMARY)
        gzipped = plain + '.gz'
        with gzip.open(gzipped, 'wt') as handle:
            handle.write(SUMMARY)

        loaded = self.load(gzipped)
        self.assertEqual(loaded, self.load(plain))
        self.assertEqual(loaded['gtdb_type_designation_ncbi_taxa'],
                         [('GB_GCA_000000002.1', 'not type material'),
                          ('RS_GCF_000000001.1', 'type strain of species')])


if __name__ == '__main__':
    unittest.main()
