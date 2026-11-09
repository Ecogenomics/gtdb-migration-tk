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

"""Offline unit tests for metadata_export_manager.py: export_metadata.

'gtdb metadata export' assigned a representative to every genome without one,
and committed, before it exported; wrote the rows in whatever order the view
gave; and left a truncated file where an export failed part way. export_metadata
only reads, in one REPEATABLE READ transaction, warns of the genomes without a
representative, orders the rows by accession, and writes a file only once it is
complete. The database is stood in for by a fake cursor whose COPY writes the
text it is given: what is tested is the statement sent, and what becomes of the
file."""

import gzip
import logging
import os
import shutil
import tempfile
import unittest
from unittest import mock

from gtdb_migration_tk import metadata_export_manager as E

VIEW_COLUMNS = ['id', 'accession', 'formatted_accession', 'gtdb_representative', 'gtdb_taxonomy']
EXPORT = ('accession\tformatted_accession\tgtdb_representative\tgtdb_taxonomy\n'
          'GB_GCA_000000002.1\tG000000002\tf\tnone\n'
          'RS_GCF_000000001.1\tG000000001\tt\td__Bacteria;s__A\n')


class FakeCursor(object):
    """A cursor over a database of three genomes, one of them without a gtdb_representative.

    Its COPY writes EXPORT, or writes part of it and raises where fail_copy.
    """

    def __init__(self, no_representative=('GB_GCA_000000003.1',), fail_copy=False):
        self.no_representative = list(no_representative)
        self.fail_copy = fail_copy
        self.statements = []
        self.description = None
        self.result = []
        self.rowcount = -1

    def execute(self, sql, params=None):
        self.statements.append(sql)
        if sql.startswith('SELECT * FROM metadata_view LIMIT 0'):
            self.description = [(name,) for name in VIEW_COLUMNS]
            self.result = []
        elif sql == E.NUM_GENOMES:
            self.result = [(3,)]
        elif sql == E.NO_REPRESENTATIVE:
            self.result = [(genome,) for genome in self.no_representative]

    def fetchone(self):
        return self.result[0]

    def fetchall(self):
        return list(self.result)

    def copy_expert(self, sql, handle):
        self.statements.append(sql)
        if self.fail_copy:
            handle.write(EXPORT[:40])
            raise RuntimeError('the connection was lost')
        handle.write(EXPORT)
        self.rowcount = 2


class FakeConnection(object):
    def __init__(self):
        self.commits = 0
        self.rollbacks = 0

    def commit(self):
        self.commits += 1

    def rollback(self):
        self.rollbacks += 1


class ExportingMetadata(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='metadata_export_manager_test.')
        self.addCleanup(shutil.rmtree, self.dir, True)
        logging.getLogger('timestamp').addHandler(logging.NullHandler())
        patcher = mock.patch.object(E, 'log_directory', return_value=self.dir)
        patcher.start()
        self.addCleanup(patcher.stop)

    def manager(self, cursor=None):
        manager = E.MetadataExportManager.__new__(E.MetadataExportManager)
        manager.logger = logging.getLogger('timestamp')
        manager.temp_cur = cursor or FakeCursor()
        manager.temp_con = FakeConnection()
        return manager

    def export(self, name='metadata.tsv', out_format='tab', cursor=None):
        manager = self.manager(cursor)
        path = os.path.join(self.dir, name)
        with self.assertLogs('timestamp', level='INFO') as logged:
            exported = manager.export_metadata(path, out_format)
        return manager, path, exported, logged.records

    def copy_statement(self, manager):
        (copy,) = [sql for sql in manager.temp_cur.statements if sql.startswith('COPY')]
        return copy

    def test_every_column_of_the_view_but_its_id_is_exported_ordered_by_accession(self):
        manager, path, exported, _ = self.export()
        self.assertEqual(self.copy_statement(manager),
                         'COPY (SELECT "accession", "formatted_accession", "gtdb_representative", '
                         '"gtdb_taxonomy" FROM metadata_view ORDER BY "accession") TO STDOUT WITH '
                         "(FORMAT csv, HEADER, DELIMITER E'\\t', NULL 'none')")
        with open(path) as handle:
            self.assertEqual(handle.read(), EXPORT)
        self.assertEqual(exported, 2)

    def test_tab_is_the_default_and_csv_writes_null_empty_as_gtdb_did(self):
        self.assertEqual(E.DEFAULT_FORMAT, 'tab')
        manager, _, _, _ = self.export(out_format='csv')
        self.assertTrue(self.copy_statement(manager).endswith('TO STDOUT WITH (FORMAT csv, HEADER)'))

    def test_an_unknown_format_or_a_view_without_accession_is_refused(self):
        with self.assertRaises(ValueError):
            E.export_query(VIEW_COLUMNS, 'xlsx')
        with self.assertRaises(ValueError):
            E.export_query(['id', 'formatted_accession'], 'tab')

    def test_nothing_is_written_to_the_database_and_the_transaction_is_ended(self):
        manager, _, _, _ = self.export()
        writes = [sql for sql in manager.temp_cur.statements
                  if sql.split()[0].upper() in ('UPDATE', 'INSERT', 'DELETE', 'CREATE', 'DROP', 'ALTER', 'TRUNCATE')]
        self.assertEqual(writes, [])
        self.assertEqual((manager.temp_con.commits, manager.temp_con.rollbacks), (0, 1))

    def test_the_connection_is_read_only_and_repeatable_read(self):
        # the count of genomes without a representative is of the snapshot exported
        with mock.patch.object(E, 'GenomeDatabaseConnectionFTPUpdate') as connection:
            E.MetadataExportManager({'service': 'gtdb_r237'})
        connection.return_value.conn.set_session.assert_called_once_with(
            isolation_level=E.ISOLATION_LEVEL_REPEATABLE_READ, readonly=True)

    def test_a_genome_without_a_representative_is_warned_of_listed_and_still_exported(self):
        # gtdb metadata export gave it a representative by AAI, and committed
        _, path, _, records = self.export()
        warnings = [r.getMessage() for r in records if r.levelno == logging.WARNING]
        self.assertEqual(len(warnings), 1)
        self.assertIn('1 of the 3 genomes of the database have no gtdb_representative', warnings[0])
        self.assertIn('GB_GCA_000000003.1', warnings[0])
        with open(os.path.join(self.dir, E.NO_REPRESENTATIVE_NAME)) as handle:
            self.assertEqual(handle.read(), 'genome_id\nGB_GCA_000000003.1\n')
        self.assertTrue(os.path.exists(path))

    def test_every_genome_with_a_representative_is_said_with_no_warning(self):
        _, _, _, records = self.export(cursor=FakeCursor(no_representative=()))
        self.assertFalse([r for r in records if r.levelno >= logging.WARNING])
        self.assertIn('Every one of the 3 genomes of the database has a gtdb_representative.',
                      [r.getMessage() for r in records])
        self.assertFalse(os.path.exists(os.path.join(self.dir, E.NO_REPRESENTATIVE_NAME)))

    def test_a_gz_output_is_gzipped_and_the_same_bytes_each_time(self):
        _, path, _, _ = self.export(name='metadata.tsv.gz')
        with gzip.open(path, 'rt') as handle:
            self.assertEqual(handle.read(), EXPORT)
        with open(path, 'rb') as handle:
            first = handle.read()
        self.export(name='metadata.tsv.gz')
        with open(path, 'rb') as handle:
            self.assertEqual(handle.read(), first)

    def test_the_output_directory_is_made(self):
        _, path, _, _ = self.export(name=os.path.join('exports', 'r237', 'metadata.tsv'))
        self.assertTrue(os.path.isfile(path))

    def test_a_failed_export_leaves_no_file_and_an_earlier_export_as_it_was(self):
        path = os.path.join(self.dir, 'metadata.tsv')
        with open(path, 'w') as handle:
            handle.write('an earlier export\n')
        manager = self.manager(FakeCursor(fail_copy=True))
        with self.assertRaises(RuntimeError):
            manager.export_metadata(path, 'tab')

        with open(path) as handle:
            self.assertEqual(handle.read(), 'an earlier export\n')
        self.assertEqual(sorted(os.listdir(self.dir)), ['export_no_gtdb_representative.tsv', 'metadata.tsv'])
        self.assertEqual(manager.temp_con.rollbacks, 1)

    def test_the_counts_are_logged(self):
        _, path, _, records = self.export()
        self.assertIn('Exported the metadata of 2 genomes, 4 columns of metadata_view, as tab, to {}'.format(path),
                      ' '.join(r.getMessage() for r in records))


class TheCommandLine(unittest.TestCase):
    def parse(self, *args):
        from gtdb_migration_tk import __main__ as main_module
        return main_module.get_main_parser().parse_args(['export_metadata'] + list(args))

    def test_the_command_takes_the_database_an_output_file_and_a_log(self):
        options = self.parse('--db_service', 'gtdb_r237', '-o', 'metadata.tsv.gz', '-l', 'export.log')
        self.assertEqual((options.db_service, options.output_file, options.log, options.out_format),
                         ('gtdb_r237', 'metadata.tsv.gz', 'export.log', 'tab'))
        self.assertEqual(self.parse('--db_service', 'gtdb_r237', '-o', 'm.csv', '-l', 'e.log',
                                    '--format', 'csv').out_format, 'csv')

    def test_the_output_file_and_log_are_required_and_the_format_one_of_two(self):
        for args in (['--db_service', 'gtdb_r237', '-l', 'e.log'],
                     ['--db_service', 'gtdb_r237', '-o', 'm.tsv'],
                     ['--db_service', 'gtdb_r237', '-o', 'm.tsv', '-l', 'e.log', '--format', 'xlsx']):
            with mock.patch('sys.stderr'), self.assertRaises(SystemExit) as ended:
                self.parse(*args)
            self.assertEqual(ended.exception.code, 2, args)

    def test_the_command_exports_through_the_manager(self):
        from gtdb_migration_tk import main as main_py
        options = self.parse('--db_service', 'gtdb_r237', '-o', 'metadata.tsv.gz', '-l', 'export.log',
                             '--format', 'csv')
        with mock.patch.object(main_py, 'MetadataExportManager') as manager:
            main_py.OptionsParser().parse_options(options)
        manager.assert_called_once_with({'service': 'gtdb_r237'})
        manager.return_value.export_metadata.assert_called_once_with('metadata.tsv.gz', 'csv')

    def test_the_help_names_the_command(self):
        from gtdb_migration_tk import __main__ as main_module
        with mock.patch('builtins.print') as printed:
            main_module.print_help()
        text = '\n'.join(str(call.args[0]) for call in printed.call_args_list if call.args)
        self.assertIn('export_metadata -> Export the metadata of every genome in the database', text)


if __name__ == '__main__':
    unittest.main()
