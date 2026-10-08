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
the test being of what is read and handed to the importer.

And each metadata command is one transaction: committed once, when it is done,
or rolled back. They committed as they went, so a run that failed part way left
some fields of the new release and some of the old, or a field set to NULL for
every genome with its new values never written. None of them asks a question
on the terminal, which a run under nohup could not answer, but update_metadata_db
given --genome_list without --do_not_null_field: it then removes the metadata of
every genome and writes it again only for those listed.

update_metadata_db loads only genomes the database holds unless given a list,
reads a table once, in chunks, and refuses a table with a malformed row, a
genome named twice or an INT field that is not a whole number."""

import gzip
import logging
import os
import shutil
import tempfile
import unittest
from collections import defaultdict
from unittest import mock

from gtdb_migration_tk import metadata_database_manager as M
from gtdb_migration_tk import ncbi_genome_category as GENOME_CATEGORY
from gtdb_migration_tk import ncbi_strain_summary as NCBI_STRAINS
from gtdb_migration_tk.gtdb_lite.gtdb_importer import SKIP, UnknownGenomesError

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
            importer.return_value.genomes.return_value = {'GCF_000000001.1', 'GCA_000000002.1'}
            manager.process_metadata_files(None, do_not_null_field=True, table_file=metadata_file,
                                           table_file_desc=self.description)
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
                         [('GCA_000000002.1', 'not type material'),
                          ('GCF_000000001.1', 'type strain of species')])


class FakeCursor(object):
    """A cursor recording each statement, over the rows of metadata_view."""

    def __init__(self, accessions=(), fail_on=None, genomes=()):
        self.accessions = list(accessions)
        self.fail_on = fail_on
        self.genomes = list(genomes)
        self.statements = []
        self.upserts = []
        self.copied = []
        self.result = []
        self.rowcount = 0

    def copy_expert(self, sql, handle):
        self.statements.append(sql)
        self.copied.append(sorted(line for line in handle.read().splitlines() if line))

    def resets(self):
        """(statement, the genomes it was told were written) of each reset_unwritten()."""
        resets = [sql for sql in self.statements if sql.startswith('UPDATE') and 'NULL' in sql]
        return list(zip(resets, self.copied))

    def written(self):
        """(table, field) -> sorted [(genome, value), ...], every chunk of each joined."""
        fields = defaultdict(list)
        for table, field, _type, genomes, values in self.upserts:
            fields[(table, field)].extend(zip(genomes, values))
        return {key: sorted(rows) for key, rows in fields.items()}

    def execute(self, sql, params=None):
        self.statements.append(sql)
        if self.fail_on and self.fail_on in sql:
            raise RuntimeError('the server refused ' + self.fail_on)
        if 'FROM metadata_view' in sql:
            self.result = [(accession,) for accession in self.accessions]
        elif sql == 'SELECT id_at_source FROM genomes':
            self.result = [(genome,) for genome in self.genomes]
        elif sql.startswith('SELECT upsert('):
            self.upserts.append(params)

    def fetchall(self):
        return list(self.result)

    def executemany(self, sql, rows):
        self.statements.append(sql)
        self.rows = list(rows)

    def __iter__(self):
        return iter(self.result)


class FakeConnection(object):
    def __init__(self):
        self.commits = 0
        self.rollbacks = 0

    def commit(self):
        self.commits += 1

    def rollback(self):
        self.rollbacks += 1


class OneTransaction(unittest.TestCase):
    """A manager built without a database, its cursor and connection fakes."""

    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='metadata_database_manager_test.')
        self.addCleanup(shutil.rmtree, self.dir, True)
        logging.getLogger('timestamp').addHandler(logging.NullHandler())
        # no command may ask a question: a run under nohup has no one to answer
        patcher = mock.patch('builtins.input', side_effect=AssertionError('asked a question'))
        patcher.start()
        self.addCleanup(patcher.stop)

    def manager(self, cls=M.MetadataDatabaseManager, cursor=None):
        manager = cls.__new__(cls)
        manager.logger = logging.getLogger('timestamp')
        manager.temp_cur = cursor or FakeCursor()
        manager.temp_con = FakeConnection()
        return manager

    def write(self, name, text):
        path = os.path.join(self.dir, name)
        with open(path, 'w') as handle:
            handle.write(text)
        return path


class LoadingMetadataInOneTransaction(OneTransaction):
    def load(self, manager, **kwargs):
        table = self.write('metadata_type_material.tsv', SUMMARY)
        description = self.write('desc.tsv', DESCRIPTION)
        manager.process_metadata_files(None, table_file=table, table_file_desc=description, **kwargs)

    @staticmethod
    def importer():
        patcher = mock.patch.object(M, 'GTDBImporter')
        importer = patcher.start()
        importer.return_value.genomes.return_value = {'GCF_000000001.1', 'GCA_000000002.1'}
        return patcher, importer

    def test_fields_are_set_to_null_and_written_in_one_commit(self):
        manager = self.manager()
        patcher, _ = self.importer()
        self.addCleanup(patcher.stop)
        self.load(manager)

        # each field is set to NULL for the genomes it holds a value of that the
        # table did not give one: GB_GCA_000000002.1 has no lpsn_priority_year
        self.assertEqual(manager.temp_cur.resets(), [
            (M.RESET_UNWRITTEN.format(table='metadata_type_material', field='gtdb_type_designation_ncbi_taxa'),
             ['GCA_000000002.1', 'GCF_000000001.1']),
            (M.RESET_UNWRITTEN.format(table='metadata_type_material', field='lpsn_priority_year'),
             ['GCF_000000001.1'])])
        self.assertEqual((manager.temp_con.commits, manager.temp_con.rollbacks), (1, 0))

    def test_do_not_null_field_sets_nothing_to_null(self):
        manager = self.manager()
        patcher, _ = self.importer()
        self.addCleanup(patcher.stop)
        self.load(manager, do_not_null_field=True)

        self.assertFalse([sql for sql in manager.temp_cur.statements if 'NULL' in sql])

    def test_a_field_that_cannot_be_written_rolls_back_the_nulls_before_it(self):
        # the NULLs were committed first, leaving the field NULL for every genome
        manager = self.manager()
        patcher, importer = self.importer()
        self.addCleanup(patcher.stop)
        importer.return_value.import_metadata_to_db.side_effect = [None, RuntimeError('refused')]
        with self.assertRaises(RuntimeError):
            self.load(manager)

        self.assertEqual((manager.temp_con.commits, manager.temp_con.rollbacks), (0, 1))

    def test_a_table_that_is_not_standard_ends_the_run_with_a_rollback(self):
        folder = os.path.join(self.dir, 'tables')
        os.makedirs(folder)
        with open(os.path.join(folder, 'not_a_standard_table.tsv'), 'w') as handle:
            handle.write('accession\tx\n')
        manager = self.manager()
        manager.description_table = {'metadata_gene.tsv': ['metadata_gene.desc.tsv']}
        with self.assertRaisesRegex(M.MetadataTableError, 'not_a_standard_table.tsv'):
            manager.process_metadata_files(None, table_folder=folder)

        self.assertEqual(manager.temp_cur.statements, [])
        self.assertEqual((manager.temp_con.commits, manager.temp_con.rollbacks), (0, 1))


NCBI_TABLE = ('genome_id\tncbi_taxid\tncbi_organism_name\tncbi_contig_l50\n'
              'RS_GCF_000000001.1\t562\tEscherichia coli\t3\n'
              'GB_GCA_000000002.1\t1280.0\tStaphylococcus aureus\t1\n'
              'GB_GCA_999999999.1\t9606\tnot in the database\t2\n')
NCBI_DESCRIPTION = ('ncbi_taxid\tNCBI taxonomy identifier.\tINT\tmetadata_ncbi\n'
                    'ncbi_organism_name\tName of organism.\tTEXT\tmetadata_ncbi\n')
HELD = ('GCF_000000001.1', 'GCA_000000002.1')


class LoadingOnlyWhatTheDatabaseHolds(OneTransaction):
    """A table loaded through the real importer, the database a cursor holding HELD."""

    def setUp(self):
        super().setUp()
        # the importer lists genomes it refuses beside the log
        handler = logging.FileHandler(os.path.join(self.dir, 'run.log'))
        logging.getLogger('timestamp').addHandler(handler)
        self.addCleanup(logging.getLogger('timestamp').removeHandler, handler)
        self.addCleanup(handler.close)

    def load(self, table=NCBI_TABLE, description=NCBI_DESCRIPTION, genome_list=None, genomes=HELD, **kwargs):
        manager = self.manager(cursor=FakeCursor(genomes=genomes))
        table_file = self.write('ncbi_assembly_metadata.tsv', table)
        description_file = self.write('desc.tsv', description)
        with self.assertLogs('timestamp', level='INFO') as logged:
            manager.process_metadata_files(genome_list, table_file=table_file, table_file_desc=description_file,
                                           **kwargs)
        return manager, [record.getMessage() for record in logged.records]

    def test_without_a_genome_list_genomes_the_database_does_not_hold_are_skipped_not_refused(self):
        # parse_ncbi_assemblies writes every genome of NCBI's summaries
        manager, messages = self.load()

        self.assertEqual(manager.temp_cur.written(), {
            ('metadata_ncbi', 'ncbi_taxid'): [('GCA_000000002.1', '1280'), ('GCF_000000001.1', '562')],
            ('metadata_ncbi', 'ncbi_organism_name'): [('GCA_000000002.1', 'Staphylococcus aureus'),
                                                      ('GCF_000000001.1', 'Escherichia coli')]})
        self.assertIn('2 loaded, 1 not in the genomes table and skipped', '\n'.join(messages))
        self.assertEqual((manager.temp_con.commits, manager.temp_con.rollbacks), (1, 0))

    def test_with_a_genome_list_only_its_genomes_are_loaded(self):
        genome_list = self.write('genomes.tsv', 'accession\tx\nGB_GCA_000000002.1\tx\n')
        manager, _ = self.load(genome_list=genome_list, do_not_null_field=True)

        self.assertEqual(manager.temp_cur.written()[('metadata_ncbi', 'ncbi_taxid')],
                         [('GCA_000000002.1', '1280')])

    def test_a_listed_genome_the_database_does_not_hold_is_still_refused(self):
        genome_list = self.write('genomes.tsv', 'GCA_999999999.1\nGCF_000000001.1\n')
        manager = self.manager(cursor=FakeCursor(genomes=HELD))
        table = self.write('ncbi_assembly_metadata.tsv', NCBI_TABLE)
        description = self.write('desc.tsv', NCBI_DESCRIPTION)
        with self.assertRaisesRegex(UnknownGenomesError, 'GCA_999999999.1'):
            manager.process_metadata_files(genome_list, do_not_null_field=True, table_file=table,
                                           table_file_desc=description)
        self.assertEqual((manager.temp_con.commits, manager.temp_con.rollbacks), (0, 1))

    def test_columns_in_no_description_are_named_in_the_log(self):
        _, messages = self.load()
        self.assertIn('1 column(s) of', '\n'.join(messages))
        self.assertIn('are in no description and are not loaded: ncbi_contig_l50.', '\n'.join(messages))

    def test_a_whole_number_written_as_a_float_is_written_as_an_integer(self):
        # '0.0' passed float(value)'s truth test as false and went to the INT field as it was
        table = 'genome_id\tncbi_taxid\nGCF_000000001.1\t0.0\nGCA_000000002.1\t-7\n'
        manager, _ = self.load(table=table)

        self.assertEqual(manager.temp_cur.written()[('metadata_ncbi', 'ncbi_taxid')],
                         [('GCA_000000002.1', '-7'), ('GCF_000000001.1', '0')])

    def test_a_whole_number_too_big_for_an_int_field_is_left_unwritten_and_listed_beside_the_log(self):
        # r237's GCA_964261755.1, 9.5 Gbp of MAGs as one assembly, failed upsert()'s cast to INT
        # and rolled back the run with a traceback
        table = ('genome_id\tncbi_taxid\nGCF_000000001.1\t9528631298\nGCA_000000002.1\t2147483647\n'
                 'GCA_000000003.1\t-2147483649\n')
        with mock.patch.object(M, 'log_directory', return_value=self.dir):
            manager, messages = self.load(table=table, genomes=HELD + ('GCA_000000003.1',))

        self.assertEqual(manager.temp_cur.written()[('metadata_ncbi', 'ncbi_taxid')],
                         [('GCA_000000002.1', '2147483647')])
        # each is left without a value, so set to NULL where it held one
        self.assertEqual([written for sql, written in manager.temp_cur.resets() if 'ncbi_taxid' in sql],
                         [['GCA_000000002.1']])
        self.assertTrue(any('gives 2 value(s) too big for an INT field, which are not written' in m
                            and 'GCF_000000001.1 ncbi_taxid 9528631298' in m for m in messages))
        with open(os.path.join(self.dir, 'int_out_of_range.ncbi_assembly_metadata.tsv')) as handle:
            self.assertEqual(handle.read().splitlines(), ['genome_id\tfield\tvalue',
                                                          'GCF_000000001.1\tncbi_taxid\t9528631298',
                                                          'GCA_000000003.1\tncbi_taxid\t-2147483649'])
        self.assertEqual((manager.temp_con.commits, manager.temp_con.rollbacks), (1, 0))

    def refused(self, table, pattern, description=NCBI_DESCRIPTION):
        manager = self.manager(cursor=FakeCursor(genomes=HELD))
        table_file = self.write('ncbi_assembly_metadata.tsv', table)
        description_file = self.write('desc.tsv', description)
        with self.assertRaisesRegex(M.MetadataTableError, pattern):
            manager.process_metadata_files(None, table_file=table_file, table_file_desc=description_file)
        self.assertEqual(manager.temp_cur.upserts, [])
        self.assertEqual((manager.temp_con.commits, manager.temp_con.rollbacks), (0, 1))

    def test_an_int_field_that_is_not_a_whole_number_refuses_the_table(self):
        # 12.5 was written as 12
        self.refused('genome_id\tncbi_taxid\nGCF_000000001.1\t12.5\n', "line 2 gives GCF_000000001.1 ncbi_taxid the value '12.5'")
        self.refused('genome_id\tncbi_taxid\nGCF_000000001.1\tna\n', "'na', which is not a whole number")

    def test_a_short_or_long_row_refuses_the_table_naming_its_line(self):
        # a short row left its last fields unset; a long one raised an IndexError
        self.refused('genome_id\tncbi_taxid\tncbi_organism_name\nGCF_000000001.1\t562\n',
                     r'line 2 \(GCF_000000001.1\) has 2 column\(s\) where its header has 3')
        self.refused('genome_id\tncbi_taxid\nGCF_000000001.1\t562\tx\n',
                     r'line 2 \(GCF_000000001.1\) has 3 column\(s\) where its header has 2')

    def test_a_genome_named_twice_refuses_the_table(self):
        # the last row was kept, whichever run it was of
        self.refused('genome_id\tncbi_taxid\nGCF_000000001.1\t562\nRS_GCF_000000001.1\t563\n',
                     'names GCF_000000001.1 more than once, again on line 3')

    def test_a_description_line_short_of_a_table_refuses_the_run(self):
        self.refused('genome_id\tncbi_taxid\nGCF_000000001.1\t562\n', 'line 1 has 3 column',
                     description='ncbi_taxid\tNCBI taxonomy identifier.\tINT\n')

    def test_a_field_is_set_to_null_after_it_is_written_and_only_for_genomes_it_was_not_written_for(self):
        # setting it to NULL for every genome first rewrote every row of the table
        # once more, upsert() rewriting the rows of every genome given a value
        manager, messages = self.load()

        statements = manager.temp_cur.statements
        upserts = [i for i, sql in enumerate(statements) if sql.startswith('SELECT upsert(')]
        resets = [i for i, sql in enumerate(statements) if sql.startswith('UPDATE') and 'NULL' in sql]
        self.assertEqual(len(resets), 2)
        self.assertLess(max(upserts), min(resets))
        self.assertTrue(all('IS NOT NULL' in statements[i] and 'NOT EXISTS' in statements[i] for i in resets))
        self.assertEqual([written for _, written in manager.temp_cur.resets()],
                         [['GCA_000000002.1', 'GCF_000000001.1']] * 2)
        self.assertIn('Set metadata_ncbi.ncbi_taxid to NULL for 0 genome(s) holding a value this table did '
                      'not give them.', messages)

    def test_a_table_is_written_in_chunks_and_the_chunks_are_the_whole_table(self):
        genomes = ['GCF_{:09d}.1'.format(i) for i in range(25)]
        table = 'genome_id\tncbi_taxid\n' + ''.join('{}\t{}\n'.format(gid, i) for i, gid in enumerate(genomes))
        with mock.patch.object(M, 'CHUNK_GENOMES', 10):
            chunked, _ = self.load(table=table, genomes=genomes)
        whole, _ = self.load(table=table, genomes=genomes)

        self.assertEqual(len(chunked.temp_cur.upserts), 3)
        self.assertEqual(len(whole.temp_cur.upserts), 1)
        self.assertEqual(chunked.temp_cur.written(), whole.temp_cur.written())
        self.assertEqual((chunked.temp_con.commits, chunked.temp_con.rollbacks), (1, 0))


class ChoosingTheTables(OneTransaction):
    def test_neither_both_or_a_table_without_its_description_is_refused_before_anything(self):
        manager = self.manager()
        table = self.write('t.tsv', 'genome_id\tx\n')
        for kwargs, pattern in (({}, 'not neither'),
                                ({'table_folder': self.dir, 'table_file': table}, 'not both'),
                                ({'table_file': table}, 'without --metadata_table_desc')):
            with self.assertRaisesRegex(M.MetadataTableError, pattern):
                manager.process_metadata_files(None, **kwargs)
        self.assertEqual(manager.temp_cur.statements, [])

    def test_every_table_of_a_folder_this_command_does_not_know_is_named(self):
        folder = os.path.join(self.dir, 'tables')
        os.makedirs(folder)
        for name in ('metadata_gene.tsv', 'one.tsv', 'two.tsv'):
            open(os.path.join(folder, name), 'w').close()
        manager = self.manager()
        manager.description_table = {'metadata_gene.tsv': ['metadata_gene.desc.tsv']}

        with self.assertRaisesRegex(M.MetadataTableError, 'holds 2 table.*: one.tsv, two.tsv'):
            manager.process_metadata_files(None, table_folder=folder)

    def gzipped_folder(self, *names):
        folder = os.path.join(self.dir, 'tables')
        os.makedirs(folder)
        for name in names:
            path = os.path.join(folder, name)
            text = ('genome_id\tncbi_strain_identifiers\tncbi_type_material_designation\n'
                    'GCF_000000001.1\tK-12\tassembly from type material\n')
            with (gzip.open(path, 'wt') if name.endswith('.gz') else open(path, 'w')) as handle:
                handle.write(text)
        manager = self.manager(cursor=FakeCursor(genomes=HELD))
        with mock.patch.object(M.GenomeDatabaseConnectionFTPUpdate, 'GenomeDatabaseConnectionFTPUpdate'):
            manager.description_table = M.MetadataDatabaseManager({}).description_table
        return folder, manager

    def test_a_folder_of_gzipped_tables_is_loaded_as_one_of_plain_ones_is(self):
        # create_tables, parse_ncbi_assemblies, parse_ncbi_dir and ncbi_strains write them gzipped
        folder, manager = self.gzipped_folder('strain_summary_file.tsv.gz')
        manager.process_metadata_files(None, table_folder=folder)

        self.assertEqual(manager.temp_cur.written()[('metadata_ncbi', 'ncbi_strain_identifiers')],
                         [('GCF_000000001.1', 'K-12')])
        self.assertEqual((manager.temp_con.commits, manager.temp_con.rollbacks), (1, 0))

    def test_a_table_both_gzipped_and_not_is_refused_before_anything(self):
        # each would be loaded, the second over the first
        folder, manager = self.gzipped_folder('strain_summary_file.tsv', 'strain_summary_file.tsv.gz')
        with self.assertRaisesRegex(M.MetadataTableError, 'both gzipped and not: strain_summary_file.tsv'):
            manager.process_metadata_files(None, table_folder=folder)
        self.assertEqual(manager.temp_cur.statements, [])

    def test_a_gzipped_table_it_does_not_know_is_named_as_it_is(self):
        folder, manager = self.gzipped_folder('not_a_table.tsv.gz')
        with self.assertRaisesRegex(M.MetadataTableError, 'does not know: not_a_table.tsv.gz'):
            manager.process_metadata_files(None, table_folder=folder)

    def test_the_substrains_ncbi_strains_writes_beside_its_table_are_passed_over_not_refused(self):
        # r237 wrote ncbi_strains into the folder update_metadata_db loads
        folder, manager = self.gzipped_folder('strain_summary_file.tsv.gz')
        with open(os.path.join(folder, NCBI_STRAINS.SUBSTRAINS_NAME), 'w') as handle:
            handle.write('genome_id\tstrain_id\nGCF_000000001.1\tstrain=K-12 substr. MG1655\n')
        manager.process_metadata_files(None, table_folder=folder)

        self.assertIn(NCBI_STRAINS.SUBSTRAINS_NAME, M.NOT_TABLES)
        self.assertEqual(manager.temp_cur.written()[('metadata_ncbi', 'ncbi_strain_identifiers')],
                         [('GCF_000000001.1', 'K-12')])

    def test_the_genome_category_table_is_known_and_its_evidence_passed_over(self):
        # parse_ncbi_genome_category wrote to whatever file -o named, loaded by hand
        folder = os.path.join(self.dir, 'tables')
        os.makedirs(folder)
        with gzip.open(os.path.join(folder, GENOME_CATEGORY.CATEGORY_TABLE_NAME + '.gz'), 'wt') as handle:
            handle.write('genome_id\tncbi_genome_category\tsource\n'
                         'GCF_000000001.1\tderived from single cell\tGBFF file\n')
        with open(os.path.join(folder, GENOME_CATEGORY.EVIDENCE_NAME), 'w') as handle:
            handle.write('genome_id\tevidence\nGCF_000000001.1\t/note="single cell"\n')
        manager = self.manager(cursor=FakeCursor(genomes=HELD))
        with mock.patch.object(M.GenomeDatabaseConnectionFTPUpdate, 'GenomeDatabaseConnectionFTPUpdate'):
            manager.description_table = M.MetadataDatabaseManager({}).description_table

        manager.process_metadata_files(None, table_folder=folder)

        self.assertIn(GENOME_CATEGORY.EVIDENCE_NAME, M.NOT_TABLES)
        self.assertEqual(manager.temp_cur.written(), {('metadata_ncbi', 'ncbi_genome_category'):
                                                      [('GCF_000000001.1', 'derived from single cell')]})

    def test_the_strain_summary_is_read_once_against_both_its_descriptions(self):
        # it was read, and its fields set to NULL, once for each description
        folder = os.path.join(self.dir, 'tables')
        os.makedirs(folder)
        with open(os.path.join(folder, 'strain_summary_file.tsv'), 'w') as handle:
            handle.write('genome_id\tOrganism name\tncbi_strain_identifiers\tncbi_type_material_designation\n'
                         'GCF_000000001.1\tEscherichia coli\tK-12\tassembly from type material\n')
        manager = self.manager(cursor=FakeCursor(genomes=HELD))
        with mock.patch.object(M.GenomeDatabaseConnectionFTPUpdate, 'GenomeDatabaseConnectionFTPUpdate'):
            manager.description_table = M.MetadataDatabaseManager({}).description_table
        opened = []
        real_open_text = M.open_text
        with mock.patch.object(M, 'open_text', side_effect=lambda path: opened.append(path) or real_open_text(path)):
            manager.process_metadata_files(None, table_folder=folder)

        self.assertEqual(len(opened), 1)
        self.assertEqual(sorted(manager.temp_cur.written()), [('metadata_ncbi', 'ncbi_strain_identifiers'),
                                                             ('metadata_taxonomy', 'ncbi_type_material_designation')])
        nulls = [sql for sql in manager.temp_cur.statements if 'NULL' in sql]
        self.assertEqual(len(nulls), 2)


class TheDescriptions(unittest.TestCase):
    """The description files update_metadata_db loads a table against."""

    def descriptions(self):
        with mock.patch.object(M.GenomeDatabaseConnectionFTPUpdate, 'GenomeDatabaseConnectionFTPUpdate'):
            table = M.MetadataDatabaseManager({}).description_table
        folder = os.path.join(os.path.dirname(M.__file__), 'data_files', 'table_description')
        return {name: os.path.join(folder, name) for files in table.values() for name in files}, folder

    def test_every_description_a_known_table_is_loaded_against_exists(self):
        used, _ = self.descriptions()
        self.assertEqual([name for name, path in used.items() if not os.path.exists(path)], [])

    def test_no_description_names_a_greengenes_field(self):
        # the ssu_gg_* fields were dropped from metadata_rna in 0.1.58, and
        # nothing has written them since rna_silva took over
        _, folder = self.descriptions()
        for name in sorted(os.listdir(folder)):
            with open(os.path.join(folder, name)) as handle:
                fields = [line.split('\t')[0] for line in handle if line.strip()]
            self.assertEqual([field for field in fields if '_gg_' in field or field.startswith('gg_')], [], name)


class AskingBeforeAPartialLoad(unittest.TestCase):
    """--genome_list without --do_not_null_field removes every genome's metadata and writes the list's."""

    def test_yes_proceeds(self):
        for answer in ('y', 'Yes', ' YES '):
            with mock.patch('builtins.input', return_value=answer) as asked:
                M.confirm_partial_load('genomes.tsv')
            self.assertIn('genomes.tsv', asked.call_args.args[0])
            self.assertIn('every genome in the database', asked.call_args.args[0])

    def test_anything_else_ends_the_run_with_nothing_changed(self):
        for answer in ('n', '', 'maybe'):
            with mock.patch('builtins.input', return_value=answer), self.assertRaises(SystemExit) as ended:
                M.confirm_partial_load('genomes.tsv')
            self.assertIn('nothing was changed', str(ended.exception.code))

    def test_no_terminal_to_answer_ends_the_run_rather_than_waiting(self):
        with mock.patch('builtins.input', side_effect=EOFError), self.assertRaises(SystemExit) as ended:
            M.confirm_partial_load('genomes.tsv')
        self.assertIn('no answer', str(ended.exception.code))

    def test_the_command_asks_only_for_a_genome_list_without_do_not_null_field(self):
        from gtdb_migration_tk import __main__ as main_module
        from gtdb_migration_tk import main as main_py
        base = ['update_metadata_db', '--db_service', 'gtdb_r237', '-l', 'run.log', '-i', 'tables']
        for extra, asks in (([], False), (['--genome_list', 'g.tsv'], True),
                            (['--genome_list', 'g.tsv', '--do_not_null_field'], False)):
            options = main_module.get_main_parser().parse_args(base + extra)
            with mock.patch.object(main_py, 'confirm_partial_load') as confirm, \
                    mock.patch.object(main_py, 'MetadataDatabaseManager'):
                main_py.OptionsParser().parse_options(options)
            self.assertEqual(confirm.called, asks, extra)

    def test_the_genome_list_help_says_what_the_list_is(self):
        from gtdb_migration_tk import __main__ as main_module
        parser = main_module.get_main_parser()
        sub = [action for action in parser._actions if hasattr(action, 'choices') and action.choices
               and 'update_metadata_db' in action.choices][0].choices['update_metadata_db']
        helps = {action.dest: action.help for action in sub._actions}
        self.assertEqual(helps['genome_list'], 'Only process genomes in this list (e.g. metadata file exported from GTDB)')


class UpdatingRepresentatives(OneTransaction):
    VIEW = ['RS_GCF_000000001.1', 'GB_GCA_000000002.1', 'GB_GCA_000000003.1']

    def clusters(self, *rows):
        return self.write('clusters.tsv', 'Representative\tClustered genomes\n'
                          + ''.join('\t'.join(row) + '\n' for row in rows))

    def test_a_genome_of_the_cluster_file_not_in_the_database_is_refused_before_any_null(self):
        # the NULLs were committed before a KeyError on such a genome, leaving
        # the database with no representatives at all
        manager = self.manager(cursor=FakeCursor(self.VIEW))
        clusters = self.clusters(('G000000001', 'G000000002,G000000099'))
        with mock.patch.object(M, 'GTDBImporter') as importer:
            with self.assertRaises(UnknownGenomesError) as raised:
                manager.update_reps(clusters)

        self.assertIn('G000000099', str(raised.exception))
        self.assertFalse([sql for sql in manager.temp_cur.statements if 'NULL' in sql])
        importer.return_value.import_metadata_to_db.assert_not_called()
        self.assertEqual((manager.temp_con.commits, manager.temp_con.rollbacks), (0, 1))

    def test_the_representatives_are_written_in_one_commit(self):
        manager = self.manager(cursor=FakeCursor(self.VIEW))
        clusters = self.clusters(('G000000001', 'G000000002'), ('G000000003',))
        with mock.patch.object(M, 'GTDBImporter') as importer:
            manager.update_reps(clusters)

        written = {call.args[1]: sorted(call.args[3])
                   for call in importer.return_value.import_metadata_to_db.call_args_list}
        self.assertEqual(written['gtdb_genome_representative'],
                         [('GB_GCA_000000002.1', 'RS_GCF_000000001.1'),
                          ('GB_GCA_000000003.1', 'GB_GCA_000000003.1'),
                          ('RS_GCF_000000001.1', 'RS_GCF_000000001.1')])
        self.assertEqual(written['gtdb_representative'],
                         [('GB_GCA_000000002.1', 'False'), ('GB_GCA_000000003.1', 'True'),
                          ('RS_GCF_000000001.1', 'True')])
        self.assertEqual((manager.temp_con.commits, manager.temp_con.rollbacks), (1, 0))


class ReplacingTheSurveillanceGenomes(OneTransaction):
    def test_the_list_replaces_the_table_in_one_commit_blanks_and_repeats_dropped(self):
        manager = self.manager()
        listed = self.write('surveillance.txt', 'G000000001\n\nG000000002\nG000000001\n')
        manager.add_surveillance_genomes(listed)

        self.assertEqual(manager.temp_cur.statements[0], 'TRUNCATE survey_genomes')
        self.assertEqual(manager.temp_cur.rows, [('G000000001',), ('G000000002',)])
        self.assertEqual((manager.temp_con.commits, manager.temp_con.rollbacks), (1, 0))

    def test_a_list_that_cannot_be_inserted_leaves_the_table_as_it_was(self):
        manager = self.manager()
        manager.temp_cur.executemany = mock.Mock(side_effect=RuntimeError('refused'))
        with self.assertRaises(RuntimeError):
            manager.add_surveillance_genomes(self.write('surveillance.txt', 'G000000001\n'))

        self.assertEqual((manager.temp_con.commits, manager.temp_con.rollbacks), (0, 1))


class UpdatingNCBITaxonomy(OneTransaction):
    """update_ncbi_tax_db, through the real importer, the database a cursor holding HELD."""

    NAMES = ('RS_GCF_000000001.1\tEscherichia coli\n'
             'GB_GCA_000000002.1\tStaphylococcus aureus\n'
             'GB_GCA_999999999.1\tnot in the database\n')
    FILTERED = ('GCF_000000001.1\td__Bacteria; p__Pseudomonadota;\n'
                'GCA_999999999.1\td__Bacteria\n')
    UNFILTERED = ('GCF_000000001.1\td__Bacteria;x__Pseudomonadati;p__Pseudomonadota\n'
                  'GCA_000000002.1\td__Bacteria;p__Bacillota\n')

    def setUp(self):
        super().setUp()
        handler = logging.FileHandler(os.path.join(self.dir, 'run.log'))
        logging.getLogger('timestamp').addHandler(handler)
        self.addCleanup(logging.getLogger('timestamp').removeHandler, handler)
        self.addCleanup(handler.close)

    def run_update(self, names=NAMES, filtered=FILTERED, unfiltered=UNFILTERED, genome_list=None,
                   genomes=HELD, **kwargs):
        manager = self.manager(cls=M.NCBITaxDatabaseManager, cursor=FakeCursor(genomes=genomes))
        files = [self.write(name, text) for name, text in (('names.tsv', names), ('filtered.tsv', filtered),
                                                           ('unfiltered.tsv', unfiltered))]
        out_dir = os.path.join(self.dir, 'out')
        os.makedirs(out_dir, exist_ok=True)
        with self.assertLogs('timestamp', level='INFO') as logged:
            manager.update_ncbi_tax_db(*files, genome_list, out_dir, **kwargs)
        with open(os.path.join(out_dir, M.NCBI_TAX_MISSING_NAME)) as handle:
            missing = [line.split('\t') for line in handle.read().splitlines()]
        return manager, missing, logged.records

    def test_the_genomes_the_database_holds_are_written_and_the_rest_skipped_in_one_commit(self):
        # NCBI's files cover every assembly NCBI holds
        manager, _, _ = self.run_update()

        self.assertEqual(manager.temp_cur.written(), {
            ('metadata_ncbi', 'ncbi_organism_name'): [('GCA_000000002.1', 'Staphylococcus aureus'),
                                                      ('GCF_000000001.1', 'Escherichia coli')],
            ('metadata_taxonomy', 'ncbi_taxonomy'): [('GCF_000000001.1', 'd__Bacteria;p__Pseudomonadota')],
            ('metadata_taxonomy', 'ncbi_taxonomy_unfiltered'): [
                ('GCA_000000002.1', 'd__Bacteria;p__Bacillota'),
                ('GCF_000000001.1', 'd__Bacteria;x__Pseudomonadati;p__Pseudomonadota')]})
        self.assertEqual(len([sql for sql in manager.temp_cur.statements if 'NULL' in sql]), 3)
        self.assertEqual((manager.temp_con.commits, manager.temp_con.rollbacks), (1, 0))

    def test_each_field_is_set_to_null_only_for_genomes_its_file_gave_no_value(self):
        manager, _, _ = self.run_update()

        self.assertEqual(manager.temp_cur.resets(), [
            (M.RESET_UNWRITTEN.format(table='metadata_ncbi', field='ncbi_organism_name'),
             ['GCA_000000002.1', 'GCF_000000001.1']),
            (M.RESET_UNWRITTEN.format(table='metadata_taxonomy', field='ncbi_taxonomy'),
             ['GCF_000000001.1']),
            (M.RESET_UNWRITTEN.format(table='metadata_taxonomy', field='ncbi_taxonomy_unfiltered'),
             ['GCA_000000002.1', 'GCF_000000001.1'])])
        statements = manager.temp_cur.statements
        self.assertLess(max(i for i, sql in enumerate(statements) if sql.startswith('SELECT upsert(')),
                        max(i for i, sql in enumerate(statements) if sql.startswith('UPDATE')))

    def test_do_not_null_field_sets_nothing_to_null(self):
        manager, _, _ = self.run_update(do_not_null_field=True)
        self.assertEqual(manager.temp_cur.resets(), [])

    def test_each_genome_without_one_of_the_three_is_written_to_the_error_file_and_warned_of(self):
        _, missing, records = self.run_update()

        self.assertEqual(missing, [['genome_id', 'ncbi_organism_name', 'ncbi_taxonomy', 'ncbi_taxonomy_unfiltered'],
                                   ['GCA_000000002.1', '', 'missing', '']])
        warnings = [r.getMessage() for r in records if r.levelno == logging.WARNING]
        self.assertEqual(len(warnings), 1)
        self.assertIn('1 genome(s) of the genomes table have no ncbi_taxonomy, e.g. GCA_000000002.1', warnings[0])

    def test_an_empty_value_is_missing_not_written(self):
        # an empty taxonomy made Taxonomy().read() fail
        names = 'RS_GCF_000000001.1\t\nGB_GCA_000000002.1\tStaphylococcus aureus\n'
        filtered = 'GCF_000000001.1\t\nGCA_000000002.1\td__Bacteria\n'
        manager, missing, _ = self.run_update(names=names, filtered=filtered)

        self.assertNotIn(('GCF_000000001.1', ''), manager.temp_cur.written()[('metadata_ncbi', 'ncbi_organism_name')])
        self.assertEqual(missing[1:], [['GCF_000000001.1', 'missing', 'missing', '']])

    def test_a_genome_list_matches_every_file_whichever_way_it_names_a_genome(self):
        # the organism names were matched as written (RS_GCF_...) and the
        # taxonomies with a prefix added, so no list matched both
        for listed in ('GCF_000000001.1\n', 'RS_GCF_000000001.1\n'):
            genome_list = self.write('genomes.tsv', listed)
            manager, missing, _ = self.run_update(genome_list=genome_list, do_not_null_field=True)

            self.assertEqual(sorted(manager.temp_cur.written()), [
                ('metadata_ncbi', 'ncbi_organism_name'), ('metadata_taxonomy', 'ncbi_taxonomy'),
                ('metadata_taxonomy', 'ncbi_taxonomy_unfiltered')])
            self.assertTrue(all(rows == [('GCF_000000001.1', rows[0][1])]
                                for rows in manager.temp_cur.written().values()))
            self.assertEqual(missing[1:], [])

    def test_a_genome_named_twice_refuses_the_run(self):
        manager = self.manager(cls=M.NCBITaxDatabaseManager, cursor=FakeCursor(genomes=HELD))
        files = [self.write('names.tsv', self.NAMES),
                 self.write('filtered.tsv', 'GCF_000000001.1\td__A\nRS_GCF_000000001.1\td__B\n'),
                 self.write('unfiltered.tsv', self.UNFILTERED)]
        with self.assertRaisesRegex(M.MetadataTableError, 'names GCF_000000001.1 more than once, again on line 2'):
            manager.update_ncbi_tax_db(*files, None, self.dir)
        self.assertEqual((manager.temp_con.commits, manager.temp_con.rollbacks), (0, 1))

    def test_a_file_is_written_in_chunks_that_are_the_whole_file(self):
        genomes = ['GCF_{:09d}.1'.format(i) for i in range(25)]
        names = ''.join('{}\tname {}\n'.format(gid, i) for i, gid in enumerate(genomes))
        taxonomy = ''.join('{}\td__Bacteria\n'.format(gid) for gid in genomes)
        with mock.patch.object(M, 'CHUNK_GENOMES', 10):
            chunked, _, _ = self.run_update(names, taxonomy, taxonomy, genomes=genomes)
        whole, _, _ = self.run_update(names, taxonomy, taxonomy, genomes=genomes)

        self.assertEqual(len(chunked.temp_cur.upserts), 9)
        self.assertEqual(chunked.temp_cur.written(), whole.temp_cur.written())


class TheNCBITaxonomyCommandLine(unittest.TestCase):
    def options(self, *extra):
        from gtdb_migration_tk import __main__ as main_module
        return main_module.get_main_parser().parse_args(
            ['update_ncbi_tax_db', '--db_service', 'gtdb_r237', '-n', 'names.tsv', '--filtered', 'f.tsv',
             '--unfiltered', 'u.tsv', '-l', 'run.log'] + list(extra))

    def test_the_command_takes_the_organism_names_as_n_and_requires_an_out_dir(self):
        options = self.options('-o', 'out')
        self.assertEqual((options.organism_names, options.output_dir), ('names.tsv', 'out'))
        with mock.patch('sys.stderr'), self.assertRaises(SystemExit) as ended:
            self.options()
        self.assertEqual(ended.exception.code, 2)

    def test_update_ncbitax_db_is_no_longer_a_command(self):
        from gtdb_migration_tk import __main__ as main_module
        with mock.patch('sys.stderr'), self.assertRaises(SystemExit) as ended:
            main_module.get_main_parser().parse_args(
                ['update_ncbitax_db', '--db_service', 'gtdb_r237', '-o', 'names.tsv', '--filtered', 'f.tsv',
                 '--unfiltered', 'u.tsv', '-l', 'run.log'])
        self.assertEqual(ended.exception.code, 2)

    def test_the_command_asks_only_for_a_genome_list_without_do_not_null_field(self):
        from gtdb_migration_tk import main as main_py
        out_dir = tempfile.mkdtemp(prefix='metadata_database_manager_test.')
        self.addCleanup(shutil.rmtree, out_dir, True)
        for extra, asks in (([], False), (['--genome_list', 'g.tsv'], True),
                            (['--genome_list', 'g.tsv', '--do_not_null_field'], False)):
            options = self.options('-o', out_dir, *extra)
            with mock.patch.object(main_py, 'confirm_partial_load') as confirm, \
                    mock.patch.object(main_py, 'check_file_exists'), \
                    mock.patch.object(main_py, 'NCBITaxDatabaseManager') as manager:
                main_py.OptionsParser().parse_options(options)
            self.assertEqual(confirm.called, asks, extra)
            if asks:
                self.assertEqual(confirm.call_args.args, ('g.tsv', 'update_ncbi_tax_db'))
            manager.return_value.update_ncbi_tax_db.assert_called_once_with(
                'names.tsv', 'f.tsv', 'u.tsv', options.genome_list, out_dir, options.do_not_null_field)

    def test_the_prompt_names_the_command_asking(self):
        with mock.patch('builtins.input', return_value='n'), self.assertRaises(SystemExit) as ended:
            M.confirm_partial_load('g.tsv', 'update_ncbi_tax_db')
        self.assertTrue(str(ended.exception.code).startswith('update_ncbi_tax_db:'))


class DecidingTypeDesignations(unittest.TestCase):
    def test_a_species_valid_under_the_seqcode_is_type_strain_with_seqcode_a_source(self):
        designations, sources = M.type_designation_changes(
            [(1, 'Valid', None, 'not type material', 'NCBI;LPSN'),
             (2, 'Valid (Seqcode)', None, None, None)])

        self.assertEqual(designations, {1: M.TYPE_STRAIN_OF_SPECIES, 2: M.TYPE_STRAIN_OF_SPECIES})
        self.assertEqual(sources, {1: 'NCBI;LPSN;Seqcode', 2: 'Seqcode'})

    def test_the_sources_keep_their_order_and_seqcode_is_not_added_twice(self):
        # they went through a set, which ordered them afresh on every run
        _, sources = M.type_designation_changes([(1, 'Valid', None, None, 'LPSN;Seqcode;NCBI')])
        self.assertEqual(sources, {1: 'LPSN;Seqcode;NCBI'})

    def test_a_metagenome_not_used_as_type_is_not_used_as_type_whatever_the_seqcode(self):
        designations, sources = M.type_designation_changes(
            [(1, 'Valid', 'derived from metagenome; not used as type', None, None),
             (2, None, 'derived from metagenome', None, None),
             (3, 'Invalid', None, None, None)])

        self.assertEqual(designations, {1: M.NOT_USED_AS_TYPE})
        self.assertEqual(sources, {1: 'Seqcode'})


if __name__ == '__main__':
    unittest.main()
