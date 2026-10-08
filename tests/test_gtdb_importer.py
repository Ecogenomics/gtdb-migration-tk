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

"""Offline unit tests for gtdb_lite/gtdb_importer.py -- how one metadata field
is written, through the database's upsert().

What is tested is what the importer hands the database and what it does when
the database cannot take it: a failure that is printed and passed over leaves a
transaction PostgreSQL then rolls back on commit() without a word, and a genome
the database does not hold fails upsert() for every genome of the field. The
database is stood in for by FakeCursor; nothing is reached.
"""

import logging
import os
import shutil
import tempfile
import unittest

from gtdb_migration_tk.gtdb_lite import gtdb_importer as I


class FakeCursor(object):
    """A cursor over a database holding the given genomes, recording each statement."""

    def __init__(self, genomes=(), fail_on=None):
        self.genomes = list(genomes)
        self.fail_on = fail_on
        self.statements = []
        self.result = []

    def execute(self, sql, params=None):
        self.statements.append((sql, params))
        if self.fail_on and self.fail_on in sql:
            raise RuntimeError('the server refused ' + self.fail_on)
        if sql.startswith('SELECT id_at_source FROM genomes'):
            self.result = [(gid,) for gid in self.genomes]

    def fetchall(self):
        return list(self.result)

    def copy_expert(self, sql, handle):
        self.statements.append((sql, handle.read()))
        if self.fail_on and self.fail_on in sql:
            raise RuntimeError('the server refused ' + self.fail_on)

    def sql(self):
        return [sql for sql, _ in self.statements]

    def upserts(self):
        return [params for sql, params in self.statements if 'upsert(' in sql]


class TempDirCase(unittest.TestCase):
    """A run logged to a directory of its own, where files beside the log go."""

    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='gtdb_importer_test.')
        self.addCleanup(shutil.rmtree, self.dir, True)
        logger = logging.getLogger('timestamp')
        handler = logging.FileHandler(os.path.join(self.dir, 'run.log'))
        logger.addHandler(handler)
        self.addCleanup(logger.removeHandler, handler)
        self.addCleanup(handler.close)


class WritingAField(TempDirCase):
    def test_upsert_is_handed_each_genome_as_the_database_names_it(self):
        cur = FakeCursor(['GCA_000000001.1', 'GCF_000000002.1'])
        I.GTDBImporter(cur).import_metadata_to_db(
            'metadata_genes', 'checkm_completeness', 'FLOAT',
            [('GB_GCA_000000001.1', '99.1'), ('RS_GCF_000000002.1', '87.0')])

        self.assertEqual(cur.upserts(), [('metadata_genes', 'checkm_completeness', 'FLOAT',
                                          ['GCA_000000001.1', 'GCF_000000002.1'],
                                          ['99.1', '87.0'])])

    def test_an_error_writing_is_raised_rather_than_printed(self):
        # printed and passed over, it left a transaction PostgreSQL rolls back on
        # commit() without raising: the command finished as though it had written
        cur = FakeCursor(['GCA_000000001.1'], fail_on='upsert(')
        with self.assertRaises(RuntimeError):
            I.GTDBImporter(cur).import_metadata_to_db(
                'metadata_genes', 'checkm_completeness', 'FLOAT',
                [('GB_GCA_000000001.1', 'not a number')])

    def test_the_genomes_of_the_database_are_read_once_for_every_field(self):
        cur = FakeCursor(['GCA_000000001.1'])
        importer = I.GTDBImporter(cur)
        for field in ('gtdb_domain', 'gtdb_phylum'):
            importer.import_metadata_to_db('metadata_taxonomy', field, 'TEXT',
                                           [('GB_GCA_000000001.1', 'd__Bacteria')])

        reads = [sql for sql, _ in cur.statements if sql.startswith('SELECT id_at_source')]
        self.assertEqual(len(reads), 1)
        self.assertEqual(len(cur.upserts()), 2)

    def test_a_field_with_no_genomes_to_write_calls_no_upsert(self):
        cur = FakeCursor(['GCA_000000001.1'])
        I.GTDBImporter(cur).import_metadata_to_db('metadata_genes', 'checkm_completeness',
                                                  'FLOAT', [])
        self.assertEqual(cur.upserts(), [])


class GenomesTheDatabaseDoesNotHold(TempDirCase):
    DATA = [('GB_GCA_000000001.1', 'x'), ('GB_GCA_000000009.1', 'y'), ('RS_GCF_000000008.1', 'z')]

    def test_a_release_input_naming_one_is_refused_and_nothing_written(self):
        # upsert() would fail the whole field for one such genome
        cur = FakeCursor(['GCA_000000001.1'])
        with self.assertRaises(I.UnknownGenomesError) as raised:
            I.GTDBImporter(cur).import_metadata_to_db('metadata_genes', 'checkm_completeness',
                                                      'FLOAT', self.DATA)

        self.assertEqual(cur.upserts(), [])
        message = str(raised.exception)
        self.assertIn('2 genome(s) to be written to metadata_genes.checkm_completeness', message)
        self.assertIn('GB_GCA_000000009.1', message)

    def test_every_one_refused_is_listed_beside_the_log(self):
        cur = FakeCursor(['GCA_000000001.1'])
        with self.assertRaises(I.UnknownGenomesError):
            I.GTDBImporter(cur).import_metadata_to_db('metadata_genes', 'checkm_completeness',
                                                      'FLOAT', self.DATA)

        path = os.path.join(self.dir, 'unknown_genomes.metadata_genes.checkm_completeness.tsv')
        with open(path) as handle:
            self.assertEqual(handle.read(),
                             'genome_id\nGB_GCA_000000009.1\nRS_GCF_000000008.1\n')

    def test_an_ncbi_wide_input_has_them_skipped_and_the_rest_written(self):
        cur = FakeCursor(['GCA_000000001.1'])
        skipped = I.GTDBImporter(cur).import_metadata_to_db(
            'metadata_ncbi', 'ncbi_organism_name', 'TEXT', self.DATA, unknown=I.SKIP)

        self.assertEqual(skipped, 2)
        self.assertEqual(cur.upserts(), [('metadata_ncbi', 'ncbi_organism_name', 'TEXT',
                                          ['GCA_000000001.1'], ['x'])])

    def test_a_genome_named_without_its_gtdb_prefix_is_written(self):
        # everything before the first '_' was taken off, whatever it was:
        # GCA_000000001.1 became 000000001.1, which no genome is
        cur = FakeCursor(['GCA_000000001.1'])
        I.GTDBImporter(cur).import_metadata_to_db('metadata_genes', 'checkm_completeness',
                                                  'FLOAT', [('GCA_000000001.1', '99.1')])
        self.assertEqual(cur.upserts()[0][3], ['GCA_000000001.1'])

    def test_only_a_leading_gb_or_rs_is_taken_off(self):
        self.assertEqual(I.id_at_source('GB_GCA_000000001.1'), 'GCA_000000001.1')
        self.assertEqual(I.id_at_source('RS_GCF_000000002.1'), 'GCF_000000002.1')
        self.assertEqual(I.id_at_source('GCF_000000002.1'), 'GCF_000000002.1')
        self.assertEqual(I.id_at_source('U_12345'), 'U_12345')
        self.assertIsNone(I.id_at_source(None))

    def test_a_missing_genome_id_is_one_the_database_does_not_hold(self):
        cur = FakeCursor(['GCA_000000001.1'])
        with self.assertRaises(I.UnknownGenomesError):
            I.GTDBImporter(cur).import_metadata_to_db('metadata_taxonomy', 'gtdb_domain',
                                                      'TEXT', [(None, 'd__Bacteria')])


class WritingATablesFieldsTogether(TempDirCase):
    """update_metadata_db wrote each field through its own upsert(), rewriting every row once a field."""

    FIELDS = [('ncbi_taxid', 'INT'), ('ncbi_organism_name', 'TEXT')]

    def test_the_fields_are_copied_in_and_every_row_written_by_one_update(self):
        cur = FakeCursor(['GCA_000000001.1', 'GCF_000000002.1'])
        I.GTDBImporter(cur).import_fields_to_db(
            'metadata_ncbi', self.FIELDS, [('GB_GCA_000000001.1', ['562', 'Escherichia coli']),
                                           ('GCF_000000002.1', [None, 'a\tb\\c'])])

        sql = cur.sql()
        self.assertEqual(sql[1], 'DROP TABLE IF EXISTS gtdb_new_values')
        self.assertEqual(sql[2], 'CREATE TEMPORARY TABLE gtdb_new_values (id_at_source TEXT, "ncbi_taxid" INT, '
                                 '"ncbi_organism_name" TEXT) ON COMMIT DROP')
        copy, copied = cur.statements[3]
        self.assertEqual(copy, 'COPY gtdb_new_values (id_at_source, "ncbi_taxid", "ncbi_organism_name") FROM STDIN')
        # the prefix taken off, a value not given NULL, a tab and a backslash escaped
        self.assertEqual(copied, 'GCA_000000001.1\t562\tEscherichia coli\n'
                                 'GCF_000000002.1\t\\N\ta\\tb\\\\c\n')
        updates = [s for s in sql if s.startswith('UPDATE')]
        self.assertEqual(len(updates), 1)
        self.assertIn('SET "ncbi_taxid" = COALESCE(n."ncbi_taxid", m."ncbi_taxid"), '
                      '"ncbi_organism_name" = COALESCE(n."ncbi_organism_name", m."ncbi_organism_name")', updates[0])
        self.assertIn('LOCK TABLE "metadata_ncbi" IN EXCLUSIVE MODE', sql)
        inserts = [s for s in sql if s.startswith('INSERT')]
        self.assertEqual(len(inserts), 1)
        self.assertIn('WHERE NOT EXISTS (SELECT 1 FROM "metadata_ncbi" m WHERE m.id = g.id)', inserts[0])
        self.assertEqual(cur.upserts(), [])

    def test_a_genome_the_database_does_not_hold_is_refused_before_anything_is_written(self):
        cur = FakeCursor(['GCA_000000001.1'])
        with self.assertRaises(I.UnknownGenomesError) as raised:
            I.GTDBImporter(cur).import_fields_to_db(
                'metadata_ncbi', self.FIELDS, [('GCA_000000001.1', ['562', 'x']), ('GCA_000000099.1', ['1', None])])

        self.assertIn('GCA_000000099.1', str(raised.exception))
        self.assertEqual(cur.sql(), ['SELECT id_at_source FROM genomes'])
        # listed for the field it was given a value of
        self.assertTrue(os.path.exists(os.path.join(self.dir, 'unknown_genomes.metadata_ncbi.ncbi_taxid.tsv')))
        self.assertFalse(os.path.exists(os.path.join(self.dir, 'unknown_genomes.metadata_ncbi.ncbi_organism_name.tsv')))

    def test_skipped_genomes_leave_the_rest_written(self):
        cur = FakeCursor(['GCA_000000001.1'])
        skipped = I.GTDBImporter(cur).import_fields_to_db(
            'metadata_ncbi', self.FIELDS, [('GCA_000000001.1', ['562', 'x']), ('GCA_000000099.1', ['1', 'y'])],
            unknown=I.SKIP)
        self.assertEqual(skipped, 1)
        self.assertEqual(cur.statements[3][1], 'GCA_000000001.1\t562\tx\n')

    def test_no_genome_to_write_writes_nothing(self):
        cur = FakeCursor(['GCA_000000001.1'])
        I.GTDBImporter(cur).import_fields_to_db('metadata_ncbi', self.FIELDS, [])
        self.assertEqual(cur.sql(), ['SELECT id_at_source FROM genomes'])

    def test_a_name_or_type_that_is_not_plain_is_refused(self):
        cur = FakeCursor(['GCA_000000001.1'])
        for table, fields in (('metadata_ncbi; DROP TABLE genomes', self.FIELDS),
                              ('metadata_ncbi', [('ncbi_taxid"', 'INT')]),
                              ('metadata_ncbi', [('ncbi_taxid', 'INT); DROP TABLE genomes; --')])):
            with self.subTest(table=table, fields=fields), self.assertRaises(ValueError):
                I.GTDBImporter(cur).import_fields_to_db(table, fields, [('GCA_000000001.1', ['1'])])


if __name__ == '__main__':
    unittest.main()
