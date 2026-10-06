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
on the terminal, which a run under nohup could not answer."""

import gzip
import logging
import os
import shutil
import tempfile
import unittest
from unittest import mock

from gtdb_migration_tk import metadata_database_manager as M
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


class FakeCursor(object):
    """A cursor recording each statement, over the rows of metadata_view."""

    def __init__(self, accessions=(), fail_on=None):
        self.accessions = list(accessions)
        self.fail_on = fail_on
        self.statements = []
        self.result = []

    def execute(self, sql, params=None):
        self.statements.append(sql)
        if self.fail_on and self.fail_on in sql:
            raise RuntimeError('the server refused ' + self.fail_on)
        if 'FROM metadata_view' in sql:
            self.result = [(accession,) for accession in self.accessions]

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

    def test_fields_are_set_to_null_and_written_in_one_commit(self):
        manager = self.manager()
        with mock.patch.object(M, 'GTDBImporter'):
            self.load(manager)

        nulls = [sql for sql in manager.temp_cur.statements if 'NULL' in sql]
        self.assertEqual(sorted(nulls),
                         ['UPDATE metadata_type_material SET gtdb_type_designation_ncbi_taxa = NULL',
                          'UPDATE metadata_type_material SET lpsn_priority_year = NULL'])
        self.assertEqual((manager.temp_con.commits, manager.temp_con.rollbacks), (1, 0))

    def test_do_not_null_field_sets_nothing_to_null(self):
        manager = self.manager()
        with mock.patch.object(M, 'GTDBImporter'):
            self.load(manager, do_not_null_field=True)

        self.assertFalse([sql for sql in manager.temp_cur.statements if 'NULL' in sql])

    def test_a_field_that_cannot_be_written_rolls_back_the_nulls_before_it(self):
        # the NULLs were committed first, leaving the field NULL for every genome
        manager = self.manager()
        with mock.patch.object(M, 'GTDBImporter') as importer:
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
        with self.assertRaises(SystemExit) as raised:
            manager.process_metadata_files(None, table_folder=folder)

        self.assertNotEqual(raised.exception.code, 0)
        self.assertEqual((manager.temp_con.commits, manager.temp_con.rollbacks), (0, 1))


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
    def test_genomes_not_in_the_database_are_skipped_in_one_commit(self):
        # NCBI's organism names and taxonomies cover every assembly NCBI holds
        manager = self.manager(cls=M.NCBITaxDatabaseManager)
        names = self.write('names.tsv', 'GB_GCA_000000002.1\tEscherichia coli\n')
        taxonomy = self.write('tax.tsv', 'GB_GCA_000000002.1\td__Bacteria;p__;c__;o__;f__;g__;s__\n')
        with mock.patch.object(M, 'GTDBImporter') as importer:
            manager.update_ncbitax_db(names, taxonomy, taxonomy, None)

        calls = importer.return_value.import_metadata_to_db.call_args_list
        self.assertEqual([call.args[1] for call in calls],
                         ['ncbi_organism_name', 'ncbi_taxonomy', 'ncbi_taxonomy_unfiltered'])
        self.assertTrue(all(call.kwargs['unknown'] == SKIP for call in calls))
        self.assertEqual(len([sql for sql in manager.temp_cur.statements if 'NULL' in sql]), 3)
        self.assertEqual((manager.temp_con.commits, manager.temp_con.rollbacks), (1, 0))


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
