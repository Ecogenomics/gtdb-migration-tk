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

"""Offline unit tests for the database commands of propagate_taxonomy.py:
update_propagated_tax, add_taxonomy_to_database and set_gtdb_domain.

Each writes the taxonomy rank by rank, and committed each rank as it went, so a
run that failed part way left genomes with a domain of the new release and a
species of the old. Each is now one transaction. The database is stood in for
by fakes, and the importer by a mock: what is tested is what is written, and
when it is committed.
"""

import logging
import os
import shutil
import tempfile
import unittest
from unittest import mock

from gtdb_migration_tk import propagate_taxonomy as P
from gtdb_migration_tk.gtdb_lite.gtdb_importer import UnknownGenomesError

TAXONOMY = 'd__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus;s__Bacillus subtilis'


class FakeCursor(object):
    """A cursor answering each statement with the next of the results given it."""

    def __init__(self, fetchone=(), fetchall=()):
        self.to_fetchone = list(fetchone)
        self.to_fetchall = list(fetchall)
        self.statements = []

    def execute(self, sql, params=None):
        self.statements.append(sql)

    def executemany(self, sql, rows):
        self.statements.append(sql)
        self.rows = list(rows)

    def fetchone(self):
        return self.to_fetchone.pop(0)

    def fetchall(self):
        return self.to_fetchall.pop(0)


class FakeConnection(object):
    def __init__(self):
        self.commits = 0
        self.rollbacks = 0

    def commit(self):
        self.commits += 1

    def rollback(self):
        self.rollbacks += 1


class TaxonomyCase(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='propagate_taxonomy_test.')
        self.addCleanup(shutil.rmtree, self.dir, True)
        logging.getLogger('timestamp').addHandler(logging.NullHandler())

    def propagate(self, cursor=None):
        """A Propagate whose database is fakes."""
        propagate = P.Propagate()
        propagate.temp_cur = cursor or FakeCursor()
        propagate.temp_con = FakeConnection()
        return propagate

    def write(self, name, text):
        path = os.path.join(self.dir, name)
        with open(path, 'w') as handle:
            handle.write(text)
        return path

    def metadata(self):
        return self.write('metadata.tsv', 'accession\tformatted_accession\tgtdb_taxonomy\n'
                          'RS_GCF_000000001.1\tG000000001\t{}\n'.format(TAXONOMY))


class AddingTheTaxonomy(TaxonomyCase):
    def test_a_canonical_id_the_metadata_does_not_map_is_refused_before_anything_is_written(self):
        # it became None, which failed the rank's write and was passed over
        taxonomy = self.write('taxonomy.tsv', 'G000000001\t{0}\nG000000099\t{0}\n'.format(TAXONOMY))
        propagate = self.propagate()
        with mock.patch.object(P, 'GTDBImporter') as importer:
            with self.assertRaises(UnknownGenomesError) as raised:
                propagate.add_taxonomy_to_database(taxonomy, self.metadata(), True)

        self.assertIn('G000000099', str(raised.exception))
        importer.return_value.import_metadata_to_db.assert_not_called()
        self.assertEqual((propagate.temp_con.commits, propagate.temp_con.rollbacks), (0, 1))

    def test_every_rank_is_written_for_the_genome_a_canonical_id_names_in_one_commit(self):
        taxonomy = self.write('taxonomy.tsv', 'G000000001\t{}\n'.format(TAXONOMY))
        propagate = self.propagate()
        with mock.patch.object(P, 'GTDBImporter') as importer:
            propagate.add_taxonomy_to_database(taxonomy, self.metadata(), False)

        calls = importer.return_value.import_metadata_to_db.call_args_list
        self.assertEqual([call.args[1] for call in calls],
                         ['gtdb_domain', 'gtdb_phylum', 'gtdb_class', 'gtdb_order',
                          'gtdb_family', 'gtdb_genus', 'gtdb_species'])
        self.assertEqual(calls[0].args[3], [('RS_GCF_000000001.1', 'd__Bacteria')])
        self.assertEqual((propagate.temp_con.commits, propagate.temp_con.rollbacks), (1, 0))


class UpdatingThePropagatedTaxonomy(TaxonomyCase):
    def test_the_truncated_taxonomy_the_ranks_and_the_representatives_are_one_commit(self):
        # the truncated taxonomy was committed first, rank by rank
        taxonomy = self.write('taxonomy.tsv', 'RS_GCF_000000001.1\t{}\n'.format(TAXONOMY))
        reps = self.write('reps.tsv', 'RS_GCF_000000001.1\tTrue\n')
        propagate = self.propagate()
        with mock.patch.object(P, 'GTDBImporter') as importer:
            propagate.add_propagated_taxonomy(taxonomy, self.metadata(), None, True, reps)

        self.assertEqual(len(importer.return_value.import_metadata_to_db.call_args_list), 15)
        self.assertEqual((propagate.temp_con.commits, propagate.temp_con.rollbacks), (1, 0))

    def test_a_rank_that_cannot_be_written_rolls_back_every_rank_before_it(self):
        taxonomy = self.write('taxonomy.tsv', 'RS_GCF_000000001.1\t{}\n'.format(TAXONOMY))
        reps = self.write('reps.tsv', 'RS_GCF_000000001.1\tTrue\n')
        propagate = self.propagate()
        with mock.patch.object(P, 'GTDBImporter') as importer:
            importer.return_value.import_metadata_to_db.side_effect = (
                [None] * 4 + [RuntimeError('refused')])
            with self.assertRaises(RuntimeError):
                propagate.add_propagated_taxonomy(taxonomy, self.metadata(), None, False, reps)

        self.assertEqual((propagate.temp_con.commits, propagate.temp_con.rollbacks), (0, 1))


class SettingTheGTDBDomain(TaxonomyCase):
    def test_an_ncbi_domain_without_its_prefix_ends_the_run_non_zero_with_a_rollback(self):
        # sys.exit() with no code ended it with exit status 0
        cursor = FakeCursor(fetchone=[(120,), (53,)],
                            fetchall=[[(7, 'RS_GCF_000000001.1', 'Bacteria;p__Bacillota')]])
        propagate = self.propagate(cursor)
        with self.assertRaises(SystemExit) as raised:
            propagate.set_gtdb_domain()

        self.assertEqual(raised.exception.code, 1)
        self.assertEqual((propagate.temp_con.commits, propagate.temp_con.rollbacks), (0, 1))


if __name__ == '__main__':
    unittest.main()
