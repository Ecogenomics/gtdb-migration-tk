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

"""Offline unit tests for propagate_taxonomy.py: propagate_gtdb_taxonomy, and the
database commands update_propagated_tax, add_taxonomy_to_database and
set_gtdb_domain.

propagate_gtdb_taxonomy read the new release's genomes from a metadata file
exported by hand, matched a genome to every predecessor it found, counted the
genomes that moved to RefSeq as moved to GenBank, stopped at the first genome
whose taxonomy differed and failed with a KeyError on one the export gave none.

Each writes the taxonomy rank by rank, and committed each rank as it went, so a
run that failed part way left genomes with a domain of the new release and a
species of the old. Each is now one transaction. The database is stood in for
by fakes, and the importer by a mock: what is tested is what is written, and
when it is committed.
"""

import contextlib
import gzip
import io
import logging
import os
import shutil
import tempfile
import unittest
from unittest import mock

from gtdb_migration_tk import __main__ as main_module
from gtdb_migration_tk import main as main_py
from gtdb_migration_tk import propagate_taxonomy as P
from gtdb_migration_tk.gtdb_lite.gtdb_importer import UnknownGenomesError

TAXONOMY = 'd__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus;s__Bacillus subtilis'
ARCHAEON = 'd__Archaea;p__Thermoproteota;c__Thermoprotei;o__Sulfolobales;f__Sulfolobaceae;g__Sulfolobus;s__Sulfolobus acidocaldarius'
OTHER = 'd__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus;s__Bacillus velezensis'
UNSET = (None,) * 7
DOMAIN_ONLY = ('d__Bacteria',) + (None,) * 6
TRUNCATED = ('d__Bacteria', 'p__', 'c__', 'o__', 'f__', 'g__', 's__')


def ranks(taxonomy):
    """The seven ranks of a taxonomy string, as metadata_taxonomy holds them."""
    return tuple(taxonomy.split(';'))


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

    def write_gz(self, name, text):
        path = os.path.join(self.dir, name)
        with gzip.open(path, 'wt') as handle:
            handle.write(text)
        return path

    def propagated(self, taxonomy_rows, reps_rows):
        """A propagate_gtdb_taxonomy output directory."""
        self.write_gz(P.TAXONOMY_NAME, ''.join('{}\t{}\n'.format(*row) for row in taxonomy_rows))
        self.write_gz(P.REPS_NAME, ''.join('{}\t{}\n'.format(*row) for row in reps_rows))
        return self.dir

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
    def update(self, in_dir, cleared=0):
        """Run update_propagated_tax with the importer and reset_unwritten() mocked."""
        propagate = self.propagate()
        with mock.patch.object(P, 'GTDBImporter') as importer, \
                mock.patch.object(P, 'reset_unwritten', return_value=cleared) as reset:
            propagate.add_propagated_taxonomy(in_dir)
        return propagate, importer.return_value.import_metadata_to_db.call_args_list, reset.call_args_list

    def test_the_ranks_the_ranks_below_the_domain_cleared_and_the_representatives_are_one_commit(self):
        # the truncated taxonomy was committed first, rank by rank, from a hand export (-m)
        in_dir = self.propagated([('RS_GCF_000000001.1', TAXONOMY)],
                                 [('RS_GCF_000000001.1', 'True'), ('GB_GCA_000000002.1', 'False')])
        propagate, calls, resets = self.update(in_dir)

        self.assertEqual([call.args[1] for call in calls],
                         ['gtdb_domain', 'gtdb_phylum', 'gtdb_class', 'gtdb_order', 'gtdb_family', 'gtdb_genus',
                          'gtdb_species', 'gtdb_representative'])
        self.assertEqual(calls[6].args[3], [('RS_GCF_000000001.1', 's__Bacillus subtilis')])
        self.assertEqual(calls[7].args[2:4], ('BOOLEAN', [('RS_GCF_000000001.1', 'True'),
                                                          ('GB_GCA_000000002.1', 'False')]))
        self.assertEqual((propagate.temp_con.commits, propagate.temp_con.rollbacks), (1, 0))

    def test_a_genome_the_taxonomy_does_not_name_keeps_its_domain_and_loses_the_ranks_below_it(self):
        # --truncate_taxonomy did that for every genome of a hand export, before writing
        in_dir = self.propagated([('RS_GCF_000000001.1', TAXONOMY)], [('RS_GCF_000000001.1', 'True')])
        with self.assertLogs('timestamp', level='WARNING') as logged:
            _, _, resets = self.update(in_dir, cleared=3)

        self.assertEqual([call.args[1:] for call in resets],
                         [('metadata_taxonomy', 'gtdb_' + rank, ['GCF_000000001.1'])
                          for rank in ('phylum', 'class', 'order', 'family', 'genus', 'species')])
        self.assertEqual(len(logged.records), 6)
        self.assertIn('Set gtdb_phylum to NULL for 3 genome(s) the propagated taxonomy does not name, in ',
                      logged.records[0].getMessage())

    def test_each_step_has_a_line_in_the_log_with_how_long_it_took(self):
        # an 11 minute run logged nothing between reading the file and the end of the taxonomy
        in_dir = self.propagated([('RS_GCF_000000001.1', TAXONOMY)], [('RS_GCF_000000001.1', 'True')])
        with self.assertLogs('timestamp', level='INFO') as logged:
            self.update(in_dir)

        steps = [r.getMessage() for r in logged.records if r.getMessage().endswith(' s.')]
        self.assertEqual(len(steps), 14)
        self.assertTrue(steps[0].startswith('Wrote gtdb_domain for 1 genomes in '))
        self.assertTrue(steps[1].startswith('Wrote gtdb_phylum for 1 genomes in '))
        self.assertTrue(steps[2].startswith('Set gtdb_phylum to NULL for 0 genome(s)'))
        self.assertTrue(steps[-1].startswith('Wrote whether each of 1 genomes is a representative, 1 of them, in '))

    def test_a_rank_that_cannot_be_written_rolls_back_every_rank_before_it(self):
        in_dir = self.propagated([('RS_GCF_000000001.1', TAXONOMY)], [('RS_GCF_000000001.1', 'True')])
        propagate = self.propagate()
        with mock.patch.object(P, 'GTDBImporter') as importer, mock.patch.object(P, 'reset_unwritten', return_value=0):
            importer.return_value.import_metadata_to_db.side_effect = [None] * 4 + [RuntimeError('refused')]
            with self.assertRaises(RuntimeError):
                propagate.add_propagated_taxonomy(in_dir)
        self.assertEqual((propagate.temp_con.commits, propagate.temp_con.rollbacks), (0, 1))


class DomainCursor(FakeCursor):
    """A database of genomes without a domain and an aligned_markers of a few blocks.

    genomes are (id, GTDB name, id_at_source, NCBI taxonomy); ranges map the first
    block of each range of aligned_markers to its rows (genome id, bac120 found,
    ar53 found).
    """

    def __init__(self, genomes, ranges, blocks):
        super().__init__()
        self.genomes, self.ranges, self.blocks = genomes, ranges, blocks
        self.ranges_read = []

    def execute(self, sql, params=None):
        self.statements.append(sql)
        if sql.startswith('SELECT set_id, count(*) FROM marker_set_contents'):
            self.result = [(1, 120), (19, 53)]
        elif sql == P.MISSING_DOMAIN:
            self.result = list(self.genomes)
        elif sql == P.ALIGNED_MARKERS_BLOCKS:
            self.result = [(self.blocks, 8192)]
        elif sql == P.MARKERS_FOUND:
            low, high = params[-2:]
            self.ranges_read.append((low, high))
            self.result = self.ranges.get(int(low[1:].split(',')[0]), [])

    def fetchone(self):
        return self.result[0]

    def fetchall(self):
        return list(self.result)


class SettingTheGTDBDomain(TaxonomyCase):
    def set_domain(self, genomes, ranges=None, blocks=1):
        """set_gtdb_domain over genomes (id, GTDB name, id_at_source, NCBI taxonomy).

        @return: (propagate, importer mock, log messages, what reached stdout or stderr).
        """
        cursor = DomainCursor(genomes, ranges or {}, blocks)
        propagate = self.propagate(cursor)
        console = io.StringIO()
        with mock.patch.object(P, 'GTDBImporter') as importer, contextlib.redirect_stdout(console), \
                contextlib.redirect_stderr(console), \
                self.assertLogs('timestamp', level='INFO') as logged:
            try:
                propagate.set_gtdb_domain(self.dir)
            except SystemExit as exc:
                propagate.exit_code = exc.code
        return propagate, importer, [r.getMessage() for r in logged.records], console.getvalue()

    def table(self, name):
        with open(os.path.join(self.dir, name)) as handle:
            return [line.split('\t') for line in handle.read().splitlines()]

    def test_the_markers_decide_the_domain_and_ncbi_does_where_they_are_too_few(self):
        _, importer, _, _ = self.set_domain([
            (1, 'RS_GCF_000000001.1', 'GCF_000000001.1', TAXONOMY),
            (2, 'GB_GCA_000000002.1', 'GCA_000000002.1', ARCHAEON),   # 10% and 11.3%: the larger share
            (3, 'GB_GCA_000000003.1', 'GCA_000000003.1', TAXONOMY),   # the markers over NCBI
            (4, 'GB_GCA_000000004.1', 'GCA_000000004.1', ARCHAEON),   # 9.2% and 9.4%: NCBI's
            (5, 'GB_GCA_000000005.1', 'GCA_000000005.1', TAXONOMY)],  # not aligned: NCBI's
            {0: [(1, 118, 4), (2, 12, 6), (3, 0, 50), (4, 11, 5)]})
        importer.return_value.import_metadata_to_db.assert_called_once_with(
            'metadata_taxonomy', 'gtdb_domain', 'TEXT',
            [('GCF_000000001.1', 'd__Bacteria'), ('GCA_000000002.1', 'd__Archaea'),
             ('GCA_000000003.1', 'd__Archaea'), ('GCA_000000004.1', 'd__Archaea'),
             ('GCA_000000005.1', 'd__Bacteria')])

    def test_a_disagreement_is_listed_in_the_output_directory_and_counted_in_one_line_and_nothing_is_printed(self):
        # each was a WARNING of its own, and each genome the markers could not decide a print()
        propagate, _, messages, stdout = self.set_domain([
            (1, 'RS_GCF_000000001.1', 'GCF_000000001.1', TAXONOMY),
            (3, 'GB_GCA_000000003.1', 'GCA_000000003.1', TAXONOMY),
            (5, 'GB_GCA_000000005.1', 'GCA_000000005.1', TAXONOMY)],
            {0: [(1, 120, 0), (3, 0, 50)]})

        self.assertEqual(stdout, '')
        self.assertEqual(propagate.temp_con.commits, 1)
        self.assertEqual(self.table(P.DOMAIN_DISAGREEMENTS_NAME), [
            list(P.DOMAIN_DISAGREEMENTS_HEADER),
            ['GB_GCA_000000003.1', 'd__Bacteria', 'd__Archaea', '0.00', '94.34']])
        self.assertEqual(self.table(P.DOMAIN_FROM_NCBI_NAME), [
            list(P.DOMAIN_FROM_NCBI_HEADER), ['GB_GCA_000000005.1', 'd__Bacteria', '0.00', '0.00']])
        about_disagreements = [m for m in messages if 'not their NCBI domain' in m]
        self.assertEqual(about_disagreements, ['1 genome(s) were given a GTDB domain by their markers that is not '
                                               'their NCBI domain; each is listed in {}.'.format(
                                                   os.path.join(self.dir, P.DOMAIN_DISAGREEMENTS_NAME))])
        self.assertFalse(any('GCA_000000003.1' in m for m in messages))

    def test_with_no_disagreement_the_list_is_its_header_and_the_line_says_none(self):
        _, _, messages, _ = self.set_domain([(1, 'RS_GCF_000000001.1', 'GCF_000000001.1', TAXONOMY)],
                                            {0: [(1, 120, 0)]})
        self.assertEqual(self.table(P.DOMAIN_DISAGREEMENTS_NAME), [list(P.DOMAIN_DISAGREEMENTS_HEADER)])
        self.assertTrue(any(m.startswith('0 genome(s) were given a GTDB domain') for m in messages))

    def test_aligned_markers_is_read_once_in_ranges_of_its_blocks_and_a_genome_across_two_is_summed(self):
        # two queries a genome, 586,000 over r237, and then one query of which nothing said how far it was
        propagate, importer, messages, _ = self.set_domain(
            [(1, 'RS_GCF_000000001.1', 'GCF_000000001.1', ARCHAEON)],
            {0: [(1, 3, 2), (9, 120, 0)], 250: [(1, 3, 2)], 500: [(1, 0, 2)]}, blocks=1000)

        cursor = propagate.temp_cur
        self.assertEqual(cursor.statements[0], 'SET TRANSACTION ISOLATION LEVEL REPEATABLE READ')
        self.assertEqual(len(cursor.ranges_read), P.ALIGNED_MARKER_PARTS)
        bounds = [int(t[1:].split(',')[0]) for pair in cursor.ranges_read for t in pair]
        self.assertEqual((bounds[0], bounds[-1]), (0, 1000))
        self.assertTrue(all(bounds[i] == bounds[i + 1] for i in range(1, len(bounds) - 1, 2)))
        # 6 of 120 and 6 of 53, read in three ranges; genome 9 is not one without a domain
        importer.return_value.import_metadata_to_db.assert_called_once_with(
            'metadata_taxonomy', 'gtdb_domain', 'TEXT', [('GCF_000000001.1', 'd__Archaea')])
        self.assertEqual([m for m in messages if m.startswith('Read ') and '% of aligned_markers' in m],
                         ['Read {}% of aligned_markers in 0.0 min.'.format(p) for p in range(10, 100, 10)])

    def test_a_table_smaller_than_the_parts_is_read_a_block_a_part(self):
        propagate, _, _, _ = self.set_domain([(1, 'RS_GCF_000000001.1', 'GCF_000000001.1', TAXONOMY)],
                                             {0: [(1, 120, 0)]}, blocks=3)
        self.assertEqual(propagate.temp_cur.ranges_read, [('(0,0)', '(1,0)'), ('(1,0)', '(2,0)'), ('(2,0)', '(3,0)')])

    def test_with_no_genome_without_a_domain_aligned_markers_is_not_read(self):
        propagate, importer, _, _ = self.set_domain([])
        self.assertEqual(propagate.temp_cur.ranges_read, [])
        importer.return_value.import_metadata_to_db.assert_not_called()
        self.assertEqual(self.table(P.DOMAIN_DISAGREEMENTS_NAME), [list(P.DOMAIN_DISAGREEMENTS_HEADER)])

    def test_an_ncbi_domain_without_its_prefix_ends_the_run_non_zero_before_reading_and_each_is_listed(self):
        # sys.exit() with no code ended it with exit status 0, at the first
        propagate, importer, messages, stdout = self.set_domain([
            (1, 'RS_GCF_000000001.1', 'GCF_000000001.1', 'Bacteria;p__Bacillota'),
            (2, 'RS_GCF_000000002.1', 'GCF_000000002.1', TAXONOMY)])

        self.assertEqual(propagate.exit_code, 1)
        self.assertEqual((propagate.temp_con.commits, propagate.temp_con.rollbacks), (0, 1))
        self.assertEqual(propagate.temp_cur.ranges_read, [])
        importer.return_value.import_metadata_to_db.assert_not_called()
        self.assertEqual(self.table(P.NCBI_DOMAIN_ERRORS_NAME), [
            list(P.NCBI_DOMAIN_ERRORS_HEADER), ['RS_GCF_000000001.1', 'Bacteria']])
        self.assertTrue(any(m.startswith('1 genome(s) have an NCBI domain without its d__ prefix') for m in messages))
        self.assertEqual(stdout, '')


class PropagatingCase(TaxonomyCase):
    def setUp(self):
        super().setUp()
        self.out_dir = os.path.join(self.dir, 'out')
        os.makedirs(self.out_dir)
        self.log_dir = os.path.join(self.dir, 'logs')
        os.makedirs(self.log_dir)
        patcher = mock.patch.object(P, 'log_directory', return_value=self.log_dir)
        patcher.start()
        self.addCleanup(patcher.stop)

    def previous(self, name, rows, compress=False, header=('accession', 'ncbi_genbank_assembly_accession',
                                                           'gtdb_representative', 'gtdb_taxonomy')):
        """A GTDB metadata file: rows are (accession, representative 't'/'f', taxonomy)."""
        lines = ['\t'.join(header)]
        for accession, rep, taxonomy in rows:
            values = {'accession': accession, 'ncbi_genbank_assembly_accession': 'x',
                      'gtdb_representative': rep, 'gtdb_taxonomy': taxonomy}
            lines.append('\t'.join(values[column] for column in header))
        text = '\n'.join(lines) + '\n'
        return self.write_gz(name, text) if compress else self.write(name, text)

    def run_propagation(self, previous_files, database_rows):
        """Propagate against a database of rows (accession, seven ranks)."""
        propagate = self.propagate(FakeCursor(fetchall=[[(row[0],) + tuple(row[1]) for row in database_rows]]))
        with self.assertLogs('timestamp', level='INFO') as logged:
            propagate.propagate_taxonomy(previous_files, self.out_dir)
        return self.read_out(P.TAXONOMY_NAME), self.read_out(P.REPS_NAME), [r.getMessage() for r in logged.records]

    def read_out(self, name):
        path = os.path.join(self.out_dir, name)
        with open(path, 'rb') as handle:
            self.assertEqual(handle.read(2), b'\x1f\x8b')
        with gzip.open(path, 'rt') as handle:
            return [line.split('\t') for line in handle.read().splitlines()]


class PropagatingTheTaxonomy(PropagatingCase):
    def test_every_metadata_file_is_read_and_the_two_tables_are_written_gzipped_under_their_names(self):
        # a release's metadata is an archaeal and a bacterial file; it took one
        archaea = self.previous('ar53_metadata_r232.tsv.gz', [('GB_GCA_000000001.1', 't', ARCHAEON)],
                                compress=True)
        bacteria = self.previous('bac120_metadata_r232.tsv', [('RS_GCF_000000002.1', 'f', TAXONOMY)],
                                 header=('gtdb_taxonomy', 'gtdb_representative', 'accession'))
        taxonomy, reps, _ = self.run_propagation(
            [archaea, bacteria],
            [('RS_GCF_000000002.1', ranks(TAXONOMY)), ('GB_GCA_000000001.1', ranks(ARCHAEON)),
             ('GB_GCA_000000003.1', UNSET)])

        self.assertEqual(taxonomy, [['GB_GCA_000000001.1', ARCHAEON], ['RS_GCF_000000002.1', TAXONOMY]])
        self.assertEqual(reps, [['GB_GCA_000000001.1', 'True'], ['GB_GCA_000000003.1', 'False'],
                                ['RS_GCF_000000002.1', 'False']])
        self.assertEqual(sorted(os.listdir(self.out_dir)), sorted([P.TAXONOMY_NAME, P.REPS_NAME]))

    def test_a_new_version_or_database_inherits_from_its_one_predecessor_and_each_move_is_counted_as_itself(self):
        # the moves were counted the other way round, and a genome took every predecessor it found
        previous = self.previous('bac120.tsv', [('GB_GCA_000000001.1', 't', TAXONOMY),
                                               ('RS_GCF_000000002.1', 'f', OTHER),
                                               ('GB_GCA_000000003.2', 'f', ARCHAEON),
                                               ('GB_GCA_000000004.1', 't', TAXONOMY)])
        taxonomy, reps, messages = self.run_propagation(
            [previous], [('RS_GCF_000000001.1', UNSET),      # GenBank to RefSeq
                         ('GB_GCA_000000002.1', DOMAIN_ONLY),  # RefSeq to GenBank
                         ('GB_GCA_000000003.1', UNSET),      # an earlier version
                         ('RS_GCF_000000004.3', UNSET)])     # to RefSeq, and a new version

        self.assertEqual(taxonomy, [['GB_GCA_000000002.1', OTHER], ['GB_GCA_000000003.1', ARCHAEON],
                                    ['RS_GCF_000000001.1', TAXONOMY], ['RS_GCF_000000004.3', TAXONOMY]])
        self.assertEqual(reps, [['GB_GCA_000000002.1', 'False'], ['GB_GCA_000000003.1', 'False'],
                                ['RS_GCF_000000001.1', 'True'], ['RS_GCF_000000004.3', 'True']])
        self.assertIn('2 genomes of the previous release moved from GenBank to RefSeq.', messages)
        self.assertIn('1 genomes of the previous release moved from RefSeq to GenBank.', messages)
        self.assertIn('1 genomes of the previous release have a new version.', messages)

    def test_a_genome_not_in_the_database_is_counted_with_the_representatives_among_them(self):
        previous = self.previous('bac120.tsv', [('GB_GCA_000000001.1', 't', TAXONOMY),
                                               ('GB_GCA_000000002.1', 'f', OTHER)])
        taxonomy, reps, messages = self.run_propagation([previous], [('GB_GCA_000000002.1', UNSET)])

        self.assertEqual(taxonomy, [['GB_GCA_000000002.1', OTHER]])
        self.assertIn('1 genomes of the previous release are not in the database, 1 of them representatives, '
                      'e.g.: GB_GCA_000000001.1.', messages)

    def test_a_user_genome_inherits_only_from_itself(self):
        previous = self.previous('bac120.tsv', [('U_10', 'f', TAXONOMY)])
        taxonomy, _, _ = self.run_propagation([previous], [('U_10', UNSET), ('U_11', UNSET)])
        self.assertEqual(taxonomy, [['U_10', TAXONOMY]])


class TheDatabasesTaxonomy(PropagatingCase):
    def test_a_taxonomy_not_yet_set_is_not_a_mismatch(self):
        # one the export gave none for failed the run with a KeyError
        previous = self.previous('bac120.tsv', [('GB_GCA_000000001.1', 'f', TAXONOMY),
                                               ('GB_GCA_000000002.1', 'f', OTHER),
                                               ('GB_GCA_000000003.1', 'f', OTHER + ';'),
                                               ('GB_GCA_000000004.1', 'f', TAXONOMY)])
        # the bare prefixes are what update_propagated_tax --truncate_taxonomy writes
        taxonomy, _, _ = self.run_propagation(
            [previous], [('GB_GCA_000000001.1', UNSET), ('GB_GCA_000000002.1', DOMAIN_ONLY),
                         ('GB_GCA_000000003.1', ranks(OTHER)), ('GB_GCA_000000004.1', TRUNCATED)])
        self.assertEqual([row[0] for row in taxonomy], ['GB_GCA_000000001.1', 'GB_GCA_000000002.1',
                                                        'GB_GCA_000000003.1', 'GB_GCA_000000004.1'])

    def test_every_genome_whose_taxonomy_is_not_the_previous_releases_is_listed_and_nothing_written(self):
        # the run stopped at the first, printing it
        previous = self.previous('bac120.tsv', [('GB_GCA_000000001.1', 'f', TAXONOMY),
                                               ('GB_GCA_000000002.1', 'f', TAXONOMY),
                                               ('GB_GCA_000000003.1', 'f', TAXONOMY)])
        propagate = self.propagate(FakeCursor(fetchall=[[
            ('GB_GCA_000000001.1',) + ranks(OTHER), ('GB_GCA_000000002.1',) + ranks(TAXONOMY),
            ('GB_GCA_000000003.1',) + ranks(ARCHAEON)]]))

        with self.assertRaisesRegex(P.PropagationError, '2 genome.*GB_GCA_000000001.1, GB_GCA_000000003.1'):
            propagate.propagate_taxonomy([previous], self.out_dir)

        self.assertEqual(os.listdir(self.out_dir), [])
        with open(os.path.join(self.log_dir, P.MISMATCHES_NAME)) as handle:
            self.assertEqual(handle.read().splitlines(), [
                'genome_id\tprevious_release\tdatabase',
                'GB_GCA_000000001.1\t{}\t{}'.format(TAXONOMY, OTHER),
                'GB_GCA_000000003.1\t{}\t{}'.format(TAXONOMY, ARCHAEON)])

    def test_the_genomes_are_read_from_the_database_with_their_seven_ranks(self):
        cursor = FakeCursor(fetchall=[[('GB_GCA_000000001.1',) + UNSET]])
        propagate = self.propagate(cursor)
        propagate.current_release()
        self.assertEqual(cursor.statements, [P.CURRENT_GENOMES])
        for column in ('gtdb_domain', 'gtdb_species', 'external_id_prefix', 'LEFT JOIN metadata_taxonomy'):
            self.assertIn(column, P.CURRENT_GENOMES)


class OnePredecessorEach(PropagatingCase):
    def test_a_genome_named_twice_by_the_previous_release_is_refused(self):
        archaea = self.previous('ar53.tsv', [('GB_GCA_000000001.1', 'f', ARCHAEON)])
        bacteria = self.previous('bac120.tsv', [('GB_GCA_000000001.1', 'f', TAXONOMY)])
        with self.assertRaisesRegex(P.PropagationError, 'named more than once.*GB_GCA_000000001.1'):
            self.propagate().propagate_taxonomy([archaea, bacteria], self.out_dir)

    def test_a_canonical_accession_held_twice_by_either_side_is_refused(self):
        twice = self.previous('bac120.tsv', [('GB_GCA_000000001.1', 'f', TAXONOMY),
                                            ('RS_GCF_000000001.2', 'f', TAXONOMY)])
        with self.assertRaisesRegex(P.PropagationError, 'G000000001: GB_GCA_000000001.1, RS_GCF_000000001.2'):
            self.propagate().propagate_taxonomy([twice], self.out_dir)

        once = self.previous('ar53.tsv', [('GB_GCA_000000001.1', 'f', TAXONOMY)])
        propagate = self.propagate(FakeCursor(fetchall=[[('GB_GCA_000000001.2',) + UNSET,
                                                         ('RS_GCF_000000001.2',) + UNSET]]))
        with self.assertRaisesRegex(P.PropagationError, 'database.*G000000001'):
            propagate.propagate_taxonomy([once], self.out_dir)
        self.assertEqual(os.listdir(self.out_dir), [])


class TheCommandLine(TaxonomyCase):
    def test_propagate_gtdb_taxonomy_names_its_database_takes_several_previous_files_and_an_output_dir(self):
        archaea, bacteria = self.write('ar53.tsv', ''), self.write('bac120.tsv', '')
        out_dir = os.path.join(self.dir, 'propagated')
        options = main_module.get_main_parser().parse_args(
            ['propagate_gtdb_taxonomy', '--db_service', 'gtdb_r237', '--gtdb_metadata_prev', archaea, bacteria,
             '-o', out_dir, '-l', os.path.join(self.dir, 'run.log')])
        with mock.patch.object(main_py, 'Propagate') as propagate:
            main_py.OptionsParser().parse_options(options)
        propagate.assert_called_once_with({'service': 'gtdb_r237'})
        propagate.return_value.propagate_taxonomy.assert_called_once_with([archaea, bacteria], out_dir)
        self.assertTrue(os.path.isdir(out_dir))

    def test_a_propagation_error_ends_the_run_exiting_1(self):
        options = main_module.get_main_parser().parse_args(
            ['propagate_gtdb_taxonomy', '--db_service', 'gtdb_r237', '--gtdb_metadata_prev',
             self.write('bac120.tsv', ''), '-o', self.dir, '-l', os.path.join(self.dir, 'run.log')])
        with mock.patch.object(main_py, 'Propagate') as propagate, \
                self.assertLogs('timestamp', level='ERROR'), self.assertRaises(SystemExit) as ended:
            propagate.return_value.propagate_taxonomy.side_effect = P.PropagationError('no')
            main_py.OptionsParser().parse_options(options)
        self.assertEqual(ended.exception.code, 1)

    def test_update_propagated_tax_reads_the_directory_it_is_given_as_i(self):
        options = main_module.get_main_parser().parse_args(
            ['update_propagated_tax', '--db_service', 'gtdb_r237', '-i', self.dir,
             '-l', os.path.join(self.dir, 'run.log')])
        with mock.patch.object(main_py, 'Propagate') as propagate:
            main_py.OptionsParser().parse_options(options)
        propagate.return_value.add_propagated_taxonomy.assert_called_once_with(self.dir)

    def test_set_gtdb_domain_requires_an_output_directory_and_is_handed_it(self):
        out_dir = os.path.join(self.dir, 'domain')
        argv = ['set_gtdb_domain', '--db_service', 'gtdb_r237', '-l', os.path.join(self.dir, 'run.log')]
        with mock.patch('sys.stderr'), self.assertRaises(SystemExit) as ended:
            main_module.get_main_parser().parse_args(argv)
        self.assertEqual(ended.exception.code, 2)

        options = main_module.get_main_parser().parse_args(argv + ['-o', out_dir])
        with mock.patch.object(main_py, 'Propagate') as propagate:
            main_py.OptionsParser().parse_options(options)
        propagate.return_value.set_gtdb_domain.assert_called_once_with(out_dir)
        self.assertTrue(os.path.isdir(out_dir))

    def test_the_options_the_files_replaced_are_no_longer_accepted(self):
        for argv in (['propagate_gtdb_taxonomy', '--db_service', 'x', '--gtdb_metadata_prev', 'p.tsv',
                      '--gtdb_metadata_cur', 'c.tsv', '-o', self.dir, '-l', 'run.log'],
                     ['propagate_gtdb_taxonomy', '--db_service', 'x', '--gtdb_metadata_prev', 'p.tsv',
                      '-t', 't.tsv', '--rep_file', 'r.tsv', '-l', 'run.log'],
                     ['update_propagated_tax', '--db_service', 'x', '-t', 't.tsv', '--rep_file', 'r.tsv',
                      '-l', 'run.log'],
                     ['update_propagated_tax', '--db_service', 'x', '-i', self.dir, '-m', 'm.tsv',
                      '-l', 'run.log'],
                     ['update_propagated_tax', '--db_service', 'x', '-i', self.dir, '--genome_list', 'g.tsv',
                      '-l', 'run.log'],
                     ['update_propagated_tax', '--db_service', 'x', '-i', self.dir, '--truncate_taxonomy',
                      '-l', 'run.log']):
            with self.subTest(argv=argv[0]), mock.patch('sys.stderr'), self.assertRaises(SystemExit) as ended:
                main_module.get_main_parser().parse_args(argv)
            self.assertEqual(ended.exception.code, 2)


if __name__ == '__main__':
    unittest.main()
