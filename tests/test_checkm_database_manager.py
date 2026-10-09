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

"""Offline unit tests for checkm_database_manager.py -- update_checkm_db, which
writes the CheckM estimates of a release to metadata_genes.

What is tested is what the command decides to write from checkm's release
files (plan_checkm_import(), which reads no database) and what it hands the
importer. The command wrote only the genomes of a table exported from the
database (--metadata), printing every other one as skipped, so a genome the
database did not hold, or an export of the wrong release, lost its estimates in
a run that succeeded. And a genome checkm left out kept the estimates of the
sequences it had before the update, as though they were its own. And a genome
checkm never planned had no estimates, and nothing said so. The database is
stood in for by a cursor over the genomes it holds; nothing is reached.
"""

import gzip
import logging
import os
import shutil
import tempfile
import unittest
from unittest import mock

from gtdb_migration_tk import __main__ as main_module
from gtdb_migration_tk import checkm_database_manager as C
from gtdb_migration_tk import main as main_py
from gtdb_migration_tk.gtdb_lite.gtdb_importer import UnknownGenomesError

# the columns of qa joined with tree_qa -o 2, as checkm writes checkm.profiles.tsv
PROFILE_HEADER = ['Bin Id', 'Marker lineage', '# genomes', '# markers', '# marker sets',
                  '0', '1', '2', '3', '4', '5+', 'Completeness', 'Contamination',
                  'Strain heterogeneity', '# unique markers (of 43)', '# multi-copy',
                  'Taxonomy (contained)', 'Taxonomy (sister lineage)', 'GC', 'Genome size (Mbp)']
QA_HEADER = ['Bin Id', 'Marker lineage', '# genomes', '# markers', '# marker sets',
             '0', '1', '2', '3', '4', '5+', 'Completeness', 'Contamination', 'Strain heterogeneity']

GENOMES = ('GCA_000000001.1', 'GCF_000000002.1')


def profile_row(accession, completeness='98.5', contamination='1.2', sh='50.0',
                lineage='k__Bacteria (UID203)'):
    return [accession + '_protein', lineage, '5449', '104', '58', '0', '104', '0', '0', '0', '0',
            completeness, contamination, sh, '43', '0', 'k__Bacteria', 'p__Firmicutes', '50.1', '3.2']


def qa_row(accession, sh='0.0'):
    return [accession + '_protein', 'k__Bacteria (UID203)', '5449', '104', '58',
            '0', '104', '0', '0', '0', '0', '98.5', '1.2', sh]


# the estimates a genome holds in metadata_genes, one value for each of C.CHECKM_FIELDS
HELD = (99.47, 0.5, 0.0, 'o__Sphingomonadales (UID3310)', 77, 474, 300, 0.0)


class FakeCursor(object):
    """A cursor over a database holding the given genomes, recording each statement.

    estimates maps a genome to (has_changed, the CheckM estimates it holds);
    missing is (genome, has_changed, added on its last update) of each NCBI
    genome with no completeness once the update is made.
    """

    def __init__(self, genomes=GENOMES, estimates=None, missing=()):
        self.genomes = list(genomes)
        self.estimates = dict(estimates or {})
        self.missing = sorted(missing)
        self.statements = []
        self.result = []

    def execute(self, sql, params=None):
        self.statements.append((sql, params))
        if sql.startswith('SELECT id_at_source FROM genomes'):
            self.result = [(gid,) for gid in self.genomes]
        elif sql.startswith('SELECT g.id, g.id_at_source, g.has_changed'):
            self.result = [(self.genomes.index(gid), gid, has_changed) + tuple(values)
                           for gid, (has_changed, values) in sorted(self.estimates.items())
                           if gid in params[0] and any(v is not None for v in values)]
        elif sql.startswith('SELECT g.id_at_source, g.has_changed'):
            self.result = list(self.missing)

    def fetchall(self):
        return list(self.result)

    def upserts(self):
        return {params[1]: (params[2], sorted(zip(params[3], params[4])))
                for sql, params in self.statements if 'upsert(' in sql}

    def cleared(self):
        """The ids of the rows of metadata_genes set to NULL."""

        return [params[0] for sql, params in self.statements if sql.startswith('UPDATE metadata_genes')]


class FakeConnection(object):
    def __init__(self):
        self.commits = 0
        self.rollbacks = 0

    def commit(self):
        self.commits += 1

    def rollback(self):
        self.rollbacks += 1


class TempDirCase(unittest.TestCase):
    """A release's CheckM files, written to a directory of the test's own."""

    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='checkm_database_manager_test.')
        self.addCleanup(shutil.rmtree, self.dir, True)
        # genomes the importer refuses are listed beside the log
        logger = logging.getLogger('timestamp')
        handler = logging.FileHandler(os.path.join(self.dir, 'run.log'))
        logger.addHandler(handler)
        self.addCleanup(logger.removeHandler, handler)
        self.addCleanup(handler.close)
        # what the run warns of, kept beside the log file rather than in place of it
        self.warnings = []
        warnings = logging.Handler(level=logging.WARNING)
        warnings.emit = lambda record: self.warnings.append(record.getMessage())
        logger.addHandler(warnings)
        self.addCleanup(logger.removeHandler, warnings)

    def table(self, name, header, rows, compress=False):
        path = os.path.join(self.dir, name)
        text = ''.join('\t'.join(row) + '\n' for row in [header] + list(rows))
        with (gzip.open(path, 'wt') if compress else open(path, 'w')) as handle:
            handle.write(text)
        return path

    def release_files(self, profile_rows=None, qa_rows=None, not_assessed=()):
        profile = self.table('checkm.profiles.tsv', PROFILE_HEADER,
                             profile_rows if profile_rows is not None
                             else [profile_row(GENOMES[0]), profile_row(GENOMES[1], '87.0', '3.4')])
        qa = self.table('checkm.qa_sh100.tsv', QA_HEADER,
                        qa_rows if qa_rows is not None
                        else [qa_row(GENOMES[0], '12.5'), qa_row(GENOMES[1], '0.0')])
        left_out = self.table('checkm_not_assessed.tsv', ['genome_id', 'reason', 'detail'],
                              [(gid, 'genome_too_large', '2067143420') for gid in not_assessed])
        return profile, qa, left_out

    def manager(self, cursor):
        """A manager built without a database, its cursor and connection fakes."""

        manager = C.CheckMDatabaseManager.__new__(C.CheckMDatabaseManager)
        manager.logger = logging.getLogger('timestamp')
        manager.temp_cur = cursor
        manager.temp_con = FakeConnection()
        return manager


class PlanningTheImport(TempDirCase):
    def plan(self, *args, **kwargs):
        fields = C.plan_checkm_import(*self.release_files(*args, **kwargs)).fields
        return {field: (data_type, rows) for field, data_type, rows in fields}

    def test_every_field_is_written_for_every_genome_named_by_its_accession(self):
        genomes, fields, _not_assessed = C.plan_checkm_import(*self.release_files())

        self.assertEqual(genomes, list(GENOMES))
        self.assertEqual([(field, data_type) for field, data_type, _rows in fields],
                         [('checkm_completeness', 'FLOAT'),
                          ('checkm_contamination', 'FLOAT'),
                          ('checkm_strain_heterogeneity', 'FLOAT'),
                          ('checkm_marker_lineage', 'TEXT'),
                          ('checkm_genome_count', 'INT'),
                          ('checkm_marker_count', 'INT'),
                          ('checkm_marker_set_count', 'INT'),
                          ('checkm_strain_heterogeneity_100', 'FLOAT')])
        plan = {field: rows for field, _type, rows in fields}
        self.assertEqual(plan['checkm_completeness'],
                         [('GCA_000000001.1', '98.5'), ('GCF_000000002.1', '87.0')])
        self.assertEqual(plan['checkm_marker_lineage'],
                         [('GCA_000000001.1', 'k__Bacteria (UID203)'),
                          ('GCF_000000002.1', 'k__Bacteria (UID203)')])

    def test_strain_heterogeneity_at_100_percent_is_taken_from_the_sh100_table_not_the_profile(self):
        plan = self.plan()

        self.assertEqual(plan['checkm_strain_heterogeneity'][1],
                         [('GCA_000000001.1', '50.0'), ('GCF_000000002.1', '50.0')])
        self.assertEqual(plan['checkm_strain_heterogeneity_100'][1],
                         [('GCA_000000001.1', '12.5'), ('GCF_000000002.1', '0.0')])

    def test_columns_are_found_by_name_wherever_they_are(self):
        # a profile whose columns are in another order is read the same
        order = list(reversed(range(1, len(PROFILE_HEADER))))
        header = [PROFILE_HEADER[0]] + [PROFILE_HEADER[i] for i in order]
        rows = [[row[0]] + [row[i] for i in order]
                for row in (profile_row(GENOMES[0]), profile_row(GENOMES[1], '87.0', '3.4'))]
        profile = self.table('reordered.tsv', header, rows)
        _profile, qa, left_out = self.release_files()

        self.assertEqual(C.plan_checkm_import(profile, qa, left_out),
                         C.plan_checkm_import(*self.release_files()))

    def test_a_gzipped_table_is_read_as_a_plain_one_is(self):
        plain = self.release_files()
        profile = self.table('checkm.profiles.tsv.gz', PROFILE_HEADER,
                             [profile_row(GENOMES[0]), profile_row(GENOMES[1], '87.0', '3.4')],
                             compress=True)

        self.assertEqual(C.plan_checkm_import(profile, *plain[1:]), C.plan_checkm_import(*plain))

    def test_a_table_missing_a_column_is_refused_naming_it(self):
        header = [column for column in PROFILE_HEADER if column != 'Contamination']
        rows = [[value for column, value in zip(PROFILE_HEADER, profile_row(gid)) if column != 'Contamination']
                for gid in GENOMES]
        _profile, qa, left_out = self.release_files()
        profile = self.table('checkm.profiles.tsv', header, rows)

        with self.assertRaisesRegex(C.CheckMTableError, "'Contamination'"):
            C.plan_checkm_import(profile, qa, left_out)

    def test_a_genome_named_twice_is_refused_rather_than_one_value_chosen(self):
        # two runs concatenated would otherwise write whichever row came last
        with self.assertRaisesRegex(C.CheckMTableError, 'GCA_000000001.1'):
            self.plan(profile_rows=[profile_row(GENOMES[0], '98.5'), profile_row(GENOMES[1]),
                                    profile_row(GENOMES[0], '60.0')])

    def test_tables_of_different_genomes_are_refused_naming_the_genomes(self):
        # a genome's two strain heterogeneities would be of different releases
        with self.assertRaisesRegex(C.CheckMTableError, 'GCF_000000002.1') as raised:
            self.plan(qa_rows=[qa_row(GENOMES[0]), qa_row('GCA_000000003.1')])
        self.assertIn('GCA_000000003.1', str(raised.exception))

    def test_a_table_with_no_genomes_is_refused(self):
        with self.assertRaisesRegex(C.CheckMTableError, 'holds no genomes'):
            self.plan(profile_rows=[], qa_rows=[])

    def test_a_short_row_is_refused_naming_its_line(self):
        with self.assertRaisesRegex(C.CheckMTableError, 'line 3'):
            self.plan(profile_rows=[profile_row(GENOMES[0]), profile_row(GENOMES[1])[:5]])

    def test_the_genomes_checkm_left_out_are_read_from_its_not_assessed_file(self):
        plan = C.plan_checkm_import(*self.release_files(not_assessed=['GCA_977065575.1']))

        self.assertEqual(plan.not_assessed, ['GCA_977065575.1'])

    def test_a_not_assessed_file_naming_no_genome_leaves_none_out(self):
        # checkm writes the file, its header alone, where it left no genome out
        self.assertEqual(C.plan_checkm_import(*self.release_files()).not_assessed, [])

    def test_a_genome_both_assessed_and_not_is_refused(self):
        with self.assertRaisesRegex(C.CheckMTableError, 'GCF_000000002.1'):
            self.plan(not_assessed=[GENOMES[1]])


class WritingToTheDatabase(TempDirCase):
    def test_every_field_is_upserted_to_metadata_genes_in_one_commit(self):
        cursor = FakeCursor()
        manager = self.manager(cursor)
        manager.add_checkm_to_db(*self.release_files(), self.dir)

        upserts = cursor.upserts()
        self.assertEqual(len(upserts), 8)
        self.assertTrue(all(params[0] == C.CHECKM_TABLE
                            for sql, params in cursor.statements if 'upsert(' in sql))
        self.assertEqual(upserts['checkm_contamination'],
                         ('FLOAT', [('GCA_000000001.1', '1.2'), ('GCF_000000002.1', '3.4')]))
        self.assertEqual((manager.temp_con.commits, manager.temp_con.rollbacks), (1, 0))

    def test_a_genome_the_database_does_not_hold_is_refused_and_nothing_is_written(self):
        # --metadata dropped such a genome with a line on the console, and the
        # run succeeded without its estimates
        cursor = FakeCursor(genomes=[GENOMES[0]])
        manager = self.manager(cursor)

        with self.assertRaisesRegex(UnknownGenomesError, 'GCF_000000002.1'):
            manager.add_checkm_to_db(*self.release_files(), self.dir)
        self.assertEqual(cursor.upserts(), {})
        self.assertEqual((manager.temp_con.commits, manager.temp_con.rollbacks), (0, 1))

    def test_a_genome_checkm_left_out_has_the_estimates_of_its_earlier_sequences_cleared(self):
        # update_db kept the row of a new version or changed sequences, and with
        # it the estimates of what the genome was; checkm could not assess it
        stale = 'GCA_977065575.1'
        cursor = FakeCursor(genomes=GENOMES + (stale,), estimates={stale: (True, HELD)})
        manager = self.manager(cursor)
        manager.add_checkm_to_db(*self.release_files(not_assessed=[stale]), self.dir)

        self.assertEqual(cursor.cleared(), [[2]])
        update = [sql for sql, _params in cursor.statements if sql.startswith('UPDATE metadata_genes')][0]
        for field in C.CHECKM_FIELDS:
            self.assertIn('{} = NULL'.format(field), update)
        self.assertIn('Cleared the CheckM estimates of 1 genome(s)', self.warnings[0])
        self.assertIn(stale, self.warnings[0])
        with open(os.path.join(self.dir, C.CLEARED_NAME)) as handle:
            self.assertEqual(handle.read().splitlines(),
                             ['\t'.join(('genome_id',) + C.CHECKM_FIELDS),
                              '\t'.join([stale] + [str(v) for v in HELD])])
        self.assertEqual((manager.temp_con.commits, manager.temp_con.rollbacks), (1, 0))

    def test_a_genome_left_out_that_holds_no_estimates_has_nothing_cleared(self):
        # a genome new to the database has no row of metadata_genes
        new = 'GCA_977065575.1'
        cursor = FakeCursor(genomes=GENOMES + (new,), estimates={new: (True, (None,) * 8)})
        manager = self.manager(cursor)
        manager.add_checkm_to_db(*self.release_files(not_assessed=[new]), self.dir)

        self.assertEqual(cursor.cleared(), [])
        with open(os.path.join(self.dir, C.CLEARED_NAME)) as handle:
            self.assertEqual(handle.read().splitlines(), ['\t'.join(('genome_id',) + C.CHECKM_FIELDS)])

    def test_a_genome_left_out_that_this_update_did_not_change_is_warned_of_and_kept(self):
        # its estimates are of the sequences it has: the file is of another release
        kept = 'GCA_000000003.1'
        cursor = FakeCursor(genomes=GENOMES + (kept,), estimates={kept: (False, HELD)})
        manager = self.manager(cursor)
        manager.add_checkm_to_db(*self.release_files(not_assessed=[kept]), self.dir)

        self.assertEqual(cursor.cleared(), [])
        self.assertIn('not changed by this update', self.warnings[0])
        self.assertIn(kept, self.warnings[0])

    def test_a_changed_genome_absent_from_the_not_assessed_file_is_not_cleared(self):
        # CheckM run again over a few genomes and loaded on its own clears
        # nothing of the genomes the release's run assessed
        versioned = 'GCF_000000009.2'
        cursor = FakeCursor(genomes=GENOMES + (versioned,), estimates={versioned: (True, HELD)})
        manager = self.manager(cursor)
        manager.add_checkm_to_db(*self.release_files(), self.dir)

        self.assertEqual(cursor.cleared(), [])

    def test_estimates_are_cleared_only_once_every_field_is_written(self):
        # a genome the database does not hold refuses the run before anything is cleared
        stale = 'GCA_977065575.1'
        cursor = FakeCursor(genomes=(GENOMES[0], stale), estimates={stale: (True, HELD)})
        manager = self.manager(cursor)
        with self.assertRaises(UnknownGenomesError):
            manager.add_checkm_to_db(*self.release_files(not_assessed=[stale]), self.dir)

        self.assertEqual(cursor.cleared(), [])
        self.assertEqual((manager.temp_con.commits, manager.temp_con.rollbacks), (0, 1))

    def test_tables_that_cannot_be_written_are_refused_before_the_database_is_asked(self):
        cursor = FakeCursor()
        manager = self.manager(cursor)

        with self.assertRaises(C.CheckMTableError):
            manager.add_checkm_to_db(*self.release_files(qa_rows=[qa_row(GENOMES[0])]), self.dir)
        self.assertEqual(cursor.statements, [])
        self.assertEqual((manager.temp_con.commits, manager.temp_con.rollbacks), (0, 1))


class GenomesWithNoEstimates(TempDirCase):
    """A genome checkm never planned is in neither its tables nor its not-assessed file."""

    def missing_file(self):
        with open(os.path.join(self.dir, C.MISSING_NAME)) as handle:
            return [line.split('\t') for line in handle.read().splitlines()]

    def test_a_genome_with_no_estimates_that_checkm_did_not_leave_out_is_warned_of_and_written(self):
        # r237: added by update_db as unchanged, so never planned by checkm
        unplanned = 'GCA_001341675.1'
        cursor = FakeCursor(genomes=GENOMES + (unplanned,), missing=[(unplanned, True, True)])
        manager = self.manager(cursor)
        manager.add_checkm_to_db(*self.release_files(), self.dir)

        self.assertEqual(self.missing_file(), [list(C.MISSING_HEADER), [unplanned, C.STATUS_NEW]])
        self.assertIn('1 genome(s) have no CheckM estimates', self.warnings[0])
        self.assertIn('1 new', self.warnings[0])
        self.assertIn(unplanned, self.warnings[0])
        self.assertEqual((manager.temp_con.commits, manager.temp_con.rollbacks), (1, 0))

    def test_a_genome_checkm_left_out_is_not_warned_of_again(self):
        # checkm_not_assessed.tsv already says why it has none
        left_out = 'GCA_977065575.1'
        cursor = FakeCursor(genomes=GENOMES + (left_out,), missing=[(left_out, True, True)])
        manager = self.manager(cursor)
        manager.add_checkm_to_db(*self.release_files(not_assessed=[left_out]), self.dir)

        self.assertEqual(self.missing_file(), [list(C.MISSING_HEADER)])
        self.assertEqual(self.warnings, [])

    def test_each_genome_says_whether_this_update_made_it_new_or_updated_it(self):
        cursor = FakeCursor(missing=[('GCA_000000005.1', True, True),     # added, or a new version
                                     ('GCF_000000006.1', True, False),    # sequences changed
                                     ('GCA_000000007.1', False, False)])  # left as it was
        manager = self.manager(cursor)
        manager.add_checkm_to_db(*self.release_files(), self.dir)

        self.assertEqual(self.missing_file()[1:], [['GCA_000000005.1', C.STATUS_NEW],
                                                   ['GCA_000000007.1', C.STATUS_UNCHANGED],
                                                   ['GCF_000000006.1', C.STATUS_UPDATED]])
        self.assertIn('1 new, 1 updated, 1 unchanged', self.warnings[0])

    def test_genomes_with_no_estimates_are_asked_for_once_every_field_is_written_and_cleared(self):
        # what the database will hold, not what it held before the update
        stale = 'GCA_977065575.1'
        cursor = FakeCursor(genomes=GENOMES + (stale,), estimates={stale: (True, HELD)})
        manager = self.manager(cursor)
        manager.add_checkm_to_db(*self.release_files(not_assessed=[stale]), self.dir)

        statements = [sql for sql, _params in cursor.statements]
        asked = [i for i, sql in enumerate(statements) if sql.startswith('SELECT g.id_at_source, g.has_changed')]
        written = [i for i, sql in enumerate(statements)
                   if 'upsert(' in sql or sql.startswith('UPDATE metadata_genes')]
        self.assertEqual(len(asked), 1)
        self.assertGreater(asked[0], max(written))

    def test_only_ncbi_genomes_are_asked_for(self):
        # user genomes are not CheckM'd by the release
        cursor = FakeCursor()
        manager = self.manager(cursor)
        manager.add_checkm_to_db(*self.release_files(), self.dir)

        params = [params for sql, params in cursor.statements
                  if sql.startswith('SELECT g.id_at_source, g.has_changed')][0]
        self.assertEqual(sorted(params[0]), ['GenBank', 'RefSeq'])

    def test_the_file_is_written_with_its_header_alone_where_no_genome_is_missing(self):
        manager = self.manager(FakeCursor())
        manager.add_checkm_to_db(*self.release_files(), self.dir)

        self.assertEqual(self.missing_file(), [list(C.MISSING_HEADER)])
        self.assertEqual(self.warnings, [])


class TheCommandLine(TempDirCase):
    def test_update_checkm_db_takes_checkms_release_files_and_no_metadata_file(self):
        profile, qa, left_out = self.release_files()
        out_dir = os.path.join(self.dir, 'update_checkm_db')
        options = main_module.get_main_parser().parse_args(
            ['update_checkm_db', '--db_service', 'gtdb_r237', '-f', profile, '-q', qa,
             '-n', left_out, '-o', out_dir, '-l', os.path.join(self.dir, 'run.log')])
        self.assertFalse(hasattr(options, 'metadata'))

        with mock.patch.object(main_py, 'CheckMDatabaseManager') as manager:
            main_py.OptionsParser().parse_options(options)
        manager.assert_called_once_with({'service': 'gtdb_r237'})
        manager.return_value.add_checkm_to_db.assert_called_once_with(profile, qa, left_out, out_dir)
        self.assertTrue(os.path.isdir(out_dir))

    def test_an_out_dir_is_required(self):
        # it is where the genomes with no estimates are written
        profile, qa, left_out = self.release_files()
        with mock.patch('sys.stderr'), self.assertRaises(SystemExit):
            main_module.get_main_parser().parse_args(
                ['update_checkm_db', '--db_service', 'gtdb_r237', '-f', profile, '-q', qa,
                 '-n', left_out, '-l', os.path.join(self.dir, 'run.log')])

    def test_a_metadata_file_is_no_longer_accepted(self):
        profile, qa, left_out = self.release_files()
        with mock.patch('sys.stderr'), self.assertRaises(SystemExit):
            main_module.get_main_parser().parse_args(
                ['update_checkm_db', '--db_service', 'gtdb_r237', '-f', profile, '-q', qa,
                 '-n', left_out, '-o', self.dir, '-m', 'metadata.tsv',
                 '-l', os.path.join(self.dir, 'run.log')])



# ------------------------------------------------------------- update_checkm2_db

# the columns of CheckM2's quality_report.tsv as checkm2 gathers it for the release
QUALITY_REPORT_HEADER = ['Name', 'Completeness', 'Contamination', 'Completeness_Model_Used',
                         'Additional_Notes']
SPECIFIC = 'Neural Network (Specific Model)'
GENERAL = 'Gradient Boost (General Model)'

# the CheckM2 estimates a genome holds in metadata_genes, one value for each of C.CHECKM2_FIELDS
HELD2 = (99.1, 0.4, 'Specific')


def report_row(accession, completeness='100.0', contamination='2.76', model=SPECIFIC):
    return [accession, completeness, contamination, model, 'None']


class CheckM2Case(TempDirCase):
    def checkm2_files(self, rows=None, not_assessed=(), header=QUALITY_REPORT_HEADER, compress=True):
        rows = rows if rows is not None else [report_row(GENOMES[0]),
                                              report_row(GENOMES[1], '87.0', '3.4', GENERAL)]
        by_name = [dict(zip(QUALITY_REPORT_HEADER, row)) for row in rows]
        report = self.table('checkm2.quality_report.tsv.gz' if compress else 'checkm2.quality_report.tsv',
                            header, [[row[column] for column in header] for row in by_name],
                            compress=compress)
        left_out = self.table('checkm2_not_assessed.tsv', ['genome_id', 'reason', 'detail'],
                              [(gid, 'no_protein_file', '') for gid in not_assessed])
        return report, left_out


class PlanningTheCheckM2Import(CheckM2Case):
    def test_completeness_contamination_and_the_model_are_written_for_every_genome(self):
        genomes, fields, not_assessed = C.plan_checkm2_import(*self.checkm2_files())

        self.assertEqual(genomes, list(GENOMES))
        self.assertEqual(not_assessed, [])
        self.assertEqual(fields, [
            ('checkm2_completeness', 'FLOAT', [('GCA_000000001.1', '100.0'), ('GCF_000000002.1', '87.0')]),
            ('checkm2_contamination', 'FLOAT', [('GCA_000000001.1', '2.76'), ('GCF_000000002.1', '3.4')]),
            ('checkm2_model', 'TEXT', [('GCA_000000001.1', 'Specific'), ('GCF_000000002.1', 'General')])])

    def test_columns_are_found_by_name_wherever_they_are(self):
        reordered = ['Name', 'Additional_Notes', 'Completeness_Model_Used', 'Contamination', 'Completeness']
        self.assertEqual(C.plan_checkm2_import(*self.checkm2_files(header=reordered)),
                         C.plan_checkm2_import(*self.checkm2_files()))

    def test_a_plain_report_is_read_as_a_gzipped_one_is(self):
        self.assertEqual(C.plan_checkm2_import(*self.checkm2_files(compress=False)),
                         C.plan_checkm2_import(*self.checkm2_files()))

    def test_a_model_naming_neither_general_nor_specific_is_refused_naming_it(self):
        # the database would otherwise hold a value nothing reading it expects
        with self.assertRaisesRegex(C.CheckMTableError, 'Transformer \\(Universal Model\\)') as raised:
            C.plan_checkm2_import(*self.checkm2_files(rows=[
                report_row(GENOMES[0]), report_row(GENOMES[1], model='Transformer (Universal Model)')]))
        self.assertIn('GCF_000000002.1', str(raised.exception))

    def test_a_report_missing_a_column_is_refused_naming_it(self):
        header = [column for column in QUALITY_REPORT_HEADER if column != 'Completeness_Model_Used']
        with self.assertRaisesRegex(C.CheckMTableError, "'Completeness_Model_Used'"):
            C.plan_checkm2_import(*self.checkm2_files(header=header))

    def test_a_genome_both_assessed_and_not_is_refused(self):
        with self.assertRaisesRegex(C.CheckMTableError, 'CheckM2 tables'):
            C.plan_checkm2_import(*self.checkm2_files(not_assessed=[GENOMES[1]]))


class WritingCheckM2ToTheDatabase(CheckM2Case):
    def manager(self, cursor):
        manager = super().manager(cursor)
        manager.__class__ = C.CheckM2DatabaseManager
        return manager

    def test_the_three_fields_are_upserted_to_metadata_genes_in_one_commit(self):
        cursor = FakeCursor()
        manager = self.manager(cursor)
        manager.add_checkm2_to_db(*self.checkm2_files(), self.dir)

        upserts = cursor.upserts()
        self.assertEqual(sorted(upserts), sorted(C.CHECKM2_FIELDS))
        self.assertEqual(upserts['checkm2_model'],
                         ('TEXT', [('GCA_000000001.1', 'Specific'), ('GCF_000000002.1', 'General')]))
        self.assertEqual((manager.temp_con.commits, manager.temp_con.rollbacks), (1, 0))

    def test_a_genome_checkm2_left_out_has_its_earlier_checkm2_estimates_cleared(self):
        stale = 'GCA_977065575.1'
        cursor = FakeCursor(genomes=GENOMES + (stale,), estimates={stale: (True, HELD2)})
        manager = self.manager(cursor)
        manager.add_checkm2_to_db(*self.checkm2_files(not_assessed=[stale]), self.dir)

        self.assertEqual(cursor.cleared(), [[2]])
        update = [sql for sql, _params in cursor.statements if sql.startswith('UPDATE metadata_genes')][0]
        self.assertEqual(update.count('= NULL'), len(C.CHECKM2_FIELDS))
        for field in C.CHECKM2_FIELDS:
            self.assertIn('{} = NULL'.format(field), update)
        self.assertNotIn('checkm_completeness', update)
        self.assertIn('Cleared the CheckM2 estimates of 1 genome(s) checkm2 left out', self.warnings[0])
        with open(os.path.join(self.dir, C.CHECKM2_CLEARED_NAME)) as handle:
            self.assertEqual(handle.read().splitlines(),
                             ['\t'.join(('genome_id',) + C.CHECKM2_FIELDS),
                              '\t'.join([stale] + [str(v) for v in HELD2])])

    def test_genomes_with_no_checkm2_estimates_are_written_to_their_own_file(self):
        unplanned = 'GCA_001341675.1'
        cursor = FakeCursor(genomes=GENOMES + (unplanned,), missing=[(unplanned, True, True)])
        manager = self.manager(cursor)
        manager.add_checkm2_to_db(*self.checkm2_files(), self.dir)

        asked = [sql for sql, _params in cursor.statements
                 if sql.startswith('SELECT g.id_at_source, g.has_changed')][0]
        self.assertIn('m.checkm2_completeness IS NULL', asked)
        with open(os.path.join(self.dir, C.CHECKM2_MISSING_NAME)) as handle:
            self.assertEqual(handle.read().splitlines(), ['genome_id\tstatus', unplanned + '\tnew'])
        self.assertIn('have no CheckM2 estimates', self.warnings[0])
        self.assertFalse(os.path.exists(os.path.join(self.dir, C.MISSING_NAME)))

    def test_a_report_that_cannot_be_written_is_refused_before_the_database_is_asked(self):
        cursor = FakeCursor()
        manager = self.manager(cursor)
        with self.assertRaises(C.CheckMTableError):
            manager.add_checkm2_to_db(*self.checkm2_files(rows=[report_row(GENOMES[0], model='?')]), self.dir)
        self.assertEqual(cursor.statements, [])
        self.assertEqual((manager.temp_con.commits, manager.temp_con.rollbacks), (0, 1))


class TheCheckM2CommandLine(CheckM2Case):
    def test_update_checkm2_db_takes_checkm2s_release_files(self):
        report, left_out = self.checkm2_files()
        out_dir = os.path.join(self.dir, 'update_checkm2_db')
        options = main_module.get_main_parser().parse_args(
            ['update_checkm2_db', '--db_service', 'gtdb_r237', '-f', report, '-n', left_out,
             '-o', out_dir, '-l', os.path.join(self.dir, 'run.log')])

        with mock.patch.object(main_py, 'CheckM2DatabaseManager') as manager:
            main_py.OptionsParser().parse_options(options)
        manager.assert_called_once_with({'service': 'gtdb_r237'})
        manager.return_value.add_checkm2_to_db.assert_called_once_with(report, left_out, out_dir)
        self.assertTrue(os.path.isdir(out_dir))

    def test_the_profile_files_are_f_for_both_commands_never_c(self):
        # -c is --cpus throughout the toolkit, and -p is --password in every *_db command
        profile, qa, left_out = self.release_files()
        report, left_out2 = self.checkm2_files()
        common = ['--db_service', 'gtdb_r237', '-o', self.dir, '-l', os.path.join(self.dir, 'run.log')]
        parser = main_module.get_main_parser()

        options = parser.parse_args(['update_checkm_db', '--checkm_profile_file', profile, '-q', qa,
                                     '-n', left_out] + common)
        self.assertEqual(options.checkm_profile_file, profile)
        options = parser.parse_args(['update_checkm2_db', '--checkm2_profile_file', report,
                                     '-n', left_out2] + common)
        self.assertEqual(options.checkm2_profile_file, report)
        for argv in (['update_checkm_db', '-c', profile, '-q', qa, '-n', left_out],
                     ['update_checkm2_db', '-c', report, '-n', left_out2],
                     ['update_checkm2_db', '-r', report, '-n', left_out2]):
            with mock.patch('sys.stderr'), self.assertRaises(SystemExit) as ended:
                parser.parse_args(argv + common)
            self.assertEqual(ended.exception.code, 2)

    def test_join_checkm_join_checkm2_and_prepare_checkm2_are_gone(self):
        # checkm and checkm2 batch, run and gather CheckM and CheckM2 themselves, and
        # the update_*_db commands load what they write
        for command in ('join_checkm', 'join_checkm2', 'prepare_checkm2'):
            with mock.patch('sys.stdout'), mock.patch('sys.stderr') as stderr, \
                    self.assertRaises(SystemExit) as ended:
                main_module.get_main_parser().parse_args([command, '-h'])
            self.assertEqual(ended.exception.code, 2)
            self.assertIn('invalid choice', ''.join(call.args[0] for call in stderr.write.call_args_list))


if __name__ == '__main__':
    unittest.main()
