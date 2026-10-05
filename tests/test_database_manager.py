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

"""update_db decides every genome from report.log and genome_dirs.tsv. What it
decides, where it says a file is, and what it hashes are tested here without a
database; the statements are tested as they are handed to psycopg2.
"""

import gzip
import hashlib
import logging
import os
import shutil
import tempfile
import unittest
from datetime import datetime
from unittest import mock

from gtdb_migration_tk import database_manager as M
from gtdb_migration_tk import update_genomes as UG
from gtdb_migration_tk.update_genomes import (STATUS_FASTA_CHANGED,
                                             STATUS_FASTA_UNCHANGED, STATUS_NEW,
                                             STATUS_REMOVED,
                                             STATUS_SEQUENCES_UNCHANGED,
                                             STATUS_TO_CURATE)

SOURCES = {'GCF': 2, 'GCA': 3}


class TempDirCase(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='database_manager_test.')
        self.addCleanup(shutil.rmtree, self.dir, True)
        self.release = os.path.join(self.dir, 'release237')
        logging.getLogger('timestamp').addHandler(logging.NullHandler())

    def genome_dir(self, accession):
        """Where a release keeps a genome, as update_genomes lays it out."""
        database = 'refseq' if accession.startswith('GCF') else 'genbank'
        digits = accession[4:13]
        return os.path.join(self.release, database, accession[:3], digits[0:3], digits[3:6],
                            digits[6:9], accession + '_ASM1v1')

    def make_genome(self, accession, genomic=b'>c1\nACGT\n', proteins=b'>g1\nMK\n'):
        gdir = self.genome_dir(accession)
        os.makedirs(os.path.join(gdir, 'prodigal'))
        with gzip.open(os.path.join(gdir, accession + '_ASM1v1_genomic.fna.gz'), 'wb') as handle:
            handle.write(genomic)
        if proteins is not None:
            with gzip.open(os.path.join(gdir, 'prodigal', accession + '_protein.faa.gz'), 'wb') as handle:
                handle.write(proteins)
        return gdir

    def location(self, accession, genes=False):
        rel = os.path.relpath(self.genome_dir(accession), os.path.join(self.release, 'genbank'
                              if accession.startswith('GCA') else 'refseq'))
        if genes:
            return '{}/prodigal/{}_protein.faa.gz'.format(rel, accession)
        return '{}/{}_ASM1v1_genomic.fna.gz'.format(rel, accession)

    def row(self, accession, id=1, fasta_hash='old', genes_hash='old', genes=True):
        return M.DatabaseGenome(id, accession, SOURCES[accession[:3]],
                                self.location(accession), fasta_hash,
                                self.location(accession, genes=True) if genes else None,
                                genes_hash if genes else None)

    def plan(self, outcomes, database=(), rehash_all=False, no_genomic_fasta=frozenset()):
        dirs = {acc: self.genome_dir(acc) for acc, outcome in outcomes.items()
                if outcome in UG.STATUS_IN_RELEASE}
        decisions = M.plan_update(outcomes, dirs, {r.name: r for r in database}, rehash_all,
                                  no_genomic_fasta)
        return {d.accession: d for d in decisions}


def content_sha1(data):
    return hashlib.sha1(data).hexdigest()


# ------------------------------------------------------------ deciding

class DecidingEachGenome(TempDirCase):
    """Each genome is decided from its outcome and whether the database holds it."""

    def test_a_new_genome_is_added_with_both_files_hashed(self):
        d = self.plan({'GCA_000000001.1': STATUS_NEW})['GCA_000000001.1']

        self.assertEqual(d.action, M.ACTION_ADDED)
        self.assertTrue(d.hash_genomic and d.hash_genes)

    def test_a_new_version_takes_over_its_predecessors_row(self):
        old = self.row('GCA_000000001.2', id=7)
        plan = self.plan({'GCA_000000001.3': STATUS_NEW, 'GCA_000000001.2': STATUS_REMOVED}, [old])

        self.assertEqual(plan['GCA_000000001.3'].action, M.ACTION_VERSIONED)
        self.assertEqual(plan['GCA_000000001.3'].row.id, 7)
        self.assertEqual(plan['GCA_000000001.2'].action, M.ACTION_REPLACED)

    def test_a_predecessor_still_in_the_release_is_kept_and_the_new_version_added(self):
        old = self.row('GCA_000000001.1')
        plan = self.plan({'GCA_000000001.2': STATUS_NEW,
                          'GCA_000000001.1': STATUS_FASTA_UNCHANGED}, [old])

        self.assertEqual(plan['GCA_000000001.2'].action, M.ACTION_ADDED)
        self.assertEqual(plan['GCA_000000001.1'].action, M.ACTION_UPDATED)

    def test_a_refseq_genome_is_never_taken_for_a_version_of_a_genbank_one(self):
        old = self.row('GCA_000000001.1')
        plan = self.plan({'GCF_000000001.2': STATUS_NEW, 'GCA_000000001.1': STATUS_REMOVED}, [old])

        self.assertEqual(plan['GCF_000000001.2'].action, M.ACTION_ADDED)
        self.assertEqual(plan['GCA_000000001.1'].action, M.ACTION_DELETED)

    def test_a_genome_whose_sequences_changed_under_its_accession_is_said_to_have(self):
        d = self.plan({'GCA_000000001.1': STATUS_FASTA_CHANGED},
                      [self.row('GCA_000000001.1')])['GCA_000000001.1']

        self.assertEqual(d.action, M.ACTION_SEQUENCES_CHANGED)
        self.assertTrue(d.hash_genomic and d.hash_genes)

    def test_rewritten_deflines_rehash_the_genomic_fasta_and_not_the_proteins(self):
        d = self.plan({'GCA_000000001.1': STATUS_SEQUENCES_UNCHANGED},
                      [self.row('GCA_000000001.1')])['GCA_000000001.1']

        self.assertEqual(d.action, M.ACTION_UPDATED)
        self.assertEqual((d.hash_genomic, d.hash_genes), (True, False))

    def test_an_unchanged_genome_is_hashed_only_when_every_genome_is(self):
        outcomes = {'GCA_000000001.1': STATUS_FASTA_UNCHANGED}
        database = [self.row('GCA_000000001.1')]

        d = self.plan(outcomes, database)['GCA_000000001.1']
        self.assertEqual((d.hash_genomic, d.hash_genes), (False, False))
        d = self.plan(outcomes, database, rehash_all=True)['GCA_000000001.1']
        self.assertEqual((d.hash_genomic, d.hash_genes), (True, True))

    def test_a_protein_file_the_database_has_no_hash_for_is_hashed(self):
        d = self.plan({'GCA_000000001.1': STATUS_FASTA_UNCHANGED},
                      [self.row('GCA_000000001.1', genes=False)])['GCA_000000001.1']

        self.assertEqual((d.hash_genomic, d.hash_genes), (False, True))

    def test_a_genome_of_the_release_the_database_lacks_is_added_whatever_its_outcome(self):
        for outcome in (STATUS_FASTA_UNCHANGED, STATUS_SEQUENCES_UNCHANGED, STATUS_FASTA_CHANGED):
            d = self.plan({'GCA_000000001.1': outcome})['GCA_000000001.1']
            self.assertEqual(d.action, M.ACTION_ADDED, outcome)

    def test_a_removed_genome_is_deleted(self):
        d = self.plan({'GCA_000000001.1': STATUS_REMOVED},
                      [self.row('GCA_000000001.1')])['GCA_000000001.1']

        self.assertEqual(d.action, M.ACTION_DELETED)

    def test_a_genome_to_curate_is_deleted(self):
        # it was never copied into the release, so its row names a file the
        # release does not hold
        d = self.plan({'GCF_030217685.1': STATUS_TO_CURATE},
                      [self.row('GCF_030217685.1')])['GCF_030217685.1']

        self.assertEqual(d.action, M.ACTION_DELETED)

    def test_a_genome_leaving_that_the_database_never_held_needs_nothing(self):
        d = self.plan({'GCA_000000001.1': STATUS_REMOVED})['GCA_000000001.1']

        self.assertEqual(d.action, M.ACTION_NOT_IN_DATABASE)

    def test_a_genome_of_the_database_the_report_does_not_name_is_left_alone(self):
        plan = self.plan({'GCA_000000001.1': STATUS_NEW}, [self.row('GCA_000000009.1')])

        self.assertEqual(plan['GCA_000000009.1'].action, M.ACTION_NOT_IN_REPORT)

    def test_a_new_genome_ncbi_published_without_a_genomic_fasta_is_not_added(self):
        # GCA_056491145.1 in r237: its directory holds NCBI's reports alone
        d = self.plan({'GCA_056491145.1': STATUS_NEW},
                      no_genomic_fasta={'GCA_056491145.1'})['GCA_056491145.1']

        self.assertEqual(d.action, M.ACTION_NO_GENOMIC_FASTA)
        self.assertFalse(d.hash_genomic or d.hash_genes)

    def test_a_new_version_without_a_genomic_fasta_does_not_take_over_its_predecessor(self):
        plan = self.plan({'GCA_000000001.3': STATUS_NEW, 'GCA_000000001.2': STATUS_REMOVED},
                         [self.row('GCA_000000001.2', id=7)],
                         no_genomic_fasta={'GCA_000000001.3'})

        self.assertEqual(plan['GCA_000000001.3'].action, M.ACTION_NO_GENOMIC_FASTA)
        self.assertEqual(plan['GCA_000000001.2'].action, M.ACTION_DELETED)

    def test_every_outcome_update_genomes_writes_is_decided(self):
        statuses = {value for name, value in vars(UG).items()
                    if name.startswith('STATUS_') and isinstance(value, str)}
        self.assertEqual(statuses - M.KNOWN_OUTCOMES, set())

    def test_an_outcome_it_does_not_know_is_refused(self):
        with self.assertRaises(M.UpdateDbError):
            self.plan({'GCA_000000001.1': 'modified'})

    def test_a_report_and_genome_dirs_of_different_releases_are_refused(self):
        with self.assertRaises(M.UpdateDbError):
            M.plan_update({'GCA_000000001.1': STATUS_NEW},
                          {'GCA_000000002.1': self.genome_dir('GCA_000000002.1')}, {})

    def test_an_accession_of_neither_database_is_refused(self):
        with self.assertRaises(M.UpdateDbError):
            M.plan_update({'XYZ_000000001.1': STATUS_REMOVED}, {}, {})


# ------------------------------------------------------------ files

class WhereTheDatabaseSaysAFileIs(TempDirCase):
    def test_a_file_is_recorded_below_the_release_from_the_prefix_down(self):
        gdir = self.genome_dir('GCA_024206075.2')
        path = os.path.join(gdir, 'GCA_024206075.2_ASM1v1_genomic.fna.gz')

        self.assertEqual(M.database_location(path, gdir, 'GCA_024206075.2'),
                         'GCA/024/206/075/GCA_024206075.2_ASM1v1/GCA_024206075.2_ASM1v1_genomic.fna.gz')

    def test_a_directory_not_laid_out_as_ncbi_lays_it_out_is_refused(self):
        gdir = os.path.join(self.dir, 'flat', 'GCA_024206075.2_ASM1v1')
        with self.assertRaises(M.UpdateDbError):
            M.database_location(os.path.join(gdir, 'x'), gdir, 'GCA_024206075.2')


class TheHash(TempDirCase):
    def test_it_is_sha1_of_the_content_not_of_the_gzip(self):
        path = os.path.join(self.dir, 'a.gz')
        with gzip.open(path, 'wb') as handle:
            handle.write(b'>g1\nMK\n')

        self.assertEqual(M.content_hash(path), content_sha1(b'>g1\nMK\n'))

    def test_the_same_content_gzipped_differently_hashes_the_same(self):
        paths, gzipped = [], []
        for mtime, level in ((0, 9), (12345, 1)):
            path = os.path.join(self.dir, '{}.gz'.format(mtime))
            gzipped.append(gzip.compress(b'>g1\nMKLV\n' * 100, compresslevel=level, mtime=mtime))
            with open(path, 'wb') as handle:
                handle.write(gzipped[-1])
            paths.append(path)

        self.assertNotEqual(gzipped[0], gzipped[1])
        self.assertEqual(M.content_hash(paths[0]), M.content_hash(paths[1]))

    def test_a_file_that_is_not_there_has_none(self):
        self.assertIsNone(M.content_hash(os.path.join(self.dir, 'gone.gz')))


class TheFilesRecorded(TempDirCase):
    """What a genome's row says once its files are hashed."""

    def manager(self):
        return M.DatabaseManager('host', 'user', 'pw', 'db', cpus=1)

    def finish(self, outcomes, database=(), rehash_all=False):
        decisions = list(self.plan(outcomes, database, rehash_all).values())
        manager = self.manager()
        hashes = manager.hash_files(decisions)
        return {d.accession: d for d in manager.finish_decisions(decisions, hashes)}, hashes

    def test_an_unchanged_genome_whose_row_agrees_is_left_unchanged(self):
        self.make_genome('GCA_000000001.1')
        finished, _ = self.finish({'GCA_000000001.1': STATUS_FASTA_UNCHANGED},
                                  [self.row('GCA_000000001.1')])

        self.assertEqual(finished['GCA_000000001.1'].action, M.ACTION_UNCHANGED)

    def test_a_stale_hash_is_put_right_by_rehashing_and_named(self):
        self.make_genome('GCA_000000001.1')
        finished, _ = self.finish({'GCA_000000001.1': STATUS_FASTA_UNCHANGED},
                                  [self.row('GCA_000000001.1', genes_hash=content_sha1(b'>g1\nMK\n'))],
                                  rehash_all=True)

        self.assertEqual(finished['GCA_000000001.1'].action, M.ACTION_UPDATED)
        self.assertEqual(finished['GCA_000000001.1'].detail, 'fasta_file_sha256')

    def test_a_genome_whose_directory_moved_is_updated_with_its_new_path(self):
        self.make_genome('GCA_000000001.1')
        row = self.row('GCA_000000001.1')._replace(fasta_location='GCA/elsewhere_genomic.fna.gz')
        finished, hashes = self.finish({'GCA_000000001.1': STATUS_FASTA_UNCHANGED}, [row])

        d = finished['GCA_000000001.1']
        self.assertEqual(d.action, M.ACTION_UPDATED)
        self.assertEqual(M.genome_files(d, hashes).fasta_location, self.location('GCA_000000001.1'))

    def test_a_genome_with_no_protein_file_is_recorded_with_none(self):
        self.make_genome('GCA_000000001.1', proteins=None)
        finished, hashes = self.finish({'GCA_000000001.1': STATUS_NEW})

        d = finished['GCA_000000001.1']
        files = M.genome_files(d, hashes)
        self.assertEqual((files.genes_location, files.genes_hash), (None, None))
        self.assertIn('no protein file', d.detail)

    def test_the_genomes_without_a_genomic_fasta_are_found(self):
        self.make_genome('GCA_000000001.1')
        os.makedirs(self.genome_dir('GCA_000000002.1'))
        dirs = {acc: self.genome_dir(acc) for acc in ('GCA_000000001.1', 'GCA_000000002.1')}

        self.assertEqual(M.missing_genomic_fasta(dirs, dirs, threads=2), {'GCA_000000002.1'})

    def test_a_genome_the_database_holds_that_lost_its_genomic_fasta_is_refused(self):
        # a release that lost a file, not a genome NCBI never published
        os.makedirs(self.genome_dir('GCA_000000001.1'))
        with self.assertRaises(M.UpdateDbError) as caught:
            self.finish({'GCA_000000001.1': STATUS_SEQUENCES_UNCHANGED},
                        [self.row('GCA_000000001.1')])
        self.assertIn('lost the file', str(caught.exception))


# ------------------------------------------------------------ the statements

class FakeCursor(object):
    def __init__(self, rows=()):
        self.statements = []
        self.rows = list(rows)

    def execute(self, sql, params=None):
        self.statements.append((' '.join(sql.split()), params))

    def fetchall(self):
        return self.rows


class WritingTheDatabase(TempDirCase):
    """What apply() hands psycopg2, in the order it hands it."""

    DATE = datetime(2026, 10, 6)

    def apply(self, outcomes, database=()):
        for accession, outcome in outcomes.items():
            if outcome in UG.STATUS_IN_RELEASE:
                self.make_genome(accession)
        manager = M.DatabaseManager('host', 'user', 'pw', 'db', cpus=1)
        decisions = list(self.plan(outcomes, database).values())
        hashes = manager.hash_files(decisions)
        decisions = manager.finish_decisions(decisions, hashes)

        cur, values = FakeCursor(), []
        with mock.patch.object(M, 'execute_values',
                               lambda cur, sql, rows, **kw: values.append((' '.join(sql.split()), rows))):
            manager.apply(cur, decisions, hashes, SOURCES, self.DATE)
        return cur.statements, values

    def test_has_changed_is_cleared_on_every_ncbi_genome_before_anything_else(self):
        statements, _ = self.apply({'GCA_000000001.1': STATUS_NEW})

        sql, params = statements[0]
        self.assertTrue(sql.startswith('UPDATE genomes SET has_changed = FALSE'))
        self.assertEqual(sorted(params[0]), [2, 3])

    def test_a_genome_is_added_under_the_source_of_its_prefix_with_has_changed(self):
        _, values = self.apply({'GCA_000000001.1': STATUS_NEW, 'GCF_000000002.1': STATUS_NEW})

        sql, rows = [v for v in values if v[0].startswith('INSERT')][0]
        by_name = {r[0]: r for r in rows}
        self.assertEqual(by_name['GCA_000000001.1'][6], 3)
        self.assertEqual(by_name['GCF_000000002.1'][6], 2)
        self.assertEqual(by_name['GCA_000000001.1'][9], True)
        self.assertEqual(by_name['GCA_000000001.1'][13], 'G000000001')
        self.assertEqual(by_name['GCA_000000001.1'][5], content_sha1(b'>c1\nACGT\n'))

    def test_a_new_version_renames_its_predecessors_row_and_loses_its_markers(self):
        statements, values = self.apply({'GCA_000000001.3': STATUS_NEW,
                                         'GCA_000000001.2': STATUS_REMOVED},
                                        [self.row('GCA_000000001.2', id=7)])

        sql, rows = [v for v in values if 'SET name = v.name' in v[0]][0]
        self.assertEqual(rows[0][:2], (7, 'GCA_000000001.3'))
        self.assertIn(('DELETE FROM aligned_markers WHERE genome_id = ANY(%s)', ([7],)), statements)
        self.assertFalse([s for s in statements if s[0].startswith('DELETE FROM genomes')])

    def test_a_removed_genome_is_deleted_by_id(self):
        statements, _ = self.apply({'GCA_000000001.1': STATUS_REMOVED},
                                   [self.row('GCA_000000001.1', id=42)])

        self.assertIn(('DELETE FROM genomes WHERE id = ANY(%s)', ([42],)), statements)

    def test_a_genome_whose_sequences_changed_loses_its_markers_and_is_flagged(self):
        statements, values = self.apply({'GCA_000000001.1': STATUS_FASTA_CHANGED},
                                        [self.row('GCA_000000001.1', id=5)])

        sql, rows = [v for v in values if 'v.has_changed' in v[0]][0]
        self.assertEqual((rows[0][0], rows[0][-1]), (5, True))
        self.assertIn(('DELETE FROM aligned_markers WHERE genome_id = ANY(%s)', ([5],)), statements)

    def test_rewritten_deflines_update_the_hash_without_flagging_the_genome(self):
        _, values = self.apply({'GCA_000000001.1': STATUS_SEQUENCES_UNCHANGED},
                               [self.row('GCA_000000001.1', id=5)])

        sql, rows = [v for v in values if 'v.has_changed' in v[0]][0]
        self.assertEqual((rows[0][0], rows[0][2], rows[0][-1]),
                         (5, content_sha1(b'>c1\nACGT\n'), False))

    def test_a_genome_not_added_for_want_of_a_genomic_fasta_is_not_written(self):
        os.makedirs(self.genome_dir('GCA_000000001.1'))
        manager = M.DatabaseManager('host', 'user', 'pw', 'db', cpus=1)
        decisions = list(self.plan({'GCA_000000001.1': STATUS_NEW},
                                   no_genomic_fasta={'GCA_000000001.1'}).values())
        decisions = manager.finish_decisions(decisions, manager.hash_files(decisions))

        cur, values = FakeCursor(), []
        with mock.patch.object(M, 'execute_values',
                               lambda cur, sql, rows, **kw: values.append((sql, rows))):
            manager.apply(cur, decisions, {}, SOURCES, self.DATE)
        self.assertEqual(len(cur.statements), 1)
        self.assertEqual(values, [])

    def test_an_unchanged_genome_is_not_written(self):
        statements, values = self.apply(
            {'GCA_000000001.1': STATUS_FASTA_UNCHANGED},
            [self.row('GCA_000000001.1', genes_hash=content_sha1(b'>g1\nMK\n'))])

        self.assertEqual(len(statements), 1)
        self.assertEqual(values, [])


# ------------------------------------------------------------ stopped part way

class ConnectingToTheDatabase(unittest.TestCase):
    def test_the_server_aborts_a_transaction_its_client_abandoned(self):
        # a machine reset mid-transaction leaves it idle; without the timeout the
        # server keeps its locks until TCP gives up, and a run started again waits
        with mock.patch.object(M.psycopg2, 'connect') as connect:
            M.DatabaseManager('host', 'user', 'pw', 'db').connect()

        kwargs = connect.call_args.kwargs
        self.assertIn('idle_in_transaction_session_timeout=' + M.IDLE_IN_TRANSACTION_TIMEOUT,
                      kwargs['options'])
        self.assertEqual(kwargs['keepalives'], 1)
        self.assertEqual(kwargs['password'], 'pw')


class RunningTheUpdateAgain(TempDirCase):
    """A second run of an update that committed does what the first did, and no more."""

    DATE = datetime(2026, 10, 6)
    GENOMIC, PROTEINS = b'>c1\nACGT\n', b'>g1\nMK\n'

    def held(self, accession, id=1, date_added=datetime(2024, 9, 14)):
        """A genome the database holds as an earlier run of this update left it."""
        return self.row(accession, id=id, fasta_hash=content_sha1(self.GENOMIC),
                        genes_hash=content_sha1(self.PROTEINS))._replace(date_added=date_added)

    def apply(self, outcomes, database):
        for accession, outcome in outcomes.items():
            if outcome in UG.STATUS_IN_RELEASE:
                self.make_genome(accession, self.GENOMIC, self.PROTEINS)
        manager = M.DatabaseManager('host', 'user', 'pw', 'db', cpus=1)
        decisions = list(self.plan(outcomes, database).values())
        hashes = manager.hash_files(decisions)
        decisions = manager.finish_decisions(decisions, hashes)
        cur = FakeCursor()
        with mock.patch.object(M, 'execute_values', lambda *a, **kw: None):
            manager.apply(cur, decisions, hashes, SOURCES, self.DATE)
        return {d.accession: d for d in decisions}, cur.statements

    def flagged(self, statements):
        for sql, params in statements:
            if sql == 'UPDATE genomes SET has_changed = TRUE WHERE id = ANY(%s)':
                return set(params[0])
        return set()

    def test_a_new_genome_an_earlier_run_added_is_already_updated_and_keeps_its_markers(self):
        decisions, statements = self.apply({'GCA_000000001.1': STATUS_NEW},
                                           [self.held('GCA_000000001.1', id=9, date_added=self.DATE)])

        self.assertEqual(decisions['GCA_000000001.1'].action, M.ACTION_ALREADY_UPDATED)
        self.assertFalse([s for s in statements if 'aligned_markers' in s[0]])
        self.assertEqual(self.flagged(statements), {9})

    def test_a_new_version_an_earlier_run_put_in_place_is_already_updated(self):
        decisions, _ = self.apply({'GCA_000000001.3': STATUS_NEW, 'GCA_000000001.2': STATUS_REMOVED},
                                  [self.held('GCA_000000001.3', id=7, date_added=self.DATE)])

        self.assertEqual(decisions['GCA_000000001.3'].action, M.ACTION_ALREADY_UPDATED)
        self.assertEqual(decisions['GCA_000000001.2'].action, M.ACTION_NOT_IN_DATABASE)

    def test_a_genome_whose_sequences_really_changed_is_still_said_to_have(self):
        # the database holds other files for it: this run is the first to see them
        decisions, statements = self.apply({'GCA_000000001.1': STATUS_FASTA_CHANGED},
                                           [self.row('GCA_000000001.1', id=5)])

        self.assertEqual(decisions['GCA_000000001.1'].action, M.ACTION_SEQUENCES_CHANGED)
        self.assertIn(('DELETE FROM aligned_markers WHERE genome_id = ANY(%s)', ([5],)), statements)

    def test_a_genome_an_earlier_run_added_from_the_release_stays_flagged(self):
        # one of r237's 15: in the release, missing from the database, added with
        # this update's date and so has_changed on every run of it
        decisions, statements = self.apply({'GCA_000000001.1': STATUS_FASTA_UNCHANGED},
                                           [self.held('GCA_000000001.1', id=3, date_added=self.DATE)])

        self.assertEqual(decisions['GCA_000000001.1'].action, M.ACTION_UNCHANGED)
        self.assertEqual(self.flagged(statements), {3})

    def test_a_genome_an_earlier_release_added_is_not_flagged(self):
        _, statements = self.apply({'GCA_000000001.1': STATUS_FASTA_UNCHANGED},
                                   [self.held('GCA_000000001.1', id=3)])

        self.assertEqual(self.flagged(statements), set())


class TheListsAffected(TempDirCase):
    """Which curated lists a deleted genome was in, kept across runs."""

    def write(self, run_started, dry_run):
        manager = M.DatabaseManager('host', 'user', 'pw', 'db')
        decisions = [M.Decision('GCF_000806395.1', STATUS_REMOVED, M.ACTION_DELETED,
                                self.row('GCF_000806395.1', id=11))]
        manager.write_lists_affected(FakeCursor([(11, 1014, 'Genomes for MS')]),
                                     decisions, self.dir, run_started, dry_run)

    def rows(self):
        with open(os.path.join(self.dir, M.LISTS_AFFECTED_NAME)) as handle:
            return [line.rstrip('\n').split('\t') for line in handle]

    def test_each_run_appends_its_rows_with_when_it_began_and_whether_it_was_dry(self):
        self.write('2026-10-06T09:00:00', True)
        self.write('2026-10-07T09:00:00', False)

        self.assertEqual(self.rows(), [
            list(M.LISTS_AFFECTED_HEADER),
            ['GCF_000806395.1', '1014', 'Genomes for MS', '2026-10-06T09:00:00', 'True'],
            ['GCF_000806395.1', '1014', 'Genomes for MS', '2026-10-07T09:00:00', 'False']])

    def test_a_file_of_other_columns_is_refused_before_anything_is_hashed(self):
        with open(os.path.join(self.dir, M.LISTS_AFFECTED_NAME), 'w') as handle:
            handle.write('genome_id\tlist_id\tlist_name\n')

        with self.assertRaises(M.UpdateDbError):
            M.DatabaseManager('host', 'user', 'pw', 'db').check_lists_affected(self.dir)


class TheHashCache(TempDirCase):
    """The hashes one run made, read back by the next."""

    def setUp(self):
        super().setUp()
        self.cache_file = os.path.join(self.dir, M.HASH_CACHE_NAME)
        self.make_genome('GCA_000000001.1')
        self.decisions = list(self.plan({'GCA_000000001.1': STATUS_NEW}).values())
        self.manager = M.DatabaseManager('host', 'user', 'pw', 'db', cpus=1)

    def hashes(self):
        return self.manager.hash_files(self.decisions, self.cache_file)['GCA_000000001.1']

    def cached(self, digest):
        """Make every hash of the cache say `digest`, as if an earlier run made it."""
        cache = M.HashCache(self.cache_file)
        cache.load()
        cache.entries = {path: (size, mtime, digest) for path, (size, mtime, _) in cache.entries.items()}
        cache.rewrite()

    def test_the_hashes_made_are_written_to_the_cache(self):
        self.hashes()

        cache = M.HashCache(self.cache_file)
        self.assertFalse(cache.load())
        self.assertEqual(len(cache.entries), 2)

    def test_a_run_started_again_takes_the_hashes_of_the_last(self):
        self.hashes()
        self.cached('from an earlier run')

        self.assertEqual(self.hashes(), ('from an earlier run', 'from an earlier run'))

    def test_a_file_rewritten_since_is_hashed_again(self):
        self.hashes()
        self.cached('from an earlier run')
        genomic = UG.genomic_fasta(self.genome_dir('GCA_000000001.1'))
        with gzip.open(genomic, 'wb') as handle:
            handle.write(b'>c1\nTTTT\n')
        os.utime(genomic, ns=(1, 1))

        self.assertEqual(self.hashes(), (content_sha1(b'>c1\nTTTT\n'), 'from an earlier run'))

    def test_a_cache_a_crash_cut_short_keeps_what_was_flushed_and_is_made_whole(self):
        # as a run reset mid-hashing leaves it: flushed, and no gzip trailer
        with M.HashCache(self.cache_file) as cache:
            cache.add('/a', 1, 1, 'aaa')
            cache.flush()
            with open(self.cache_file, 'rb') as handle:
                cut_short = handle.read()
        with open(self.cache_file, 'wb') as handle:
            handle.write(cut_short)

        cache = M.HashCache(self.cache_file)
        self.assertTrue(cache.load())
        self.assertEqual(cache.entries, {'/a': (1, 1, 'aaa')})
        # written again whole, so what the next run appends can be read
        self.assertFalse(M.HashCache(self.cache_file).load())


if __name__ == '__main__':
    unittest.main()
