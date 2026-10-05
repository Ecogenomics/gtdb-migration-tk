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

    def plan(self, outcomes, database=(), rehash_all=False):
        dirs = {acc: self.genome_dir(acc) for acc, outcome in outcomes.items()
                if outcome in UG.STATUS_IN_RELEASE}
        decisions = M.plan_update(outcomes, dirs, {r.name: r for r in database}, rehash_all)
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

    def test_a_genome_of_the_release_with_no_genomic_fasta_is_refused(self):
        os.makedirs(self.genome_dir('GCA_000000001.1'))
        with self.assertRaises(M.UpdateDbError):
            self.finish({'GCA_000000001.1': STATUS_NEW})


# ------------------------------------------------------------ the statements

class FakeCursor(object):
    def __init__(self):
        self.statements = []

    def execute(self, sql, params=None):
        self.statements.append((' '.join(sql.split()), params))


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

    def test_an_unchanged_genome_is_not_written(self):
        statements, values = self.apply(
            {'GCA_000000001.1': STATUS_FASTA_UNCHANGED},
            [self.row('GCA_000000001.1', genes_hash=content_sha1(b'>g1\nMK\n'))])

        self.assertEqual(len(statements), 1)
        self.assertEqual(values, [])


if __name__ == '__main__':
    unittest.main()
