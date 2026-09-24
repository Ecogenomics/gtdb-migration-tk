#!/usr/bin/env python3
"""Offline unit tests for trans_table.py -- gTranslate itself is never run.

What is tested here is everything decided before the subprocess starts and
everything made of what comes back: which file of a genome directory is handed
over, which genomes are left out, the command line built from the options, and
the comparison, the conflicts and the CheckM2 estimates written afterwards.
Running gTranslate is gTranslate's business; getting the wrong genomes to it, or
the right ones under the wrong names, is this module's.

The batches those genomes arrive in are not tested here. Cutting a release into
batches, claiming one, and the canaries and leases that let several machines share
an output directory belong to batching.py, and are in tests/test_batching.py --
under a layout naming no real command, so that they stay tests of the mechanism.
"""

import gzip
import os
import shutil
import tempfile
import unittest
from unittest import mock

from gtdb_migration_tk import trans_table as G


# What the stubbed gTranslate and CheckM2 say they are. Neither is installed for
# the tests, and every GTranslate() asks both, so the question is answered here
# for the whole module rather than by each test that builds one.
VERSIONS = {G.GTRANSLATE_BIN: 'gtranslate: version 0.0.4', G.CHECKM2_BIN: '1.1.0'}


def setUpModule():
    global _record_program_version
    _record_program_version = G.record_program_version
    G.record_program_version = VERSIONS.__getitem__


def tearDownModule():
    G.record_program_version = _record_program_version


def conflict_row(**fields):
    """A conflict row of whatever width CONFLICT_HEADER currently is.

    Built by name so that a column added to the file does not silently turn every
    hand-written fixture into a row too short to be read back.
    """

    row = {'genome_id': 'GCA_1.1', 'gtranslate_tt': '25', 'ncbi_tt': '11',
           'checkm_tt': '4', 'checkm_conflict': 'False',
           'coding_density_4': '90.1', 'coding_density_11': '64.2',
           'gc_percent': '35.33882', 'n50': '5269725', 'genome_size': '5269725',
           'ncbi_taxonomy': 'd__Bacteria'}
    row.update(fields)

    return [row[column] for column in G.CONFLICT_HEADER]


class TempDirCase(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='trans_table_test.')

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def genome_dir(self, assembly, fasta=True, empty=False):
        """A genome directory as a release holds it, named for its assembly."""
        path = os.path.join(self.dir, assembly)
        os.makedirs(path)
        if fasta:
            with open(os.path.join(path, assembly + '_genomic.fna.gz'), 'wb') as handle:
                if not empty:
                    handle.write(b'>contig\nACGT\n')
        return path


# ------------------------------------------------------------- naming the FASTA

# ------------------------------------------------------------- reading the release

# ------------------------------------------------------------- what is asked about

class CheckBatchFastasTests(TempDirCase):
    """The check belongs to the batch that is about to run, not to the plan."""

    def batch(self, *genomes):
        """A batch directory holding the batchfile the plan would have cut."""
        batch_dir = os.path.join(self.dir, 'batch_000001')
        os.makedirs(batch_dir)
        G.write_batchfile([(G.genomic_fasta(path), accession)
                           for accession, path in genomes],
                          os.path.join(batch_dir, G.BATCHFILE_NAME), compress=True)
        return batch_dir

    def test_gtranslate_is_handed_a_plain_copy_and_never_the_plan(self):
        """gTranslate reads a batchfile with a plain open(), and the plan is
        gzipped; handing it the plan would end every batch before it started."""
        batch_dir = self.batch(('GCF_1.1', self.genome_dir('GCF_1.1_ASM1')),
                               ('GCF_2.1', self.genome_dir('GCF_2.1_ASM2')))

        batchfile, present, missing = G.check_batch_fastas(batch_dir)

        self.assertEqual(batchfile,
                         os.path.join(batch_dir, G.PRESENT_BATCHFILE_NAME))
        self.assertEqual(len(present), 2)
        self.assertEqual(missing, [])
        with open(batchfile, 'rb') as handle:
            self.assertNotEqual(handle.read(2), b'\x1f\x8b')

    def test_the_plan_itself_is_gzipped(self):
        batch_dir = self.batch(('GCF_1.1', self.genome_dir('GCF_1.1_ASM1')))
        with open(os.path.join(batch_dir, G.BATCHFILE_NAME), 'rb') as handle:
            self.assertEqual(handle.read(2), b'\x1f\x8b')

    def test_nothing_is_recorded_as_missing_where_nothing_is(self):
        batch_dir = self.batch(('GCF_1.1', self.genome_dir('GCF_1.1_ASM1')))
        G.check_batch_fastas(batch_dir)
        self.assertFalse(os.path.exists(os.path.join(batch_dir, G.MISSING_NAME)))

    def test_a_missing_genome_is_left_out_of_what_gtranslate_is_handed(self):
        batch_dir = self.batch(('GCF_1.1', self.genome_dir('GCF_1.1_ASM1')),
                               ('GCA_2.1', self.genome_dir('GCA_2.1_ASM2', fasta=False)),
                               ('GCF_3.1', self.genome_dir('GCF_3.1_ASM3')))

        batchfile, present, missing = G.check_batch_fastas(batch_dir)

        self.assertEqual(batchfile,
                         os.path.join(batch_dir, G.PRESENT_BATCHFILE_NAME))
        self.assertEqual([accession for _, accession in G.read_batchfile(batchfile)],
                         ['GCF_1.1', 'GCF_3.1'])
        self.assertEqual(missing, ['GCA_2.1'])

    def test_the_genome_left_out_is_named_in_the_batch_directory(self):
        batch_dir = self.batch(('GCF_1.1', self.genome_dir('GCF_1.1_ASM1')),
                               ('GCA_2.1', self.genome_dir('GCA_2.1_ASM2', fasta=False)))

        G.check_batch_fastas(batch_dir)

        with open(os.path.join(batch_dir, G.MISSING_NAME)) as handle:
            self.assertEqual(handle.read().split(), ['GCA_2.1'])

    def test_the_batchfile_the_plan_cut_is_left_as_the_record_of_the_batch(self):
        # the comparison reads it to know which genomes the batch was, and a
        # filtered copy must not become that record
        batch_dir = self.batch(('GCF_1.1', self.genome_dir('GCF_1.1_ASM1')),
                               ('GCA_2.1', self.genome_dir('GCA_2.1_ASM2', fasta=False)))

        G.check_batch_fastas(batch_dir)

        self.assertEqual([accession for _, accession in G.read_batchfile(
            G.batchfile_path(batch_dir, G.LAYOUT))], ['GCF_1.1', 'GCA_2.1'])


# ------------------------------------------------------------- the command line

class DetectTableCommandTests(unittest.TestCase):
    """An option not given is left off, so gTranslate's defaults stay its own."""

    def test_batchfile_is_used_rather_than_genome_dir(self):
        cmd = G.detect_table_command('/out/batch.tsv', '/out')
        self.assertIn('--batchfile', cmd)
        self.assertNotIn('--genome_dir', cmd)

    def test_cpus_is_passed_as_a_string(self):
        cmd = G.detect_table_command('/out/batch.tsv', '/out', cpus=16)
        self.assertEqual(cmd[cmd.index('--cpus') + 1], '16')

    def test_unset_options_are_absent(self):
        cmd = G.detect_table_command('/out/batch.tsv', '/out')
        for flag in ('--tmpdir', '--prefix', '--custom_model_path',
                     '--force', '--keep_called_genes'):
            self.assertNotIn(flag, cmd)

    def test_flags_appear_only_when_asked_for(self):
        cmd = G.detect_table_command('/out/batch.tsv', '/out',
                                     tmp_dir='/scratch', force=True,
                                     keep_called_genes=True, prefix='r232',
                                     custom_model_path='/srv/db/gtranslate/models')
        self.assertEqual(cmd[cmd.index('--tmpdir') + 1], '/scratch')
        self.assertEqual(cmd[cmd.index('--prefix') + 1], 'r232')
        self.assertEqual(cmd[cmd.index('--custom_model_path') + 1],
                         '/srv/db/gtranslate/models')
        self.assertIn('--force', cmd)
        self.assertIn('--keep_called_genes', cmd)

    def test_subcommand_is_detect_table(self):
        cmd = G.detect_table_command('/out/batch.tsv', '/out')
        self.assertEqual(cmd[:2], [G.GTRANSLATE_BIN, 'detect_table'])


# ------------------------------------------------------------- the manager itself

class ManagerTests(unittest.TestCase):
    """gtranslate, checkm2 and prodigal are checked for when the manager is built."""

    def setUp(self):
        # the real check exits the process, and neither tool is wanted offline
        self._check, G.check_dependencies = G.check_dependencies, lambda *a, **k: True

    def tearDown(self):
        G.check_dependencies = self._check

    def test_every_third_party_tool_is_checked_for(self):
        """gTranslate calls Prodigal and its package does not depend on it, and
        CheckM2 is run from an environment of its own, so neither arrives with
        this package and a run that cannot find one should say so before it
        spends hours finding out."""
        asked = []
        G.check_dependencies = lambda programs, *a, **k: asked.extend(programs)
        G.GTranslate()
        self.assertEqual(sorted(asked), ['checkm2', 'gtranslate', 'prodigal'])

    def test_batch_size_default_matches_the_command_line(self):
        """The two defaults are written out separately and must not drift apart."""
        self.assertEqual(G.GTranslate().batch_size, G.DEFAULT_BATCH_SIZE)
        self.assertEqual(G.DEFAULT_BATCH_SIZE, 10000)


if __name__ == '__main__':
    unittest.main()


# ------------------------------------------------------------- planning the batches

class PlanBatchesTests(TempDirCase):
    """The plan is what several machines agree on, so it must not move."""

    def rows(self, *accessions):
        return [('/rel/x/{}_ASM1_genomic.fna.gz'.format(a), a) for a in accessions]

    def test_release_is_cut_into_batches_of_the_given_size(self):
        batches = G.create_batches(self.rows(*['GCF_{}.1'.format(i) for i in range(5)]),
                                   2, self.dir, G.LAYOUT)
        self.assertEqual(len(batches), 3)
        self.assertEqual(len(G.read_batchfile(
            os.path.join(batches[0], G.BATCHFILE_NAME))), 2)
        self.assertEqual(len(G.read_batchfile(
            os.path.join(batches[-1], G.BATCHFILE_NAME))), 1)

    def test_batches_are_numbered_in_order(self):
        batches = G.create_batches(self.rows('GCF_1.1', 'GCF_2.1'), 1, self.dir, G.LAYOUT)
        self.assertEqual([os.path.basename(b) for b in batches],
                         ['batch_000001', 'batch_000002'])

    def test_an_existing_plan_is_found_and_reused(self):
        G.create_batches(self.rows('GCF_1.1', 'GCF_2.1'), 1, self.dir, G.LAYOUT)
        self.assertEqual(len(G.batch_dir_names(self.dir, G.LAYOUT)), 2)

    def test_a_directory_without_a_batchfile_is_not_a_batch(self):
        os.makedirs(os.path.join(self.dir, 'batch_000001'))
        self.assertEqual(G.batch_dir_names(self.dir, G.LAYOUT), [])

    def test_batchfile_survives_a_round_trip(self):
        batches = G.create_batches(self.rows('GCF_1.1'), 10, self.dir, G.LAYOUT)
        self.assertEqual(G.read_batchfile(os.path.join(batches[0], G.BATCHFILE_NAME)),
                         self.rows('GCF_1.1'))

    def test_the_plan_is_cut_from_the_genome_dirs_file_without_touching_a_genome(self):
        # the plan must be written at once rather than after hours of stat calls
        # over NFS, so a genome whose FASTA is not there is still planned; whether
        # the file exists is the running batch's question, not the plan's
        out_dir = os.path.join(self.dir, 'out')
        os.makedirs(out_dir)
        genome_dirs = os.path.join(self.dir, 'genome_dirs.tsv')
        with open(genome_dirs, 'w') as handle:
            handle.write('GCF_2.1\t{}\tG2\n'.format(
                self.genome_dir('GCF_2.1_ASM2', fasta=False)))
            handle.write('GCF_1.1\t/gone/GCF_1.1_ASM1\tG1\n')

        with mock.patch.object(G, 'check_dependencies', lambda *a, **k: True):
            batches = G.GTranslate(batch_size=10).plan_batches(genome_dirs, out_dir)

        # sorted by accession, and both of them there
        self.assertEqual([accession for _, accession in G.read_batchfile(
            os.path.join(batches[0], G.BATCHFILE_NAME))], ['GCF_1.1', 'GCF_2.1'])

    def test_a_genome_dirs_file_naming_nothing_is_an_error(self):
        out_dir = os.path.join(self.dir, 'out')
        os.makedirs(out_dir)
        empty = os.path.join(self.dir, 'empty.tsv')
        open(empty, 'w').close()

        with mock.patch.object(G, 'check_dependencies', lambda *a, **k: True):
            with self.assertRaises(RuntimeError):
                G.GTranslate().plan_batches(empty, out_dir)


# ------------------------------------------------------------- the canaries

# ------------------------------------------------------------- the lease

# ------------------------------------------------------------- work already done

class PredictedTests(TempDirCase):
    """gTranslate is the hours of a batch; it is not run twice over one batch."""

    def batch(self, summary=True, canary=True):
        path = os.path.join(self.dir, 'batch_000001')
        os.makedirs(path)
        if summary:
            open(os.path.join(path, G.summary_name()), 'w').close()
        if canary:
            G.mark_predicted(path, genomes=10)
        return path

    def test_a_batch_gtranslate_finished_is_not_predicted_again(self):
        self.assertTrue(G.already_predicted(self.batch(), G.summary_name()))

    def test_the_canary_alone_is_not_taken_for_results(self):
        """A batch whose summary is not there has nothing for the comparison to read."""
        self.assertFalse(G.already_predicted(self.batch(summary=False),
                                             G.summary_name()))

    def test_results_without_the_canary_are_not_assumed_final(self):
        """gTranslate writes as it goes; only its exit code says the batch is done."""
        self.assertFalse(G.already_predicted(self.batch(canary=False),
                                             G.summary_name()))

    def test_what_was_predicted_is_recorded_with_the_canary(self):
        batch = self.batch()
        self.assertEqual(
            G.read_canary(os.path.join(batch, G.PREDICTED_CANARY))['genomes'], '10')

    def test_a_predicted_batch_is_still_unfinished(self):
        """Nothing is finished until it has been compared and said so."""
        self.assertEqual(G.batch_state(self.batch()), G.STATE_PENDING)


# ------------------------------------------------------------- genomes dropped

class NoPredictionTests(TempDirCase):
    """A genome --force drops is named, not discovered by prodigal months later."""

    def batch(self):
        path = os.path.join(self.dir, 'batch_000001')
        os.makedirs(path)
        return path

    def test_a_genome_with_no_prediction_is_named(self):
        batch = self.batch()
        missing = G.report_no_prediction(batch, ['G1', 'G2', 'G3'], ['G1', 'G3'])
        self.assertEqual(missing, ['G2'])
        self.assertEqual(
            open(os.path.join(batch, G.NO_PREDICTION_NAME)).read().split(), ['G2'])

    def test_nothing_is_written_when_every_genome_was_predicted(self):
        """The file being there at all says a batch lost genomes."""
        batch = self.batch()
        self.assertEqual(G.report_no_prediction(batch, ['G1'], ['G1']), [])
        self.assertFalse(os.path.exists(os.path.join(batch, G.NO_PREDICTION_NAME)))

    def test_the_order_the_genomes_were_given_in_is_kept(self):
        batch = self.batch()
        self.assertEqual(G.report_no_prediction(batch, ['G3', 'G1', 'G2'], []),
                         ['G3', 'G1', 'G2'])


# ------------------------------------------------------------- the batch's own log

# ------------------------------------------------------------- the run over a release

class RunTests(TempDirCase):
    """What a run does to a batch it finds part-done, and to one it is stopped in."""

    def setUp(self):
        super().setUp()
        self._check, G.check_dependencies = G.check_dependencies, lambda *a, **k: True
        self.out_dir = os.path.join(self.dir, 'out')
        self.genome_dirs = os.path.join(self.dir, 'genome_dirs.tsv')
        with open(self.genome_dirs, 'w') as handle:
            handle.write('GCA_000001.1\t{}\tG000001\n'.format(
                self.genome_dir('GCA_000001.1_ASM1')))
        self.taxonomy = os.path.join(self.dir, 'taxonomy.tsv')
        open(self.taxonomy, 'w').close()

    def tearDown(self):
        G.check_dependencies = self._check
        super().tearDown()

    def manager(self):
        return G.GTranslate(batch_size=1)

    def run_one(self, manager):
        return manager.run(self.genome_dirs, self.taxonomy, self.out_dir)

    def comparison(self, manager, batch_dir, taxonomy):
        """Stands in for compare_batch, leaving what the run aggregates."""
        G.write_table([], os.path.join(batch_dir, G.CONFLICT_NAME), G.CONFLICT_HEADER)
        G.write_table([], os.path.join(batch_dir, G.COMPARISON_NAME),
                      header=G.COMPARISON_HEADER, compress=True)
        with open(os.path.join(batch_dir, G.summary_name()), 'w') as handle:
            handle.write('genome_id\n')
        return G.ComparisonCounts(compared=1, conflicts=0, no_ncbi_table=0)

    def test_a_batch_gtranslate_already_finished_is_only_compared(self):
        """The prediction is hours and the comparison is seconds; a lost machine
        between the two must not cost the hours."""
        manager = self.manager()
        manager.plan_batches(self.genome_dirs, self.out_dir)
        batch = os.path.join(self.out_dir, 'batch_000001')
        open(os.path.join(batch, G.summary_name()), 'w').close()
        G.mark_predicted(batch)

        with mock.patch.object(G.GTranslate, 'run_gtranslate') as predict, \
                mock.patch.object(G.GTranslate, 'compare_batch', autospec=True,
                                  side_effect=self.comparison) as compare:
            self.run_one(manager)

        self.assertFalse(predict.called)
        self.assertTrue(compare.called)
        self.assertEqual(G.batch_state(batch), G.STATE_SUCCESS)

    def test_a_batch_nothing_has_been_done_to_is_predicted(self):
        manager = self.manager()
        with mock.patch.object(G.GTranslate, 'run_gtranslate') as predict, \
                mock.patch.object(G.GTranslate, 'compare_batch', autospec=True,
                                  side_effect=self.comparison):
            self.run_one(manager)

        self.assertTrue(predict.called)

    def test_a_run_interrupted_in_a_batch_gives_the_batch_back(self):
        """Ctrl-C is not a failure of the batch, and the next run should not wait
        out a lease on a machine that has already stopped."""
        manager = self.manager()
        with mock.patch.object(G.GTranslate, 'run_gtranslate',
                               side_effect=KeyboardInterrupt):
            self.assertRaises(KeyboardInterrupt, self.run_one, manager)

        batch = os.path.join(self.out_dir, 'batch_000001')
        self.assertEqual(G.batch_state(batch), G.STATE_PENDING)

    def test_what_happened_to_a_batch_is_written_in_the_batch(self):
        manager = self.manager()
        with mock.patch.object(G.GTranslate, 'run_gtranslate',
                               side_effect=RuntimeError('gtranslate returned exit code 1.')):
            self.run_one_expecting_failure(manager)

        batch = os.path.join(self.out_dir, 'batch_000001')
        log = open(os.path.join(batch, G.BATCH_LOG_NAME)).read()
        self.assertIn('exit code 1', log)
        self.assertEqual(G.batch_state(batch), G.STATE_FAILED)

    def run_one_expecting_failure(self, manager):
        """A run ending with a failed batch says so, having done the others."""
        self.assertFalse(self.run_one(manager))

    def test_a_failed_batch_is_reported_rather_than_raised(self):
        """The batch says why it failed and is retried by the next run; a
        traceback out of days of work over five machines adds nothing."""
        manager = self.manager()
        with mock.patch.object(G.GTranslate, 'run_gtranslate',
                               side_effect=RuntimeError('gtranslate returned exit code 1.')):
            with self.assertLogs(manager.logger, level='ERROR') as logged:
                finished = self.run_one(manager)

        self.assertFalse(finished)
        self.assertTrue(any('1 batch(es) failed' in message
                            for message in logged.output))


# ------------------------------------------------------------- the comparison

class ComparisonTests(TempDirCase):
    """The comparison holds every genome the two both called; the conflicts are
    that file filtered."""

    def genome_dir(self, accession, ncbi_table=None):
        assembly = '{}_ASM1'.format(accession)
        path = os.path.join(self.dir, assembly)
        os.makedirs(path)
        if ncbi_table is not None:
            with gzip.open(os.path.join(path, assembly + '_genomic.gff.gz'), 'wt') as handle:
                handle.write('##gff-version 3\n')
                handle.write('c\tRefSeq\tCDS\t1\t9\t.\t+\t0\t'
                             'ID=cds1;product=x;transl_table={}\n'.format(ncbi_table))
        return path

    def summary(self, *rows):
        path = os.path.join(self.dir, 'gtranslate.translation_table_summary.tsv')
        with open(path, 'w') as handle:
            handle.write('user_genome\tbest_tln_table\tcoding_density_4\t'
                         'coding_density_11\tconfidence\n')
            for accession, table in rows:
                handle.write('{}\t{}\t90.1\t64.2\t1.0\n'.format(accession, table))
        return path

    def test_a_genome_agreeing_with_ncbi_is_in_the_comparison_not_the_conflicts(self):
        """The rate the two differ at, and whether the genomes they differ about
        are unlike the ones they agree about, need the agreements present. The
        conflict file is the few hundred rows worth looking at."""
        path = self.genome_dir('GCF_1.1', ncbi_table=11)
        rows, compared, _ = G.comparison_rows(
            G.read_translation_table_summary(self.summary(('GCF_1.1', '11'))),
            {'GCF_1.1': path}, {})
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0][G.COMPARISON_HEADER.index('ncbi_conflict')], 'False')
        self.assertEqual(G.conflicts_from_comparison(rows), [])
        self.assertEqual(compared, 1)

    def test_a_genome_the_two_differ_about_is_marked_and_kept(self):
        path = self.genome_dir('GCF_1.1', ncbi_table=11)
        rows, compared, _ = G.comparison_rows(
            G.read_translation_table_summary(self.summary(('GCF_1.1', '4'))),
            {'GCF_1.1': path}, {})
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0][G.COMPARISON_HEADER.index('gtranslate_tt')], '4')
        self.assertEqual(rows[0][G.COMPARISON_HEADER.index('ncbi_tt')], '11')
        self.assertEqual(rows[0][G.COMPARISON_HEADER.index('ncbi_conflict')], 'True')
        self.assertEqual(len(G.conflicts_from_comparison(rows)), 1)
        self.assertEqual(compared, 1)

    def test_the_conflict_file_says_nothing_about_whether_a_row_is_a_conflict(self):
        """Every row of it is one, so the column would say True on every line.
        It is in the comparison, where it tells the rows apart."""
        self.assertNotIn('ncbi_conflict', G.CONFLICT_HEADER)
        self.assertIn('ncbi_conflict', G.COMPARISON_HEADER)

        path = os.path.join(self.dir, 'conflicts.tsv')
        G.write_table([conflict_row()], path, G.CONFLICT_HEADER)
        with open(path) as handle:
            header, row = handle.read().splitlines()
        self.assertEqual(header.split('\t'), list(G.CONFLICT_HEADER))
        self.assertEqual(len(row.split('\t')), len(G.CONFLICT_HEADER))

    def test_a_batch_with_nothing_to_report_still_writes_the_file(self):
        """The release file is every batch's concatenated, so a batch that
        conflicted nowhere has to leave a header behind."""
        path = os.path.join(self.dir, 'none.tsv')
        G.write_table([], path, G.CONFLICT_HEADER)
        self.assertEqual(open(path).read().splitlines(), ['\t'.join(G.CONFLICT_HEADER)])

    def test_genome_ncbi_declares_no_table_for_is_left_out(self):
        path = self.genome_dir('GCF_1.1')
        rows, compared, no_table = G.comparison_rows(
            G.read_translation_table_summary(self.summary(('GCF_1.1', '11'))), {'GCF_1.1': path}, {})
        self.assertEqual(rows, [])
        self.assertEqual(compared, 0)
        self.assertEqual(no_table, 1)

    def test_coding_densities_and_lineage_are_carried(self):
        path = self.genome_dir('GCF_1.1', ncbi_table=11)
        rows, _, _ = G.comparison_rows(
            G.read_translation_table_summary(self.summary(('GCF_1.1', '4'))),
            {'GCF_1.1': path},
            {'GCF_1.1': 'd__Bacteria;p__Pseudomonadota;s__Escherichia coli'})
        self.assertEqual(rows[0][G.COMPARISON_HEADER.index('coding_density_4')], '90.1')
        self.assertEqual(rows[0][G.COMPARISON_HEADER.index('coding_density_11')], '64.2')
        self.assertIn('d__Bacteria', rows[0][G.COMPARISON_HEADER.index('ncbi_taxonomy')])

    def test_lineage_is_found_through_the_canonical_accession(self):
        """A GenBank genome takes the lineage held against its RefSeq counterpart."""
        path = self.genome_dir('GCA_005435135.1', ncbi_table=11)
        rows, _, _ = G.comparison_rows(
            G.read_translation_table_summary(self.summary(('GCA_005435135.1', '4'))),
            {'GCA_005435135.1': path},
            G.read_taxonomy(self.write_taxonomy('GCF_005435135.1\td__Bacteria;s__X\n')))
        self.assertEqual(rows[0][G.COMPARISON_HEADER.index('ncbi_taxonomy')], 'd__Bacteria;s__X')

    def test_genome_missing_from_the_taxonomy_is_still_reported(self):
        path = self.genome_dir('GCF_1.1', ncbi_table=4)
        rows, _, _ = G.comparison_rows(
            G.read_translation_table_summary(self.summary(('GCF_1.1', '11'))),
            {'GCF_1.1': path}, {})
        self.assertEqual(rows[0][G.COMPARISON_HEADER.index('ncbi_taxonomy')], 'na')

    def write_taxonomy(self, text):
        path = os.path.join(self.dir, 'taxonomy.tsv')
        with open(path, 'w') as handle:
            handle.write(text)
        return path

    def test_checkm_table_is_reported_beside_the_others(self):
        """The density rule alone cannot express 25, which is the point of it."""
        path = self.genome_dir('GCF_1.1', ncbi_table=11)
        rows, _, _ = G.comparison_rows(
            G.read_translation_table_summary(self.summary(('GCF_1.1', '25'))),
            {'GCF_1.1': path}, {})
        row = rows[0]
        self.assertEqual(row[G.COMPARISON_HEADER.index('gtranslate_tt')], '25')
        self.assertEqual(row[G.COMPARISON_HEADER.index('checkm_tt')], '4')

    def test_checkm_conflict_is_true_where_the_two_disagree_about_recoding(self):
        """11 against 4, and 4 or 25 against 11, are calls of different kinds:
        genes called under the wrong one break at every TGA."""
        for predicted, checkm in (('11', 4), ('4', 11), ('25', 11)):
            with self.subTest(gtranslate=predicted, checkm=checkm):
                self.assertTrue(G.checkm_conflict(predicted, checkm))

    def test_checkm_conflict_is_false_where_the_rule_cannot_say_otherwise(self):
        """The density rule picks between 4 and 11 alone, so 4 is the only
        recoding it can return; 25 against 4 is as close to agreement as it
        gets, which is the whole reason gTranslate is run."""
        for predicted, checkm in (('11', 11), ('4', 4), ('25', 4)):
            with self.subTest(gtranslate=predicted, checkm=checkm):
                self.assertFalse(G.checkm_conflict(predicted, checkm))

    def test_checkm_conflict_is_false_where_either_table_is_unknown(self):
        self.assertFalse(G.checkm_conflict('11', None))
        self.assertFalse(G.checkm_conflict('', 4))
        self.assertFalse(G.checkm_conflict('na', 11))

    def test_checkm_conflict_is_reported_beside_the_checkm_table(self):
        """It qualifies checkm_tt, so it is read next to it rather than hunted
        for at the end of the row."""
        # in the conflict file checkm_conflict follows checkm_tt; in the
        # comparison ncbi_conflict comes between them, the two verdicts sitting
        # together because that file is where they differ
        self.assertEqual(G.CONFLICT_HEADER.index('checkm_conflict'),
                         G.CONFLICT_HEADER.index('checkm_tt') + 1)
        self.assertEqual(G.COMPARISON_HEADER.index('checkm_conflict'),
                         G.COMPARISON_HEADER.index('ncbi_conflict') + 1)

        path = self.genome_dir('GCF_1.1', ncbi_table=4)
        rows, _, _ = G.comparison_rows(
            G.read_translation_table_summary(self.summary(('GCF_1.1', '11'))),
            {'GCF_1.1': path}, {})
        self.assertEqual(rows[0][G.COMPARISON_HEADER.index('checkm_tt')], '4')
        self.assertEqual(rows[0][G.COMPARISON_HEADER.index('checkm_conflict')], 'True')

    def test_checkm_conflict_is_false_on_a_row_the_density_rule_agrees_with(self):
        path = self.genome_dir('GCF_1.1', ncbi_table=11)
        rows, _, _ = G.comparison_rows(
            G.read_translation_table_summary(self.summary(('GCF_1.1', '25'))),
            {'GCF_1.1': path}, {})
        self.assertEqual(rows[0][G.COMPARISON_HEADER.index('checkm_tt')], '4')
        self.assertEqual(rows[0][G.COMPARISON_HEADER.index('checkm_conflict')], 'False')

    def test_checkm_table_is_11_where_the_densities_are_close(self):
        path = self.genome_dir('GCF_1.1', ncbi_table=4)
        summary = os.path.join(self.dir, 'close.tsv')
        with open(summary, 'w') as handle:
            handle.write('user_genome\tbest_tln_table\tcoding_density_4\t'
                         'coding_density_11\n')
            handle.write('GCF_1.1\t11\t86.24689\t86.64953\n')
        rows, _, _ = G.comparison_rows(
            G.read_translation_table_summary(summary), {'GCF_1.1': path}, {})
        self.assertEqual(rows[0][G.COMPARISON_HEADER.index('checkm_tt')], '11')


class AggregationTests(TempDirCase):
    """The release file is the whole release or absent, never a part of it."""

    def comparison(self, name, *rows):
        path = os.path.join(self.dir, name)
        G.write_table(rows, path, G.CONFLICT_HEADER)
        return path

    def test_headers_are_not_repeated(self):
        first = self.comparison('a.tsv', ('GCF_1.1', '4', '11', '4', 'True', '', '', 'na'))
        second = self.comparison('b.tsv', ('GCF_2.1', '25', '11', '4', 'False', '', '', 'na'))
        out = os.path.join(self.dir, 'all.tsv')
        self.assertEqual(G.concatenate([first, second], out), 2)
        with open(out) as handle:
            lines = handle.read().splitlines()
        self.assertEqual(len(lines), 3)
        self.assertEqual(lines[0].split('\t')[0], 'genome_id')


# ------------------------------------------- how much of the release NCBI annotates

class DisagreementRateTests(unittest.TestCase):
    """The rate is of the genomes that could be compared, not of the release."""

    def test_the_denominator_is_the_genomes_ncbi_declares_a_table_for(self):
        """Counting the unannotated genomes in would make the number a measure of
        how much of the release NCBI has annotated, not of how often the two
        differ."""
        self.assertAlmostEqual(G.disagreement_rate(1, 4), 25.0)

    def test_a_release_nothing_could_be_compared_in_has_a_rate_of_zero(self):
        """A batch of genomes NCBI has annotated none of divides by nothing, and
        the log line is still written."""
        self.assertEqual(G.disagreement_rate(0, 0), 0.0)


class BatchCountsTests(TempDirCase):
    """The release counts come from the batches' canaries, which is the only place
    a batch run on another machine leaves them."""

    def batch(self, name, **fields):
        path = os.path.join(self.dir, name)
        os.makedirs(path)
        G.finish_batch(path, **fields)
        return path

    def test_the_counts_of_every_batch_are_added_up(self):
        batches = [self.batch('batch_000001', compared=10, no_ncbi_table=2),
                   self.batch('batch_000002', compared=7, no_ncbi_table=3)]
        self.assertEqual(G.batch_counts(batches), (17, 5))

    def test_a_batch_finished_before_the_field_existed_leaves_that_total_unknown(self):
        """r237 was predicted before no_ncbi_table was recorded. Its release line
        still reports the genomes compared and the rate, and says nothing about
        the genomes NCBI declares no table for rather than reporting a total that
        is short by every batch predicted then."""
        batches = [self.batch('batch_000001', compared=10, no_ncbi_table=2),
                   self.batch('batch_000002', compared=7)]
        compared, no_ncbi_table = G.batch_counts(batches)
        self.assertEqual(compared, 17)
        self.assertIsNone(no_ncbi_table)


# --------------------------------------------- the quality of a conflicting genome

class ConflictTableTests(unittest.TestCase):
    """A conflicting genome is asked about under BOTH of the tables in dispute."""

    def row(self, accession, gtranslate_tt, ncbi_tt):
        return conflict_row(genome_id=accession, gtranslate_tt=gtranslate_tt,
                            ncbi_tt=ncbi_tt)

    def test_a_genome_is_run_under_each_of_the_two_tables_it_disputes(self):
        """Completeness under a table is the evidence about that table; one run
        would say how good the genome is, two say which table makes it look like
        a genome at all."""
        tables = G.conflict_tables([self.row('GCA_1.1', '25', '11')])
        self.assertEqual(tables, {25: ['GCA_1.1'], 11: ['GCA_1.1']})

    def test_the_genomes_disputing_one_table_are_gathered_into_one_run(self):
        """CheckM2 searches the whole DIAMOND database once per run whatever the
        run holds, so the cost is the number of runs and barely the number of
        genomes."""
        tables = G.conflict_tables([self.row('GCA_2.1', '4', '11'),
                                    self.row('GCA_1.1', '25', '11')])
        self.assertEqual(tables[11], ['GCA_1.1', 'GCA_2.1'])
        self.assertEqual(sorted(tables), [4, 11, 25])

    def test_a_table_that_is_not_a_number_is_not_run(self):
        """There is nothing to force Prodigal to, and the genome keeps its row."""
        tables = G.conflict_tables([self.row('GCA_1.1', '25', 'na')])
        self.assertEqual(tables, {25: ['GCA_1.1']})


class StageCheckM2InputTests(TempDirCase):
    """CheckM2 names a result for the file it read, so the file is named for the
    accession the conflict row is keyed by."""

    def test_a_genome_is_linked_under_its_accession(self):
        """NCBI names the FASTA for the assembly, so the genome would otherwise
        come back as GCA_000238995.1_ASM23899v1_genomic and join to nothing."""
        path = self.genome_dir('GCA_1.1_ASM1')
        fasta = os.path.join(path, 'GCA_1.1_ASM1_genomic.fna.gz')
        staged = G.stage_checkm2_input(['GCA_1.1'], {'GCA_1.1': fasta},
                                       os.path.join(self.dir, 'input'))
        self.assertEqual([os.path.basename(link) for link in staged],
                         ['GCA_1.1.fna.gz'])
        self.assertEqual(os.path.realpath(staged[0]), os.path.realpath(fasta))

    def test_a_genome_with_no_fasta_is_left_out_rather_than_linked_to_nothing(self):
        """The row stays and reports na; a dangling link would fail the whole run
        of that table."""
        staged = G.stage_checkm2_input(
            ['GCA_1.1'], {'GCA_1.1': os.path.join(self.dir, 'gone.fna.gz')},
            os.path.join(self.dir, 'input'))
        self.assertEqual(staged, [])

    def test_staging_again_replaces_the_link_rather_than_failing(self):
        """A table whose run failed is staged again by the next run."""
        path = self.genome_dir('GCA_1.1_ASM1')
        fasta = os.path.join(path, 'GCA_1.1_ASM1_genomic.fna.gz')
        staging = os.path.join(self.dir, 'input')
        G.stage_checkm2_input(['GCA_1.1'], {'GCA_1.1': fasta}, staging)
        staged = G.stage_checkm2_input(['GCA_1.1'], {'GCA_1.1': fasta}, staging)
        self.assertEqual(len(staged), 1)


class CheckM2CommandTests(unittest.TestCase):
    """The table is forced; letting CheckM2 choose would ask a question that is
    already answered."""

    def command(self, **kwargs):
        options = dict(fastas=['/x/GCA_1.1.fna.gz'], table=25,
                       out_dir='/out/table_25', threads=8)
        options.update(kwargs)
        return G.checkm2_command(**options)

    def test_the_table_is_forced(self):
        """CheckM2 left to itself picks between 4 and 11 by coding density, which
        is the rule checkm_tt already reports and cannot express 25 at all."""
        cmd = self.command()
        self.assertEqual(cmd[cmd.index('--ttable') + 1], '25')

    def test_the_genomes_come_last_because_input_takes_the_rest_of_the_line(self):
        """--input is nargs='+', so anything after it is read as a genome."""
        cmd = self.command(fastas=['/x/a.fna.gz', '/x/b.fna.gz'])
        self.assertEqual(cmd[-3:], ['--input', '/x/a.fna.gz', '/x/b.fna.gz'])

    def test_the_output_directory_is_the_one_the_report_is_read_from(self):
        cmd = self.command()
        self.assertEqual(cmd[cmd.index('--output-directory') + 1], '/out/table_25')


class ReadCheckM2ReportTests(TempDirCase):
    """The report is read by column name, as the NCBI tables are."""

    def report(self, header, *rows):
        path = os.path.join(self.dir, G.CHECKM2_REPORT)
        with open(path, 'w') as handle:
            handle.write('\t'.join(header) + '\n')
            for row in rows:
                handle.write('\t'.join(row) + '\n')
        return path

    def test_the_columns_are_found_by_name_and_not_by_position(self):
        """CheckM2 writes a dozen columns and has added to them between releases;
        the two wanted sit in the middle of the rest."""
        path = self.report(('Name', 'Translation_Table_Used', 'Contamination',
                            'Completeness'),
                           ('GCA_1.1', '25', '0.17', '94.3'))
        self.assertEqual(G.read_checkm2_report(path), {'GCA_1.1': ('94.3', '0.17')})

    def test_a_run_that_wrote_no_report_is_not_an_exception(self):
        """A failed run costs four columns and not the release."""
        self.assertEqual(G.read_checkm2_report(os.path.join(self.dir, 'gone.tsv')), {})

    def test_a_report_with_no_recognisable_header_is_read_as_nothing(self):
        path = self.report(('something', 'else'), ('a', 'b'))
        self.assertEqual(G.read_checkm2_report(path), {})


class AnnotateConflictsTests(unittest.TestCase):
    """Each genome carries the estimate made under each of its two tables."""

    def row(self, accession='GCA_1.1', gtranslate_tt='25', ncbi_tt='11'):
        return conflict_row(genome_id=accession, gtranslate_tt=gtranslate_tt,
                            ncbi_tt=ncbi_tt)

    def annotated(self, rows, quality):
        return G.annotate_conflicts(rows, quality)[0]

    def field(self, row, column):
        return row[G.CONFLICT_HEADER_CHECKM2.index(column)]

    def test_each_estimate_goes_under_the_table_it_was_made_for(self):
        """The columns pair by name with gtranslate_tt and ncbi_tt, and mixing
        them up would reverse what the table says about the conflict."""
        row = self.annotated([self.row()],
                             {25: {'GCA_1.1': ('94.3', '0.17')},
                              11: {'GCA_1.1': ('51.0', '16.4')}})
        self.assertEqual(self.field(row, 'cm2_completeness_gtranslate_tt'), '94.3')
        self.assertEqual(self.field(row, 'cm2_contamination_gtranslate_tt'), '0.17')
        self.assertEqual(self.field(row, 'cm2_completeness_ncbi_tt'), '51.0')
        self.assertEqual(self.field(row, 'cm2_contamination_ncbi_tt'), '16.4')

    def test_the_lineage_stays_the_last_field_of_the_row(self):
        """It is the longest field and the columns of the comparison belong
        together; a row read by eye is unreadable otherwise."""
        row = self.annotated([self.row()], {})
        self.assertEqual(row[-1], 'd__Bacteria')
        self.assertEqual(len(row), len(G.CONFLICT_HEADER_CHECKM2))

    def test_a_genome_checkm2_returned_nothing_for_keeps_its_row(self):
        """The row is the conflict, and the conflict is there whether or not its
        quality could be estimated."""
        row = self.annotated([self.row()], {25: {}, 11: {}})
        self.assertEqual(self.field(row, 'cm2_completeness_gtranslate_tt'), G.NCBI_NA)
        self.assertEqual(self.field(row, 'cm2_completeness_ncbi_tt'), G.NCBI_NA)
        self.assertEqual(row[0], 'GCA_1.1')

    def test_a_genome_estimated_under_one_table_only_keeps_that_one(self):
        """One table's run failing does not cost the other table's answer."""
        row = self.annotated([self.row()], {25: {'GCA_1.1': ('94.3', '0.17')}})
        self.assertEqual(self.field(row, 'cm2_completeness_gtranslate_tt'), '94.3')
        self.assertEqual(self.field(row, 'cm2_completeness_ncbi_tt'), G.NCBI_NA)


class ConflictHeaderTests(unittest.TestCase):
    """The release header is derived from the batch header and must stay so."""

    def test_the_release_header_is_the_batch_header_plus_the_checkm2_columns(self):
        """A batch file has no CheckM2 columns: the runs are made once for the
        release, over the genomes every batch together found."""
        self.assertEqual(len(G.CONFLICT_HEADER_CHECKM2),
                         len(G.CONFLICT_HEADER) + len(G.CHECKM2_COLUMNS))
        for column in G.CONFLICT_HEADER:
            self.assertIn(column, G.CONFLICT_HEADER_CHECKM2)

    def test_the_columns_the_batch_writes_keep_their_positions(self):
        """concatenate() takes the header from the first batch file, so the
        release file starts as a batch file and is annotated afterwards."""
        self.assertEqual(G.CONFLICT_HEADER_CHECKM2[:len(G.CONFLICT_HEADER) - 1],
                         G.CONFLICT_HEADER[:-1])


class ReadConflictsTests(TempDirCase):
    """A batch that conflicted about nothing still writes its header."""

    def test_a_file_of_nothing_but_a_header_holds_no_conflicts(self):
        path = os.path.join(self.dir, G.CONFLICT_NAME)
        G.write_table([], path, G.CONFLICT_HEADER)
        self.assertEqual(G.read_conflicts(path), [])

    def test_a_row_is_read_back_as_it_was_written(self):
        path = os.path.join(self.dir, G.CONFLICT_NAME)
        row = conflict_row()
        G.write_table([row], path, G.CONFLICT_HEADER)
        self.assertEqual(G.read_conflicts(path), [row])

    def test_a_file_already_annotated_is_read_as_the_row_it_was_made_from(self):
        """The command is run again over a finished output directory, so the
        release file it reads is one it has already annotated; the columns are
        taken by name so the CheckM2 ones are simply not among them."""
        path = os.path.join(self.dir, G.CONFLICT_NAME)
        row = conflict_row()
        # built by the annotator rather than by hand, so that a column added to
        # CONFLICT_HEADER_CHECKM2 cannot quietly turn this into a row too short
        # to be read back at all
        annotated = G.annotate_conflicts([row],
                                         {25: {'GCA_1.1': ('94.3', '0.17')},
                                          11: {'GCA_1.1': ('51.0', '16.4')}})
        G.write_table(annotated, path, header=G.CONFLICT_HEADER_CHECKM2)
        self.assertEqual(len(annotated[0]), len(G.CONFLICT_HEADER_CHECKM2))
        self.assertEqual(G.read_conflicts(path), [row])

    def test_a_file_that_is_not_a_conflict_file_is_refused(self):
        """It is written by this command and read by it, so one that does not
        look like one is a bug, and annotating it would put one genome's quality
        against another genome's conflict."""
        path = os.path.join(self.dir, G.CONFLICT_NAME)
        with open(path, 'w') as handle:
            handle.write('genome_id\nGCA_1.1\n')
        self.assertRaises(G.BadConflictFile, G.read_conflicts, path)


class ReleaseFastasTests(TempDirCase):
    """Where a genome is, is taken from the batchfiles for the reason the
    comparison takes it from them: they say what the release was RUN on."""

    def batch(self, name, *rows):
        path = os.path.join(self.dir, name)
        os.makedirs(path)
        G.write_batchfile(rows, os.path.join(path, G.BATCHFILE_NAME))
        return path

    def test_a_genome_is_found_in_whichever_batch_holds_it(self):
        batches = [self.batch('batch_000001', ('/m/a.fna.gz', 'GCA_1.1')),
                   self.batch('batch_000002', ('/m/b.fna.gz', 'GCA_2.1'))]
        self.assertEqual(G.release_fastas(batches, ['GCA_2.1']),
                         {'GCA_2.1': '/m/b.fna.gz'})

    def test_only_the_genomes_asked_about_are_kept(self):
        """The batchfiles are the whole release and the conflicts are a few
        hundred of it."""
        batches = [self.batch('batch_000001', ('/m/a.fna.gz', 'GCA_1.1'),
                              ('/m/b.fna.gz', 'GCA_2.1'))]
        self.assertEqual(list(G.release_fastas(batches, ['GCA_1.1'])), ['GCA_1.1'])

    def test_a_batch_with_no_batchfile_is_skipped_rather_than_raising(self):
        batches = [self.batch('batch_000001', ('/m/a.fna.gz', 'GCA_1.1')),
                   os.path.join(self.dir, 'batch_000002')]
        os.makedirs(batches[1])
        self.assertEqual(G.release_fastas(batches, ['GCA_1.1']),
                         {'GCA_1.1': '/m/a.fna.gz'})


class EstimateConflictQualityTests(TempDirCase):
    """CheckM2 runs once for the release, and what it says lands in the release
    conflict file."""

    def setUp(self):
        super().setUp()
        self._check, G.check_dependencies = G.check_dependencies, lambda *a, **k: True
        self.out_dir = os.path.join(self.dir, 'out')
        os.makedirs(self.out_dir)

    def tearDown(self):
        G.check_dependencies = self._check
        super().tearDown()

    def conflicts(self, *rows):
        G.write_table(rows, os.path.join(self.out_dir, G.CONFLICT_NAME), G.CONFLICT_HEADER)

    def release(self):
        with open(os.path.join(self.out_dir, G.CONFLICT_NAME)) as handle:
            return [line.rstrip('\n').split('\t') for line in handle]

    def row(self, accession='GCA_1.1', gtranslate_tt='25', ncbi_tt='11'):
        return conflict_row(genome_id=accession, gtranslate_tt=gtranslate_tt,
                            ncbi_tt=ncbi_tt)

    def test_the_release_file_gains_the_checkm2_columns(self):
        self.conflicts(self.row())
        with mock.patch.object(G.GTranslate, 'checkm2_table', autospec=True,
                               return_value={'GCA_1.1': ('94.3', '0.17')}):
            G.GTranslate().estimate_conflict_quality([], self.out_dir)

        lines = self.release()
        self.assertEqual(tuple(lines[0]), G.CONFLICT_HEADER_CHECKM2)
        self.assertEqual(lines[1][G.CONFLICT_HEADER_CHECKM2.index(
            'cm2_completeness_gtranslate_tt')], '94.3')

    def test_one_run_is_made_for_each_table_in_dispute_and_no_more(self):
        """Two genomes disputing 25 against 11 are three runs between them at
        most, and are two: one per table."""
        self.conflicts(self.row('GCA_1.1'), self.row('GCA_2.1'))
        with mock.patch.object(G.GTranslate, 'checkm2_table', autospec=True,
                               return_value={}) as run:
            G.GTranslate().estimate_conflict_quality([], self.out_dir)

        self.assertEqual(sorted(call.args[4] for call in run.call_args_list),
                         [11, 25])

    def test_a_file_that_cannot_be_read_costs_the_columns_and_not_the_run(self):
        """A run of five machines over days should not end in a traceback for
        want of four columns."""
        with open(os.path.join(self.out_dir, G.CONFLICT_NAME), 'w') as handle:
            handle.write('genome_id\nGCA_1.1\n')
        with mock.patch.object(G.GTranslate, 'checkm2_table', autospec=True) as run:
            G.GTranslate().estimate_conflict_quality([], self.out_dir)

        self.assertFalse(run.called)

    def test_a_release_that_conflicted_about_nothing_runs_checkm2_not_at_all(self):
        self.conflicts()
        with mock.patch.object(G.GTranslate, 'checkm2_table', autospec=True) as run:
            G.GTranslate().estimate_conflict_quality([], self.out_dir)

        self.assertFalse(run.called)

    def test_annotating_twice_does_not_double_the_columns(self):
        """The command is run again over a finished output directory to pick this
        up at all, so the release file it reads is one it has already written."""
        self.conflicts(self.row())
        with mock.patch.object(G.GTranslate, 'checkm2_table', autospec=True,
                               return_value={'GCA_1.1': ('94.3', '0.17')}):
            G.GTranslate().estimate_conflict_quality([], self.out_dir)
            G.GTranslate().estimate_conflict_quality([], self.out_dir)

        self.assertEqual(tuple(self.release()[0]), G.CONFLICT_HEADER_CHECKM2)
        self.assertEqual(len(self.release()[1]), len(G.CONFLICT_HEADER_CHECKM2))


class CheckM2RunTests(TempDirCase):
    """A run already made is not made again, and a run that fails costs the
    release four columns and nothing else."""

    def setUp(self):
        super().setUp()
        self._check, G.check_dependencies = G.check_dependencies, lambda *a, **k: True
        self.checkm2_dir = os.path.join(self.dir, 'checkm2')

    def tearDown(self):
        G.check_dependencies = self._check
        super().tearDown()

    def report(self, table, *rows):
        path = os.path.join(self.checkm2_dir, G.CHECKM2_TABLE_DIR.format(table))
        os.makedirs(path)
        with open(os.path.join(path, G.CHECKM2_REPORT), 'w') as handle:
            handle.write('Name\tCompleteness\tContamination\n')
            for row in rows:
                handle.write('\t'.join(row) + '\n')

    def test_a_table_already_run_is_read_rather_than_run_again(self):
        """The release is finished by whichever machine happens to be last, and
        every run over a finished output directory reaches this."""
        self.report(25, ('GCA_1.1', '94.3', '0.17'))
        with mock.patch.object(G.subprocess, 'run') as run:
            quality = G.GTranslate().checkm2_table(
                ['GCA_1.1'], {}, self.checkm2_dir, 25)

        self.assertFalse(run.called)
        self.assertEqual(quality, {'GCA_1.1': ('94.3', '0.17')})

    def test_a_table_no_genome_could_be_staged_for_is_not_run(self):
        """A run with no genomes in it is a run that fails."""
        with mock.patch.object(G.subprocess, 'run') as run:
            quality = G.GTranslate().checkm2_table(
                ['GCA_1.1'], {'GCA_1.1': '/gone.fna.gz'}, self.checkm2_dir, 25)

        self.assertFalse(run.called)
        self.assertEqual(quality, {})

    def test_a_failed_run_gives_back_nothing_rather_than_raising(self):
        """The conflicts are found, written and counted before this runs, and
        they are what the command is for."""
        path = self.genome_dir('GCA_1.1_ASM1')
        fasta = os.path.join(path, 'GCA_1.1_ASM1_genomic.fna.gz')
        with mock.patch.object(G.subprocess, 'run',
                               return_value=mock.Mock(returncode=1)):
            quality = G.GTranslate().checkm2_table(
                ['GCA_1.1'], {'GCA_1.1': fasta}, self.checkm2_dir, 25)

        self.assertEqual(quality, {})

    def test_the_staged_genomes_are_not_written_under_the_run_directory(self):
        """--force empties the output directory before CheckM2 starts, so
        staging under it would delete the genomes about to be read."""
        path = self.genome_dir('GCA_1.1_ASM1')
        fasta = os.path.join(path, 'GCA_1.1_ASM1_genomic.fna.gz')
        seen = {}

        def record(cmd, **kwargs):
            seen['input'] = cmd[cmd.index('--input') + 1:]
            seen['out'] = cmd[cmd.index('--output-directory') + 1]
            return mock.Mock(returncode=1)

        with mock.patch.object(G.subprocess, 'run', side_effect=record):
            G.GTranslate().checkm2_table(['GCA_1.1'], {'GCA_1.1': fasta},
                                         self.checkm2_dir, 25)

        for staged in seen['input']:
            self.assertFalse(staged.startswith(seen['out'] + os.sep))


class ProgramVersionTests(TempDirCase):
    """What predicted a batch, and what estimated the quality of the conflicts,
    are recorded beside what they made."""

    def setUp(self):
        super().setUp()
        self._check, G.check_dependencies = G.check_dependencies, lambda *a, **k: True
        self.out_dir = os.path.join(self.dir, 'out')
        genome_dirs = os.path.join(self.dir, 'genome_dirs.tsv')
        with open(genome_dirs, 'w') as handle:
            handle.write('GCA_000001.1\t{}\tG000001\n'.format(
                self.genome_dir('GCA_000001.1_ASM1')))
        self.batch = G.GTranslate().plan_batches(genome_dirs, self.out_dir)[0]

    def tearDown(self):
        G.check_dependencies = self._check
        super().tearDown()

    def predict(self, returncode):
        with mock.patch.object(G.subprocess, 'run',
                               return_value=mock.Mock(returncode=returncode)):
            try:
                G.GTranslate().run_gtranslate(self.batch)
            except (RuntimeError, OSError):
                # a failed run raises, and a run that succeeded goes on to read
                # the summary the stub never wrote -- after the version is down
                pass

    def test_a_predicted_batch_records_the_gtranslate_that_predicted_it(self):
        self.predict(0)

        with open(os.path.join(self.batch, 'gtranslate.version')) as handle:
            self.assertEqual(handle.read(), VERSIONS[G.GTRANSLATE_BIN] + '\n')

    def test_a_batch_gtranslate_failed_on_records_no_version(self):
        self.predict(1)

        self.assertFalse(os.path.exists(os.path.join(self.batch, 'gtranslate.version')))


# ------------------------------------------------------------- standard GTDB QC

class PassesQCTests(unittest.TestCase):
    """All three conditions have to hold, and a genome with no estimate is not a
    genome that failed."""

    def test_a_good_genome_passes(self):
        self.assertIs(G.passes_qc('94.3', '0.17'), True)

    def test_a_genome_half_there_or_less_fails_however_clean(self):
        """The threshold is exclusive: exactly 50% complete does not pass."""
        self.assertIs(G.passes_qc('50.0', '0.0'), False)
        self.assertIs(G.passes_qc('50.01', '0.0'), True)

    def test_a_genome_at_ten_percent_contamination_fails_however_complete(self):
        """Exclusive too, and it is the condition the score alone would miss: 100
        complete against 10 contaminated scores 50 and is refused twice over."""
        self.assertIs(G.passes_qc('100.0', '10.0'), False)
        self.assertIs(G.passes_qc('100.0', '9.99'), True)

    def test_the_score_refuses_a_genome_the_other_two_conditions_would_keep(self):
        """60% complete at 2.5% contamination is more than half there and barely
        contaminated, and scores 47.5. This is why the score is asked for."""
        self.assertIs(G.passes_qc('60.0', '2.5'), False)
        self.assertIs(G.passes_qc('60.0', '1.9'), True)

    def test_contamination_is_charged_at_five_times_the_weight(self):
        self.assertEqual(G.quality_score(94.3, 0.17), 94.3 - 5 * 0.17)

    def test_a_genome_with_no_estimate_is_unknown_rather_than_failed(self):
        """na is not False: a genome whose FASTA is missing was not looked at and
        found wanting."""
        self.assertIsNone(G.passes_qc(G.NCBI_NA, G.NCBI_NA))
        self.assertIsNone(G.passes_qc('94.3', G.NCBI_NA))
        self.assertIsNone(G.passes_qc('', ''))


class QCColumnTests(unittest.TestCase):
    """The verdict sits beside the numbers it was reached from."""

    def row(self):
        return conflict_row()

    def field(self, row, column):
        return row[G.CONFLICT_HEADER_CHECKM2.index(column)]

    def test_each_table_is_judged_on_its_own_estimate(self):
        """A genome can pass under one table and fail under the other; that it
        does is the whole reason for running both."""
        row = G.annotate_conflicts([self.row()],
                                   {25: {'GCA_1.1': ('94.3', '0.17')},
                                    11: {'GCA_1.1': ('51.0', '16.4')}})[0]
        self.assertEqual(self.field(row, 'pass_qc_gtranslate_tt'), 'True')
        self.assertEqual(self.field(row, 'pass_qc_ncbi_tt'), 'False')

    def test_a_verdict_follows_the_numbers_it_was_reached_from(self):
        """Read by eye the row is two answers to one question, not four numbers
        and two verdicts at the end."""
        header = G.CONFLICT_HEADER_CHECKM2
        self.assertEqual(header.index('pass_qc_gtranslate_tt'),
                         header.index('cm2_contamination_gtranslate_tt') + 1)
        self.assertEqual(header.index('pass_qc_ncbi_tt'),
                         header.index('cm2_contamination_ncbi_tt') + 1)

    def test_a_genome_checkm2_returned_nothing_for_is_not_reported_as_failing(self):
        row = G.annotate_conflicts([self.row()], {})[0]
        self.assertEqual(self.field(row, 'pass_qc_gtranslate_tt'), G.NCBI_NA)
        self.assertEqual(self.field(row, 'pass_qc_ncbi_tt'), G.NCBI_NA)


class RemoveCheckM2DirTests(TempDirCase):
    """The working directory goes once what it was for is in the conflict file."""

    def setUp(self):
        super().setUp()
        self._check, G.check_dependencies = G.check_dependencies, lambda *a, **k: True
        self.out_dir = os.path.join(self.dir, 'out')
        os.makedirs(self.out_dir)
        self.checkm2_dir = os.path.join(self.out_dir, G.CHECKM2_DIR)

    def tearDown(self):
        G.check_dependencies = self._check
        super().tearDown()

    def conflicts(self, *rows):
        G.write_table(rows, os.path.join(self.out_dir, G.CONFLICT_NAME), G.CONFLICT_HEADER)

    def row(self, accession='GCA_1.1'):
        return conflict_row(genome_id=accession)

    def estimate(self, quality):
        """Run the stage with CheckM2 standing in, leaving a directory behind."""
        def stand_in(manager, accessions, fastas, checkm2_dir, table):
            os.makedirs(os.path.join(checkm2_dir, G.CHECKM2_TABLE_DIR.format(table)),
                        exist_ok=True)
            return quality.get(table, {})

        with mock.patch.object(G.GTranslate, 'checkm2_table', autospec=True,
                               side_effect=stand_in):
            G.GTranslate().estimate_conflict_quality([], self.out_dir)

    def test_the_directory_is_gone_once_every_table_has_been_run(self):
        self.conflicts(self.row())
        self.estimate({25: {'GCA_1.1': ('94.3', '0.17')},
                       11: {'GCA_1.1': ('51.0', '16.4')}})
        self.assertFalse(os.path.exists(self.checkm2_dir))

    def test_the_estimates_are_in_the_table_before_the_directory_goes(self):
        """Removing it is the last thing done, so nothing is swept up that has
        not landed."""
        self.conflicts(self.row())
        self.estimate({25: {'GCA_1.1': ('94.3', '0.17')},
                       11: {'GCA_1.1': ('51.0', '16.4')}})
        with open(os.path.join(self.out_dir, G.CONFLICT_NAME)) as handle:
            rows = handle.read().splitlines()
        self.assertEqual(rows[1].split('\t')[
            G.CONFLICT_HEADER_CHECKM2.index('cm2_completeness_gtranslate_tt')], '94.3')

    def test_a_table_that_produced_nothing_keeps_the_directory(self):
        """The tables that succeeded are then read rather than run again, which
        is what makes retrying a failed table cheap."""
        self.conflicts(self.row())
        self.estimate({25: {'GCA_1.1': ('94.3', '0.17')}, 11: {}})
        self.assertTrue(os.path.isdir(self.checkm2_dir))

    def test_a_release_that_conflicted_about_nothing_removes_nothing(self):
        """There is no directory to remove, and rmtree on a path that was never
        made should not be how that is found out."""
        self.conflicts()
        self.estimate({})
        self.assertFalse(os.path.exists(self.checkm2_dir))

    def test_removing_a_directory_that_is_not_there_is_not_an_error(self):
        G.remove_checkm2_dir(os.path.join(self.dir, 'never_made'))


# ------------------------------------------------- the release summary is gzipped

class ReleaseSummaryNameTests(unittest.TestCase):
    """A batch's summary is gTranslate's output; the release's is ours."""

    def test_the_release_summary_is_the_batch_name_gzipped(self):
        self.assertEqual(G.release_summary_name(),
                         G.summary_name() + G.GZIP_EXT)

    def test_the_prefix_is_carried_through(self):
        """--prefix names gTranslate's output, and the release file is named for
        the same run."""
        self.assertEqual(G.release_summary_name('run7'),
                         'run7.translation_table_summary.tsv.gz')


class CompressedConcatenateTests(TempDirCase):
    """The summary is a row per genome of the release; the conflicts are a few
    hundred rows meant to be looked at."""

    def batch_summary(self, name, *rows):
        path = os.path.join(self.dir, name)
        with open(path, 'w') as handle:
            handle.write('user_genome\tbest_tln_table\n')
            for row in rows:
                handle.write('\t'.join(row) + '\n')
        return path

    def test_a_compressed_join_is_gzip_and_reads_back_as_the_rows(self):
        first = self.batch_summary('a.tsv', ('GCA_1.1', '25'))
        second = self.batch_summary('b.tsv', ('GCA_2.1', '11'))
        out = os.path.join(self.dir, 'all.tsv.gz')
        self.assertEqual(G.concatenate([first, second], out, compress=True), 2)

        with open(out, 'rb') as handle:
            self.assertEqual(handle.read(2), b'\x1f\x8b')
        self.assertEqual(sorted(G.read_translation_table_summary(out)),
                         ['GCA_1.1', 'GCA_2.1'])

    def test_the_header_is_still_written_once(self):
        first = self.batch_summary('a.tsv', ('GCA_1.1', '25'))
        second = self.batch_summary('b.tsv', ('GCA_2.1', '11'))
        out = os.path.join(self.dir, 'all.tsv.gz')
        G.concatenate([first, second], out, compress=True)
        with gzip.open(out, 'rt') as handle:
            lines = handle.read().splitlines()
        self.assertEqual(len(lines), 3)
        self.assertEqual(lines[0].split('\t')[0], 'user_genome')

    def test_an_uncompressed_join_is_still_plain_text(self):
        """ncbi_tt_conflict.tsv is read by eye and is not compressed."""
        first = self.batch_summary('a.tsv', ('GCA_1.1', '25'))
        out = os.path.join(self.dir, 'all.tsv')
        G.concatenate([first], out)
        with open(out, 'rb') as handle:
            self.assertNotEqual(handle.read(2), b'\x1f\x8b')


class AggregateSummaryTests(TempDirCase):
    """What the release directory holds afterwards."""

    def setUp(self):
        super().setUp()
        self._check, G.check_dependencies = G.check_dependencies, lambda *a, **k: True
        self.out_dir = os.path.join(self.dir, 'out')
        self.batch = os.path.join(self.out_dir, 'batch_000001')
        os.makedirs(self.batch)
        G.write_table([], os.path.join(self.batch, G.CONFLICT_NAME), G.CONFLICT_HEADER)
        G.write_table([], os.path.join(self.batch, G.COMPARISON_NAME),
                      header=G.COMPARISON_HEADER, compress=True)
        G.write_batchfile([('/m/a.fna.gz', 'GCA_1.1')],
                          os.path.join(self.batch, G.BATCHFILE_NAME), compress=True)
        with open(os.path.join(self.batch, G.summary_name()), 'w') as handle:
            handle.write('user_genome\tbest_tln_table\nGCA_1.1\t25\n')
        G.finish_batch(self.batch, compared=1, conflicts=0, no_ncbi_table=0)

    def tearDown(self):
        G.check_dependencies = self._check
        super().tearDown()

    def aggregate(self):
        with mock.patch.object(G.GTranslate, 'checkm2_table', autospec=True,
                               return_value={}):
            G.GTranslate().aggregate([self.batch], self.out_dir)

    def test_the_release_summary_is_written_gzipped(self):
        self.aggregate()
        path = os.path.join(self.out_dir, G.release_summary_name())
        self.assertTrue(os.path.exists(path))
        self.assertEqual(G.read_translation_table_summary(path)['GCA_1.1'
                                                                ]['best_tln_table'], '25')

    def test_the_batch_summary_is_left_as_gtranslate_wrote_it(self):
        """It is gTranslate's output and not this command's to compress."""
        self.aggregate()
        with open(os.path.join(self.batch, G.summary_name()), 'rb') as handle:
            self.assertNotEqual(handle.read(2), b'\x1f\x8b')

    def test_an_uncompressed_summary_left_by_an_older_run_is_removed(self):
        """Two files a genome apart of which one is stale is how a release gets
        called under the wrong tables."""
        stale = os.path.join(self.out_dir, G.summary_name())
        with open(stale, 'w') as handle:
            handle.write('user_genome\tbest_tln_table\nGCA_1.1\t11\n')
        self.aggregate()
        self.assertFalse(os.path.exists(stale))
        self.assertTrue(os.path.exists(
            os.path.join(self.out_dir, G.release_summary_name())))


# ------------------------------------------------ the comparison of every genome

class ComparisonFileTests(TempDirCase):
    """Every genome the two both called, with the conflicts derived from it."""

    def genome_dir(self, accession, ncbi_table):
        assembly = '{}_ASM1'.format(accession)
        path = os.path.join(self.dir, assembly)
        os.makedirs(path)
        with gzip.open(os.path.join(path, assembly + '_genomic.gff.gz'), 'wt') as handle:
            handle.write('##gff-version 3\n')
            handle.write('c\tRefSeq\tCDS\t1\t9\t.\t+\t0\t'
                         'ID=cds1;product=x;transl_table={}\n'.format(ncbi_table))
        return path

    def summary(self, *rows):
        """A summary with the statistics gTranslate actually reports."""
        path = os.path.join(self.dir, 'gtranslate.translation_table_summary.tsv')
        with open(path, 'w') as handle:
            handle.write('user_genome\tbest_tln_table\tcoding_density_4\t'
                         'coding_density_11\tgc_percent\tn50\tgenome_size\n')
            for accession, table in rows:
                handle.write('{}\t{}\t90.1\t64.2\t35.33882\t5269725\t5269725\n'.format(
                    accession, table))
        return path

    def compare(self, *genomes):
        summary = self.summary(*[(accession, table) for accession, table, _ in genomes])
        dirs = {accession: self.genome_dir(accession, ncbi)
                for accession, _, ncbi in genomes}
        return G.comparison_rows(
            G.read_translation_table_summary(summary), dirs, {})

    def field(self, row, column):
        return row[G.COMPARISON_HEADER.index(column)]

    def test_the_statistics_are_carried_from_the_summary(self):
        """gTranslate measured them of the FASTA it read; measuring them here
        would be reading 1.3M genomes again to learn what is already known."""
        rows, _, _ = self.compare(('GCF_1.1', '4', 11))
        self.assertEqual(self.field(rows[0], 'gc_percent'), '35.33882')
        self.assertEqual(self.field(rows[0], 'n50'), '5269725')
        self.assertEqual(self.field(rows[0], 'genome_size'), '5269725')

    def test_a_summary_without_the_statistics_reports_them_as_unknown(self):
        """An older gTranslate reported fewer columns, and an empty field in a
        TSV does not say which of absent and zero it means."""
        summary = os.path.join(self.dir, 'old.tsv')
        with open(summary, 'w') as handle:
            handle.write('user_genome\tbest_tln_table\tcoding_density_4\t'
                         'coding_density_11\n')
            handle.write('GCF_1.1\t4\t90.1\t64.2\n')
        rows, _, _ = G.comparison_rows(
            G.read_translation_table_summary(summary),
            {'GCF_1.1': self.genome_dir('GCF_1.1', 11)}, {})
        self.assertEqual(self.field(rows[0], 'gc_percent'), G.NCBI_NA)
        self.assertEqual(self.field(rows[0], 'n50'), G.NCBI_NA)
        self.assertEqual(self.field(rows[0], 'genome_size'), G.NCBI_NA)

    def test_agreements_and_disagreements_are_both_kept(self):
        rows, compared, _ = self.compare(('GCF_1.1', '11', 11), ('GCF_2.1', '4', 11))
        self.assertEqual(compared, 2)
        self.assertEqual([self.field(row, 'ncbi_conflict') for row in rows],
                         ['False', 'True'])

    def test_the_conflicts_are_the_comparison_filtered(self):
        """One walk over the genomes, not two: the GFF of every genome has
        already been read once to make the comparison."""
        rows, _, _ = self.compare(('GCF_1.1', '11', 11), ('GCF_2.1', '4', 11))
        conflicts = G.conflicts_from_comparison(rows)
        self.assertEqual([row[0] for row in conflicts], ['GCF_2.1'])

    def test_a_conflict_keeps_every_column_the_comparison_gave_it(self):
        rows, _, _ = self.compare(('GCF_2.1', '4', 11))
        conflict = G.conflicts_from_comparison(rows)[0]
        self.assertEqual(len(conflict), len(G.CONFLICT_HEADER))
        for column in ('gc_percent', 'n50', 'genome_size', 'checkm_conflict',
                       'coding_density_4', 'ncbi_taxonomy'):
            self.assertEqual(conflict[G.CONFLICT_HEADER.index(column)],
                             self.field(rows[0], column))

    def test_every_conflict_column_comes_from_the_comparison(self):
        """conflicts_from_comparison() takes them by name, so a column added to
        one file and not the other is a failure here and not in a release."""
        for column in G.CONFLICT_HEADER:
            self.assertIn(column, G.COMPARISON_HEADER)

    def test_a_genome_ncbi_declares_no_table_for_is_in_neither_file(self):
        """It has nothing to be compared against and no ncbi_conflict to report."""
        assembly = 'GCF_9.1_ASM1'
        path = os.path.join(self.dir, assembly)
        os.makedirs(path)
        rows, compared, no_table = G.comparison_rows(
            G.read_translation_table_summary(self.summary(('GCF_9.1', '11'))),
            {'GCF_9.1': path}, {})
        self.assertEqual(rows, [])
        self.assertEqual(compared, 0)
        self.assertEqual(no_table, 1)

    def test_the_batch_comparison_is_written_gzipped(self):
        """A row per genome compared, where the conflicts are a few hundred."""
        path = os.path.join(self.dir, G.COMPARISON_NAME)
        rows, _, _ = self.compare(('GCF_1.1', '4', 11))
        G.write_table(rows, path, header=G.COMPARISON_HEADER, compress=True)
        with open(path, 'rb') as handle:
            self.assertEqual(handle.read(2), b'\x1f\x8b')
        with gzip.open(path, 'rt') as handle:
            lines = handle.read().splitlines()
        self.assertEqual(lines[0].split('\t'), list(G.COMPARISON_HEADER))
        self.assertEqual(len(lines[1].split('\t')), len(G.COMPARISON_HEADER))

    def test_the_name_says_it_is_gzipped(self):
        self.assertTrue(G.COMPARISON_NAME.endswith('.tsv' + G.GZIP_EXT))


class ComparisonAggregationTests(TempDirCase):
    """The release comparison is every batch's, joined."""

    def setUp(self):
        super().setUp()
        self._check, G.check_dependencies = G.check_dependencies, lambda *a, **k: True
        self.out_dir = os.path.join(self.dir, 'out')
        self.batches = []
        for index, accession in enumerate(('GCA_1.1', 'GCA_2.1'), start=1):
            batch = os.path.join(self.out_dir, 'batch_{:06d}'.format(index))
            os.makedirs(batch)
            row = [accession, '25', '11', '4', 'True', 'False', '90.1', '64.2',
                   '35.3', '5269725', '5269725', 'd__Bacteria']
            G.write_table([row], os.path.join(batch, G.COMPARISON_NAME),
                          header=G.COMPARISON_HEADER, compress=True)
            G.write_table([], os.path.join(batch, G.CONFLICT_NAME), G.CONFLICT_HEADER)
            G.write_batchfile([('/m/{}.fna.gz'.format(accession), accession)],
                              os.path.join(batch, G.BATCHFILE_NAME), compress=True)
            with open(os.path.join(batch, G.summary_name()), 'w') as handle:
                handle.write('user_genome\tbest_tln_table\n{}\t25\n'.format(accession))
            G.finish_batch(batch, compared=1, conflicts=1, no_ncbi_table=0)
            self.batches.append(batch)

    def tearDown(self):
        G.check_dependencies = self._check
        super().tearDown()

    def test_the_release_comparison_joins_every_batch_under_one_header(self):
        with mock.patch.object(G.GTranslate, 'checkm2_table', autospec=True,
                               return_value={}):
            G.GTranslate().aggregate(self.batches, self.out_dir)

        path = os.path.join(self.out_dir, G.COMPARISON_NAME)
        with gzip.open(path, 'rt') as handle:
            lines = handle.read().splitlines()
        self.assertEqual(lines[0].split('\t'), list(G.COMPARISON_HEADER))
        self.assertEqual([line.split('\t')[0] for line in lines[1:]],
                         ['GCA_1.1', 'GCA_2.1'])

    def test_a_gzipped_batch_file_is_joined_without_being_decompressed_first(self):
        """concatenate() reads what it is given, which for the comparison is
        gzip and for the conflicts is text."""
        out = os.path.join(self.dir, 'joined.tsv.gz')
        self.assertEqual(G.concatenate(
            [os.path.join(batch, G.COMPARISON_NAME) for batch in self.batches],
            out, compress=True), 2)


# ------------------------------------------------------- the batch plan is gzipped

class BatchfileCleanupTests(TempDirCase):
    """The uncompressed copy exists only while the batch is being worked on."""

    def setUp(self):
        super().setUp()
        self._check, G.check_dependencies = G.check_dependencies, lambda *a, **k: True
        self.batch_dir = os.path.join(self.dir, 'batch_000001')
        os.makedirs(self.batch_dir)
        path = self.genome_dir('GCF_1.1_ASM1')
        G.write_batchfile([(G.genomic_fasta(path), 'GCF_1.1')],
                          os.path.join(self.batch_dir, G.BATCHFILE_NAME),
                          compress=True)

    def tearDown(self):
        G.check_dependencies = self._check
        super().tearDown()

    def predict(self, returncode=0):
        manager = G.GTranslate(prefix=None)
        with open(os.path.join(self.batch_dir, G.summary_name()), 'w') as handle:
            handle.write('user_genome\tbest_tln_table\nGCF_1.1\t11\n')
        with mock.patch.object(G.subprocess, 'run',
                               return_value=mock.Mock(returncode=returncode)):
            manager.run_gtranslate(self.batch_dir)

    def test_the_plain_copy_is_gone_once_gtranslate_has_read_it(self):
        self.predict()
        self.assertFalse(os.path.exists(
            os.path.join(self.batch_dir, G.PRESENT_BATCHFILE_NAME)))

    def test_the_compressed_plan_is_what_the_batch_keeps(self):
        self.predict()
        self.assertTrue(os.path.exists(
            os.path.join(self.batch_dir, G.BATCHFILE_NAME)))
        self.assertEqual([accession for _, accession in G.read_batchfile(
            G.batchfile_path(self.batch_dir, G.LAYOUT))], ['GCF_1.1'])

    def test_a_batch_gtranslate_failed_on_keeps_the_copy_it_was_given(self):
        """The batch is retried, and what it was handed is what the retry looks
        at to see what went in."""
        self.assertRaises(RuntimeError, self.predict, returncode=1)
        self.assertTrue(os.path.exists(
            os.path.join(self.batch_dir, G.PRESENT_BATCHFILE_NAME)))


# --------------------------------------- the genomes the release has no table for

class NoPredictionReleaseTests(TempDirCase):
    """The one file that says which genomes the release has no answer for."""

    def setUp(self):
        super().setUp()
        self._check, G.check_dependencies = G.check_dependencies, lambda *a, **k: True
        self.out_dir = os.path.join(self.dir, 'out')

    def tearDown(self):
        G.check_dependencies = self._check
        super().tearDown()

    def batch(self, index, planned, predicted, no_fasta=()):
        batch_dir = os.path.join(self.out_dir, 'batch_{:06d}'.format(index))
        os.makedirs(batch_dir)
        G.write_batchfile([('/m/{}.fna.gz'.format(a), a) for a in planned],
                          os.path.join(batch_dir, G.BATCHFILE_NAME), compress=True)
        with open(os.path.join(batch_dir, G.summary_name()), 'w') as handle:
            handle.write('user_genome\tbest_tln_table\n')
            for accession in predicted:
                handle.write('{}\t11\n'.format(accession))
        if no_fasta:
            with open(os.path.join(batch_dir, G.MISSING_NAME), 'w') as handle:
                for accession in no_fasta:
                    handle.write('{}\n'.format(accession))
        return batch_dir

    def test_a_genome_gtranslate_returned_nothing_for_is_named(self):
        batch = self.batch(1, ['GCA_1.1', 'GCA_2.1'], ['GCA_1.1'])
        self.assertEqual(G.no_prediction_rows([batch]),
                         [('GCA_2.1', G.REASON_NO_PREDICTION)])

    def test_a_genome_with_no_fasta_is_reported_as_that_and_not_as_a_failure(self):
        """gTranslate was never given it, so calling it a prediction failure
        would send whoever reads this looking at the wrong thing."""
        batch = self.batch(1, ['GCA_1.1', 'GCA_2.1'], ['GCA_1.1'],
                           no_fasta=['GCA_2.1'])
        self.assertEqual(G.no_prediction_rows([batch]),
                         [('GCA_2.1', G.REASON_NO_FASTA)])

    def test_the_two_kinds_are_both_reported_and_told_apart(self):
        batch = self.batch(1, ['GCA_1.1', 'GCA_2.1', 'GCA_3.1'], ['GCA_1.1'],
                           no_fasta=['GCA_3.1'])
        self.assertEqual(G.no_prediction_rows([batch]),
                         [('GCA_2.1', G.REASON_NO_PREDICTION),
                          ('GCA_3.1', G.REASON_NO_FASTA)])

    def test_every_batch_of_the_release_is_gathered_in_accession_order(self):
        """They were in 135 directories, to be found by whoever thought to look."""
        batches = [self.batch(1, ['GCA_9.1', 'GCA_1.1'], ['GCA_1.1']),
                   self.batch(2, ['GCA_5.1'], [])]
        self.assertEqual([row[0] for row in G.no_prediction_rows(batches)],
                         ['GCA_5.1', 'GCA_9.1'])

    def test_a_batch_that_never_wrote_no_prediction_tsv_is_still_accounted_for(self):
        """It is worked out from what the batch was asked against what it
        answered, not read from the file the predicting run wrote -- a batch
        already predicted is never predicted again, so that file may not be
        there."""
        batch = self.batch(1, ['GCA_1.1', 'GCA_2.1'], ['GCA_1.1'])
        self.assertFalse(os.path.exists(os.path.join(batch, G.NO_PREDICTION_NAME)))
        self.assertEqual(len(G.no_prediction_rows([batch])), 1)

    def test_a_release_with_nothing_missing_yields_no_rows(self):
        batch = self.batch(1, ['GCA_1.1'], ['GCA_1.1'])
        self.assertEqual(G.no_prediction_rows([batch]), [])

    def test_the_file_is_written_even_when_there_is_nothing_to_report(self):
        """A release with nothing missing says so, rather than leaving the
        question open."""
        self.batch(1, ['GCA_1.1'], ['GCA_1.1'])
        G.GTranslate().report_no_prediction(
            [os.path.join(self.out_dir, 'batch_000001')], self.out_dir)
        path = os.path.join(self.out_dir, G.NO_PREDICTION_RELEASE_NAME)
        self.assertEqual(open(path).read().splitlines(),
                         ['\t'.join(G.NO_PREDICTION_HEADER)])

    def test_the_file_names_the_genomes_and_why(self):
        batch = self.batch(1, ['GCA_1.1', 'GCA_2.1', 'GCA_3.1'], ['GCA_1.1'],
                           no_fasta=['GCA_3.1'])
        G.GTranslate().report_no_prediction([batch], self.out_dir)
        with open(os.path.join(self.out_dir, G.NO_PREDICTION_RELEASE_NAME)) as handle:
            lines = handle.read().splitlines()
        self.assertEqual(lines[0].split('\t'), list(G.NO_PREDICTION_HEADER))
        self.assertEqual(lines[1:], ['GCA_2.1\t' + G.REASON_NO_PREDICTION,
                                     'GCA_3.1\t' + G.REASON_NO_FASTA])

    def test_the_reasons_are_tallied_for_the_log(self):
        rows = [('GCA_1.1', G.REASON_NO_PREDICTION),
                ('GCA_2.1', G.REASON_NO_PREDICTION),
                ('GCA_3.1', G.REASON_NO_FASTA)]
        self.assertEqual(G.tally_reasons(rows),
                         {G.REASON_NO_PREDICTION: 2, G.REASON_NO_FASTA: 1})

    def test_a_batch_with_no_plan_is_skipped_rather_than_raising(self):
        empty = os.path.join(self.out_dir, 'batch_000009')
        os.makedirs(empty)
        self.assertEqual(G.no_prediction_rows([empty]), [])
