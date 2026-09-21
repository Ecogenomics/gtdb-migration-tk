#!/usr/bin/env python3
"""Offline unit tests for prodigal_manager.py -- Prodigal itself is never run.

Run with the interpreter that has tqdm:

    /opt/centos7/sw/miniconda3/envs/python3/bin/python -m unittest -v tests.test_prodigal_manager

The vendored Prodigal wrapper is replaced by a stub that records the translation
tables it was handed and writes the files run_prodigal() then files away, so what
is tested here is the bookkeeping around the gene calling: which table each genome
is called under, where that table came from, and that a release the predictions do
not cover is refused before any genes are called rather than hours into a run.
"""

import collections
import gzip
import os
import shutil
import tempfile
import unittest
from unittest import mock

from gtdb_migration_tk import batching as B
from gtdb_migration_tk import prodigal_manager as P


Summary = collections.namedtuple(
    'Summary', 'best_translation_table coding_density_4 coding_density_11')


class StubProdigal:
    """Stands in for the vendored wrapper, which would call the real Prodigal."""

    last_tasks = None

    def __init__(self, cpus=1):
        self.cpus = cpus

    def run(self, tasks):
        """Write what the real wrapper writes: the results where the task says."""
        StubProdigal.last_tasks = list(tasks)

        summary_stats = {}
        for task in tasks:
            for path in (task.aa_gene_file, task.nt_gene_file, task.gff_file):
                os.makedirs(os.path.dirname(path), exist_ok=True)
                with gzip.open(path, 'wt') as handle:
                    handle.write('>gene\nMA\n')
            if task.checksum_file:
                # of the UNCOMPRESSED bytes, as the wrapper writes it and as
                # prodigal_parser() reads it back to decide a genome can be skipped
                with open(task.checksum_file, 'w') as handle:
                    handle.write(P.sha256_rb(
                        gzip.GzipFile(fileobj=open(task.aa_gene_file, 'rb'))) + '\n')
            summary_stats[task.genome_id] = Summary(task.translation_table, -1, -1)

        return summary_stats

    @classmethod
    def genomes_called(cls):
        """The genomes the wrapper was handed, in the order it was given them."""
        return [task.genome_id for task in cls.last_tasks or []]

    @classmethod
    def table_for(cls, accession):
        """The table the wrapper was handed for one genome."""
        for task in cls.last_tasks or []:
            if task.genome_id == accession:
                return task.translation_table
        raise KeyError(accession)


class TempDirCase(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='prodigal_manager_test.')
        self.tmp_dir = os.path.join(self.dir, 'tmp')
        os.makedirs(self.tmp_dir)
        # --out_dir holds the batches and the state of the run, never the genes
        self.out_dir = os.path.join(self.dir, 'out')

        # neither the real Prodigal nor the check for it is wanted offline
        self._prodigal, P.Prodigal = P.Prodigal, StubProdigal
        self._check, P.check_dependencies = P.check_dependencies, lambda *a, **k: True
        StubProdigal.last_tasks = None

    def tearDown(self):
        P.Prodigal = self._prodigal
        P.check_dependencies = self._check
        shutil.rmtree(self.dir, ignore_errors=True)

    def genome(self, accession):
        """A genome directory as a release holds it, named for its assembly."""
        assembly = '{}_ASM1'.format(accession)
        gpath = os.path.join(self.dir, assembly)
        os.makedirs(gpath)
        with gzip.open(os.path.join(gpath, assembly + '_genomic.fna.gz'), 'wt') as handle:
            handle.write('>contig\nACGT\n')
        return accession, gpath

    def genome_dirs(self, *genomes):
        path = os.path.join(self.dir, 'genome_dirs.tsv')
        with open(path, 'w') as handle:
            for accession, gpath in genomes:
                handle.write('{}\t{}\tG{}\n'.format(accession, gpath, accession[4:13]))
        return path

    def summary(self, *rows):
        """The translation table summary gTranslate writes."""
        path = os.path.join(self.dir, 'gtranslate.translation_table_summary.tsv')
        with open(path, 'w') as handle:
            handle.write('user_genome\tbest_tln_table\tcoding_density_4\t'
                         'coding_density_11\tconfidence\n')
            for accession, table in rows:
                handle.write('{}\t{}\t90.1\t64.2\t1.0\n'.format(accession, table))
        return path

    def override(self, *rows):
        path = os.path.join(self.dir, 'tt_override.tsv')
        with open(path, 'w') as handle:
            handle.write('genome_id\ttranslation_table\n')
            for accession, table in rows:
                handle.write('{}\t{}\n'.format(accession, table))
        return path

    def written_table(self, gpath):
        path = os.path.join(gpath, 'prodigal', 'prodigal_translation_table.tsv')
        with open(path) as handle:
            return handle.read()


# ------------------------------------------------------ reading the tables

class ReadTranslationTablesTests(TempDirCase):
    """A correction replaces a prediction; it is not weighed against it."""

    def test_predictions_are_read_from_the_summary(self):
        tables, sources = P.read_translation_tables(self.summary(('GCF_1.1', '25')))
        self.assertEqual(tables['GCF_1.1'], 25)
        self.assertEqual(sources['GCF_1.1'], P.SOURCE_PREDICTED)

    def test_an_override_replaces_the_prediction(self):
        tables, sources = P.read_translation_tables(
            self.summary(('GCF_1.1', '11')), self.override(('GCF_1.1', '4')))
        self.assertEqual(tables['GCF_1.1'], 4)
        self.assertEqual(sources['GCF_1.1'], P.SOURCE_OVERRIDE)

    def test_a_genome_not_corrected_keeps_its_prediction(self):
        tables, sources = P.read_translation_tables(
            self.summary(('GCF_1.1', '11'), ('GCF_2.1', '4')), self.override(('GCF_1.1', '25')))
        self.assertEqual(tables['GCF_2.1'], 4)
        self.assertEqual(sources['GCF_2.1'], P.SOURCE_PREDICTED)

    def test_a_correction_can_add_a_genome_the_summary_lacks(self):
        """A genome no classifier could answer for is callable by hand."""
        tables, sources = P.read_translation_tables(
            self.summary(('GCF_1.1', '11')), self.override(('GCF_9.1', '4')))
        self.assertEqual(tables['GCF_9.1'], 4)
        self.assertEqual(sources['GCF_9.1'], P.SOURCE_OVERRIDE)

    def test_a_gzipped_release_summary_reads_the_same_as_a_plain_one(self):
        """trans_table writes the release summary gzipped -- 116 MB of text for
        r237 against 34 MB -- and gTranslate writes a batch's plain."""
        path = os.path.join(self.dir, 'gtranslate.translation_table_summary.tsv.gz')
        with gzip.open(path, 'wt') as handle:
            handle.write('user_genome\tbest_tln_table\tconfidence\n')
            handle.write('GCF_1.1\t25\t1.0\n')
        tables, sources = P.read_translation_tables(path)
        self.assertEqual(tables['GCF_1.1'], 25)
        self.assertEqual(sources['GCF_1.1'], P.SOURCE_PREDICTED)

    def test_a_summary_is_read_by_what_it_is_and_not_by_what_it_is_called(self):
        """A release summary gunzipped, or renamed on the way to another machine,
        should not be the reason a release is called under the wrong tables."""
        path = os.path.join(self.dir, 'named_plain_but_gzipped.tsv')
        with gzip.open(path, 'wt') as handle:
            handle.write('user_genome\tbest_tln_table\n')
            handle.write('GCF_1.1\t4\n')
        tables, _ = P.read_translation_tables(path)
        self.assertEqual(tables['GCF_1.1'], 4)

    def test_summary_columns_are_read_by_name(self):
        """gTranslate has reordered its summary before; position must not matter."""
        path = os.path.join(self.dir, 'reordered.tsv')
        with open(path, 'w') as handle:
            handle.write('confidence\tbest_tln_table\tuser_genome\n')
            handle.write('1.0\t25\tGCF_1.1\n')
        tables, _ = P.read_translation_tables(path)
        self.assertEqual(tables['GCF_1.1'], 25)


# ------------------------------------------------------ a release not covered

class MissingTableTests(TempDirCase):
    """A genome with nothing to call its genes under is left, not guessed at."""

    def test_a_genome_without_a_table_is_not_called(self):
        """Prodigal choosing a table by coding density is the very thing handing
        it the summary is there to prevent."""
        one, two = self.genome('GCF_000000001.1'), self.genome('GCF_000000002.1')
        P.ProdigalManager(self.tmp_dir).run(
            self.genome_dirs(one, two), self.summary(('GCF_000000001.1', '11')),
                self.out_dir)
        self.assertNotIn('GCF_000000002.1', StubProdigal.genomes_called())

    def test_the_genomes_that_do_have_a_table_are_still_called(self):
        """Eight genomes of r237 have no prediction; the other 1.35M should not
        wait on them."""
        one, two = self.genome('GCF_000000001.1'), self.genome('GCF_000000002.1')
        P.ProdigalManager(self.tmp_dir).run(
            self.genome_dirs(one, two), self.summary(('GCF_000000001.1', '11')),
                self.out_dir)
        self.assertEqual(StubProdigal.genomes_called(), ['GCF_000000001.1'])
        self.assertEqual(StubProdigal.table_for('GCF_000000001.1'), 11)

    def test_a_summary_covering_no_genome_of_the_release_stops_the_run(self):
        """That is the wrong file, not a release with a few unpredictable
        genomes, and carrying on would call nothing and report that it had
        finished."""
        one = self.genome('GCF_000000001.1')
        with self.assertRaises(RuntimeError) as raised:
            P.ProdigalManager(self.tmp_dir).run(
                self.genome_dirs(one), self.summary(('GCF_000000009.1', '11')),
                self.out_dir)
        self.assertIn('another release', str(raised.exception))
        self.assertIsNone(StubProdigal.last_tasks)

    def test_the_run_reports_that_it_finished(self):
        """The genomes with no table are a known handful, not a failure of the
        run, and the 1.35M that were called are called."""
        one, two = self.genome('GCF_000000001.1'), self.genome('GCF_000000002.1')
        self.assertTrue(P.ProdigalManager(self.tmp_dir).run(
            self.genome_dirs(one, two), self.summary(('GCF_000000001.1', '11')),
                self.out_dir))

    def test_the_genomes_left_uncalled_are_named_in_a_file_of_the_release(self):
        """Silently calling fewer genomes than the release holds is how a gap
        reaches the next command. A batch may leave thousands, so the log says how
        many and of what kind and the file says which."""
        one, two = self.genome('GCF_000000001.1'), self.genome('GCF_000000002.1')
        with self.assertLogs('timestamp', level='WARNING') as caught:
            P.ProdigalManager(self.tmp_dir).run(
                self.genome_dirs(one, two), self.summary(('GCF_000000001.1', '11')),
                self.out_dir)

        warning = '\n'.join(caught.output)
        self.assertIn(P.REASON_NO_TABLE, warning)
        self.assertIn(P.NOT_CALLED_RELEASE_NAME, warning)

        with open(os.path.join(self.out_dir, P.NOT_CALLED_RELEASE_NAME)) as handle:
            rows = handle.read().splitlines()
        self.assertEqual(rows[0].split('\t'), list(P.NOT_CALLED_HEADER))
        self.assertEqual(rows[1:], ['GCF_000000002.1\t' + P.REASON_NO_TABLE])

    def test_a_genome_left_uncalled_has_no_results_written_for_it(self):
        one, two = self.genome('GCF_000000001.1'), self.genome('GCF_000000002.1')
        P.ProdigalManager(self.tmp_dir).run(
            self.genome_dirs(one, two), self.summary(('GCF_000000001.1', '11')),
                self.out_dir)
        self.assertFalse(os.path.exists(os.path.join(two[1], 'prodigal')))

    def test_an_override_can_fill_the_gap(self):
        one = self.genome('GCF_000000001.1')
        P.ProdigalManager(self.tmp_dir).run(
            self.genome_dirs(one),
            self.summary(('GCF_000000009.1', '11')),
            self.out_dir,
            self.override(('GCF_000000001.1', '4')))
        self.assertEqual(StubProdigal.table_for('GCF_000000001.1'), 4)


# ------------------------------------------------------ what Prodigal is told

class TranslationTableHandoverTests(TempDirCase):
    """Prodigal indexes this dictionary for every genome it is handed."""

    def test_the_predicted_table_is_passed_through(self):
        one = self.genome('GCF_000000001.1')
        P.ProdigalManager(self.tmp_dir).run(
            self.genome_dirs(one), self.summary(('GCF_000000001.1', '25')),
                self.out_dir)
        self.assertEqual(StubProdigal.table_for('GCF_000000001.1'), 25)

    def test_a_corrected_table_is_what_prodigal_is_given(self):
        one = self.genome('GCF_000000001.1')
        P.ProdigalManager(self.tmp_dir).run(
            self.genome_dirs(one),
            self.summary(('GCF_000000001.1', '11')),
            self.out_dir,
            self.override(('GCF_000000001.1', '4')))
        self.assertEqual(StubProdigal.table_for('GCF_000000001.1'), 4)


# ------------------------------------------------------ what each genome records

class TranslationTableFileTests(TempDirCase):
    """The file is written per genome, so it must describe THAT genome."""

    def test_a_predicted_genome_says_so(self):
        one = self.genome('GCF_000000001.1')
        P.ProdigalManager(self.tmp_dir).run(
            self.genome_dirs(one), self.summary(('GCF_000000001.1', '11')),
                self.out_dir)
        written = self.written_table(one[1])
        self.assertIn('best_translation_table\t11', written)
        self.assertIn(P.SOURCE_PREDICTED, written)

    def test_a_corrected_genome_says_it_was_corrected(self):
        one = self.genome('GCF_000000001.1')
        P.ProdigalManager(self.tmp_dir).run(
            self.genome_dirs(one),
            self.summary(('GCF_000000001.1', '11')),
            self.out_dir,
            self.override(('GCF_000000001.1', '4')))
        written = self.written_table(one[1])
        self.assertIn('best_translation_table\t4', written)
        self.assertIn(P.SOURCE_OVERRIDE, written)

    def test_one_genome_does_not_decide_what_is_written_for_another(self):
        one, two = self.genome('GCF_000000001.1'), self.genome('GCF_000000002.1')
        P.ProdigalManager(self.tmp_dir).run(
            self.genome_dirs(one, two),
            self.summary(('GCF_000000001.1', '11'), ('GCF_000000002.1', '4')),
            self.out_dir,
            self.override(('GCF_000000002.1', '25')))
        self.assertIn(P.SOURCE_PREDICTED, self.written_table(one[1]))
        self.assertIn(P.SOURCE_OVERRIDE, self.written_table(two[1]))

    def test_coding_densities_are_no_longer_written(self):
        """Only one table is run now, so the wrapper never measures the other."""
        one = self.genome('GCF_000000001.1')
        P.ProdigalManager(self.tmp_dir).run(
            self.genome_dirs(one), self.summary(('GCF_000000001.1', '11')),
                self.out_dir)
        self.assertNotIn('coding_density', self.written_table(one[1]))


if __name__ == '__main__':
    unittest.main()


# ------------------------------------------------------ what the run reports

class RunCountsTests(TempDirCase):
    """How much of a run is work and how much was already done."""

    def test_a_first_run_reports_every_genome_as_requiring_calling(self):
        one, two = self.genome('GCF_000000001.1'), self.genome('GCF_000000002.1')
        with self.assertLogs('timestamp', level='INFO') as captured:
            P.ProdigalManager(self.tmp_dir).run(
                self.genome_dirs(one, two),
                self.summary(('GCF_000000001.1', '11'), ('GCF_000000002.1', '4')),
                self.out_dir)
        self.assertIn('2 genome(s) require gene calling; 0 already have valid',
                      '\n'.join(captured.output))

    def test_a_rerun_skips_the_batch_rather_than_rechecking_its_genomes(self):
        """A finished batch is not looked into again. That is what batching buys
        over the per-genome checksum alone: a rerun does not re-read the proteins
        of 1.35M genomes to find out they are still there."""
        one, two = self.genome('GCF_000000001.1'), self.genome('GCF_000000002.1')
        paths = self.genome_dirs(one, two)
        summary = self.summary(('GCF_000000001.1', '11'), ('GCF_000000002.1', '4'))
        P.ProdigalManager(self.tmp_dir).run(paths, summary, self.out_dir)

        StubProdigal.last_tasks = None
        with self.assertLogs('timestamp', level='INFO') as captured:
            P.ProdigalManager(self.tmp_dir).run(paths, summary, self.out_dir)

        self.assertIn('already finished, skipping', '\n'.join(captured.output))
        self.assertIsNone(StubProdigal.last_tasks)

    def test_a_genome_already_called_is_skipped_within_an_unfinished_batch(self):
        """The checksum still decides genome by genome, which is what a batch
        retried after a failure leans on."""
        one, two = self.genome('GCF_000000001.1'), self.genome('GCF_000000002.1')
        paths = self.genome_dirs(one, two)
        summary = self.summary(('GCF_000000001.1', '11'), ('GCF_000000002.1', '4'))
        P.ProdigalManager(self.tmp_dir).run(paths, summary, self.out_dir)

        # as a batch left FAILED by an earlier run would be found
        batch = os.path.join(self.out_dir, 'batch_000001')
        os.unlink(os.path.join(batch, 'SUCCESS'))

        with self.assertLogs('timestamp', level='INFO') as captured:
            P.ProdigalManager(self.tmp_dir).run(paths, summary, self.out_dir)
        self.assertIn('0 genome(s) require gene calling; 2 already have valid',
                      '\n'.join(captured.output))


# ------------------------------------------ batches, and sharing them between machines

class BatchingTests(TempDirCase):
    """--out_dir holds the state of the run; the genes go to the genome directories."""

    def release(self, *accessions):
        genomes = [self.genome(a) for a in accessions]
        return (self.genome_dirs(*genomes),
                self.summary(*[(a, '11') for a in accessions]),
                genomes)

    def manager(self, **kwargs):
        return P.ProdigalManager(self.tmp_dir, batch_size=kwargs.pop('batch_size', 1),
                                 **kwargs)

    def test_the_release_is_cut_into_batches_under_the_output_directory(self):
        paths, summary, _ = self.release('GCF_000000001.1', 'GCF_000000002.1')
        self.manager().run(paths, summary, self.out_dir)
        self.assertEqual(sorted(d for d in os.listdir(self.out_dir)
                                if d.startswith('batch_')),
                         ['batch_000001', 'batch_000002'])

    def test_no_called_genes_are_written_to_the_output_directory(self):
        """The genes go where the rest of the toolkit looks for them, which is
        the genome's own directory, and that is why two machines on different
        batches never write to the same place."""
        paths, summary, genomes = self.release('GCF_000000001.1')
        self.manager().run(paths, summary, self.out_dir)

        written = []
        for root, _, files in os.walk(self.out_dir):
            written.extend(files)
        self.assertEqual([f for f in written if f.endswith('.faa.gz')], [])
        self.assertTrue(os.path.exists(os.path.join(
            genomes[0][1], 'prodigal', 'GCF_000000001.1_protein.faa.gz')))

    def test_a_finished_batch_carries_its_canary(self):
        paths, summary, _ = self.release('GCF_000000001.1')
        self.manager().run(paths, summary, self.out_dir)
        batch = os.path.join(self.out_dir, 'batch_000001')
        self.assertEqual(B.batch_state(batch), B.STATE_SUCCESS)
        self.assertFalse(os.path.exists(os.path.join(batch, B.RUNNING_CANARY)))

    def test_a_batch_another_machine_holds_is_left_alone(self):
        """This is what lets several machines share one --out_dir: each takes the
        batches no other machine holds."""
        paths, summary, _ = self.release('GCF_000000001.1', 'GCF_000000002.1')
        manager = self.manager()
        batches = manager.plan_batches(paths, self.out_dir, {'GCF_000000001.1': 11})

        # a claim held by a live process, which is what another machine still
        # working on the batch looks like from here
        B.claim_batch(batches[0])
        manager.run(paths, summary, self.out_dir)

        self.assertEqual(B.batch_state(batches[0]), B.STATE_RUNNING)
        self.assertEqual(B.batch_state(batches[1]), B.STATE_SUCCESS)

    def test_a_batch_that_failed_is_marked_and_the_run_says_so(self):
        """A batch lost is a batch, not the release."""
        paths, summary, _ = self.release('GCF_000000001.1', 'GCF_000000002.1')
        with mock.patch.object(P.ProdigalManager, 'call_batch',
                               side_effect=RuntimeError('prodigal fell over')):
            self.assertFalse(self.manager().run(paths, summary, self.out_dir))
        batch = os.path.join(self.out_dir, 'batch_000001')
        self.assertEqual(B.batch_state(batch), B.STATE_FAILED)

    def test_a_failed_batch_is_retried_by_the_next_run(self):
        paths, summary, _ = self.release('GCF_000000001.1')
        with mock.patch.object(P.ProdigalManager, 'call_batch',
                               side_effect=RuntimeError('prodigal fell over')):
            self.manager().run(paths, summary, self.out_dir)

        self.assertTrue(self.manager().run(paths, summary, self.out_dir))
        self.assertEqual(B.batch_state(os.path.join(self.out_dir, 'batch_000001')),
                         B.STATE_SUCCESS)

    def test_an_interrupted_batch_is_given_back_rather_than_left_claimed(self):
        """Ctrl-C is not a failure of the batch, and the next run should not wait
        out a lease on a machine that has already stopped."""
        paths, summary, _ = self.release('GCF_000000001.1')
        with mock.patch.object(P.ProdigalManager, 'call_batch',
                               side_effect=KeyboardInterrupt):
            self.assertRaises(KeyboardInterrupt,
                              self.manager().run, paths, summary, self.out_dir)
        self.assertEqual(B.batch_state(os.path.join(self.out_dir, 'batch_000001')),
                         B.STATE_PENDING)

    def test_what_happened_to_a_batch_is_written_in_the_batch(self):
        """Five machines appending to one --log over NFS leave none of it."""
        paths, summary, _ = self.release('GCF_000000001.1')
        # assertLogs lifts the logger to INFO, which logger_setup() does in a run
        with self.assertLogs('timestamp', level='INFO'):
            self.manager().run(paths, summary, self.out_dir)
        log = os.path.join(self.out_dir, 'batch_000001', P.BATCH_LOG_NAME)
        self.assertIn('starting', open(log).read())

    def test_the_plan_is_reused_rather_than_cut_again(self):
        """Partitioning a release again that has gained a genome would move
        genomes between batches that are already finished."""
        paths, summary, _ = self.release('GCF_000000001.1', 'GCF_000000002.1')
        self.manager().run(paths, summary, self.out_dir)
        before = sorted(os.listdir(os.path.join(self.out_dir, 'batch_000001')))

        self.manager(batch_size=10).run(paths, summary, self.out_dir)
        self.assertEqual(sorted(d for d in os.listdir(self.out_dir)
                                if d.startswith('batch_')),
                         ['batch_000001', 'batch_000002'])
        self.assertEqual(sorted(os.listdir(os.path.join(self.out_dir, 'batch_000001'))),
                         before)

    def test_all_genomes_does_the_finished_batches_again(self):
        """Otherwise it would discard nothing and call nothing, having skipped
        every batch it was asked to redo."""
        paths, summary, _ = self.release('GCF_000000001.1')
        self.manager().run(paths, summary, self.out_dir)

        StubProdigal.last_tasks = None
        self.manager().run(paths, summary, self.out_dir, all_genomes=True)
        self.assertEqual(StubProdigal.genomes_called(), ['GCF_000000001.1'])

    def test_the_release_file_waits_for_every_batch(self):
        """The file at the top is the whole release or absent, never a part of it
        that reads like the whole."""
        paths, summary, _ = self.release('GCF_000000001.1', 'GCF_000000002.1')
        manager = self.manager()
        batches = manager.plan_batches(paths, self.out_dir, {'GCF_000000001.1': 11})
        B.claim_batch(batches[1])
        manager.run(paths, summary, self.out_dir)
        self.assertFalse(os.path.exists(
            os.path.join(self.out_dir, P.NOT_CALLED_RELEASE_NAME)))

    def test_a_genome_with_no_fasta_is_recorded_rather_than_called(self):
        one = self.genome('GCF_000000001.1')
        two = ('GCF_000000002.1', os.path.join(self.dir, 'GCF_000000002.1_ASM1'))
        os.makedirs(two[1])
        paths = self.genome_dirs(one, two)
        summary = self.summary(('GCF_000000001.1', '11'), ('GCF_000000002.1', '11'))

        self.manager().run(paths, summary, self.out_dir)

        with open(os.path.join(self.out_dir, P.NOT_CALLED_RELEASE_NAME)) as handle:
            rows = handle.read().splitlines()[1:]
        self.assertEqual(rows, ['GCF_000000002.1\t' + P.REASON_NO_FASTA])
