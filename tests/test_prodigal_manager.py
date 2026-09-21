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
            self.genome_dirs(one, two), self.summary(('GCF_000000001.1', '11')))
        self.assertNotIn('GCF_000000002.1', StubProdigal.genomes_called())

    def test_the_genomes_that_do_have_a_table_are_still_called(self):
        """Eight genomes of r237 have no prediction; the other 1.35M should not
        wait on them."""
        one, two = self.genome('GCF_000000001.1'), self.genome('GCF_000000002.1')
        P.ProdigalManager(self.tmp_dir).run(
            self.genome_dirs(one, two), self.summary(('GCF_000000001.1', '11')))
        self.assertEqual(StubProdigal.genomes_called(), ['GCF_000000001.1'])
        self.assertEqual(StubProdigal.table_for('GCF_000000001.1'), 11)

    def test_a_summary_covering_no_genome_of_the_release_stops_the_run(self):
        """That is the wrong file, not a release with a few unpredictable
        genomes, and carrying on would call nothing and report that it had
        finished."""
        one = self.genome('GCF_000000001.1')
        with self.assertRaises(RuntimeError) as raised:
            P.ProdigalManager(self.tmp_dir).run(
                self.genome_dirs(one), self.summary(('GCF_000000009.1', '11')))
        self.assertIn('another release', str(raised.exception))
        self.assertIsNone(StubProdigal.last_tasks)

    def test_the_run_reports_that_it_finished(self):
        """The genomes with no table are a known handful, not a failure of the
        run, and the 1.35M that were called are called."""
        one, two = self.genome('GCF_000000001.1'), self.genome('GCF_000000002.1')
        self.assertTrue(P.ProdigalManager(self.tmp_dir).run(
            self.genome_dirs(one, two), self.summary(('GCF_000000001.1', '11'))))

    def test_the_genomes_left_uncalled_are_named_in_the_log(self):
        """Silently calling fewer genomes than the release holds is how a gap
        reaches the next command."""
        one, two = self.genome('GCF_000000001.1'), self.genome('GCF_000000002.1')
        with self.assertLogs('timestamp', level='WARNING') as caught:
            P.ProdigalManager(self.tmp_dir).run(
                self.genome_dirs(one, two), self.summary(('GCF_000000001.1', '11')))
        warning = '\n'.join(caught.output)
        self.assertIn('GCF_000000002.1', warning)
        self.assertIn('gtranslate_no_prediction.tsv', warning)

    def test_a_genome_left_uncalled_has_no_results_written_for_it(self):
        one, two = self.genome('GCF_000000001.1'), self.genome('GCF_000000002.1')
        P.ProdigalManager(self.tmp_dir).run(
            self.genome_dirs(one, two), self.summary(('GCF_000000001.1', '11')))
        self.assertFalse(os.path.exists(os.path.join(two[1], 'prodigal')))

    def test_an_override_can_fill_the_gap(self):
        one = self.genome('GCF_000000001.1')
        P.ProdigalManager(self.tmp_dir).run(
            self.genome_dirs(one),
            self.summary(('GCF_000000009.1', '11')),
            self.override(('GCF_000000001.1', '4')))
        self.assertEqual(StubProdigal.table_for('GCF_000000001.1'), 4)


# ------------------------------------------------------ what Prodigal is told

class TranslationTableHandoverTests(TempDirCase):
    """Prodigal indexes this dictionary for every genome it is handed."""

    def test_the_predicted_table_is_passed_through(self):
        one = self.genome('GCF_000000001.1')
        P.ProdigalManager(self.tmp_dir).run(
            self.genome_dirs(one), self.summary(('GCF_000000001.1', '25')))
        self.assertEqual(StubProdigal.table_for('GCF_000000001.1'), 25)

    def test_a_corrected_table_is_what_prodigal_is_given(self):
        one = self.genome('GCF_000000001.1')
        P.ProdigalManager(self.tmp_dir).run(
            self.genome_dirs(one),
            self.summary(('GCF_000000001.1', '11')),
            self.override(('GCF_000000001.1', '4')))
        self.assertEqual(StubProdigal.table_for('GCF_000000001.1'), 4)


# ------------------------------------------------------ what each genome records

class TranslationTableFileTests(TempDirCase):
    """The file is written per genome, so it must describe THAT genome."""

    def test_a_predicted_genome_says_so(self):
        one = self.genome('GCF_000000001.1')
        P.ProdigalManager(self.tmp_dir).run(
            self.genome_dirs(one), self.summary(('GCF_000000001.1', '11')))
        written = self.written_table(one[1])
        self.assertIn('best_translation_table\t11', written)
        self.assertIn(P.SOURCE_PREDICTED, written)

    def test_a_corrected_genome_says_it_was_corrected(self):
        one = self.genome('GCF_000000001.1')
        P.ProdigalManager(self.tmp_dir).run(
            self.genome_dirs(one),
            self.summary(('GCF_000000001.1', '11')),
            self.override(('GCF_000000001.1', '4')))
        written = self.written_table(one[1])
        self.assertIn('best_translation_table\t4', written)
        self.assertIn(P.SOURCE_OVERRIDE, written)

    def test_one_genome_does_not_decide_what_is_written_for_another(self):
        one, two = self.genome('GCF_000000001.1'), self.genome('GCF_000000002.1')
        P.ProdigalManager(self.tmp_dir).run(
            self.genome_dirs(one, two),
            self.summary(('GCF_000000001.1', '11'), ('GCF_000000002.1', '4')),
            self.override(('GCF_000000002.1', '25')))
        self.assertIn(P.SOURCE_PREDICTED, self.written_table(one[1]))
        self.assertIn(P.SOURCE_OVERRIDE, self.written_table(two[1]))

    def test_coding_densities_are_no_longer_written(self):
        """Only one table is run now, so the wrapper never measures the other."""
        one = self.genome('GCF_000000001.1')
        P.ProdigalManager(self.tmp_dir).run(
            self.genome_dirs(one), self.summary(('GCF_000000001.1', '11')))
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
                self.summary(('GCF_000000001.1', '11'), ('GCF_000000002.1', '4')))
        self.assertIn('2 genome(s) require gene calling; 0 already have valid',
                      '\n'.join(captured.output))

    def test_a_rerun_reports_the_genomes_it_skipped(self):
        """The checksum written by the first run is what vouches for the second."""
        one, two = self.genome('GCF_000000001.1'), self.genome('GCF_000000002.1')
        paths = self.genome_dirs(one, two)
        summary = self.summary(('GCF_000000001.1', '11'), ('GCF_000000002.1', '4'))
        P.ProdigalManager(self.tmp_dir).run(paths, summary)
        with self.assertLogs('timestamp', level='INFO') as captured:
            P.ProdigalManager(self.tmp_dir).run(paths, summary)
        self.assertIn('0 genome(s) require gene calling; 2 already have valid',
                      '\n'.join(captured.output))
