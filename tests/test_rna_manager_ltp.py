#!/usr/bin/env python3
"""Offline unit tests for rna_manager_ltp.py -- blastn is never run.

Run with the interpreter that has tqdm:

    /opt/centos7/sw/miniconda3/envs/gtdb_migration_tk-r237/bin/python -m unittest discover -s tests -p test_rna_manager_ltp.py

genometk_lite's RNA is replaced by a stub whose classify() writes the files the
real one writes, so what is tested here is the bookkeeping: which genomes are
classified, which have nothing to classify, which are named as unclassified and
why, what lands in the genome directory and in what order, and what a batch
records about itself. The stub is a module-level class because the workers run
in forked processes.
"""

import logging
import os
import pickle
import shutil
import tempfile
import unittest
from unittest import mock

from gtdb_migration_tk import batching as B
from gtdb_migration_tk import config
from gtdb_migration_tk import rna_manager_ltp as L
from gtdb_migration_tk import rna_manager_silva as S
from gtdb_migration_tk.biolib_lite.external.blast import BlastError
from gtdb_migration_tk.utils import common as C


# A 16S gene whose sequence says this is one the stubbed blastn fails on.
REFUSE = 'REFUSE'

BLASTN_VERSION = 'blastn: 2.16.0+'
LTP_VERSION = '10_2024'
SSU_VERSION = '138.2'


class StubRNA(object):
    """Stands in for genometk_lite's RNA, writing what classify() writes."""

    def __init__(self, rna_name, domain, cpus):
        self.rna_name = rna_name

    def classify(self, seq_file, db, taxonomy_file, output_dir):
        with open(seq_file) as handle:
            if REFUSE in handle.read():
                raise BlastError('blastn failed with status 2: database error')
        with open(os.path.join(output_dir, 'ssu.blastn.tsv'), 'w') as handle:
            handle.write('contig_1\t1500\tLTP_1\n')
        with open(os.path.join(output_dir, 'ssu.taxonomy.tsv'), 'w') as handle:
            handle.write('query_id\ttaxonomy\ncontig_1\td__Bacteria\n')


class QuietTqdm(object):
    def __init__(self, iterable=None, **kwargs):
        self.iterable = iterable

    def __iter__(self):
        return iter(() if self.iterable is None else self.iterable)

    def __enter__(self):
        return self

    def __exit__(self, *exc_info):
        return False

    def update(self, n=1):
        pass

    def close(self):
        pass


class InProcessPool(object):
    """Runs the pool's work in this process, weighing each task if asked to."""

    sizes = None

    def __init__(self, processes=None):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def imap_unordered(self, func, items):
        for item in items:
            if InProcessPool.sizes is not None:
                InProcessPool.sizes.append(len(pickle.dumps((func, item))))
            yield func(item)


class TempDirCase(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='rna_ltp_test.')
        self.out_dir = os.path.join(self.dir, 'out')

        for module in (L, B, C):
            patch = mock.patch.object(module, 'tqdm', QuietTqdm)
            patch.start()
            self.addCleanup(patch.stop)

        for name, value in (('record_program_version', mock.Mock(return_value=BLASTN_VERSION)),
                            ('check_dependencies', mock.Mock()),
                            ('RNA', StubRNA)):
            patch = mock.patch.object(L, name, value)
            patch.start()
            self.addCleanup(patch.stop)

        logger = logging.getLogger('timestamp')
        self.addCleanup(logger.setLevel, logger.level)
        logger.setLevel(logging.CRITICAL)

        self.ltp_dir = os.path.join(self.dir, 'ltp')
        release = os.path.join(self.ltp_dir, LTP_VERSION)
        os.makedirs(release)
        for name in ('ltp_{}.fna', 'ltp_{}_taxonomy.tsv'):
            with open(os.path.join(release, name.format(LTP_VERSION)), 'w') as handle:
                handle.write('x\n')

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    # ------------------------------------------------------------- the release

    def genome_dir(self, accession, ssu='gene', classified=False):
        """A genome directory after rna_silva.

        ssu is 'gene' where rna_silva found a 16S gene, 'refuse' for one blastn
        fails on, 'none' where it searched and found none, and None where it has
        not searched the genome yet.
        """
        gpath = os.path.join(self.dir, '{}_ASM{}v1'.format(accession, accession[4:10]))
        silva = os.path.join(gpath, L.SSU_RESULTS_DIR_FORMAT.format(SSU_VERSION))
        os.makedirs(gpath)

        if ssu is not None:
            os.makedirs(silva)
            with open(os.path.join(silva, L.SSU_CANARY), 'w') as handle:
                handle.write('done.\n')
            if ssu != 'none':
                with open(os.path.join(silva, L.SSU_FASTA), 'w') as handle:
                    handle.write('>contig_1\n{}\n'.format(REFUSE if ssu == 'refuse' else 'ACGT'))

        if classified:
            ltp = self.results(gpath)
            os.makedirs(ltp)
            for name in (L.LTP_CANARY, 'ssu.taxonomy.tsv', 'ssu.stale_from_an_old_run.tsv'):
                with open(os.path.join(ltp, name), 'w') as handle:
                    handle.write('from an earlier run\n')

        return gpath

    def genome_dirs_file(self, genomes):
        path = os.path.join(self.dir, 'genome_dirs.tsv')
        with open(path, 'w') as handle:
            for accession, gpath in genomes:
                handle.write('{}\t{}\t{}\n'.format(accession, gpath, 'G' + accession[4:-2]))
        return path

    def manager(self, **kwargs):
        return L.RnaManagerLTP(LTP_VERSION, SSU_VERSION, self.ltp_dir, cpus=1,
                               tmp_dir=os.path.join(self.dir, 'tmp'), **kwargs)

    def run_ltp(self, genomes, all_genomes=False, remove_prior=False,
                batch_size=B.DEFAULT_BATCH_SIZE):
        manager = self.manager(batch_size=batch_size)
        ok = manager.run(self.genome_dirs_file(genomes), self.out_dir,
                         all_genomes, remove_prior)
        return manager, ok

    # ------------------------------------------------------------- reading back

    def results(self, gpath):
        return os.path.join(gpath, L.LTP_RESULTS_DIR_FORMAT.format(LTP_VERSION))

    def classified(self, gpath):
        """Whether this run classified the genome."""
        path = os.path.join(self.results(gpath), 'ssu.blastn.tsv')
        return os.path.exists(path)

    def run_dir(self):
        return os.path.join(self.out_dir, 'rna_ltp_{}-silva_{}'.format(LTP_VERSION, SSU_VERSION))

    def batches(self):
        return B.batch_dir_names(self.run_dir(), L.LAYOUT)

    def release_report(self):
        with open(os.path.join(self.run_dir(), L.NOT_CLASSIFIED_RELEASE_NAME)) as handle:
            header = handle.readline().rstrip('\n').split('\t')
            self.assertEqual(tuple(header), L.NOT_CLASSIFIED_HEADER)
            return sorted(tuple(line.rstrip('\n').split('\t')) for line in handle)


# ------------------------------------------------------- what a genome has to classify

class WhatAGenomeHasToClassify(TempDirCase):
    """rna_silva's ssu.fna, and rna_silva's canary for why there is none."""

    def test_a_genome_with_a_16s_gene_is_classified(self):
        gpath = self.genome_dir('GCA_000001.1')

        self.run_ltp([('GCA_000001.1', gpath)])

        self.assertTrue(self.classified(gpath))
        self.assertEqual(self.release_report(), [])

    def test_a_genome_rna_silva_found_no_gene_in_is_counted_not_named(self):
        """Nothing is wrong with it; there is nothing to classify."""
        gpath = self.genome_dir('GCA_000001.1', ssu='none')

        _, ok = self.run_ltp([('GCA_000001.1', gpath)])

        self.assertTrue(ok)
        self.assertEqual(self.release_report(), [])
        canary = B.read_canary(os.path.join(self.batches()[0], B.SUCCESS_CANARY))
        self.assertEqual(canary['no_ssu_gene'], '1')

    def test_a_genome_rna_silva_has_not_searched_is_named(self):
        """Or the order the two commands ran in reads as genomes with no 16S."""
        gpath = self.genome_dir('GCA_000001.1', ssu=None)

        self.run_ltp([('GCA_000001.1', gpath)])

        self.assertEqual(self.release_report(),
                         [('GCA_000001.1', L.REASON_SSU_NOT_IDENTIFIED)])

    def test_the_batchfile_names_the_file_classified_not_the_genome(self):
        gpath = self.genome_dir('GCA_000001.1')

        self.run_ltp([('GCA_000001.1', gpath)])

        rows = B.read_batchfile(B.batchfile_path(self.batches()[0], L.LAYOUT))
        self.assertEqual(rows, [(L.ssu_fasta(SSU_VERSION, 'GCA_000001.1', gpath),
                                 'GCA_000001.1')])

    def test_no_domain_file_is_asked_for(self):
        """classify() never read the domain the command looked up for it."""
        self.assertNotIn('gtdb_domain_file',
                         L.RnaManagerLTP.__init__.__code__.co_varnames)


class WhereTheGenesAreRead(unittest.TestCase):
    """From where rna_silva writes them, which is the directory the release
    carries across -- until 0.1.36 rna_silva wrote to a -testing directory and
    this command classified the old workflow's genes."""

    def test_it_reads_the_directory_rna_silva_writes(self):
        self.assertEqual(L.SSU_RESULTS_DIR_FORMAT, S.RESULTS_DIR_FORMAT)

    def test_that_is_the_directory_the_release_carries_across(self):
        self.assertIn(S.RESULTS_DIR_FORMAT.format(config.SILVA_VERSION),
                      config.GTDB_DERIVED_DIRS_TO_COPY)


# ---------------------------------------------------------- what decides the work

class WhatDecidesTheWork(TempDirCase):

    def test_a_genome_with_an_ltp_canary_is_skipped(self):
        gpath = self.genome_dir('GCA_000001.1', classified=True)

        self.run_ltp([('GCA_000001.1', gpath)])

        self.assertFalse(self.classified(gpath))

    def test_all_classifies_it_again(self):
        gpath = self.genome_dir('GCA_000001.1', classified=True)

        self.run_ltp([('GCA_000001.1', gpath)], all_genomes=True)

        self.assertTrue(self.classified(gpath))

    def test_a_genome_classified_again_keeps_nothing_of_its_last_classification(self):
        gpath = self.genome_dir('GCA_000001.1', classified=True)

        self.run_ltp([('GCA_000001.1', gpath)], all_genomes=True)

        self.assertFalse(os.path.exists(os.path.join(
            self.results(gpath), 'ssu.stale_from_an_old_run.tsv')))

    def test_the_canary_records_both_versions(self):
        gpath = self.genome_dir('GCA_000001.1')

        self.run_ltp([('GCA_000001.1', gpath)])

        with open(os.path.join(self.results(gpath), L.LTP_CANARY)) as handle:
            text = handle.read()
        self.assertIn('Silva version:{}.'.format(SSU_VERSION), text)
        self.assertIn('LTP version:{}.'.format(LTP_VERSION), text)

    def test_remove_does_not_fall_over_on_a_genome_never_classified(self):
        gpath = self.genome_dir('GCA_000001.1')

        _, ok = self.run_ltp([('GCA_000001.1', gpath)], remove_prior=True)

        self.assertTrue(ok)
        self.assertTrue(self.classified(gpath))

    def test_the_canary_is_copied_after_everything_it_vouches_for(self):
        gpath = self.genome_dir('GCA_000001.1')
        copied = []
        real_copy = shutil.copy

        def recording_copy(source, target):
            copied.append(os.path.basename(target))
            return real_copy(source, target)

        manager = self.manager()
        with mock.patch.object(L.mp, 'Pool', InProcessPool), \
                mock.patch.object(L.shutil, 'copy', recording_copy):
            manager.run(self.genome_dirs_file([('GCA_000001.1', gpath)]), self.out_dir)

        self.assertEqual(copied.index(L.LTP_CANARY), len(copied) - 1)
        self.assertIn('ssu.taxonomy.tsv', copied)


class WhatThePoolIsHanded(TempDirCase):

    MOST_A_TASK_SHOULD_WEIGH = 10000

    def test_a_task_is_the_genome_and_not_the_instance(self):
        genomes = [('GCA_00000{}.1'.format(i),
                    self.genome_dir('GCA_00000{}.1'.format(i))) for i in (1, 2)]
        manager = self.manager()
        manager.ballast = ['x' * 100] * 1000

        InProcessPool.sizes = []
        self.addCleanup(setattr, InProcessPool, 'sizes', None)
        with mock.patch.object(L.mp, 'Pool', InProcessPool):
            ok = manager.run(self.genome_dirs_file(genomes), self.out_dir)

        self.assertTrue(ok)
        self.assertEqual(len(InProcessPool.sizes), 4)
        self.assertLess(max(InProcessPool.sizes), self.MOST_A_TASK_SHOULD_WEIGH)


# ------------------------------------------------------------- blastn failing

class WhenBlastnFails(TempDirCase):
    """BLAST's wrapper raises now, where it returned as though it had succeeded."""

    def test_the_genome_is_named_and_the_batch_still_succeeds(self):
        whole = self.genome_dir('GCA_000001.1')
        bad = self.genome_dir('GCA_000004.1', ssu='refuse')

        _, ok = self.run_ltp([('GCA_000001.1', whole), ('GCA_000004.1', bad)])

        self.assertTrue(ok)
        self.assertEqual(self.release_report(),
                         [('GCA_000004.1', L.REASON_BLASTN_FAILED)])

    def test_a_genome_blastn_failed_on_gets_no_canary(self):
        """Which is what wrote it down as classified with no hits."""
        whole = self.genome_dir('GCA_000001.1')
        bad = self.genome_dir('GCA_000004.1', ssu='refuse')

        self.run_ltp([('GCA_000001.1', whole), ('GCA_000004.1', bad)])

        self.assertFalse(os.path.exists(os.path.join(self.results(bad), L.LTP_CANARY)))
        self.assertFalse(os.path.exists(C.version_file(self.results(bad), L.BLASTN)))

    def test_a_batch_in_which_every_genome_failed_is_failed(self):
        bad = self.genome_dir('GCA_000004.1', ssu='refuse')
        worse = self.genome_dir('GCA_000005.1', ssu='refuse')

        _, ok = self.run_ltp([('GCA_000004.1', bad), ('GCA_000005.1', worse)])

        self.assertFalse(ok)
        self.assertEqual(B.batch_state(self.batches()[0]), B.STATE_FAILED)


# --------------------------------------------------------- batches and the release

class TheBatchesAndTheRelease(TempDirCase):

    def setUp(self):
        super().setUp()
        self.genomes = [(accession, self.genome_dir(accession))
                        for accession in ('GCA_000001.1', 'GCA_000002.1',
                                          'GCA_000003.1', 'GCA_000004.1')]

    def test_the_release_is_cut_into_batches_of_batch_size(self):
        self.run_ltp(self.genomes, batch_size=2)

        self.assertEqual(len(self.batches()), 2)

    def test_each_batch_records_what_it_came_to(self):
        self.run_ltp(self.genomes, batch_size=2)

        canary = B.read_canary(os.path.join(self.batches()[0], B.SUCCESS_CANARY))
        self.assertEqual(canary['classified'], '2')

    def test_a_scanned_genome_records_the_blastn_that_classified_it(self):
        self.run_ltp(self.genomes, batch_size=2)

        with open(C.version_file(self.results(self.genomes[0][1]), L.BLASTN)) as handle:
            self.assertEqual(handle.read().strip(), BLASTN_VERSION)

    def test_no_results_are_written_to_the_out_dir(self):
        self.run_ltp(self.genomes, batch_size=2)

        top = sorted(name for name in os.listdir(self.run_dir())
                     if not name.startswith(B.BATCH_DIR_PREFIX))
        self.assertEqual(top, [L.NOT_CLASSIFIED_RELEASE_NAME])

    def test_a_finished_batch_is_skipped_by_a_later_run(self):
        self.run_ltp(self.genomes, batch_size=2)
        shutil.rmtree(self.results(self.genomes[0][1]))

        self.run_ltp(self.genomes, batch_size=2)

        self.assertFalse(self.classified(self.genomes[0][1]))

    def test_a_batch_held_by_another_machine_is_left_alone(self):
        self.run_ltp(self.genomes, batch_size=2)
        batch = self.batches()[0]
        os.remove(os.path.join(batch, B.SUCCESS_CANARY))
        with open(os.path.join(batch, B.RUNNING_CANARY), 'w') as handle:
            handle.write(B.canary_payload(host='another-machine', pid='1'))

        self.run_ltp(self.genomes, batch_size=2)

        self.assertEqual(B.batch_state(batch), B.STATE_RUNNING)


if __name__ == '__main__':
    unittest.main()
