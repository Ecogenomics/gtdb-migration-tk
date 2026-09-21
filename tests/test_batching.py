#!/usr/bin/env python3
"""Offline unit tests for batching.py -- the batches of a release, and sharing
them between machines.

What is tested here is the mechanism and not either command's use of it: cutting
a release into batches, the batchfiles, who holds which batch and for how long,
and where what happened to a batch is written. trans_table predicting a table and
prodigal calling genes are those commands' business; two machines both believing
they hold one batch, or a release being partitioned again with its batches
already finished, is this module's.

The layout used here names no real command. That is the point: the mechanism is
told which files a batch keeps, and a test that passed only under trans_table's
names would be testing trans_table.
"""

import logging
import os
import shutil
import tempfile
import time
import unittest
from unittest import mock

from gtdb_migration_tk import batching as B
from gtdb_migration_tk.ncbi_utils import genomic_fasta


# A command of no name at all, so that nothing here leans on what trans_table or
# prodigal happen to call their files. The older batchfile is in the layout
# because a directory planned by an earlier version must still read as planned.
BATCHFILE_NAME = 'example_batchfile.tsv.gz'
LEGACY_BATCHFILE_NAME = 'example_batchfile.tsv'
BATCH_LOG_NAME = 'example.log'
LAYOUT = B.BatchLayout(batchfiles=(BATCHFILE_NAME, LEGACY_BATCHFILE_NAME),
                       log=BATCH_LOG_NAME)


class TempDirCase(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='batching_test.')

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


class LayoutTests(unittest.TestCase):
    """What differs between the commands is a name, and travels in the layout."""

    def test_a_command_names_the_files_its_batches_keep(self):
        from gtdb_migration_tk import prodigal_manager, trans_table
        self.assertNotEqual(trans_table.LAYOUT, prodigal_manager.LAYOUT)
        for layout in (trans_table.LAYOUT, prodigal_manager.LAYOUT):
            self.assertTrue(layout.batchfiles)
            self.assertTrue(layout.log)

    def test_the_batchfiles_are_tried_in_the_order_given(self):
        """The first is what a new batch is written under; the rest are names an
        earlier version used and are still read."""
        self.assertEqual(LAYOUT.batchfiles[0], BATCHFILE_NAME)


class ReadGenomeDirsTests(TempDirCase):
    """The genome_dirs file is the lingua franca; extra columns may be appended."""

    def write_genome_dirs(self, text):
        path = os.path.join(self.dir, 'genome_dirs.tsv')
        with open(path, 'w') as handle:
            handle.write(text)
        return path

    def test_accession_and_directory_are_read(self):
        path = self.write_genome_dirs(
            'GCF_000001405.39\t/rel/refseq/GCF/000/001/405/x\tG000001405\n')
        self.assertEqual(B.read_genome_dirs(path),
                         [('GCF_000001405.39', '/rel/refseq/GCF/000/001/405/x')])

    def test_further_columns_are_ignored(self):
        path = self.write_genome_dirs('GCA_1.1\t/rel/x\tG1\tsomething\telse\n')
        self.assertEqual(B.read_genome_dirs(path), [('GCA_1.1', '/rel/x')])

    def test_blank_lines_are_skipped(self):
        path = self.write_genome_dirs('GCA_1.1\t/rel/x\tG1\n\n')
        self.assertEqual(len(B.read_genome_dirs(path)), 1)


class FastaSizeTests(TempDirCase):
    """One stat, where two calls asked the file server the same question twice."""

    def test_the_size_of_a_fasta_is_returned(self):
        path = genomic_fasta(self.genome_dir('GCF_1.1_ASM1'))
        self.assertEqual(B.fasta_size(path), os.path.getsize(path))

    def test_a_file_that_is_not_there_has_no_size_rather_than_raising(self):
        self.assertEqual(B.fasta_size(os.path.join(self.dir, 'no_such.fna.gz')), 0)

    def test_an_empty_fasta_has_no_size(self):
        path = genomic_fasta(self.genome_dir('GCA_2.1_ASM2', empty=True))
        self.assertEqual(B.fasta_size(path), 0)


class SplitByFastaTests(TempDirCase):
    """A genome the release does not hold is left out rather than taken along.

    gTranslate refuses a whole batch over one missing path; prodigal records the
    genome as having no FASTA and calls the rest. Neither wants it in the list.
    """

    def rows(self, *genomes):
        """(FASTA path, accession) as a batchfile names them."""
        return [(genomic_fasta(path), accession) for accession, path in genomes]

    def test_genome_with_a_fasta_is_kept(self):
        rows = self.rows(('GCF_000001405.39',
                          self.genome_dir('GCF_000001405.39_GRCh38.p13')))
        present, missing = B.split_by_fasta(rows)
        self.assertEqual(missing, [])
        self.assertEqual(present, rows)

    def test_genome_without_a_fasta_is_left_out_and_named(self):
        rows = self.rows(('GCA_1.1', self.genome_dir('GCA_1.1_ASM1', fasta=False)))
        present, missing = B.split_by_fasta(rows)
        self.assertEqual(present, [])
        self.assertEqual(missing, ['GCA_1.1'])

    def test_empty_fasta_is_left_out(self):
        rows = self.rows(('GCA_2.1', self.genome_dir('GCA_2.1_ASM2', empty=True)))
        present, missing = B.split_by_fasta(rows)
        self.assertEqual(present, [])
        self.assertEqual(missing, ['GCA_2.1'])

    def test_the_rest_of_the_batch_survives_one_missing_genome(self):
        rows = self.rows(('GCF_1.1', self.genome_dir('GCF_1.1_ASM1')),
                         ('GCA_2.1', self.genome_dir('GCA_2.1_ASM2', fasta=False)))
        present, missing = B.split_by_fasta(rows)
        self.assertEqual(len(present), 1)
        self.assertEqual(missing, ['GCA_2.1'])

    def test_the_order_given_is_the_order_returned_whatever_the_thread_count(self):
        # it is the order the genomes are then worked on in, and the answer must
        # not depend on how the pool happened to be scheduled
        genomes = []
        for i in range(1, 60):
            # every third genome has no FASTA, so both lists are interleaved
            genomes.append(('GCF_%09d.1' % i,
                            self.genome_dir('GCF_%09d.1_ASM%dv1' % (i, i),
                                            fasta=bool(i % 3))))
        rows = self.rows(*genomes)

        serial = B.split_by_fasta(rows, threads=1)
        parallel = B.split_by_fasta(rows, threads=16)

        self.assertEqual(serial, parallel)
        self.assertEqual([accession for _, accession in serial[0]],
                         [accession for accession, _ in genomes
                          if int(accession.split('_')[1].split('.')[0]) % 3])
        self.assertEqual(len(serial[1]), 19)

    def test_more_genomes_than_one_chunk_are_all_answered(self):
        # the pool is fed a chunk at a time, and a genome must not be lost at the seam
        rows = self.rows(*[('GCF_%09d.1' % i,
                            self.genome_dir('GCF_%09d.1_ASM%dv1' % (i, i)))
                           for i in range(1, 12)])

        with mock.patch.object(B, 'STAT_CHUNK', 4):
            present, missing = B.split_by_fasta(rows, threads=3)

        self.assertEqual(missing, [])
        self.assertEqual(present, rows)


class WriteBatchfileTests(TempDirCase):
    """FASTA first, genome ID second; the ID is the accession the genome_dirs file
    names, so a result can be matched back without canonicalising anything."""

    def test_columns_are_fasta_then_accession(self):
        batchfile = os.path.join(self.dir, 'batch.tsv')
        B.write_batchfile([('/rel/x/GCF_1.1_ASM1_genomic.fna.gz', 'GCF_1.1')], batchfile)
        with open(batchfile) as handle:
            self.assertEqual(handle.read(),
                             '/rel/x/GCF_1.1_ASM1_genomic.fna.gz\tGCF_1.1\n')


class BatchfilePathTests(TempDirCase):
    """A directory an earlier version planned is still a planned directory."""

    def batch(self, name, rows=(('/m/a.fna.gz', 'GCA_1.1'),)):
        batch_dir = os.path.join(self.dir, 'batch_000001')
        if not os.path.isdir(batch_dir):
            os.makedirs(batch_dir)
        B.write_batchfile(rows, os.path.join(batch_dir, name),
                          compress=name.endswith(B.GZIP_EXT))
        return batch_dir

    def test_the_gzipped_plan_is_found(self):
        batch_dir = self.batch(BATCHFILE_NAME)
        self.assertEqual(B.batchfile_path(batch_dir, LAYOUT),
                         os.path.join(batch_dir, BATCHFILE_NAME))

    def test_a_plan_an_earlier_version_wrote_is_found(self):
        """r237 was planned before the plan was compressed, and the command is
        run again over a finished output directory to pick up later work."""
        batch_dir = self.batch(LEGACY_BATCHFILE_NAME)
        self.assertEqual(B.batchfile_path(batch_dir, LAYOUT),
                         os.path.join(batch_dir, LEGACY_BATCHFILE_NAME))
        self.assertEqual(B.read_batchfile(B.batchfile_path(batch_dir, LAYOUT)),
                         [('/m/a.fna.gz', 'GCA_1.1')])

    def test_a_batch_planned_by_an_earlier_version_is_not_planned_again(self):
        """Repartitioning a release whose batches are already done would move
        genomes between batches that have finished."""
        out_dir = os.path.join(self.dir, 'out')
        batch_dir = os.path.join(out_dir, 'batch_000001')
        os.makedirs(batch_dir)
        B.write_batchfile([('/m/a.fna.gz', 'GCA_1.1')],
                          os.path.join(batch_dir, LEGACY_BATCHFILE_NAME))
        self.assertEqual(B.batch_dir_names(out_dir, LAYOUT), [batch_dir])

    def test_the_gzipped_plan_wins_where_a_directory_holds_both(self):
        """A directory part-way through an upgrade reads the current one."""
        batch_dir = self.batch(LEGACY_BATCHFILE_NAME,
                               rows=(('/m/old.fna.gz', 'GCA_OLD.1'),))
        self.batch(BATCHFILE_NAME, rows=(('/m/new.fna.gz', 'GCA_NEW.1'),))
        self.assertEqual(B.read_batchfile(B.batchfile_path(batch_dir, LAYOUT)),
                         [('/m/new.fna.gz', 'GCA_NEW.1')])

    def test_a_batch_with_no_plan_names_where_one_would_go(self):
        """So that a caller's error names the file it was looking for."""
        empty = os.path.join(self.dir, 'unplanned')
        os.makedirs(empty)
        self.assertEqual(B.batchfile_path(empty, LAYOUT),
                         os.path.join(empty, BATCHFILE_NAME))


class ClaimTests(TempDirCase):
    """Two machines must never both take one batch, and a reset must not lose one."""

    def batch(self):
        path = os.path.join(self.dir, 'batch_000001')
        os.makedirs(path)
        return path

    def test_an_unclaimed_batch_is_claimed(self):
        batch = self.batch()
        self.assertTrue(B.claim_batch(batch))
        self.assertEqual(B.batch_state(batch), B.STATE_RUNNING)

    def test_a_batch_claimed_by_a_live_process_is_not_taken(self):
        batch = self.batch()
        B.claim_batch(batch)
        self.assertFalse(B.claim_batch(batch))

    def test_a_claim_of_a_dead_process_on_this_host_is_reclaimed(self):
        """What a machine reset leaves behind; nothing else would ever run it."""
        batch = self.batch()
        with open(os.path.join(batch, B.RUNNING_CANARY), 'w') as handle:
            handle.write(B.canary_payload())
        # a PID that cannot be running, recorded against this host
        text = open(os.path.join(batch, B.RUNNING_CANARY)).read()
        with open(os.path.join(batch, B.RUNNING_CANARY), 'w') as handle:
            handle.write(text.replace('pid\t{}'.format(os.getpid()), 'pid\t2147483646'))
        self.assertTrue(B.claim_batch(batch))

    def test_a_claim_of_another_host_is_left_alone(self):
        batch = self.batch()
        with open(os.path.join(batch, B.RUNNING_CANARY), 'w') as handle:
            handle.write('host\tsome-other-machine\npid\t1\ntime\tnow\n')
        self.assertFalse(B.claim_batch(batch))

    def test_reclaim_takes_another_hosts_claim(self):
        batch = self.batch()
        with open(os.path.join(batch, B.RUNNING_CANARY), 'w') as handle:
            handle.write('host\tsome-other-machine\npid\t1\ntime\tnow\n')
        self.assertTrue(B.claim_batch(batch, reclaim=True))

    def test_a_finished_batch_reports_success_and_gives_up_the_claim(self):
        batch = self.batch()
        B.claim_batch(batch)
        B.finish_batch(batch, compared=7)
        self.assertEqual(B.batch_state(batch), B.STATE_SUCCESS)
        self.assertFalse(os.path.exists(os.path.join(batch, B.RUNNING_CANARY)))
        self.assertEqual(B.read_canary(os.path.join(batch, B.SUCCESS_CANARY))['compared'], '7')

    def test_a_failed_batch_gives_up_the_claim_so_it_is_retried(self):
        batch = self.batch()
        B.claim_batch(batch)
        B.fail_batch(batch, 'the work returned exit code 1.')
        self.assertEqual(B.batch_state(batch), B.STATE_FAILED)
        self.assertTrue(B.claim_batch(batch))

    def test_claiming_a_failed_batch_clears_the_failure(self):
        batch = self.batch()
        B.fail_batch(batch, 'whatever')
        B.claim_batch(batch)
        self.assertFalse(os.path.exists(os.path.join(batch, B.FAILED_CANARY)))

    def test_what_a_failed_batch_said_survives_the_retry(self):
        """A batch failing the same way each time is read by what it wrote."""
        batch = self.batch()
        B.fail_batch(batch, 'the work returned exit code 1.')
        B.claim_batch(batch)
        kept = [name for name in os.listdir(batch)
                if name.startswith(B.FAILED_CANARY + '.')]
        self.assertEqual(len(kept), 1)
        self.assertIn('exit code 1',
                      B.read_canary(os.path.join(batch, kept[0]))['reason'])

    def test_success_outranks_running(self):
        batch = self.batch()
        B.claim_batch(batch)
        open(os.path.join(batch, B.SUCCESS_CANARY), 'w').close()
        self.assertEqual(B.batch_state(batch), B.STATE_SUCCESS)


class LeaseTests(TempDirCase):
    """A claim is held by saying so, not by having said so once."""

    def batch(self):
        path = os.path.join(self.dir, 'batch_000001')
        os.makedirs(path)
        return path

    def other_machine_claims(self, batch, age=0.0):
        """A RUNNING file of another host, last touched age seconds ago."""
        running = os.path.join(batch, B.RUNNING_CANARY)
        with open(running, 'w') as handle:
            handle.write('host\tsome-other-machine\npid\t1\ntime\tnow\n')
        touched = B.server_time(batch) - age
        os.utime(running, (touched, touched))
        return running

    def test_a_claim_that_has_stopped_being_touched_is_taken(self):
        """The machine holding it reset, wedged or was killed; nothing else says so."""
        batch = self.batch()
        self.other_machine_claims(batch, age=3 * B.CLAIM_LEASE_SECONDS)
        self.assertTrue(B.claim_batch(batch))

    def test_a_claim_still_being_touched_is_left_alone_however_old_the_run(self):
        """A batch of 10,000 genomes runs for hours and stays its machine's."""
        batch = self.batch()
        self.other_machine_claims(batch, age=1.0)
        self.assertFalse(B.claim_batch(batch))

    def test_the_lease_is_what_says_when_a_claim_has_gone_quiet(self):
        batch = self.batch()
        running = self.other_machine_claims(batch, age=600.0)
        self.assertFalse(B.stale_claim(running, lease=3600))
        self.assertTrue(B.stale_claim(running, lease=60))

    def test_a_heartbeat_keeps_a_claim_from_expiring(self):
        batch = self.batch()
        running = self.other_machine_claims(batch, age=600.0)
        with B.Heartbeat(running, interval=0.05):
            time.sleep(0.3)
            self.assertFalse(B.stale_claim(running, lease=60))

    def test_the_heartbeat_stops_with_the_batch(self):
        """A claim outlives the process holding it by one lease and no longer."""
        batch = self.batch()
        running = self.other_machine_claims(batch)
        with B.Heartbeat(running, interval=0.05) as beat:
            time.sleep(0.1)
        self.assertTrue(beat.stop.is_set())
        self.assertFalse(beat.thread.is_alive())

    def test_the_age_of_a_claim_is_the_time_since_it_was_touched(self):
        batch = self.batch()
        running = self.other_machine_claims(batch, age=1800.0)
        self.assertAlmostEqual(B.claim_age(running), 1800.0, delta=30)

    def test_a_claim_that_has_gone_has_no_age_rather_than_raising(self):
        self.assertIsNone(B.claim_age(os.path.join(self.dir, 'nothing')))

    def test_the_clock_a_lease_is_measured_against_is_the_file_servers(self):
        """Not this machine's: the machines sharing an --out_dir have a clock each."""
        with mock.patch.object(B.time, 'time', return_value=0.0):
            self.assertGreater(B.server_time(self.dir), 1e9)

    def test_an_interrupted_batch_hands_its_claim_straight_back(self):
        batch = self.batch()
        B.claim_batch(batch)
        B.release_claim(batch)
        self.assertEqual(B.batch_state(batch), B.STATE_PENDING)
        self.assertTrue(B.claim_batch(batch))


class BatchLogTests(TempDirCase):
    """Several machines share an --out_dir; no two of them share a log file."""

    def test_what_happens_to_a_batch_is_written_in_the_batch(self):
        batch = os.path.join(self.dir, 'batch_000001')
        os.makedirs(batch)
        logger = logging.getLogger('batching_test_batch_log')
        with B.batch_log(batch, logger, LAYOUT):
            logger.error('batch_000001: failed for a reason')
        self.assertIn('failed for a reason',
                      open(os.path.join(batch, BATCH_LOG_NAME)).read())

    def test_the_log_is_let_go_of_when_the_batch_is(self):
        """A run works through many batches and must not hold a handle on each."""
        batch = os.path.join(self.dir, 'batch_000001')
        os.makedirs(batch)
        logger = logging.getLogger('batching_test_batch_log_handles')
        before = len(logger.handlers)
        with B.batch_log(batch, logger, LAYOUT):
            self.assertEqual(len(logger.handlers), before + 1)
        self.assertEqual(len(logger.handlers), before)

if __name__ == '__main__':
    unittest.main()
