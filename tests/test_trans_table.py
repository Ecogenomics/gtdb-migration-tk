#!/usr/bin/env python3
"""Offline unit tests for trans_table.py -- gTranslate itself is never run.

Run with the interpreter that has tqdm:

    /opt/centos7/sw/miniconda3/envs/python3/bin/python -m unittest -v tests.test_trans_table

What is tested here is everything decided before the subprocess starts: which file
of a genome directory is handed over, which genomes are left out, and the command
line built from the options. Running gTranslate is gTranslate's business; getting
the wrong genomes to it, or the right ones under the wrong names, is this module's.
"""

import gzip
import logging
import os
import shutil
import tempfile
import time
import unittest
from unittest import mock

from gtdb_migration_tk import trans_table as G


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

class GenomicFastaTests(TempDirCase):
    """The genes of a genome end in the same suffix as the genome itself."""

    def test_fasta_is_named_for_the_assembly_directory(self):
        path = G.genomic_fasta('/rel/genbank/GCA/047/639/395/GCA_047639395.1_ASM4763939v1')
        self.assertEqual(os.path.basename(path),
                         'GCA_047639395.1_ASM4763939v1_genomic.fna.gz')

    def test_trailing_separator_does_not_change_the_name(self):
        without = G.genomic_fasta('/rel/refseq/GCF/000/001/405/GCF_000001405.39_GRCh38.p13')
        with_sep = G.genomic_fasta('/rel/refseq/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/')
        self.assertEqual(without, with_sep)


# ------------------------------------------------------------- reading the release

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
        self.assertEqual(G.read_genome_dirs(path),
                         [('GCF_000001405.39', '/rel/refseq/GCF/000/001/405/x')])

    def test_further_columns_are_ignored(self):
        path = self.write_genome_dirs('GCA_1.1\t/rel/x\tG1\tsomething\telse\n')
        self.assertEqual(G.read_genome_dirs(path), [('GCA_1.1', '/rel/x')])

    def test_blank_lines_are_skipped(self):
        path = self.write_genome_dirs('GCA_1.1\t/rel/x\tG1\n\n')
        self.assertEqual(len(G.read_genome_dirs(path)), 1)


# ------------------------------------------------------------- what is asked about

class SplitByFastaTests(TempDirCase):
    """gTranslate refuses a whole batch over one missing path, so one is left out here."""

    def rows(self, *genomes):
        """(FASTA path, accession) as a batchfile names them."""
        return [(G.genomic_fasta(path), accession) for accession, path in genomes]

    def test_genome_with_a_fasta_is_kept(self):
        rows = self.rows(('GCF_000001405.39',
                          self.genome_dir('GCF_000001405.39_GRCh38.p13')))
        present, missing = G.split_by_fasta(rows)
        self.assertEqual(missing, [])
        self.assertEqual(present, rows)

    def test_genome_without_a_fasta_is_left_out_and_named(self):
        rows = self.rows(('GCA_1.1', self.genome_dir('GCA_1.1_ASM1', fasta=False)))
        present, missing = G.split_by_fasta(rows)
        self.assertEqual(present, [])
        self.assertEqual(missing, ['GCA_1.1'])

    def test_empty_fasta_is_left_out(self):
        rows = self.rows(('GCA_2.1', self.genome_dir('GCA_2.1_ASM2', empty=True)))
        present, missing = G.split_by_fasta(rows)
        self.assertEqual(present, [])
        self.assertEqual(missing, ['GCA_2.1'])

    def test_the_rest_of_the_batch_survives_one_missing_genome(self):
        rows = self.rows(('GCF_1.1', self.genome_dir('GCF_1.1_ASM1')),
                         ('GCA_2.1', self.genome_dir('GCA_2.1_ASM2', fasta=False)))
        present, missing = G.split_by_fasta(rows)
        self.assertEqual(len(present), 1)
        self.assertEqual(missing, ['GCA_2.1'])

    def test_the_order_given_is_the_order_returned_whatever_the_thread_count(self):
        # it is the order gTranslate is handed the genomes in, and the answer must
        # not depend on how the pool happened to be scheduled
        genomes = []
        for i in range(1, 60):
            # every third genome has no FASTA, so both lists are interleaved
            genomes.append(('GCF_%09d.1' % i,
                            self.genome_dir('GCF_%09d.1_ASM%dv1' % (i, i),
                                            fasta=bool(i % 3))))
        rows = self.rows(*genomes)

        serial = G.split_by_fasta(rows, threads=1)
        parallel = G.split_by_fasta(rows, threads=16)

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

        with mock.patch.object(G, 'STAT_CHUNK', 4):
            present, missing = G.split_by_fasta(rows, threads=3)

        self.assertEqual(missing, [])
        self.assertEqual(present, rows)


class CheckBatchFastasTests(TempDirCase):
    """The check belongs to the batch that is about to run, not to the plan."""

    def batch(self, *genomes):
        """A batch directory holding the batchfile the plan would have cut."""
        batch_dir = os.path.join(self.dir, 'batch_000001')
        os.makedirs(batch_dir)
        G.write_batchfile([(G.genomic_fasta(path), accession)
                           for accession, path in genomes],
                          os.path.join(batch_dir, G.BATCHFILE_NAME))
        return batch_dir

    def test_a_whole_batch_is_handed_over_as_it_stands(self):
        # nothing is written where nothing is wrong, and gTranslate reads the
        # batchfile the plan cut
        batch_dir = self.batch(('GCF_1.1', self.genome_dir('GCF_1.1_ASM1')),
                               ('GCF_2.1', self.genome_dir('GCF_2.1_ASM2')))

        batchfile, present, missing = G.check_batch_fastas(batch_dir)

        self.assertEqual(batchfile, os.path.join(batch_dir, G.BATCHFILE_NAME))
        self.assertEqual(len(present), 2)
        self.assertEqual(missing, [])
        self.assertEqual(os.listdir(batch_dir), [G.BATCHFILE_NAME])

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
            os.path.join(batch_dir, G.BATCHFILE_NAME))], ['GCF_1.1', 'GCA_2.1'])


class FastaSizeTests(TempDirCase):
    """One stat, where two calls asked the file server the same question twice."""

    def test_the_size_of_a_fasta_is_returned(self):
        path = G.genomic_fasta(self.genome_dir('GCF_1.1_ASM1'))
        self.assertEqual(G.fasta_size(path), os.path.getsize(path))

    def test_a_file_that_is_not_there_has_no_size_rather_than_raising(self):
        self.assertEqual(G.fasta_size(os.path.join(self.dir, 'no_such.fna.gz')), 0)

    def test_an_empty_fasta_has_no_size(self):
        path = G.genomic_fasta(self.genome_dir('GCA_2.1_ASM2', empty=True))
        self.assertEqual(G.fasta_size(path), 0)


class WriteBatchfileTests(TempDirCase):
    """gTranslate reads FASTA first, genome ID second; the ID is the accession."""

    def test_columns_are_fasta_then_accession(self):
        batchfile = os.path.join(self.dir, 'batch.tsv')
        G.write_batchfile([('/rel/x/GCF_1.1_ASM1_genomic.fna.gz', 'GCF_1.1')], batchfile)
        with open(batchfile) as handle:
            self.assertEqual(handle.read(),
                             '/rel/x/GCF_1.1_ASM1_genomic.fna.gz\tGCF_1.1\n')


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
                                   2, self.dir)
        self.assertEqual(len(batches), 3)
        self.assertEqual(len(G.read_batchfile(
            os.path.join(batches[0], G.BATCHFILE_NAME))), 2)
        self.assertEqual(len(G.read_batchfile(
            os.path.join(batches[-1], G.BATCHFILE_NAME))), 1)

    def test_batches_are_numbered_in_order(self):
        batches = G.create_batches(self.rows('GCF_1.1', 'GCF_2.1'), 1, self.dir)
        self.assertEqual([os.path.basename(b) for b in batches],
                         ['batch_000001', 'batch_000002'])

    def test_an_existing_plan_is_found_and_reused(self):
        G.create_batches(self.rows('GCF_1.1', 'GCF_2.1'), 1, self.dir)
        self.assertEqual(len(G.batch_dir_names(self.dir)), 2)

    def test_a_directory_without_a_batchfile_is_not_a_batch(self):
        os.makedirs(os.path.join(self.dir, 'batch_000001'))
        self.assertEqual(G.batch_dir_names(self.dir), [])

    def test_batchfile_survives_a_round_trip(self):
        batches = G.create_batches(self.rows('GCF_1.1'), 10, self.dir)
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

class ClaimTests(TempDirCase):
    """Two machines must never both take one batch, and a reset must not lose one."""

    def batch(self):
        path = os.path.join(self.dir, 'batch_000001')
        os.makedirs(path)
        return path

    def test_an_unclaimed_batch_is_claimed(self):
        batch = self.batch()
        self.assertTrue(G.claim_batch(batch))
        self.assertEqual(G.batch_state(batch), G.STATE_RUNNING)

    def test_a_batch_claimed_by_a_live_process_is_not_taken(self):
        batch = self.batch()
        G.claim_batch(batch)
        self.assertFalse(G.claim_batch(batch))

    def test_a_claim_of_a_dead_process_on_this_host_is_reclaimed(self):
        """What a machine reset leaves behind; nothing else would ever run it."""
        batch = self.batch()
        with open(os.path.join(batch, G.RUNNING_CANARY), 'w') as handle:
            handle.write(G.canary_payload())
        # a PID that cannot be running, recorded against this host
        text = open(os.path.join(batch, G.RUNNING_CANARY)).read()
        with open(os.path.join(batch, G.RUNNING_CANARY), 'w') as handle:
            handle.write(text.replace('pid\t{}'.format(os.getpid()), 'pid\t2147483646'))
        self.assertTrue(G.claim_batch(batch))

    def test_a_claim_of_another_host_is_left_alone(self):
        batch = self.batch()
        with open(os.path.join(batch, G.RUNNING_CANARY), 'w') as handle:
            handle.write('host\tsome-other-machine\npid\t1\ntime\tnow\n')
        self.assertFalse(G.claim_batch(batch))

    def test_reclaim_takes_another_hosts_claim(self):
        batch = self.batch()
        with open(os.path.join(batch, G.RUNNING_CANARY), 'w') as handle:
            handle.write('host\tsome-other-machine\npid\t1\ntime\tnow\n')
        self.assertTrue(G.claim_batch(batch, reclaim=True))

    def test_a_finished_batch_reports_success_and_gives_up_the_claim(self):
        batch = self.batch()
        G.claim_batch(batch)
        G.finish_batch(batch, compared=7)
        self.assertEqual(G.batch_state(batch), G.STATE_SUCCESS)
        self.assertFalse(os.path.exists(os.path.join(batch, G.RUNNING_CANARY)))
        self.assertEqual(G.read_canary(os.path.join(batch, G.SUCCESS_CANARY))['compared'], '7')

    def test_a_failed_batch_gives_up_the_claim_so_it_is_retried(self):
        batch = self.batch()
        G.claim_batch(batch)
        G.fail_batch(batch, 'gtranslate returned exit code 1.')
        self.assertEqual(G.batch_state(batch), G.STATE_FAILED)
        self.assertTrue(G.claim_batch(batch))

    def test_claiming_a_failed_batch_clears_the_failure(self):
        batch = self.batch()
        G.fail_batch(batch, 'whatever')
        G.claim_batch(batch)
        self.assertFalse(os.path.exists(os.path.join(batch, G.FAILED_CANARY)))

    def test_what_a_failed_batch_said_survives_the_retry(self):
        """A batch failing the same way each time is read by what it wrote."""
        batch = self.batch()
        G.fail_batch(batch, 'gtranslate returned exit code 1.')
        G.claim_batch(batch)
        kept = [name for name in os.listdir(batch)
                if name.startswith(G.FAILED_CANARY + '.')]
        self.assertEqual(len(kept), 1)
        self.assertIn('exit code 1',
                      G.read_canary(os.path.join(batch, kept[0]))['reason'])

    def test_success_outranks_running(self):
        batch = self.batch()
        G.claim_batch(batch)
        open(os.path.join(batch, G.SUCCESS_CANARY), 'w').close()
        self.assertEqual(G.batch_state(batch), G.STATE_SUCCESS)


# ------------------------------------------------------------- the lease

class LeaseTests(TempDirCase):
    """A claim is held by saying so, not by having said so once."""

    def batch(self):
        path = os.path.join(self.dir, 'batch_000001')
        os.makedirs(path)
        return path

    def other_machine_claims(self, batch, age=0.0):
        """A RUNNING file of another host, last touched age seconds ago."""
        running = os.path.join(batch, G.RUNNING_CANARY)
        with open(running, 'w') as handle:
            handle.write('host\tsome-other-machine\npid\t1\ntime\tnow\n')
        touched = G.server_time(batch) - age
        os.utime(running, (touched, touched))
        return running

    def test_a_claim_that_has_stopped_being_touched_is_taken(self):
        """The machine holding it reset, wedged or was killed; nothing else says so."""
        batch = self.batch()
        self.other_machine_claims(batch, age=3 * G.CLAIM_LEASE_SECONDS)
        self.assertTrue(G.claim_batch(batch))

    def test_a_claim_still_being_touched_is_left_alone_however_old_the_run(self):
        """A batch of 10,000 genomes runs for hours and stays its machine's."""
        batch = self.batch()
        self.other_machine_claims(batch, age=1.0)
        self.assertFalse(G.claim_batch(batch))

    def test_the_lease_is_what_says_when_a_claim_has_gone_quiet(self):
        batch = self.batch()
        running = self.other_machine_claims(batch, age=600.0)
        self.assertFalse(G.stale_claim(running, lease=3600))
        self.assertTrue(G.stale_claim(running, lease=60))

    def test_a_heartbeat_keeps_a_claim_from_expiring(self):
        batch = self.batch()
        running = self.other_machine_claims(batch, age=600.0)
        with G.Heartbeat(running, interval=0.05):
            time.sleep(0.3)
            self.assertFalse(G.stale_claim(running, lease=60))

    def test_the_heartbeat_stops_with_the_batch(self):
        """A claim outlives the process holding it by one lease and no longer."""
        batch = self.batch()
        running = self.other_machine_claims(batch)
        with G.Heartbeat(running, interval=0.05) as beat:
            time.sleep(0.1)
        self.assertTrue(beat.stop.is_set())
        self.assertFalse(beat.thread.is_alive())

    def test_the_age_of_a_claim_is_the_time_since_it_was_touched(self):
        batch = self.batch()
        running = self.other_machine_claims(batch, age=1800.0)
        self.assertAlmostEqual(G.claim_age(running), 1800.0, delta=30)

    def test_a_claim_that_has_gone_has_no_age_rather_than_raising(self):
        self.assertIsNone(G.claim_age(os.path.join(self.dir, 'nothing')))

    def test_the_clock_a_lease_is_measured_against_is_the_file_servers(self):
        """Not this machine's: the machines sharing an --out_dir have a clock each."""
        with mock.patch.object(G.time, 'time', return_value=0.0):
            self.assertGreater(G.server_time(self.dir), 1e9)

    def test_an_interrupted_batch_hands_its_claim_straight_back(self):
        batch = self.batch()
        G.claim_batch(batch)
        G.release_claim(batch)
        self.assertEqual(G.batch_state(batch), G.STATE_PENDING)
        self.assertTrue(G.claim_batch(batch))


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

class BatchLogTests(TempDirCase):
    """Several machines share an --out_dir; no two of them share a log file."""

    def test_what_happens_to_a_batch_is_written_in_the_batch(self):
        batch = os.path.join(self.dir, 'batch_000001')
        os.makedirs(batch)
        logger = logging.getLogger('trans_table_test_batch_log')
        with G.batch_log(batch, logger):
            logger.error('batch_000001: failed for a reason')
        self.assertIn('failed for a reason',
                      open(os.path.join(batch, G.BATCH_LOG_NAME)).read())

    def test_the_log_is_let_go_of_when_the_batch_is(self):
        """A run works through many batches and must not hold a handle on each."""
        batch = os.path.join(self.dir, 'batch_000001')
        os.makedirs(batch)
        logger = logging.getLogger('trans_table_test_batch_log_handles')
        before = len(logger.handlers)
        with G.batch_log(batch, logger):
            self.assertEqual(len(logger.handlers), before + 1)
        self.assertEqual(len(logger.handlers), before)


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
        G.write_conflicts([], os.path.join(batch_dir, G.CONFLICT_NAME))
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
    """A conflict is a genome GTDB would call genes for under a table NCBI rejects,
    and a conflict is the only thing the file holds."""

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

    def test_a_genome_agreeing_with_ncbi_is_counted_and_not_reported(self):
        """Agreement is nearly every genome of a release; the file is the
        disagreements, and the count is what says how many were looked at."""
        path = self.genome_dir('GCF_1.1', ncbi_table=11)
        rows, compared, _ = G.conflict_rows(
            G.read_translation_table_summary(self.summary(('GCF_1.1', '11'))),
            {'GCF_1.1': path}, {})
        self.assertEqual(rows, [])
        self.assertEqual(compared, 1)

    def test_differing_tables_conflict(self):
        path = self.genome_dir('GCF_1.1', ncbi_table=11)
        rows, compared, _ = G.conflict_rows(
            G.read_translation_table_summary(self.summary(('GCF_1.1', '4'))),
            {'GCF_1.1': path}, {})
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0][G.CONFLICT_HEADER.index('gtranslate_tt')], '4')
        self.assertEqual(rows[0][G.CONFLICT_HEADER.index('ncbi_tt')], '11')
        self.assertEqual(compared, 1)

    def test_the_file_says_nothing_about_whether_a_row_is_a_conflict(self):
        """Every row of it is one, so a result column would say 'conflict' and
        nothing else on every line of the release."""
        self.assertNotIn('result', G.CONFLICT_HEADER)
        path = os.path.join(self.dir, 'conflicts.tsv')
        G.write_conflicts([('GCF_1.1', '4', '11', '4', 'True', '90.1', '64.2', 'na')], path)
        with open(path) as handle:
            header, row = handle.read().splitlines()
        self.assertEqual(header.split('\t'), list(G.CONFLICT_HEADER))
        self.assertEqual(len(row.split('\t')), len(G.CONFLICT_HEADER))

    def test_a_batch_with_nothing_to_report_still_writes_the_file(self):
        """The release file is every batch's concatenated, so a batch that
        conflicted nowhere has to leave a header behind."""
        path = os.path.join(self.dir, 'none.tsv')
        G.write_conflicts([], path)
        self.assertEqual(open(path).read().splitlines(), ['\t'.join(G.CONFLICT_HEADER)])

    def test_genome_ncbi_declares_no_table_for_is_left_out(self):
        path = self.genome_dir('GCF_1.1')
        rows, compared, no_table = G.conflict_rows(
            G.read_translation_table_summary(self.summary(('GCF_1.1', '11'))), {'GCF_1.1': path}, {})
        self.assertEqual(rows, [])
        self.assertEqual(compared, 0)
        self.assertEqual(no_table, 1)

    def test_coding_densities_and_lineage_are_carried(self):
        path = self.genome_dir('GCF_1.1', ncbi_table=11)
        rows, _, _ = G.conflict_rows(
            G.read_translation_table_summary(self.summary(('GCF_1.1', '4'))),
            {'GCF_1.1': path},
            {'GCF_1.1': 'd__Bacteria;p__Pseudomonadota;s__Escherichia coli'})
        self.assertEqual(rows[0][G.CONFLICT_HEADER.index('coding_density_4')], '90.1')
        self.assertEqual(rows[0][G.CONFLICT_HEADER.index('coding_density_11')], '64.2')
        self.assertIn('d__Bacteria', rows[0][G.CONFLICT_HEADER.index('ncbi_taxonomy')])

    def test_lineage_is_found_through_the_canonical_accession(self):
        """A GenBank genome takes the lineage held against its RefSeq counterpart."""
        path = self.genome_dir('GCA_005435135.1', ncbi_table=11)
        rows, _, _ = G.conflict_rows(
            G.read_translation_table_summary(self.summary(('GCA_005435135.1', '4'))),
            {'GCA_005435135.1': path},
            G.read_taxonomy(self.write_taxonomy('GCF_005435135.1\td__Bacteria;s__X\n')))
        self.assertEqual(rows[0][G.CONFLICT_HEADER.index('ncbi_taxonomy')], 'd__Bacteria;s__X')

    def test_genome_missing_from_the_taxonomy_is_still_reported(self):
        path = self.genome_dir('GCF_1.1', ncbi_table=4)
        rows, _, _ = G.conflict_rows(
            G.read_translation_table_summary(self.summary(('GCF_1.1', '11'))),
            {'GCF_1.1': path}, {})
        self.assertEqual(rows[0][G.CONFLICT_HEADER.index('ncbi_taxonomy')], 'na')

    def write_taxonomy(self, text):
        path = os.path.join(self.dir, 'taxonomy.tsv')
        with open(path, 'w') as handle:
            handle.write(text)
        return path

    def test_checkm_table_is_reported_beside_the_others(self):
        """The density rule alone cannot express 25, which is the point of it."""
        path = self.genome_dir('GCF_1.1', ncbi_table=11)
        rows, _, _ = G.conflict_rows(
            G.read_translation_table_summary(self.summary(('GCF_1.1', '25'))),
            {'GCF_1.1': path}, {})
        row = rows[0]
        self.assertEqual(row[G.CONFLICT_HEADER.index('gtranslate_tt')], '25')
        self.assertEqual(row[G.CONFLICT_HEADER.index('checkm_tt')], '4')

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
        self.assertEqual(G.CONFLICT_HEADER.index('checkm_conflict'),
                         G.CONFLICT_HEADER.index('checkm_tt') + 1)

        path = self.genome_dir('GCF_1.1', ncbi_table=4)
        rows, _, _ = G.conflict_rows(
            G.read_translation_table_summary(self.summary(('GCF_1.1', '11'))),
            {'GCF_1.1': path}, {})
        self.assertEqual(rows[0][G.CONFLICT_HEADER.index('checkm_tt')], '4')
        self.assertEqual(rows[0][G.CONFLICT_HEADER.index('checkm_conflict')], 'True')

    def test_checkm_conflict_is_false_on_a_row_the_density_rule_agrees_with(self):
        path = self.genome_dir('GCF_1.1', ncbi_table=11)
        rows, _, _ = G.conflict_rows(
            G.read_translation_table_summary(self.summary(('GCF_1.1', '25'))),
            {'GCF_1.1': path}, {})
        self.assertEqual(rows[0][G.CONFLICT_HEADER.index('checkm_tt')], '4')
        self.assertEqual(rows[0][G.CONFLICT_HEADER.index('checkm_conflict')], 'False')

    def test_checkm_table_is_11_where_the_densities_are_close(self):
        path = self.genome_dir('GCF_1.1', ncbi_table=4)
        summary = os.path.join(self.dir, 'close.tsv')
        with open(summary, 'w') as handle:
            handle.write('user_genome\tbest_tln_table\tcoding_density_4\t'
                         'coding_density_11\n')
            handle.write('GCF_1.1\t11\t86.24689\t86.64953\n')
        rows, _, _ = G.conflict_rows(
            G.read_translation_table_summary(summary), {'GCF_1.1': path}, {})
        self.assertEqual(rows[0][G.CONFLICT_HEADER.index('checkm_tt')], '11')


class AggregationTests(TempDirCase):
    """The release file is the whole release or absent, never a part of it."""

    def comparison(self, name, *rows):
        path = os.path.join(self.dir, name)
        G.write_conflicts(rows, path)
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
        return (accession, gtranslate_tt, ncbi_tt, '4', 'False', '90.1', '64.2', 'na')

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
        return [accession, gtranslate_tt, ncbi_tt, '4', 'False', '90.1', '64.2',
                'd__Bacteria']

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
        G.write_conflicts([], path)
        self.assertEqual(G.read_conflicts(path), [])

    def test_a_row_is_read_back_as_it_was_written(self):
        path = os.path.join(self.dir, G.CONFLICT_NAME)
        row = ('GCA_1.1', '25', '11', '4', 'False', '90.1', '64.2', 'd__Bacteria')
        G.write_conflicts([row], path)
        self.assertEqual(G.read_conflicts(path), [list(row)])

    def test_a_file_already_annotated_is_read_as_the_row_it_was_made_from(self):
        """The command is run again over a finished output directory, so the
        release file it reads is one it has already annotated; the columns are
        taken by name so the CheckM2 ones are simply not among them."""
        path = os.path.join(self.dir, G.CONFLICT_NAME)
        row = ('GCA_1.1', '25', '11', '4', 'False', '90.1', '64.2', 'd__Bacteria')
        # built by the annotator rather than by hand, so that a column added to
        # CONFLICT_HEADER_CHECKM2 cannot quietly turn this into a row too short
        # to be read back at all
        annotated = G.annotate_conflicts([list(row)],
                                         {25: {'GCA_1.1': ('94.3', '0.17')},
                                          11: {'GCA_1.1': ('51.0', '16.4')}})
        G.write_conflicts(annotated, path, header=G.CONFLICT_HEADER_CHECKM2)
        self.assertEqual(len(annotated[0]), len(G.CONFLICT_HEADER_CHECKM2))
        self.assertEqual(G.read_conflicts(path), [list(row)])

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
        G.write_conflicts(rows, os.path.join(self.out_dir, G.CONFLICT_NAME))

    def release(self):
        with open(os.path.join(self.out_dir, G.CONFLICT_NAME)) as handle:
            return [line.rstrip('\n').split('\t') for line in handle]

    def row(self, accession='GCA_1.1', gtranslate_tt='25', ncbi_tt='11'):
        return (accession, gtranslate_tt, ncbi_tt, '4', 'False', '90.1', '64.2',
                'd__Bacteria')

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
        return ['GCA_1.1', '25', '11', '4', 'False', '90.1', '64.2', 'd__Bacteria']

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
