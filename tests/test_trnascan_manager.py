#!/usr/bin/env python3
"""Offline unit tests for trnascan_manager.py -- tRNAscan-SE itself is never run.

Run with the interpreter that has tqdm:

    /opt/centos7/sw/miniconda3/envs/gtdb_migration_tk-r237/bin/python -m unittest discover -s tests -p test_trnascan_manager.py

subprocess.Popen is replaced by a stub that writes the files the real
tRNAscan-SE writes and records the model flag it was given, so what is tested
here is the bookkeeping around the scan: which genomes are handed over, which
model each is scanned with, which are reported as having got nothing, and what a
batch records about itself. The stub is a module-level class because the workers
run in forked processes.
"""

import gzip
import logging
import os
import shutil
import tempfile
import unittest
from unittest import mock

from gtdb_migration_tk import batching as B
from gtdb_migration_tk import trnascan_manager as T


# What the stub was asked to do, written into the log file tRNAscan-SE would
# write, since a forked worker cannot hand an object back to the test.
FLAG_LINE = 'model flag: {}\n'

# A genome whose sequences say this is one the stub refuses, standing in for the
# genomes the real tRNAscan-SE cannot process.
REFUSE = 'REFUSE'


class StubPopen(object):
    """Stands in for the tRNAscan-SE process."""

    def __init__(self, command, stdout=None, stderr=None):
        self.command = command
        flag = command[1]
        args = {command[i]: command[i + 1]
                for i in range(len(command)) if command[i] in ('-o', '-m', '-l')}
        fasta = command[-1]

        with open(fasta) as handle:
            refused = REFUSE in handle.read()

        self.returncode = 3 if refused else 0
        if refused:
            return

        with open(args['-o'], 'w') as handle:
            handle.write('Name\ttRNA#\tBegin\tEnd\tType\n')
            handle.write('{}\t1\t10\t80\tAla\n'.format(os.path.basename(fasta)))
        with open(args['-m'], 'w') as handle:
            handle.write('tRNAs decoding Standard 20 AA:\t1\n')
        with open(args['-l'], 'w') as handle:
            handle.write(FLAG_LINE.format(flag))

    def communicate(self):
        return (b'', b'refused' if self.returncode else b'')


class QuietTqdm(object):
    """tqdm with the drawing taken out: iterated over, and used as a context
    manager with a bar to update, as batching.py and this module use it."""

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


class TempDirCase(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='trnascan_test.')
        self.out_dir = os.path.join(self.dir, 'out')

        # the progress bars and the per-genome warnings do not belong in the
        # output of a test run; assertLogs() raises the logger back up
        for module in (T, B):
            patch = mock.patch.object(module, 'tqdm', QuietTqdm)
            patch.start()
            self.addCleanup(patch.stop)

        logger = logging.getLogger('timestamp')
        self.addCleanup(logger.setLevel, logger.level)
        logger.setLevel(logging.CRITICAL)

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    # ------------------------------------------------------------- the release

    def summary(self, name, accessions, compress=True):
        """An NCBI assembly summary file, gzipped as a release holds it."""
        path = os.path.join(self.dir, name)
        opener = gzip.open if compress else open
        with opener(path, 'wt') as handle:
            handle.write('#   See ftp://ftp.ncbi.nlm.nih.gov/genomes/README.txt\n')
            handle.write('#assembly_accession\tbioproject\tftp_path\n')
            for accession in accessions:
                handle.write('{}\tPRJNA1\thttps://x\n'.format(accession))
        return path

    def summaries(self, archaea=(), bacteria=(), compress=True):
        """The four files the command takes, archaea and bacteria from GenBank.

        Named .gz only when they are gzipped: open_summary() decides by the
        name, which is how NCBI's own uncompressed files are read.
        """
        ext = '.txt.gz' if compress else '.txt'
        return dict(
            gbk_arc_assembly_file=self.summary('arc_gbk' + ext, archaea, compress),
            gbk_bac_assembly_file=self.summary('bac_gbk' + ext, bacteria, compress),
            rfq_arc_assembly_file=self.summary('arc_rfq' + ext, (), compress),
            rfq_bac_assembly_file=self.summary('bac_rfq' + ext, (), compress))

    def genome_dir(self, accession, fasta=True, refuse=False, trna=None):
        """A genome directory as a release holds it.

        trna is None for a genome with no tRNAs, 'valid' for one whose results
        are vouched for, 'stale' for one whose checksum disagrees, and 'nosum'
        for one whose checksum file was never written.
        """
        assembly = '{}_ASM{}v1'.format(accession, accession[4:10])
        gpath = os.path.join(self.dir, assembly)
        os.makedirs(gpath, exist_ok=True)

        if fasta:
            path = os.path.join(gpath, assembly + '_genomic.fna.gz')
            with gzip.open(path, 'wt') as handle:
                handle.write('>contig_1\n{}\n'.format(REFUSE if refuse else 'ACGT' * 20))

        if trna:
            trna_dir = os.path.join(gpath, T.TRNA_DIR)
            os.makedirs(trna_dir, exist_ok=True)
            table = os.path.join(trna_dir, accession + T.TRNA_EXT)
            with open(table, 'w') as handle:
                handle.write('Name\ttRNA#\n')
            if trna != 'nosum':
                from gtdb_migration_tk.biolib_lite.checksum import sha256
                with open(table + T.CHECKSUM_EXT, 'w') as handle:
                    handle.write('{}\n'.format(
                        sha256(table) if trna == 'valid' else 'not-the-checksum'))

        return gpath

    def genome_dirs_file(self, genomes):
        path = os.path.join(self.dir, 'genome_dirs.tsv')
        with open(path, 'w') as handle:
            for accession, gpath in genomes:
                handle.write('{}\t{}\t{}\n'.format(accession, gpath, 'G' + accession[4:-2]))
        return path

    # --------------------------------------------------------------- the run

    def scanner(self, archaea=(), bacteria=(), compress=True, **kwargs):
        with mock.patch.object(T, 'check_dependencies'):
            return T.tRNAScan(cpus=1, tmp_dir=os.path.join(self.dir, 'tmp'),
                              **dict(self.summaries(archaea, bacteria, compress),
                                     **kwargs))

    def run_trnascan(self, genomes, archaea=(), bacteria=(), all_genomes=False,
                     batch_size=B.DEFAULT_BATCH_SIZE):
        scanner = self.scanner(archaea, bacteria, batch_size=batch_size)
        with mock.patch.object(T.subprocess, 'Popen', StubPopen):
            ok = scanner.run(self.genome_dirs_file(genomes), self.out_dir, all_genomes)
        return scanner, ok

    # ------------------------------------------------------------- reading back

    def not_scanned(self, path):
        """A not_scanned table, as (genome_id, reason) pairs."""
        with open(path) as handle:
            header = handle.readline().rstrip('\n').split('\t')
            self.assertEqual(tuple(header), T.NOT_SCANNED_HEADER)
            return sorted(tuple(line.rstrip('\n').split('\t')) for line in handle)

    def release_report(self):
        return self.not_scanned(os.path.join(self.out_dir, T.NOT_SCANNED_RELEASE_NAME))

    def model_flag(self, gpath, accession):
        """The flag tRNAscan-SE was given for this genome, or None if never run."""
        log = os.path.join(gpath, T.TRNA_DIR, accession + T.TRNA_LOG_EXT)
        if not os.path.exists(log):
            return None
        with open(log) as handle:
            return handle.read().strip().split()[-1]

    def batches(self):
        return B.batch_dir_names(self.out_dir, T.LAYOUT)


# ------------------------------------------------- reading the assembly summaries

class TheDomainOfEachGenome(TempDirCase):
    """tRNAscan-SE searches with a bacterial or an archaeal model, so this decides
    the answer, and it comes out of files the module could not previously open."""

    def test_the_gzipped_summaries_a_release_holds_are_read(self):
        scanner = self.scanner(archaea=['GCA_000002.1'], bacteria=['GCA_000001.1'])

        self.assertEqual(scanner.domains, {'GCA_000002.1': T.DOMAIN_ARCHAEA,
                                           'GCA_000001.1': T.DOMAIN_BACTERIA})

    def test_an_uncompressed_summary_is_read_too(self):
        scanner = self.scanner(archaea=['GCA_000002.1'], compress=False)

        self.assertEqual(scanner.domains, {'GCA_000002.1': T.DOMAIN_ARCHAEA})

    def test_the_header_row_is_not_taken_for_a_genome(self):
        """NCBI writes two comment lines above the data, not one."""
        scanner = self.scanner(bacteria=['GCA_000001.1'])

        self.assertEqual(list(scanner.domains), ['GCA_000001.1'])

    def test_an_archaeon_is_scanned_with_the_archaeal_model(self):
        scanner = self.scanner(archaea=['GCA_000002.1'])

        self.assertEqual(scanner.domain_flag('GCA_000002.1'), T.ARCHAEAL_FLAG)

    def test_a_genome_in_none_of_the_files_is_scanned_as_a_bacterium(self):
        scanner = self.scanner(bacteria=['GCA_000001.1'])

        self.assertEqual(scanner.domain_flag('GCA_000009.1'), T.BACTERIAL_FLAG)

    def test_the_model_each_genome_was_scanned_with(self):
        arc = self.genome_dir('GCA_000002.1')
        bac = self.genome_dir('GCA_000001.1')
        unknown = self.genome_dir('GCA_000009.1')

        self.run_trnascan([('GCA_000002.1', arc), ('GCA_000001.1', bac),
                           ('GCA_000009.1', unknown)],
                          archaea=['GCA_000002.1'], bacteria=['GCA_000001.1'])

        self.assertEqual(self.model_flag(arc, 'GCA_000002.1'), T.ARCHAEAL_FLAG)
        self.assertEqual(self.model_flag(bac, 'GCA_000001.1'), T.BACTERIAL_FLAG)
        self.assertEqual(self.model_flag(unknown, 'GCA_000009.1'), T.BACTERIAL_FLAG)

    def test_a_genome_of_no_known_domain_is_warned_about(self):
        unknown = self.genome_dir('GCA_000009.1')

        with self.assertLogs('timestamp', level='WARNING') as captured:
            self.run_trnascan([('GCA_000009.1', unknown)])

        self.assertTrue(any('none of the NCBI assembly summary files' in r.getMessage()
                            for r in captured.records), captured.records)


# ---------------------------------------------------------- what decides the work

class WhatDecidesTheWork(TempDirCase):
    """The checksum beside a genome's own results, and not a report."""

    def decide(self, accession, **kwargs):
        gpath = self.genome_dir(accession, **kwargs)
        scanner = self.scanner()
        job = T.TrnaJob(accession,
                        os.path.join(gpath, os.path.basename(gpath) + '_genomic.fna.gz'),
                        T.BACTERIAL_FLAG)
        return scanner.trnascan_parser(job)

    def test_a_genome_with_no_trnas_is_scanned(self):
        self.assertIsNotNone(self.decide('GCA_000001.1'))

    def test_a_genome_whose_trnas_are_vouched_for_is_skipped(self):
        self.assertIsNone(self.decide('GCA_000001.1', trna='valid'))

    def test_a_genome_whose_checksum_disagrees_is_scanned_again(self):
        """A table written by a run interrupted partway through it."""
        self.assertIsNotNone(self.decide('GCA_000001.1', trna='stale'))

    def test_a_table_with_no_checksum_beside_it_is_scanned_again(self):
        self.assertIsNotNone(self.decide('GCA_000001.1', trna='nosum'))

    def test_a_genome_already_scanned_is_not_handed_to_trnascan(self):
        gpath = self.genome_dir('GCA_000001.1', trna='valid')

        self.run_trnascan([('GCA_000001.1', gpath)])

        self.assertIsNone(self.model_flag(gpath, 'GCA_000001.1'))

    def test_all_scans_it_again(self):
        gpath = self.genome_dir('GCA_000001.1', trna='valid')

        self.run_trnascan([('GCA_000001.1', gpath)], all_genomes=True)

        self.assertIsNotNone(self.model_flag(gpath, 'GCA_000001.1'))


# ------------------------------------------------------ a genome missing its FASTA

class AGenomeWithNoSequences(TempDirCase):
    """Named and left, rather than taking the other ten thousand of the batch."""

    def test_it_is_named_in_the_release_report(self):
        whole = self.genome_dir('GCA_000001.1')
        gone = self.genome_dir('GCA_000003.1', fasta=False)

        self.run_trnascan([('GCA_000001.1', whole), ('GCA_000003.1', gone)])

        self.assertEqual(self.release_report(),
                         [('GCA_000003.1', T.REASON_NO_GENOMIC_FASTA)])

    def test_the_other_genomes_of_the_batch_are_still_scanned(self):
        whole = self.genome_dir('GCA_000001.1')
        gone = self.genome_dir('GCA_000003.1', fasta=False)

        self.run_trnascan([('GCA_000001.1', whole), ('GCA_000003.1', gone)])

        self.assertIsNotNone(self.model_flag(whole, 'GCA_000001.1'))

    def test_all_does_not_hand_it_over_either(self):
        """--all used to queue a genome whose FASTA is not there, and the scan of
        a file that does not exist took the whole run down."""
        whole = self.genome_dir('GCA_000001.1')
        gone = self.genome_dir('GCA_000003.1', fasta=False)

        _, ok = self.run_trnascan([('GCA_000001.1', whole), ('GCA_000003.1', gone)],
                                  all_genomes=True)

        self.assertTrue(ok)
        self.assertEqual(self.release_report(),
                         [('GCA_000003.1', T.REASON_NO_GENOMIC_FASTA)])


# ------------------------------------------------------------ a scan that fails

class WhenTheScanFails(TempDirCase):
    """One genome tRNAscan-SE refuses is one genome, not a batch."""

    def test_the_genome_is_named_and_the_batch_still_succeeds(self):
        whole = self.genome_dir('GCA_000001.1')
        bad = self.genome_dir('GCA_000004.1', refuse=True)

        _, ok = self.run_trnascan([('GCA_000001.1', whole), ('GCA_000004.1', bad)])

        self.assertTrue(ok)
        self.assertEqual(self.release_report(),
                         [('GCA_000004.1', T.REASON_TRNASCAN_FAILED)])
        self.assertEqual(B.batch_state(self.batches()[0]), B.STATE_SUCCESS)

    def test_a_batch_in_which_every_genome_failed_is_failed(self):
        """That is not a batch of difficult genomes; it is tRNAscan-SE not working."""
        bad = self.genome_dir('GCA_000004.1', refuse=True)
        worse = self.genome_dir('GCA_000005.1', refuse=True)

        _, ok = self.run_trnascan([('GCA_000004.1', bad), ('GCA_000005.1', worse)])

        self.assertFalse(ok)
        self.assertEqual(B.batch_state(self.batches()[0]), B.STATE_FAILED)

    def test_a_lone_failing_genome_does_not_fail_the_batch(self):
        """A batch whose last unscanned genome is a bad one would otherwise fail
        identically on every retry, which is what naming it exists to avoid."""
        done = self.genome_dir('GCA_000001.1', trna='valid')
        bad = self.genome_dir('GCA_000004.1', refuse=True)

        _, ok = self.run_trnascan([('GCA_000001.1', done), ('GCA_000004.1', bad)])

        self.assertTrue(ok)
        self.assertEqual(self.release_report(),
                         [('GCA_000004.1', T.REASON_TRNASCAN_FAILED)])


# --------------------------------------------------------- batches and the release

class TheBatchesAndTheRelease(TempDirCase):
    """--out_dir holds the state of the run; the tRNAs go to the genomes."""

    def setUp(self):
        super().setUp()
        self.genomes = [(accession, self.genome_dir(accession))
                        for accession in ('GCA_000001.1', 'GCA_000002.1',
                                          'GCA_000003.1', 'GCA_000004.1')]

    def test_the_release_is_cut_into_batches_of_batch_size(self):
        self.run_trnascan(self.genomes, batch_size=2)

        self.assertEqual(len(self.batches()), 2)

    def test_each_batch_records_what_it_came_to(self):
        self.run_trnascan(self.genomes, batch_size=2)

        canary = B.read_canary(os.path.join(self.batches()[0], B.SUCCESS_CANARY))
        self.assertEqual(canary['scanned'], '2')

    def test_no_trnas_are_written_to_the_out_dir(self):
        self.run_trnascan(self.genomes, batch_size=2)

        top = sorted(name for name in os.listdir(self.out_dir)
                     if not name.startswith(B.BATCH_DIR_PREFIX))
        self.assertEqual(top, [T.NOT_SCANNED_RELEASE_NAME])

    def test_the_release_report_is_written_even_with_nothing_in_it(self):
        self.run_trnascan(self.genomes, batch_size=2)

        self.assertEqual(self.release_report(), [])

    def test_a_finished_batch_is_skipped_by_a_later_run(self):
        self.run_trnascan(self.genomes, batch_size=2)
        for accession, gpath in self.genomes:
            os.remove(os.path.join(gpath, T.TRNA_DIR, accession + T.TRNA_LOG_EXT))

        self.run_trnascan(self.genomes, batch_size=2)

        self.assertIsNone(self.model_flag(*reversed(self.genomes[0])))

    def test_rerunning_a_finished_release_does_not_fall_over(self):
        """It used to divide by the number of genomes left to scan, which is zero."""
        self.run_trnascan(self.genomes, batch_size=2)
        for batch in self.batches():
            os.remove(os.path.join(batch, B.SUCCESS_CANARY))

        _, ok = self.run_trnascan(self.genomes, batch_size=2)

        self.assertTrue(ok)

    def test_a_batch_held_by_another_machine_is_left_alone(self):
        self.run_trnascan(self.genomes, batch_size=2)
        batch = self.batches()[0]
        os.remove(os.path.join(batch, B.SUCCESS_CANARY))
        with open(os.path.join(batch, B.RUNNING_CANARY), 'w') as handle:
            handle.write(B.canary_payload(host='another-machine', pid='1'))

        self.run_trnascan(self.genomes, batch_size=2)

        self.assertEqual(B.batch_state(batch), B.STATE_RUNNING)


if __name__ == '__main__':
    unittest.main()
