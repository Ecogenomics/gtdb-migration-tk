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
import pickle
import shutil
import tempfile
import unittest
from unittest import mock

from gtdb_migration_tk import batching as B
from gtdb_migration_tk import trnascan_manager as T
from gtdb_migration_tk.ncbi_utils import assembly_stats
from gtdb_migration_tk.utils import common as C


# What the stub was asked to do, written into the log file tRNAscan-SE would
# write, since a forked worker cannot hand an object back to the test.
FLAG_LINE = 'model flag: {}\n'

# A genome whose sequences say this is one the stub refuses, standing in for the
# genomes the real tRNAscan-SE cannot process.
REFUSE = 'REFUSE'


# What the stubbed tRNAscan-SE says it is.
VERSION = 'tRNAscan-SE 2.0.12 (Nov 2022)'


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
        for module in (T, B, C):
            patch = mock.patch.object(module, 'tqdm', QuietTqdm)
            patch.start()
            self.addCleanup(patch.stop)

        # tRNAscan-SE is not installed for the tests, so it cannot be asked
        patch = mock.patch.object(T, 'record_program_version', return_value=VERSION)
        patch.start()
        self.addCleanup(patch.stop)

        logger = logging.getLogger('timestamp')
        self.addCleanup(logger.setLevel, logger.level)
        logger.setLevel(logging.CRITICAL)

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    # ------------------------------------------------------------- the release

    def domain_file(self, predicted):
        """The GTDB domain report, keyed by GTDB's own genome ids.

        predicted maps an accession to 'ar', 'bac' or None, None being a genome
        the marker genes gave no prediction for.
        """
        path = os.path.join(self.dir, 'gtdb_domain_report.tsv')
        with open(path, 'w') as handle:
            handle.write('Genome Id\tGenome size\t{}\n'.format(T.DOMAIN_FILE_DOMAIN))
            for accession, domain in predicted.items():
                prefix = 'RS_' if accession.startswith('GCF') else 'GB_'
                handle.write('{}{}\t4000000\t{}\n'.format(
                    prefix, accession,
                    {'ar': 'd__Archaea', 'bac': 'd__Bacteria'}.get(domain,
                                                                   T.NO_PREDICTION)))
        return path

    def taxonomy_file(self, lineages):
        """The standardised NCBI taxonomy: accession and lineage, no header."""
        path = os.path.join(self.dir, 'ncbi_taxonomy.tsv')
        with open(path, 'w') as handle:
            for accession, domain in lineages.items():
                handle.write('{}\td__{};p__Whatever;c__;o__;f__;g__;s__\n'.format(
                    accession, 'Archaea' if domain == 'ar' else 'Bacteria'))
        return path

    def inputs(self, predicted=None, taxonomy=None):
        """The two files the command reads the domain of each genome from."""
        return dict(gtdb_domain_file=self.domain_file(predicted or {}),
                    taxonomy_file=self.taxonomy_file(taxonomy or {}))

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

    def scanner(self, predicted=None, taxonomy=None, **kwargs):
        with mock.patch.object(T, 'check_dependencies'):
            return T.tRNAScan(cpus=1, tmp_dir=os.path.join(self.dir, 'tmp'),
                              **dict(self.inputs(predicted, taxonomy), **kwargs))

    def run_trnascan(self, genomes, predicted=None, taxonomy=None,
                     all_genomes=False, batch_size=B.DEFAULT_BATCH_SIZE):
        scanner = self.scanner(predicted, taxonomy, batch_size=batch_size)
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


# ------------------------------------------------------- where the domain comes from

class TheDomainOfEachGenome(TempDirCase):
    """tRNAscan-SE searches with a bacterial or an archaeal model and the two give
    different answers, so this decides what the genome's tRNAs come out as."""

    def test_gtdb_s_own_prediction_is_used(self):
        scanner = self.scanner(predicted={'GCA_000002.1': 'ar',
                                          'GCA_000001.1': 'bac'})

        self.assertEqual(scanner.domain_flag('GCA_000002.1'), T.ARCHAEAL_FLAG)
        self.assertEqual(scanner.domain_flag('GCA_000001.1'), T.BACTERIAL_FLAG)

    def test_the_gtdb_id_prefix_is_stripped(self):
        """GTDB names a genome GB_GCA_000002.1 where a genome_dirs file says
        GCA_000002.1."""
        scanner = self.scanner(predicted={'GCA_000002.1': 'ar'})

        self.assertIn('GCA_000002.1', scanner.domains)

    def test_the_ncbi_taxonomy_answers_for_a_genome_gtdb_could_not_predict(self):
        """'None' is what the column holds when the markers gave no answer."""
        scanner = self.scanner(predicted={'GCA_000002.1': None},
                               taxonomy={'GCA_000002.1': 'ar'})

        self.assertEqual(scanner.domain_flag('GCA_000002.1'), T.ARCHAEAL_FLAG)

    def test_the_ncbi_taxonomy_answers_for_a_genome_the_domain_file_omits(self):
        scanner = self.scanner(taxonomy={'GCA_000002.1': 'ar'})

        self.assertEqual(scanner.domain_flag('GCA_000002.1'), T.ARCHAEAL_FLAG)

    def test_gtdb_s_prediction_wins_over_the_ncbi_taxonomy(self):
        """Which is the point of preferring it: it is made from the genome, and
        catches one NCBI has filed under the wrong domain."""
        scanner = self.scanner(predicted={'GCA_000002.1': 'ar'},
                               taxonomy={'GCA_000002.1': 'bac'})

        self.assertEqual(scanner.domain_flag('GCA_000002.1'), T.ARCHAEAL_FLAG)

    def test_a_genbank_genome_finds_the_taxonomy_of_its_refseq_counterpart(self):
        scanner = self.scanner(taxonomy={'GCF_000002.1': 'ar'})

        self.assertEqual(scanner.domain_flag('GCA_000002.1'), T.ARCHAEAL_FLAG)

    def test_a_genome_neither_file_answers_for_is_scanned_as_a_bacterium(self):
        scanner = self.scanner(predicted={'GCA_000001.1': 'bac'})

        self.assertEqual(scanner.domain_flag('GCA_000009.1'), T.BACTERIAL_FLAG)

    def test_the_model_each_genome_was_scanned_with(self):
        arc = self.genome_dir('GCA_000002.1')
        bac = self.genome_dir('GCA_000001.1')
        unknown = self.genome_dir('GCA_000009.1')

        self.run_trnascan([('GCA_000002.1', arc), ('GCA_000001.1', bac),
                           ('GCA_000009.1', unknown)],
                          predicted={'GCA_000002.1': 'ar', 'GCA_000001.1': 'bac'})

        self.assertEqual(self.model_flag(arc, 'GCA_000002.1'), T.ARCHAEAL_FLAG)
        self.assertEqual(self.model_flag(bac, 'GCA_000001.1'), T.BACTERIAL_FLAG)
        self.assertEqual(self.model_flag(unknown, 'GCA_000009.1'), T.BACTERIAL_FLAG)

    def test_a_genome_of_no_known_domain_is_warned_about(self):
        unknown = self.genome_dir('GCA_000009.1')

        with self.assertLogs('timestamp', level='WARNING') as captured:
            self.run_trnascan([('GCA_000009.1', unknown)])

        self.assertTrue(any('no domain in either' in record.getMessage()
                            for record in captured.records), captured.records)


# ---------------------------------------------------------- what decides the work

class WhatDecidesTheWork(TempDirCase):
    """The checksum beside a genome's own results, and not a report."""

    def decide(self, accession, **kwargs):
        gpath = self.genome_dir(accession, **kwargs)
        job = T.TrnaJob(accession,
                        os.path.join(gpath, os.path.basename(gpath) + '_genomic.fna.gz'),
                        T.BACTERIAL_FLAG)
        return T.trnascan_parser(job)

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


class WhatThePoolIsHanded(TempDirCase):
    """Each genome, and not the release's domain table with it.

    The pool pickles what it is handed once per genome, one at a time in the
    parent. Handed a bound method, it pickled the instance and the 2.7M-entry
    domain table in it -- 47 MB and 0.6 s a genome -- and a batch spent hours
    'Checking tRNAs' that is seconds of work.
    """

    # a genome's task is its function, its accession, two paths and a flag
    MOST_A_TASK_SHOULD_WEIGH = 10000

    def test_neither_pass_hands_a_worker_the_domain_table(self):
        genomes = [('GCA_00000{}.1'.format(i),
                    self.genome_dir('GCA_00000{}.1'.format(i))) for i in (1, 2)]
        scanner = self.scanner()
        # a table of a size that would weigh far more than a task, were it sent
        scanner.domains = {'GCA_{:09d}.1'.format(i): T.DOMAIN_BACTERIA
                           for i in range(100000)}

        sizes = []

        class RecordingPool(object):
            """Runs the work in this process, weighing what each genome sends."""

            def __init__(self, processes=None):
                pass

            def __enter__(self):
                return self

            def __exit__(self, *exc):
                return False

            def imap_unordered(self, func, items):
                for item in items:
                    sizes.append(len(pickle.dumps((func, item))))
                    yield func(item)

        with mock.patch.object(T.mp, 'Pool', RecordingPool), \
                mock.patch.object(T.subprocess, 'Popen', StubPopen):
            ok = scanner.run(self.genome_dirs_file(genomes), self.out_dir)

        self.assertTrue(ok)
        # both passes, both genomes: the check, and the scan it decided on
        self.assertEqual(len(sizes), 4)
        self.assertLess(max(sizes), self.MOST_A_TASK_SHOULD_WEIGH)


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


# -------------------------------------------------- what the tRNAs were made with

class TheVersionThatScannedAGenome(TempDirCase):
    """Recorded beside the tRNAs, since they outlive the run that made them."""

    def version_file(self, gpath):
        return os.path.join(gpath, T.TRNA_DIR, 'trnascan-se.version')

    def test_a_scanned_genome_records_the_version_that_scanned_it(self):
        gpath = self.genome_dir('GCA_000001.1')

        self.run_trnascan([('GCA_000001.1', gpath)])

        with open(self.version_file(gpath)) as handle:
            self.assertEqual(handle.read(), VERSION + '\n')

    def test_the_version_is_asked_of_trnascan_se_once_for_the_run(self):
        """Not once per genome: a million-odd questions with one answer."""
        genomes = [('GCA_00000{}.1'.format(i),
                    self.genome_dir('GCA_00000{}.1'.format(i))) for i in (1, 2, 3)]

        with mock.patch.object(T, 'record_program_version',
                               return_value=VERSION) as asked:
            scanner = self.scanner()
        with mock.patch.object(T.subprocess, 'Popen', StubPopen):
            scanner.run(self.genome_dirs_file(genomes), self.out_dir)

        asked.assert_called_once_with('tRNAscan-SE')

    def test_a_genome_trnascan_se_failed_on_records_no_version(self):
        """Nothing was made, so there is nothing for a version to vouch for."""
        whole = self.genome_dir('GCA_000001.1')
        bad = self.genome_dir('GCA_000004.1', refuse=True)

        self.run_trnascan([('GCA_000001.1', whole), ('GCA_000004.1', bad)])

        self.assertFalse(os.path.exists(self.version_file(bad)))

    def test_a_genome_skipped_keeps_the_version_it_was_scanned_with(self):
        """Its tRNAs were not made by this run, so this run's version would say
        something untrue about them."""
        gpath = self.genome_dir('GCA_000001.1', trna='valid')
        with open(self.version_file(gpath), 'w') as handle:
            handle.write('tRNAscan-SE 2.0.9 (July 2021)\n')

        self.run_trnascan([('GCA_000001.1', gpath)])

        with open(self.version_file(gpath)) as handle:
            self.assertEqual(handle.read(), 'tRNAscan-SE 2.0.9 (July 2021)\n')


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


# ------------------------------------------------------------ a genome too large

# GCA_964261755.1, a faecal metagenome deposited as one genome
METAGENOME_BASES = 9528631298


def state_size(gpath, bases):
    """Write the assembly statistics NCBI publishes beside a genome, stating its size."""
    with open(assembly_stats(gpath), 'w') as handle:
        handle.write('all\tall\tall\tall\ttotal-length\t{}\n'.format(bases))


class AGenomeTooLargeToScan(TempDirCase):
    """Named and left, rather than holding its batch for as long as it takes."""

    def run_with_limit(self, genomes, max_genome_size=C.DEFAULT_MAX_GENOME_SIZE,
                       all_genomes=False):
        scanner = self.scanner(max_genome_size=max_genome_size)
        with mock.patch.object(T.subprocess, 'Popen', StubPopen):
            return scanner.run(self.genome_dirs_file(genomes), self.out_dir, all_genomes)

    def test_it_is_named_in_the_release_report_and_not_scanned(self):
        small = self.genome_dir('GCA_000001.1')
        large = self.genome_dir('GCA_000003.1')
        state_size(large, METAGENOME_BASES)

        ok = self.run_with_limit([('GCA_000001.1', small), ('GCA_000003.1', large)])

        self.assertTrue(ok)
        self.assertEqual(self.release_report(),
                         [('GCA_000003.1', T.REASON_GENOME_TOO_LARGE)])
        self.assertIsNone(self.model_flag(large, 'GCA_000003.1'))
        self.assertIsNotNone(self.model_flag(small, 'GCA_000001.1'))

    def test_the_limit_is_the_one_given(self):
        large = self.genome_dir('GCA_000003.1')
        state_size(large, METAGENOME_BASES)

        self.run_with_limit([('GCA_000003.1', large)], max_genome_size=10000)

        self.assertIsNotNone(self.model_flag(large, 'GCA_000003.1'))
        self.assertEqual(self.release_report(), [])

    def test_the_default_limit_is_100_mbp(self):
        self.assertEqual(self.scanner().max_genome_bases, 100 * 1000 * 1000)

    def test_a_genome_whose_trnas_are_already_there_is_not_named(self):
        """It has not been left out of anything."""
        large = self.genome_dir('GCA_000003.1', trna='valid')
        state_size(large, METAGENOME_BASES)

        self.run_with_limit([('GCA_000003.1', large)])

        self.assertEqual(self.release_report(), [])

    def test_all_does_not_hand_it_over_either(self):
        large = self.genome_dir('GCA_000003.1', trna='valid')
        state_size(large, METAGENOME_BASES)

        self.run_with_limit([('GCA_000003.1', large)], all_genomes=True)

        self.assertEqual(self.release_report(),
                         [('GCA_000003.1', T.REASON_GENOME_TOO_LARGE)])
        self.assertIsNone(self.model_flag(large, 'GCA_000003.1'))

    def test_the_batch_counts_it_as_not_scanned(self):
        large = self.genome_dir('GCA_000003.1')
        state_size(large, METAGENOME_BASES)

        self.run_with_limit([('GCA_000003.1', large)])

        canary = B.read_canary(os.path.join(self.batches()[0], B.SUCCESS_CANARY))
        self.assertEqual((canary['scanned'], canary['not_scanned']), ('0', '1'))


if __name__ == '__main__':
    unittest.main()
