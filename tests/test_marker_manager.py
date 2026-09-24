#!/usr/bin/env python3
"""Offline unit tests for marker_manager.py -- no HMM is ever searched.

Run with the interpreter that has tqdm:

    /opt/centos7/sw/miniconda3/envs/python3/bin/python -m unittest -v tests.test_marker_manager

What is tested here is everything hmmsearch decides before an HMM is searched:
how the release is cut into batches and shared between machines, which genomes of
a batch still need their markers found, and what is handed to the workers.
search_markers() is replaced by a recorder, so the real run_hmmsearch() runs over
real batch directories and what it would have searched can be read back.

What a search coming to nothing costs is tested here too, with no HMM searched:
that --hmm_db_path is refused at the start where it is not the HMMs --db will be
searched against, and that a batch whose workers died is failed rather than
finished. r237 was searched against a directory hmmsearch could not read, and all
135 batches were marked SUCCESS with no marker table written by any of them.

The batching machinery itself is tested in tests/test_batching.py. What is tested
here is this command's use of it: that a Pfam run and a TIGRFAM run of one output
directory do not read each other's canaries, and that the work list can hold
nothing a worker would mistake for the end of the work.
"""

import gzip
import os
import shutil
import subprocess
import tempfile
import unittest
from unittest import mock

from gtdb_migration_tk import batching as B
from gtdb_migration_tk import marker_manager as M


SUFFIX = '33.1_lite'
MARKER_DIR = 'pfam_' + SUFFIX
MARKER_EXT = '_pfam_{}.tsv'.format(SUFFIX)
TIGR_SUFFIX = '15.0_lite'


def dying_worker(queue_in, queue_out, dir_suffix):
    """A worker whose search cannot run, as all 96 of them could not in r237."""

    raise RuntimeError('hmmsearch exited 1: the HMM file appears to be empty')


def exiting_worker(queue_in, queue_out, dir_suffix):
    """A worker that leaves by sys.exit(), which is how PfamScan gives up."""

    raise SystemExit(1)


def quiet_worker(queue_in, queue_out, dir_suffix):
    """A worker that does the work and ends as one should."""

    while True:
        item = queue_in.get()
        if item is None:
            break
        queue_out.put(item)


# What the stubbed HMMER says it is.
VERSION = 'HMMER 3.4 (Aug 2023)'


class TempDirCase(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='marker_manager_test.')
        self.out_dir = os.path.join(self.dir, 'out')
        os.makedirs(self.out_dir)
        self.searched = []

        # HMMER is not installed for the tests, so it cannot be asked
        patch = mock.patch.object(M, 'record_program_version', return_value=VERSION)
        patch.start()
        self.addCleanup(patch.stop)

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def hmm_db(self, db='pfam', magic=b'HMMER3/f [3.1b2]\n'):
        """An HMM database of the shape --db asks for, holding nothing to search.

        marker_setup() checks --hmm_db_path before a batch is claimed, so a test
        that runs the command needs the directory Pfam is given or the file
        TIGRFAM is, and both are read for their first five bytes alone.

        @return: what --hmm_db_path would be.
        """

        root = os.path.join(self.dir, 'hmms')
        if not os.path.isdir(root):
            os.makedirs(root)

        name = M.PFAM_LIBRARY if db == 'pfam' else 'tigrfam.hmm'
        library = os.path.join(root, name)
        with open(library, 'wb') as handle:
            handle.write(magic)

        return root if db == 'pfam' else library

    def manager(self, cpus=1, batch_size=10000, tmp_dir=None, **kwargs):
        # the real one checks for prodigal and hmmsearch on PATH and exits without them
        with mock.patch.object(M, 'check_dependencies'):
            return M.MarkerManager(tmp_dir=tmp_dir or self.dir, cpus=cpus,
                                   batch_size=batch_size, **kwargs)

    def genome(self, gid, proteins=b'>gene\nMA\n', annotated=False, checksum=True,
               marker_dir=MARKER_DIR, marker_ext=MARKER_EXT):
        """A genome directory, optionally already carrying its marker table.

        Parameters
        ----------
        gid : str
            Accession, which names the files within.
        proteins : bytes or None
            Contents of the protein FASTA, gzipped; None writes no protein file
            at all, and b'' writes a zero byte one, which is what split_by_fasta()
            means by missing -- it stats the file, and a gzip of nothing is still
            twenty-odd bytes of header.
        annotated : bool
            Whether a marker table is already there.
        checksum : bool or str
            True writes the table's .sha256 correctly, False writes none at all,
            and a string is written as the digest -- a table that disagrees with
            its checksum rather than one with none, which are different states on
            disk and the same answer.

        @return: the genome directory.
        """

        path = os.path.join(self.dir, 'release', gid)
        prodigal_dir = os.path.join(path, 'prodigal')
        os.makedirs(prodigal_dir)

        if proteins == b'':
            open(os.path.join(prodigal_dir, gid + '_protein.faa.gz'), 'wb').close()
        elif proteins is not None:
            with gzip.open(os.path.join(prodigal_dir, gid + '_protein.faa.gz'), 'wb') as handle:
                handle.write(proteins)

        if annotated:
            table_dir = os.path.join(prodigal_dir, marker_dir)
            os.makedirs(table_dir)
            table = os.path.join(table_dir, gid + marker_ext)
            with gzip.open(table + '.gz', 'wb') as handle:
                handle.write(b'# hits\n')
            if checksum is True:
                # of the UNCOMPRESSED bytes, as marker_parser() reads it back
                with open(table + '.gz', 'rb') as raw:
                    digest = M.sha256_rb(gzip.GzipFile(fileobj=raw))
                with open(table + '.sha256', 'w') as handle:
                    handle.write(digest + '\n')
            elif checksum:
                with open(table + '.sha256', 'w') as handle:
                    handle.write(checksum + '\n')

        return path

    def inputs(self, genomes, outcomes=None):
        """The two files run_hmmsearch() reads.

        @return: (genome_dirs file, report file).
        """

        outcomes = outcomes or {gid: 'new' for gid in genomes}
        dirs_file = os.path.join(self.dir, 'genome_dirs.tsv')
        with open(dirs_file, 'w') as handle:
            for gid, path in genomes.items():
                handle.write('{}\t{}\tG{}\n'.format(gid, path, gid[4:]))

        report = os.path.join(self.dir, 'report.log')
        with open(report, 'w') as handle:
            for gid, outcome in outcomes.items():
                handle.write('{}\t{}\n'.format(gid, outcome))

        return dirs_file, report

    def run_hmmsearch(self, genomes, outcomes=None, db='pfam', suffix=SUFFIX,
                      manager=None, searcher=None, **kwargs):
        """Run the real run_hmmsearch() with the searching replaced by a recorder.

        @return: what run_hmmsearch() returned.
        """

        dirs_file, report = self.inputs(genomes, outcomes)
        manager = manager or self.manager()

        def record(manager, genome_files, worker, dir_suffix):
            self.searched.append(list(genome_files))

        with mock.patch.object(M.MarkerManager, 'search_markers',
                               searcher or record):
            return manager.run_hmmsearch(dirs_file, report, db, suffix,
                                         self.hmm_db(db), self.out_dir, **kwargs)

    def state_dir(self, marker_dir=MARKER_DIR):
        return os.path.join(self.out_dir, marker_dir)

    def batches(self, marker_dir=MARKER_DIR):
        return B.batch_dir_names(self.state_dir(marker_dir), M.LAYOUT)

    def queued(self):
        """Every genome handed to the searching, across all batches."""
        return sorted(path for batch in self.searched for path in batch)


class TheBatchesOfARun(TempDirCase):
    """Where a run's batches live, and what their batchfiles name."""

    def test_the_batches_are_under_the_marker_directory_of_the_output_directory(self):
        # so that --out_dir is one directory for a release rather than one per
        # database, which five servers would each have to be told the right one of
        genomes = {'GCF_000000001.1': self.genome('GCF_000000001.1')}

        self.run_hmmsearch(genomes)

        self.assertTrue(os.path.isdir(self.state_dir()), os.listdir(self.out_dir))
        self.assertEqual([os.path.basename(b) for b in self.batches()],
                         ['batch_000001'])

    def test_a_batchfile_names_each_genomes_proteins_not_its_genomic_fasta(self):
        # this command reads the proteins prodigal called; the genomic FASTA is
        # what trans_table and prodigal read
        path = self.genome('GCF_000000001.1')

        self.run_hmmsearch({'GCF_000000001.1': path})

        rows = B.read_batchfile(B.batchfile_path(self.batches()[0], M.LAYOUT))
        self.assertEqual(rows, [(os.path.join(
            path, 'prodigal', 'GCF_000000001.1_protein.faa.gz'), 'GCF_000000001.1')])

    def test_the_release_is_cut_into_batches_of_the_size_asked_for(self):
        genomes = {'GCF_00000000{}.1'.format(n): self.genome('GCF_00000000{}.1'.format(n))
                   for n in range(1, 6)}

        self.run_hmmsearch(genomes, manager=self.manager(batch_size=2))

        self.assertEqual(len(self.batches()), 3)

    def test_a_plan_already_there_is_used_rather_than_made_again(self):
        # another machine is working from it, and repartitioning a release that
        # has since gained a genome would move genomes between finished batches
        genomes = {'GCF_00000000{}.1'.format(n): self.genome('GCF_00000000{}.1'.format(n))
                   for n in range(1, 6)}
        self.run_hmmsearch(genomes, manager=self.manager(batch_size=2))

        self.run_hmmsearch(genomes, manager=self.manager(batch_size=5))

        self.assertEqual(len(self.batches()), 3)


class TwoDatabasesOneOutputDirectory(TempDirCase):
    """A Pfam run and a TIGRFAM run of one --out_dir are different work.

    They cover the same genomes and finish at different times, so a batch that is
    done for one is not done for the other. Sharing batch directories would have
    the SUCCESS canary of a Pfam batch tell a TIGRFAM run it had nothing to do.
    """

    def test_each_database_gets_its_own_batches(self):
        genomes = {'GCF_000000001.1': self.genome('GCF_000000001.1')}

        self.run_hmmsearch(genomes, db='pfam', suffix=SUFFIX)
        self.run_hmmsearch(genomes, db='tigrfam', suffix=TIGR_SUFFIX)

        self.assertEqual(sorted(os.listdir(self.out_dir)),
                         ['pfam_' + SUFFIX, 'tigrfam_' + TIGR_SUFFIX])

    def test_a_finished_pfam_batch_does_not_finish_the_tigrfam_one(self):
        genomes = {'GCF_000000001.1': self.genome('GCF_000000001.1')}
        self.run_hmmsearch(genomes, db='pfam', suffix=SUFFIX)
        self.assertEqual(len(self.queued()), 1)

        self.searched = []
        self.run_hmmsearch(genomes, db='tigrfam', suffix=TIGR_SUFFIX)

        self.assertEqual(len(self.queued()), 1)

    def test_one_version_of_a_database_does_not_finish_another(self):
        # annotating against a new Pfam release is new work over the same genomes
        genomes = {'GCF_000000001.1': self.genome('GCF_000000001.1')}
        self.run_hmmsearch(genomes, db='pfam', suffix=SUFFIX)

        self.searched = []
        self.run_hmmsearch(genomes, db='pfam', suffix='37.0_lite')

        self.assertEqual(len(self.queued()), 1)
        self.assertIn('pfam_37.0_lite', os.listdir(self.out_dir))


class SharingTheReleaseBetweenMachines(TempDirCase):
    """What a second run over the same output directory does."""

    def held_batch(self, genomes, batch_size=1):
        """Plan the batches and let another machine hold the first.

        @return: (genome_dirs file, report file).
        """

        dirs_file, report = self.inputs(genomes)
        B.plan_batches(dirs_file, self.state_dir(), batch_size, M.LAYOUT,
                       self.manager().logger, genome_file=M.protein_fasta)
        with open(os.path.join(self.batches()[0], B.RUNNING_CANARY), 'w') as handle:
            handle.write('host\tsome-other-machine\npid\t1\ntime\tnow\n')
        return dirs_file, report

    def test_a_finished_batch_is_not_done_again(self):
        genomes = {'GCF_000000001.1': self.genome('GCF_000000001.1')}
        self.run_hmmsearch(genomes)
        self.assertEqual(B.batch_state(self.batches()[0]), B.STATE_SUCCESS)

        self.searched = []
        self.run_hmmsearch(genomes)

        self.assertEqual(self.searched, [])

    def test_a_batch_another_machine_holds_is_left_to_it(self):
        genomes = {'GCF_00000000{}.1'.format(n): self.genome('GCF_00000000{}.1'.format(n))
                   for n in range(1, 3)}
        self.held_batch(genomes)

        self.run_hmmsearch(genomes, manager=self.manager(batch_size=1))

        # the held batch was skipped and the other was done
        self.assertEqual(len(self.searched), 1)
        self.assertEqual(B.batch_state(self.batches()[0]), B.STATE_RUNNING)
        self.assertEqual(B.batch_state(self.batches()[1]), B.STATE_SUCCESS)

    def test_a_batch_whose_search_died_is_failed_and_not_called_a_success(self):
        # the searching used to swallow every exception, which would have a batch
        # that searched nothing be written a SUCCESS no machine looks behind
        genomes = {'GCF_000000001.1': self.genome('GCF_000000001.1')}

        def explode(manager, genome_files, worker, dir_suffix):
            raise RuntimeError('hmmsearch died')

        finished = self.run_hmmsearch(genomes, searcher=explode)

        self.assertFalse(finished)
        self.assertEqual(B.batch_state(self.batches()[0]), B.STATE_FAILED)

    def test_a_failed_batch_is_retried_by_the_next_run(self):
        genomes = {'GCF_000000001.1': self.genome('GCF_000000001.1')}

        def explode(manager, genome_files, worker, dir_suffix):
            raise RuntimeError('hmmsearch died')

        self.run_hmmsearch(genomes, searcher=explode)
        self.assertTrue(self.run_hmmsearch(genomes))

        self.assertEqual(B.batch_state(self.batches()[0]), B.STATE_SUCCESS)
        self.assertEqual(len(self.queued()), 1)

    def test_a_run_that_finished_every_batch_says_so(self):
        genomes = {'GCF_000000001.1': self.genome('GCF_000000001.1')}

        self.assertTrue(self.run_hmmsearch(genomes))


class WhatReachesTheWorkers(TempDirCase):
    """The work list holds no genome a worker would mistake for the end of it.

    search_markers() ends its queue with one None per worker, and every worker
    breaks on the first None it draws. A None among the genomes is a second stop
    signal: the worker that takes it exits with the rest of the batch still
    queued, and does so silently.
    """

    def test_a_genome_needing_markers_is_searched(self):
        path = self.genome('GCF_000000001.1')

        self.run_hmmsearch({'GCF_000000001.1': path})

        self.assertEqual(self.queued(), [os.path.join(
            path, 'prodigal', 'GCF_000000001.1_protein.faa.gz')])

    def test_a_genome_already_annotated_is_left_out_of_the_work_list(self):
        genomes = {'GCF_000000001.1': self.genome('GCF_000000001.1'),
                   'GCF_000000002.1': self.genome('GCF_000000002.1', annotated=True)}

        self.run_hmmsearch(genomes)

        self.assertEqual(len(self.queued()), 1)
        self.assertNotIn(None, self.queued())

    def test_a_genome_with_no_proteins_is_left_out_of_the_work_list(self):
        genomes = {'GCF_000000001.1': self.genome('GCF_000000001.1'),
                   'GCF_000000002.1': self.genome('GCF_000000002.1', proteins=None)}

        self.run_hmmsearch(genomes)

        self.assertEqual(len(self.queued()), 1)
        self.assertNotIn(None, self.queued())

    def test_a_genome_with_empty_proteins_is_left_out_of_the_work_list(self):
        genomes = {'GCF_000000001.1': self.genome('GCF_000000001.1'),
                   'GCF_000000002.1': self.genome('GCF_000000002.1', proteins=b'')}

        self.run_hmmsearch(genomes)

        self.assertEqual(len(self.queued()), 1)
        self.assertNotIn(None, self.queued())

    def test_a_batch_with_nothing_to_do_hands_over_an_empty_work_list(self):
        genomes = {'GCF_000000001.1': self.genome('GCF_000000001.1', annotated=True),
                   'GCF_000000002.1': self.genome('GCF_000000002.1', proteins=None)}

        self.run_hmmsearch(genomes)

        self.assertEqual(self.searched, [[]])

    def test_all_genomes_searches_what_is_already_annotated(self):
        genomes = {'GCF_000000001.1': self.genome('GCF_000000001.1', annotated=True)}

        self.run_hmmsearch(genomes, all_genomes=True)

        self.assertEqual(len(self.queued()), 1)

    def test_all_genomes_does_not_search_a_genome_with_no_proteins(self):
        # there is nothing to search whatever the run was told; --all discards
        # results, it does not conjure proteins
        genomes = {'GCF_000000001.1': self.genome('GCF_000000001.1', proteins=None)}

        self.run_hmmsearch(genomes, all_genomes=True)

        self.assertEqual(self.queued(), [])

    def test_all_genomes_does_a_batch_that_already_finished(self):
        genomes = {'GCF_000000001.1': self.genome('GCF_000000001.1')}
        self.run_hmmsearch(genomes)

        self.searched = []
        self.run_hmmsearch(genomes, all_genomes=True)

        self.assertEqual(len(self.queued()), 1)


class WhichGenomesMarkerParserSkips(TempDirCase):
    """marker_parser() decides one genome, and says so with one value.

    Skipping used to be said two ways -- the string 'null' for a genome already
    annotated, and falling off the end as None for one with no protein file --
    and only the first was filtered out of the work list.
    """

    def parse(self, gid, consider=True, all_genomes=False, **kwargs):
        path = self.genome(gid, **kwargs)
        proteins = M.protein_fasta(gid, path)
        job = (gid, proteins, MARKER_DIR, MARKER_EXT,
               {gid} if consider else set(), 'Pfam', all_genomes)
        return self.manager().marker_parser(job)

    def test_a_genome_to_annotate_gives_its_protein_file(self):
        result = self.parse('GCF_000000001.1')

        self.assertTrue(result.endswith('GCF_000000001.1_protein.faa.gz'), result)
        self.assertTrue(os.path.exists(result))

    def test_an_annotated_genome_with_a_matching_checksum_is_skipped(self):
        self.assertIsNone(self.parse('GCF_000000001.1', annotated=True))

    def test_an_annotated_genome_with_no_checksum_is_annotated_again(self):
        # nothing accounts for the table, so it is not results anyone can use
        self.assertIsNotNone(self.parse('GCF_000000001.1', annotated=True,
                                        checksum=False))

    def test_an_annotated_genome_whose_checksum_does_not_match_is_annotated_again(self):
        self.assertIsNotNone(self.parse('GCF_000000001.1', annotated=True,
                                        checksum='0' * 64))

    def test_a_genome_the_report_does_not_name_is_still_annotated(self):
        # the marker table decides the work; the report only says whether the two
        # agree, and a genome with no table is searched either way
        self.assertIsNotNone(self.parse('GCF_000000001.1', consider=False))

    def test_all_genomes_ignores_a_table_that_is_already_there(self):
        self.assertIsNotNone(self.parse('GCF_000000001.1', annotated=True,
                                        all_genomes=True))

    def test_skipping_is_one_value_so_nothing_can_be_dropped_by_halves(self):
        # the contract search_batch()'s filter rests on: there is no second skip
        # value for it to miss
        self.assertIsNone(self.parse('GCF_000000001.1', annotated=True))
        self.assertIsNone(self.parse('GCF_000000002.1', annotated=True,
                                     consider=False))


class WhatTheReportIsCrossCheckedAgainst(TempDirCase):
    """The report does not decide the work; it says whether disk agrees with it.

    A genome is in the report's genomes_to_regenerate() when the release did NOT
    carry its derived data across, so the report's claim is that it has no marker
    table yet. The marker table decides whether it is searched either way; where
    the two disagree the log says so, and it says which way round the
    disagreement is rather than asserting one of them.
    """

    def warnings(self, gid, consider, **kwargs):
        path = self.genome(gid, **kwargs)
        job = (gid, M.protein_fasta(gid, path), MARKER_DIR, MARKER_EXT,
               {gid} if consider else set(), 'Pfam', False)
        manager = self.manager()
        with self.assertLogs('timestamp', level='WARNING') as caught:
            manager.marker_parser(job)
            # assertLogs fails an empty block, so every case logs at least this
            manager.logger.warning('end of case')
        return ' '.join(caught.output)

    def test_an_annotated_genome_the_release_calls_new_is_flagged(self):
        # top left: the release says it carried nothing across, and the table is
        # there anyway. Expected of a batch being resumed, so it is said and the
        # genome is still skipped
        said = self.warnings('GCF_000000001.1', consider=True, annotated=True)

        self.assertIn('marked as new or modified, but already has Pfam annotations', said)
        self.assertIn('being skipped', said)

    def test_an_annotated_genome_the_release_expects_is_not_flagged(self):
        said = self.warnings('GCF_000000001.1', consider=False, annotated=True)

        self.assertEqual(said.count('WARNING'), 1)      # the marker of the block

    def test_an_unannotated_genome_the_release_expects_annotated_is_flagged(self):
        # bottom right: the release says this genome kept its derived data, and
        # the marker table is not there. The one cell worth a warning every time
        said = self.warnings('GCF_000000001.1', consider=False)

        self.assertIn('has no Pfam annotations, but is also not marked for processing', said)

    def test_an_unannotated_genome_the_release_calls_new_is_not_flagged(self):
        said = self.warnings('GCF_000000001.1', consider=True)

        self.assertEqual(said.count('WARNING'), 1)

    def test_an_unvouched_table_the_release_calls_new_says_so(self):
        # this used to say the genome "was not marked for reannotation" without
        # looking, which is false for exactly the genome it most often describes:
        # one whose annotation run was interrupted partway through writing it
        said = self.warnings('GCF_000000001.1', consider=True,
                             annotated=True, checksum=False)

        self.assertIn('no valid checksum, and is marked as new or modified', said)
        self.assertNotIn('not marked for reannotation', said)

    def test_an_unvouched_table_the_release_does_not_call_new_says_so(self):
        said = self.warnings('GCF_000000001.1', consider=False,
                             annotated=True, checksum=False)

        self.assertIn('no valid checksum, though it is not marked for reannotation', said)

    def test_a_table_that_disagrees_with_its_checksum_reads_the_same_way(self):
        # an absent .sha256 and one that disagrees are different states on disk
        # and the same answer: neither shows the annotations to be right
        said = self.warnings('GCF_000000001.1', consider=True,
                             annotated=True, checksum='0' * 64)

        self.assertIn('no valid checksum, and is marked as new or modified', said)
        self.assertIn('will be reannotated', said)


class TheGenomesThatGotNoMarkers(TempDirCase):
    """not_searched.tsv, per batch and gathered for the release."""

    def test_a_genome_with_no_proteins_is_named_in_the_batch(self):
        genomes = {'GCF_000000001.1': self.genome('GCF_000000001.1'),
                   'GCF_000000002.1': self.genome('GCF_000000002.1', proteins=None)}

        self.run_hmmsearch(genomes)

        with open(os.path.join(self.batches()[0], M.NOT_SEARCHED_NAME)) as handle:
            rows = [line.rstrip('\n').split('\t') for line in handle]
        self.assertEqual(rows, [list(M.NOT_SEARCHED_HEADER),
                                ['GCF_000000002.1', M.REASON_NO_PROTEINS]])

    def test_the_release_file_is_written_once_every_batch_is_done(self):
        genomes = {'GCF_00000000{}.1'.format(n): self.genome(
            'GCF_00000000{}.1'.format(n), proteins=None if n == 2 else b'>g\nMA\n')
            for n in range(1, 4)}

        self.run_hmmsearch(genomes, manager=self.manager(batch_size=1))

        path = os.path.join(self.state_dir(), M.NOT_SEARCHED_RELEASE_NAME)
        with open(path) as handle:
            rows = [line.rstrip('\n').split('\t') for line in handle]
        self.assertEqual(rows, [list(M.NOT_SEARCHED_HEADER),
                                ['GCF_000000002.1', M.REASON_NO_PROTEINS]])

    def test_the_release_file_is_absent_while_a_batch_is_unfinished(self):
        # either the whole release or nothing: a part of it reads like the whole
        genomes = {'GCF_00000000{}.1'.format(n): self.genome('GCF_00000000{}.1'.format(n))
                   for n in range(1, 3)}
        dirs_file, report = self.inputs(genomes)
        B.plan_batches(dirs_file, self.state_dir(), 1, M.LAYOUT,
                       self.manager().logger, genome_file=M.protein_fasta)
        with open(os.path.join(self.batches()[0], B.RUNNING_CANARY), 'w') as handle:
            handle.write('host\tsome-other-machine\npid\t1\ntime\tnow\n')

        self.run_hmmsearch(genomes, manager=self.manager(batch_size=1))

        self.assertFalse(os.path.exists(
            os.path.join(self.state_dir(), M.NOT_SEARCHED_RELEASE_NAME)))


class WhatABatchRecordsOfItself(TempDirCase):
    """The counts in the SUCCESS canary, which the release totals are added from."""

    def test_a_finished_batch_records_what_it_came_to(self):
        genomes = {'GCF_000000001.1': self.genome('GCF_000000001.1'),
                   'GCF_000000002.1': self.genome('GCF_000000002.1', annotated=True),
                   'GCF_000000003.1': self.genome('GCF_000000003.1', proteins=None)}

        self.run_hmmsearch(genomes)

        canary = B.read_canary(os.path.join(self.batches()[0], B.SUCCESS_CANARY))
        self.assertEqual(canary['searched'], '1')
        self.assertEqual(canary['already_searched'], '1')
        self.assertEqual(canary['not_searched'], '1')


class WhereTheScratchCopiesGo(TempDirCase):
    """--tmp_dir, which the workers ignored.

    Each genome's proteins are decompressed before they are searched, one
    directory per genome in flight. Both workers called tempfile.mkdtemp() with
    no dir, so every copy went to /tmp whatever --tmp_dir said -- on a machine
    running -c 40 that is forty decompressed proteomes in /tmp at a time, and a
    scratch directory chosen to keep them off it did nothing.
    """

    class Stop(Exception):
        """Raised in place of making the directory, to stop the worker there."""

    def scratch_dir_of(self, worker_name):
        """Where the named worker asks for its scratch directory.

        @return: the `dir` it passed to tempfile.mkdtemp().
        """

        manager = self.manager()
        asked = {}

        def record(dir=None, **kwargs):
            asked['dir'] = dir
            raise self.Stop()

        queue_in = mock.Mock()
        queue_in.get.return_value = os.path.join(
            self.genome('GCF_000000001.1'), 'prodigal', 'GCF_000000001.1_protein.faa.gz')
        worker = getattr(manager, '_MarkerManager' + worker_name)

        with mock.patch.object(M.tempfile, 'mkdtemp', record):
            with self.assertRaises(self.Stop):
                worker(queue_in, mock.Mock(), SUFFIX)

        return asked['dir']

    def test_the_pfam_worker_unzips_under_the_tmp_dir_it_was_given(self):
        self.assertEqual(self.scratch_dir_of('__pfam_worker'), self.dir)

    def test_the_tigrfam_worker_unzips_under_the_tmp_dir_it_was_given(self):
        self.assertEqual(self.scratch_dir_of('__tigrfam_worker'), self.dir)

    def test_a_tmp_dir_that_is_not_there_is_made_before_any_batch_is_claimed(self):
        # met once, at the start, rather than once per genome inside a batch this
        # machine has already taken and would then fail
        scratch = os.path.join(self.dir, 'scratch', 'deeper')

        self.manager(tmp_dir=scratch)

        self.assertTrue(os.path.isdir(scratch))


class ChoosingTheMarkerDatabase(TempDirCase):
    def test_pfam_and_tigrfam_name_their_own_directories_and_files(self):
        pfam = self.manager().marker_setup('pfam', SUFFIX, self.hmm_db('pfam'))
        tigr = self.manager().marker_setup('tigrfam', TIGR_SUFFIX,
                                           self.hmm_db('tigrfam'))

        self.assertEqual((pfam.marker_dir, pfam.extension, pfam.name),
                         ('pfam_33.1_lite', '_pfam_33.1_lite.tsv', 'Pfam'))
        self.assertEqual((tigr.marker_dir, tigr.extension, tigr.name),
                         ('tigrfam_15.0_lite', '_tigrfam_15.0_lite.tsv', 'Tigrfam'))

    def test_an_unknown_database_is_refused_rather_than_failing_later(self):
        # it used to leave the marker directory unbound and fail further in with
        # a NameError naming nothing
        with self.assertRaises(ValueError):
            self.manager().marker_setup('panther', SUFFIX, self.hmm_db('pfam'))


class TheHmmDatabaseIsCheckedBeforeTheSearch(TempDirCase):
    """--hmm_db_path, which means a directory for one --db and a file for the other.

    Pfam is handed the directory Pfam-A.hmm sits in and TIGRFAM the HMM file
    itself, and the r237 run gave TIGRFAM the directory. hmmsearch read it as a
    file that "appears to be empty" and wrote no marker table; the run met that as
    a FileNotFoundError in a worker, one call later, having already claimed and
    finished 135 batches.
    """

    def test_a_pfam_directory_holding_the_library_is_taken(self):
        M.check_hmm_db('pfam', self.hmm_db('pfam'))

    def test_a_tigrfam_hmm_file_is_taken(self):
        M.check_hmm_db('tigrfam', self.hmm_db('tigrfam'))

    def test_the_directory_pfam_wants_is_refused_for_tigrfam(self):
        directory = os.path.dirname(self.hmm_db('tigrfam'))

        with self.assertRaises(M.BadHmmDatabase) as raised:
            M.check_hmm_db('tigrfam', directory)

        # and says which file in it was meant, that being the whole of the mistake
        self.assertIn('tigrfam.hmm', str(raised.exception))

    def test_the_file_tigrfam_wants_is_refused_for_pfam(self):
        with self.assertRaises(M.BadHmmDatabase):
            M.check_hmm_db('pfam', self.hmm_db('tigrfam'))

    def test_a_directory_with_no_library_in_it_is_refused_for_pfam(self):
        empty = os.path.join(self.dir, 'empty')
        os.makedirs(empty)

        with self.assertRaises(M.BadHmmDatabase) as raised:
            M.check_hmm_db('pfam', empty)

        self.assertIn(M.PFAM_LIBRARY, str(raised.exception))

    def test_a_path_that_is_not_there_is_refused(self):
        with self.assertRaises(M.BadHmmDatabase):
            M.check_hmm_db('tigrfam', os.path.join(self.dir, 'nothing.hmm'))

        with self.assertRaises(M.BadHmmDatabase):
            M.check_hmm_db('pfam', os.path.join(self.dir, 'nothing'))

    def test_a_file_that_is_not_an_hmm_library_is_refused(self):
        # the format line is read rather than the name: a FASTA called tigrfam.hmm
        # is refused, and hmmsearch would have said the same thing a batch later
        not_hmms = self.hmm_db('tigrfam', magic=b'>gene\nMAKV\n')

        with self.assertRaises(M.BadHmmDatabase):
            M.check_hmm_db('tigrfam', not_hmms)

    def test_the_run_claims_no_batch_when_the_hmms_are_wrong(self):
        # the point of checking in marker_setup(): nothing is planned, claimed or
        # finished, so no other machine is told the release was searched
        genomes = {'GCF_000000001.1': self.genome('GCF_000000001.1')}
        dirs_file, report = self.inputs(genomes)
        directory = os.path.dirname(self.hmm_db('tigrfam'))

        with self.assertRaises(M.BadHmmDatabase):
            self.manager().run_hmmsearch(dirs_file, report, 'tigrfam', TIGR_SUFFIX,
                                         directory, self.out_dir)

        self.assertFalse(os.path.isdir(self.state_dir('tigrfam_' + TIGR_SUFFIX)))


class WhenTheSearchProcessesDie(TempDirCase):
    """A batch is finished on what came back, not on what was handed out.

    search_markers() started the workers, joined them and returned, whatever they
    had done: BatchCounts.searched is the length of the work list. So a batch in
    which every worker died on its first genome was written a SUCCESS canary
    saying searched=10,000, and no later run would look behind it.
    """

    def test_a_worker_that_raised_fails_the_search(self):
        with self.assertRaises(RuntimeError) as raised:
            self.manager().search_markers(['/no/such/genome_protein.faa.gz'],
                                          dying_worker, SUFFIX)

        self.assertIn('1 of 1 search process(es) ended in error',
                      str(raised.exception))

    def test_a_worker_that_exited_fails_the_search(self):
        # PfamScan gives up with sys.exit() rather than an exception, and that is
        # a worker gone just the same
        with self.assertRaises(RuntimeError):
            self.manager().search_markers(['/no/such/genome_protein.faa.gz'],
                                          exiting_worker, SUFFIX)

    def test_workers_that_did_the_work_are_not_called_a_failure(self):
        self.manager(cpus=2).search_markers(
            ['/one_protein.faa.gz', '/two_protein.faa.gz'], quiet_worker, SUFFIX)

    def test_a_search_with_nothing_to_do_is_not_called_a_failure(self):
        self.manager().search_markers([], quiet_worker, SUFFIX)

    def test_the_batch_is_failed_rather_than_finished(self):
        genomes = {'GCF_000000001.1': self.genome('GCF_000000001.1')}

        def search(manager, genome_files, worker, dir_suffix):
            M.MarkerManager.search_markers(manager, genome_files, dying_worker,
                                           dir_suffix)

        finished = self.run_hmmsearch(genomes, searcher=search)

        self.assertFalse(finished)
        batch = self.batches()[0]
        self.assertEqual(B.batch_state(batch), B.STATE_FAILED)
        self.assertFalse(os.path.exists(os.path.join(batch, B.SUCCESS_CANARY)))


class TheVersionThatSearchedAGenome(TempDirCase):
    """hmmsearch.version goes into the marker directory with the marker table,
    since the table outlives the run that made it."""

    def search_one(self, returncode=0):
        """Run the TIGRFAM worker over one genome, hmmsearch replaced by a stub
        that writes what the real one writes.

        @return: the genome's marker directory.
        """

        manager = self.manager()
        manager.tigrfam_hmms = self.hmm_db('tigrfam')
        manager.hmmer_version = VERSION
        gene_file = os.path.join(self.genome('GCF_000000001.1'), 'prodigal',
                                 'GCF_000000001.1_protein.faa.gz')

        def hmmsearch(cmd, **kwargs):
            if returncode == 0:
                for flag in ('-o', '--tblout'):
                    with open(cmd[cmd.index(flag) + 1], 'w') as handle:
                        handle.write('# nothing found\n')
            return subprocess.CompletedProcess(cmd, returncode, stderr='refused')

        queue_in = mock.Mock()
        queue_in.get.side_effect = [gene_file, None]

        with mock.patch.object(M.subprocess, 'run', hmmsearch):
            try:
                manager._MarkerManager__tigrfam_worker(queue_in, mock.Mock(),
                                                       TIGR_SUFFIX)
            except RuntimeError:
                pass

        return os.path.join(os.path.dirname(gene_file),
                            'tigrfam_{}'.format(TIGR_SUFFIX))

    def test_a_searched_genome_records_the_version_that_searched_it(self):
        marker_dir = self.search_one()

        with open(os.path.join(marker_dir, 'hmmsearch.version')) as handle:
            self.assertEqual(handle.read(), VERSION + '\n')

    def test_a_search_that_failed_records_no_version(self):
        marker_dir = self.search_one(returncode=1)

        self.assertFalse(os.path.exists(os.path.join(marker_dir, 'hmmsearch.version')))

    def test_both_databases_ask_hmmsearch_which_is_what_both_run(self):
        """PfamScan runs hmmsearch too, whatever its name suggests."""
        for db in ('pfam', 'tigrfam'):
            manager = self.manager()
            with mock.patch.object(M, 'record_program_version',
                                   return_value=VERSION) as asked:
                manager.marker_setup(db, SUFFIX, self.hmm_db(db))
            asked.assert_called_once_with('hmmsearch')
            self.assertEqual(manager.hmmer_version, VERSION)


class WhatTheTigrfamWorkerDoesWithAFailedSearch(TempDirCase):
    """os.system() returned the exit status to nobody.

    A search that failed was met one call later, in _tigr_top_hit(), as a
    FileNotFoundError on the marker table that was never written -- a traceback
    naming the missing file and not the reason there was none.
    """

    def worker_on(self, completed):
        """Run the TIGRFAM worker over one genome with hmmsearch replaced.

        @return: nothing; what the worker raised is what is being tested.
        """

        manager = self.manager()
        manager.tigrfam_hmms = self.hmm_db('tigrfam')
        gene_file = os.path.join(self.genome('GCF_000000001.1'), 'prodigal',
                                 'GCF_000000001.1_protein.faa.gz')

        queue_in = mock.Mock()
        queue_in.get.return_value = gene_file

        with mock.patch.object(M.subprocess, 'run', return_value=completed):
            manager._MarkerManager__tigrfam_worker(queue_in, mock.Mock(),
                                                   TIGR_SUFFIX)

    def test_a_search_that_failed_raises_what_hmmsearch_said(self):
        failed = subprocess.CompletedProcess(
            [], 1, stderr='Error: File format problem in trying to open HMM file '
                          '/srv/db/gtdb/marker_genes/hmms. File exists, but appears '
                          'to be empty?\n')

        with self.assertRaises(RuntimeError) as raised:
            self.worker_on(failed)

        message = str(raised.exception)
        self.assertIn('exited 1', message)
        self.assertIn('File format problem', message)
        # and not the FileNotFoundError on the table the search never wrote
        self.assertNotIn('_tigrfam_15.0_lite.tsv', message)

    def test_a_search_that_said_nothing_still_raises(self):
        with self.assertRaises(RuntimeError):
            self.worker_on(subprocess.CompletedProcess([], 1, stderr=''))

    def test_a_search_that_ran_is_not_raised_on(self):
        # it gets as far as reading the table, which a search that ran would have
        # written; the point is that it is reached at all
        with self.assertRaises(FileNotFoundError):
            self.worker_on(subprocess.CompletedProcess([], 0, stderr=''))


if __name__ == '__main__':
    unittest.main()
