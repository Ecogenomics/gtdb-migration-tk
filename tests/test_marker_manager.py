#!/usr/bin/env python3
"""Offline unit tests for marker_manager.py -- no HMM is ever searched.

Run with the interpreter that has tqdm:

    /opt/centos7/sw/miniconda3/envs/python3/bin/python -m unittest -v tests.test_marker_manager

What is tested here is everything hmmsearch decides before an HMM is searched:
how the release is cut into batches and shared between machines, which genomes of
a batch still need their markers found, and what is handed to the workers.
search_markers() is replaced by a recorder, so the real run_hmmsearch() runs over
real batch directories and what it would have searched can be read back.

The batching machinery itself is tested in tests/test_batching.py. What is tested
here is this command's use of it: that a Pfam run and a TIGRFAM run of one output
directory do not read each other's canaries, and that the work list can hold
nothing a worker would mistake for the end of the work.
"""

import gzip
import os
import shutil
import tempfile
import unittest
from unittest import mock

from gtdb_migration_tk import batching as B
from gtdb_migration_tk import marker_manager as M


SUFFIX = '33.1_lite'
MARKER_DIR = 'pfam_' + SUFFIX
MARKER_EXT = '_pfam_{}.tsv'.format(SUFFIX)
TIGR_SUFFIX = '15.0_lite'


class TempDirCase(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='marker_manager_test.')
        self.out_dir = os.path.join(self.dir, 'out')
        os.makedirs(self.out_dir)
        self.searched = []

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def manager(self, cpus=1, batch_size=10000, **kwargs):
        # the real one checks for prodigal and hmmsearch on PATH and exits without them
        with mock.patch.object(M, 'check_dependencies'):
            return M.MarkerManager(tmp_dir=self.dir, cpus=cpus,
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
                                         '/nonexistent/hmms', self.out_dir, **kwargs)

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


class ChoosingTheMarkerDatabase(TempDirCase):
    def test_pfam_and_tigrfam_name_their_own_directories_and_files(self):
        pfam = self.manager().marker_setup('pfam', SUFFIX, '/hmms')
        tigr = self.manager().marker_setup('tigrfam', TIGR_SUFFIX, '/hmms')

        self.assertEqual((pfam.marker_dir, pfam.extension, pfam.name),
                         ('pfam_33.1_lite', '_pfam_33.1_lite.tsv', 'Pfam'))
        self.assertEqual((tigr.marker_dir, tigr.extension, tigr.name),
                         ('tigrfam_15.0_lite', '_tigrfam_15.0_lite.tsv', 'Tigrfam'))

    def test_an_unknown_database_is_refused_rather_than_failing_later(self):
        # it used to leave the marker directory unbound and fail further in with
        # a NameError naming nothing
        with self.assertRaises(ValueError):
            self.manager().marker_setup('panther', SUFFIX, '/hmms')


if __name__ == '__main__':
    unittest.main()
