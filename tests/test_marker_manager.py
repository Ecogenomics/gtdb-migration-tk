#!/usr/bin/env python3
"""Offline unit tests for marker_manager.py -- no HMM is ever searched.

Run with the interpreter that has tqdm:

    /opt/centos7/sw/miniconda3/envs/python3/bin/python -m unittest -v tests.test_marker_manager

What is tested here is what run_hmmsearch() decides before any HMM is searched:
which genomes of a release still need their markers found, and what reaches the
queue the workers draw from. The multiprocessing is replaced by stubs that run
the real marker_parser() serially and record what was queued, because the bug
this file exists for is not in either worker -- it is in what the work list is
allowed to contain.
"""

import gzip
import os
import shutil
import tempfile
import types
import unittest
from unittest import mock

from gtdb_migration_tk import marker_manager as M


SUFFIX = '33.1_lite'
MARKER_DIR = 'pfam_' + SUFFIX
MARKER_EXT = '_pfam_{}.tsv'.format(SUFFIX)


class RecordingQueue:
    """Records what was put on it, in order, rather than crossing a process."""

    def __init__(self):
        self.items = []

    def put(self, item):
        self.items.append(item)

    def get(self, block=True, timeout=None):
        return self.items.pop(0)


class SerialPool:
    """Runs marker_parser() in this process, so the decisions under test are real."""

    def __init__(self, processes=1):
        self.processes = processes

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def imap_unordered(self, func, iterable):
        return [func(item) for item in iterable]


class NoopProcess:
    """A worker that is never started: the queue it would drain is the evidence."""

    started = []

    def __init__(self, target=None, args=()):
        self.target = target
        self.args = args

    def start(self):
        NoopProcess.started.append(self)

    def join(self):
        pass

    def terminate(self):
        pass


class TempDirCase(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='marker_manager_test.')
        self.queues = []
        NoopProcess.started = []

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def manager(self, cpus=1):
        # the real one checks for prodigal and hmmsearch on PATH and exits without them
        with mock.patch.object(M, 'check_dependencies'):
            return M.MarkerManager(tmp_dir=self.dir, cpus=cpus)

    def genome(self, gid, protein=b'>gene\nMA\n', annotated=False, checksum=True):
        """A genome directory, optionally already carrying its Pfam results.

        Parameters
        ----------
        gid : str
            Accession, which names the files within.
        protein : bytes or None
            Contents of the protein FASTA, gzipped; None writes no protein file
            at all, and b'' writes a zero byte one, which is what the code means
            by empty -- it stats the file on disk, and a gzip of nothing is still
            twenty-odd bytes of header.
        annotated : bool
            Whether a marker table is already there.
        checksum : bool
            Whether that table's .sha256 is there and correct.

        @return: the genome directory.
        """

        path = os.path.join(self.dir, gid)
        prodigal_dir = os.path.join(path, 'prodigal')
        os.makedirs(prodigal_dir)

        if protein == b'':
            open(os.path.join(prodigal_dir, gid + '_protein.faa.gz'), 'wb').close()
        elif protein is not None:
            with gzip.open(os.path.join(prodigal_dir, gid + '_protein.faa.gz'), 'wb') as handle:
                handle.write(protein)

        if annotated:
            marker_dir = os.path.join(prodigal_dir, MARKER_DIR)
            os.makedirs(marker_dir)
            table = os.path.join(marker_dir, gid + MARKER_EXT)
            with gzip.open(table + '.gz', 'wb') as handle:
                handle.write(b'# hits\n')
            if checksum:
                # of the UNCOMPRESSED bytes, as marker_parser() reads it back
                with open(table + '.gz', 'rb') as raw:
                    digest = M.sha256_rb(gzip.GzipFile(fileobj=raw))
                with open(table + '.sha256', 'w') as handle:
                    handle.write(digest + '\n')

        return path

    def inputs(self, genomes, outcomes=None):
        """The two files run_hmmsearch() reads.

        Parameters
        ----------
        genomes : dict
            Accession to genome directory.
        outcomes : dict
            Accession to its report outcome; the default is 'new' for every one.

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

    def run_hmmsearch(self, genomes, outcomes=None, cpus=1):
        """Run the real run_hmmsearch() with the multiprocessing stubbed out.

        @return: (what was put on the work queue, what the progress bar was told).
        """

        dirs_file, report = self.inputs(genomes, outcomes)

        def make_queue():
            self.queues.append(RecordingQueue())
            return self.queues[-1]

        stub = types.SimpleNamespace(Queue=make_queue, Pool=SerialPool, Process=NoopProcess)
        with mock.patch.object(M, 'mp', stub):
            self.manager(cpus=cpus).run_hmmsearch(
                dirs_file, report, 'pfam', SUFFIX, '/nonexistent/hmms')

        worker_queue = self.queues[0]
        progress = [p for p in NoopProcess.started if p.target.__name__.endswith('__progress')]
        return worker_queue.items, progress[0].args[0]


class WhatReachesTheWorkers(TempDirCase):
    """The work queue ends with one None per worker, and holds none anywhere else.

    Every worker breaks on the first None it draws, so a None among the genomes is
    a second stop signal: the worker that takes it exits with the rest of the
    release still queued, and does so silently. marker_parser() returns None for
    every genome it skips, which is why nothing it returns can be queued untested.
    """

    def test_a_genome_needing_markers_is_queued(self):
        genomes = {'GCF_000000001.1': self.genome('GCF_000000001.1')}

        queued, _ = self.run_hmmsearch(genomes)

        self.assertEqual(queued[:-1], [os.path.join(
            genomes['GCF_000000001.1'], 'prodigal', 'GCF_000000001.1_protein.faa.gz')])

    def test_the_only_none_is_the_one_that_stops_each_worker(self):
        genomes = {'GCF_000000001.1': self.genome('GCF_000000001.1')}

        queued, _ = self.run_hmmsearch(genomes, cpus=3)

        self.assertEqual(queued[-3:], [None, None, None])
        self.assertNotIn(None, queued[:-3])

    def test_a_genome_with_no_protein_file_does_not_stop_a_worker(self):
        # it cannot be searched, so it is skipped -- but skipped by being left out
        # of the queue, not by being put on it as the value that ends the queue
        genomes = {'GCF_000000001.1': self.genome('GCF_000000001.1'),
                   'GCF_000000002.1': self.genome('GCF_000000002.1', protein=None)}

        queued, _ = self.run_hmmsearch(genomes)

        self.assertEqual(len(queued), 2)              # one genome, one sentinel
        self.assertIsNone(queued[-1])
        self.assertNotIn(None, queued[:-1])

    def test_a_genome_with_an_empty_protein_file_does_not_stop_a_worker(self):
        genomes = {'GCF_000000001.1': self.genome('GCF_000000001.1'),
                   'GCF_000000002.1': self.genome('GCF_000000002.1', protein=b'')}

        queued, _ = self.run_hmmsearch(genomes)

        self.assertEqual(len(queued), 2)
        self.assertNotIn(None, queued[:-1])

    def test_a_genome_already_annotated_does_not_stop_a_worker_either(self):
        # the other reason marker_parser() skips a genome; it used to be told apart
        # from the two above, which is how those two came to be queued
        genomes = {'GCF_000000001.1': self.genome('GCF_000000001.1'),
                   'GCF_000000002.1': self.genome('GCF_000000002.1', annotated=True)}

        queued, _ = self.run_hmmsearch(genomes)

        self.assertEqual(len(queued), 2)
        self.assertNotIn(None, queued[:-1])

    def test_every_genome_skipped_leaves_only_the_sentinels(self):
        # nothing to search is a run that ends, not a run that hangs or crashes
        genomes = {'GCF_000000001.1': self.genome('GCF_000000001.1', protein=None),
                   'GCF_000000002.1': self.genome('GCF_000000002.1', annotated=True)}

        queued, _ = self.run_hmmsearch(genomes, cpus=2)

        self.assertEqual(queued, [None, None])

    def test_the_progress_bar_counts_the_genomes_that_will_be_searched(self):
        # the denominator used to include the skipped genomes, so a run that did
        # everything asked of it still reported less than 100%
        genomes = {'GCF_000000001.1': self.genome('GCF_000000001.1'),
                   'GCF_000000002.1': self.genome('GCF_000000002.1', protein=None),
                   'GCF_000000003.1': self.genome('GCF_000000003.1', annotated=True)}

        queued, denominator = self.run_hmmsearch(genomes)

        self.assertEqual(denominator, 1)
        self.assertEqual(denominator, len(queued) - 1)


class WhichGenomesMarkerParserSkips(TempDirCase):
    """marker_parser() decides one genome, and says so with one value.

    'null' used to mean 'already annotated' while a genome with no protein file
    fell off the end of the function as None. run_hmmsearch() dropped the first
    and queued the second.
    """

    def parse(self, gid, consider=True, **kwargs):
        path = self.genome(gid, **kwargs)
        job = (gid, path, MARKER_DIR, MARKER_EXT, {gid} if consider else set(), 'Pfam')
        return self.manager().marker_parser(job)

    def test_a_genome_to_annotate_gives_its_protein_file(self):
        result = self.parse('GCF_000000001.1')

        self.assertTrue(result.endswith('GCF_000000001.1_protein.faa.gz'), result)
        self.assertTrue(os.path.exists(result))

    def test_an_annotated_genome_with_a_matching_checksum_is_skipped(self):
        self.assertIsNone(self.parse('GCF_000000001.1', annotated=True))

    def test_a_genome_with_no_protein_file_is_skipped(self):
        self.assertIsNone(self.parse('GCF_000000001.1', protein=None))

    def test_a_genome_with_an_empty_protein_file_is_skipped(self):
        self.assertIsNone(self.parse('GCF_000000001.1', protein=b''))

    def test_an_annotated_genome_whose_checksum_does_not_match_is_annotated_again(self):
        # a table that does not match its own checksum is not results anyone can use
        result = self.parse('GCF_000000001.1', annotated=True, checksum=False)

        self.assertIsNotNone(result)
        self.assertTrue(result.endswith('_protein.faa.gz'), result)

    def test_skipping_is_one_value_so_nothing_can_be_dropped_by_halves(self):
        # the contract run_hmmsearch()'s filter rests on: there is no second skip
        # value for it to miss
        skipped = [self.parse('GCF_00000000%d.1' % n, **case)
                   for n, case in enumerate((dict(annotated=True),
                                             dict(protein=None),
                                             dict(protein=b'')), start=1)]

        self.assertEqual(skipped, [None, None, None])


if __name__ == '__main__':
    unittest.main()
