#!/usr/bin/env python3
"""Offline unit tests for the vendored Prodigal wrapper -- Prodigal is never run.

What is tested here is how the wrapper FILES a genome's results. The gene calling
is Prodigal's business; leaving a half-written protein file in a genome directory
of the release, where nothing sweeps it up and the next command reads it, is the
wrapper's.
"""

import gzip
import hashlib
import os
import shutil
import tempfile
import unittest
from unittest import mock

from gtdb_migration_tk.biolib_lite.external import prodigal as W


class TempDirCase(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='prodigal_wrapper_test.')

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def source(self, text='>gene\nMAAA\n'):
        path = os.path.join(self.dir, 'genes.faa')
        with open(path, 'w') as handle:
            handle.write(text)
        return path

    def destination(self, name='GCF_1.1_protein.faa.gz'):
        genome = os.path.join(self.dir, 'genome', 'prodigal')
        os.makedirs(genome, exist_ok=True)
        return os.path.join(genome, name)


class CompressToTests(TempDirCase):
    """The results are the release, not scratch: a reader sees a whole file."""

    def test_the_bytes_arrive_gzipped_and_unchanged(self):
        source, destination = self.source(), self.destination()
        W.compress_to(source, destination)
        with gzip.open(destination, 'rt') as handle:
            self.assertEqual(handle.read(), '>gene\nMAAA\n')

    def test_the_digest_is_of_the_uncompressed_bytes(self):
        """It is what vouches for the genes, not for the gzip container, which
        recompressing or copying about can change on its own."""
        source, destination = self.source(), self.destination()
        digest = W.compress_to(source, destination)
        self.assertEqual(digest, hashlib.sha1(b'>gene\nMAAA\n').hexdigest())

    def test_the_destination_is_never_opened_for_writing(self):
        """Two machines can meet on one genome -- --reclaim takes a claim
        deliberately -- and writing into the destination left them interleaving
        into one file and a corrupt gzip behind. Nothing opens it: it is moved
        onto."""
        source, destination = self.source(), self.destination()
        opened = []

        real_open = gzip.open

        def watching(path, *args, **kwargs):
            opened.append(path)
            return real_open(path, *args, **kwargs)

        with mock.patch.object(W.gzip, 'open', watching):
            W.compress_to(source, destination)

        self.assertNotIn(destination, opened)
        self.assertEqual(len(opened), 1)
        self.assertTrue(os.path.exists(destination))

    def test_a_write_that_fails_leaves_no_file_behind(self):
        """A half-written file beside a genome is swept up by nothing: this
        directory is the release and not a scratch directory. The failure is put
        where the staged file already exists, which is the only place the
        cleanup is reached from."""
        source, destination = self.source(), self.destination()

        # inside the copy, once the gzip is open and bytes have been written:
        # the old code had the destination itself half written by now
        broken = mock.Mock()
        broken.update.side_effect = RuntimeError('disk full')
        with mock.patch.object(W.hashlib, 'sha1', return_value=broken):
            self.assertRaises(RuntimeError, W.compress_to, source, destination)

        self.assertFalse(os.path.exists(destination))
        self.assertEqual(os.listdir(os.path.dirname(destination)), [])

    def test_a_move_that_fails_leaves_no_file_behind(self):
        """The last step, where the staged file is complete and the destination
        is not yet there."""
        source, destination = self.source(), self.destination()

        with mock.patch.object(W.os, 'replace',
                               side_effect=OSError('read-only filesystem')):
            self.assertRaises(OSError, W.compress_to, source, destination)

        self.assertFalse(os.path.exists(destination))
        self.assertEqual(os.listdir(os.path.dirname(destination)), [])

    def test_a_write_that_fails_leaves_an_earlier_result_untouched(self):
        """A genome that had proteins keeps them: better the ones that are
        vouched for than nothing at all."""
        source, destination = self.source(), self.destination()
        W.compress_to(source, destination)

        with mock.patch.object(W.os, 'replace',
                               side_effect=OSError('read-only filesystem')):
            self.assertRaises(OSError, W.compress_to,
                              self.source('>gene\nMBBB\n'), destination)

        with gzip.open(destination, 'rt') as handle:
            self.assertEqual(handle.read(), '>gene\nMAAA\n')

    def test_a_second_write_replaces_the_first(self):
        """--all_genomes calls a genome again over results that are already
        there."""
        destination = self.destination()
        W.compress_to(self.source('>gene\nMAAA\n'), destination)
        W.compress_to(self.source('>gene\nMCCC\n'), destination)
        with gzip.open(destination, 'rt') as handle:
            self.assertEqual(handle.read(), '>gene\nMCCC\n')

    def test_the_staged_file_is_beside_the_destination(self):
        """So the move is a rename within one filesystem and not a copy of
        every protein file of a release."""
        source, destination = self.source(), self.destination()
        staged = []

        real_mkstemp = tempfile.mkstemp

        def watching(*args, **kwargs):
            handle, path = real_mkstemp(*args, **kwargs)
            staged.append(path)
            return handle, path

        with mock.patch.object(W.tempfile, 'mkstemp', watching):
            W.compress_to(source, destination)

        self.assertEqual(len(staged), 1)
        self.assertEqual(os.path.dirname(staged[0]),
                         os.path.dirname(destination))


if __name__ == '__main__':
    unittest.main()
