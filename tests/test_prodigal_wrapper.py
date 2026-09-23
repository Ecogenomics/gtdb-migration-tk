#!/usr/bin/env python3
"""Offline unit tests for the vendored Prodigal wrapper -- Prodigal is never run.

What is tested here is how the wrapper FILES a genome's results. The gene calling
is Prodigal's business; leaving a half-written protein file in a genome directory
of the release, where nothing sweeps it up and the next command reads it, is the
wrapper's.

The meta mode fallback is tested against a Prodigal of a few lines of shell, put
on PATH for the test: what matters is what the wrapper does with an exit status
and with a mode, and a stub that answers on its argv tests that where a mocked
subprocess would only test the mock.
"""

import gzip
import hashlib
import os
import shutil
import stat
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


class MetaFallbackTests(TempDirCase):
    """Single mode trains on the genome; meta mode answers for what it refuses.

    Eleven genomes of one 2014 submission of N-rich actinomycetes, 9 to 13 Mb
    each, were called by neither mode in r237: Prodigal exits 52 on them saying
    "saw too many regions of N's", and the empty protein file it leaves behind was
    carried forward from release to release. Meta mode calls eight to twelve
    thousand genes for each of them.
    """

    def prodigal(self, single_exit=0, meta_exit=0, complaint="too many N's",
                 read_stdin=True):
        """A Prodigal of a few lines of shell, first on PATH for this test.

        It records its argv for each run, answers on the mode it was given, and
        writes a GFF to stdout where it succeeds.

        Parameters
        ----------
        single_exit : int
            What it exits with under -p single.
        meta_exit : int
            What it exits with under -p meta.
        complaint : str
            What it says on stderr where it exits non-zero.
        read_stdin : bool
            Whether it reads the genome at all. A Prodigal that refuses one
            without reading it closes stdin while the wrapper is still writing.

        @return: the file its argv are recorded in, one run per line.
        """

        argv_log = os.path.join(self.dir, 'argv.log')
        script = os.path.join(self.dir, 'bin', 'prodigal')
        os.makedirs(os.path.dirname(script), exist_ok=True)
        with open(script, 'w') as handle:
            handle.write(
                '#!/bin/bash\n'
                'echo "$@" >> {log}\n'
                'pwd -P > {cwd}\n'
                'mode=single; for a in "$@"; do [ "$a" = meta ] && mode=meta; done\n'
                '{read_in}\n'
                'if [ $mode = single ]; then code={single}; else code={meta}; fi\n'
                'if [ $code -ne 0 ]; then echo "Error: {complaint}" >&2; exit $code; fi\n'
                'echo "##gff-version 3"\n'
                'for f in "$@"; do case $prev in -a|-d) : > "$f";; esac; prev=$f; done\n'
                'exit 0\n'.format(
                    log=argv_log, cwd=os.path.join(self.dir, 'cwd'),
                    single=single_exit, meta=meta_exit,
                    complaint=complaint,
                    read_in='cat > /dev/null' if read_stdin else 'true'))
        os.chmod(script, os.stat(script).st_mode | stat.S_IEXEC)

        self.addCleanup(os.environ.__setitem__, 'PATH', os.environ['PATH'])
        os.environ['PATH'] = os.path.dirname(script) + os.pathsep + os.environ['PATH']
        return argv_log

    def genome(self, bases=b'ACGT' * 1000):
        path = os.path.join(self.dir, 'genome.fna.gz')
        with gzip.open(path, 'wb') as handle:
            handle.write(b'>contig\n' + bases + b'\n')
        return path

    def call(self, table=11, mode='single'):
        """Run the wrapper's one attempt-and-fall-back over the stub.

        @return: what run_prodigal() returned.
        """

        cmd = ['prodigal', '-m', '-p', mode, '-q', '-f', 'gff', '-g', str(table),
               '-a', os.path.join(self.dir, 'genes.faa'),
               '-d', os.path.join(self.dir, 'genes.fna')]
        return W.run_prodigal(cmd, self.genome(),
                              os.path.join(self.dir, 'genes.gff'), self.dir)

    def runs(self, argv_log):
        with open(argv_log) as handle:
            return handle.read().splitlines()

    def test_a_genome_called_in_single_mode_reports_no_fallback(self):
        argv_log = self.prodigal()

        self.assertIsNone(self.call())
        self.assertEqual(len(self.runs(argv_log)), 1)

    def test_a_genome_single_mode_refuses_is_called_in_meta_mode(self):
        argv_log = self.prodigal(single_exit=52)

        complaint = self.call()

        self.assertIn("too many N's", complaint)
        runs = self.runs(argv_log)
        self.assertEqual(len(runs), 2)
        self.assertIn('-p single', runs[0])
        self.assertIn('-p meta', runs[1])

    def test_the_fallback_keeps_the_translation_table(self):
        # the table is what the genome is called under and is not in question:
        # gTranslate predicted it from the genome, and the mode is what failed
        argv_log = self.prodigal(single_exit=52)

        self.call(table=4)

        for run in self.runs(argv_log):
            self.assertIn('-g 4', run)

    def test_a_genome_neither_mode_calls_raises_with_what_both_said(self):
        self.prodigal(single_exit=52, meta_exit=55)

        with self.assertRaises(RuntimeError) as raised:
            self.call()

        message = str(raised.exception)
        self.assertIn('single mode with exit code 52', message)
        self.assertIn('meta mode with exit code 55', message)

    def test_a_genome_already_in_meta_mode_is_not_tried_again(self):
        # a genome too small to train on is called in meta mode from the start,
        # and there is nothing to fall back to
        argv_log = self.prodigal(meta_exit=55)

        with self.assertRaises(RuntimeError) as raised:
            self.call(mode='meta')

        self.assertEqual(len(self.runs(argv_log)), 1)
        self.assertIn('exit code 55', str(raised.exception))

    def test_prodigal_closing_stdin_early_is_not_an_error_of_its_own(self):
        """A refusal that does not read the genome leaves the wrapper writing into
        a closed pipe. BrokenPipeError is not RuntimeError, so it would come out
        of the pool and fail the batch rather than the genome."""
        self.prodigal(single_exit=52, meta_exit=55, read_stdin=False)

        with self.assertRaises(RuntimeError):
            self.call()

    def test_a_genome_the_fallback_saves_is_not_raised_on(self):
        self.prodigal(single_exit=52, meta_exit=0, read_stdin=False)

        self.assertIn("too many N's", self.call())

    def test_prodigal_runs_in_the_scratch_directory(self):
        """It spools stdin into tmp.prodigal.stdin.<pid> in the CURRENT directory
        and leaves it there when it exits non-zero: a copy of the genome,
        uncompressed, per failure, in the working directory of whoever started the
        run. Here it goes when the genome's scratch directory does."""
        self.prodigal(single_exit=52)
        here = os.getcwd()

        self.call()

        with open(os.path.join(self.dir, 'cwd')) as handle:
            self.assertEqual(handle.read().strip(), os.path.realpath(self.dir))
        self.assertEqual(os.getcwd(), here)


if __name__ == '__main__':
    unittest.main()
