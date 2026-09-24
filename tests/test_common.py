#!/usr/bin/env python3
"""Offline unit tests for utils/common.py -- how an external program's version
is asked for and recorded.

Every program the toolkit runs is recorded in the run's log and in a
<program>.version file beside what it made, because the results outlive the
run: a genome's tRNAs, proteins and marker tables are carried across from
release to release while its sequences are unchanged, and only the file can
then say what made them. The programs are stood in for by a few lines of shell
put first on PATH, so what is tested is what the helper does with what a
program prints and how it exits, not any one build of it.
"""

import logging
import os
import shutil
import stat
import tempfile
import unittest

from gtdb_migration_tk.utils import common as C


# What each program printed when asked, on the machine building r237. A pattern
# that stops matching the program it was written for would fail every run of
# that command before it started, so each one is held against the real thing.
REAL_OUTPUT = {
    'tRNAscan-SE': ('\ntRNAscan-SE 2.0.13 (Jul 2026)\nCopyright (C) 2022 Patricia '
                    'Chan and Todd Lowe\n', 'tRNAscan-SE 2.0.13 (Jul 2026)'),
    'prodigal': ('\nProdigal V2.6.3: February, 2016\n\n',
                 'Prodigal V2.6.3: February, 2016'),
    'hmmsearch': ('# hmmsearch :: search profile(s) against a sequence database\n'
                  '# HMMER 3.4 (Aug 2023); http://hmmer.org/\n',
                  'HMMER 3.4 (Aug 2023)'),
    'nhmmer': ('# nhmmer :: search a DNA model, alignment, or sequence against a '
               'DNA database\n# HMMER 3.4 (Aug 2023); http://hmmer.org/\n',
               'HMMER 3.4 (Aug 2023)'),
    'blastn': ('blastn: 2.12.0+\n Package: blast 2.12.0, build Jul 13 2021\n',
               'blastn: 2.12.0+'),
    'makeblastdb': ('makeblastdb: 2.12.0+\n Package: blast 2.12.0, build Jul 13 '
                    '2021\n', 'makeblastdb: 2.12.0+'),
    'gtranslate': ('gtranslate: version 0.0.4 Copyright 2025 Pierre-Alain Chaumeil '
                   'and Donovan Parks\n', 'gtranslate: version 0.0.4'),
    'checkm2': ('1.1.0\n', '1.1.0'),
    'checkm': ('\n                ...::: CheckM v1.0.18 :::...\n\n  Lineage-specific '
               'marker set:\n', 'CheckM v1.0.18'),
    'busco': ('BUSCO 5.4.7\n', 'BUSCO 5.4.7'),
}


class TempDirCase(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='common_test.')
        self.bin = os.path.join(self.dir, 'bin')
        os.makedirs(self.bin)
        self.addCleanup(os.environ.__setitem__, 'PATH', os.environ['PATH'])
        os.environ['PATH'] = self.bin + os.pathsep + os.environ['PATH']

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def program(self, name, stdout='', stderr='', status=0):
        """A program of a few lines of shell, first on PATH for this test."""
        path = os.path.join(self.bin, name)
        with open(path, 'w') as handle:
            handle.write('#!/bin/sh\n')
            handle.write("printf '%s' '{}'\n".format(stdout))
            handle.write("printf '%s' '{}' >&2\n".format(stderr))
            handle.write('exit {}\n'.format(status))
        os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR)


class AskingAProgramItsVersion(TempDirCase):
    def test_a_version_printed_on_stdout_is_read(self):
        self.program('prodigal', stdout=REAL_OUTPUT['prodigal'][0])

        self.assertEqual(C.program_version('prodigal'),
                         'Prodigal V2.6.3: February, 2016')

    def test_a_version_printed_on_stderr_is_read(self):
        """tRNAscan-SE -h prints its usage on stdout and its version on stderr."""
        self.program('tRNAscan-SE', stdout='Usage: tRNAscan-SE [-options]',
                     stderr=REAL_OUTPUT['tRNAscan-SE'][0])

        self.assertEqual(C.program_version('tRNAscan-SE'),
                         'tRNAscan-SE 2.0.13 (Jul 2026)')

    def test_a_version_printed_by_a_program_that_exits_non_zero_is_still_read(self):
        """A version printed is a version, whatever the program made of the
        rest of what it was asked."""
        self.program('tRNAscan-SE', stderr=REAL_OUTPUT['tRNAscan-SE'][0], status=1)

        self.assertEqual(C.program_version('tRNAscan-SE'),
                         'tRNAscan-SE 2.0.13 (Jul 2026)')

    def test_a_program_that_is_not_there_is_an_error_naming_it(self):
        os.environ['PATH'] = self.bin

        with self.assertRaises(RuntimeError) as raised:
            C.program_version('busco')

        self.assertIn('busco --version', str(raised.exception))

    def test_a_program_that_states_no_version_is_an_error_rather_than_a_guess(self):
        """A version file saying something other than what ran is worse than
        none."""
        self.program('busco', stdout='usage: busco [-h] ...')

        with self.assertRaises(RuntimeError) as raised:
            C.program_version('busco')

        self.assertIn('usage: busco', str(raised.exception))

    def test_every_pattern_matches_what_its_program_really_prints(self):
        self.assertEqual(set(REAL_OUTPUT), set(C.VERSION_QUERIES))
        for program, (printed, version) in REAL_OUTPUT.items():
            with self.subTest(program=program):
                self.program(program, stdout=printed)
                self.assertEqual(C.program_version(program), version)


class RecordingTheVersion(TempDirCase):
    def test_the_version_is_said_in_the_log_of_the_run(self):
        self.program('prodigal', stdout=REAL_OUTPUT['prodigal'][0])

        with self.assertLogs('timestamp', logging.INFO) as logged:
            version = C.record_program_version('prodigal')

        self.assertEqual(version, 'Prodigal V2.6.3: February, 2016')
        self.assertIn('prodigal: Prodigal V2.6.3: February, 2016', logged.output[0])

    def test_the_version_file_is_named_for_the_program_in_lower_case(self):
        self.assertEqual(C.version_file('/g/trna', 'tRNAscan-SE'),
                         os.path.join('/g/trna', 'trnascan-se.version'))

    def test_the_version_file_holds_the_version_on_one_line(self):
        C.write_version_file(self.dir, 'tRNAscan-SE', 'tRNAscan-SE 2.0.13 (Jul 2026)')

        with open(os.path.join(self.dir, 'trnascan-se.version')) as handle:
            self.assertEqual(handle.read(), 'tRNAscan-SE 2.0.13 (Jul 2026)\n')


if __name__ == '__main__':
    unittest.main()
