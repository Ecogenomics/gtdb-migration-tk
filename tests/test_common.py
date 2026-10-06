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

import argparse
import json
import logging
import os
import shutil
import stat
import tempfile
import unittest
from unittest import mock

import psycopg2

from gtdb_migration_tk.database_configuration.GenomeDatabaseConnectionFTPUpdate import \
    GenomeDatabaseConnectionFTPUpdate
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

    def conda_env(self, packages):
        """A conda environment holding a program, with packages in conda-meta."""
        prefix = os.path.join(self.dir, 'env')
        os.makedirs(os.path.join(prefix, 'bin'))
        os.makedirs(os.path.join(prefix, 'conda-meta'))
        executable = os.path.join(prefix, 'bin', 'pplacer.exe')
        open(executable, 'w').close()
        os.symlink('pplacer.exe', os.path.join(prefix, 'bin', 'pplacer'))
        for name, version in packages:
            with open(os.path.join(prefix, 'conda-meta',
                                   '{}-{}-h0_0.json'.format(name, version)), 'w') as handle:
                json.dump({'name': name, 'version': version}, handle)
        return os.path.join(prefix, 'bin', 'pplacer')

    def test_a_conda_package_version_is_read_from_its_environments_conda_meta(self):
        executable = self.conda_env([('pplacer', '1.1.alpha22')])

        self.assertEqual(C.conda_package_version(executable, 'pplacer'), '1.1.alpha22')

    def test_a_package_whose_name_begins_with_anothers_is_not_taken_for_it(self):
        executable = self.conda_env([('pplacer-extra', '9.9'), ('pplacer', '1.1.alpha20')])

        self.assertEqual(C.conda_package_version(executable, 'pplacer'), '1.1.alpha20')

    def test_an_executable_outside_a_conda_environment_has_no_package_version(self):
        executable = self.conda_env([])
        shutil.rmtree(os.path.join(self.dir, 'env', 'conda-meta'))

        self.assertIsNone(C.conda_package_version(executable, 'pplacer'))

    def test_an_environment_without_the_package_has_no_version_of_it(self):
        executable = self.conda_env([('guppy', '1.0')])

        self.assertIsNone(C.conda_package_version(executable, 'pplacer'))

    def test_the_version_file_holds_the_version_on_one_line(self):
        C.write_version_file(self.dir, 'tRNAscan-SE', 'tRNAscan-SE 2.0.13 (Jul 2026)')

        with open(os.path.join(self.dir, 'trnascan-se.version')) as handle:
            self.assertEqual(handle.read(), 'tRNAscan-SE 2.0.13 (Jul 2026)\n')


def database_options(**given):
    """The options __database_setup() gives a command, as argparse sets them."""
    options = dict(db_service=None, hostname=None, user=None, db=None, password=None)
    options.update(given)
    return argparse.Namespace(**options)


class NamingTheDatabase(unittest.TestCase):
    def setUp(self):
        patcher = mock.patch.dict(os.environ)
        patcher.start()
        self.addCleanup(patcher.stop)
        os.environ.pop('PGSERVICE', None)

    def test_a_service_alone_names_the_database(self):
        self.assertEqual(C.database_keywords(database_options(db_service='gtdb_r237')),
                         {'service': 'gtdb_r237'})

    def test_host_user_and_database_name_it_with_no_password(self):
        # the password is then read from ~/.pgpass
        keywords = C.database_keywords(database_options(hostname='h', user='u', db='d'))
        self.assertEqual(keywords, {'host': 'h', 'user': 'u', 'dbname': 'd'})

    def test_options_given_with_a_service_are_handed_over_beside_it(self):
        # libpq takes a keyword given over the service's value for it
        keywords = C.database_keywords(database_options(db_service='gtdb_r237', db='gtdb_r237_test'))
        self.assertEqual(keywords, {'service': 'gtdb_r237', 'dbname': 'gtdb_r237_test'})

    def test_pgservice_in_the_environment_names_the_database(self):
        os.environ['PGSERVICE'] = 'gtdb_r237'
        self.assertEqual(C.database_keywords(database_options()), {})

    def test_nothing_given_is_refused(self):
        # libpq's defaults are a socket on this machine and a user named for
        # whoever runs the command, which is no database the toolkit is used on
        with self.assertRaises(C.DatabaseOptionsError):
            C.database_keywords(database_options())

    def test_a_host_without_a_user_and_database_is_refused(self):
        with self.assertRaises(C.DatabaseOptionsError):
            C.database_keywords(database_options(hostname='h', password='pw'))


class HandingTheKeywordsOver(unittest.TestCase):
    PASSWORD = 'a pass@word/with \'quotes\''

    def connect_with(self, opener):
        """The keywords psycopg2.connect() is called with, nothing being reached."""
        with mock.patch.object(psycopg2, 'connect', side_effect=ConnectionRefusedError) as connect:
            with self.assertRaises(ConnectionRefusedError):
                opener()
        return connect.call_args.kwargs

    def test_the_connection_wrapper_hands_a_password_over_unchanged(self):
        # it was written into a connection string unquoted, where a space ended it
        database = {'host': 'h', 'user': 'u', 'dbname': 'd', 'password': self.PASSWORD}
        kwargs = self.connect_with(GenomeDatabaseConnectionFTPUpdate(database).MakePostgresConnection)
        self.assertEqual(kwargs, database)

    def test_the_sqlalchemy_engine_hands_psycopg2_the_keywords(self):
        # it was given a URL, which held no service, fixed the port at 5432, and
        # read a password holding '@' or '/' as part of the host
        database = {'service': 'gtdb_r237', 'password': self.PASSWORD}
        kwargs = self.connect_with(C.database_engine(database).connect)
        self.assertEqual(kwargs, database)


if __name__ == '__main__':
    unittest.main()
