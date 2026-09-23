#!/usr/bin/env python3
"""Offline unit tests for biolib_lite/logger.py -- where a run's log goes.

logger_setup() is called once by __main__ before anything is dispatched, and a
second time where the first could not open the log it was given. Both loggers are
named, so the second call is handed the SAME objects, and everything here is about
that: a run whose --log named a file where a directory was wanted printed itself
double from beginning to end, and said nothing about the log it had not opened.
"""

import io
import logging
import os
import shutil
import sys
import tempfile
import unittest

from gtdb_migration_tk.__main__ import FALLBACK_LOG, log_candidates
from gtdb_migration_tk.biolib_lite.logger import logger_setup


class LoggerSetupCase(unittest.TestCase):
    """Both loggers are global, so what a test does to them is put back."""

    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='logger_test.')
        self.saved = {name: list(logging.getLogger(name).handlers)
                      for name in ('timestamp', 'no_timestamp')}
        for name in self.saved:
            logging.getLogger(name).handlers = []

        self.stdout, sys.stdout = sys.stdout, io.StringIO()

    def tearDown(self):
        for logger in (logging.getLogger('timestamp'),
                       logging.getLogger('no_timestamp')):
            for handler in list(logger.handlers):
                logger.removeHandler(handler)
                handler.close()

        for name, handlers in self.saved.items():
            logging.getLogger(name).handlers = handlers

        sys.stdout = self.stdout
        shutil.rmtree(self.dir, ignore_errors=True)

    def setup(self, log_dir=None, log_file='run.log', silent=False):
        logger_setup(self.dir if log_dir is None else log_dir, log_file,
                     'GTDB Migration Tk', 'gtdb_migration_tk', '0.0.0', silent)

    def console(self):
        return sys.stdout.getvalue()

    def log_file(self, name='run.log'):
        with open(os.path.join(self.dir, name)) as handle:
            return handle.read()


class WhatASecondCallDoes(LoggerSetupCase):
    """__main__ calls logger_setup() again where the first call could not open
    the log. The handlers of the first were left on the logger, and every line of
    the run reached the console once per handler."""

    def test_a_second_call_does_not_double_the_console(self):
        self.setup()
        self.setup(log_dir=os.path.join(self.dir, 'second'), log_file='other.log')

        logging.getLogger('timestamp').info('one line')

        self.assertEqual(self.console().count('one line'), 1)

    def test_a_second_call_leaves_one_handler_of_each_kind(self):
        self.setup()
        self.setup()

        # a console handler and a file handler, on each of the two loggers
        for name in ('timestamp', 'no_timestamp'):
            self.assertEqual(len(logging.getLogger(name).handlers), 2)

    def test_the_log_of_the_first_call_stops_being_written_to(self):
        # it is not the log of this run; the second call is what said where that
        # goes, and a line in both would have a reader believe either
        self.setup()
        self.setup(log_file='second.log')

        logging.getLogger('timestamp').info('after the second call')

        self.assertNotIn('after the second call', self.log_file('run.log'))
        self.assertIn('after the second call', self.log_file('second.log'))

    def test_the_second_calls_log_holds_each_line_once(self):
        self.setup()
        self.setup(log_file='second.log')

        logging.getLogger('timestamp').info('one line')

        self.assertEqual(self.log_file('second.log').count('one line'), 1)

    def test_the_no_timestamp_logger_is_not_doubled_either(self):
        self.setup()
        self.setup()

        logging.getLogger('no_timestamp').info('plain line')

        self.assertEqual(self.console().count('plain line'), 1)

    def test_a_silent_second_call_does_not_leave_a_talking_handler(self):
        # --silent lifts the console handler to ERROR; the handler of a first,
        # not silent call would have gone on printing everything
        self.setup()
        self.setup(silent=True)

        logging.getLogger('timestamp').info('should not be seen')

        self.assertNotIn('should not be seen', self.console())


class WhenTheLogCannotBeOpened(LoggerSetupCase):
    """__main__ falls back to ./gtdb_migration_tk.log, and it can only do that if
    it is told. The failure has to come out of here as an exception."""

    def test_a_log_directory_that_is_a_file_raises(self):
        # -l logs/run.log where ./logs is a file, which is one keystroke from
        # -l logs and what the r237 patch run met
        not_a_dir = os.path.join(self.dir, 'logs')
        with open(not_a_dir, 'w') as handle:
            handle.write('written by an earlier run\n')

        with self.assertRaises(OSError):
            self.setup(log_dir=not_a_dir)

    def test_the_run_can_still_be_logged_afterwards(self):
        not_a_dir = os.path.join(self.dir, 'logs')
        with open(not_a_dir, 'w') as handle:
            handle.write('written by an earlier run\n')
        with self.assertRaises(OSError):
            self.setup(log_dir=not_a_dir)

        # what __main__ does next, and the line it writes about it
        self.setup(log_file='gtdb_migration_tk.log')
        logging.getLogger('timestamp').warning('logged here instead')

        self.assertEqual(self.console().count('logged here instead'), 1)
        self.assertIn('logged here instead',
                      self.log_file('gtdb_migration_tk.log'))

    def test_no_log_directory_leaves_the_console_alone(self):
        # a falsy directory is "no log file", which is what ncbi_genome_sync run
        # as a script relies on
        self.setup(log_dir='')

        logging.getLogger('timestamp').info('console only')

        self.assertIn('console only', self.console())
        self.assertFalse(os.path.exists(os.path.join(self.dir, 'run.log')))


class WhereARunIsLogged(unittest.TestCase):
    """The order __main__ tries, which decides where a log nobody asked for goes.

    A run whose --log could not be opened was logged to the directory the command
    happened to be run from. The rest of what a run produces goes under its
    --out_dir, and that is where someone looks for its log.
    """

    def test_the_log_that_was_asked_for_comes_first(self):
        self.assertEqual(log_candidates('logs/run.log', 'out')[0],
                         ('logs', 'run.log'))

    def test_the_output_directory_is_tried_before_the_current_one(self):
        candidates = log_candidates('logs/run.log', 'out')

        self.assertEqual(candidates[1], ('out', FALLBACK_LOG))
        self.assertEqual(candidates[2], ('.', FALLBACK_LOG))

    def test_a_command_with_no_log_flag_is_logged_under_its_output_directory(self):
        self.assertEqual(log_candidates(None, 'out')[0], ('out', FALLBACK_LOG))

    def test_a_command_with_no_output_directory_falls_to_the_current_one(self):
        self.assertEqual(log_candidates('run.log', None),
                         [('.', 'run.log'), ('.', FALLBACK_LOG),
                          (None, FALLBACK_LOG)])

    def test_a_bare_filename_is_logged_to_the_current_directory(self):
        # dirname('sync.log') is '', which logger_setup() reads as "no log file"
        self.assertEqual(log_candidates('sync.log', None)[0], ('.', 'sync.log'))

    def test_the_console_is_the_last_resort_and_cannot_fail(self):
        for candidates in (log_candidates('logs/run.log', 'out'),
                           log_candidates(None, None)):
            self.assertEqual(candidates[-1], (None, FALLBACK_LOG))

    def test_the_same_place_is_not_tried_twice(self):
        # -l gtdb_migration_tk.log with --out_dir . names one file three ways
        self.assertEqual(log_candidates(FALLBACK_LOG, '.'),
                         [('.', FALLBACK_LOG), (None, FALLBACK_LOG)])


if __name__ == '__main__':
    unittest.main()
