#!/usr/bin/env python3
"""Offline unit tests for lpsn.py -- no LPSN website, no HTML.

Run with the interpreter that has pandas, bs4, requests and sqlalchemy:

    /opt/centos7/sw/miniconda3/envs/python3/bin/python -m unittest -v tests.test_lpsn

What is tested here is the bookkeeping of pull_html: a page whose name
duplicates a validly published one is skipped, recorded in <rank>_skipped.lst
rather than the failed list, and that file exists only when there is something
in it. Nothing is downloaded: the skip is decided before any request is made.
"""

import io
import os
import shutil
import tempfile
import unittest

from gtdb_migration_tk.lpsn import LPSN


class TempDirCase(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='lpsn_test.')
        self.lpsn = LPSN(skip_taxa_per_letter_dl=True, lpsn_output_dir=self.dir)

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def skipped_file(self, rank):
        return os.path.join(self.dir, '{}_skipped.lst'.format(rank))


class SkippedFile(TempDirCase):
    def test_written_only_when_something_was_skipped(self):
        self.assertIsNone(self.lpsn.write_skipped('genus', []))
        self.assertFalse(os.path.exists(self.skipped_file('genus')))

    def test_lists_each_skipped_page_with_its_name_and_reason(self):
        skipped = [('"Alcaligenes faecalis subsp. faecalis"',
                    'https://lpsn.dsmz.de/subspecies/alcaligenes-faecalis-faecalis-1')]

        path = self.lpsn.write_skipped('subspecies', skipped)

        self.assertEqual(path, self.skipped_file('subspecies'))
        with open(path) as handle:
            self.assertEqual(handle.read(),
                             '"Alcaligenes faecalis subsp. faecalis"\t'
                             'https://lpsn.dsmz.de/subspecies/alcaligenes-faecalis-faecalis-1\t'
                             'duplicate_name\n')

    def test_a_file_from_an_earlier_run_is_removed_when_nothing_is_skipped_now(self):
        # the file's presence is the signal, so a stale one would misreport this run
        with open(self.skipped_file('genus'), 'w') as handle:
            handle.write('old\n')

        self.lpsn.write_skipped('genus', [])

        self.assertFalse(os.path.exists(self.skipped_file('genus')))


class DuplicateNames(TempDirCase):
    def test_a_duplicate_name_is_skipped_not_failed(self):
        failed = io.StringIO()
        skipped = []
        valid_names = ['Alcaligenes faecalis subsp. faecalis']   # the unquoted twin

        n = self.lpsn.download_rank_name(
            '"Alcaligenes faecalis subsp. faecalis"',
            'https://lpsn.dsmz.de/subspecies/alcaligenes-faecalis-faecalis-1',
            os.path.join(self.dir, 'alcaligenes-faecalis-faecalis-1'),
            valid_names, failed, skipped, num_already_dl=0)

        self.assertEqual(skipped, [('"Alcaligenes faecalis subsp. faecalis"',
                                    'https://lpsn.dsmz.de/subspecies/alcaligenes-faecalis-faecalis-1')])
        self.assertEqual(failed.getvalue(), '')
        self.assertEqual(n, 0)
        # and nothing was fetched or written for it
        self.assertEqual(os.listdir(self.dir), [])

    def test_an_already_downloaded_valid_page_is_counted_not_refetched(self):
        out_file = os.path.join(self.dir, 'alcaligenes-faecalis-faecalis')
        open(out_file, 'w').close()
        failed = io.StringIO()
        skipped = []

        n = self.lpsn.download_rank_name(
            'Alcaligenes faecalis subsp. faecalis',
            'https://lpsn.dsmz.de/subspecies/alcaligenes-faecalis-faecalis',
            out_file, ['Alcaligenes faecalis subsp. faecalis'], failed, skipped,
            num_already_dl=0)

        self.assertEqual(n, 1)
        self.assertEqual(skipped, [])
        self.assertEqual(failed.getvalue(), '')


if __name__ == '__main__':
    unittest.main()
