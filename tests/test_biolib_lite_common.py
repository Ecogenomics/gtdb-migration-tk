#!/usr/bin/env python3
"""Offline unit tests for biolib_lite/common.py.

get_num_lines() sizes the progress bar of most commands that walk a genome_dirs
file, so a file it cannot count ends the command before the command can say
anything about the file itself.
"""

import os
import shutil
import stat
import tempfile
import unittest

from gtdb_migration_tk.biolib_lite.common import get_num_lines


class CountingLines(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='biolib_lite_common_test.')
        self.addCleanup(shutil.rmtree, self.dir, True)

    def file_of(self, text, mode=None):
        path = os.path.join(self.dir, 'genome_dirs.tsv')
        with open(path, 'w') as handle:
            handle.write(text)
        if mode is not None:
            os.chmod(path, mode)
            self.addCleanup(os.chmod, path, stat.S_IRUSR | stat.S_IWUSR)
        return path

    def test_lines_are_counted(self):
        self.assertEqual(get_num_lines(self.file_of('a\tb\nc\td\n')), 2)

    def test_a_last_line_with_no_newline_is_counted(self):
        self.assertEqual(get_num_lines(self.file_of('a\tb\nc\td')), 2)

    def test_an_empty_file_has_no_lines(self):
        # mmap refuses a file of no bytes, and the ValueError reached the caller
        self.assertEqual(get_num_lines(self.file_of('')), 0)

    def test_a_read_only_file_is_counted(self):
        # it was opened for writing, which a read-only release file refuses
        self.assertEqual(get_num_lines(self.file_of('a\tb\n', mode=stat.S_IRUSR)), 1)


if __name__ == '__main__':
    unittest.main()
