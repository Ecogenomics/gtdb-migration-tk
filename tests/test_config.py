#!/usr/bin/env python3
"""Offline unit tests for config.py.

A marker version left behind in one place produces genome directories whose
symlinks point at annotation files that were never written, and nothing fails
loudly when it happens. These tests assert that every name in config.py derives
from the two version constants, so bumping a version cannot leave part of the
vocabulary on the old release.
"""

import unittest

from gtdb_migration_tk import config


class DerivedNames(unittest.TestCase):
    def test_every_name_carries_the_configured_version(self):
        for name, expected in (('PFAM_MARKER_DIR', config.PFAM_VERSION),
                               ('PFAM_EXT', config.PFAM_VERSION),
                               ('PFAM_TOPHIT_EXT', config.PFAM_VERSION),
                               ('TIGRFAM_MARKER_DIR', config.TIGRFAM_VERSION),
                               ('TIGRFAM_EXT', config.TIGRFAM_VERSION),
                               ('TIGRFAM_TOPHIT_EXT', config.TIGRFAM_VERSION),
                               ('TIGRFAM_OUT_EXT', config.TIGRFAM_VERSION)):
            self.assertIn(expected, getattr(config, name), name)

    def test_symlink_names_carry_no_version(self):
        # downstream code opens a genome's hits without knowing the release,
        # so only the symlink target may change when a version is bumped
        for name in ('PFAM_SYMLINK_EXT', 'PFAM_TOPHIT_SYMLINK_EXT',
                     'TIGRFAM_SYMLINK_EXT', 'TIGRFAM_TOPHIT_SYMLINK_EXT',
                     'TIGRFAM_OUT_SYMLINK_EXT'):
            value = getattr(config, name)
            self.assertNotIn(config.PFAM_VERSION, value, name)
            self.assertNotIn(config.TIGRFAM_VERSION, value, name)

    def test_hmmer_extensions_track_both_versions(self):
        joined = ''.join(config.HMMER_EXTS_TO_GZIP)
        self.assertIn(config.PFAM_VERSION, joined)
        self.assertIn(config.TIGRFAM_VERSION, joined)

    def test_a_version_bump_reaches_every_derived_name(self):
        # execute config.py with the two version literals rewritten: every name
        # that follows must change with them, or it repeats a literal instead of
        # deriving from the constant, which is the mistake this guards against
        source = open(config.__file__).read()
        source = source.replace("PFAM_VERSION = '33.1'", "PFAM_VERSION = '37.0'")
        source = source.replace("TIGRFAM_VERSION = '15.0'", "TIGRFAM_VERSION = '16.0'")

        bumped = {}
        exec(compile(source, config.__file__, 'exec'), bumped)

        self.assertEqual(bumped['PFAM_MARKER_DIR'], 'pfam_37.0_lite')
        self.assertEqual(bumped['PFAM_EXT'], '_pfam_37.0_lite.tsv.gz')
        self.assertEqual(bumped['PFAM_TOPHIT_EXT'], '_pfam_37.0_lite_tophit.tsv.gz')
        self.assertEqual(bumped['TIGRFAM_MARKER_DIR'], 'tigrfam_16.0_lite')
        self.assertEqual(bumped['TIGRFAM_OUT_EXT'], '_tigrfam_16.0_lite.out.gz')
        self.assertEqual(bumped['HMMER_EXTS_TO_GZIP'],
                         ('_pfam_37.0.tsv', '_pfam_37.0_tophit.tsv',
                          '_tigrfam_16.0.out', '_tigrfam_16.0.tsv',
                          '_tigrfam_16.0_tophit.tsv'))

        # no name may still mention the superseded releases
        for name, value in bumped.items():
            if name.startswith('_'):
                continue
            text = ''.join(value) if isinstance(value, tuple) else str(value)
            self.assertNotIn('33.1', text, name)
            self.assertNotIn('15.0', text, name)


class CurrentValues(unittest.TestCase):
    """The names in use for the current release, as a guard against typos."""

    def test_names_match_the_directories_on_disk(self):
        self.assertEqual(config.PFAM_MARKER_DIR, 'pfam_33.1_lite')
        self.assertEqual(config.TIGRFAM_MARKER_DIR, 'tigrfam_15.0_lite')
        self.assertEqual(config.HMMER_EXTS_TO_GZIP,
                         ('_pfam_33.1.tsv', '_pfam_33.1_tophit.tsv',
                          '_tigrfam_15.0.out', '_tigrfam_15.0.tsv',
                          '_tigrfam_15.0_tophit.tsv'))


if __name__ == '__main__':
    unittest.main()
