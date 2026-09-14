#!/usr/bin/env python3
"""Offline unit tests for config.py.

A database version left behind in one place makes one command write a directory
another command does not read, and nothing fails loudly when it happens. These
tests assert that every derived name in config.py follows from the version
constants, so bumping a version cannot leave part of the vocabulary on the old
release.
"""

import unittest

from gtdb_migration_tk import config


class DerivedNames(unittest.TestCase):
    def test_marker_folder_suffixes_carry_the_configured_versions(self):
        # hmmsearch and top_hit default to these, so they decide which
        # pfam_*/tigrfam_* directory a release is annotated into
        self.assertIn(config.PFAM_VERSION, config.MARKER_FOLDER_SUFFIX['pfam'])
        self.assertIn(config.TIGRFAM_VERSION, config.MARKER_FOLDER_SUFFIX['tigrfam'])

    def test_marker_folder_suffixes_name_the_lite_sets(self):
        # GTDB searches the reduced marker sets; marker_manager keys the
        # version-free symlink names off this suffix
        for suffix in config.MARKER_FOLDER_SUFFIX.values():
            self.assertTrue(suffix.endswith('_lite'), suffix)

    def test_marker_folder_suffixes_are_keyed_by_the_db_option(self):
        self.assertEqual(sorted(config.MARKER_FOLDER_SUFFIX), ['pfam', 'tigrfam'])

    def test_derived_directories_carry_the_configured_rrna_versions(self):
        # a genome directory may hold the results of several SILVA or LTP
        # releases side by side, so the wrong version here copies the wrong one
        self.assertIn('rna_silva_{}'.format(config.SILVA_VERSION),
                      config.GTDB_DERIVED_DIRS_TO_COPY)
        self.assertIn('rna_ltp_{}'.format(config.LTP_VERSION),
                      config.GTDB_DERIVED_DIRS_TO_COPY)

    def test_unversioned_derived_directories_are_listed(self):
        # prodigal and trna results are not named for a database release, but
        # are carried across with the rest and must not be dropped from the list
        for name in ('prodigal', 'trna'):
            self.assertIn(name, config.GTDB_DERIVED_DIRS_TO_COPY)

    def test_a_version_bump_reaches_every_derived_name(self):
        # execute config.py with the version literals rewritten: every name that
        # follows must change with them, or it repeats a literal instead of
        # deriving from the constant, which is the mistake this guards against
        source = open(config.__file__).read()
        source = source.replace("PFAM_VERSION = '33.1'", "PFAM_VERSION = '37.0'")
        source = source.replace("TIGRFAM_VERSION = '15.0'", "TIGRFAM_VERSION = '16.0'")
        source = source.replace("SILVA_VERSION = '138.2'", "SILVA_VERSION = '140.0'")
        source = source.replace("LTP_VERSION = '10_2024'", "LTP_VERSION = '06_2026'")

        bumped = {}
        exec(compile(source, config.__file__, 'exec'), bumped)

        self.assertEqual(bumped['MARKER_FOLDER_SUFFIX'],
                         {'pfam': '37.0_lite', 'tigrfam': '16.0_lite'})
        self.assertEqual(bumped['GTDB_DERIVED_DIRS_TO_COPY'],
                         ('prodigal', 'rna_silva_140.0', 'trna', 'rna_ltp_06_2026'))

        # no name may still mention the superseded releases
        for name, value in bumped.items():
            if name.startswith('_'):
                continue
            if isinstance(value, dict):
                text = ''.join(value.values())
            elif isinstance(value, tuple):
                text = ''.join(value)
            else:
                text = str(value)
            for superseded in ('33.1', '15.0', '138.2', '10_2024'):
                self.assertNotIn(superseded, text, name)


class CurrentValues(unittest.TestCase):
    """The names in use for the current release, as a guard against typos."""

    def test_names_match_the_directories_on_disk(self):
        self.assertEqual(config.MARKER_FOLDER_SUFFIX,
                         {'pfam': '33.1_lite', 'tigrfam': '15.0_lite'})
        self.assertEqual(config.GTDB_DERIVED_DIRS_TO_COPY,
                         ('prodigal', 'rna_silva_138.2', 'trna', 'rna_ltp_10_2024'))


if __name__ == '__main__':
    unittest.main()
