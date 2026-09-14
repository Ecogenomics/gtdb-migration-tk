#!/usr/bin/env python3
"""Offline unit tests for ncbi_ftp_manager_tools.py -- no FTP site, no real genomes.

Run with the interpreter that has tqdm:

    /opt/centos7/sw/miniconda3/envs/python3/bin/python -m unittest -v tests.test_ncbi_ftp_manager_tools

What is tested here is how a genome held by both the previous release and the
NCBI mirror is carried into the new release: which side each file comes from,
and what the report says about it. The genome directories are a handful of
small files standing in for the real ones; only their names and their entries
in md5checksums.txt matter to the code under test.
"""

import io
import os
import shutil
import tempfile
import unittest

from gtdb_migration_tk import config
from gtdb_migration_tk import ncbi_ftp_manager_tools as T


ACCESSION = 'GCF_000000001.1'
ASSEMBLY = ACCESSION + '_ASM1v1'
FASTA = ASSEMBLY + '_genomic.fna.gz'
REPORT = ASSEMBLY + '_assembly_report.txt'

# MD5s as NCBI would publish them; the FASTA contents are never hashed, so
# any 32 hex digits will do
MD5_A = 'a' * 32
MD5_B = 'b' * 32
MD5_R = 'c' * 32


def manifest(fasta_md5, *extra):
    """An md5checksums.txt naming the FASTA, the report, and any extra lines."""
    lines = ['{}  ./{}'.format(fasta_md5, FASTA),
             '{}  ./{}'.format(MD5_R, REPORT)] + list(extra)
    return '\n'.join(lines) + '\n'


class TempDirCase(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='ncbi_ftp_manager_tools_test.')
        self.report = io.StringIO()
        self.review = io.StringIO()

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def tools(self, dry_run=False):
        return T.FTPTools(self.report, self.review, dry_run)

    def genome_dir(self, name, fasta_md5, derived=()):
        """A genome directory holding NCBI's files and, optionally, derived data."""
        path = os.path.join(self.dir, name, ASSEMBLY)
        os.makedirs(path)
        for filename, text in ((FASTA, 'fasta'), (REPORT, 'report'),
                               (T.MD5_MANIFEST, manifest(fasta_md5))):
            with open(os.path.join(path, filename), 'w') as handle:
                handle.write(text)
        for subdir in derived:
            os.makedirs(os.path.join(path, subdir))
            with open(os.path.join(path, subdir, ACCESSION + '_derived.tsv'), 'w') as handle:
                handle.write(subdir)
        return path

    def target(self):
        return os.path.join(self.dir, 'new_release', ASSEMBLY)

    def status(self, row):
        accession, status = row.rstrip('\n').split('\t')
        self.assertEqual(accession, ACCESSION)
        return status


class ComparisonOutcome(TempDirCase):
    def test_unchanged_fasta_carries_the_derived_data_across(self):
        prev = self.genome_dir('previous', MD5_A, derived=config.GTDB_DERIVED_DIRS_TO_COPY)
        ftp = self.genome_dir('mirror', MD5_A)

        row = self.tools().compare_genome_directories(prev, ftp, self.target(), ACCESSION)

        self.assertEqual(self.status(row), T.STATUS_FASTA_UNCHANGED)
        for derived in config.GTDB_DERIVED_DIRS_TO_COPY:
            copied = os.path.join(self.target(), derived, ACCESSION + '_derived.tsv')
            self.assertTrue(os.path.isfile(copied), derived)
        self.assertEqual(self.review.getvalue(), '')

    def test_changed_fasta_leaves_the_derived_data_behind(self):
        prev = self.genome_dir('previous', MD5_A, derived=config.GTDB_DERIVED_DIRS_TO_COPY)
        ftp = self.genome_dir('mirror', MD5_B)

        row = self.tools().compare_genome_directories(prev, ftp, self.target(), ACCESSION)

        self.assertEqual(self.status(row), T.STATUS_FASTA_CHANGED)
        for derived in config.GTDB_DERIVED_DIRS_TO_COPY:
            self.assertFalse(os.path.exists(os.path.join(self.target(), derived)), derived)

    def test_every_ncbi_file_is_taken_from_the_mirror_whether_or_not_the_fasta_changed(self):
        # the manifest included: the new release must describe the files it
        # actually holds, and those are the mirror's
        for name, md5 in (('same', MD5_A), ('different', MD5_B)):
            with self.subTest(name):
                prev = self.genome_dir('previous_' + name, MD5_A)
                ftp = self.genome_dir('mirror_' + name, md5)
                target = os.path.join(self.dir, 'release_' + name, ASSEMBLY)

                self.tools().compare_genome_directories(prev, ftp, target, ACCESSION)

                self.assertEqual(sorted(os.listdir(target)),
                                 sorted([FASTA, REPORT, T.MD5_MANIFEST]))
                with open(os.path.join(target, T.MD5_MANIFEST)) as handle:
                    self.assertEqual(handle.read(), manifest(md5))

    def test_derived_data_of_the_previous_release_is_never_modified(self):
        prev = self.genome_dir('previous', MD5_A, derived=config.GTDB_DERIVED_DIRS_TO_COPY)
        ftp = self.genome_dir('mirror', MD5_A)
        before = sorted(os.walk(prev))

        self.tools().compare_genome_directories(prev, ftp, self.target(), ACCESSION)

        self.assertEqual(sorted(os.walk(prev)), before)

    def test_a_missing_derived_directory_is_reported_for_review_and_the_rest_still_copied(self):
        # a genome may not have had every step run on it in the previous release
        present = config.GTDB_DERIVED_DIRS_TO_COPY[:-1]
        missing = config.GTDB_DERIVED_DIRS_TO_COPY[-1]
        prev = self.genome_dir('previous', MD5_A, derived=present)
        ftp = self.genome_dir('mirror', MD5_A)

        row = self.tools().compare_genome_directories(prev, ftp, self.target(), ACCESSION)

        self.assertEqual(self.status(row), T.STATUS_FASTA_UNCHANGED)
        for derived in present:
            self.assertTrue(os.path.isdir(os.path.join(self.target(), derived)), derived)
        self.assertFalse(os.path.exists(os.path.join(self.target(), missing)))
        review = self.review.getvalue()
        self.assertIn(ACCESSION, review)
        self.assertIn(missing, review)

    def test_a_rerun_replaces_a_target_left_by_an_earlier_run(self):
        # derived data copied when the FASTA was unchanged must not survive a
        # rerun made after the FASTA changed
        prev = self.genome_dir('previous', MD5_A, derived=config.GTDB_DERIVED_DIRS_TO_COPY)
        ftp = self.genome_dir('mirror', MD5_A)
        self.tools().compare_genome_directories(prev, ftp, self.target(), ACCESSION)
        self.assertTrue(os.path.isdir(os.path.join(self.target(), 'prodigal')))

        with open(os.path.join(ftp, T.MD5_MANIFEST), 'w') as handle:
            handle.write(manifest(MD5_B))
        row = self.tools().compare_genome_directories(prev, ftp, self.target(), ACCESSION)

        self.assertEqual(self.status(row), T.STATUS_FASTA_CHANGED)
        self.assertFalse(os.path.exists(os.path.join(self.target(), 'prodigal')))


class DryRun(TempDirCase):
    def test_dry_run_reports_the_real_outcome_but_copies_nothing(self):
        # the report of a dry run must be the report the real run would write
        for name, md5, expected in (('same', MD5_A, T.STATUS_FASTA_UNCHANGED),
                                    ('different', MD5_B, T.STATUS_FASTA_CHANGED)):
            with self.subTest(name):
                prev = self.genome_dir('previous_' + name, MD5_A,
                                       derived=config.GTDB_DERIVED_DIRS_TO_COPY)
                ftp = self.genome_dir('mirror_' + name, md5)
                target = os.path.join(self.dir, 'release_' + name, ASSEMBLY)

                row = self.tools(dry_run=True).compare_genome_directories(prev, ftp, target, ACCESSION)

                self.assertEqual(self.status(row), expected)
                self.assertFalse(os.path.exists(target))

    def test_dry_run_still_reports_a_missing_derived_directory(self):
        present = config.GTDB_DERIVED_DIRS_TO_COPY[:-1]
        missing = config.GTDB_DERIVED_DIRS_TO_COPY[-1]
        prev = self.genome_dir('previous', MD5_A, derived=present)
        ftp = self.genome_dir('mirror', MD5_A)

        self.tools(dry_run=True).compare_genome_directories(prev, ftp, self.target(), ACCESSION)

        self.assertIn(missing, self.review.getvalue())
        self.assertFalse(os.path.exists(self.target()))

    def test_dry_run_leaves_an_existing_target_untouched(self):
        # a real run replaces the target; a dry run must not even do that
        prev = self.genome_dir('previous', MD5_A)
        ftp = self.genome_dir('mirror', MD5_B)
        os.makedirs(self.target())
        stale = os.path.join(self.target(), 'stale.txt')
        open(stale, 'w').close()

        self.tools(dry_run=True).compare_genome_directories(prev, ftp, self.target(), ACCESSION)

        self.assertEqual(os.listdir(self.target()), ['stale.txt'])


class ReadingTheManifest(TempDirCase):
    def write_manifest(self, text, name=ASSEMBLY):
        path = os.path.join(self.dir, name)
        os.makedirs(path)
        with open(os.path.join(path, T.MD5_MANIFEST), 'w') as handle:
            handle.write(text)
        return path

    def test_the_genomic_fasta_md5_is_read_from_the_manifest(self):
        path = self.write_manifest(manifest(MD5_A))
        self.assertEqual(self.tools().genomic_fasta_md5(path), MD5_A)

    def test_the_entry_is_found_by_the_name_the_directory_gives_it(self):
        # the name is built from the directory, not searched for by suffix, so
        # the derived FASTA files NCBI also lists cannot be picked up instead
        path = self.write_manifest(manifest(
            MD5_A,
            '{}  ./{}_cds_from_genomic.fna.gz'.format(MD5_B, ASSEMBLY),
            '{}  ./{}_rna_from_genomic.fna.gz'.format(MD5_B, ASSEMBLY)))
        self.assertEqual(self.tools().genomic_fasta_md5(path), MD5_A)

    def test_a_manifest_without_the_genomic_fasta_is_an_error(self):
        path = self.write_manifest('{}  ./{}\n'.format(MD5_R, REPORT))
        with self.assertRaises(ValueError):
            self.tools().genomic_fasta_md5(path)

    def test_a_directory_named_for_another_assembly_is_an_error(self):
        # the manifest describes ASSEMBLY, but the directory says otherwise
        path = self.write_manifest(manifest(MD5_A), name='GCF_000000002.1_ASM2v1')
        with self.assertRaises(ValueError):
            self.tools().genomic_fasta_md5(path)

    def test_a_missing_manifest_is_an_error(self):
        with self.assertRaises(OSError):
            self.tools().genomic_fasta_md5(os.path.join(self.dir, 'no_such_genome'))


if __name__ == '__main__':
    unittest.main()
