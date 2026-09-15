#!/usr/bin/env python3
"""Offline unit tests for update_genomes.py -- no FTP site, no genome files.

Run with the interpreter that has tqdm:

    /opt/centos7/sw/miniconda3/envs/python3/bin/python -m unittest -v tests.test_update_genomes

What is tested here is the update of a release from the mirror: which genomes
the previous release loses, gains and keeps, and how a genome held by both the
previous release and the NCBI mirror is carried into the new release -- which
side each file comes from, and what the report says about it. The genome
directories are a handful of small files standing in for the real ones; only
their names and their entries in md5checksums.txt matter to the code under test.
"""

import io
import os
import shutil
import tempfile
import unittest

from gtdb_migration_tk import config
from gtdb_migration_tk import ncbi_utils as U
from gtdb_migration_tk import update_genomes as UG


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
        self.dir = tempfile.mkdtemp(prefix='update_genomes_test.')
        self.report = io.StringIO()
        self.review = io.StringIO()

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def write(self, name, text):
        path = os.path.join(self.dir, name)
        with open(path, 'w') as handle:
            handle.write(text)
        return path

    def tools(self, dry_run=False):
        return UG.FTPTools(self.report, self.review, dry_run)

    def genome_dir(self, name, fasta_md5, derived=()):
        """A genome directory holding NCBI's files and, optionally, derived data."""
        path = os.path.join(self.dir, name, ASSEMBLY)
        os.makedirs(path)
        for filename, text in ((FASTA, 'fasta'), (REPORT, 'report'),
                               (UG.MD5_MANIFEST, manifest(fasta_md5))):
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


# ----------------------------------------------------------------- release comparison

class GenomeManagerTests(TempDirCase):
    def setUp(self):
        super().setUp()
        self.manager = UG.GenomeManager(U.REFSEQ_PREFIX, self.dir)
        self.new_genomes = {'G1': '/ftp/g1', 'G2': '/ftp/g2', 'G3': '/ftp/g3_new'}
        self.old_genomes = {'G2': '/gtdb/g2', 'G3': '/gtdb/g3_old', 'G4': '/gtdb/g4'}

    def test_genomes_to_remove_are_those_gone_from_ncbi(self):
        self.assertEqual(self.manager.generate_genomes_to_remove(self.new_genomes, self.old_genomes),
                         {'G4': '/gtdb/g4'})

    def test_genomes_to_add_are_those_new_to_ncbi(self):
        self.assertEqual(self.manager.generate_genomes_to_add(self.new_genomes, self.old_genomes),
                         {'G1': '/ftp/g1'})

    def test_genomes_to_compare_are_those_in_both_releases(self):
        # G3 is in both, but at a different path, so it must be compared and not ignored
        self.assertEqual(sorted(self.manager.generate_genomes_to_compare(self.new_genomes, self.old_genomes)),
                         ['G2', 'G3'])

    def test_empty_releases_are_handled(self):
        self.assertEqual(self.manager.generate_genomes_to_remove({}, {}), {})
        self.assertEqual(self.manager.generate_genomes_to_add({}, {}), {})
        self.assertEqual(self.manager.generate_genomes_to_compare({}, {}), [])

    def test_load_genome_dirs_keeps_only_the_database_of_the_manager(self):
        # one genome_dirs file describes a whole release, both databases together
        path = self.write('dirs.tsv',
                          'GCF_000000001.1\t/ftp/g1\n'
                          'GCF_000000002.1\t/ftp/g2\n'
                          'GCA_000000003.1\t/ftp/g3\n')
        self.assertEqual(self.manager.load_genome_dirs(path),
                         {'GCF_000000001.1': '/ftp/g1', 'GCF_000000002.1': '/ftp/g2'})
        self.assertEqual(UG.GenomeManager(U.GENBANK_PREFIX, self.dir).load_genome_dirs(path),
                         {'GCA_000000003.1': '/ftp/g3'})

    def test_reports_are_named_for_the_database(self):
        # the RefSeq and GenBank runs of a release share one output directory
        refseq = UG.GenomeManager(U.REFSEQ_PREFIX, self.dir)
        genbank = UG.GenomeManager(U.GENBANK_PREFIX, self.dir)
        self.assertEqual(os.path.basename(refseq.report_file()), 'report_gcf.log')
        self.assertEqual(os.path.basename(refseq.review_file()), 'gcf_to_review.log')
        self.assertEqual(os.path.basename(genbank.report_file()), 'report_gca.log')
        self.assertEqual(os.path.basename(genbank.review_file()), 'gca_to_review.log')

    def test_construction_writes_nothing(self):
        # reports are opened by run_comparison, so a manager can be built without touching disk
        UG.GenomeManager(U.REFSEQ_PREFIX, os.path.join(self.dir, 'does_not_exist'))
        self.assertEqual(os.listdir(self.dir), [])


class RunComparisonTests(TempDirCase):
    """The whole update of one database, run dry: compared and reported, nothing copied."""

    def genome(self, root, assembly, fasta_md5):
        # the two files the comparison reads: a manifest naming the genomic FASTA
        path = os.path.join(root, assembly)
        os.makedirs(path)
        with open(os.path.join(path, 'md5checksums.txt'), 'w') as handle:
            handle.write('{}  ./{}_genomic.fna.gz\n'.format(fasta_md5, assembly))
        return path

    def test_dry_run_reports_every_genome_of_the_database_and_no_other(self):
        ftp = os.path.join(self.dir, 'mirror')
        gtdb = os.path.join(self.dir, 'previous')
        shared_old = self.genome(gtdb, 'GCF_000000001.1_ASM1v1', 'a' * 32)
        shared_new = self.genome(os.path.join(ftp, 'all', 'GCF', '000', '000', '001'),
                                 'GCF_000000001.1_ASM1v1', 'b' * 32)
        old = self.write('old_dirs.tsv',
                         'GCF_000000001.1\t{}\n'.format(shared_old) +           # shared, FASTA changed
                         'GCF_000000002.1\t/gtdb/GCF_000000002.1_ASM2v1\n'        # gone from NCBI
                         'GCA_000000004.1\t/gtdb/GCA_000000004.1_ASM4v1\n')       # other database
        new = self.write('ftp_dirs.tsv',
                         'GCF_000000001.1\t{}\n'.format(shared_new) +
                         'GCF_000000003.1\t{0}/all/GCF/000/000/003/GCF_000000003.1_ASM3v1\n'
                         'GCA_000000004.1\t{0}/all/GCA/000/000/004/GCA_000000004.1_ASM4v1\n'.format(ftp))
        out = os.path.join(self.dir, 'release')
        os.mkdir(out)

        UG.GenomeManager(U.REFSEQ_PREFIX, out, dry_run=True).run_comparison(ftp, new, old)

        with open(os.path.join(out, 'report_gcf.log')) as handle:
            rows = sorted(line.rstrip('\n').split('\t') for line in handle)
        self.assertEqual(rows, [['GCF_000000001.1', 'genomic FASTA file changed'],
                                ['GCF_000000002.1', 'removed'],
                                ['GCF_000000003.1', 'new']])
        # a dry run compares but copies nothing: only the two reports appear
        self.assertEqual(sorted(os.listdir(out)), ['gcf_to_review.log', 'report_gcf.log'])


# ----------------------------------------------------------- carrying a genome across

class ComparisonOutcome(TempDirCase):
    def test_unchanged_fasta_carries_the_derived_data_across(self):
        prev = self.genome_dir('previous', MD5_A, derived=config.GTDB_DERIVED_DIRS_TO_COPY)
        ftp = self.genome_dir('mirror', MD5_A)

        row = self.tools().compare_genome_directories(prev, ftp, self.target(), ACCESSION)

        self.assertEqual(self.status(row), UG.STATUS_FASTA_UNCHANGED)
        for derived in config.GTDB_DERIVED_DIRS_TO_COPY:
            copied = os.path.join(self.target(), derived, ACCESSION + '_derived.tsv')
            self.assertTrue(os.path.isfile(copied), derived)
        self.assertEqual(self.review.getvalue(), '')

    def test_changed_fasta_leaves_the_derived_data_behind(self):
        prev = self.genome_dir('previous', MD5_A, derived=config.GTDB_DERIVED_DIRS_TO_COPY)
        ftp = self.genome_dir('mirror', MD5_B)

        row = self.tools().compare_genome_directories(prev, ftp, self.target(), ACCESSION)

        self.assertEqual(self.status(row), UG.STATUS_FASTA_CHANGED)
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
                                 sorted([FASTA, REPORT, UG.MD5_MANIFEST]))
                with open(os.path.join(target, UG.MD5_MANIFEST)) as handle:
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

        self.assertEqual(self.status(row), UG.STATUS_FASTA_UNCHANGED)
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

        with open(os.path.join(ftp, UG.MD5_MANIFEST), 'w') as handle:
            handle.write(manifest(MD5_B))
        row = self.tools().compare_genome_directories(prev, ftp, self.target(), ACCESSION)

        self.assertEqual(self.status(row), UG.STATUS_FASTA_CHANGED)
        self.assertFalse(os.path.exists(os.path.join(self.target(), 'prodigal')))


class DryRun(TempDirCase):
    def test_dry_run_reports_the_real_outcome_but_copies_nothing(self):
        # the report of a dry run must be the report the real run would write
        for name, md5, expected in (('same', MD5_A, UG.STATUS_FASTA_UNCHANGED),
                                    ('different', MD5_B, UG.STATUS_FASTA_CHANGED)):
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
        with open(os.path.join(path, UG.MD5_MANIFEST), 'w') as handle:
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
