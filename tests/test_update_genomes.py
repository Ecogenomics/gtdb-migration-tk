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
import gzip
import os
import shutil
import tempfile
import unittest
import queue
import contextlib

from gtdb_migration_tk import config
from gtdb_migration_tk.biolib_lite.common import canonical_gid
from gtdb_migration_tk import ncbi_utils as U
from gtdb_migration_tk import update_genomes as UG


ACCESSION = 'GCF_000000001.1'
ASSEMBLY = ACCESSION + '_ASM1v1'
CANONICAL = 'G000000001'

# where the mirror holds ACCESSION and where the new release is to hold it: NCBI
# nests a genome by archive and accession digits under all/, and a release keeps
# that nesting with all/ replaced by the database
MIRROR_RELDIR = os.path.join('all', 'GCF', '000', '000', '001')
RELEASE_RELDIR = os.path.join('refseq', 'GCF', '000', '000', '001')
FASTA = ASSEMBLY + '_genomic.fna.gz'
REPORT = ASSEMBLY + '_assembly_report.txt'

# MD5s as NCBI would publish them; the FASTA contents are never hashed, so
# any 32 hex digits will do
MD5_A = 'a' * 32
MD5_B = 'b' * 32
MD5_R = 'c' * 32


# A genomic FASTA as NCBI serves one: a defline of contig ID then free text,
# and the sequence wrapped. The tests vary one part at a time to say which of
# them the comparison is supposed to notice.
CONTIG = 'NZ_CP007501.1'
DESCRIPTION = 'Escherichia coli K-12, complete genome'
SEQUENCE = 'ACGTACGTACGTTTGGCCTTAAGGCCTTAA'


def fasta(sequence=SEQUENCE, contig=CONTIG, description=DESCRIPTION, wrap=10):
    """A one-contig genomic FASTA, wrapped as NCBI wraps one."""
    lines = ['>{} {}'.format(contig, description)]
    lines += [sequence[i:i + wrap] for i in range(0, len(sequence), wrap)]
    return '\n'.join(lines) + '\n'


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

    def genome_dir(self, name, fasta_md5, derived=(), sequences=None):
        """A genome directory holding NCBI's files and, optionally, derived data.

        The genomic FASTA is gzipped, as NCBI serves it: the comparison reads it
        whenever the two manifests disagree.
        """
        path = os.path.join(self.dir, name, ASSEMBLY)
        os.makedirs(path)
        with gzip.open(os.path.join(path, FASTA), 'wt') as handle:
            handle.write(fasta() if sequences is None else sequences)
        for filename, text in ((REPORT, 'report'),
                               (U.MD5_MANIFEST, manifest(fasta_md5))):
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

class UpdateGenomesTests(TempDirCase):
    def setUp(self):
        super().setUp()
        self.manager = UG.UpdateGenomes(self.dir)
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

    def test_load_genome_dirs_keeps_both_databases(self):
        # one genome_dirs file describes a whole release, both databases together,
        # and both are updated in the one pass
        path = self.write('dirs.tsv',
                          'GCF_000000001.1\t/ftp/g1\n'
                          'GCF_000000002.1\t/ftp/g2\n'
                          'GCA_000000003.1\t/ftp/g3\n')
        self.assertEqual(self.manager.load_genome_dirs(path),
                         {'GCF_000000001.1': '/ftp/g1',
                          'GCF_000000002.1': '/ftp/g2',
                          'GCA_000000003.1': '/ftp/g3'})

    def test_load_genome_dirs_reports_how_many_genomes_each_database_holds(self):
        # the breakdown the two runs of this command used to give for free
        path = self.write('dirs.tsv',
                          'GCF_000000001.1\t/ftp/g1\n'
                          'GCF_000000002.1\t/ftp/g2\n'
                          'GCA_000000003.1\t/ftp/g3\n')
        with self.assertLogs('timestamp', level='INFO') as logged:
            self.manager.load_genome_dirs(path)

        identified = [line for line in logged.output if 'identified' in line][0]
        self.assertIn('3 genomes', identified)
        self.assertIn('2 RefSeq', identified)
        self.assertIn('1 GenBank', identified)

    def test_the_three_decisions_are_reported_per_database(self):
        new = {'GCF_1.1': '/ftp/f1', 'GCA_2.1': '/ftp/a2', 'GCA_3.1': '/ftp/a3'}
        old = {'GCF_1.1': '/gtdb/f1', 'GCF_4.1': '/gtdb/f4'}
        with self.assertLogs('timestamp', level='INFO') as logged:
            self.manager.generate_genomes_to_remove(new, old)
            self.manager.generate_genomes_to_add(new, old)
            self.manager.generate_genomes_to_compare(new, old)

        remove, add, compare = [line.split(':', 2)[-1] for line in logged.output]
        self.assertIn('1 RefSeq, 0 GenBank', remove)          # GCF_4.1 alone
        self.assertIn('0 RefSeq, 2 GenBank', add)             # the two GenBank genomes
        self.assertIn('1 RefSeq, 0 GenBank', compare)         # GCF_1.1 alone

    def test_the_reports_describe_the_whole_release(self):
        manager = UG.UpdateGenomes(self.dir)
        self.assertEqual(os.path.basename(manager.report_file()), 'report.log')
        self.assertEqual(os.path.basename(manager.review_file()), 'to_review.log')

    def test_construction_writes_nothing(self):
        # reports are opened by run_comparison, so a manager can be built without touching disk
        UG.UpdateGenomes(os.path.join(self.dir, 'does_not_exist'))
        self.assertEqual(os.listdir(self.dir), [])


class ReleaseFixture:
    """A mirror and a previous release holding one genome of every outcome.

    Shared by the tests of a release built once through and of one resumed, which
    need the same release: the second is the first run twice.
    """

    def genome(self, root, assembly, fasta_md5, sequences=None):
        # the two files the comparison reads: the manifest, and the genomic FASTA
        # it names, which is read whenever two manifests disagree
        path = os.path.join(root, assembly)
        os.makedirs(path)
        with open(os.path.join(path, 'md5checksums.txt'), 'w') as handle:
            handle.write('{}  ./{}_genomic.fna.gz\n'.format(fasta_md5, assembly))
        with gzip.open(os.path.join(path, assembly + '_genomic.fna.gz'), 'wt') as handle:
            handle.write(fasta() if sequences is None else sequences)
        return path

    def release(self):
        """A release with one genome of every outcome, ready to be updated.

        @return: (ftp root, ftp_genome_dirs file, old_genome_dirs file, output dir)
        """
        ftp = os.path.join(self.dir, 'mirror')
        gtdb = os.path.join(self.dir, 'previous')

        shared = self.genome(os.path.join(ftp, 'all', 'GCF', '000', '000', '001'),
                             'GCF_000000001.1_ASM1v1', 'a' * 32)
        shared_old = self.genome(gtdb, 'GCF_000000001.1_ASM1v1', 'a' * 32)
        added = self.genome(os.path.join(ftp, 'all', 'GCA', '000', '000', '002'),
                            'GCA_000000002.1_ASM2v1', 'b' * 32)
        # its previous release directory has gone, so it cannot be compared at all
        curate = self.genome(os.path.join(ftp, 'all', 'GCA', '000', '000', '004'),
                             'GCA_000000004.1_ASM4v1', 'c' * 32)

        old = self.write('old_dirs.tsv',
                         'GCF_000000001.1\t{}\n'.format(shared_old) +
                         'GCF_000000003.1\t/gone/GCF_000000003.1_ASM3v1\n'     # removed
                         'GCA_000000004.1\t/gone/GCA_000000004.1_ASM4v1\n')    # to_curate
        new = self.write('ftp_dirs.tsv',
                         'GCF_000000001.1\t{}\n'.format(shared) +
                         'GCA_000000002.1\t{}\n'.format(added) +
                         'GCA_000000004.1\t{}\n'.format(curate))
        out = os.path.join(self.dir, 'release')
        os.mkdir(out)

        return ftp, new, old, out

class RunComparisonTests(ReleaseFixture, TempDirCase):
    """The whole update of a release, run dry: compared and reported, nothing copied."""

    def targets(self, rows):
        """The release directory of each genome named by a report row."""
        return {row.split('\t')[0]: '/release/' + row.split('\t')[0] for row in rows}

    def test_dry_run_reports_every_genome_of_both_databases(self):
        # the one run updates RefSeq and GenBank together: no genome of either is
        # left out, and each is decided on its own accession
        ftp = os.path.join(self.dir, 'mirror')
        gtdb = os.path.join(self.dir, 'previous')
        shared_old = self.genome(gtdb, 'GCF_000000001.1_ASM1v1', 'a' * 32)
        shared_new = self.genome(os.path.join(ftp, 'all', 'GCF', '000', '000', '001'),
                                 'GCF_000000001.1_ASM1v1', 'b' * 32,
                                 sequences=fasta('TTTTTTTTTT'))
        old = self.write('old_dirs.tsv',
                         'GCF_000000001.1\t{}\n'.format(shared_old) +           # shared, FASTA changed
                         'GCF_000000002.1\t/gtdb/GCF_000000002.1_ASM2v1\n'        # gone from NCBI
                         'GCA_000000005.1\t/gtdb/GCA_000000005.1_ASM5v1\n')       # gone, GenBank
        new = self.write('ftp_dirs.tsv',
                         'GCF_000000001.1\t{}\n'.format(shared_new) +
                         'GCF_000000003.1\t{0}/all/GCF/000/000/003/GCF_000000003.1_ASM3v1\n'
                         'GCA_000000004.1\t{0}/all/GCA/000/000/004/GCA_000000004.1_ASM4v1\n'.format(ftp))
        out = os.path.join(self.dir, 'release')
        os.mkdir(out)

        UG.UpdateGenomes(out, dry_run=True).run_comparison(ftp, new, old)

        with open(os.path.join(out, 'report.log')) as handle:
            rows = sorted(line.rstrip('\n').split('\t') for line in handle)
        self.assertEqual(rows, [['GCA_000000004.1', 'new'],
                                ['GCA_000000005.1', 'removed'],
                                ['GCF_000000001.1', 'genomic FASTA file changed'],
                                ['GCF_000000002.1', 'removed'],
                                ['GCF_000000003.1', 'new']])
        # a dry run compares but copies nothing: only the two reports appear
        self.assertEqual(sorted(os.listdir(out)), ['report.log', 'to_review.log'])

    def test_a_release_holds_refseq_and_genbank_in_their_own_trees(self):
        ftp, new, old, out = self.release()

        UG.UpdateGenomes(out).run_comparison(ftp, new, old)

        self.assertTrue(os.path.isdir(os.path.join(
            out, 'refseq', 'GCF', '000', '000', '001', 'GCF_000000001.1_ASM1v1')))
        self.assertTrue(os.path.isdir(os.path.join(
            out, 'genbank', 'GCA', '000', '000', '002', 'GCA_000000002.1_ASM2v1')))
        # the mirror's own all/ level is not carried into the release
        self.assertFalse(os.path.exists(os.path.join(out, 'all')))

    def test_the_run_writes_the_genome_dirs_file_of_the_release(self):
        # what every later step is pointed at, so walking the tree back with
        # list_genomes would only recover what this run already knew
        ftp, new, old, out = self.release()

        UG.UpdateGenomes(out).run_comparison(ftp, new, old)

        with open(os.path.join(out, 'genome_dirs.tsv')) as handle:
            rows = sorted(line.rstrip('\n').split('\t') for line in handle)

        self.assertEqual(rows, [
            ['GCA_000000002.1',
             os.path.join(out, 'genbank', 'GCA', '000', '000', '002', 'GCA_000000002.1_ASM2v1'),
             'G000000002'],
            ['GCF_000000001.1',
             os.path.join(out, 'refseq', 'GCF', '000', '000', '001', 'GCF_000000001.1_ASM1v1'),
             'G000000001']])

    def test_the_genome_dirs_file_names_only_directories_that_are_there(self):
        # a removed genome is not in the release, and one that could not be
        # compared was never copied: naming either would hand every later step a
        # path that does not exist
        ftp, new, old, out = self.release()

        UG.UpdateGenomes(out).run_comparison(ftp, new, old)

        with open(os.path.join(out, 'genome_dirs.tsv')) as handle:
            named = [line.split('\t') for line in handle]
        self.assertEqual(sorted(gid for gid, _, _ in named),
                         ['GCA_000000002.1', 'GCF_000000001.1'])
        for _, path, _ in named:
            self.assertTrue(os.path.isdir(path), path)

    def test_the_paths_are_absolute_however_the_output_directory_was_given(self):
        ftp, new, old, out = self.release()
        relative = os.path.relpath(out, os.getcwd())

        UG.UpdateGenomes(relative).run_comparison(ftp, new, old)

        with open(os.path.join(out, 'genome_dirs.tsv')) as handle:
            for line in handle:
                self.assertTrue(os.path.isabs(line.split('\t')[1]), line)

    def test_the_genome_dirs_file_is_reported_with_its_count_per_database(self):
        ftp, new, old, out = self.release()

        with self.assertLogs('timestamp', level='INFO') as logged:
            UG.UpdateGenomes(out).run_comparison(ftp, new, old)

        wrote = [line for line in logged.output if 'Wrote' in line]
        self.assertEqual(len(wrote), 1, logged.output)
        self.assertIn('Wrote 2 genomes', wrote[0])
        self.assertIn('1 RefSeq, 1 GenBank', wrote[0])

    def test_a_dry_run_writes_no_genome_dirs_file(self):
        # it copied no directory, so it has none to name
        ftp, new, old, out = self.release()

        UG.UpdateGenomes(out, dry_run=True).run_comparison(ftp, new, old)

        self.assertEqual(sorted(os.listdir(out)), ['report.log', 'to_review.log'])

    def test_a_genome_ncbi_reissued_unchanged_goes_through_the_whole_run(self):
        # the case this check exists for, from the mirror to the reports: NCBI
        # republished the FASTA with a new description and the same sequence
        ftp_root = os.path.join(self.dir, 'mirror')
        gtdb = os.path.join(self.dir, 'previous')
        old_dir = self.genome(gtdb, 'GCF_000000001.1_ASM1v1', 'a' * 32,
                              sequences=fasta(description='Escherichia coli K-12'))
        os.makedirs(os.path.join(old_dir, 'prodigal'))
        new_dir = self.genome(os.path.join(ftp_root, 'all', 'GCF', '000', '000', '001'),
                              'GCF_000000001.1_ASM1v1', 'b' * 32,
                              sequences=fasta(description='Escherichia coli str. K-12 MG1655'))
        old = self.write('old_dirs.tsv', 'GCF_000000001.1\t{}\n'.format(old_dir))
        new = self.write('ftp_dirs.tsv', 'GCF_000000001.1\t{}\n'.format(new_dir))
        out = os.path.join(self.dir, 'release')
        os.mkdir(out)

        with self.assertLogs('timestamp', level='INFO') as logged:
            UG.UpdateGenomes(out).run_comparison(ftp_root, new, old)

        with open(os.path.join(out, 'report.log')) as handle:
            self.assertEqual(handle.read(),
                             'GCF_000000001.1\t{}\n'.format(UG.STATUS_SEQUENCES_UNCHANGED))
        # it is in the release, so it is in the genome_dirs file
        with open(os.path.join(out, 'genome_dirs.tsv')) as handle:
            self.assertEqual(handle.read().split('\t')[0], 'GCF_000000001.1')
        # with the derived data of the previous release, which is what it was for
        self.assertTrue(os.path.isdir(os.path.join(
            out, 'refseq', 'GCF', '000', '000', '001', 'GCF_000000001.1_ASM1v1', 'prodigal')))
        summary = [line for line in logged.output if 'Compared' in line][0]
        self.assertIn('1 whose sequences were unchanged despite a differing MD5', summary)

    def test_a_dry_run_reports_how_many_genomes_changed_and_how_many_did_not(self):
        # the whole question a dry run answers: how much derived data the release
        # inherits. Two shared genomes, one whose FASTA moved and one whose did not
        ftp_root = os.path.join(self.dir, 'mirror')
        gtdb = os.path.join(self.dir, 'previous')
        mirror = os.path.join(ftp_root, 'all', 'GCF', '000', '000', '001')
        same_old = self.genome(gtdb, 'GCF_000000001.1_ASM1v1', 'a' * 32)
        same_new = self.genome(mirror, 'GCF_000000001.1_ASM1v1', 'a' * 32)
        moved_old = self.genome(gtdb, 'GCF_000000002.1_ASM2v1', 'a' * 32)
        moved_new = self.genome(mirror, 'GCF_000000002.1_ASM2v1', 'b' * 32,
                                sequences=fasta('TTTTTTTTTT'))
        old = self.write('old_dirs.tsv', 'GCF_000000001.1\t{}\nGCF_000000002.1\t{}\n'.format(
            same_old, moved_old))
        new = self.write('ftp_dirs.tsv', 'GCF_000000001.1\t{}\nGCF_000000002.1\t{}\n'.format(
            same_new, moved_new))
        out = os.path.join(self.dir, 'release')
        os.mkdir(out)

        with self.assertLogs('timestamp', level='INFO') as logged:
            UG.UpdateGenomes(out, dry_run=True).run_comparison(ftp_root, new, old)

        summary = [line for line in logged.output if 'Compared' in line]
        self.assertEqual(len(summary), 1, logged.output)
        self.assertIn('DRY RUN', summary[0])
        self.assertIn('Compared 2 shared genomes', summary[0])
        self.assertIn('1 with an unchanged genomic FASTA', summary[0])
        self.assertIn('1 changed', summary[0])
        self.assertEqual(sorted(os.listdir(out)), ['report.log', 'to_review.log'])

    def test_a_real_run_reports_the_same_counts_without_the_dry_run_wording(self):
        ftp_root = os.path.join(self.dir, 'mirror')
        gtdb = os.path.join(self.dir, 'previous')
        same_old = self.genome(gtdb, 'GCF_000000001.1_ASM1v1', 'a' * 32)
        same_new = self.genome(os.path.join(ftp_root, 'all', 'GCF', '000', '000', '001'),
                               'GCF_000000001.1_ASM1v1', 'a' * 32)
        old = self.write('old_dirs.tsv', 'GCF_000000001.1\t{}\n'.format(same_old))
        new = self.write('ftp_dirs.tsv', 'GCF_000000001.1\t{}\n'.format(same_new))
        out = os.path.join(self.dir, 'release')
        os.mkdir(out)

        with self.assertLogs('timestamp', level='INFO') as logged:
            UG.UpdateGenomes(out, dry_run=False).run_comparison(ftp_root, new, old)

        summary = [line for line in logged.output if 'Compared' in line][0]
        self.assertNotIn('DRY RUN', summary)
        self.assertIn('derived data carried across', summary)

    def test_a_genome_that_cannot_be_compared_is_counted_and_warned_about(self):
        # a manifest with no genomic FASTA entry: reported to_curate, not silently lost
        ftp_root = os.path.join(self.dir, 'mirror')
        gtdb = os.path.join(self.dir, 'previous')
        broken_old = self.genome(gtdb, 'GCF_000000001.1_ASM1v1', 'a' * 32)
        broken_new = os.path.join(ftp_root, 'all', 'GCF', '000', '000', '001',
                                  'GCF_000000001.1_ASM1v1')
        os.makedirs(broken_new)
        with open(os.path.join(broken_new, 'md5checksums.txt'), 'w') as handle:
            handle.write('{}  ./something_else.txt\n'.format('b' * 32))
        old = self.write('old_dirs.tsv', 'GCF_000000001.1\t{}\n'.format(broken_old))
        new = self.write('ftp_dirs.tsv', 'GCF_000000001.1\t{}\n'.format(broken_new))
        out = os.path.join(self.dir, 'release')
        os.mkdir(out)

        with self.assertLogs('timestamp', level='INFO') as logged:
            UG.UpdateGenomes(out, dry_run=True).run_comparison(ftp_root, new, old)

        self.assertTrue(any('1 genome(s) could not be compared' in line
                            for line in logged.output), logged.output)
        with open(os.path.join(out, 'report.log')) as handle:
            self.assertIn('to_curate;ValueError', handle.read())


    def test_a_previous_release_that_has_moved_says_so_in_one_line(self):
        # the failure that prompted this: --old_genome_dirs_file naming a tree that no
        # longer exists fails every genome identically, and the run must name the cause
        # rather than leave 16,000 identical rows to be read
        ftp_root = os.path.join(self.dir, 'mirror')
        rows_old, rows_new = [], []
        for i in (1, 2, 3):
            asm = 'GCF_00000000{}.1_ASM{}v1'.format(i, i)
            mirror = self.genome(os.path.join(ftp_root, 'all', 'GCF', '000', '000', '00%d' % i),
                                 asm, 'a' * 32)
            rows_new.append('GCF_00000000{}.1\t{}\n'.format(i, mirror))
            rows_old.append('GCF_00000000{}.1\t/gone/release232/{}\n'.format(i, asm))
        old = self.write('old_dirs.tsv', ''.join(rows_old))
        new = self.write('ftp_dirs.tsv', ''.join(rows_new))
        out = os.path.join(self.dir, 'release')
        os.mkdir(out)

        with self.assertLogs('timestamp', level='INFO') as logged:
            UG.UpdateGenomes(out, dry_run=True).run_comparison(ftp_root, new, old)

        text = '\n'.join(logged.output)
        self.assertIn('3 that could not be compared', text)      # the numbers add up
        # one line for all three, naming the kind and showing a path to look at
        self.assertIn('3 x FileNotFoundError', text)
        self.assertIn('/gone/release232/', text)


    def test_the_progress_bar_carries_the_running_counts(self):
        # the comparison of a full release runs for a long time; the bar is where the
        # three numbers are visible before the summary at the end. Driven directly:
        # in a real run the bar is drawn by the listener PROCESS, whose stderr is not
        # the parent's, so it cannot be captured from there
        rows = ['GCF_000000001.1\t{}\n'.format(UG.STATUS_FASTA_UNCHANGED),
                'GCF_000000002.1\t{}\n'.format(UG.STATUS_FASTA_UNCHANGED),
                'GCF_000000004.1\t{}\n'.format(UG.STATUS_SEQUENCES_UNCHANGED),
                'GCF_000000003.1\t{}\n'.format(UG.STATUS_FASTA_CHANGED)]
        writer, tally = queue.Queue(), queue.Queue()
        for row in rows + [None]:
            writer.put(row)

        drawn = io.StringIO()
        with contextlib.redirect_stderr(drawn):
            self.tools()._FTPTools__listener(len(rows), writer, tally, self.targets(rows))

        bar = drawn.getvalue()
        self.assertIn('unchanged=2', bar)
        self.assertIn('seqs unchanged=1', bar)
        self.assertIn('changed=1', bar)
        self.assertIn('failed=0', bar)
        # named in the order the summary names them, not sorted alphabetically
        self.assertLess(bar.rindex('unchanged=2'), bar.rindex('seqs unchanged='))
        self.assertLess(bar.rindex('seqs unchanged='), bar.rindex('changed=1'))
        self.assertLess(bar.rindex('changed=1'), bar.rindex('failed='))
        # and the bar's last figures are the ones the summary is built from
        counts = tally.get().counts
        self.assertEqual(counts[UG.STATUS_FASTA_UNCHANGED], 2)
        self.assertEqual(counts[UG.STATUS_SEQUENCES_UNCHANGED], 1)
        self.assertEqual(counts[UG.STATUS_FASTA_CHANGED], 1)

    def test_the_listener_tallies_the_outcomes_of_each_database_apart(self):
        # the one figure the two runs this command used to make gave for free, and
        # the reason the counts are kept per database now that there is one run
        rows = ['GCF_000000001.1\t{}\n'.format(UG.STATUS_FASTA_UNCHANGED),
                'GCF_000000002.1\t{}\n'.format(UG.STATUS_FASTA_CHANGED),
                'GCA_000000003.1\t{}\n'.format(UG.STATUS_FASTA_UNCHANGED),
                'GCA_000000004.1\t{}\n'.format(UG.STATUS_FASTA_UNCHANGED)]
        writer, queue_out = queue.Queue(), queue.Queue()
        for row in rows + [None]:
            writer.put(row)

        with contextlib.redirect_stderr(io.StringIO()):
            self.tools()._FTPTools__listener(len(rows), writer, queue_out, self.targets(rows))

        tally = queue_out.get()
        self.assertEqual(tally.by_database['RefSeq'],
                         {UG.STATUS_FASTA_UNCHANGED: 1, UG.STATUS_FASTA_CHANGED: 1})
        self.assertEqual(tally.by_database['GenBank'],
                         {UG.STATUS_FASTA_UNCHANGED: 2})
        # and the totals are still the two added together
        self.assertEqual(tally.counts[UG.STATUS_FASTA_UNCHANGED], 3)

    def test_the_bar_counts_a_genome_that_could_not_be_compared_as_failed(self):
        rows = ['GCF_000000001.1\t{}\n'.format(
            UG.curate_status(ValueError('no entry for x_genomic.fna.gz')))]
        writer, tally = queue.Queue(), queue.Queue()
        for row in rows + [None]:
            writer.put(row)

        drawn = io.StringIO()
        with contextlib.redirect_stderr(drawn):
            self.tools()._FTPTools__listener(len(rows), writer, tally, self.targets(rows))

        self.assertIn('failed=1', drawn.getvalue())
        self.assertIn('unchanged=0', drawn.getvalue())


# ----------------------------------------------------------------- resuming a run

class ResumeTests(ReleaseFixture, TempDirCase):
    """Continuing a run that stopped part way, from the genome_dirs file it left.

    The release built by these tests holds one genome of every outcome, so a run
    of it places two genomes, reports one removed and fails to compare one. What
    a resume must then do is settled by which of those are in the genome_dirs
    file: the two that were placed, and neither of the others.
    """

    def interrupted(self):
        """A release built once through, standing in for a run that was stopped.

        Nothing here is interrupted for real -- a run killed mid-copy is not
        something a test can arrange reliably -- but what a resume reads is the
        genome_dirs file, and a complete run leaves one of exactly the same shape.
        The tests that need a genome to be missing from it take it out.

        @return: (ftp root, ftp_genome_dirs file, old_genome_dirs file, output dir)
        """

        ftp, new, old, out = self.release()
        UG.UpdateGenomes(out).run_comparison(ftp, new, old)
        return ftp, new, old, out

    def placed(self, out):
        """Accession to genome directory, as the run's genome_dirs file names them."""
        with open(os.path.join(out, 'genome_dirs.tsv')) as handle:
            return dict(line.split('\t')[:2] for line in handle)

    def rows(self, out, report='report.log'):
        with open(os.path.join(out, report)) as handle:
            return [line.rstrip('\n').split('\t') for line in handle]

    def test_without_resume_a_run_will_not_build_over_one_already_there(self):
        # what makes the flag necessary: an interrupted run cannot simply be
        # repeated, because the first genome of the add pass is already in place
        ftp, new, old, out = self.interrupted()

        with self.assertRaises(FileExistsError):
            UG.UpdateGenomes(out).run_comparison(ftp, new, old)

    def test_a_genome_already_in_the_release_is_not_handled_again(self):
        # the whole point of the flag: a genome the genome_dirs file names is left
        # exactly as it is, neither compared nor copied over
        ftp, new, old, out = self.interrupted()
        sentinels = []
        for genome_dir in self.placed(out).values():
            sentinels.append(os.path.join(genome_dir, 'untouched'))
            open(sentinels[-1], 'w').close()

        UG.UpdateGenomes(out, resume=True).run_comparison(ftp, new, old)

        for sentinel in sentinels:
            self.assertTrue(os.path.exists(sentinel), sentinel)

    def test_the_release_is_described_once_through(self):
        # the report of a resumed run is the report of the release, not of the
        # fragment the second run happened to do, and no genome is in it twice
        ftp, new, old, out = self.interrupted()
        before = self.rows(out)

        UG.UpdateGenomes(out, resume=True).run_comparison(ftp, new, old)

        after = self.rows(out)
        self.assertEqual(sorted(row[0] for row in after),
                         sorted(row[0] for row in before))
        self.assertEqual(len(after), len(set(row[0] for row in after)))

    def test_what_was_said_about_a_genome_that_finished_is_kept(self):
        ftp, new, old, out = self.interrupted()
        before = self.rows(out, 'to_review.log')
        self.assertTrue(before, 'the fixture should leave rows to review')

        UG.UpdateGenomes(out, resume=True).run_comparison(ftp, new, old)

        self.assertEqual(self.rows(out, 'to_review.log'), before)

    def test_a_genome_the_run_could_not_place_is_tried_again(self):
        # it is not in the genome_dirs file, so its row is not carried: it is
        # compared again and described by what happens to it this time
        ftp, new, old, out = self.interrupted()

        with self.assertLogs('timestamp', level='INFO') as logged:
            UG.UpdateGenomes(out, resume=True).run_comparison(ftp, new, old)

        resuming = ' '.join(line for line in logged.output if 'Resuming:' in line)
        self.assertIn('leaving 0 to add', resuming)
        self.assertIn('leaving 1 to compare', resuming)
        curate = [row for row in self.rows(out) if row[0] == 'GCA_000000004.1']
        self.assertEqual(len(curate), 1)
        self.assertTrue(curate[0][1].startswith('to_curate;'), curate)

    def test_the_genome_dirs_file_still_names_the_whole_release(self):
        # every later step is pointed at this file, so a resumed run must leave it
        # naming the release rather than the part of it the second run did
        ftp, new, old, out = self.interrupted()
        before = self.placed(out)

        UG.UpdateGenomes(out, resume=True).run_comparison(ftp, new, old)

        self.assertEqual(self.placed(out), before)

    def test_a_directory_the_run_was_part_way_through_is_replaced(self):
        # a genome being copied when the run stopped never reached the genome_dirs
        # file, so it is placed again -- over whatever the copy had got through
        ftp, new, old, out = self.interrupted()
        added = self.placed(out)['GCA_000000002.1']
        half_written = os.path.join(added, 'half_written.tmp')
        open(half_written, 'w').close()
        with open(os.path.join(out, 'genome_dirs.tsv')) as handle:
            kept = [line for line in handle if not line.startswith('GCA_000000002.1\t')]
        with open(os.path.join(out, 'genome_dirs.tsv'), 'w') as handle:
            handle.writelines(kept)

        UG.UpdateGenomes(out, resume=True).run_comparison(ftp, new, old)

        self.assertFalse(os.path.exists(half_written))
        self.assertTrue(os.path.exists(os.path.join(added, 'md5checksums.txt')))
        self.assertIn('GCA_000000002.1', self.placed(out))

    def test_a_row_that_never_finished_being_written_is_not_trusted(self):
        # the one thing a killed run can leave here; the genome it half-names is
        # placed again, and named once when it has been
        ftp, new, old, out = self.interrupted()
        genome_dirs = os.path.join(out, 'genome_dirs.tsv')
        with open(genome_dirs) as handle:
            written = handle.readlines()
        torn = written[-1].split('\t')[0]
        with open(genome_dirs, 'w') as handle:
            handle.writelines(written[:-1] + [written[-1].rstrip('\n')])

        UG.UpdateGenomes(out, resume=True).run_comparison(ftp, new, old)

        named = [row[0] for row in self.rows(out, 'genome_dirs.tsv')]
        self.assertIn(torn, named)
        self.assertEqual(named.count(torn), 1)

    def test_a_run_with_nothing_to_resume_from_builds_the_release_from_the_start(self):
        ftp, new, old, out = self.release()

        UG.UpdateGenomes(out, resume=True).run_comparison(ftp, new, old)

        self.assertEqual(sorted(self.placed(out)),
                         ['GCA_000000002.1', 'GCF_000000001.1'])

    def test_a_dry_run_reports_what_is_left_without_disturbing_the_record(self):
        # a dry run opens no genome_dirs file, so a resume can be sized without
        # putting the record it would resume from at risk
        ftp, new, old, out = self.interrupted()
        with open(os.path.join(out, 'genome_dirs.tsv')) as handle:
            before = handle.read()

        UG.UpdateGenomes(out, dry_run=True, resume=True).run_comparison(ftp, new, old)

        with open(os.path.join(out, 'genome_dirs.tsv')) as handle:
            self.assertEqual(handle.read(), before)

    def test_a_fresh_run_resumes_from_what_it_placed(self):
        # --fresh copies the whole mirror, so it is the run with the most to lose
        ftp, new, old, out = self.release()
        UG.UpdateGenomes(out).run_fresh(ftp, new)
        before = self.placed(out)
        sentinel = os.path.join(sorted(before.values())[0], 'untouched')
        open(sentinel, 'w').close()

        UG.UpdateGenomes(out, resume=True).run_fresh(ftp, new)

        self.assertTrue(os.path.exists(sentinel))
        self.assertEqual(self.placed(out), before)


class FreshRunTests(TempDirCase):
    """A release built from the mirror alone: --fresh, with no previous release."""

    def genome(self, root, assembly, fasta_md5, derived=()):
        path = os.path.join(root, assembly)
        os.makedirs(path)
        with open(os.path.join(path, 'md5checksums.txt'), 'w') as handle:
            handle.write('{}  ./{}_genomic.fna.gz\n'.format(fasta_md5, assembly))
        with gzip.open(os.path.join(path, assembly + '_genomic.fna.gz'), 'wt') as handle:
            handle.write(fasta())
        for subdir in derived:
            os.makedirs(os.path.join(path, subdir))
        return path

    def mirror(self):
        """A mirror of one RefSeq and one GenBank genome, and somewhere to build.

        @return: (ftp root, ftp_genome_dirs file, output dir)
        """
        ftp = os.path.join(self.dir, 'mirror')
        refseq = self.genome(os.path.join(ftp, 'all', 'GCF', '000', '000', '001'),
                             'GCF_000000001.1_ASM1v1', 'a' * 32)
        genbank = self.genome(os.path.join(ftp, 'all', 'GCA', '000', '000', '002'),
                              'GCA_000000002.1_ASM2v1', 'b' * 32)
        new = self.write('ftp_dirs.tsv',
                         'GCF_000000001.1\t{}\n'.format(refseq) +
                         'GCA_000000002.1\t{}\n'.format(genbank))
        out = os.path.join(self.dir, 'release')
        os.mkdir(out)

        return ftp, new, out

    def test_every_genome_of_the_mirror_is_copied_and_reported_as_new(self):
        # the whole point of --fresh: the mirror is the release, and no genome of it
        # depends on anything the previous release did or did not hold
        ftp, new, out = self.mirror()

        UG.UpdateGenomes(out).run_fresh(ftp, new)

        with open(os.path.join(out, 'report.log')) as handle:
            rows = sorted(line.rstrip('\n').split('\t') for line in handle)
        self.assertEqual(rows, [['GCA_000000002.1', 'new'],
                                ['GCF_000000001.1', 'new']])
        self.assertTrue(os.path.isfile(os.path.join(
            out, 'refseq', 'GCF', '000', '000', '001', 'GCF_000000001.1_ASM1v1',
            'md5checksums.txt')))
        self.assertTrue(os.path.isfile(os.path.join(
            out, 'genbank', 'GCA', '000', '000', '002', 'GCA_000000002.1_ASM2v1',
            'md5checksums.txt')))

    def test_no_derived_data_is_carried_across_from_a_previous_release(self):
        # a genome the previous release holds, with its Prodigal results, is still
        # taken from the mirror alone: a fresh release regenerates everything
        ftp, new, out = self.mirror()
        previous = self.genome(os.path.join(self.dir, 'previous'),
                               'GCF_000000001.1_ASM1v1', 'a' * 32,
                               derived=('prodigal',))
        self.write('old_dirs.tsv', 'GCF_000000001.1\t{}\n'.format(previous))

        UG.UpdateGenomes(out).run_fresh(ftp, new)

        self.assertFalse(os.path.exists(os.path.join(
            out, 'refseq', 'GCF', '000', '000', '001', 'GCF_000000001.1_ASM1v1',
            'prodigal')))

    def test_nothing_is_compared_or_removed(self):
        ftp, new, out = self.mirror()

        with self.assertLogs('timestamp', level='INFO') as logged:
            UG.UpdateGenomes(out).run_fresh(ftp, new)

        text = '\n'.join(logged.output)
        self.assertIn('Identified 2 genomes to add', text)
        self.assertIn('1 RefSeq, 1 GenBank', text)
        self.assertNotIn('to compare', text)
        self.assertNotIn('to remove', text)
        # to_review.log is written and stays empty: nothing was looked for
        with open(os.path.join(out, 'to_review.log')) as handle:
            self.assertEqual(handle.read(), '')

    def test_the_run_writes_the_genome_dirs_file_of_the_release(self):
        ftp, new, out = self.mirror()

        UG.UpdateGenomes(out).run_fresh(ftp, new)

        with open(os.path.join(out, 'genome_dirs.tsv')) as handle:
            rows = sorted(line.rstrip('\n').split('\t') for line in handle)

        self.assertEqual(rows, [
            ['GCA_000000002.1',
             os.path.join(out, 'genbank', 'GCA', '000', '000', '002', 'GCA_000000002.1_ASM2v1'),
             'G000000002'],
            ['GCF_000000001.1',
             os.path.join(out, 'refseq', 'GCF', '000', '000', '001', 'GCF_000000001.1_ASM1v1'),
             'G000000001']])

    def test_a_dry_run_reports_the_release_without_building_it(self):
        ftp, new, out = self.mirror()

        UG.UpdateGenomes(out, dry_run=True).run_fresh(ftp, new)

        self.assertEqual(sorted(os.listdir(out)), ['report.log', 'to_review.log'])
        with open(os.path.join(out, 'report.log')) as handle:
            self.assertEqual(len(handle.readlines()), 2)


class SequenceDigest(TempDirCase):
    """What sequences_md5() counts as the same genome, and what it does not."""

    def digest(self, name, text):
        path = os.path.join(self.dir, name)
        with gzip.open(path, 'wt') as handle:
            handle.write(text)
        return UG.sequences_md5(path)

    def test_the_description_after_the_contig_id_is_ignored(self):
        # the change this whole check exists for: NCBI renames an organism and
        # reissues a FASTA whose sequences it has not touched
        self.assertEqual(
            self.digest('a.gz', fasta(description='Escherichia coli K-12')),
            self.digest('b.gz', fasta(description='Escherichia coli str. K-12 substr. MG1655')))

    def test_line_wrapping_is_ignored(self):
        self.assertEqual(self.digest('a.gz', fasta(wrap=10)),
                         self.digest('b.gz', fasta(wrap=70)))

    def test_soft_masking_is_ignored(self):
        # lowercase marks a repeat rather than a different base, and nothing
        # called on the sequence cares which case it was given
        self.assertEqual(self.digest('a.gz', fasta(SEQUENCE)),
                         self.digest('b.gz', fasta(SEQUENCE.lower())))

    def test_a_changed_base_is_a_changed_genome(self):
        self.assertNotEqual(self.digest('a.gz', fasta(SEQUENCE)),
                            self.digest('b.gz', fasta(SEQUENCE[:-1] + 'G')))

    def test_a_renamed_contig_is_a_changed_genome(self):
        # the derived data carried across names its contigs: Prodigal's gene calls
        # give coordinates within a contig ID that must still be there
        self.assertNotEqual(self.digest('a.gz', fasta(contig='NZ_CP007501.1')),
                            self.digest('b.gz', fasta(contig='NZ_CP044123.1')))

    def test_where_one_contig_ends_and_the_next_begins_is_part_of_the_genome(self):
        # two contigs of ACGT and TTTT are not one contig of ACGTTTTT, and gene
        # calls made on either would be wrong about the other
        one = '>c1 x\nACGTTTTT\n'
        two = '>c1 x\nACGT\n>c2 y\nTTTT\n'
        self.assertNotEqual(self.digest('a.gz', one), self.digest('b.gz', two))

    def test_contig_order_is_part_of_the_genome(self):
        first = '>c1 x\nACGT\n>c2 y\nTTTT\n'
        second = '>c2 y\nTTTT\n>c1 x\nACGT\n'
        self.assertNotEqual(self.digest('a.gz', first), self.digest('b.gz', second))


class ReleaseLayout(unittest.TestCase):
    """Where the new release puts a genome the mirror holds."""

    def test_the_database_replaces_the_all_directory_of_the_mirror(self):
        self.assertEqual(
            UG.release_genome_dir('/release',
                                  'all/GCA/047/639/395/GCA_047639395.1_ASM4763939v1',
                                  'GCA_047639395.1'),
            '/release/genbank/GCA/047/639/395/GCA_047639395.1_ASM4763939v1')
        self.assertEqual(
            UG.release_genome_dir('/release',
                                  'all/GCF/000/006/805/GCF_000006805.1_ASM680v1',
                                  'GCF_000006805.1'),
            '/release/refseq/GCF/000/006/805/GCF_000006805.1_ASM680v1')

    def test_the_same_path_comes_out_whether_or_not_the_mirror_root_included_all(self):
        # --ftp_directory may name the mirror root or the all/ directory in it
        self.assertEqual(
            UG.release_genome_dir('/release', 'GCF/000/006/805/GCF_000006805.1_ASM680v1',
                                  'GCF_000006805.1'),
            UG.release_genome_dir('/release', 'all/GCF/000/006/805/GCF_000006805.1_ASM680v1',
                                  'GCF_000006805.1'))

    def test_the_nesting_is_the_mirrors_and_is_never_rebuilt_from_the_accession(self):
        # NCBI defines the layout and the sync laid it down from NCBI's own URLs;
        # a release that reshaped it would not be what --verify checks
        self.assertEqual(
            UG.release_genome_dir('/release', 'all/GCF/123/456/789/GCF_000006805.1_ASM680v1',
                                  'GCF_000006805.1'),
            '/release/refseq/GCF/123/456/789/GCF_000006805.1_ASM680v1')

    def test_a_mirror_not_laid_out_as_the_sync_lays_it_out_is_refused(self):
        # placing the genome anyway would put it somewhere unintended, and
        # silently: every genome after it would go to the same wrong place
        with self.assertRaises(ValueError):
            UG.release_genome_dir('/release', 'GCF_000006805.1_ASM680v1', 'GCF_000006805.1')
        with self.assertRaises(ValueError):
            UG.release_genome_dir('/release', 'all/XYZ/1/2/3/XYZ_1.1_ASM1v1', 'XYZ_1.1')

    def test_a_genome_dirs_row_carries_the_canonical_accession_last(self):
        # the format list_genomes writes and every later command reads
        self.assertEqual(
            UG.genome_dirs_row('GCA_003004725.1', '/release/genbank/GCA/003/004/725/x'),
            'GCA_003004725.1\t/release/genbank/GCA/003/004/725/x\tG003004725\n')


class PerDatabaseCounts(unittest.TestCase):
    """RefSeq and GenBank are updated in one pass, but counted apart."""

    def test_a_genome_is_placed_by_its_accession_prefix(self):
        self.assertEqual(UG.database_label('GCF_000000001.1'), U.REFSEQ.label)
        self.assertEqual(UG.database_label('GCA_000000001.1'), U.GENBANK.label)

    def test_an_accession_of_neither_database_is_named_as_such(self):
        # it cannot come from a mirror, but while this command ran once per prefix
        # such a row was dropped by both runs without a word
        self.assertEqual(UG.database_label('XYZ_000000001.1'), UG.UNKNOWN_DATABASE)

    def test_both_databases_are_counted_even_when_one_has_no_genomes(self):
        # a breakdown that drops an empty database reads as though it were never
        # looked at; the counts of a release are worth being seen to add up
        self.assertEqual(UG.tally_by_database(['GCF_1.1', 'GCF_2.1']),
                         {U.REFSEQ.label: 2, U.GENBANK.label: 0})

    def test_anything_unrecognised_is_counted_only_when_it_is_there(self):
        self.assertNotIn(UG.UNKNOWN_DATABASE, UG.tally_by_database(['GCA_1.1']))
        self.assertEqual(UG.tally_by_database(['GCA_1.1', 'XYZ_2.1'])[UG.UNKNOWN_DATABASE], 1)

    def test_the_breakdown_is_written_refseq_first_with_thousands_separated(self):
        counted = UG.count_by_database(['GCF_%d.1' % i for i in range(1500)] +
                                       ['GCA_%d.1' % i for i in range(2)])
        self.assertEqual(counted, '1,500 RefSeq, 2 GenBank')

    def test_the_comparison_summary_names_each_database(self):
        tally = UG.ComparisonTally(
            counts={UG.STATUS_FASTA_UNCHANGED: 3, UG.STATUS_FASTA_CHANGED: 1},
            by_database={U.REFSEQ.label: {UG.STATUS_FASTA_UNCHANGED: 1,
                                          UG.STATUS_FASTA_CHANGED: 1},
                         U.GENBANK.label: {UG.STATUS_FASTA_UNCHANGED: 2}},
            reasons={}, examples={})

        with self.assertLogs('timestamp', level='INFO') as logged:
            UG.FTPTools(io.StringIO(), io.StringIO(), True).log_comparison(tally)

        text = '\n'.join(logged.output)
        self.assertIn('Compared 4 shared genomes', text)
        self.assertIn('RefSeq: 1 unchanged, 0 sequences unchanged, 1 changed, 0 not compared', text)
        self.assertIn('GenBank: 2 unchanged, 0 sequences unchanged, 0 changed, 0 not compared', text)

    def test_every_outcome_is_named_and_the_four_add_up_to_the_total(self):
        # giving the total and then only some of the outcomes reads as a
        # contradiction, and hides how much of a release was reissued unchanged
        tally = UG.ComparisonTally(
            counts={UG.STATUS_FASTA_UNCHANGED: 5, UG.STATUS_SEQUENCES_UNCHANGED: 3,
                    UG.STATUS_FASTA_CHANGED: 2, UG.STATUS_TO_CURATE: 1},
            by_database={U.REFSEQ.label: {UG.STATUS_SEQUENCES_UNCHANGED: 3}},
            reasons={'OSError': 1}, examples={'OSError': 'OSError: [Errno 5]'})

        with self.assertLogs('timestamp', level='INFO') as logged:
            UG.FTPTools(io.StringIO(), io.StringIO(), False).log_comparison(tally)

        summary = [line for line in logged.output if 'Compared' in line][0]
        self.assertIn('Compared 11 shared genomes', summary)          # 5 + 3 + 2 + 1
        self.assertIn('5 with an unchanged genomic FASTA', summary)
        self.assertIn('3 whose sequences were unchanged despite a differing MD5', summary)
        self.assertIn('2 changed', summary)
        self.assertIn('1 that could not be compared', summary)
        # and both outcomes that keep derived data say so
        self.assertEqual(summary.count('derived data carried across'), 2)
        self.assertIn('RefSeq: 0 unchanged, 3 sequences unchanged, 0 changed, 0 not compared',
                      '\n'.join(logged.output))

    def test_a_database_with_nothing_shared_is_still_named_in_the_summary(self):
        tally = UG.ComparisonTally(counts={UG.STATUS_FASTA_UNCHANGED: 1},
                                   by_database={U.REFSEQ.label: {UG.STATUS_FASTA_UNCHANGED: 1}},
                                   reasons={}, examples={})

        with self.assertLogs('timestamp', level='INFO') as logged:
            UG.FTPTools(io.StringIO(), io.StringIO(), True).log_comparison(tally)

        self.assertIn('GenBank: 0 unchanged, 0 sequences unchanged, 0 changed, 0 not compared',
                      '\n'.join(logged.output))


class ReportOutcome(unittest.TestCase):
    """How a report row is read back for the counts."""

    def test_the_outcome_is_the_last_column(self):
        self.assertEqual(UG.report_outcome('GCF_1.1\tgenomic FASTA file unchanged\n'),
                         UG.STATUS_FASTA_UNCHANGED)
        self.assertEqual(UG.report_outcome('GCF_1.1\tnew\n'), 'new')

    def test_every_to_curate_reason_counts_as_one_outcome(self):
        # the reason varies with what went wrong; the outcome does not
        for reason in ('ValueError: no entry', 'OSError: [Errno 5]', 'KeyError'):
            self.assertEqual(UG.report_outcome('GCF_1.1\tto_curate;%s\n' % reason),
                             UG.STATUS_TO_CURATE)

    def test_failures_group_by_exception_type_not_by_message(self):
        # every message names the genome's own path, so grouping by it would give one
        # group per genome and a summary as long as the report
        first = UG.curate_reason("GCF_1.1\tto_curate;FileNotFoundError: no such file: /a\n")
        second = UG.curate_reason("GCF_2.1\tto_curate;FileNotFoundError: no such file: /b\n")
        self.assertEqual(first[0], second[0])                    # same group
        self.assertNotEqual(first[1], second[1])                 # different example

    def test_a_curate_status_carries_the_message_not_only_the_type(self):
        status = UG.curate_status(ValueError('no entry for x_genomic.fna.gz'))
        self.assertTrue(status.startswith(UG.STATUS_TO_CURATE + ';ValueError: '))
        self.assertIn('no entry for x_genomic.fna.gz', status)

    def test_a_message_can_never_break_the_report_it_is_written_into(self):
        status = UG.curate_status(OSError('line one\nline two\tthree'))
        self.assertNotIn('\n', status)
        self.assertNotIn('\t', status)


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

    def test_a_rewritten_defline_keeps_the_derived_data(self):
        # NCBI reissues a FASTA with a new organism name and an untouched
        # sequence; its published MD5 changes, and believing that alone would
        # recompute a genome that has not changed
        prev = self.genome_dir('previous', MD5_A, derived=config.GTDB_DERIVED_DIRS_TO_COPY)
        ftp = self.genome_dir('mirror', MD5_B,
                              sequences=fasta(description='Escherichia coli str. K-12'))

        row = self.tools().compare_genome_directories(prev, ftp, self.target(), ACCESSION)

        self.assertEqual(self.status(row), UG.STATUS_SEQUENCES_UNCHANGED)
        for derived in config.GTDB_DERIVED_DIRS_TO_COPY:
            copied = os.path.join(self.target(), derived, ACCESSION + '_derived.tsv')
            self.assertTrue(os.path.isfile(copied), derived)

    def test_a_rewrapped_or_remasked_fasta_keeps_the_derived_data(self):
        prev = self.genome_dir('previous', MD5_A, derived=config.GTDB_DERIVED_DIRS_TO_COPY)
        ftp = self.genome_dir('mirror', MD5_B, sequences=fasta(SEQUENCE.lower(), wrap=70))

        row = self.tools().compare_genome_directories(prev, ftp, self.target(), ACCESSION)

        self.assertEqual(self.status(row), UG.STATUS_SEQUENCES_UNCHANGED)

    def test_a_renamed_contig_leaves_the_derived_data_behind(self):
        # the bases are the same, but the gene calls name a contig the new FASTA
        # does not have, so they cannot come across
        prev = self.genome_dir('previous', MD5_A, derived=config.GTDB_DERIVED_DIRS_TO_COPY)
        ftp = self.genome_dir('mirror', MD5_B, sequences=fasta(contig='NZ_CP044123.1'))

        row = self.tools().compare_genome_directories(prev, ftp, self.target(), ACCESSION)

        self.assertEqual(self.status(row), UG.STATUS_FASTA_CHANGED)
        for derived in config.GTDB_DERIVED_DIRS_TO_COPY:
            self.assertFalse(os.path.exists(os.path.join(self.target(), derived)), derived)

    def test_an_unchanged_md5_never_opens_the_fasta_at_all(self):
        # the cheap path must stay cheap: the check runs only where the published
        # MD5s already disagree, so a FASTA that cannot even be read is no obstacle
        prev = self.genome_dir('previous', MD5_A, derived=config.GTDB_DERIVED_DIRS_TO_COPY)
        ftp = self.genome_dir('mirror', MD5_A)
        for directory in (prev, ftp):
            with open(os.path.join(directory, FASTA), 'wb') as handle:
                handle.write(b'not gzip at all')

        row = self.tools().compare_genome_directories(prev, ftp, self.target(), ACCESSION)

        self.assertEqual(self.status(row), UG.STATUS_FASTA_UNCHANGED)

    def test_changed_fasta_leaves_the_derived_data_behind(self):
        prev = self.genome_dir('previous', MD5_A, derived=config.GTDB_DERIVED_DIRS_TO_COPY)
        ftp = self.genome_dir('mirror', MD5_B, sequences=fasta('TTTTTTTTTT'))

        row = self.tools().compare_genome_directories(prev, ftp, self.target(), ACCESSION)

        self.assertEqual(self.status(row), UG.STATUS_FASTA_CHANGED)
        for derived in config.GTDB_DERIVED_DIRS_TO_COPY:
            self.assertFalse(os.path.exists(os.path.join(self.target(), derived)), derived)

    def test_every_ncbi_file_is_taken_from_the_mirror_whether_or_not_the_fasta_changed(self):
        # the manifest included: the new release must describe the files it
        # actually holds, and those are the mirror's
        for name, md5, seqs in (('same', MD5_A, None),
                                ('different', MD5_B, fasta('TTTTTTTTTT'))):
            with self.subTest(name):
                prev = self.genome_dir('previous_' + name, MD5_A)
                ftp = self.genome_dir('mirror_' + name, md5, sequences=seqs)
                target = os.path.join(self.dir, 'release_' + name, ASSEMBLY)

                self.tools().compare_genome_directories(prev, ftp, target, ACCESSION)

                self.assertEqual(sorted(os.listdir(target)),
                                 sorted([FASTA, REPORT, U.MD5_MANIFEST]))
                with open(os.path.join(target, U.MD5_MANIFEST)) as handle:
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

        with open(os.path.join(ftp, U.MD5_MANIFEST), 'w') as handle:
            handle.write(manifest(MD5_B))
        with gzip.open(os.path.join(ftp, FASTA), 'wt') as handle:
            handle.write(fasta('TTTTTTTTTT'))
        row = self.tools().compare_genome_directories(prev, ftp, self.target(), ACCESSION)

        self.assertEqual(self.status(row), UG.STATUS_FASTA_CHANGED)
        self.assertFalse(os.path.exists(os.path.join(self.target(), 'prodigal')))


class DryRun(TempDirCase):
    def test_dry_run_reports_the_real_outcome_but_copies_nothing(self):
        # the report of a dry run must be the report the real run would write
        for name, md5, seqs, expected in (
                ('same', MD5_A, None, UG.STATUS_FASTA_UNCHANGED),
                ('headers only', MD5_B, fasta(description='renamed'),
                 UG.STATUS_SEQUENCES_UNCHANGED),
                ('different', MD5_B, fasta('TTTTTTTTTT'), UG.STATUS_FASTA_CHANGED)):
            with self.subTest(name):
                prev = self.genome_dir('previous_' + name.replace(' ', '_'), MD5_A,
                                       derived=config.GTDB_DERIVED_DIRS_TO_COPY)
                ftp = self.genome_dir('mirror_' + name.replace(' ', '_'), md5,
                                      sequences=seqs)
                target = os.path.join(self.dir, 'release_' + name.replace(' ', '_'), ASSEMBLY)

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
        with open(os.path.join(path, U.MD5_MANIFEST), 'w') as handle:
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


# -------------------------------------------------------------------- adding a genome

class AddingGenomes(TempDirCase):
    """A genome new to NCBI is the mirror's directory, whole."""

    def test_every_file_of_the_mirror_directory_is_copied(self):
        # the mirror holds only what the sync fetched, so nothing is filtered:
        # the sync's .last_synced stamp and a file the old ignore list named
        # both come across, as they do for a shared genome
        mirror = self.genome_dir(os.path.join('mirror', MIRROR_RELDIR), MD5_A)
        for extra in ('.last_synced', ASSEMBLY + '_protein.faa.gz'):
            with open(os.path.join(mirror, extra), 'w') as handle:
                handle.write('x')

        self.tools().add_genomes({ACCESSION: mirror},
                                 os.path.join(self.dir, 'mirror'),
                                 os.path.join(self.dir, 'release'))

        target = os.path.join(self.dir, 'release', RELEASE_RELDIR, ASSEMBLY)
        self.assertEqual(sorted(os.listdir(target)), sorted(os.listdir(mirror)))
        self.assertEqual(self.report.getvalue(), '{}\tnew\n'.format(ACCESSION))

    def test_a_new_genome_is_placed_under_the_database_it_belongs_to(self):
        # the mirror nests both databases under one all/; a release splits them
        mirror = self.genome_dir(os.path.join('mirror', MIRROR_RELDIR), MD5_A)

        self.tools().add_genomes({ACCESSION: mirror},
                                 os.path.join(self.dir, 'mirror'),
                                 os.path.join(self.dir, 'release'))

        self.assertTrue(os.path.isdir(
            os.path.join(self.dir, 'release', 'refseq', 'GCF', '000', '000', '001', ASSEMBLY)))
        # and nothing is left at the path the mirror's own layout would have given
        self.assertFalse(os.path.exists(os.path.join(self.dir, 'release', 'all')))

    def test_a_new_genome_is_recorded_in_the_genome_dirs_file(self):
        mirror = self.genome_dir(os.path.join('mirror', MIRROR_RELDIR), MD5_A)
        genome_dirs = io.StringIO()

        UG.FTPTools(self.report, self.review, False, genome_dirs).add_genomes(
            {ACCESSION: mirror},
            os.path.join(self.dir, 'mirror'),
            os.path.join(self.dir, 'release'))

        self.assertEqual(genome_dirs.getvalue(), '{}\t{}\t{}\n'.format(
            ACCESSION,
            os.path.join(self.dir, 'release', RELEASE_RELDIR, ASSEMBLY),
            CANONICAL))

    def test_a_dry_run_reports_the_genome_and_copies_nothing(self):
        mirror = self.genome_dir(os.path.join('mirror', MIRROR_RELDIR), MD5_A)

        self.tools(dry_run=True).add_genomes({ACCESSION: mirror},
                                             os.path.join(self.dir, 'mirror'),
                                             os.path.join(self.dir, 'release'))

        self.assertFalse(os.path.exists(os.path.join(self.dir, 'release')))
        self.assertEqual(self.report.getvalue(), '{}\tnew\n'.format(ACCESSION))

    def several(self, count=8):
        """Several genomes of the mirror, as accession to directory."""
        genomes = {}
        for i in range(1, count + 1):
            accession = 'GCF_%09d.1' % i
            reldir = os.path.join('mirror', 'all', 'GCF', '000', '%03d' % i, '000')
            path = os.path.join(self.dir, reldir, accession + '_ASM%dv1' % i)
            os.makedirs(path)
            with open(os.path.join(path, 'md5checksums.txt'), 'w') as handle:
                handle.write(accession)
            genomes[accession] = path
        return genomes

    def test_copying_across_several_cpus_places_every_genome(self):
        # the copies overlap, but every genome is placed and every one of them is
        # named in the reports: a fresh release is nothing but this loop. More
        # genomes than 2 * COPY_QUEUE_DEPTH, so they are not all submitted at once
        # and the genomes copied last are submitted as the first of them finish
        genomes = self.several(count=2 * UG.COPY_QUEUE_DEPTH * 3)
        genome_dirs = io.StringIO()

        UG.FTPTools(self.report, self.review, False, genome_dirs).add_genomes(
            genomes, os.path.join(self.dir, 'mirror'),
            os.path.join(self.dir, 'release'), cpus=2)

        placed = sorted(line.split('\t')[0] for line in
                        genome_dirs.getvalue().splitlines())
        self.assertEqual(placed, sorted(genomes))
        self.assertEqual(sorted(row.split('\t')[0] for row in
                                self.report.getvalue().splitlines()),
                         sorted(genomes))
        for i in range(1, len(genomes) + 1):
            self.assertTrue(os.path.isfile(os.path.join(
                self.dir, 'release', 'refseq', 'GCF', '000', '%03d' % i, '000',
                'GCF_%09d.1_ASM%dv1' % (i, i), 'md5checksums.txt')))

    def test_the_paths_written_are_the_paths_of_their_own_genomes(self):
        # the one way a threaded copy could go wrong quietly: a genome_dirs row
        # pairing one accession with another genome's directory
        genomes = self.several()
        genome_dirs = io.StringIO()

        UG.FTPTools(self.report, self.review, False, genome_dirs).add_genomes(
            genomes, os.path.join(self.dir, 'mirror'),
            os.path.join(self.dir, 'release'), cpus=4)

        for line in genome_dirs.getvalue().splitlines():
            gid, path, canonical = line.split('\t')
            self.assertEqual(os.path.basename(path).split('_ASM')[0], gid)
            self.assertEqual(canonical, canonical_gid(gid))
            # and it is the directory of THAT genome, as the mirror held it
            with open(os.path.join(path, 'md5checksums.txt')) as handle:
                self.assertEqual(handle.read(), gid)

    def test_a_genome_that_will_not_copy_stops_the_run(self):
        # a release quietly short of a genome is worse than one that did not
        # finish being built, so the copy's exception is raised
        genomes = self.several(count=3)
        genomes['GCF_000000009.1'] = os.path.join(self.dir, 'mirror', 'all', 'GCF',
                                                  '000', '009', '000', 'gone')

        with self.assertRaises(OSError):
            self.tools().add_genomes(genomes, os.path.join(self.dir, 'mirror'),
                                     os.path.join(self.dir, 'release'), cpus=4)

        # the report still names every genome the run handled, the failed one included
        self.assertEqual(sorted(row.split('\t')[0] for row in
                                self.report.getvalue().splitlines()),
                         sorted(genomes))
