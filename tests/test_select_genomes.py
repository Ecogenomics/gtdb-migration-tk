#!/usr/bin/env python3
"""Offline unit tests for select_genomes.py -- no FTP site, no genome files.

Run with the interpreter that has tqdm:

    /opt/centos7/sw/miniconda3/envs/python3/bin/python -m unittest -v tests.test_select_genomes

What is tested here is the decision about which genomes belong in a release. Reading
the NCBI assembly summary files is ncbi_utils.py's job and is tested in
tests/test_ncbi_utils.py; one test below still feeds this module a revised column
layout, to confirm the selection is reading those tables by column name.
"""

import gzip
import os
import shutil
import tempfile
import unittest

from gtdb_migration_tk import ncbi_utils as U
from gtdb_migration_tk import select_genomes as S


# an NCBI assembly summary: a comment block, the header, then the genomes
HEADER = ('#assembly_accession\tbioproject\tversion_status\tgbrs_paired_asm'
          '\texcluded_from_refseq\tftp_path')
COMMENT = '#   See ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt'


def summary(*rows, header=HEADER, comment=COMMENT):
    return '\n'.join([comment, header] + list(rows)) + '\n'


class TempDirCase(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='select_genomes_test.')

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def write(self, name, text):
        path = os.path.join(self.dir, name)
        with open(path, 'w') as handle:
            handle.write(text)
        return path


# ------------------------------------------------------------- multi-isolate projects

class MultiIsolateTagTests(unittest.TestCase):
    """NCBI renamed this annotation, so both spellings must exclude a genome."""

    def test_current_tag_is_recognised(self):
        self.assertTrue(S.is_multi_isolate('large multi-isolate project'))

    def test_renamed_tag_is_still_recognised(self):
        # summary files written before NCBI renamed the annotation
        self.assertTrue(S.is_multi_isolate('derived from surveillance project'))

    def test_tag_is_found_among_other_annotations(self):
        # the column holds several annotations separated by semicolons
        self.assertTrue(S.is_multi_isolate('partial;large multi-isolate project;contaminated'))

    def test_unannotated_genome_is_not_excluded(self):
        self.assertFalse(S.is_multi_isolate('na'))
        self.assertFalse(S.is_multi_isolate(''))

    def test_unrelated_annotation_does_not_exclude(self):
        self.assertFalse(S.is_multi_isolate('derived from metagenome;low contig N50'))


# ---------------------------------------------------------------------- the selection

class SelectedGenomesTests(TempDirCase):
    """Selecting the genomes of a release from the assembly summary files alone."""

    def write_summaries(self, refseq=(), genbank=()):
        """Write a RefSeq and a GenBank summary file, named as NCBI names them."""

        return [self.write('assembly_summary_refseq.txt', summary(*refseq)),
                self.write('assembly_summary_genbank.txt', summary(*genbank))]

    def table(self):
        """Read back the selected genome table, without its header."""

        with gzip.open(os.path.join(self.dir, 'gtdb_selected_genomes.tsv.gz'), 'rt') as handle:
            return [line.rstrip('\n').split('\t')
                    for line in handle if not line.startswith('#')]

    def select(self, refseq=(), genbank=()):
        """Run the selection over one RefSeq and one GenBank summary file."""

        S.SelectGenomes(self.dir).run(self.write_summaries(refseq, genbank))

        return self.table()

    def accessions(self, refseq=(), genbank=()):
        return [row[0] for row in self.select(refseq, genbank)]

    def test_every_refseq_genome_is_selected(self):
        self.assertEqual(
            self.accessions(refseq=['GCF_000000001.1\tPRJNA1\tlatest\tGCA_000000001.1\tna\tftp://a']),
            ['GCF_000000001.1'])

    def test_genbank_genome_without_a_refseq_counterpart_is_selected(self):
        self.assertEqual(
            self.accessions(genbank=['GCA_000000002.1\tPRJNA2\tlatest\tna\tna\tftp://b']),
            ['GCA_000000002.1'])

    def test_genbank_genome_covered_by_refseq_is_skipped(self):
        # the RefSeq copy is preferred, so the GenBank copy of G000000003 is redundant
        self.assertEqual(
            self.accessions(refseq=['GCF_000000003.1\tPRJNA3\tlatest\tGCA_000000003.1\tna\tftp://c'],
                            genbank=['GCA_000000003.1\tPRJNA3\tlatest\tGCF_000000003.1\tna\tftp://d']),
            ['GCF_000000003.1'])

    def test_coverage_is_decided_on_the_canonical_accession(self):
        # NCBI gives paired assemblies the same number but not the same version
        self.assertEqual(
            self.accessions(refseq=['GCF_000000004.2\tPRJNA4\tlatest\tGCA_000000004.1\tna\tftp://e'],
                            genbank=['GCA_000000004.1\tPRJNA4\tlatest\tGCF_000000004.2\tna\tftp://f']),
            ['GCF_000000004.2'])

    def test_genbank_genome_of_a_large_multi_isolate_project_is_skipped(self):
        self.assertEqual(
            self.accessions(genbank=['GCA_000000005.1\tPRJNA5\tlatest\tna\tlarge multi-isolate project\tftp://g']),
            [])

    def test_genbank_genome_paired_with_an_unlisted_refseq_assembly_is_selected(self):
        # the pairing names GCF_000000006.1, but no RefSeq row lists it, so the
        # genome would otherwise be lost from the release entirely
        self.assertEqual(
            self.accessions(genbank=['GCA_000000006.1\tPRJNA6\tlatest\tGCF_000000006.1\tna\tftp://h']),
            ['GCA_000000006.1'])

    def test_superseded_genomes_are_not_selected(self):
        self.assertEqual(
            self.accessions(refseq=['GCF_000000007.1\tPRJNA7\treplaced\tna\tna\tftp://i'],
                            genbank=['GCA_000000008.1\tPRJNA8\tsuppressed\tna\tna\tftp://j']),
            [])

    def test_refseq_genome_of_a_large_multi_isolate_project_is_not_selected(self):
        # NCBI is not expected to flag a RefSeq assembly this way; one that is
        # flagged must not enter the release by the route the GenBank rule closes
        self.assertEqual(
            self.accessions(refseq=['GCF_000000009.1\tPRJNA9\tlatest\tna\tlarge multi-isolate project\tftp://k']),
            [])

    def test_a_dropped_refseq_genome_does_not_cover_its_genbank_counterpart(self):
        # the RefSeq copy is not in the release, so the GenBank copy is needed
        self.assertEqual(
            self.accessions(refseq=['GCF_000000010.1\tPRJNA10\treplaced\tGCA_000000010.1\tna\tftp://l'],
                            genbank=['GCA_000000010.1\tPRJNA10\tlatest\tGCF_000000010.1\tna\tftp://m']),
            ['GCA_000000010.1'])

    def test_table_reports_the_ftp_path_and_paired_assembly(self):
        self.assertEqual(
            self.select(refseq=['GCF_000000011.1\tPRJNA11\tlatest\tGCA_000000011.1\tna\tftp://n']),
            [['GCF_000000011.1', 'ftp://n', 'latest', 'na',
              'GCA_000000011.1', 'na']])

    def test_an_unlisted_refseq_pairing_is_noted_on_the_row(self):
        # NCBI's summary files routinely name a RefSeq assembly they do not list;
        # the note is the only record that the genome was selected in spite of it
        row, = self.select(genbank=['GCA_000000020.1\tPRJNA20\tlatest\tGCF_000000020.1\tna\tftp://x'])
        self.assertEqual(row[:5], ['GCA_000000020.1', 'ftp://x', 'latest', 'na',
                                   'GCF_000000020.1'])
        self.assertEqual(row[-1],
                         'paired RefSeq assembly GCF_000000020.1 absent from the '
                         'assembly summary files')

    def test_an_ordinary_genome_carries_no_note(self):
        row, = self.select(genbank=['GCA_000000021.1\tPRJNA21\tlatest\tna\tna\tftp://y'])
        self.assertEqual(row[-1], 'na')

    def test_a_genome_covered_by_refseq_leaves_no_note_behind(self):
        # the pairing is correct here, so there is nothing to record
        rows = self.select(refseq=['GCF_000000022.1\tPRJNA22\tlatest\tGCA_000000022.1\tna\tftp://z'],
                           genbank=['GCA_000000022.1\tPRJNA22\tlatest\tGCF_000000022.1\tna\tftp://z2'])
        self.assertEqual([row[-1] for row in rows], ['na'])

    def test_table_is_sorted_so_releases_can_be_compared(self):
        self.assertEqual(
            self.accessions(refseq=['GCF_000000013.1\tPRJNA13\tlatest\tna\tna\tftp://p'],
                            genbank=['GCA_000000012.1\tPRJNA12\tlatest\tna\tna\tftp://o']),
            ['GCA_000000012.1', 'GCF_000000013.1'])

    def test_table_carries_a_commented_header(self):
        S.SelectGenomes(self.dir).run(self.write_summaries(
            refseq=['GCF_000000014.1\tPRJNA14\tlatest\tna\tna\tftp://q']))
        with gzip.open(os.path.join(self.dir, 'gtdb_selected_genomes.tsv.gz'), 'rt') as table:
            self.assertEqual(
                table.readline().rstrip('\n'),
                '#assembly_accession\tftp_path\tversion_status\texcluded_from_refseq'
                '\tgbrs_paired_asm\tnotes')

    def test_table_is_gzipped(self):
        # the files of a release are gzipped once they are part of GTDB, and the
        # table carries one row per genome in NCBI
        S.SelectGenomes(self.dir).run(self.write_summaries(
            refseq=['GCF_000000024.1\tPRJNA24\tlatest\tna\tna\tftp://ab']))
        path = os.path.join(self.dir, 'gtdb_selected_genomes.tsv.gz')
        with open(path, 'rb') as raw:
            self.assertEqual(raw.read(2), b'\x1f\x8b')          # gzip magic number
        self.assertFalse(os.path.exists(os.path.join(self.dir, 'gtdb_selected_genomes.tsv')))
    def test_selection_survives_a_revised_column_layout(self):
        # the tables are read by column name; NCBI has grown them from 23 fields to 38
        header = ('#assembly_accession\tftp_path\tsome_new_column\texcluded_from_refseq'
                  '\tgbrs_paired_asm\tversion_status')
        path = self.write('assembly_summary_genbank.txt',
                          summary('GCA_000000015.1\tftp://r\tx\tna\tna\tlatest',
                                  'GCA_000000016.1\tftp://s\tx\tlarge multi-isolate project\tna\tlatest',
                                  header=header))
        S.SelectGenomes(self.dir).run([path])
        self.assertEqual(self.table(),
                         [['GCA_000000015.1', 'ftp://r', 'latest', 'na', 'na', 'na']])

    def test_genomes_are_read_from_every_summary_file(self):
        # a release is described by one summary file per domain per database
        arc = self.write('assembly_summary_archaea_refseq.txt',
                         summary('GCF_000000017.1\tPRJNA17\tlatest\tna\tna\tftp://t'))
        bac = self.write('assembly_summary_bacteria_genbank.txt',
                         summary('GCA_000000018.1\tPRJNA18\tlatest\tna\tna\tftp://u'))
        S.SelectGenomes(self.dir).run([arc, bac])
        self.assertEqual([row[0] for row in self.table()],
                         ['GCA_000000018.1', 'GCF_000000017.1'])

    def test_refseq_is_read_before_genbank_regardless_of_file_order(self):
        # coverage cannot be decided until every RefSeq genome is known, so the
        # GenBank file being given first must not change the selection
        rfq, gbk = self.write_summaries(
            refseq=['GCF_000000019.1\tPRJNA19\tlatest\tGCA_000000019.1\tna\tftp://w'],
            genbank=['GCA_000000019.1\tPRJNA19\tlatest\tGCF_000000019.1\tna\tftp://v'])
        S.SelectGenomes(self.dir).run([gbk, rfq])
        self.assertEqual([row[0] for row in self.table()], ['GCF_000000019.1'])


class SelectedGenomesFeedTheSyncTests(TempDirCase):
    """The selected genome table is the input ncbi_genome_sync mirrors from.

    The two commands meet at this file, and they meet only through its column
    names: whatever select_genomes writes, the sync has to be able to read.
    """

    def selected(self, *rows):
        rfq = self.write('assembly_summary_refseq.txt', summary(*rows))
        S.SelectGenomes(self.dir).run([rfq])
        return os.path.join(self.dir, 'gtdb_selected_genomes.tsv.gz')

    def test_the_sync_reads_the_selected_genomes_table(self):
        import gtdb_migration_tk.ncbi_genome_sync as sync
        path = self.selected(
            'GCF_000000001.1\tPRJNA1\tlatest\tGCA_000000001.1\tna'
            '\thttps://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/000/001/GCF_000000001.1_x')
        genomes, skipped = sync.read_assembly_summary(path)
        self.assertEqual(skipped, [])
        self.assertEqual(len(genomes), 1)
        self.assertEqual(genomes[0].accession, 'GCF_000000001.1')

    def test_the_table_carries_what_assembly_status_is_rendered_from(self):
        # version_status and excluded_from_refseq are not needed to fetch a
        # genome, but assembly_status.txt is written from them
        import gtdb_migration_tk.ncbi_genome_sync as sync
        path = self.selected(
            'GCF_000000002.1\tPRJNA2\tlatest\tna\tderived from metagenome'
            '\thttps://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/000/002/GCF_000000002.1_y')
        genome, = sync.read_assembly_summary(path)[0]
        self.assertEqual(genome.version_status, 'latest')
        self.assertEqual(genome.excluded_from_refseq, 'derived from metagenome')

    def test_the_table_names_the_columns_the_sync_requires(self):
        # it refuses a table missing either, rather than guessing
        header = S.SELECTED_GENOMES_HEADER.lstrip('#').split('\t')
        for column in sync_required():
            self.assertIn(column, header)


def sync_required():
    import gtdb_migration_tk.ncbi_genome_sync as sync
    return sync.SYNC_COLUMNS


# ------------------------------------------------------------------- unserved genomes

class UnservedGenomeTests(TempDirCase):
    """Genomes NCBI serves no directory for, and the GenBank copies that rescue them."""

    SERVED = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/000/001/x'

    def write_summaries(self, refseq=(), genbank=()):
        return [self.write('assembly_summary_refseq.txt', summary(*refseq)),
                self.write('assembly_summary_genbank.txt', summary(*genbank))]

    def select(self, refseq=(), genbank=()):
        S.SelectGenomes(self.dir).run(self.write_summaries(refseq, genbank))
        with gzip.open(os.path.join(self.dir, 'gtdb_selected_genomes.tsv.gz'), 'rt') as handle:
            return [line.rstrip('\n').split('\t')
                    for line in handle if not line.startswith('#')]

    def accessions(self, refseq=(), genbank=()):
        return [row[0] for row in self.select(refseq, genbank)]

    def test_a_refseq_genome_with_no_ftp_path_is_not_selected(self):
        self.assertEqual(
            self.accessions(refseq=['GCF_000000001.1\tPRJNA1\tlatest\tna\tna\tna']),
            [])

    def test_a_genbank_genome_with_no_ftp_path_is_not_selected(self):
        self.assertEqual(
            self.accessions(genbank=['GCA_000000002.1\tPRJNA2\tlatest\tna\tna\tna']),
            [])

    def test_an_empty_ftp_path_column_is_treated_the_same(self):
        self.assertEqual(
            self.accessions(genbank=['GCA_000000003.1\tPRJNA3\tlatest\tna\tna\t']),
            [])

    def test_the_genbank_copy_is_selected_when_refseq_has_no_ftp_path(self):
        # the case this rule exists for: NCBI holds the genome under its GenBank
        # accession only, and dropping the unserved RefSeq row must not drop the genome
        rows = self.select(
            refseq=['GCF_000000004.1\tPRJNA4\tlatest\tGCA_000000004.1\tna\tna'],
            genbank=['GCA_000000004.1\tPRJNA4\tlatest\tGCF_000000004.1\tna\t' + self.SERVED])
        self.assertEqual([row[0] for row in rows], ['GCA_000000004.1'])
        self.assertEqual(rows[0][-1],
                         'paired RefSeq assembly GCF_000000004.1 was not selected: '
                         'NCBI serves no directory for it')

    def test_a_served_refseq_genome_still_covers_its_genbank_copy(self):
        # the rescue must not fire when RefSeq is perfectly usable
        self.assertEqual(
            self.accessions(
                refseq=['GCF_000000005.1\tPRJNA5\tlatest\tGCA_000000005.1\tna\t' + self.SERVED],
                genbank=['GCA_000000005.1\tPRJNA5\tlatest\tGCF_000000005.1\tna\t' + self.SERVED]),
            ['GCF_000000005.1'])

    def test_a_rescue_is_noted_differently_from_a_pairing_ncbi_never_listed(self):
        # both rows are GenBank copies of an apparently-paired genome, and only one of
        # them is NCBI's bookkeeping error; the table has to tell them apart
        rows = self.select(
            refseq=['GCF_000000006.1\tPRJNA6\tlatest\tGCA_000000006.1\tna\tna'],
            genbank=['GCA_000000006.1\tPRJNA6\tlatest\tGCF_000000006.1\tna\t' + self.SERVED,
                     'GCA_000000007.1\tPRJNA7\tlatest\tGCF_000000007.1\tna\t' + self.SERVED])
        notes = {row[0]: row[-1] for row in rows}
        self.assertIn('was not selected: NCBI serves no directory',
                      notes['GCA_000000006.1'])
        self.assertIn('absent from the assembly summary files', notes['GCA_000000007.1'])

    def test_a_multi_isolate_refseq_genome_also_names_its_reason(self):
        rows = self.select(
            refseq=['GCF_000000008.1\tPRJNA8\tlatest\tGCA_000000008.1'
                    '\tlarge multi-isolate project\t' + self.SERVED],
            genbank=['GCA_000000008.1\tPRJNA8\tlatest\tGCF_000000008.1\tna\t' + self.SERVED])
        self.assertEqual([row[0] for row in rows], ['GCA_000000008.1'])
        self.assertIn('large multi-isolate project', rows[0][-1])

    def test_an_unserved_genbank_genome_does_not_rescue_an_unserved_refseq_one(self):
        # neither is mirrorable, so the genome is simply not in the release
        self.assertEqual(
            self.accessions(
                refseq=['GCF_000000009.1\tPRJNA9\tlatest\tGCA_000000009.1\tna\tna'],
                genbank=['GCA_000000009.1\tPRJNA9\tlatest\tGCF_000000009.1\tna\tna']),
            [])

    def test_every_selected_row_carries_a_usable_ftp_path(self):
        # the contract the sync depends on
        rows = self.select(
            refseq=['GCF_000000010.1\tPRJNA10\tlatest\tna\tna\t' + self.SERVED,
                    'GCF_000000011.1\tPRJNA11\tlatest\tna\tna\tna'],
            genbank=['GCA_000000012.1\tPRJNA12\tlatest\tna\tna\t' + self.SERVED,
                     'GCA_000000013.1\tPRJNA13\tlatest\tna\tna\tna'])
        self.assertEqual([row[0] for row in rows],
                         ['GCA_000000012.1', 'GCF_000000010.1'])
        for row in rows:
            self.assertTrue(U.has_ftp_path(row[1]), row)


# --------------------------------------------------------- grouping the summary files

class SummaryFileGroupingTests(TempDirCase):
    """Telling a RefSeq summary file from a GenBank one by its name."""

    def group(self, *names):
        paths = [self.write(name, summary()) for name in names]
        return S.SelectGenomes(self.dir).group_by_database(paths)

    def test_files_are_grouped_by_the_database_they_are_named_for(self):
        refseq, genbank = self.group('assembly_summary_archaea_refseq.txt',
                                     'assembly_summary_bacteria_genbank.txt',
                                     'assembly_summary_bacteria_refseq.txt')
        self.assertEqual([os.path.basename(p) for p in refseq],
                         ['assembly_summary_archaea_refseq.txt',
                          'assembly_summary_bacteria_refseq.txt'])
        self.assertEqual([os.path.basename(p) for p in genbank],
                         ['assembly_summary_bacteria_genbank.txt'])

    def test_gzipped_files_are_grouped_by_the_same_suffix(self):
        # ncbi_metadata_sync stores them compressed; the database is still named
        # by the suffix, which now sits behind the .gz
        refseq, genbank = self.group('assembly_summary_archaea_refseq.txt.gz',
                                     'assembly_summary_bacteria_genbank.txt.gz')
        self.assertEqual([os.path.basename(p) for p in refseq],
                         ['assembly_summary_archaea_refseq.txt.gz'])
        self.assertEqual([os.path.basename(p) for p in genbank],
                         ['assembly_summary_bacteria_genbank.txt.gz'])

    def test_compressed_and_plain_files_can_be_mixed(self):
        # older releases hold these uncompressed
        refseq, genbank = self.group('assembly_summary_archaea_refseq.txt.gz',
                                     'assembly_summary_bacteria_refseq.txt')
        self.assertEqual(len(refseq), 2)
        self.assertEqual(genbank, [])

    def test_an_unrecognised_file_name_stops_the_run(self):
        # ignoring it would leave every genome it lists out of the release, which
        # is indistinguishable from a successful run of a smaller release
        with self.assertRaises(SystemExit):
            self.group('assembly_summary.txt')

    def test_only_one_database_is_allowed_but_reported(self):
        refseq, genbank = self.group('assembly_summary_refseq.txt')
        self.assertEqual(len(refseq), 1)
        self.assertEqual(genbank, [])


class MisfiledRowTests(TempDirCase):
    """A row whose accession disagrees with the file it was found in."""

    def test_a_genbank_row_in_a_refseq_file_is_not_treated_as_refseq(self):
        # otherwise it would silently mark the genome as covered by RefSeq and
        # suppress the GenBank copy that is in fact the only one available
        rfq = self.write('assembly_summary_refseq.txt',
                         summary('GCA_000000023.1\tPRJNA23\tlatest\tna\tna\tftp://aa'))
        gbk = self.write('assembly_summary_genbank.txt',
                         summary('GCA_000000023.1\tPRJNA23\tlatest\tna\tna\tftp://aa'))
        S.SelectGenomes(self.dir).run([rfq, gbk])
        with gzip.open(os.path.join(self.dir, 'gtdb_selected_genomes.tsv.gz'), 'rt') as handle:
            rows = [line.split('\t')[0] for line in handle if not line.startswith('#')]
        self.assertEqual(rows, ['GCA_000000023.1'])
