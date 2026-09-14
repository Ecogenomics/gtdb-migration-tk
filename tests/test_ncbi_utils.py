#!/usr/bin/env python3
"""Offline unit tests for ncbi_utils.py -- the shared reader for NCBI assembly summaries.

The contract that would break silently in production, and nowhere else, is column
lookup: NCBI has grown assembly_summary.txt from 23 fields to 38, so a reader that
addresses columns by position starts reading the wrong field whenever the table is
revised. ncbi_genome_sync then mirrors the wrong files, and ncbi_ftp_manager builds a release
from the wrong set of genomes, with nothing to indicate anything went wrong. Those
tests feed the reader tables whose columns have moved.
"""

import os
import shutil
import tempfile
import unittest

from gtdb_migration_tk import ncbi_utils as U


# an NCBI assembly summary: a comment block, the header, then the genomes
HEADER = ('#assembly_accession\tbioproject\tversion_status\tgbrs_paired_asm'
          '\texcluded_from_refseq\tftp_path')
COMMENT = '#   See ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt'


def summary(*rows, header=HEADER, comment=COMMENT):
    return '\n'.join([comment, header] + list(rows)) + '\n'


class TempDirCase(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='ncbi_utils_test.')

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def write(self, name, text):
        path = os.path.join(self.dir, name)
        with open(path, 'w') as handle:
            handle.write(text)
        return path


# -------------------------------------------------------------------- counting rows

class CountSummaryRows(TempDirCase):
    """The total a progress bar is sized from must count what the reader yields."""

    def summary(self, *rows):
        return self.write('s.txt', COMMENT + '\n' + HEADER + '\n'
                          + ''.join(row + '\n' for row in rows))

    def test_counts_genomes_not_lines(self):
        # the '#' comment and header are not genomes; counting lines instead is what
        # left the bar two short of its own total and never reaching 100%
        path = self.summary('GCA_1.1\tftp://a\tlatest', 'GCA_2.1\tftp://b\tlatest')
        self.assertEqual(U.count_summary_rows(path), 2)

    def test_a_summary_with_no_genomes_counts_zero(self):
        self.assertEqual(U.count_summary_rows(self.summary()), 0)

    def test_blank_lines_are_not_genomes(self):
        path = self.write('s.txt', COMMENT + '\n' + HEADER + '\n\nGCA_1.1\tftp://a\tlatest\n\n')
        self.assertEqual(U.count_summary_rows(path), 1)

    def test_it_counts_exactly_what_the_reader_yields(self):
        # the contract: bar total == iterations, whatever the file's shape
        path = self.write('s.txt', COMMENT + '\n' + HEADER + '\n'
                          + 'GCA_1.1\tftp://a\tlatest\n\n'
                          + COMMENT + '\n' + HEADER + '\n'      # a concatenated table
                          + 'GCA_2.1\tftp://b\tlatest\n')
        yielded = sum(1 for _ in U.read_assembly_summary(path, 'assembly_accession'))
        self.assertEqual(U.count_summary_rows(path), yielded)

    def test_a_gzipped_summary_is_counted_too(self):
        import gzip
        path = os.path.join(self.dir, 's.txt.gz')
        with gzip.open(path, 'wt') as handle:
            handle.write(COMMENT + '\n' + HEADER + '\nGCA_1.1\tftp://a\tlatest\n')
        self.assertEqual(U.count_summary_rows(path), 1)


# ------------------------------------------------------------------ column lookup

class SummaryColumns(unittest.TestCase):
    def test_header_is_located_by_content_not_position(self):
        # the comment block is '#'-prefixed too, and must not be taken for the header
        self.assertIsNone(U.summary_columns(COMMENT))
        self.assertEqual(U.summary_columns(HEADER)['assembly_accession'], 0)
        self.assertEqual(U.summary_columns(HEADER)['version_status'], 2)

    def test_a_genome_row_is_not_a_header(self):
        self.assertIsNone(U.summary_columns('GCA_000001405.1\tPRJNA31257\tlatest'))

    def test_header_is_recognised_by_either_naming_column(self):
        # the outputs ncbi_genome_sync writes carry a subset of the NCBI columns
        self.assertEqual(U.summary_columns('#ftp_path\tbioproject\tassembly_accession'),
                         {'ftp_path': 0, 'bioproject': 1, 'assembly_accession': 2})

    def test_header_written_with_a_space_after_the_hash_is_read(self):
        # NCBI wrote the marker as '# assembly_accession' until 2020, and the
        # archived summary files of earlier releases are still read
        self.assertEqual(U.summary_columns('# assembly_accession\tbioproject\tftp_path'),
                         {'assembly_accession': 0, 'bioproject': 1, 'ftp_path': 2})

    def test_an_old_style_header_satisfies_its_required_columns(self):
        # it names the accession column, so it must not be rejected as missing it
        cols = U.summary_columns('# assembly_accession\tftp_path',
                                 required=('assembly_accession',))
        self.assertEqual(cols['assembly_accession'], 0)

    def test_header_missing_a_required_column_is_rejected(self):
        # returning None would leave the caller reading the rest of the file as headerless
        with self.assertRaises(U.BadInput):
            U.summary_columns('#ftp_path\tbioproject', required=('assembly_accession',))

    def test_required_columns_present_are_accepted(self):
        cols = U.summary_columns(HEADER, required=('assembly_accession', 'ftp_path'))
        self.assertEqual(cols['ftp_path'], 5)


# ------------------------------------------------------------------ field lookup

class SummaryField(unittest.TestCase):
    COLUMNS = {'assembly_accession': 0, 'version_status': 1, 'ftp_path': 2}

    def test_absent_column_reads_as_empty(self):
        self.assertEqual(U.summary_field(['GCA_1.1'], self.COLUMNS, 'excluded_from_refseq'), '')

    def test_row_stopping_short_of_the_column_reads_as_empty(self):
        self.assertEqual(U.summary_field(['GCA_1.1'], self.COLUMNS, 'ftp_path'), '')

    def test_value_is_stripped(self):
        self.assertEqual(U.summary_field(['GCA_1.1', ' latest\n'], self.COLUMNS, 'version_status'),
                         'latest')


# ------------------------------------------------------------------ reading a table

class ReadAssemblySummary(TempDirCase):
    def read(self, text, *names):
        return list(U.read_assembly_summary(self.write('summary.txt', text), *names))

    def test_reads_named_fields(self):
        rows = self.read(summary('GCA_000001405.1\tPRJNA31257\tlatest\tGCF_000001405.40\tna\tftp://x'),
                         'assembly_accession', 'version_status')
        self.assertEqual(rows, [('GCA_000001405.1', 'latest')])

    def test_columns_are_located_by_name_when_the_table_is_revised(self):
        # NCBI inserts two columns ahead of version_status; a positional reader breaks here
        header = ('#assembly_accession\tbioproject\tbiosample\twgs_master'
                  '\tversion_status\tgbrs_paired_asm\texcluded_from_refseq')
        rows = self.read(summary('GCA_000001405.1\tPRJNA31257\tSAMN3\tABCD01\tlatest\tGCF_000001405.40\tna',
                                 header=header),
                         'assembly_accession', 'version_status', 'gbrs_paired_asm')
        self.assertEqual(rows, [('GCA_000001405.1', 'latest', 'GCF_000001405.40')])

    def test_absent_column_reads_as_empty(self):
        # older summary files do not carry excluded_from_refseq
        header = '#assembly_accession\tversion_status'
        rows = self.read(summary('GCA_000001405.1\tlatest', header=header),
                         'assembly_accession', 'excluded_from_refseq')
        self.assertEqual(rows, [('GCA_000001405.1', '')])

    def test_trailing_whitespace_is_stripped(self):
        # version_status as the final column previously compared as 'latest\n' and lost the genome
        header = '#assembly_accession\tversion_status'
        rows = self.read(summary('GCA_000001405.1\tlatest', header=header),
                         'version_status')
        self.assertEqual(rows, [('latest',)])

    def test_blank_lines_are_skipped(self):
        rows = self.read(summary('GCA_000001405.1\tPRJNA31257\tlatest\tna\tna\tftp://x', ''),
                         'assembly_accession')
        self.assertEqual(rows, [('GCA_000001405.1',)])

    def test_table_without_a_header_is_refused(self):
        # guessing at column positions is how a summary revision corrupts a release
        path = self.write('summary.txt', 'GCA_000001405.1\tPRJNA31257\tlatest\n')
        with self.assertRaises(U.BadInput) as ctx:
            list(U.read_assembly_summary(path, 'assembly_accession'))
        self.assertIn('summary.txt:1', str(ctx.exception))

    def test_refusal_is_a_value_error(self):
        # callers that predate BadInput catch ValueError
        self.assertTrue(issubclass(U.BadInput, ValueError))


class ReadSummaryRows(TempDirCase):
    def test_line_numbers_count_comments_and_blanks(self):
        # a row is reported to the user by its line in the file as shipped
        path = self.write('summary.txt', summary('GCA_1.1\tx\tlatest\tna\tna\tftp://x'))
        rows = list(U.read_summary_rows(path))
        self.assertEqual([lineno for lineno, _, _ in rows], [3])

    def test_required_columns_are_enforced_at_the_header(self):
        path = self.write('summary.txt', summary('GCA_1.1\tlatest',
                                                 header='#assembly_accession\tversion_status'))
        with self.assertRaises(U.BadInput):
            list(U.read_summary_rows(path, required=('ftp_path',)))


if __name__ == '__main__':
    unittest.main()
