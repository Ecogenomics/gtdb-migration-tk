#!/usr/bin/env python3
"""Offline unit tests for ncbi_utils.py -- the shared reader for NCBI assembly summaries.

The contract that would break silently in production, and nowhere else, is column
lookup: NCBI has grown assembly_summary.txt from 23 fields to 38, so a reader that
addresses columns by position starts reading the wrong field whenever the table is
revised. ncbi_genome_sync then mirrors the wrong files, and select_genomes builds a release
from the wrong set of genomes, with nothing to indicate anything went wrong. Those
tests feed the reader tables whose columns have moved.

Also here is what the NCBI commands know in common about NCBI's files -- the
accession prefixes, the manifest line, the test for an unserved genome -- because
ncbi_utils.py is the one module all of them may import.
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


# --------------------------------------------------------------------------- ftp_path

class FtpPathTests(unittest.TestCase):
    """A genome NCBI lists but does not serve cannot be mirrored."""

    def test_a_served_genome_has_an_ftp_path(self):
        self.assertTrue(U.has_ftp_path('https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/x'))

    def test_na_is_not_an_ftp_path(self):
        self.assertFalse(U.has_ftp_path('na'))

    def test_na_is_matched_whatever_its_case(self):
        self.assertFalse(U.has_ftp_path('NA'))

    def test_an_empty_column_is_not_an_ftp_path(self):
        # older summary files leave it empty rather than writing na
        self.assertFalse(U.has_ftp_path(''))



class FtpPathFeedsTheSyncTests(TempDirCase):
    """select_genomes keeps the rows has_ftp_path accepts; the sync must fetch
    exactly those, or the selection promises genomes the sync then refuses."""

    ROWS = {'na': 'GCF_000000001.1', 'NA': 'GCF_000000002.1', '': 'GCF_000000003.1',
            'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/000/004/x': 'GCF_000000004.1'}

    def test_the_sync_skips_exactly_the_rows_has_ftp_path_rejects(self):
        import gtdb_migration_tk.ncbi_genome_sync as sync
        table = self.write('summary.txt', summary(*(
            '{}\tPRJ\tlatest\tna\tna\t{}'.format(accession, ftp_path)
            for ftp_path, accession in self.ROWS.items())))

        genomes, skipped = sync.read_assembly_summary(table)

        self.assertEqual(sorted(g.accession for g in genomes),
                         sorted(a for f, a in self.ROWS.items() if U.has_ftp_path(f)))
        self.assertEqual(sorted(a for _, a, _ in skipped),
                         sorted(a for f, a in self.ROWS.items() if not U.has_ftp_path(f)))


# ----------------------------------------------------------------- the NCBI databases

class DatabaseTableTests(unittest.TestCase):
    """RefSeq is GCF and its summary files end _refseq.txt: said once, here."""

    def test_refseq_comes_before_genbank(self):
        # select_genomes must know what RefSeq covers before it judges GenBank
        self.assertEqual([d.name for d in U.NCBI_DATABASES], ['refseq', 'genbank'])

    def test_the_prefixes_are_the_databases_own(self):
        self.assertEqual(U.REFSEQ_PREFIX, U.REFSEQ.prefix)
        self.assertEqual(U.GENBANK_PREFIX, U.GENBANK.prefix)
        self.assertEqual({U.REFSEQ_PREFIX, U.GENBANK_PREFIX}, {'GCF', 'GCA'})


class SummaryFileNamingTests(unittest.TestCase):
    """ncbi_metadata_sync names the files and select_genomes reads the database
    back out of the name; one convention, two directions."""

    def test_a_name_reads_back_as_the_database_it_was_built_for(self):
        for database in U.NCBI_DATABASES:
            for domain in ('archaea', 'bacteria', 'fungi'):
                name = U.assembly_summary_filename(domain, database)
                self.assertIs(U.assembly_summary_database(name), database, name)

    def test_the_names_are_the_ones_gtdb_has_always_used(self):
        self.assertEqual(U.assembly_summary_filename('bacteria', U.REFSEQ),
                         'assembly_summary_bacteria_refseq.txt.gz')
        self.assertEqual(U.assembly_summary_filename('archaea', U.GENBANK),
                         'assembly_summary_archaea_genbank.txt.gz')

    def test_an_uncompressed_file_from_an_older_release_is_placed_too(self):
        self.assertIs(U.assembly_summary_database('assembly_summary_bacteria_refseq.txt'),
                      U.REFSEQ)

    def test_a_path_is_read_by_its_last_component(self):
        self.assertIs(U.assembly_summary_database('/r237/ncbi/assembly_summary_fungi_genbank.txt.gz'),
                      U.GENBANK)

    def test_a_name_saying_neither_database_is_none(self):
        # select_genomes stops on such a file rather than guessing
        for name in ('assembly_summary.txt', 'assembly_summary_bacteria.txt.gz',
                     'genbank_summary.txt', 'assembly_summary_refseq_bacteria.txt'):
            self.assertIsNone(U.assembly_summary_database(name), name)


# ----------------------------------------------------------------- the genome columns

class GenomeColumnsTests(unittest.TestCase):
    """Every table GTDB hands the sync opens with the same four columns; they are
    written down once and each table's header is built from them."""

    def test_the_columns_are_what_the_sync_needs_to_mirror_a_genome(self):
        self.assertEqual(U.GENOME_COLUMNS,
                         ('assembly_accession', 'ftp_path', 'version_status',
                          'excluded_from_refseq'))

    def test_a_table_header_is_a_comment_line_naming_the_columns(self):
        # the summary readers skip comment lines and read column names from them
        self.assertEqual(U.table_header('a', 'b'), '#a\tb')
        columns = U.summary_columns(U.table_header(*U.GENOME_COLUMNS))
        self.assertEqual(sorted(columns), sorted(U.GENOME_COLUMNS))

    def test_the_selection_opens_with_the_genome_columns(self):
        from gtdb_migration_tk import select_genomes
        self.assertTrue(select_genomes.SELECTED_GENOMES_HEADER.startswith(
            U.table_header(*U.GENOME_COLUMNS) + '\t'))

    def test_the_sync_writes_its_own_tables_with_the_genome_columns(self):
        from gtdb_migration_tk import ncbi_genome_sync as sync
        # both open with the shared columns, and each adds the one column saying what
        # went wrong -- which --retry, reading by name, ignores
        for header in (sync.BAD_HEADER, sync.FAIL_HEADER):
            self.assertTrue(header.startswith(U.table_header(*U.GENOME_COLUMNS) + '\t'))
        self.assertEqual(sync.BAD_HEADER.split('\t')[-1], 'failed_files\n')
        self.assertEqual(sync.FAIL_HEADER.split('\t')[-1], 'reason\n')
        self.assertEqual(len(sync.Genome._fields), len(U.GENOME_COLUMNS))
        self.assertEqual(tuple(sync.SYNC_COLUMNS), U.GENOME_COLUMNS[:2])


class NcbiNullTests(unittest.TestCase):
    def test_na_is_the_null_the_selection_writes_and_has_ftp_path_reads(self):
        from gtdb_migration_tk import select_genomes
        self.assertEqual(U.NCBI_NA, 'na')
        self.assertEqual(select_genomes.NO_NOTE, U.NCBI_NA)
        self.assertFalse(U.has_ftp_path(U.NCBI_NA))


# ----------------------------------------------------------------------- the manifest

class ManifestTests(unittest.TestCase):
    """md5checksums.txt is read once, for the sync and for the release update."""

    MD5 = 'a' * 32

    def read(self, text):
        return list(U.read_md5_manifest(text.splitlines()))

    def test_an_entry_is_its_md5_and_its_name_less_the_leading_dot_slash(self):
        self.assertEqual(self.read(self.MD5 + '  ./GCF_1_A_genomic.fna.gz\n'),
                         [(self.MD5, 'GCF_1_A_genomic.fna.gz')])

    def test_a_nested_entry_keeps_its_directory(self):
        # callers match the whole path, which is what keeps the subtrees out
        self.assertEqual(self.read(self.MD5 + '  ./all_assembly_versions/x.txt\n'),
                         [(self.MD5, 'all_assembly_versions/x.txt')])

    def test_lines_that_are_not_entries_are_skipped(self):
        text = '\n# a comment\nnot a checksum  ./x\n' + self.MD5 + '  ./y\n'
        self.assertEqual(self.read(text), [(self.MD5, 'y')])

    def test_manifest_order_is_kept(self):
        text = ''.join('{}  ./{}\n'.format(c * 32, n) for c, n in (('1', 'b'), ('2', 'a')))
        self.assertEqual([n for _, n in self.read(text)], ['b', 'a'])

    def test_the_manifest_and_fasta_names_are_those_the_sync_fetches(self):
        from gtdb_migration_tk import ncbi_genome_sync as sync
        self.assertIn(U.MD5_MANIFEST, sync.WANTED_EXACT)
        self.assertIn(U.GENOMIC_FASTA_EXT, sync.WANTED_SUFFIXES)


# ---------------------------------------------------------------------------- hashing

class FileMd5Tests(TempDirCase):
    """One file hash for the taxonomy dump and the mirror alike."""

    def test_the_digest_is_the_md5_of_the_whole_file(self):
        import hashlib
        path = self.write('small.bin', 'hello ncbi')
        self.assertEqual(U.file_md5(path), hashlib.md5(b'hello ncbi').hexdigest())

    def test_a_file_longer_than_a_chunk_hashes_as_one_stream(self):
        import hashlib
        payload = bytes(range(256)) * (2 * U.CHUNK // 256 + 3)
        path = os.path.join(self.dir, 'big.bin')
        with open(path, 'wb') as handle:
            handle.write(payload)
        self.assertGreater(len(payload), 2 * U.CHUNK)
        self.assertEqual(U.file_md5(path), hashlib.md5(payload).hexdigest())

    def test_the_chunk_is_one_mebibyte(self):
        # the measured flat region; the sync's memory bound is CHUNK x threads
        self.assertEqual(U.CHUNK, 1024 * 1024)
