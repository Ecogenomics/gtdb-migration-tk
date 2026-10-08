###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

"""Offline unit tests for ncbi_genome_category.py -- parse_ncbi_genome_category,
named ncbi_genome_category until 0.1.62, which writes the genomes NCBI marks as
a MAG, a SAG or an environmental genome.

It took one GenBank and one RefSeq summary where a release's are four; sent the
summaries' category of every genome of the release with each genome's work, a
16 MB pickle 1.35M times over r237; read every GenBank file whole; and caught
a failure with a bare except, writing the table without the genomes it lost and
exiting 0.
"""

import gzip
import logging
import os
import shutil
import tempfile
import unittest
from unittest import mock

from gtdb_migration_tk import __main__ as main_module
from gtdb_migration_tk import main as main_py
from gtdb_migration_tk import ncbi_genome_category as C
from gtdb_migration_tk.utils.common import open_text

SUMMARY_HEADER = ['assembly_accession', 'bioproject', 'excluded_from_refseq', 'relation_to_type_material']

HEADER = 'genome_id\tncbi_genome_category\tsource'


def gbff(source_qualifiers='', title='Direct Submission', after_source=''):
    """A GenBank file of one record: its header, source feature and a CDS."""
    return ('LOCUS       contig1\n'
            'REFERENCE   1  (bases 1 to 5000)\n'
            '  TITLE     {}\n'
            'FEATURES             Location/Qualifiers\n'
            '     source          1..5000\n'
            '                     /organism="Bacteroides fragilis"\n'
            '{}'
            '     CDS             1..300\n'
            '                     /transl_table=11\n'
            '{}'
            'ORIGIN\n'
            '        1 acgt\n'
            '//\n').format(title, source_qualifiers, after_source)


class TempDirCase(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='ncbi_genome_category_test.')
        self.addCleanup(shutil.rmtree, self.dir, True)
        logger = logging.getLogger('timestamp')
        handler = logging.NullHandler()
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
        self.addCleanup(logger.removeHandler, handler)
        self.out_dir = os.path.join(self.dir, 'out')
        os.makedirs(self.out_dir)
        self.output = os.path.join(self.out_dir, 'ncbi_genome_category.tsv.gz')
        self.evidence = os.path.join(self.out_dir, 'ncbi_genome_category_evidence.tsv')

    def summary(self, name, rows, header=SUMMARY_HEADER, compress=False):
        """An assembly summary: rows are (accession, excluded_from_refseq)."""
        path = os.path.join(self.dir, name)
        lines = ['##  See ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt',
                 '#' + '\t'.join(header)]
        for accession, excluded in rows:
            values = {'assembly_accession': accession, 'bioproject': 'PRJNA1',
                      'excluded_from_refseq': excluded, 'relation_to_type_material': 'na'}
            lines.append('\t'.join(values[column] for column in header))
        with (gzip.open(path, 'wt') if compress else open(path, 'w')) as handle:
            handle.write('\n'.join(lines) + '\n')
        return path

    def genome(self, gid, text=None, raw=None):
        """A genome directory, with a GenBank file of text (gzipped) and then raw bytes, or none."""
        assembly = gid + '_ASM1v1'
        gpath = os.path.join(self.dir, 'genomes', assembly)
        os.makedirs(gpath)
        if text is not None:
            with open(os.path.join(gpath, assembly + '_genomic.gbff.gz'), 'wb') as handle:
                handle.write(gzip.compress(text.encode() if isinstance(text, str) else text))
                if raw is not None:
                    handle.write(raw)
        return gpath

    def genome_dirs(self, *genomes):
        """genomes are (accession, directory)."""
        path = os.path.join(self.dir, 'genome_dirs.tsv')
        with open(path, 'w') as handle:
            for gid, gpath in genomes:
                handle.write('{}\t{}\tG{}\n'.format(gid, gpath, gid[4:13]))
        return path

    def run_command(self, summaries, genome_dirs, cpus=1):
        with self.assertLogs('timestamp', level='INFO') as logged:
            written = C.GenomeType(cpus).run(summaries, genome_dirs, self.out_dir)
        self.assertEqual(written, self.output)
        with open_text(self.output) as handle:
            lines = handle.read().splitlines()
        self.assertEqual(lines[0], HEADER)
        warnings = [r.getMessage() for r in logged.records if r.levelno >= logging.WARNING]
        return [line.split('\t') for line in lines[1:]], warnings


class TheSummariesDecideMost(TempDirCase):
    def test_every_summary_given_is_read_by_column_name_gzipped_or_not(self):
        # a release's summaries are four; two left a domain to the GenBank files
        archaea = self.summary('assembly_summary_archaea_genbank.txt.gz',
                               [('GCA_000000001.1', 'derived from metagenome')], compress=True)
        bacteria = self.summary('assembly_summary_bacteria_refseq.txt',
                                [('GCF_000000002.1', 'derived from single cell; partial')],
                                header=list(reversed(SUMMARY_HEADER)))
        env = self.summary('assembly_summary_bacteria_genbank.txt',
                           [('GCA_000000003.1', 'derived from environmental source')])
        genome_dirs = self.genome_dirs(*((gid, self.genome(gid)) for gid in
                                         ('GCA_000000001.1', 'GCF_000000002.1', 'GCA_000000003.1')))

        rows, warnings = self.run_command([archaea, bacteria, env], genome_dirs)

        # none has a GenBank file, and none is read: the summaries decided them
        self.assertEqual(rows, [['GCA_000000001.1', 'derived from metagenome', 'assembly report'],
                                ['GCA_000000003.1', 'derived from environmental sample', 'assembly report'],
                                ['GCF_000000002.1', 'derived from single cell', 'assembly report']])
        self.assertEqual(warnings, [])

    def test_a_derived_from_it_does_not_know_refuses_the_run_naming_each_value_before_anything(self):
        summary = self.summary('assembly_summary_bacteria_genbank.txt',
                               [('GCA_000000001.1', 'derived from a new source'),
                                ('GCA_000000002.1', 'derived from a new source'),
                                ('GCA_000000003.1', 'derived from somewhere else; partial'),
                                ('GCA_000000004.1', 'derived from surveillance project')])
        genome_dirs = self.genome_dirs(*((gid, self.genome(gid, raw=b'')) for gid in
                                         ('GCA_000000001.1', 'GCA_000000002.1', 'GCA_000000003.1',
                                          'GCA_000000004.1')))

        with mock.patch.object(C, 'read_gbff') as read:
            with self.assertRaisesRegex(C.GenomeCategoryError,
                                        '2 excluded_from_refseq value.*"derived from a new source" \\(2 genomes\\); '
                                        '"derived from somewhere else; partial" \\(1 genomes\\)'):
                C.GenomeType(1).run([summary], genome_dirs, self.out_dir)
        read.assert_not_called()
        self.assertEqual(os.listdir(self.out_dir), [])

    def test_a_surveillance_genome_has_no_category_from_the_summary_and_its_genbank_file_is_read(self):
        summary = self.summary('assembly_summary_bacteria_genbank.txt',
                               [('GCA_000000001.1', 'derived from surveillance project')])
        gpath = self.genome('GCA_000000001.1', gbff('                     /isolation_source="single cell"\n'))

        rows, _ = self.run_command([summary], self.genome_dirs(('GCA_000000001.1', gpath)))

        self.assertEqual(rows, [['GCA_000000001.1', 'derived from single cell', 'GBFF file']])


class TheGenBankFilesDecideTheRest(TempDirCase):
    def run_one(self, text, raw=None):
        summary = self.summary('assembly_summary_bacteria_genbank.txt', [('GCA_000000001.1', 'na')])
        gpath = self.genome('GCA_000000001.1', text, raw)
        return self.run_command([summary], self.genome_dirs(('GCA_000000001.1', gpath)))

    def test_each_category_is_found_where_r226_found_them_and_the_line_kept_beside_the_table(self):
        for qualifiers, title, category in (
                ('                     /isolation_source="single cell amplified by MDA"\n', 'x', 'single cell'),
                ('                     /note="This data was collected from a single cell"\n', 'x', 'single cell'),
                ('', 'Capturing single cell genomes', 'single cell'),
                ('                     /metagenome_source="soil metagenome"\n', 'x', 'metagenome'),
                ('                     /environmental_sample\n', 'x', 'environmental sample')):
            with self.subTest(qualifiers=qualifiers, title=title):
                shutil.rmtree(os.path.join(self.dir, 'genomes'), True)
                rows, warnings = self.run_one(gbff(qualifiers, title=title))
                self.assertEqual(rows, [['GCA_000000001.1', 'derived from ' + category, 'GBFF file']])
                self.assertEqual(warnings, [])
                with open(self.evidence) as handle:
                    evidence = handle.read().splitlines()
                self.assertEqual(evidence[0], 'genome_id\tevidence')
                self.assertEqual(len(evidence), 2)
                self.assertTrue(evidence[1].startswith('GCA_000000001.1\t'))

    def test_the_file_is_not_read_past_the_first_records_source_feature(self):
        # a full read meets the bytes that follow it, which are not gzip, and the
        # CDS's note after the source feature
        rows, _ = self.run_one(gbff(after_source='                     /note="single cell"\n'),
                               raw=b'not gzip at all' * 1000)
        self.assertEqual(rows, [])

    def test_a_byte_that_is_not_utf8_does_not_stop_the_run(self):
        text = gbff('                     /isolation_source="single cell"\n').encode().replace(
            b'Bacteroides', b'Bacteroid\xe9s')
        rows, _ = self.run_one(text)
        self.assertEqual(rows, [['GCA_000000001.1', 'derived from single cell', 'GBFF file']])

    def test_a_genome_marked_both_a_mag_and_a_sag_is_a_mag_and_warned_of(self):
        rows, warnings = self.run_one(gbff('                     /metagenome_source="gut metagenome"\n'
                                           '                     /isolation_source="single cell"\n'))
        self.assertEqual(rows, [['GCA_000000001.1', 'derived from metagenome', 'GBFF file']])
        self.assertEqual(len(warnings), 1)
        self.assertIn('Identified 1 genomes whose GenBank file marks them both a MAG and a SAG', warnings[0])

    def test_a_genome_with_no_genbank_file_has_no_row_and_is_warned_of(self):
        # the worker raised, and the genome left the table without a word
        rows, warnings = self.run_one(None)
        self.assertEqual(rows, [])
        self.assertEqual(len(warnings), 1)
        self.assertIn('Identified 1 genomes with a missing _genomic.gbff.gz file, e.g.: GCA_000000001.1',
                      warnings[0])


class WhatAWorkerIsHanded(TempDirCase):
    def test_a_worker_is_handed_the_genomes_line_alone_and_only_for_genomes_the_summaries_leave(self):
        # each was handed the summaries' category of every genome of the release
        summary = self.summary('assembly_summary_bacteria_genbank.txt',
                               [('GCA_000000001.1', 'derived from metagenome'), ('GCA_000000002.1', 'na')])
        genome_dirs = self.genome_dirs(('GCA_000000001.1', self.genome('GCA_000000001.1')),
                                       ('GCA_000000002.1', self.genome('GCA_000000002.1', gbff())))
        handed = []

        class SerialPool(object):
            def __init__(self, processes):
                pass

            def __enter__(self):
                return self

            def __exit__(self, *exc):
                return False

            def imap_unordered(self, func, iterable, chunksize=1):
                for item in iterable:
                    handed.append(item)
                    yield func(item)

        with mock.patch.object(C.mp, 'Pool', SerialPool):
            self.run_command([summary], genome_dirs)

        self.assertEqual(len(handed), 1)
        self.assertIsInstance(handed[0], str)
        self.assertTrue(handed[0].startswith('GCA_000000002.1\t'))


class TheTable(TempDirCase):
    def test_rows_are_in_accession_order_whatever_order_the_workers_finish_in(self):
        genomes = ['GCA_{:09d}.1'.format(i) for i in range(40)]
        summary = self.summary('assembly_summary_bacteria_genbank.txt', [(gid, 'na') for gid in genomes])
        genome_dirs = self.genome_dirs(*((gid, self.genome(gid, gbff('                     /note="single cell"\n')))
                                         for gid in reversed(genomes)))

        rows, warnings = self.run_command([summary], genome_dirs, cpus=4)

        self.assertEqual([row[0] for row in rows], genomes)
        self.assertEqual(warnings, [])

    def test_a_user_genome_is_a_mag_and_a_genome_in_no_summary_is_warned_of(self):
        summary = self.summary('assembly_summary_bacteria_genbank.txt', [])
        genome_dirs = self.genome_dirs(('U_00001', self.genome('U_00001')),
                                       ('GCA_000000001.1', self.genome('GCA_000000001.1', gbff())))

        rows, warnings = self.run_command([summary], genome_dirs)

        self.assertEqual(rows, [['U_00001', 'derived from metagenome', 'user genome']])
        self.assertEqual(len(warnings), 1)
        self.assertIn('Identified 1 genomes in none of the assembly summaries, e.g.: GCA_000000001.1', warnings[0])

    def test_a_genome_that_cannot_be_read_stops_the_run_and_leaves_the_earlier_table(self):
        # a bare except swallowed it, and the table was written without the genome
        summary = self.summary('assembly_summary_bacteria_genbank.txt',
                               [('GCA_000000001.1', 'na'), ('GCA_000000002.1', 'na')])
        with open(self.output, 'w') as handle:
            handle.write('an earlier run\n')
        genome_dirs = os.path.join(self.dir, 'genome_dirs.tsv')
        with open(genome_dirs, 'w') as handle:
            handle.write('GCA_000000001.1\t{}\tG000000001\n'.format(self.genome('GCA_000000001.1', gbff())))
            handle.write('GCA_000000002.1\n')

        with self.assertRaises(ValueError):
            C.GenomeType(2).run([summary], genome_dirs, self.out_dir)

        with open(self.output) as handle:
            self.assertEqual(handle.read(), 'an earlier run\n')
        self.assertEqual(os.listdir(self.out_dir), ['ncbi_genome_category.tsv.gz'])

    def test_the_table_is_gzipped_with_its_fixed_name_and_an_uncompressed_one_an_earlier_run_left_is_removed(self):
        # it was written to whatever file -o named, uncompressed
        open(os.path.join(self.out_dir, 'ncbi_genome_category.tsv'), 'w').close()
        summary = self.summary('assembly_summary_bacteria_genbank.txt', [('GCA_000000001.1', 'derived from metagenome')])

        self.run_command([summary], self.genome_dirs(('GCA_000000001.1', self.genome('GCA_000000001.1'))))

        with open(self.output, 'rb') as handle:
            self.assertEqual(handle.read(2), b'\x1f\x8b')
        self.assertEqual(sorted(os.listdir(self.out_dir)),
                         ['ncbi_genome_category.tsv.gz', 'ncbi_genome_category_evidence.tsv'])


class TheCommandLine(TempDirCase):
    def options(self, summaries):
        return main_module.get_main_parser().parse_args(
            ['parse_ncbi_genome_category', '-g', self.genome_dirs(), '-n'] + summaries
            + ['-o', os.path.join(self.dir, 'category'), '-l', os.path.join(self.dir, 'run.log'), '-c', '3'])

    def test_the_summaries_are_taken_as_n(self):
        summaries = [self.summary('assembly_summary_archaea_genbank.txt', []),
                     self.summary('assembly_summary_bacteria_genbank.txt', [])]
        with mock.patch.object(main_py, 'GenomeType') as command:
            main_py.OptionsParser().parse_options(self.options(summaries))
        command.assert_called_once_with(3)
        command.return_value.run.assert_called_once_with(
            summaries, os.path.join(self.dir, 'genome_dirs.tsv'), os.path.join(self.dir, 'category'))
        self.assertTrue(os.path.isdir(os.path.join(self.dir, 'category')))

    def test_a_value_it_cannot_place_ends_the_run_exiting_1(self):
        summary = self.summary('assembly_summary_bacteria_genbank.txt', [])
        with mock.patch.object(main_py.GenomeType, 'run', side_effect=C.GenomeCategoryError('no')), \
                self.assertLogs('timestamp', level='ERROR'), self.assertRaises(SystemExit) as ended:
            main_py.OptionsParser().parse_options(self.options([summary]))
        self.assertEqual(ended.exception.code, 1)

    def test_the_two_summary_arguments_are_no_longer_accepted(self):
        summary = self.summary('assembly_summary_bacteria_genbank.txt', [])
        with mock.patch('sys.stderr'), self.assertRaises(SystemExit) as ended:
            main_module.get_main_parser().parse_args(
                ['parse_ncbi_genome_category', '--genbank_assembly_summary', summary,
                 '--refseq_assembly_summary', summary, '-g', self.genome_dirs(), '-o', self.out_dir,
                 '-l', os.path.join(self.dir, 'run.log')])
        self.assertEqual(ended.exception.code, 2)

    def test_ncbi_genome_category_is_no_longer_a_command(self):
        summary = self.summary('assembly_summary_bacteria_genbank.txt', [])
        with mock.patch('sys.stderr'), self.assertRaises(SystemExit) as ended:
            main_module.get_main_parser().parse_args(
                ['ncbi_genome_category', '-g', self.genome_dirs(), '-n', summary, '-o', self.out_dir,
                 '-l', os.path.join(self.dir, 'run.log')])
        self.assertEqual(ended.exception.code, 2)


if __name__ == '__main__':
    unittest.main()
