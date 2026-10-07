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

"""Offline unit tests for ncbi_strain_summary.py -- ncbi_strains, which writes
the strain IDs and NCBI type material status of each genome of a release.

It took the assembly summaries as exactly four files, RefSeq and GenBank of
bacteria and archaea (--rb, --ra, --gb, --ga), where select_genomes, strains
type_table and parse_ncbi_assemblies take them as -n, and found its columns from
the first '#' line after skipping one.
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
from gtdb_migration_tk.ncbi_strain_summary import NCBIStrainParser
from gtdb_migration_tk.utils.common import open_text

SUMMARY_HEADER = ['assembly_accession', 'bioproject', 'infraspecific_name', 'excluded_from_refseq',
                  'relation_to_type_material']

ASSEMBLY_REPORT = """# Assembly name:  ASM584v2
# Organism name:  Escherichia coli str. K-12 substr. MG1655 (E. coli)
# Infraspecific name:  strain=K-12
# Isolate:  na
# Taxid:          511145
"""


class TempDirCase(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='ncbi_strain_summary_test.')
        self.addCleanup(shutil.rmtree, self.dir, True)
        logger = logging.getLogger('timestamp')
        handler = logging.NullHandler()
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
        self.addCleanup(logger.removeHandler, handler)

    def summary(self, name, rows, header=SUMMARY_HEADER, compress=False):
        """An assembly summary: rows are (accession, relation_to_type_material)."""
        path = os.path.join(self.dir, name)
        lines = ['##  See ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt',
                 '#' + '\t'.join(header)]
        for accession, relation in rows:
            values = {'assembly_accession': accession, 'bioproject': 'PRJNA1', 'infraspecific_name': 'na',
                      'excluded_from_refseq': 'na', 'relation_to_type_material': relation}
            lines.append('\t'.join(values[column] for column in header))
        with (gzip.open(path, 'wt') if compress else open(path, 'w')) as handle:
            handle.write('\n'.join(lines) + '\n')
        return path

    def genome(self, gid):
        assembly = gid + '_ASM584v2'
        gpath = os.path.join(self.dir, 'genomes', assembly)
        os.makedirs(gpath)
        with open(os.path.join(gpath, assembly + '_assembly_report.txt'), 'w') as handle:
            handle.write(ASSEMBLY_REPORT)
        return gpath


class ReadingTheSummaries(TempDirCase):
    def test_every_summary_given_is_read_refseq_and_genbank_alike(self):
        refseq = self.summary('assembly_summary_bacteria_refseq.txt',
                              [('GCF_000005845.2', 'assembly from type material')])
        genbank = self.summary('assembly_summary_bacteria_genbank.txt.gz',
                               [('GCA_000005845.2', 'assembly from type material'),
                                ('GCA_000000001.1', 'na')], compress=True)
        viral = self.summary('assembly_summary_viral_genbank.txt', [('GCA_000000002.1', 'na')])

        parser = NCBIStrainParser([refseq, genbank, viral], 1)

        self.assertEqual(parser.type_material, {'GCF_000005845.2': 'assembly from type material',
                                                'GCA_000005845.2': 'assembly from type material',
                                                'GCA_000000001.1': 'na',
                                                'GCA_000000002.1': 'na'})

    def test_columns_are_found_by_name_whatever_their_order(self):
        rows = [('GCF_000005845.2', 'assembly from type material')]
        plain = self.summary('a.txt', rows)
        reordered = self.summary('b.txt', rows, header=list(reversed(SUMMARY_HEADER)))

        self.assertEqual(NCBIStrainParser([reordered], 1).type_material,
                         NCBIStrainParser([plain], 1).type_material)


class WritingTheStrainSummary(TempDirCase):
    def test_each_genome_has_its_strain_ids_and_type_material_status(self):
        summary = self.summary('assembly_summary_bacteria_refseq.txt',
                               [('GCF_000005845.2', 'assembly from type material')])
        genome_dirs = os.path.join(self.dir, 'genome_dirs.tsv')
        with open(genome_dirs, 'w') as handle:
            handle.write('GCF_000005845.2\t{}\tG000005845\n'.format(self.genome('GCF_000005845.2')))

        NCBIStrainParser([summary], 1).generate_ncbi_strains_summary(genome_dirs, self.dir)

        with open_text(os.path.join(self.dir, 'strain_summary_file.tsv.gz')) as handle:
            self.assertEqual(handle.read().splitlines(), [
                'genome_id\tOrganism name\tncbi_strain_identifiers\tncbi_type_material_designation',
                'GCF_000005845.2\tEscherichia coli str. K-12 substr. MG1655\tK-12\tassembly from type material'])


    def genome_dirs(self, *lines):
        path = os.path.join(self.dir, 'genome_dirs.tsv')
        with open(path, 'w') as handle:
            handle.write(''.join(line + '\n' for line in lines))
        return path

    def run_parser(self, summaries, genome_dirs, cpus=1):
        with self.assertLogs('timestamp', level='INFO') as logged:
            NCBIStrainParser(summaries, cpus).generate_ncbi_strains_summary(genome_dirs, self.dir)
        with open_text(os.path.join(self.dir, 'strain_summary_file.tsv.gz')) as handle:
            rows = [line.split('\t') for line in handle.read().splitlines()[1:]]
        return rows, [r for r in logged.records if r.levelno >= logging.WARNING]

    def test_a_genome_in_no_summary_has_an_empty_status_not_the_text_none(self):
        # it was written as None, which update_metadata_db would have loaded
        summary = self.summary('assembly_summary_archaea_refseq.txt', [('GCF_000000009.1', 'na')])
        genome_dirs = self.genome_dirs('GCA_000005845.2\t{}\tG000005845'.format(self.genome('GCA_000005845.2')))

        rows, warnings = self.run_parser([summary], genome_dirs)

        self.assertEqual(rows, [['GCA_000005845.2', 'Escherichia coli str. K-12 substr. MG1655', 'K-12', '']])
        self.assertEqual(len(warnings), 1)
        self.assertIn('Identified 1 genomes in none of the assembly summaries, e.g.: GCA_000005845.2',
                      warnings[0].getMessage())

    def test_a_genome_without_an_assembly_report_has_its_row_and_is_warned_of(self):
        summary = self.summary('assembly_summary_bacteria_refseq.txt', [('GCF_000000003.1', 'na')])
        gpath = os.path.join(self.dir, 'genomes', 'GCF_000000003.1_ASM3v1')
        os.makedirs(gpath)

        rows, warnings = self.run_parser([summary], self.genome_dirs('GCF_000000003.1\t{}\tG000000003'.format(gpath)))

        self.assertEqual(rows, [['GCF_000000003.1', '', '', 'na']])
        self.assertEqual(len(warnings), 1)
        self.assertIn('Identified 1 genomes with a missing _assembly_report.txt file', warnings[0].getMessage())

    def test_every_genome_is_written_whatever_order_the_workers_finish_in(self):
        genomes = ['GCF_{:09d}.1'.format(i) for i in range(40)]
        summary = self.summary('assembly_summary_bacteria_refseq.txt', [(gid, 'na') for gid in genomes])
        genome_dirs = self.genome_dirs(*('{}\t{}\tG'.format(gid, self.genome(gid)) for gid in genomes))

        rows, warnings = self.run_parser([summary], genome_dirs, cpus=4)

        self.assertEqual(sorted(row[0] for row in rows), genomes)
        self.assertEqual(warnings, [])

    def test_a_genome_that_cannot_be_read_stops_the_run_and_leaves_no_table(self):
        # the failure was swallowed: the table was left empty, or short of the
        # genomes the failed worker had, and the run reported success
        summary = self.summary('assembly_summary_bacteria_refseq.txt', [('GCF_000000001.1', 'na')])
        earlier = os.path.join(self.dir, 'strain_summary_file.tsv.gz')
        with open(earlier, 'w') as handle:
            handle.write('an earlier run\n')
        genome_dirs = self.genome_dirs('GCF_000000001.1\t{}\tG000000001'.format(self.genome('GCF_000000001.1')),
                                       'GCF_000000002.1')

        with self.assertRaises(IndexError):
            NCBIStrainParser([summary], 2).generate_ncbi_strains_summary(genome_dirs, self.dir)

        with open(earlier) as handle:
            self.assertEqual(handle.read(), 'an earlier run\n')
        self.assertEqual(sorted(os.listdir(self.dir)),
                         sorted(['assembly_summary_bacteria_refseq.txt', 'genome_dirs.tsv', 'genomes',
                                 'strain_summary_file.tsv.gz']))


class WritingItGzipped(TempDirCase):
    def test_the_table_is_gzipped_and_an_uncompressed_one_an_earlier_run_left_is_removed(self):
        open(os.path.join(self.dir, 'strain_summary_file.tsv'), 'w').close()
        summary = self.summary('assembly_summary_bacteria_refseq.txt', [('GCF_000005845.2', 'na')])
        genome_dirs = os.path.join(self.dir, 'genome_dirs.tsv')
        with open(genome_dirs, 'w') as handle:
            handle.write('GCF_000005845.2\t{}\tG000005845\n'.format(self.genome('GCF_000005845.2')))

        output = NCBIStrainParser([summary], 1).generate_ncbi_strains_summary(genome_dirs, self.dir)

        self.assertEqual(output, os.path.join(self.dir, 'strain_summary_file.tsv.gz'))
        with open(output, 'rb') as handle:
            self.assertEqual(handle.read(2), b'\x1f\x8b')
        self.assertFalse(os.path.exists(os.path.join(self.dir, 'strain_summary_file.tsv')))


class TheCommandLine(TempDirCase):
    def test_ncbi_strains_takes_the_summaries_as_n(self):
        summaries = [self.summary('assembly_summary_bacteria_refseq.txt', []),
                     self.summary('assembly_summary_bacteria_genbank.txt', [])]
        genome_dirs = os.path.join(self.dir, 'genome_dirs.tsv')
        out_dir = os.path.join(self.dir, 'ncbi_strains')
        options = main_module.get_main_parser().parse_args(
            ['ncbi_strains', '-g', genome_dirs, '-n'] + summaries
            + ['-o', out_dir, '-l', os.path.join(self.dir, 'run.log'), '-c', '3'])

        with mock.patch.object(main_py, 'NCBIStrainParser') as parser:
            main_py.OptionsParser().parse_options(options)
        parser.assert_called_once_with(summaries, 3)
        parser.return_value.generate_ncbi_strains_summary.assert_called_once_with(genome_dirs, out_dir)
        self.assertTrue(os.path.isdir(out_dir))

    def test_the_four_summary_arguments_are_no_longer_accepted(self):
        summary = self.summary('assembly_summary_bacteria_refseq.txt', [])
        with mock.patch('sys.stderr'), self.assertRaises(SystemExit) as ended:
            main_module.get_main_parser().parse_args(
                ['ncbi_strains', '-g', os.path.join(self.dir, 'genome_dirs.tsv'),
                 '--gb', summary, '--ga', summary, '--rb', summary, '--ra', summary,
                 '-o', self.dir, '-l', os.path.join(self.dir, 'run.log')])
        self.assertEqual(ended.exception.code, 2)


if __name__ == '__main__':
    unittest.main()
