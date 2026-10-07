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

"""Offline unit tests for metadata_ncbi_manager.py -- parse_assemblies, which
writes the NCBI metadata of a release's genomes from the assembly summaries.

It took exactly four summaries, RefSeq and GenBank of bacteria and archaea, as
--rb, --ra, --gb and --ga, where select_genomes and strains type_table take the
summaries as -n, and found its columns from the second line of each file,
writing the header of the first: a summary whose columns were in another order
had its values written under the wrong fields.
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
from gtdb_migration_tk.metadata_ncbi_manager import NCBIMeta

# the columns of an assembly summary as NCBI publishes it today
SUMMARY_HEADER = ['assembly_accession', 'bioproject', 'biosample', 'wgs_master', 'refseq_category',
                  'taxid', 'species_taxid', 'organism_name', 'infraspecific_name', 'isolate',
                  'version_status', 'assembly_level', 'release_type', 'genome_rep', 'seq_rel_date',
                  'asm_name', 'asm_submitter', 'gbrs_paired_asm', 'paired_asm_comp', 'ftp_path',
                  'excluded_from_refseq', 'relation_to_type_material', 'asm_not_live_date']

TABLE_HEADER = ['genome_id', 'ncbi_bioproject', 'ncbi_wgs_master', 'ncbi_wgs_formatted',
                'ncbi_refseq_category', 'ncbi_species_taxid', 'ncbi_isolate', 'ncbi_version_status',
                'ncbi_seq_rel_date', 'ncbi_asm_name', 'ncbi_gbrs_paired_asm', 'ncbi_paired_asm_comp',
                'ncbi_excluded_from_refseq', 'ncbi_not_used_as_type', 'ncbi_type_material_designation']


def summary_row(accession, wgs_master='JBAFXE000000000.1', excluded='na'):
    values = {'assembly_accession': accession, 'bioproject': 'PRJNA224116', 'biosample': 'SAMN1',
              'wgs_master': wgs_master, 'refseq_category': 'na', 'taxid': '7', 'species_taxid': '7',
              'organism_name': 'Azorhizobium caulinodans', 'infraspecific_name': 'strain=ORS 571',
              'isolate': 'na', 'version_status': 'latest', 'assembly_level': 'Complete Genome',
              'release_type': 'Major', 'genome_rep': 'Full', 'seq_rel_date': '2024-02-14',
              'asm_name': 'ASM3660089v1', 'asm_submitter': 'NCBI', 'gbrs_paired_asm': 'GCA_036600895.1',
              'paired_asm_comp': 'identical', 'ftp_path': 'https://ftp.ncbi.nlm.nih.gov/genomes/all/x',
              'excluded_from_refseq': excluded, 'relation_to_type_material': 'assembly from type material',
              'asm_not_live_date': 'na'}
    return values


class TempDirCase(unittest.TestCase):
    """A release's assembly summaries and genome list, in a directory of the test's own."""

    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='metadata_ncbi_manager_test.')
        self.addCleanup(shutil.rmtree, self.dir, True)
        self.warnings = []
        logger = logging.getLogger('timestamp')
        warnings = logging.Handler(level=logging.WARNING)
        warnings.emit = lambda record: self.warnings.append(record.getMessage())
        logger.addHandler(warnings)
        self.addCleanup(logger.removeHandler, warnings)

    def summary(self, name, rows, header=SUMMARY_HEADER, compress=False):
        path = os.path.join(self.dir, name)
        lines = ['##  See ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt',
                 '#' + '\t'.join(header)]
        lines += ['\t'.join(row[column] for column in header) for row in rows]
        with (gzip.open(path, 'wt') if compress else open(path, 'w')) as handle:
            handle.write('\n'.join(lines) + '\n')
        return path

    def genome_list(self, genomes):
        path = os.path.join(self.dir, 'metadata.tsv')
        with open(path, 'w') as handle:
            handle.write('accession\tother\n')
            handle.write(''.join('{}\tx\n'.format(gid) for gid in genomes))
        return path

    def parse(self, summaries, genomes):
        output = os.path.join(self.dir, 'ncbi_assembly_metadata.tsv')
        NCBIMeta().parse_assemblies(summaries, self.genome_list(genomes), output)
        with open(output) as handle:
            return [line.split('\t') for line in handle.read().splitlines()]


class ParsingTheSummaries(TempDirCase):
    def test_every_summary_given_is_read_however_many_there_are(self):
        refseq = self.summary('assembly_summary_bacteria_refseq.txt', [summary_row('GCF_000000001.1')])
        genbank = self.summary('assembly_summary_bacteria_genbank.txt', [summary_row('GCA_000000002.1')])
        viral = self.summary('assembly_summary_viral_genbank.txt', [summary_row('GCA_000000003.1')])

        table = self.parse([refseq, genbank, viral], ['GCF_000000001.1', 'GB_GCA_000000002.1', 'GCA_000000003.1'])

        self.assertEqual(table[0], TABLE_HEADER)
        self.assertEqual([row[0] for row in table[1:]],
                         ['RS_GCF_000000001.1', 'GB_GCA_000000002.1', 'GB_GCA_000000003.1'])

    def test_columns_are_found_by_name_in_each_summary_whatever_their_order(self):
        # the header was read from the first summary and the values of every
        # other summary written in that summary's own order beneath it
        reordered = list(reversed(SUMMARY_HEADER))
        first = self.summary('assembly_summary_bacteria_refseq.txt', [summary_row('GCF_000000001.1')])
        second = self.summary('assembly_summary_archaea_refseq.txt', [summary_row('GCF_000000002.1')],
                              header=reordered)

        table = self.parse([first, second], ['GCF_000000001.1', 'GCF_000000002.1'])

        self.assertEqual(table[2][1:], table[1][1:])

    def test_a_gzipped_summary_is_read_as_a_plain_one_is(self):
        rows = [summary_row('GCF_000000001.1')]
        plain = self.parse([self.summary('assembly_summary_bacteria_refseq.txt', rows)], ['GCF_000000001.1'])
        gzipped = self.parse([self.summary('assembly_summary_bacteria_refseq.txt.gz', rows, compress=True)],
                             ['GCF_000000001.1'])

        self.assertEqual(gzipped, plain)

    def test_the_wgs_master_and_refseq_exclusion_are_written_as_before(self):
        rows = [summary_row('GCF_000000001.1'),
                summary_row('GCA_000000002.1', wgs_master='na',
                            excluded='derived from surveillance project; not used as type')]
        table = self.parse([self.summary('assembly_summary_bacteria_refseq.txt', rows)],
                           ['GCF_000000001.1', 'GCA_000000002.1'])
        fields = [dict(zip(TABLE_HEADER, row)) for row in table[1:]]

        self.assertEqual((fields[0]['ncbi_wgs_master'], fields[0]['ncbi_wgs_formatted']),
                         ('JBAFXE000000000.1', 'JBAFXE01'))
        self.assertEqual((fields[0]['ncbi_excluded_from_refseq'], fields[0]['ncbi_not_used_as_type']),
                         ('', 'False'))
        self.assertEqual((fields[1]['ncbi_wgs_master'], fields[1]['ncbi_wgs_formatted']), ('na', ''))
        self.assertEqual(fields[1]['ncbi_not_used_as_type'], 'True')

    def test_only_the_genomes_listed_are_written_and_a_listed_genome_in_no_summary_is_warned_of(self):
        summary = self.summary('assembly_summary_bacteria_refseq.txt',
                               [summary_row('GCF_000000001.1'), summary_row('GCF_000000002.1')])

        table = self.parse([summary], ['GCF_000000001.1', 'GCF_000000009.1', 'U_12345'])

        self.assertEqual([row[0] for row in table[1:]], ['RS_GCF_000000001.1'])
        self.assertEqual(len(self.warnings), 1)
        self.assertIn('RS_GCF_000000009.1', self.warnings[0])
        self.assertNotIn('U_12345', self.warnings[0])


class TheCommandLine(TempDirCase):
    def test_parse_assemblies_takes_the_summaries_as_n(self):
        summaries = [self.summary('assembly_summary_bacteria_refseq.txt', [summary_row('GCF_000000001.1')]),
                     self.summary('assembly_summary_bacteria_genbank.txt', [summary_row('GCA_000000002.1')])]
        genomes = self.genome_list(['GCF_000000001.1'])
        output = os.path.join(self.dir, 'out.tsv')
        options = main_module.get_main_parser().parse_args(
            ['parse_assemblies', '-n'] + summaries + ['-m', genomes, '-o', output])

        with mock.patch.object(main_py, 'NCBIMeta') as meta:
            main_py.OptionsParser().parse_options(options)
        meta.return_value.parse_assemblies.assert_called_once_with(summaries, genomes, output)

    def test_the_log_is_written_where_l_says(self):
        summary = self.summary('assembly_summary_bacteria_refseq.txt', [summary_row('GCF_000000001.1')])
        log = os.path.join(self.dir, 'parse_assemblies.log')
        options = main_module.get_main_parser().parse_args(
            ['parse_assemblies', '-n', summary, '-m', self.genome_list([]),
             '-o', os.path.join(self.dir, 'out.tsv'), '-l', log])

        self.assertEqual(options.log, log)
        self.assertEqual(main_module.log_candidates(options.log, getattr(options, 'output_dir', None))[0],
                         (self.dir, 'parse_assemblies.log'))

    def test_the_four_summary_arguments_are_no_longer_accepted(self):
        summary = self.summary('assembly_summary_bacteria_refseq.txt', [summary_row('GCF_000000001.1')])
        with mock.patch('sys.stderr'), self.assertRaises(SystemExit):
            main_module.get_main_parser().parse_args(
                ['parse_assemblies', '--rb', summary, '--ra', summary, '--gb', summary, '--ga', summary,
                 '-m', self.genome_list([]), '-o', os.path.join(self.dir, 'out.tsv')])


if __name__ == '__main__':
    unittest.main()
