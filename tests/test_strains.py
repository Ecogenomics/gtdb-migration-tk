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

"""strains type_table, run whole over a release of three genomes: which is type
material, and that a genome the matching fails on fails the command rather than
being written as not type material.
"""

import contextlib
import io
import logging
import os
import shutil
import tempfile
import unittest
from unittest import mock

from gtdb_migration_tk import strains as S

TYPE_GENOME = 'RS_GCF_000000001.1'
OTHER_GENOME = 'GB_GCA_000000002.1'
USER_GENOME = 'U_000001'
TAXONOMY = 'd__Bacteria;p__Pseudomonadota;g__Escherichia;s__Escherichia coli'

METADATA_COLUMNS = ('accession', 'ncbi_organism_name', 'ncbi_type_material_designation',
                    'ncbi_strain_identifiers', 'ncbi_taxonomy_unfiltered', 'ncbi_taxid')
METADATA_ROWS = (
    (TYPE_GENOME, 'Escherichia coli ATCC 11775', 'assembly from type material',
     'ATCC 11775; DSM 30083', TAXONOMY, '562'),
    (OTHER_GENOME, 'Escherichia coli K-12', 'none', 'K-12', TAXONOMY, '562'),
    (USER_GENOME, 'a user genome', 'none', 'none', 'none', 'none'),
)

GSS_HEADER = ('genus_name,sp_epithet,subsp_epithet,reference,status,authors,address,'
              'risk_grp,nomenclatural_type,record_no,record_lnk\n')
GSS_ROW = ('Escherichia,coli,,ref,"validly published under the ICNP;correct name",'
           'Castellani and Chalmers 1919,,2,ATCC 11775; DSM 30083,1,link\n')


class StrainsCase(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='strains_test.')
        self.addCleanup(shutil.rmtree, self.dir, True)
        self.out = os.path.join(self.dir, 'out')
        os.makedirs(self.out)
        logging.getLogger('timestamp').addHandler(logging.NullHandler())

        self.metadata = self.write('metadata.tsv', '\t'.join(METADATA_COLUMNS) + '\n' +
                                   ''.join('\t'.join(row) + '\n' for row in METADATA_ROWS))
        self.names = self.write('names.dmp',
                                '562\t|\tEscherichia coli\t|\t\t|\tscientific name\t|\n'
                                '561\t|\tEscherichia\t|\t\t|\tscientific name\t|\n')
        self.nodes = self.write('nodes.dmp',
                                '1\t|\t1\t|\tno rank\t|\n'
                                '561\t|\t1\t|\tgenus\t|\n'
                                '562\t|\t561\t|\tspecies\t|\n')
        self.lpsn_dir = os.path.join(self.dir, 'lpsn')
        os.makedirs(self.lpsn_dir)
        self.write(os.path.join('lpsn', 'lpsn_strains.tsv'),
                   'lpsn_strain\tco-identical strain IDs\ttype_designation\n'
                   'Escherichia coli\tATCC11775=DSM30083\tType strain\n')
        self.write(os.path.join('lpsn', 'lpsn_species.tsv'),
                   'lpsn_species\tlpsn_type_species\tlpsn_species_authority\tsource\n'
                   's__Escherichia coli\tg__Escherichia\tCastellani and Chalmers 1919\tGSS\n')
        self.gss = self.write('lpsn_gss.csv', GSS_HEADER + GSS_ROW)
        self.years = self.write('years.tsv', 'Escherichia coli\t1919\n')

    def write(self, name, text):
        path = os.path.join(self.dir, name)
        with open(path, 'w') as handle:
            handle.write(text)
        return path

    def run_type_table(self, cpus=1):
        S.Strains(self.out, cpus).generate_type_strain_table(
            self.metadata, self.names, self.nodes, self.gss, self.lpsn_dir, self.years)

    def summary(self):
        with open(os.path.join(self.out, 'gtdb_type_strain_summary.tsv')) as handle:
            header = handle.readline().rstrip('\n').split('\t')
            return {row[0]: dict(zip(header, row)) for row in
                    (line.rstrip('\n').split('\t') for line in handle)}


class DecidingTypeMaterial(StrainsCase):
    def test_the_type_strain_is_type_material_and_another_strain_is_not(self):
        for cpus in (1, 2):
            with self.subTest(cpus=cpus):
                self.run_type_table(cpus)
                summary = self.summary()

                self.assertEqual(summary[TYPE_GENOME]['gtdb_type_designation_ncbi_taxa'],
                                 'type strain of species')
                self.assertEqual(summary[TYPE_GENOME]['gtdb_type_designation_ncbi_taxa_sources'], 'LPSN')
                self.assertEqual(summary[TYPE_GENOME]['gtdb_type_species_of_genus'], 'True')
                self.assertEqual(summary[TYPE_GENOME]['lpsn_priority_year'], '1919')
                self.assertEqual(summary[OTHER_GENOME]['gtdb_type_designation_ncbi_taxa'],
                                 'not type material')

    def test_every_ncbi_genome_is_in_the_summary_and_no_user_genome(self):
        self.run_type_table(cpus=2)

        self.assertEqual(sorted(self.summary()), sorted([TYPE_GENOME, OTHER_GENOME]))

    def test_nothing_is_printed_to_the_console(self):
        # what the run has to say goes to the log, where --silent governs it
        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout):
            self.run_type_table(cpus=2)

        self.assertEqual(stdout.getvalue(), '')


class AFailureFailsTheCommand(StrainsCase):
    """A genome the matching fails on was dropped, then written as not type material."""

    def failing_on(self, gid):
        real = S.Strains.match_genome

        def match_genome(strains, genome):
            if genome == gid:
                raise S.StrainsError('cannot decide {}'.format(genome))
            return real(strains, genome)
        return mock.patch.object(S.Strains, 'match_genome', match_genome)

    def test_a_worker_that_fails_fails_the_command(self):
        for cpus in (1, 2):
            with self.subTest(cpus=cpus), self.failing_on(TYPE_GENOME):
                with self.assertRaises(S.StrainsError):
                    self.run_type_table(cpus)

    def test_a_failed_run_leaves_no_table_that_reads_as_finished(self):
        with self.failing_on(TYPE_GENOME), self.assertRaises(S.StrainsError):
            self.run_type_table(cpus=2)

        self.assertEqual(os.listdir(self.out), [])

    def test_a_contradiction_raises_rather_than_exiting_the_process(self):
        # sys.exit() in a worker ended the worker and nothing else
        with self.assertRaises(S.StrainsError):
            S.Strains().select_category_name('Escherichia coli', {}, {}, {})


class WhatLpsnCallsTheType(StrainsCase):
    """The third column of lpsn_strains.tsv, and the GSS file's say over it."""

    WEB_ONLY_GENOME = 'GB_GCA_000000003.1'

    def lpsn_strains(self, *rows):
        self.write(os.path.join('lpsn', 'lpsn_strains.tsv'),
                   'lpsn_strain\tco-identical strain IDs\ttype_designation\n' +
                   ''.join('\t'.join(row) + '\n' for row in rows))

    def add_web_only_species(self):
        """Escherichia albertii: a genome of a species on LPSN's web pages and not in the GSS file."""
        with open(self.metadata, 'a') as handle:
            handle.write('\t'.join((self.WEB_ONLY_GENOME, 'Escherichia albertii LMG 20976', 'none',
                                    'LMG 20976', TAXONOMY.replace('coli', 'albertii'), '208962')) + '\n')
        with open(self.names, 'a') as handle:
            handle.write('208962\t|\tEscherichia albertii\t|\t\t|\tscientific name\t|\n')
        with open(self.nodes, 'a') as handle:
            handle.write('208962\t|\t561\t|\tspecies\t|\n')

    def test_a_combined_designation_is_refused_rather_than_guessed_at(self):
        self.lpsn_strains(('Escherichia coli', 'ATCC11775=DSM30083', 'Type strain;Holotype'))

        with self.assertRaises(S.StrainsError) as caught:
            self.run_type_table()
        self.assertIn('Escherichia coli (Type strain;Holotype)', str(caught.exception))

    def test_a_designation_it_does_not_know_is_refused(self):
        self.lpsn_strains(('Escherichia coli', 'ATCC11775=DSM30083', 'Neotype'))

        with self.assertRaises(S.StrainsError):
            self.run_type_table()

    def test_a_gss_species_is_a_type_strain_whatever_its_web_page_says(self):
        # as 152 species of r232 are: the web page parsed was that of another LPSN
        # record of the name, not validly published
        self.lpsn_strains(('Escherichia coli', 'ATCC11775=DSM30083', 'Nomenclatural type'))
        self.run_type_table()

        self.assertEqual(self.summary()[TYPE_GENOME]['lpsn_type_designation'], 'type strain of species')

    def test_a_species_only_on_the_web_pages_keeps_their_designation(self):
        self.add_web_only_species()
        self.lpsn_strains(('Escherichia coli', 'ATCC11775=DSM30083', 'Type strain'),
                          ('Escherichia albertii', 'LMG20976', 'Nomenclatural type'))
        self.run_type_table()

        row = self.summary()[self.WEB_ONLY_GENOME]
        self.assertEqual(row['lpsn_type_designation'], 'nomenclatural type of species')
        # and GTDB's own designation folds it into a type strain of species
        self.assertEqual(row['gtdb_type_designation_ncbi_taxa'], 'type strain of species')


class ReadingTheMetadata(StrainsCase):
    def test_a_quoted_comma_in_a_comma_separated_export_stays_in_its_field(self):
        csv_file = self.write('metadata.csv', ','.join(METADATA_COLUMNS) + '\n' +
                              '{},"Escherichia coli, strain ATCC 11775",none,ATCC 11775,{},562\n'.format(
                                  TYPE_GENOME, TAXONOMY))

        metadata, taxids = S.Strains().load_metadata(csv_file)

        self.assertEqual(metadata[TYPE_GENOME]['ncbi_organism_name'],
                         'Escherichia coli, strain ATCC 11775')
        self.assertEqual(metadata[TYPE_GENOME]['ncbi_taxid'], 562)
        self.assertEqual(taxids, {562})

    def test_a_tab_separated_export_is_read_as_before(self):
        metadata, _ = S.Strains().load_metadata(self.metadata)

        self.assertEqual(metadata[TYPE_GENOME]['ncbi_standardised_strain_ids'],
                         {'ATCC11775', 'DSM30083'})
        self.assertNotIn(USER_GENOME, metadata)


if __name__ == '__main__':
    unittest.main()
