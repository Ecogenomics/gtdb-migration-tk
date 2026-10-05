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

"""strains type_table, run whole over a small release built from NCBI's files --
genome_dirs.tsv, the assembly summaries, names.dmp and nodes.dmp -- with no
database: which genome is type material, and that a genome the matching fails
on fails the command rather than being written as not type material.
"""

import contextlib
import gzip
import io
import logging
import os
import shutil
import tempfile
import unittest
from unittest import mock

from gtdb_migration_tk import strains as S

TYPE_ACCESSION, OTHER_ACCESSION = 'GCF_000000001.1', 'GCA_000000002.1'
TYPE_GENOME, OTHER_GENOME = 'RS_' + TYPE_ACCESSION, 'GB_' + OTHER_ACCESSION

SUMMARY_HEADER = ('#   See ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt\n'
                  '#assembly_accession\tbioproject\ttaxid\tspecies_taxid\torganism_name\t'
                  'infraspecific_name\tisolate\trelation_to_type_material\n')

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

        # accession -> (summary, taxid, organism name, infraspecific name,
        # isolate, relation to type material)
        self.genomes = {
            TYPE_ACCESSION: ('rb', '562', 'Escherichia coli ATCC 11775', 'strain=ATCC 11775',
                             'na', 'assembly from type material'),
            OTHER_ACCESSION: ('gb', '562', 'Escherichia coli K-12', 'strain=K-12 substr. MG1655',
                              'na', 'na'),
        }
        self.names_rows = ['562\t|\tEscherichia coli\t|\t\t|\tscientific name\t|',
                           '561\t|\tEscherichia\t|\t\t|\tscientific name\t|',
                           '562\t|\tATCC 11775\t|\t\t|\ttype material\t|',
                           '562\t|\tDSM 30083\t|\t\t|\ttype material\t|']
        self.nodes_rows = ['1\t|\t1\t|\tno rank\t|',
                           '561\t|\t1\t|\tgenus\t|',
                           '562\t|\t561\t|\tspecies\t|']

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

    def release_files(self, release=None):
        """genome_dirs.tsv of the release, and the four gzipped assembly summaries."""
        release = self.genomes if release is None else release
        genome_dirs = self.write('genome_dirs.tsv', ''.join(
            '{0}\t/release/{0}_ASM1v1\tG{1}\n'.format(acc, acc[4:13]) for acc in release))
        summaries = []
        for name in ('rb', 'ra', 'gb', 'ga'):
            path = os.path.join(self.dir, 'assembly_summary_{}.txt.gz'.format(name))
            with gzip.open(path, 'wt') as handle:
                handle.write(SUMMARY_HEADER)
                for acc, (summary, taxid, organism, infraspecific, isolate, type_material) in self.genomes.items():
                    if summary == name:
                        handle.write('\t'.join((acc, 'PRJNA1', taxid, taxid, organism, infraspecific,
                                                isolate, type_material)) + '\n')
            summaries.append(path)
        names = self.write('names.dmp', '\n'.join(self.names_rows) + '\n')
        nodes = self.write('nodes.dmp', '\n'.join(self.nodes_rows) + '\n')
        return genome_dirs, summaries, names, nodes

    def run_type_table(self, cpus=1):
        genome_dirs, summaries, names, nodes = self.release_files()
        S.Strains(self.out, cpus).generate_type_strain_table(
            genome_dirs, summaries, names, nodes, self.gss, self.lpsn_dir, self.years)

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

    def test_every_genome_of_the_release_is_in_the_summary_and_no_other(self):
        # a genome of the summaries the release does not hold is not decided
        self.genomes['GCA_000000009.1'] = ('gb', '562', 'Escherichia coli X', 'strain=X', 'na', 'na')
        genome_dirs, summaries, names, nodes = self.release_files(
            release=[TYPE_ACCESSION, OTHER_ACCESSION])
        S.Strains(self.out, 2).generate_type_strain_table(
            genome_dirs, summaries, names, nodes, self.gss, self.lpsn_dir, self.years)

        self.assertEqual(sorted(self.summary()), sorted([TYPE_GENOME, OTHER_GENOME]))

    def test_nothing_is_printed_to_the_console(self):
        # what the run has to say goes to the log, where --silent governs it
        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout):
            self.run_type_table(cpus=2)

        self.assertEqual(stdout.getvalue(), '')


class ReadingNcbisFiles(StrainsCase):
    """What type_table knows of a genome, read from NCBI's files rather than the database."""

    def genome(self):
        genome_dirs, summaries, names, nodes = self.release_files()
        strains = S.Strains(self.out)
        strains.metadata, taxids = strains.load_genomes(genome_dirs, summaries)
        rtn = strains.parse_ncbi_names_and_nodes(names, nodes, taxids)
        strains.name_species(rtn[-1])
        return strains

    def test_a_genome_is_named_by_its_accession_as_gtdb_writes_it(self):
        self.assertEqual(sorted(self.genome().metadata), sorted([TYPE_GENOME, OTHER_GENOME]))

    def test_strain_ids_are_split_from_the_infraspecific_name_and_a_substrain_dropped(self):
        metadata = self.genome().metadata

        self.assertEqual(metadata[TYPE_GENOME]['ncbi_standardised_strain_ids'], {'ATCC11775'})
        self.assertEqual(metadata[OTHER_GENOME]['ncbi_strain_ids'], 'K-12')

    def test_a_strain_of_n_a_is_no_strain(self):
        # it was split on its '/' into the strain IDs n and a
        self.genomes[OTHER_ACCESSION] = ('gb', '562', 'Escherichia coli', 'strain=n/a', 'H08', 'na')

        self.assertEqual(self.genome().metadata[OTHER_GENOME]['ncbi_strain_ids'], 'H08')

    def test_the_type_material_status_is_ncbis_null_where_it_gives_none(self):
        self.assertEqual(self.genome().metadata[OTHER_GENOME]['ncbi_type_material_designation'], 'na')

    def test_the_species_is_named_from_names_dmp_and_nodes_dmp(self):
        self.assertEqual(self.genome().get_species_name(TYPE_GENOME), 'Escherichia coli')

    def test_a_subspecies_is_named_where_ncbi_has_one(self):
        self.genomes[OTHER_ACCESSION] = ('gb', '563', 'Escherichia coli subsp. x', 'na', 'na', 'na')
        self.names_rows.append('563\t|\tEscherichia coli subsp. x\t|\t\t|\tscientific name\t|')
        self.nodes_rows.append('563\t|\t562\t|\tsubspecies\t|')

        self.assertEqual(self.genome().get_species_name(OTHER_GENOME), 'Escherichia coli subsp. x')

    def test_a_seqcode_name_is_named_without_its_code(self):
        # as ncbi_tax_manager names it
        self.names_rows[0] = '562\t|\tEscherichia coli (SeqCode)\t|\t\t|\tscientific name\t|'

        self.assertEqual(self.genome().get_species_name(TYPE_GENOME), 'Escherichia coli')

    def test_a_genome_whose_taxid_ncbi_deleted_is_not_type_material(self):
        # 1,485 genomes of r237, nearly all '<genus> sp.' placeholders
        self.genomes[OTHER_ACCESSION] = ('gb', '999999', 'Escherichia sp.', 'strain=ATCC 11775', 'na', 'na')
        self.assertIsNone(self.genome().get_species_name(OTHER_GENOME))

        self.run_type_table()
        row = self.summary()[OTHER_GENOME]
        self.assertEqual(row['ncbi_species'], '')
        self.assertEqual(row['gtdb_type_designation_ncbi_taxa'], 'not type material')

    def test_a_genome_of_the_release_in_no_summary_is_refused(self):
        genome_dirs, summaries, _, _ = self.release_files(
            release=[TYPE_ACCESSION, OTHER_ACCESSION, 'GCF_000000099.1'])

        with self.assertRaises(S.StrainsError) as caught:
            S.Strains(self.out).load_genomes(genome_dirs, summaries)
        self.assertIn('GCF_000000099.1', str(caught.exception))


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

    WEB_ONLY_ACCESSION = 'GCA_000000003.1'
    WEB_ONLY_GENOME = 'GB_' + WEB_ONLY_ACCESSION

    def lpsn_strains(self, *rows):
        self.write(os.path.join('lpsn', 'lpsn_strains.tsv'),
                   'lpsn_strain\tco-identical strain IDs\ttype_designation\n' +
                   ''.join('\t'.join(row) + '\n' for row in rows))

    def add_web_only_species(self):
        """Escherichia albertii: a genome of a species on LPSN's web pages and not in the GSS file."""
        self.genomes[self.WEB_ONLY_ACCESSION] = ('gb', '208962', 'Escherichia albertii LMG 20976',
                                                 'strain=LMG 20976', 'na', 'na')
        self.names_rows.append('208962\t|\tEscherichia albertii\t|\t\t|\tscientific name\t|')
        self.nodes_rows.append('208962\t|\t561\t|\tspecies\t|')

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


class ChoosingBetweenNames(StrainsCase):
    """A genome matching LPSN under two of its NCBI names, to the same type strain.

    The names came from a set, and the first in Python's order won, so a run
    could report either: RS_GCF_964245115.1 as Clostridium ramosum (a type strain)
    on one run and Erysipelatoclostridium ramosum (a nomenclatural type) on the next.
    """

    RAMOSA_ACCESSION = 'GCF_000000005.1'
    RAMOSA_GENOME = 'RS_' + RAMOSA_ACCESSION

    def ramosa(self, *lpsn_rows):
        """Thomasclavelia ramosa, known to NCBI also by two names LPSN holds."""
        self.genomes[self.RAMOSA_ACCESSION] = ('rb', '1547', 'Thomasclavelia ramosa',
                                               'strain=ATCC 25582', 'na', 'na')
        self.names_rows += ['1547\t|\tThomasclavelia ramosa\t|\t\t|\tscientific name\t|',
                            '1547\t|\tAclostridium ramosum\t|\t\t|\tequivalent name\t|',
                            '1547\t|\tClostridium ramosum\t|\t\t|\tequivalent name\t|']
        self.nodes_rows.append('1547\t|\t1\t|\tspecies\t|')
        self.write(os.path.join('lpsn', 'lpsn_strains.tsv'),
                   'lpsn_strain\tco-identical strain IDs\ttype_designation\n'
                   'Escherichia coli\tATCC11775=DSM30083\tType strain\n' +
                   ''.join('\t'.join(row) + '\n' for row in lpsn_rows))

    def lpsn_match(self):
        with open(os.path.join(self.out, 'lpsn_summary.tsv')) as handle:
            header = handle.readline().rstrip('\n').split('\t')
            for line in handle:
                row = dict(zip(header, line.rstrip('\n').split('\t')))
                if row['genome'] == self.RAMOSA_GENOME:
                    return row['lpsn_match_name'], row['gtdb_type_designation']

    def test_a_type_strain_outranks_a_nomenclatural_type_whichever_name_comes_first(self):
        # the nomenclatural type's name sorts first, so it is the ranking that decides
        self.ramosa(('Aclostridium ramosum', 'ATCC25582', 'Nomenclatural type'),
                    ('Clostridium ramosum', 'ATCC25582', 'Type strain'))
        self.run_type_table(cpus=2)

        self.assertEqual(self.lpsn_match(), ('Clostridium ramosum', 'type strain of species'))
        self.assertEqual(self.summary()[self.RAMOSA_GENOME]['lpsn_type_designation'],
                         'type strain of species')

    def test_two_names_as_good_as_each_other_give_the_same_one_on_every_run(self):
        self.ramosa(('Aclostridium ramosum', 'ATCC25582', 'Type strain'),
                    ('Clostridium ramosum', 'ATCC25582', 'Type strain'))
        self.run_type_table(cpus=2)

        self.assertEqual(self.lpsn_match(), ('Aclostridium ramosum', 'type strain of species'))


if __name__ == '__main__':
    unittest.main()
