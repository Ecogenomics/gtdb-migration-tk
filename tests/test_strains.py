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
                  'infraspecific_name\tisolate\trelation_to_type_material\texcluded_from_refseq\n')

# the columns of download_seqcode_data's table type_table reads, and one it does not
SEQCODE_HEADER = 'seqcode_type_material_accn\tseqcode_id\tseqcode_species_status\tseqcode_type_species_of_genus\n'
METAGENOME_NOT_TYPE = 'derived from metagenome; not used as type'

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
        # isolate, relation to type material[, excluded from RefSeq])
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
        # (accession, status, type species of genus) of seqcode_table.tsv
        self.seqcode_rows = []

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
                for acc, genome in self.genomes.items():
                    summary, taxid, organism, infraspecific, isolate, type_material = genome[:6]
                    excluded = genome[6] if len(genome) > 6 else 'na'
                    if summary == name:
                        handle.write('\t'.join((acc, 'PRJNA1', taxid, taxid, organism, infraspecific,
                                                isolate, type_material, excluded)) + '\n')
            summaries.append(path)
        names = self.write('names.dmp', '\n'.join(self.names_rows) + '\n')
        nodes = self.write('nodes.dmp', '\n'.join(self.nodes_rows) + '\n')
        return genome_dirs, summaries, names, nodes

    def seqcode_table(self):
        """download_seqcode_data's seqcode_table.tsv, of self.seqcode_rows."""
        return self.write('seqcode_table.tsv', SEQCODE_HEADER + ''.join(
            '{}\t{}\t{}\t{}\n'.format(acc, number, status, type_species)
            for number, (acc, status, type_species) in enumerate(self.seqcode_rows, start=100)))

    def run_type_table(self, cpus=1):
        genome_dirs, summaries, names, nodes = self.release_files()
        S.Strains(self.out, cpus).generate_type_strain_table(
            genome_dirs, summaries, names, nodes, self.gss, self.lpsn_dir, self.years, self.seqcode_table())

    def summary(self):
        with gzip.open(os.path.join(self.out, S.TYPE_STRAIN_SUMMARY_NAME), 'rt') as handle:
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

    def test_every_row_of_the_summary_has_a_value_for_each_column_of_its_header(self):
        # update_metadata_db refuses a table with a short row, and the header
        # named an is_from_standard column no row held
        self.run_type_table()
        with gzip.open(os.path.join(self.out, S.TYPE_STRAIN_SUMMARY_NAME), 'rt') as handle:
            header = handle.readline().rstrip('\n').split('\t')
            rows = [line.rstrip('\n').split('\t') for line in handle]

        self.assertTrue(rows)
        for row in rows:
            self.assertEqual(len(row), len(header), row[0])

    def test_every_genome_of_the_release_is_in_the_summary_and_no_other(self):
        # a genome of the summaries the release does not hold is not decided
        self.genomes['GCA_000000009.1'] = ('gb', '562', 'Escherichia coli X', 'strain=X', 'na', 'na')
        genome_dirs, summaries, names, nodes = self.release_files(
            release=[TYPE_ACCESSION, OTHER_ACCESSION])
        S.Strains(self.out, 2).generate_type_strain_table(
            genome_dirs, summaries, names, nodes, self.gss, self.lpsn_dir, self.years, self.seqcode_table())

        self.assertEqual(sorted(self.summary()), sorted([TYPE_GENOME, OTHER_GENOME]))

    def test_the_summary_is_gzipped(self):
        self.run_type_table()

        self.assertIn(S.TYPE_STRAIN_SUMMARY_NAME, os.listdir(self.out))
        self.assertNotIn('gtdb_type_strain_summary.tsv', os.listdir(self.out))
        with open(os.path.join(self.out, S.TYPE_STRAIN_SUMMARY_NAME), 'rb') as handle:
            self.assertEqual(handle.read(2), b'\x1f\x8b')

    def test_the_same_inputs_write_the_same_bytes(self):
        # the gzip header records no time and no file name
        self.run_type_table()
        with open(os.path.join(self.out, S.TYPE_STRAIN_SUMMARY_NAME), 'rb') as handle:
            first = handle.read()
        self.out = os.path.join(self.dir, 'again')
        os.makedirs(self.out)
        self.run_type_table()
        with open(os.path.join(self.out, S.TYPE_STRAIN_SUMMARY_NAME), 'rb') as handle:
            self.assertEqual(handle.read(), first)

    def test_nothing_is_printed_to_the_console(self):
        # what the run has to say goes to the log, where --silent governs it
        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout):
            self.run_type_table(cpus=2)

        self.assertEqual(stdout.getvalue(), '')


class TypeMaterialUnderTheSeqCode(StrainsCase):
    """LPSN, then the SeqCode, then NCBI's exclusion of a metagenome not used as
    type: what update_type_designation did to metadata_type_material once the
    summary was loaded, decided where the summary is written."""

    def test_a_genome_typing_a_species_valid_under_the_seqcode_is_a_type_strain_of_species(self):
        self.seqcode_rows = [(OTHER_ACCESSION, 'Valid (SeqCode)', 'True')]
        self.run_type_table()
        other = self.summary()[OTHER_GENOME]
        self.assertEqual((other['gtdb_type_designation_ncbi_taxa'], other['gtdb_type_designation_ncbi_taxa_sources'],
                          other['gtdb_type_species_of_genus'], other['lpsn_type_designation']),
                         ('type strain of species', 'SeqCode', 'True', 'not type material'))

    def test_lpsn_and_the_seqcode_agreeing_are_both_sources_joined_by_one_separator(self):
        # update_type_designation appended ';Seqcode' to type_table's '; '-joined sources
        self.seqcode_rows = [(TYPE_ACCESSION, 'Valid (SeqCode)', 'False')]
        self.run_type_table()
        typed = self.summary()[TYPE_GENOME]
        self.assertEqual(typed['gtdb_type_designation_ncbi_taxa_sources'], 'LPSN; SeqCode')
        # LPSN's type species of genus stands where the SeqCode does not say so
        self.assertEqual(typed['gtdb_type_species_of_genus'], 'True')

    def test_a_species_not_valid_under_the_seqcode_types_nothing(self):
        self.seqcode_rows = [(OTHER_ACCESSION, '', 'True')]
        self.run_type_table()
        other = self.summary()[OTHER_GENOME]
        self.assertEqual((other['gtdb_type_designation_ncbi_taxa'], other['gtdb_type_species_of_genus']),
                         ('not type material', 'False'))

    def test_ncbis_metagenome_not_used_as_type_overrides_lpsn_and_the_seqcode_and_says_so(self):
        # update_type_designation left the sources and the type species of genus as they were
        self.genomes[TYPE_ACCESSION] += (METAGENOME_NOT_TYPE,)
        self.seqcode_rows = [(TYPE_ACCESSION, 'Valid (SeqCode)', 'True')]
        self.run_type_table()
        typed = self.summary()[TYPE_GENOME]
        self.assertEqual((typed['gtdb_type_designation_ncbi_taxa'], typed['gtdb_type_designation_ncbi_taxa_sources'],
                          typed['gtdb_type_species_of_genus'], typed['lpsn_type_designation']),
                         ('not used as type', '', 'False', 'type strain of species'))
        self.assertEqual(typed['gtdb_type_designation_notes'],
                         "NCBI's 'derived from metagenome; not used as type' rule overrides LPSN and SeqCode.")

    def test_the_ncbi_rule_needs_both_and_applies_to_a_genome_neither_types(self):
        self.genomes[OTHER_ACCESSION] += ('contaminated; ' + METAGENOME_NOT_TYPE,)
        self.genomes[TYPE_ACCESSION] += ('derived from metagenome',)
        self.run_type_table()
        summary = self.summary()
        self.assertEqual((summary[OTHER_GENOME]['gtdb_type_designation_ncbi_taxa'],
                          summary[OTHER_GENOME]['gtdb_type_designation_notes']),
                         ('not used as type', "NCBI's 'derived from metagenome; not used as type' rule applies."))
        self.assertEqual((summary[TYPE_GENOME]['gtdb_type_designation_ncbi_taxa'],
                          summary[TYPE_GENOME]['gtdb_type_designation_notes']), ('type strain of species', ''))

    def test_the_notes_are_the_last_column_and_every_row_has_one(self):
        # columns are appended, never reordered
        self.run_type_table()
        with gzip.open(os.path.join(self.out, S.TYPE_STRAIN_SUMMARY_NAME), 'rt') as handle:
            header = handle.readline().rstrip('\n').split('\t')
            rows = [line.rstrip('\n').split('\t') for line in handle]
        self.assertEqual(header[-2:], ['gtdb_type_species_of_genus', 'gtdb_type_designation_notes'])
        self.assertTrue(all(len(row) == len(header) for row in rows))

    def test_a_seqcode_genome_the_release_does_not_hold_is_warned_of_and_passed_over(self):
        self.seqcode_rows = [('GCA_000000099.1', 'Valid (SeqCode)', 'True')]
        with self.assertLogs('timestamp', level='WARNING') as logged:
            self.run_type_table()
        self.assertNotIn('GB_GCA_000000099.1', self.summary())
        self.assertTrue(any('1 genomes of the SeqCode table not in the release' in m for m in logged.output))

    def test_a_seqcode_table_without_its_columns_is_refused(self):
        genome_dirs, summaries, names, nodes = self.release_files()
        table = self.write('not_seqcode.tsv', 'genome\tstatus\n{}\tValid\n'.format(OTHER_ACCESSION))
        with self.assertRaisesRegex(S.StrainsError, 'seqcode_type_material_accn'):
            S.Strains(self.out).generate_type_strain_table(
                genome_dirs, summaries, names, nodes, self.gss, self.lpsn_dir, self.years, table)

    def test_a_summary_without_excluded_from_refseq_is_refused(self):
        # read as empty, it would pass over NCBI's exclusion of every genome unsaid
        from gtdb_migration_tk.ncbi_utils import BadInput
        genome_dirs, summaries, _, _ = self.release_files()
        with gzip.open(summaries[0], 'rt') as handle:
            text = handle.read().replace('\texcluded_from_refseq', '')
        with gzip.open(summaries[0], 'wt') as handle:
            handle.write(text)
        with self.assertRaisesRegex(BadInput, 'excluded_from_refseq'):
            S.Strains(self.out).load_genomes(genome_dirs, summaries)

    def test_the_columns_read_are_those_download_seqcode_data_writes(self):
        from gtdb_migration_tk import seqcode_manager
        for column in (S.SEQCODE_GENOME, S.SEQCODE_STATUS, S.SEQCODE_TYPE_SPECIES_OF_GENUS):
            self.assertIn(column, seqcode_manager.TABLE_HEADER)


class TheTypeTableCommandLine(StrainsCase):
    ARGV = ['strains', 'type_table', '-g', 'genome_dirs.tsv', '-n', 'a.txt', '--ncbi_names', 'names.dmp',
            '--ncbi_nodes', 'nodes.dmp', '--lpsn_gss_file', 'gss.csv', '--lpsn_dir', 'lpsn',
            '--year_table', 'years.tsv', '-o', 'out']

    def test_the_seqcode_table_is_required_and_handed_to_the_command(self):
        from gtdb_migration_tk import __main__ as main_module
        from gtdb_migration_tk import main as main_py
        with mock.patch('sys.stderr'), self.assertRaises(SystemExit) as ended:
            main_module.get_main_parser().parse_args(self.ARGV)
        self.assertEqual(ended.exception.code, 2)

        options = main_module.get_main_parser().parse_args(self.ARGV + ['--seqcode_table', 'seqcode_table.tsv'])
        with mock.patch.object(main_py, 'Strains') as strains, mock.patch.object(main_py, 'check_file_exists'):
            main_py.OptionsParser().parse_options(options)
        strains.return_value.generate_type_strain_table.assert_called_once_with(
            'genome_dirs.tsv', ['a.txt'], 'names.dmp', 'nodes.dmp', 'gss.csv', 'lpsn', 'years.tsv',
            'seqcode_table.tsv')


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


class TheWarnings(StrainsCase):
    """One WARNING per kind at the end, and every warning in type_table_warnings.tsv."""

    def setUp(self):
        super().setUp()
        self.records = []
        handler = logging.Handler(logging.WARNING)
        handler.emit = self.records.append
        logging.getLogger('timestamp').addHandler(handler)
        self.addCleanup(logging.getLogger('timestamp').removeHandler, handler)

    def warn_of_everything(self):
        # a strain under four spellings NCBI does not consider type
        self.names_rows += ['562\t|\t{0}\t|\t{0} <not considered type>\t|\ttype material\t|'.format(s)
                            for s in ('ATCC 23581', 'ATCC-23581', 'ATCC:23581', 'ATCC23581')]
        # a genome whose taxid NCBI deleted
        self.genomes['GCA_000000007.1'] = ('gb', '999999', 'Escherichia sp.', 'na', 'na', 'na')
        # a genome under a subspecies-rank strain name
        self.genomes['GCA_000000008.1'] = ('gb', '600', 'Escherichia coli X1', 'na', 'na', 'na')
        self.names_rows.append('600\t|\tEscherichia coli X1\t|\t\t|\tscientific name\t|')
        self.nodes_rows.append('600\t|\t562\t|\tsubspecies\t|')
        # a genus LPSN gives two type species
        self.write(os.path.join('lpsn', 'lpsn_species.tsv'),
                   'lpsn_species\tlpsn_type_species\tlpsn_species_authority\tsource\n'
                   's__Escherichia coli\tg__Escherichia\tCastellani and Chalmers 1919\tGSS\n'
                   's__Escherichia albertii\tg__Escherichia\tHuys et al. 2003\tGSS\n')

    def warnings_tsv(self):
        with open(os.path.join(self.out, S.WARNINGS_NAME)) as handle:
            header = handle.readline().rstrip('\n').split('\t')
            return header, [dict(zip(header, line.rstrip('\n').split('\t'))) for line in handle]

    def test_each_kind_of_warning_is_one_line_with_its_count_and_three_examples(self):
        self.warn_of_everything()
        self.run_type_table(cpus=2)

        lines = [r.getMessage() for r in self.records]
        self.assertEqual(len(lines), 4, lines)
        not_type = [l for l in lines if 'not considered type' in l][0]
        self.assertTrue(not_type.startswith('4 strain IDs names.dmp lists'), not_type)
        self.assertEqual(not_type.count('ATCC'), 3)
        self.assertIn(S.WARNINGS_NAME, not_type)
        self.assertTrue([l for l in lines if l.startswith('1 genomes with no NCBI species')])
        self.assertTrue([l for l in lines if 'Escherichia coli X1' in l])
        self.assertTrue([l for l in lines if l.startswith('1 genera LPSN gives')])

    def test_every_warning_is_in_the_tsv_with_its_meaning_and_its_data(self):
        self.warn_of_everything()
        self.run_type_table()

        header, rows = self.warnings_tsv()
        self.assertEqual(header, list(S.WARNINGS_HEADER))
        kinds = [row['warning_type'] for row in rows]
        self.assertEqual(kinds.count(S.NOT_CONSIDERED_TYPE), 4)
        for row in rows:
            self.assertEqual(row['description'], S.WARNING_KINDS[row['warning_type']][1])

        not_type = [r for r in rows if r['warning_type'] == S.NOT_CONSIDERED_TYPE][0]
        self.assertEqual(not_type['warning'], 'Ignoring ATCC 23581 as it is not considered type material.')
        self.assertEqual(not_type['extra_data'], 'taxid=562; unique_name=ATCC 23581 <not considered type>')
        unnamed = [r for r in rows if r['warning_type'] == S.NO_NCBI_SPECIES][0]
        self.assertIn('taxid=999999; in_nodes_dmp=no', unnamed['extra_data'])
        genus = [r for r in rows if r['warning_type'] == S.MULTIPLE_TYPE_SPECIES][0]
        self.assertIn('type_species=Escherichia coli, Escherichia albertii', genus['extra_data'])

    def test_a_subspecies_name_is_warned_of_once_a_genome(self):
        # it was logged wherever the genome's species name was asked for
        self.warn_of_everything()
        self.run_type_table()

        _, rows = self.warnings_tsv()
        self.assertEqual([r['warning_type'] for r in rows].count(S.SUBSPECIES_WITHOUT_SUBSP), 1)

    def test_a_run_with_nothing_to_warn_of_says_so(self):
        self.run_type_table()

        self.assertEqual(self.records, [])
        self.assertEqual(self.warnings_tsv(), (list(S.WARNINGS_HEADER), []))


if __name__ == '__main__':
    unittest.main()
