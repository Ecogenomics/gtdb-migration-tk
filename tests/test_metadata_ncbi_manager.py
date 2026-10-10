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

"""Offline unit tests for metadata_ncbi_manager.py -- parse_ncbi_assemblies, which
writes the NCBI metadata of every genome of the assembly summaries.

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
from gtdb_migration_tk import metadata_database_manager
from gtdb_migration_tk import metadata_ncbi_manager as N
from gtdb_migration_tk.metadata_ncbi_manager import NCBI_ASSEMBLY_TABLE, NCBIMeta
from gtdb_migration_tk.utils.common import open_text

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

    def parse(self, summaries):
        output = NCBIMeta().parse_assemblies(summaries, self.dir)
        self.assertEqual(output, os.path.join(self.dir, NCBI_ASSEMBLY_TABLE + '.gz'))
        with open_text(output) as handle:
            return [line.split('\t') for line in handle.read().splitlines()]


class ParsingTheSummaries(TempDirCase):
    def test_every_summary_given_is_read_however_many_there_are(self):
        refseq = self.summary('assembly_summary_bacteria_refseq.txt', [summary_row('GCF_000000001.1')])
        genbank = self.summary('assembly_summary_bacteria_genbank.txt', [summary_row('GCA_000000002.1')])
        viral = self.summary('assembly_summary_viral_genbank.txt', [summary_row('GCA_000000003.1')])

        table = self.parse([refseq, genbank, viral])

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

        table = self.parse([first, second])

        self.assertEqual(table[2][1:], table[1][1:])

    def test_a_gzipped_summary_is_read_as_a_plain_one_is(self):
        rows = [summary_row('GCF_000000001.1')]
        plain = self.parse([self.summary('assembly_summary_bacteria_refseq.txt', rows)])
        gzipped = self.parse([self.summary('assembly_summary_bacteria_refseq.txt.gz', rows, compress=True)])

        self.assertEqual(gzipped, plain)

    def test_the_wgs_master_and_refseq_exclusion_are_written_as_before(self):
        rows = [summary_row('GCF_000000001.1'),
                summary_row('GCA_000000002.1', wgs_master='na',
                            excluded='derived from surveillance project; not used as type')]
        table = self.parse([self.summary('assembly_summary_bacteria_refseq.txt', rows)])
        fields = [dict(zip(TABLE_HEADER, row)) for row in table[1:]]

        self.assertEqual((fields[0]['ncbi_wgs_master'], fields[0]['ncbi_wgs_formatted']),
                         ('JBAFXE000000000.1', 'JBAFXE01'))
        self.assertEqual((fields[0]['ncbi_excluded_from_refseq'], fields[0]['ncbi_not_used_as_type']),
                         ('', 'False'))
        self.assertEqual((fields[1]['ncbi_wgs_master'], fields[1]['ncbi_wgs_formatted']), ('na', ''))
        self.assertEqual(fields[1]['ncbi_not_used_as_type'], 'True')

    def test_every_genome_of_the_summaries_is_written_in_their_order(self):
        # which of them the database loads is update_metadata_db's to decide
        summary = self.summary('assembly_summary_bacteria_refseq.txt',
                               [summary_row('GCF_000000002.1'), summary_row('GCF_000000001.1')])

        table = self.parse([summary])

        self.assertEqual([row[0] for row in table[1:]], ['RS_GCF_000000002.1', 'RS_GCF_000000001.1'])
        self.assertEqual(self.warnings, [])

    def test_a_genome_in_two_summaries_is_written_once_and_warned_of(self):
        # a summary given twice gave each of its genomes two rows
        summary = self.summary('assembly_summary_bacteria_refseq.txt', [summary_row('GCF_000000001.1')])

        table = self.parse([summary, summary])

        self.assertEqual([row[0] for row in table[1:]], ['RS_GCF_000000001.1'])
        self.assertEqual(len(self.warnings), 1)
        self.assertIn('1 genome(s) are in more than one assembly summary', self.warnings[0])
        self.assertIn('RS_GCF_000000001.1 again in ' + summary, self.warnings[0])


class WhatUpdateMetadataDbLoads(TempDirCase):
    def test_the_table_is_named_as_update_metadata_db_loads_it_against_its_own_fields(self):
        # --input_folder refuses a table of any other name, and loads
        # ncbi_assembly_metadata.tsv against the fields of parse_ncbi_dir
        with mock.patch.object(metadata_database_manager.GenomeDatabaseConnectionFTPUpdate,
                               'GenomeDatabaseConnectionFTPUpdate'):
            manager = metadata_database_manager.MetadataDatabaseManager({})
        descriptions = manager.description_table[NCBI_ASSEMBLY_TABLE]
        self.assertEqual(len(descriptions), 1)
        path = os.path.join(os.path.dirname(metadata_database_manager.__file__),
                            'data_files', 'table_description', descriptions[0])
        with open(path) as handle:
            described = {line.split('\t')[0] for line in handle if line.strip()}

        self.assertEqual(set(TABLE_HEADER[1:]), described)


class TheCommandLine(TempDirCase):
    def test_parse_ncbi_assemblies_takes_the_summaries_as_n(self):
        summaries = [self.summary('assembly_summary_bacteria_refseq.txt', [summary_row('GCF_000000001.1')]),
                     self.summary('assembly_summary_bacteria_genbank.txt', [summary_row('GCA_000000002.1')])]
        out_dir = os.path.join(self.dir, 'parse_ncbi_assemblies')
        options = main_module.get_main_parser().parse_args(
            ['parse_ncbi_assemblies', '-n'] + summaries + ['-o', out_dir,
                                                          '-l', os.path.join(self.dir, 'run.log')])

        with mock.patch.object(main_py, 'NCBIMeta') as meta:
            main_py.OptionsParser().parse_options(options)
        meta.return_value.parse_assemblies.assert_called_once_with(summaries, out_dir)
        self.assertTrue(os.path.isdir(out_dir))

    def test_an_out_dir_and_a_log_are_required(self):
        summary = self.summary('assembly_summary_bacteria_refseq.txt', [summary_row('GCF_000000001.1')])
        complete = ['parse_ncbi_assemblies', '-n', summary,
                    '-o', os.path.join(self.dir, 'out'), '-l', os.path.join(self.dir, 'run.log')]
        main_module.get_main_parser().parse_args(complete)
        for flag in ('-o', '-l'):
            at = complete.index(flag)
            with mock.patch('sys.stderr'), self.assertRaises(SystemExit) as ended:
                main_module.get_main_parser().parse_args(complete[:at] + complete[at + 2:])
            self.assertEqual(ended.exception.code, 2)

    def test_the_log_is_written_where_l_says(self):
        summary = self.summary('assembly_summary_bacteria_refseq.txt', [summary_row('GCF_000000001.1')])
        log = os.path.join(self.dir, 'parse_ncbi_assemblies.log')
        options = main_module.get_main_parser().parse_args(
            ['parse_ncbi_assemblies', '-n', summary,
             '-o', os.path.join(self.dir, 'out'), '-l', log])

        self.assertEqual(options.log, log)
        self.assertEqual(main_module.log_candidates(options.log, options.output_dir)[0],
                         (self.dir, 'parse_ncbi_assemblies.log'))

    def test_the_command_is_parse_ncbi_assemblies_not_parse_assemblies(self):
        summary = self.summary('assembly_summary_bacteria_refseq.txt', [summary_row('GCF_000000001.1')])
        with mock.patch('sys.stderr'), self.assertRaises(SystemExit) as ended:
            main_module.get_main_parser().parse_args(
                ['parse_assemblies', '-n', summary,
                 '-o', os.path.join(self.dir, 'out'), '-l', os.path.join(self.dir, 'run.log')])
        self.assertEqual(ended.exception.code, 2)

    def test_a_genome_list_is_no_longer_taken(self):
        summary = self.summary('assembly_summary_bacteria_refseq.txt', [summary_row('GCF_000000001.1')])
        with mock.patch('sys.stderr'), self.assertRaises(SystemExit) as ended:
            main_module.get_main_parser().parse_args(
                ['parse_ncbi_assemblies', '-n', summary, '-m', os.path.join(self.dir, 'metadata.tsv'),
                 '-o', os.path.join(self.dir, 'out'), '-l', os.path.join(self.dir, 'run.log')])
        self.assertEqual(ended.exception.code, 2)

    def test_the_four_summary_arguments_are_no_longer_accepted(self):
        summary = self.summary('assembly_summary_bacteria_refseq.txt', [summary_row('GCF_000000001.1')])
        with mock.patch('sys.stderr'), self.assertRaises(SystemExit):
            main_module.get_main_parser().parse_args(
                ['parse_ncbi_assemblies', '--rb', summary, '--ra', summary, '--gb', summary, '--ga', summary,
                 '-o', os.path.join(self.dir, 'out'),
                 '-l', os.path.join(self.dir, 'run.log')])



# ---------------------------------------------------------------- parse_ncbi_dir

ASSEMBLY_STATS = """# Assembly Statistics Report
# Assembly name:  ASM15701v1
# Organism name:  Bacteroides fragilis 3_1_12 (CFB group bacteria)
# Taxid:          457424
# BioSample:      SAMN02463689
# Submitter:      Broad Institute
# Date:           2009-02-05
# Assembly type:  na
# Release type:   major
# Assembly level: Scaffold
# Genome representation: full
# GenBank assembly accession: GCA_000157015.1
# RefSeq assembly and GenBank assemblies identical: yes
#
# Statistic Types
# Statistic\tDescription
# total-length\tTotal sequence length
#
# Sequence-type Description
all\tall\tall\tall\ttotal-length\t5530115
all\tall\tall\tall\tcontig-count\t104
all\tall\tall\tall\tscaffold-N50\t748878
"""

GFF = ('contig1\tGenbank\tCDS\t1\t300\t.\t+\t0\tID=cds-1;protein_id=WP_1.1;transl_table=11\n'
       'contig1\tGenbank\tCDS\t400\t900\t.\t+\t0\tID=cds-2;protein_id=WP_2.1;transl_table=11\n'
       'contig1\tGenbank\ttRNA\t1000\t1075\t.\t+\t.\tID=rna-1\n'
       'contig1\tGenbank\trRNA\t2000\t3500\t.\t+\t.\tID=rna-2;product=16S ribosomal RNA\n')

GBFF = """LOCUS       contig1
FEATURES             Location/Qualifiers
     source          1..5000
                     /organism="Bacteroides fragilis"
                     /isolation_source="human gut"
                     /geo_loc_name="Australia: Brisbane"
                     /lat_lon="27.47 S 153.02 E"
     CDS             1..300
                     /transl_table=11
//
"""


class WritingGzippedTables(TempDirCase):
    def test_parse_ncbi_assemblies_writes_its_table_gzipped_removing_an_uncompressed_one(self):
        open(os.path.join(self.dir, NCBI_ASSEMBLY_TABLE), 'w').close()
        summary = self.summary('assembly_summary_bacteria_refseq.txt', [summary_row('GCF_000000001.1')])
        output = NCBIMeta().parse_assemblies([summary], self.dir)

        with open(output, 'rb') as handle:
            self.assertEqual(handle.read(2), b'\x1f\x8b')
        self.assertFalse(os.path.exists(os.path.join(self.dir, NCBI_ASSEMBLY_TABLE)))

    def test_parse_ncbi_dir_writes_its_table_gzipped_removing_an_uncompressed_one(self):
        open(os.path.join(self.dir, N.NCBI_DIR_TABLE), 'w').close()
        genome_dirs = os.path.join(self.dir, 'genome_dirs.tsv')
        open(genome_dirs, 'w').close()
        output = N.NCBIMetaDir(1).parse_ncbi_dir(genome_dirs, self.dir)

        with open(output, 'rb') as handle:
            self.assertEqual(handle.read(2), b'\x1f\x8b')
        self.assertFalse(os.path.exists(os.path.join(self.dir, N.NCBI_DIR_TABLE)))


class ParsingTheNcbiDirectories(TempDirCase):
    """parse_ncbi_dir gives every genome a row, empty where an NCBI file is missing."""

    def genome(self, gid, stats=True, gff=True, gbff=True, proteins=False):
        assembly = gid + '_ASM1v1'
        gpath = os.path.join(self.dir, 'genomes', assembly)
        os.makedirs(gpath)
        if stats:
            with open(os.path.join(gpath, assembly + '_assembly_stats.txt'), 'w') as handle:
                handle.write(ASSEMBLY_STATS)
        if gff:
            with gzip.open(os.path.join(gpath, assembly + '_genomic.gff.gz'), 'wt') as handle:
                handle.write(GFF)
        if gbff:
            with gzip.open(os.path.join(gpath, assembly + '_genomic.gbff.gz'), 'wt') as handle:
                handle.write(GBFF)
        if proteins:
            os.makedirs(os.path.join(gpath, 'prodigal'))
            with gzip.open(os.path.join(gpath, 'prodigal', gid + '_protein.faa.gz'), 'wt') as handle:
                handle.write('>p1\nM\n')
        return gid, gpath

    def parse(self, genomes, cpus=2):
        genome_dirs = os.path.join(self.dir, 'genome_dirs.tsv')
        with open(genome_dirs, 'w') as handle:
            for gid, gpath in genomes:
                handle.write('{}\t{}\tG{}\n'.format(gid, gpath, gid[4:-2]))
        with self.assertLogs('timestamp', level='INFO') as logged:
            output = N.NCBIMetaDir(cpus).parse_ncbi_dir(genome_dirs, self.dir)
        self.assertEqual(output, os.path.join(self.dir, N.NCBI_DIR_TABLE + '.gz'))
        with open_text(output) as handle:
            lines = handle.read().splitlines()
        header = lines[0].split('\t')
        rows = {line.split('\t')[0]: dict(zip(header, line.split('\t'))) for line in lines[1:]}
        for line in lines[1:]:
            self.assertEqual(len(line.split('\t')), len(header), line)
        return rows, [record for record in logged.records if 'Identified' in record.getMessage()]

    def test_a_genome_with_every_file_has_them_all_read(self):
        rows, summaries = self.parse([self.genome('GCA_000000001.1', proteins=True)])

        row = rows['GCA_000000001.1']
        self.assertEqual(row['ncbi_assembly_name'], 'ASM15701v1')
        self.assertEqual(row['ncbi_total_length'], '5530115')
        self.assertEqual((row['ncbi_cds_count'], row['ncbi_trna_count'], row['ncbi_ssu_count']), ('2', '1', '1'))
        self.assertEqual(row['ncbi_translation_table'], '11')
        self.assertEqual(row['ncbi_isolation_source'], 'human gut')
        self.assertEqual((row['ncbi_country'], row['ncbi_protein_count']), ('Australia: Brisbane', '2'))
        self.assertEqual(summaries, [])

    def test_a_genome_without_proteins_is_given_its_row(self):
        # nothing is read from prodigal's proteins; a genome without them had no row
        with_proteins = self.parse([self.genome('GCA_000000001.1', proteins=True)])[0]['GCA_000000001.1']
        shutil.rmtree(os.path.join(self.dir, 'genomes'))
        without = self.parse([self.genome('GCA_000000001.1', proteins=False)])[0]['GCA_000000001.1']

        self.assertEqual(without, with_proteins)

    def test_a_genome_without_assembly_statistics_has_a_full_row_empty_in_those_fields(self):
        rows, summaries = self.parse([self.genome('GCA_000000002.1', stats=False)])

        row = rows['GCA_000000002.1']
        self.assertEqual((row['ncbi_assembly_name'], row['ncbi_total_length']), ('', ''))
        self.assertEqual(row['ncbi_cds_count'], '2')
        self.assertEqual(row['ncbi_translation_table'], '11')
        self.assertEqual([(r.levelname, r.getMessage().split(';')[0]) for r in summaries],
                         [('WARNING', 'Identified 1 genomes with a missing _assembly_stats.txt file, '
                                      'e.g.: GCA_000000002.1')])

    def test_a_missing_gbff_is_one_warning_naming_the_genomes(self):
        rows, summaries = self.parse([self.genome('GCA_000000003.1', gbff=False),
                                      self.genome('GCA_000000004.1', gbff=False),
                                      self.genome('GCA_000000005.1')])

        # the translation table is the GFF's, the source qualifiers the GenBank file's
        self.assertEqual(rows['GCA_000000003.1']['ncbi_isolation_source'], '')
        self.assertEqual(rows['GCA_000000003.1']['ncbi_translation_table'], '11')
        self.assertEqual(rows['GCA_000000005.1']['ncbi_isolation_source'], 'human gut')
        self.assertEqual(len(summaries), 1)
        self.assertEqual(summaries[0].levelname, 'WARNING')
        self.assertIn('Identified 2 genomes with a missing _genomic.gbff.gz file, e.g.: '
                      'GCA_000000003.1, GCA_000000004.1', summaries[0].getMessage())

    def test_an_unannotated_genome_is_told_of_not_warned_of(self):
        # NCBI publishes no GFF for many GenBank assemblies
        rows, summaries = self.parse([self.genome('GCA_000000006.1', gff=False)])

        self.assertEqual(rows['GCA_000000006.1']['ncbi_cds_count'], '')
        # an unannotated genome declares no translation table
        self.assertEqual(rows['GCA_000000006.1']['ncbi_translation_table'], '')
        self.assertEqual([r.levelname for r in summaries], ['INFO'])
        self.assertIn('missing _genomic.gff.gz file', summaries[0].getMessage())

    def test_every_genome_is_written_whatever_the_order_the_workers_finish_in(self):
        genomes = [self.genome('GCA_{:09d}.1'.format(i), gff=i % 2 == 0) for i in range(30)]
        rows, _ = self.parse(genomes, cpus=4)

        self.assertEqual(sorted(rows), sorted(gid for gid, _ in genomes))


SECOND_RECORD = """LOCUS       contig2
FEATURES             Location/Qualifiers
     source          1..800
                     /isolation_source="Human gut"
                     /lat_lon="0.00 N 0.00 E"
//
"""


class ReadingAsLittleAsItCan(TempDirCase):
    """parse_ncbi_dir reads a GenBank file's first record, as far as its source feature."""

    def gbff(self, *members):
        """A gzipped GenBank file of the members given, each a str (gzipped) or bytes (as they are)."""
        path = os.path.join(self.dir, 'genome_genomic.gbff.gz')
        with open(path, 'wb') as handle:
            for member in members:
                handle.write(gzip.compress(member.encode()) if isinstance(member, str) else member)
        return path

    def test_the_source_qualifiers_are_the_first_records(self):
        # the last record's were taken; the records of r237 disagree in case alone
        values = N.NCBIMetaDir()._parse_gbff(self.gbff(GBFF + SECOND_RECORD))
        fields = dict(zip(N.NCBIMetaDir().gbff_fields, values))

        self.assertEqual((fields['isolation_source'], fields['lat_lon']), ('human gut', '27.47 S 153.02 E'))
        self.assertEqual(fields['translation_table'], '')

    def test_the_file_is_not_read_past_the_first_records_source_feature(self):
        # a full read meets the bytes that follow it, which are not gzip
        path = self.gbff(GBFF, b'not gzip at all' * 1000)
        with self.assertRaises(Exception):
            with gzip.open(path, 'rt') as handle:
                handle.read()

        fields = dict(zip(N.NCBIMetaDir().gbff_fields, N.NCBIMetaDir()._parse_gbff(path)))
        self.assertEqual(fields['isolation_source'], 'human gut')

    def test_a_record_with_no_source_feature_is_read_to_its_end_and_no_further(self):
        path = self.gbff('LOCUS       contig1\nFEATURES             Location/Qualifiers\n//\n',
                         b'not gzip at all' * 1000)
        self.assertEqual(set(N.NCBIMetaDir()._parse_gbff(path)), {''})

    def test_the_country_is_the_geo_loc_name_or_the_older_country_qualifier(self):
        # read as geo_loc_name, a column no description named, so ncbi_country held
        # only values from before the toolkit: 0.0% of r237's new genomes
        fields = dict(zip(N.NCBIMetaDir().gbff_fields, N.NCBIMetaDir()._parse_gbff(self.gbff(GBFF))))
        self.assertEqual(fields['country'], 'Australia: Brisbane')

        older = GBFF.replace('/geo_loc_name="Australia: Brisbane"', '/country="USA: Los Angeles"')
        fields = dict(zip(N.NCBIMetaDir().gbff_fields, N.NCBIMetaDir()._parse_gbff(self.gbff(older))))
        self.assertEqual(fields['country'], 'USA: Los Angeles')

    def test_the_protein_count_is_the_cds_features_with_a_protein_a_split_one_once(self):
        # NCBI's 'CDSs (with protein)': over four r237 genomes it is what the
        # GenBank file's annotation summary says, where the CDS lines are more
        path = os.path.join(self.dir, 'genome_genomic.gff.gz')
        with gzip.open(path, 'wt') as handle:
            handle.write('##gff-version 3\n'
                         'c1\tRefSeq\tCDS\t1\t300\t.\t+\t0\tID=cds-A;protein_id=WP_1.1\n'
                         'c1\tRefSeq\tCDS\t400\t600\t.\t+\t0\tID=cds-B;protein_id=WP_2.1\n'
                         'c1\tRefSeq\tCDS\t600\t900\t.\t+\t0\tID=cds-B;protein_id=WP_2.1\n'
                         'c1\tRefSeq\tCDS\t1000\t1300\t.\t+\t0\tID=cds-C;pseudo=true\n'
                         'c1\tRefSeq\tCDS\t1400\t1700\t.\t+\t0\tID=cds-D;protein_id=WP_1.1\n')
        counts = dict(zip(N.NCBIMetaDir().gff_fields, N.NCBIMetaDir()._parse_gff(path)[0]))

        # a protein two CDSs share is two; the pseudogene none; the split CDS one
        self.assertEqual((counts['protein_count'], counts['cds_count']), (3, 5))

    def test_the_gff_gives_its_first_translation_table_and_builds_no_coding_mask_unasked(self):
        path = os.path.join(self.dir, 'genome_genomic.gff.gz')
        with gzip.open(path, 'wt') as handle:
            handle.write(GFF.replace('transl_table=11\n', 'transl_table=4\n', 1))
        parser = N.GenericFeatureParser(path)

        self.assertEqual(parser.translation_table, 4)
        self.assertEqual(parser.coding_mask, {})
        self.assertEqual(parser.total_coding_bases(), 300 + 501)
        self.assertEqual(N.NCBIMetaDir()._parse_gff(path)[1], 4)


class TheNcbiDirCommandLine(TempDirCase):
    def test_parse_ncbi_dir_writes_to_an_out_dir_and_requires_a_log(self):
        genome_dirs = os.path.join(self.dir, 'genome_dirs.tsv')
        open(genome_dirs, 'w').close()
        out_dir = os.path.join(self.dir, 'parse_ncbi_dir')
        complete = ['parse_ncbi_dir', '-g', genome_dirs, '-o', out_dir,
                    '-l', os.path.join(self.dir, 'run.log'), '--cpus', '3']
        options = main_module.get_main_parser().parse_args(complete)

        with mock.patch.object(main_py, 'NCBIMetaDir') as meta:
            main_py.OptionsParser().parse_options(options)
        meta.assert_called_once_with(3)
        meta.return_value.parse_ncbi_dir.assert_called_once_with(genome_dirs, out_dir)
        self.assertTrue(os.path.isdir(out_dir))

        for flag in ('-o', '-l'):
            at = complete.index(flag)
            with mock.patch('sys.stderr'), self.assertRaises(SystemExit) as ended:
                main_module.get_main_parser().parse_args(complete[:at] + complete[at + 2:])
            self.assertEqual(ended.exception.code, 2)

    def test_the_table_is_named_as_update_metadata_db_loads_it_and_not_as_parse_ncbi_assemblies_names_its(self):
        with mock.patch.object(metadata_database_manager.GenomeDatabaseConnectionFTPUpdate,
                               'GenomeDatabaseConnectionFTPUpdate'):
            manager = metadata_database_manager.MetadataDatabaseManager({})

        self.assertEqual(manager.description_table[N.NCBI_DIR_TABLE], ['metadata_ncbi_assembly.desc.tsv'])
        self.assertNotEqual(N.NCBI_DIR_TABLE, NCBI_ASSEMBLY_TABLE)

    def test_every_column_is_loaded_but_the_three_the_database_has_no_field_for(self):
        # ncbi_isolation_source and ncbi_lat_lon were read from the GenBank file
        # and never loaded, the description not naming them, and the country was
        # read as ncbi_geo_loc_name until 0.1.78; the database has no column for
        # these three (checked against gtdb_r237_dev)
        path = os.path.join(os.path.dirname(metadata_database_manager.__file__),
                            'data_files', 'table_description', 'metadata_ncbi_assembly.desc.tsv')
        with open(path) as handle:
            described = {line.split('\t')[0] for line in handle if line.strip()}
        written = set(N.NCBIMetaDir().header().rstrip('\n').split('\t')[1:])

        self.assertEqual(written - described, {'ncbi_contig_l50', 'ncbi_component_count',
                                               'ncbi_metagenome_source'})
        self.assertLessEqual({'ncbi_isolation_source', 'ncbi_lat_lon', 'ncbi_country', 'ncbi_protein_count'},
                             described)

    def test_the_organism_name_is_neither_written_nor_described_update_ncbi_tax_db_alone_writes_it(self):
        # the assembly report's name was loaded from this table, after update_ncbi_tax_db
        # had written the taxonomy's, and replaced it for 13,538 genomes of r237
        self.assertNotIn('ncbi_organism_name', N.NCBIMetaDir().header().rstrip('\n').split('\t'))

        directory = os.path.join(os.path.dirname(metadata_database_manager.__file__),
                                 'data_files', 'table_description')
        owned = {field for _, field in metadata_database_manager.NCBI_TAX_FIELDS}
        for name in sorted(os.listdir(directory)):
            with open(os.path.join(directory, name)) as handle:
                described = {line.split('\t')[0] for line in handle if line.strip()}
            self.assertEqual(described & owned, set(), name)


if __name__ == '__main__':
    unittest.main()
