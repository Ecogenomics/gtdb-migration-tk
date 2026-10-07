#!/usr/bin/env python3
"""Offline unit tests for metadata_manager.py -- no release and no mirror.

Run with the interpreter that has tqdm and numpy:

    /opt/centos7/sw/miniconda3/envs/gtdb_migration_tk-r237/bin/python -m unittest -v tests.test_metadata_manager

The metadata of a release is generated while the gene calling of its last genomes
is still finishing, so what these cover is what happens to a genome whose files
are not all there. That used to end the command -- check_file_exists() calls
sys.exit(), inside a pool worker, on the first genome missing anything -- and a
run of a million genomes would be thrown away for one straggler. The genomes are
real enough to be calculated over: a few hundred bases of FASTA and a GFF with one
CDS in it, so the calculators run rather than being stubbed out.
"""

import contextlib
import gzip
import io
import logging
import os
import shutil
import sys
import tempfile
import unittest
from unittest import mock

from gtdb_migration_tk import metadata_manager as M


class TempDirCase(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='metadata_manager_test.')

        # the command draws a progress bar and warns once per missing file.
        # Neither belongs in the output of a test run, and assertLogs() raises
        # the logger back up for the tests that are about what it says.
        patch = mock.patch.object(M, 'tqdm', lambda iterable, **kwargs: iterable)
        patch.start()
        self.addCleanup(patch.stop)

        logger = logging.getLogger('timestamp')
        self.addCleanup(logger.setLevel, logger.level)
        logger.setLevel(logging.CRITICAL)

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def genome_dir(self, gid, assembly, fasta=True, gff=True):
        """A genome directory as a release holds it, with or without its files."""
        gpath = os.path.join(self.dir, assembly)
        os.makedirs(os.path.join(gpath, 'prodigal'))

        if fasta:
            path = os.path.join(gpath, assembly + '_genomic.fna.gz')
            with gzip.open(path, 'wt') as handle:
                handle.write('>contig_1 an assembly of one contig\n')
                handle.write('ATGCGCGCATGCATGCATGC' * 20 + '\n')

        if gff:
            path = os.path.join(gpath, 'prodigal', gid + '_protein.gff.gz')
            with gzip.open(path, 'wt') as handle:
                handle.write('##gff-version 3\n')
                handle.write('contig_1\tProdigal\tCDS\t1\t99\t.\t+\t0\tID=1_1\n')

        return gpath

    def genome_dirs_file(self, genomes):
        """The genome_dirs file of a release: accession, directory, canonical."""
        path = os.path.join(self.dir, 'genome_dirs.tsv')
        with open(path, 'w') as handle:
            for gid, gpath in genomes:
                handle.write('{}\t{}\t{}\n'.format(gid, gpath, 'G' + gid[4:-2]))
        return path

    def run_metadata(self, genomes, cpus=1):
        """Generate the metadata of a release, as the command does."""
        out_dir = os.path.join(self.dir, 'out')
        manager = M.MetadataManager(cpus)
        manager.generate_metadata(self.genome_dirs_file(genomes), out_dir)
        return out_dir

    def report_rows(self, out_dir):
        """The missing files report, as (genome_id, missing) pairs."""
        path = os.path.join(out_dir, M.MISSING_FILES_NAME)
        with open(path) as handle:
            header = handle.readline().rstrip('\n').split('\t')
            self.assertEqual(tuple(header), M.MISSING_FILES_HEADER)
            return [tuple(line.rstrip('\n').split('\t')[:2]) for line in handle]

    def wrote(self, gpath, name):
        return os.path.exists(os.path.join(gpath, name))


# --------------------------------------------- a genome that has both its files

class AWholeGenome(TempDirCase):
    """Nothing about the reporting stops an ordinary genome being processed."""

    def test_both_metadata_files_and_their_descriptions_are_written(self):
        gpath = self.genome_dir('GCA_000001.1', 'GCA_000001.1_ASM1v1')

        self.run_metadata([('GCA_000001.1', gpath)])

        for name in ('metadata.genome_nt.tsv', 'metadata.genome_nt.desc.tsv',
                     'metadata.genome_gene.tsv', 'metadata.genome_gene.desc.tsv'):
            self.assertTrue(self.wrote(gpath, name), name + ' was not written')

    def test_it_is_not_named_in_the_report(self):
        gpath = self.genome_dir('GCA_000001.1', 'GCA_000001.1_ASM1v1')

        out_dir = self.run_metadata([('GCA_000001.1', gpath)])

        self.assertEqual(self.report_rows(out_dir), [])


# ------------------------------------------------- a genome missing its GFF only

class AGenomeWhoseGenesAreNotCalledYet(TempDirCase):
    """The nucleotide half needs only the FASTA, and is worth having on its own.

    create_metadata_tables() reads metadata.genome_nt.tsv and
    metadata.genome_gene.tsv independently, so half a genome is a row in one
    table rather than a broken genome, and it is not calculated again when
    prodigal catches up.
    """

    def setUp(self):
        super().setUp()
        self.gpath = self.genome_dir('GCA_000002.1', 'GCA_000002.1_ASM2v1', gff=False)
        self.out_dir = self.run_metadata([('GCA_000002.1', self.gpath)])

    def test_the_nucleotide_metadata_is_written(self):
        self.assertTrue(self.wrote(self.gpath, 'metadata.genome_nt.tsv'))

    def test_the_gene_metadata_is_not(self):
        self.assertFalse(self.wrote(self.gpath, 'metadata.genome_gene.tsv'))

    def test_the_report_says_which_file_was_missing(self):
        self.assertEqual(self.report_rows(self.out_dir),
                         [('GCA_000002.1', M.MISSING_PROTEIN_GFF)])

    def test_the_report_gives_the_path_that_was_looked_for(self):
        with open(os.path.join(self.out_dir, M.MISSING_FILES_NAME)) as handle:
            handle.readline()
            path = handle.readline().rstrip('\n').split('\t')[2]

        self.assertEqual(path, os.path.join(self.gpath, 'prodigal',
                                            'GCA_000002.1_protein.gff.gz'))


# --------------------------------------- a genome whose GFF has gone since a run

class AGenomeWhoseGenesHaveGone(TempDirCase):
    """Gene metadata of an earlier run does not outlive the GFF it came from.

    create_metadata_tables() reads metadata.genome_gene.tsv wherever it is. r237
    called 15 genomes as empty proteomes, this command gave each a protein count
    of 0, and the patch that called them again left 3 with no GFF and the 0, to
    be loaded as though Prodigal had found no genes.
    """

    def setUp(self):
        super().setUp()
        self.gpath = self.genome_dir('GCA_000005.1', 'GCA_000005.1_ASM5v1')
        self.run_metadata([('GCA_000005.1', self.gpath)])
        os.remove(os.path.join(self.gpath, 'prodigal', 'GCA_000005.1_protein.gff.gz'))

    def test_the_gene_metadata_of_the_earlier_run_is_removed(self):
        self.run_metadata([('GCA_000005.1', self.gpath)])

        for name in M.GENE_METADATA_FILES:
            self.assertFalse(self.wrote(self.gpath, name), name + ' was left')

    def test_the_nucleotide_metadata_is_still_written(self):
        self.run_metadata([('GCA_000005.1', self.gpath)])

        self.assertTrue(self.wrote(self.gpath, 'metadata.genome_nt.tsv'))

    def test_create_tables_gives_it_no_gene_row(self):
        self.run_metadata([('GCA_000005.1', self.gpath)])
        table = M.MetadataTable('138.2')
        # the file removed is the one create_tables reads
        self.assertEqual(table.metadata_gene_file, M.GENE_METADATA_FILES[0])
        self.assertIsNone(table._read_field_table(
            'GCA_000005.1', os.path.join(self.gpath, table.metadata_gene_file)))

    def test_the_removal_is_warned_of_with_the_missing_gff_and_counted(self):
        logger = logging.getLogger('timestamp')
        logger.setLevel(logging.WARNING)
        with self.assertLogs('timestamp', level='WARNING') as logged:
            self.run_metadata([('GCA_000005.1', self.gpath)])

        self.assertTrue(any('GCA_000005.1 has no called genes (GFF)' in line
                            and 'removed the gene metadata an earlier run wrote' in line
                            for line in logged.output), logged.output)
        self.assertTrue(any('for 1 genome(s) that no longer have called genes' in line
                            for line in logged.output), logged.output)

    def test_a_genome_that_never_had_gene_metadata_is_not_said_to_have_lost_it(self):
        gpath = self.genome_dir('GCA_000006.1', 'GCA_000006.1_ASM6v1', gff=False)
        logger = logging.getLogger('timestamp')
        logger.setLevel(logging.WARNING)
        with self.assertLogs('timestamp', level='WARNING') as logged:
            self.run_metadata([('GCA_000006.1', gpath)])

        self.assertFalse(any('removed' in line.lower() for line in logged.output), logged.output)

    def test_a_genome_with_no_sequences_keeps_its_gene_metadata(self):
        # with no FASTA the genome directory is left as it was found: the mirror
        # has lost a file, which is not the genes having gone
        os.remove(os.path.join(self.gpath, 'GCA_000005.1_ASM5v1_genomic.fna.gz'))
        self.run_metadata([('GCA_000005.1', self.gpath)])

        self.assertTrue(self.wrote(self.gpath, M.GENE_METADATA_FILES[0]))


# ----------------------------------------------- a genome missing its sequences

class AGenomeWithNoSequences(TempDirCase):
    """Without the FASTA there is nothing to calculate, the gene metadata included.

    The gene metadata is a coding density, which needs the genome size as well as
    the GFF, so a genome with called genes and no sequences yields neither file.
    """

    def setUp(self):
        super().setUp()
        self.gpath = self.genome_dir('GCA_000003.1', 'GCA_000003.1_ASM3v1', fasta=False)
        self.out_dir = self.run_metadata([('GCA_000003.1', self.gpath)])

    def test_neither_metadata_file_is_written(self):
        self.assertFalse(self.wrote(self.gpath, 'metadata.genome_nt.tsv'))
        self.assertFalse(self.wrote(self.gpath, 'metadata.genome_gene.tsv'))

    def test_the_report_says_the_fasta_was_missing(self):
        self.assertEqual(self.report_rows(self.out_dir),
                         [('GCA_000003.1', M.MISSING_GENOMIC_FASTA)])

    def test_a_genome_it_cannot_process_keeps_the_log_of_the_run_that_could(self):
        """The old genometk.log is removed by the run that replaces it, not before."""
        gpath = self.genome_dir('GCA_000004.1', 'GCA_000004.1_ASM4v1', fasta=False)
        log = os.path.join(gpath, 'genometk.log')
        open(log, 'w').close()

        self.run_metadata([('GCA_000004.1', gpath)])

        self.assertTrue(os.path.exists(log))


# ------------------------------------------------------------ the release as one

class WhatTheRunReports(TempDirCase):
    """One release of stragglers and whole genomes, as a run in progress has."""

    def setUp(self):
        super().setUp()
        self.whole = self.genome_dir('GCA_000001.1', 'GCA_000001.1_ASM1v1')
        self.no_gff = self.genome_dir('GCA_000002.1', 'GCA_000002.1_ASM2v1', gff=False)
        self.no_fasta = self.genome_dir('GCA_000003.1', 'GCA_000003.1_ASM3v1', fasta=False)
        self.nothing = self.genome_dir('GCA_000004.1', 'GCA_000004.1_ASM4v1',
                                       fasta=False, gff=False)
        self.out_dir = self.run_metadata(
            [('GCA_000001.1', self.whole), ('GCA_000002.1', self.no_gff),
             ('GCA_000003.1', self.no_fasta), ('GCA_000004.1', self.nothing)],
            cpus=2)

    def test_one_genome_missing_a_file_does_not_cost_the_release(self):
        """What the exit in check_file_exists() used to do to a run of a million."""
        self.assertTrue(self.wrote(self.whole, 'metadata.genome_gene.tsv'))

    def test_every_missing_file_is_named_once(self):
        self.assertEqual(sorted(self.report_rows(self.out_dir)),
                         [('GCA_000002.1', M.MISSING_PROTEIN_GFF),
                          ('GCA_000003.1', M.MISSING_GENOMIC_FASTA),
                          ('GCA_000004.1', M.MISSING_GENOMIC_FASTA),
                          ('GCA_000004.1', M.MISSING_PROTEIN_GFF)])

    def test_a_genome_missing_both_is_one_row_per_file(self):
        rows = [row for row in self.report_rows(self.out_dir) if row[0] == 'GCA_000004.1']

        self.assertEqual(len(rows), 2)


class TheReportIsWrittenEitherWay(TempDirCase):
    """A release with nothing missing says so, rather than leaving no file."""

    def test_a_release_with_nothing_missing_still_gets_a_report(self):
        gpath = self.genome_dir('GCA_000001.1', 'GCA_000001.1_ASM1v1')

        out_dir = self.run_metadata([('GCA_000001.1', gpath)])

        self.assertTrue(os.path.exists(os.path.join(out_dir, M.MISSING_FILES_NAME)))
        self.assertEqual(self.report_rows(out_dir), [])

    def test_the_out_dir_is_made_if_it_is_not_there(self):
        gpath = self.genome_dir('GCA_000001.1', 'GCA_000001.1_ASM1v1')
        out_dir = os.path.join(self.dir, 'made', 'by', 'the', 'run')

        M.MetadataManager(1).generate_metadata(
            self.genome_dirs_file([('GCA_000001.1', gpath)]), out_dir)

        self.assertTrue(os.path.exists(os.path.join(out_dir, M.MISSING_FILES_NAME)))

    def test_no_metadata_is_written_to_the_out_dir(self):
        """--out_dir holds the account of the run; the metadata goes to the genomes."""
        gpath = self.genome_dir('GCA_000001.1', 'GCA_000001.1_ASM1v1')

        out_dir = self.run_metadata([('GCA_000001.1', gpath)])

        self.assertEqual(sorted(os.listdir(out_dir)), [M.MISSING_FILES_NAME])


# ---------------------------------------------------------------- and the log

class WhatTheLogSays(TempDirCase):
    """The report is the list; the log is what is watched while the run goes on."""

    def test_a_missing_file_is_warned_about_as_it_is_met(self):
        gpath = self.genome_dir('GCA_000002.1', 'GCA_000002.1_ASM2v1', gff=False)

        with self.assertLogs('timestamp', level='WARNING') as captured:
            self.run_metadata([('GCA_000002.1', gpath)])

        warned = [record.getMessage() for record in captured.records
                  if record.levelno == logging.WARNING]
        self.assertTrue(any('GCA_000002.1' in message and 'GCA_000002.1_protein.gff.gz' in message
                            for message in warned), warned)

    def test_the_run_ends_by_saying_how_many_genomes_were_short_a_file(self):
        whole = self.genome_dir('GCA_000001.1', 'GCA_000001.1_ASM1v1')
        no_gff = self.genome_dir('GCA_000002.1', 'GCA_000002.1_ASM2v1', gff=False)

        with self.assertLogs('timestamp', level='WARNING') as captured:
            self.run_metadata([('GCA_000001.1', whole), ('GCA_000002.1', no_gff)])

        self.assertIn('1 of 2 genomes were missing a file',
                      captured.records[-1].getMessage())


# ------------------------------------------------------------ create_tables

TRNA_STATS = ('tRNAs decoding Standard 20 AA:              40\n'
              'Selenocysteine tRNAs (TCA):                 0\n')


class WhatCreateTablesLogs(TempDirCase):
    """create_tables passes over a genome without a file, and the log says how many.

    A genome missing a table's file is given no row in it and nothing else says
    so, which across a release is how genomic_metadata, rna_silva or trnascan
    not getting to every genome would go unseen until the database was loaded.
    """

    def gathered_genome(self, gid, nt=True, gene=True, trna=True):
        """A genome directory holding what earlier commands wrote into it."""
        gpath = os.path.join(self.dir, gid)
        os.makedirs(os.path.join(gpath, 'trna'))
        if nt:
            with open(os.path.join(gpath, 'metadata.genome_nt.tsv'), 'w') as handle:
                handle.write('gc_percentage\t50.0\ngenome_size\t400\n')
        if gene:
            with open(os.path.join(gpath, 'metadata.genome_gene.tsv'), 'w') as handle:
                handle.write('protein_count\t1\n')
        if trna:
            with open(os.path.join(gpath, 'trna', gid + '_trna_stats.tsv'), 'w') as handle:
                handle.write(TRNA_STATS)
        return gpath

    def create_tables(self, genomes):
        out_dir = os.path.join(self.dir, 'tables')
        logger = logging.getLogger('timestamp')
        logger.setLevel(logging.INFO)
        with self.assertLogs('timestamp', level='INFO') as captured:
            M.MetadataTable('138.2').create_metadata_tables(
                self.genome_dirs_file(genomes), out_dir)
        return out_dir, [record.getMessage() for record in captured.records]

    def line_for(self, messages, table):
        lines = [message for message in messages if message.strip().startswith(table + ':')]
        self.assertEqual(len(lines), 1, messages)
        return lines[0]

    def test_each_table_is_given_the_genomes_it_has_a_row_for(self):
        whole = self.gathered_genome('GCA_000001.1')
        no_gene = self.gathered_genome('GCA_000002.1', gene=False)
        out_dir, messages = self.create_tables([('GCA_000001.1', whole), ('GCA_000002.1', no_gene)])

        self.assertIn('Rows written for 2 genomes:', messages)
        self.assertEqual(self.line_for(messages, M.NT_TABLE).strip(), 'metadata_nt.tsv: 2 with a row.')
        self.assertEqual(self.line_for(messages, M.GENE_TABLE).strip(),
                         'metadata_gene.tsv: 1 with a row; 1 had no metadata.genome_gene.tsv.')

        # and the count is of the rows the table holds, header aside
        with open(os.path.join(out_dir, M.GENE_TABLE)) as handle:
            self.assertEqual(len(handle.readlines()) - 1, 1)

    def test_a_genome_missing_a_file_is_counted_against_the_table_read_from_it(self):
        no_trna = self.gathered_genome('GCA_000003.1', trna=False)
        _, messages = self.create_tables([('GCA_000003.1', no_trna)])

        self.assertIn('1 had no ' + os.path.join('trna', '<gid>_trna_stats.tsv'),
                      self.line_for(messages, M.TRNA_TABLE))
        self.assertIn('1 had no ' + os.path.join('rna_silva_138.2', 'ssu.taxonomy.tsv'),
                      self.line_for(messages, M.taxonomy_table('ssu_silva')))

    def test_every_table_with_rows_only_for_some_genomes_is_named(self):
        _, messages = self.create_tables([('GCA_000001.1', self.gathered_genome('GCA_000001.1'))])

        for table in (M.NT_TABLE, M.GENE_TABLE, M.taxonomy_table('ssu_gg'),
                      M.taxonomy_table('ssu_silva'), M.taxonomy_table('lsu_silva_23s'),
                      M.LSU_5S_TABLE, M.TRNA_TABLE):
            self.line_for(messages, table)

    def test_a_file_with_nothing_to_report_is_not_called_missing(self):
        gpath = self.gathered_genome('GCA_000004.1')
        silva = os.path.join(gpath, 'rna_silva_138.2')
        os.makedirs(silva)
        with open(os.path.join(silva, 'ssu.taxonomy.tsv'), 'w') as handle:
            handle.write('query_id\ttaxonomy\tlength\n')
        _, messages = self.create_tables([('GCA_000004.1', gpath)])

        self.assertEqual(self.line_for(messages, M.taxonomy_table('ssu_silva')).strip(),
                         'metadata_ssu_silva.tsv: 0 with a row; 1 had nothing to report.')

    def test_a_genome_with_no_nucleotide_metadata_is_warned_of(self):
        # every genome has a FASTA, so one without this is a genome
        # genomic_metadata did not get to
        whole = self.gathered_genome('GCA_000001.1')
        no_nt = self.gathered_genome('GCA_000002.1', nt=False)
        out_dir = os.path.join(self.dir, 'tables')
        logger = logging.getLogger('timestamp')
        logger.setLevel(logging.INFO)
        with self.assertLogs('timestamp', level='INFO') as captured:
            M.MetadataTable('138.2').create_metadata_tables(
                self.genome_dirs_file([('GCA_000001.1', whole), ('GCA_000002.1', no_nt)]), out_dir)

        warned = [record.getMessage().strip() for record in captured.records
                  if record.levelno == logging.WARNING]
        self.assertEqual(warned, ['metadata_nt.tsv: 1 with a row; 1 had no metadata.genome_nt.tsv.'])

    def test_a_table_rightly_missing_for_some_genomes_is_not_warned_of(self):
        # no proteins called, no tRNA or rRNA gene found: not a step left undone
        gpath = self.gathered_genome('GCA_000003.1', gene=False, trna=False)
        logger = logging.getLogger('timestamp')
        logger.setLevel(logging.INFO)
        with self.assertLogs('timestamp', level='INFO') as captured:
            M.MetadataTable('138.2').create_metadata_tables(
                self.genome_dirs_file([('GCA_000003.1', gpath)]), os.path.join(self.dir, 'tables'))

        self.assertFalse([record for record in captured.records if record.levelno >= logging.WARNING])


class AnEmptyGenomeDirsFile(TempDirCase):
    """A genome_dirs file naming no genomes ends the command before anything is written.

    Ten tables of no rows would replace the release's in --out_dir, for
    update_metadata_db to load.
    """

    def test_it_is_refused_and_no_tables_are_written(self):
        out_dir = os.path.join(self.dir, 'tables')
        with self.assertRaises(M.EmptyGenomeDirs):
            M.MetadataTable('138.2').create_metadata_tables(self.genome_dirs_file([]), out_dir)
        self.assertFalse(os.path.exists(out_dir))

    def test_an_earlier_releases_tables_are_left_as_they_were(self):
        out_dir = os.path.join(self.dir, 'tables')
        os.makedirs(out_dir)
        table = os.path.join(out_dir, M.NT_TABLE)
        with open(table, 'w') as handle:
            handle.write('genome_id\tgc_percentage\nGCA_000001.1\t50.0\n')

        with self.assertRaises(M.EmptyGenomeDirs):
            M.MetadataTable('138.2').create_metadata_tables(self.genome_dirs_file([]), out_dir)
        with open(table) as handle:
            self.assertEqual(handle.read(), 'genome_id\tgc_percentage\nGCA_000001.1\t50.0\n')

    def test_the_command_says_so_and_exits_1(self):
        from gtdb_migration_tk import __main__ as main_module
        for name in ('timestamp', 'no_timestamp'):
            self.addCleanup(self.drop_handlers, logging.getLogger(name))
        genome_dirs = self.genome_dirs_file([])
        log = os.path.join(self.dir, 'create_tables.log')
        argv = ['gtdb_migration_tk', 'create_tables', '-g', genome_dirs,
                '-o', os.path.join(self.dir, 'tables'), '-v', '138.2', '-l', log, '--silent']

        with mock.patch.object(sys, 'argv', argv), \
                contextlib.redirect_stdout(io.StringIO()), \
                contextlib.redirect_stderr(io.StringIO()), \
                self.assertRaises(SystemExit) as ended:
            main_module.main()

        self.assertEqual(ended.exception.code, 1)
        with open(log) as handle:
            self.assertIn('ERROR: {} names no genomes'.format(genome_dirs), handle.read())

    @staticmethod
    def drop_handlers(logger):
        for handler in list(logger.handlers):
            logger.removeHandler(handler)
            handler.close()


class CreateTablesCommandLine(TempDirCase):
    def test_the_log_is_written_where_l_says(self):
        from gtdb_migration_tk import __main__ as main_module
        log = os.path.join(self.dir, 'create_tables.log')
        out_dir = os.path.join(self.dir, 'tables')
        options = main_module.get_main_parser().parse_args(
            ['create_tables', '-g', self.genome_dirs_file([]), '-o', out_dir,
             '-v', '138.2', '-l', log])

        self.assertEqual(main_module.log_candidates(options.log, options.output_dir)[0],
                         (self.dir, 'create_tables.log'))

    def test_a_run_naming_no_log_is_refused(self):
        from gtdb_migration_tk import __main__ as main_module
        out_dir = os.path.join(self.dir, 'tables')
        with mock.patch('sys.stderr'), self.assertRaises(SystemExit) as ended:
            main_module.get_main_parser().parse_args(
                ['create_tables', '-g', self.genome_dirs_file([]), '-o', out_dir, '-v', '138.2'])

        self.assertEqual(ended.exception.code, 2)


    def test_cpus_defaults_to_what_the_file_server_was_measured_to_take_and_is_passed_on(self):
        from gtdb_migration_tk import __main__ as main_module
        from gtdb_migration_tk import main as main_py
        out_dir = os.path.join(self.dir, 'tables')
        argv = ['create_tables', '-g', self.genome_dirs_file([]), '-o', out_dir, '-v', '138.2',
                '-l', os.path.join(self.dir, 'create_tables.log')]
        options = main_module.get_main_parser().parse_args(argv)
        self.assertEqual(options.cpus, M.CREATE_TABLES_THREADS)

        options = main_module.get_main_parser().parse_args(argv + ['--cpus', '3'])
        with mock.patch.object(main_py, 'MetadataTable') as table:
            main_py.OptionsParser().parse_options(options)
        table.return_value.create_metadata_tables.assert_called_once_with(
            options.gtdb_genome_path_file, out_dir, 3)


class ReadingOnThreads(TempDirCase):
    """create_tables reads genomes on threads and writes them in the order given.

    Each table is headed by the first genome that has the file it is read
    from, so a table written as genomes finished would be headed, and ordered,
    by whichever thread came back first.
    """

    def genome(self, gid, gene=True, ssu=True):
        gpath = WhatCreateTablesLogs.gathered_genome(self, gid, gene=gene)
        if ssu:
            silva = os.path.join(gpath, 'rna_silva_138.2')
            os.makedirs(silva)
            with open(os.path.join(silva, 'ssu.taxonomy.tsv'), 'w') as handle:
                handle.write('query_id\ttaxonomy\tlength\n{0}_1\td__Bacteria\t1500\n'.format(gid))
            with open(os.path.join(silva, 'ssu.fna'), 'w') as handle:
                handle.write('>{0}_1\nACGTACGT\n'.format(gid))
            with open(os.path.join(silva, 'ssu.hmm_summary.tsv'), 'w') as handle:
                handle.write('Sequence Id\tSequence length\n{0}_1\t5000\n'.format(gid))
        return gpath

    def tables(self, genome_dirs, cpus):
        out_dir = os.path.join(self.dir, 'tables_{}'.format(cpus))
        with self.assertLogs('timestamp', level='INFO'):
            M.MetadataTable('138.2').create_metadata_tables(genome_dirs, out_dir, cpus)
        written = {}
        for name in sorted(os.listdir(out_dir)):
            with open(os.path.join(out_dir, name)) as handle:
                written[name] = handle.read()
        return written

    def test_the_tables_are_the_same_whatever_cpus_is(self):
        # the first genome has no gene table and no rRNA genes, so those
        # tables are headed by a later one
        genomes = [('GCA_{:06d}.1'.format(i), self.genome('GCA_{:06d}.1'.format(i),
                                                          gene=i > 0, ssu=i % 3 == 1))
                   for i in range(40)]
        genome_dirs = self.genome_dirs_file(genomes)

        serial = self.tables(genome_dirs, 1)
        self.assertEqual(self.tables(genome_dirs, 8), serial)
        self.assertEqual(serial[M.GENE_TABLE].splitlines()[:2],
                         ['genome_id\tprotein_count', 'GCA_000001.1\t1'])
        self.assertEqual([line.split('\t')[0] for line in serial[M.NT_TABLE].splitlines()[1:]],
                         [gid for gid, _ in genomes])
        self.assertEqual(serial[M.taxonomy_table('ssu_silva')].splitlines()[:2],
                         ['genome_id\tssu_query_id\tssu_silva_taxonomy\tssu_length\tssu_sequence\tssu_contig_len',
                          'GCA_000001.1\tGCA_000001.1_1\td__Bacteria\t1500\tACGTACGT\t5000'])

    def test_results_come_back_in_the_order_given_however_they_finish(self):
        import time
        # the first items take longest, so later ones finish first
        results = list(M.ordered_map(lambda i: time.sleep((10 - i) / 1000.0) or i, range(10), 4))

        self.assertEqual(results, list(range(10)))

    def test_no_more_than_threads_times_depth_genomes_are_in_hand_at_once(self):
        drawn = []

        def items():
            for i in range(100):
                drawn.append(i)
                yield i

        in_hand = []
        for yielded, _ in enumerate(M.ordered_map(lambda i: i, items(), 3, depth=2), start=1):
            in_hand.append(len(drawn) - yielded)

        self.assertLessEqual(max(in_hand), 3 * 2)

    def test_a_genome_that_cannot_be_read_stops_the_run(self):
        def read(i):
            if i == 5:
                raise ValueError('genome 5')
            return i

        with self.assertRaisesRegex(ValueError, 'genome 5'):
            list(M.ordered_map(read, range(50), 4))


if __name__ == '__main__':
    unittest.main()
