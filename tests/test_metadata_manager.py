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

import gzip
import logging
import os
import shutil
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


if __name__ == '__main__':
    unittest.main()
