#!/usr/bin/env python3
"""Offline unit tests for directory_manager.py -- no mirror, no network.

Run with the interpreter that has tqdm:

    /opt/centos7/sw/miniconda3/envs/python3/bin/python -m unittest -v tests.test_directory_manager

What is tested is the genome_dirs file: the lingua franca the rest of the toolkit reads,
one line per genome of the release saying where it is held. The contract that would
otherwise break silently is that the file describes THIS release -- a tree may be a GTDB
release accumulated over several cycles, and a genome_dirs file quietly carrying a
previous release's genomes would feed them to every command downstream.
"""

import gzip
import logging
import os
import shutil
import tempfile
import unittest

from gtdb_migration_tk.directory_manager import DirectoryManager


SERVED = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/000/001/x'
SELECTION_HEADER = ('#assembly_accession\tftp_path\tversion_status\texcluded_from_refseq'
                    '\tgbrs_paired_asm\tnotes')


class TempDirCase(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='directory_manager_test.')
        self.warnings = []
        self.manager = DirectoryManager()
        self.manager.logger = type('Log', (), {
            'info': lambda _s, *a: None,
            'warning': lambda _s, msg, *a: self.warnings.append(msg.format(*a) if a else msg),
        })()

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def path(self, *parts):
        return os.path.join(self.dir, *parts)

    def selection(self, *accessions, gzipped=False):
        """The table select_genomes writes, naming the given genomes."""
        text = SELECTION_HEADER + '\n' + ''.join(
            '%s\t%s\tlatest\tna\tna\tna\n' % (acc, SERVED) for acc in accessions)
        path = self.path('sel.tsv.gz' if gzipped else 'sel.tsv')
        opener = gzip.open if gzipped else open
        with opener(path, 'wt') as handle:
            handle.write(text)
        return path

    def genome(self, archive, triplets, leaf):
        """A genome directory in the tree layout: <tree>/<archive>/NNN/NNN/NNN/<leaf>."""
        path = self.path('tree', archive, *(triplets + (leaf,)))
        os.makedirs(path)
        open(os.path.join(path, 'md5checksums.txt'), 'w').close()
        return path

    def run_list(self, selection):
        out = self.path('genome_dirs.tsv')
        self.manager.generate_genome_dir_file(self.path('tree'), out, selection, cpus=2)
        with open(out) as handle:
            return [line.rstrip('\n').split('\t') for line in handle]


class ReadSelectedGenomes(TempDirCase):
    def test_reads_the_accessions_of_the_release(self):
        path = self.selection('GCF_000000001.1', 'GCA_000000002.1')
        self.assertEqual(self.manager.read_selected_genomes(path),
                         {'GCF_000000001.1', 'GCA_000000002.1'})

    def test_reads_a_gzipped_table(self):
        # select_genomes writes gtdb_selected_genomes.tsv.gz
        path = self.selection('GCF_000000001.1', gzipped=True)
        self.assertEqual(self.manager.read_selected_genomes(path), {'GCF_000000001.1'})

    def test_an_empty_release_is_an_empty_set(self):
        self.assertEqual(self.manager.read_selected_genomes(self.selection()), set())


class GenerateGenomeDirFile(TempDirCase):
    def test_writes_accession_path_and_canonical_accession(self):
        held = self.genome('GCF', ('000', '000', '001'), 'GCF_000000001.1_ASM1v1')
        rows = self.run_list(self.selection('GCF_000000001.1'))
        self.assertEqual(rows, [['GCF_000000001.1', held, 'G000000001']])

    def test_finds_genomes_in_both_archives(self):
        self.genome('GCA', ('000', '000', '002'), 'GCA_000000002.1_ASM2v1')
        self.genome('GCF', ('000', '000', '001'), 'GCF_000000001.1_ASM1v1')
        rows = self.run_list(self.selection('GCF_000000001.1', 'GCA_000000002.1'))
        self.assertEqual(sorted(row[0] for row in rows),
                         ['GCA_000000002.1', 'GCF_000000001.1'])

    def test_a_tree_with_only_one_archive_is_normal(self):
        self.genome('GCF', ('000', '000', '001'), 'GCF_000000001.1_ASM1v1')
        self.assertEqual(len(self.run_list(self.selection('GCF_000000001.1'))), 1)

    def test_a_genome_the_release_does_not_list_is_passed_over(self):
        # the tree may be a GTDB release accumulated over several cycles; the
        # genome_dirs file describes THIS release
        self.genome('GCF', ('000', '000', '001'), 'GCF_000000001.1_ASM1v1')
        self.genome('GCF', ('000', '000', '009'), 'GCF_000000009.1_OLD')
        rows = self.run_list(self.selection('GCF_000000001.1'))
        self.assertEqual([row[0] for row in rows], ['GCF_000000001.1'])

    def test_a_superseded_version_does_not_match(self):
        # accessions carry their version, so a genome revised since the tree was
        # built is not written under the version the release asks for
        self.genome('GCF', ('000', '000', '001'), 'GCF_000000001.1_ASM1v1')
        rows = self.run_list(self.selection('GCF_000000001.2'))
        self.assertEqual(rows, [])

    def test_an_assembly_name_with_underscores_still_yields_the_accession(self):
        held = self.genome('GCA', ('000', '001', '405'), 'GCA_000001405.28_my_odd_name_v2')
        rows = self.run_list(self.selection('GCA_000001405.28'))
        self.assertEqual(rows, [['GCA_000001405.28', held, 'G000001405']])

    def test_a_genome_of_the_release_with_no_directory_is_counted(self):
        # not a verification -- that is ncbi_genome_sync --verify -- but a file
        # quietly short of the release it describes must not pass unremarked
        self.genome('GCF', ('000', '000', '001'), 'GCF_000000001.1_ASM1v1')
        rows = self.run_list(self.selection('GCF_000000001.1', 'GCF_000000002.1'))
        self.assertEqual(len(rows), 1)
        self.assertTrue(any('1 genomes of the release have no directory' in w
                            for w in self.warnings), self.warnings)

    def test_a_complete_tree_warns_about_nothing(self):
        self.genome('GCF', ('000', '000', '001'), 'GCF_000000001.1_ASM1v1')
        self.run_list(self.selection('GCF_000000001.1'))
        self.assertEqual(self.warnings, [])

    def test_no_missing_or_extra_report_is_written_any_more(self):
        # verification moved to ncbi_genome_sync
        self.genome('GCF', ('000', '000', '001'), 'GCF_000000001.1_ASM1v1')
        self.run_list(self.selection('GCF_000000001.1', 'GCF_000000002.1'))
        for suffix in ('-missing', '-extra'):
            self.assertFalse(os.path.exists(self.path('genome_dirs.tsv' + suffix)), suffix)

    def test_paths_written_are_absolute(self):
        # downstream commands open them from their own working directory
        self.genome('GCF', ('000', '000', '001'), 'GCF_000000001.1_ASM1v1')
        rows = self.run_list(self.selection('GCF_000000001.1'))
        self.assertTrue(os.path.isabs(rows[0][1]), rows[0][1])
