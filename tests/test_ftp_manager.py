#!/usr/bin/env python3
"""Offline unit tests for ftp_manager.py -- no FTP site, no genome files.

Run with the interpreter that has tqdm:

    /opt/centos7/sw/miniconda3/envs/python3/bin/python -m unittest -v tests.test_ftp_manager

What is tested here is the decision about which genomes belong in a release. Reading
the NCBI assembly summary files is ncbi_utils.py's job and is tested in
tests/test_ncbi_utils.py; one test below still feeds this module a revised column
layout, to confirm the selection is reading those tables by column name.
"""

import os
import shutil
import tempfile
import unittest

from gtdb_migration_tk import ftp_manager as F


# an NCBI assembly summary: a comment block, the header, then the genomes
HEADER = ('#assembly_accession\tbioproject\tversion_status\tgbrs_paired_asm'
          '\texcluded_from_refseq\tftp_path')
COMMENT = '#   See ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt'


def summary(*rows, header=HEADER, comment=COMMENT):
    return '\n'.join([comment, header] + list(rows)) + '\n'


class TempDirCase(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='ftp_manager_test.')

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def write(self, name, text):
        path = os.path.join(self.dir, name)
        with open(path, 'w') as handle:
            handle.write(text)
        return path


# ------------------------------------------------------------------ release comparison

class GenericDatabaseManagerTests(TempDirCase):
    def setUp(self):
        super().setUp()
        self.manager = F.GenericDatabaseManager()
        self.new_genomes = {'G1': '/ftp/g1', 'G2': '/ftp/g2', 'G3': '/ftp/g3_new'}
        self.old_genomes = {'G2': '/gtdb/g2', 'G3': '/gtdb/g3_old', 'G4': '/gtdb/g4'}

    def test_genomes_to_remove_are_those_gone_from_ncbi(self):
        self.assertEqual(self.manager.generate_genomes_to_remove(self.new_genomes, self.old_genomes),
                         {'G4': '/gtdb/g4'})

    def test_genomes_to_add_are_those_new_to_ncbi(self):
        self.assertEqual(self.manager.generate_genomes_to_add(self.new_genomes, self.old_genomes),
                         {'G1': '/ftp/g1'})

    def test_genomes_to_compare_are_those_in_both_releases(self):
        # G3 is in both, but at a different path, so it must be compared and not ignored
        self.assertEqual(sorted(self.manager.generate_genomes_to_compare(self.new_genomes, self.old_genomes)),
                         ['G2', 'G3'])

    def test_empty_releases_are_handled(self):
        self.assertEqual(self.manager.generate_genomes_to_remove({}, {}), {})
        self.assertEqual(self.manager.generate_genomes_to_add({}, {}), {})
        self.assertEqual(self.manager.generate_genomes_to_compare({}, {}), [])

    def test_load_previous_records(self):
        path = self.write('old_dirs.tsv', 'GCF_000001405.40\t/gtdb/g1\nGCA_000002305.1\t/gtdb/g2\n')
        self.assertEqual(self.manager.load_previous_records(path),
                         {'GCF_000001405.40': '/gtdb/g1', 'GCA_000002305.1': '/gtdb/g2'})

    def test_load_ftp_records_keeps_only_selected_genomes_of_the_right_database(self):
        path = self.write('ftp_dirs.tsv',
                          'GCF_000000001.1\t/ftp/keep\n'       # selected
                          'GCF_000000002.1\t/ftp/drop\n'       # not selected
                          'GCA_000000003.1\t/ftp/wrong_db\n')  # GenBank, not RefSeq
        selected = {'GCF_000000001.1', 'GCA_000000003.1'}
        self.assertEqual(self.manager.load_ftp_records(path, 'GCF', selected),
                         {'GCF_000000001.1': '/ftp/keep'})


class RefSeqManagerTests(TempDirCase):
    def test_only_latest_assemblies_are_selected(self):
        path = self.write('summary.txt',
                          summary('GCF_000000001.1\tPRJNA1\tlatest\tGCA_000000001.1\tna\tftp://a',
                                  'GCF_000000002.1\tPRJNA2\treplaced\tGCA_000000002.1\tna\tftp://b',
                                  'GCF_000000003.1\tPRJNA3\tsuppressed\tGCA_000000003.1\tna\tftp://c'))
        manager = F.RefSeqManager(self.dir)
        self.assertEqual(manager.parse_assembly_summary(path), ['GCF_000000001.1'])

    def test_selection_survives_a_revised_column_layout(self):
        header = '#assembly_accession\tbioproject\tbiosample\tversion_status'
        path = self.write('summary.txt',
                          summary('GCF_000000001.1\tPRJNA1\tSAMN1\tlatest',
                                  'GCF_000000002.1\tPRJNA2\tSAMN2\treplaced',
                                  header=header))
        manager = F.RefSeqManager(self.dir)
        self.assertEqual(manager.parse_assembly_summary(path), ['GCF_000000001.1'])

    def test_construction_writes_nothing(self):
        # reports are opened by run_comparison, so a manager can be built without touching disk
        F.RefSeqManager(os.path.join(self.dir, 'does_not_exist'))
        F.GenBankManager(os.path.join(self.dir, 'does_not_exist'))
        self.assertEqual(os.listdir(self.dir), [])


class RefSeqDirectoryCompleteness(TempDirCase):
    """A GenBank genome is rescued when its RefSeq counterpart holds no assembly."""

    def refseq_dir(self, *file_names):
        genome_dir = os.path.join(self.dir, 'GCF_000000009.1_ASM9v1')
        os.mkdir(genome_dir)
        for name in file_names:
            open(os.path.join(genome_dir, name), 'w').close()
        return genome_dir

    def select_against(self, genome_dir):
        refseq_dirs = self.write('refseq_dirs.tsv', 'GCF_000000009.1\t{}\n'.format(genome_dir))
        arc = self.write('arc.txt',
                         summary('GCA_000000009.1\tPRJNA9\tlatest\tGCF_000000009.1\tna\tftp://b'))
        bac = self.write('bac.txt', summary())
        manager = F.GenBankManager(self.dir)
        with open(os.path.join(self.dir, 'gca_selection.log'), 'w') as handle:
            manager.select_gca = handle
            return manager.select_genbank_genomes(arc, bac, refseq_dirs)

    def test_complete_refseq_directory_means_the_genbank_copy_is_not_needed(self):
        genome_dir = self.refseq_dir('GCF_000000009.1_ASM9v1_genomic.fna.gz')
        self.assertEqual(self.select_against(genome_dir), [])

    def test_derived_files_alone_do_not_count_as_an_assembly(self):
        # '*_genomic.fna.gz' also matches these two, so a glob reports the directory
        # as complete and the GenBank genome is never rescued
        genome_dir = self.refseq_dir('GCF_000000009.1_ASM9v1_cds_from_genomic.fna.gz',
                                     'GCF_000000009.1_ASM9v1_rna_from_genomic.fna.gz',
                                     'GCF_000000009.1_ASM9v1_protein.faa.gz')
        self.assertEqual(self.select_against(genome_dir), ['GCA_000000009.1'])

    def test_empty_refseq_directory_rescues_the_genbank_genome(self):
        self.assertEqual(self.select_against(self.refseq_dir()), ['GCA_000000009.1'])

    def test_missing_refseq_directory_rescues_the_genbank_genome(self):
        self.assertEqual(self.select_against(os.path.join(self.dir, 'gone')),
                         ['GCA_000000009.1'])


class GenBankManagerTests(TempDirCase):
    def select(self, *rows):
        # a complete RefSeq directory for GCF_000000009.1, so genomes paired with it
        # are covered by RefSeq and need no GenBank copy
        genome_dir = os.path.join(self.dir, 'GCF_000000009.1_ASM9v1')
        os.mkdir(genome_dir)
        open(os.path.join(genome_dir, 'GCF_000000009.1_ASM9v1_genomic.fna.gz'), 'w').close()
        refseq_dirs = self.write('refseq_dirs.tsv',
                                 'GCF_000000009.1\t{}\n'.format(genome_dir))
        arc = self.write('arc.txt', summary(*rows))
        bac = self.write('bac.txt', summary())
        manager = F.GenBankManager(self.dir)
        with open(os.path.join(self.dir, 'gca_selection.log'), 'w') as handle:
            manager.select_gca = handle
            return manager.select_genbank_genomes(arc, bac, refseq_dirs), manager

    def test_genome_without_a_refseq_counterpart_is_selected(self):
        selected, _ = self.select('GCA_000000001.1\tPRJNA1\tlatest\tna\tna\tftp://a')
        self.assertEqual(selected, ['GCA_000000001.1'])

    def test_genome_already_covered_by_refseq_is_skipped(self):
        # G000000009 is in the RefSeq genome directories, so the GenBank copy is redundant
        selected, _ = self.select('GCA_000000009.1\tPRJNA9\tlatest\tGCF_000000009.1\tna\tftp://b')
        self.assertEqual(selected, [])

    def test_surveillance_genomes_are_skipped(self):
        selected, _ = self.select('GCA_000000002.1\tPRJNA2\tlatest\tna\tsurveillance\tftp://c')
        self.assertEqual(selected, [])

    def test_superseded_genomes_are_skipped(self):
        selected, _ = self.select('GCA_000000003.1\tPRJNA3\treplaced\tna\tna\tftp://d')
        self.assertEqual(selected, [])

    def test_domain_is_recorded_for_selected_genomes(self):
        _, manager = self.select('GCA_000000001.1\tPRJNA1\tlatest\tna\tna\tftp://a')
        self.assertEqual(manager.genome_domain_dict['GCA_000000001.1'], F.ARCHAEA)

    def test_populate_genomes_dict_keys_on_the_canonical_accession(self):
        path = self.write('refseq_dirs.tsv', 'GCF_000001405.40\t/gtdb/g1\n')
        manager = F.GenBankManager(self.dir)
        self.assertEqual(manager._populate_genomes_dict(path), {'G000001405': '/gtdb/g1'})

    def test_paired_accessions_reduce_to_the_same_key(self):
        # matching a GenBank genome to its RefSeq counterpart rests on this
        path = self.write('dirs.tsv', 'GCA_000001405.1\t/gtdb/gca\n')
        manager = F.GenBankManager(self.dir)
        self.assertEqual(list(manager._populate_genomes_dict(path)), ['G000001405'])


if __name__ == '__main__':
    unittest.main()
