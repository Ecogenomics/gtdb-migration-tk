#!/usr/bin/env python3
"""Offline unit tests for trans_table.py -- gTranslate itself is never run.

Run with the interpreter that has tqdm:

    /opt/centos7/sw/miniconda3/envs/python3/bin/python -m unittest -v tests.test_trans_table

What is tested here is everything decided before the subprocess starts: which file
of a genome directory is handed over, which genomes are left out, and the command
line built from the options. Running gTranslate is gTranslate's business; getting
the wrong genomes to it, or the right ones under the wrong names, is this module's.
"""

import gzip
import os
import shutil
import tempfile
import unittest

from gtdb_migration_tk import trans_table as G


class TempDirCase(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='trans_table_test.')

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def genome_dir(self, assembly, fasta=True, empty=False):
        """A genome directory as a release holds it, named for its assembly."""
        path = os.path.join(self.dir, assembly)
        os.makedirs(path)
        if fasta:
            with open(os.path.join(path, assembly + '_genomic.fna.gz'), 'wb') as handle:
                if not empty:
                    handle.write(b'>contig\nACGT\n')
        return path


# ------------------------------------------------------------- naming the FASTA

class GenomicFastaTests(TempDirCase):
    """The genes of a genome end in the same suffix as the genome itself."""

    def test_fasta_is_named_for_the_assembly_directory(self):
        path = G.genomic_fasta('/rel/genbank/GCA/047/639/395/GCA_047639395.1_ASM4763939v1')
        self.assertEqual(os.path.basename(path),
                         'GCA_047639395.1_ASM4763939v1_genomic.fna.gz')

    def test_trailing_separator_does_not_change_the_name(self):
        without = G.genomic_fasta('/rel/refseq/GCF/000/001/405/GCF_000001405.39_GRCh38.p13')
        with_sep = G.genomic_fasta('/rel/refseq/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/')
        self.assertEqual(without, with_sep)


# ------------------------------------------------------------- reading the release

class ReadGenomeDirsTests(TempDirCase):
    """The genome_dirs file is the lingua franca; extra columns may be appended."""

    def write_genome_dirs(self, text):
        path = os.path.join(self.dir, 'genome_dirs.tsv')
        with open(path, 'w') as handle:
            handle.write(text)
        return path

    def test_accession_and_directory_are_read(self):
        path = self.write_genome_dirs(
            'GCF_000001405.39\t/rel/refseq/GCF/000/001/405/x\tG000001405\n')
        self.assertEqual(G.read_genome_dirs(path),
                         [('GCF_000001405.39', '/rel/refseq/GCF/000/001/405/x')])

    def test_further_columns_are_ignored(self):
        path = self.write_genome_dirs('GCA_1.1\t/rel/x\tG1\tsomething\telse\n')
        self.assertEqual(G.read_genome_dirs(path), [('GCA_1.1', '/rel/x')])

    def test_blank_lines_are_skipped(self):
        path = self.write_genome_dirs('GCA_1.1\t/rel/x\tG1\n\n')
        self.assertEqual(len(G.read_genome_dirs(path)), 1)


# ------------------------------------------------------------- what is asked about

class BatchfileRowsTests(TempDirCase):
    """A genome with no FASTA must be named here, not fail mid-run."""

    def test_genome_with_a_fasta_is_included(self):
        path = self.genome_dir('GCF_000001405.39_GRCh38.p13')
        rows, missing = G.batchfile_rows([('GCF_000001405.39', path)])
        self.assertEqual(missing, [])
        self.assertEqual(rows[0][1], 'GCF_000001405.39')
        self.assertTrue(rows[0][0].endswith('_genomic.fna.gz'))

    def test_genome_without_a_fasta_is_left_out_and_named(self):
        path = self.genome_dir('GCA_1.1_ASM1', fasta=False)
        rows, missing = G.batchfile_rows([('GCA_1.1', path)])
        self.assertEqual(rows, [])
        self.assertEqual(missing, ['GCA_1.1'])

    def test_empty_fasta_is_left_out(self):
        path = self.genome_dir('GCA_2.1_ASM2', empty=True)
        rows, missing = G.batchfile_rows([('GCA_2.1', path)])
        self.assertEqual(rows, [])
        self.assertEqual(missing, ['GCA_2.1'])

    def test_the_rest_of_the_release_survives_one_missing_genome(self):
        good = self.genome_dir('GCF_1.1_ASM1')
        bad = self.genome_dir('GCA_2.1_ASM2', fasta=False)
        rows, missing = G.batchfile_rows([('GCF_1.1', good), ('GCA_2.1', bad)])
        self.assertEqual(len(rows), 1)
        self.assertEqual(missing, ['GCA_2.1'])


class WriteBatchfileTests(TempDirCase):
    """gTranslate reads FASTA first, genome ID second; the ID is the accession."""

    def test_columns_are_fasta_then_accession(self):
        batchfile = os.path.join(self.dir, 'batch.tsv')
        G.write_batchfile([('/rel/x/GCF_1.1_ASM1_genomic.fna.gz', 'GCF_1.1')], batchfile)
        with open(batchfile) as handle:
            self.assertEqual(handle.read(),
                             '/rel/x/GCF_1.1_ASM1_genomic.fna.gz\tGCF_1.1\n')


# ------------------------------------------------------------- the command line

class DetectTableCommandTests(unittest.TestCase):
    """An option not given is left off, so gTranslate's defaults stay its own."""

    def test_batchfile_is_used_rather_than_genome_dir(self):
        cmd = G.detect_table_command('/out/batch.tsv', '/out')
        self.assertIn('--batchfile', cmd)
        self.assertNotIn('--genome_dir', cmd)

    def test_cpus_is_passed_as_a_string(self):
        cmd = G.detect_table_command('/out/batch.tsv', '/out', cpus=16)
        self.assertEqual(cmd[cmd.index('--cpus') + 1], '16')

    def test_unset_options_are_absent(self):
        cmd = G.detect_table_command('/out/batch.tsv', '/out')
        for flag in ('--tmpdir', '--prefix', '--custom_model_path',
                     '--force', '--keep_called_genes'):
            self.assertNotIn(flag, cmd)

    def test_flags_appear_only_when_asked_for(self):
        cmd = G.detect_table_command('/out/batch.tsv', '/out',
                                     tmp_dir='/scratch', force=True,
                                     keep_called_genes=True, prefix='r232',
                                     custom_model_path='/srv/db/gtranslate/models')
        self.assertEqual(cmd[cmd.index('--tmpdir') + 1], '/scratch')
        self.assertEqual(cmd[cmd.index('--prefix') + 1], 'r232')
        self.assertEqual(cmd[cmd.index('--custom_model_path') + 1],
                         '/srv/db/gtranslate/models')
        self.assertIn('--force', cmd)
        self.assertIn('--keep_called_genes', cmd)

    def test_subcommand_is_detect_table(self):
        cmd = G.detect_table_command('/out/batch.tsv', '/out')
        self.assertEqual(cmd[:2], [G.GTRANSLATE_BIN, 'detect_table'])


# ------------------------------------------------------------- the manager itself

class ManagerTests(unittest.TestCase):
    """gtranslate and prodigal are checked for when the manager is built."""

    def setUp(self):
        # the real check exits the process, and neither tool is wanted offline
        self._check, G.check_dependencies = G.check_dependencies, lambda *a, **k: True

    def tearDown(self):
        G.check_dependencies = self._check

    def test_both_third_party_tools_are_checked_for(self):
        """gTranslate calls Prodigal, and its package does not depend on it."""
        asked = []
        G.check_dependencies = lambda programs, *a, **k: asked.extend(programs)
        G.GTranslate()
        self.assertEqual(sorted(asked), ['gtranslate', 'prodigal'])

    def test_batch_size_default_matches_the_command_line(self):
        """The two defaults are written out separately and must not drift apart."""
        self.assertEqual(G.GTranslate().batch_size, G.DEFAULT_BATCH_SIZE)
        self.assertEqual(G.DEFAULT_BATCH_SIZE, 10000)


if __name__ == '__main__':
    unittest.main()


# ------------------------------------------------------------- planning the batches

class PlanBatchesTests(TempDirCase):
    """The plan is what several machines agree on, so it must not move."""

    def rows(self, *accessions):
        return [('/rel/x/{}_ASM1_genomic.fna.gz'.format(a), a) for a in accessions]

    def test_release_is_cut_into_batches_of_the_given_size(self):
        batches = G.create_batches(self.rows(*['GCF_{}.1'.format(i) for i in range(5)]),
                                   2, self.dir)
        self.assertEqual(len(batches), 3)
        self.assertEqual(len(G.read_batchfile(
            os.path.join(batches[0], G.BATCHFILE_NAME))), 2)
        self.assertEqual(len(G.read_batchfile(
            os.path.join(batches[-1], G.BATCHFILE_NAME))), 1)

    def test_batches_are_numbered_in_order(self):
        batches = G.create_batches(self.rows('GCF_1.1', 'GCF_2.1'), 1, self.dir)
        self.assertEqual([os.path.basename(b) for b in batches],
                         ['batch_000001', 'batch_000002'])

    def test_an_existing_plan_is_found_and_reused(self):
        G.create_batches(self.rows('GCF_1.1', 'GCF_2.1'), 1, self.dir)
        self.assertEqual(len(G.batch_dir_names(self.dir)), 2)

    def test_a_directory_without_a_batchfile_is_not_a_batch(self):
        os.makedirs(os.path.join(self.dir, 'batch_000001'))
        self.assertEqual(G.batch_dir_names(self.dir), [])

    def test_batchfile_survives_a_round_trip(self):
        batches = G.create_batches(self.rows('GCF_1.1'), 10, self.dir)
        self.assertEqual(G.read_batchfile(os.path.join(batches[0], G.BATCHFILE_NAME)),
                         self.rows('GCF_1.1'))


# ------------------------------------------------------------- the canaries

class ClaimTests(TempDirCase):
    """Two machines must never both take one batch, and a reset must not lose one."""

    def batch(self):
        path = os.path.join(self.dir, 'batch_000001')
        os.makedirs(path)
        return path

    def test_an_unclaimed_batch_is_claimed(self):
        batch = self.batch()
        self.assertTrue(G.claim_batch(batch))
        self.assertEqual(G.batch_state(batch), G.STATE_RUNNING)

    def test_a_batch_claimed_by_a_live_process_is_not_taken(self):
        batch = self.batch()
        G.claim_batch(batch)
        self.assertFalse(G.claim_batch(batch))

    def test_a_claim_of_a_dead_process_on_this_host_is_reclaimed(self):
        """What a machine reset leaves behind; nothing else would ever run it."""
        batch = self.batch()
        with open(os.path.join(batch, G.RUNNING_CANARY), 'w') as handle:
            handle.write(G.canary_payload())
        # a PID that cannot be running, recorded against this host
        text = open(os.path.join(batch, G.RUNNING_CANARY)).read()
        with open(os.path.join(batch, G.RUNNING_CANARY), 'w') as handle:
            handle.write(text.replace('pid\t{}'.format(os.getpid()), 'pid\t2147483646'))
        self.assertTrue(G.claim_batch(batch))

    def test_a_claim_of_another_host_is_left_alone(self):
        batch = self.batch()
        with open(os.path.join(batch, G.RUNNING_CANARY), 'w') as handle:
            handle.write('host\tsome-other-machine\npid\t1\ntime\tnow\n')
        self.assertFalse(G.claim_batch(batch))

    def test_reclaim_takes_another_hosts_claim(self):
        batch = self.batch()
        with open(os.path.join(batch, G.RUNNING_CANARY), 'w') as handle:
            handle.write('host\tsome-other-machine\npid\t1\ntime\tnow\n')
        self.assertTrue(G.claim_batch(batch, reclaim=True))

    def test_a_finished_batch_reports_success_and_gives_up_the_claim(self):
        batch = self.batch()
        G.claim_batch(batch)
        G.finish_batch(batch, compared=7)
        self.assertEqual(G.batch_state(batch), G.STATE_SUCCESS)
        self.assertFalse(os.path.exists(os.path.join(batch, G.RUNNING_CANARY)))
        self.assertEqual(G.read_canary(os.path.join(batch, G.SUCCESS_CANARY))['compared'], '7')

    def test_a_failed_batch_gives_up_the_claim_so_it_is_retried(self):
        batch = self.batch()
        G.claim_batch(batch)
        G.fail_batch(batch, 'gtranslate returned exit code 1.')
        self.assertEqual(G.batch_state(batch), G.STATE_FAILED)
        self.assertTrue(G.claim_batch(batch))

    def test_claiming_a_failed_batch_clears_the_failure(self):
        batch = self.batch()
        G.fail_batch(batch, 'whatever')
        G.claim_batch(batch)
        self.assertFalse(os.path.exists(os.path.join(batch, G.FAILED_CANARY)))

    def test_success_outranks_running(self):
        batch = self.batch()
        G.claim_batch(batch)
        open(os.path.join(batch, G.SUCCESS_CANARY), 'w').close()
        self.assertEqual(G.batch_state(batch), G.STATE_SUCCESS)


# ------------------------------------------------------------- the comparison

class ComparisonTests(TempDirCase):
    """A conflict is a genome GTDB would call genes for under a table NCBI rejects."""

    def genome_dir(self, accession, ncbi_table=None):
        assembly = '{}_ASM1'.format(accession)
        path = os.path.join(self.dir, assembly)
        os.makedirs(path)
        if ncbi_table is not None:
            with gzip.open(os.path.join(path, assembly + '_genomic.gff.gz'), 'wt') as handle:
                handle.write('##gff-version 3\n')
                handle.write('c\tRefSeq\tCDS\t1\t9\t.\t+\t0\t'
                             'ID=cds1;product=x;transl_table={}\n'.format(ncbi_table))
        return path

    def summary(self, *rows):
        path = os.path.join(self.dir, 'gtranslate.translation_table_summary.tsv')
        with open(path, 'w') as handle:
            handle.write('user_genome\tbest_tln_table\tcoding_density_4\t'
                         'coding_density_11\tconfidence\n')
            for accession, table in rows:
                handle.write('{}\t{}\t90.1\t64.2\t1.0\n'.format(accession, table))
        return path

    def test_matching_tables_agree(self):
        path = self.genome_dir('GCF_1.1', ncbi_table=11)
        rows, _ = G.comparison_rows(G.read_translation_table_summary(self.summary(('GCF_1.1', '11'))),
                                    {'GCF_1.1': path}, {})
        self.assertEqual(rows[0][G.COMPARISON_HEADER.index('result')], G.AGREE)

    def test_differing_tables_conflict(self):
        path = self.genome_dir('GCF_1.1', ncbi_table=11)
        rows, _ = G.comparison_rows(G.read_translation_table_summary(self.summary(('GCF_1.1', '4'))),
                                    {'GCF_1.1': path}, {})
        self.assertEqual(rows[0][G.COMPARISON_HEADER.index('result')], G.CONFLICT)

    def test_genome_ncbi_declares_no_table_for_is_left_out(self):
        path = self.genome_dir('GCF_1.1')
        rows, no_table = G.comparison_rows(
            G.read_translation_table_summary(self.summary(('GCF_1.1', '11'))), {'GCF_1.1': path}, {})
        self.assertEqual(rows, [])
        self.assertEqual(no_table, 1)

    def test_coding_densities_and_lineage_are_carried(self):
        path = self.genome_dir('GCF_1.1', ncbi_table=11)
        rows, _ = G.comparison_rows(
            G.read_translation_table_summary(self.summary(('GCF_1.1', '11'))),
            {'GCF_1.1': path},
            {'GCF_1.1': 'd__Bacteria;p__Pseudomonadota;s__Escherichia coli'})
        self.assertEqual(rows[0][G.COMPARISON_HEADER.index('coding_density_4')], '90.1')
        self.assertEqual(rows[0][G.COMPARISON_HEADER.index('coding_density_11')], '64.2')
        self.assertIn('d__Bacteria', rows[0][G.COMPARISON_HEADER.index('ncbi_taxonomy')])

    def test_lineage_is_found_through_the_canonical_accession(self):
        """A GenBank genome takes the lineage held against its RefSeq counterpart."""
        path = self.genome_dir('GCA_005435135.1', ncbi_table=11)
        rows, _ = G.comparison_rows(
            G.read_translation_table_summary(self.summary(('GCA_005435135.1', '11'))),
            {'GCA_005435135.1': path},
            G.read_taxonomy(self.write_taxonomy('GCF_005435135.1\td__Bacteria;s__X\n')))
        self.assertEqual(rows[0][G.COMPARISON_HEADER.index('ncbi_taxonomy')], 'd__Bacteria;s__X')

    def test_genome_missing_from_the_taxonomy_is_still_compared(self):
        path = self.genome_dir('GCF_1.1', ncbi_table=4)
        rows, _ = G.comparison_rows(G.read_translation_table_summary(self.summary(('GCF_1.1', '4'))),
                                    {'GCF_1.1': path}, {})
        self.assertEqual(rows[0][G.COMPARISON_HEADER.index('ncbi_taxonomy')], 'na')

    def write_taxonomy(self, text):
        path = os.path.join(self.dir, 'taxonomy.tsv')
        with open(path, 'w') as handle:
            handle.write(text)
        return path

    def test_checkm_table_is_reported_beside_the_others(self):
        """The density rule alone cannot express 25, which is the point of it."""
        path = self.genome_dir('GCF_1.1', ncbi_table=11)
        rows, _ = G.comparison_rows(
            G.read_translation_table_summary(self.summary(('GCF_1.1', '25'))),
            {'GCF_1.1': path}, {})
        row = rows[0]
        self.assertEqual(row[G.COMPARISON_HEADER.index('gtranslate_tt')], '25')
        self.assertEqual(row[G.COMPARISON_HEADER.index('checkm_tt')], '4')

    def test_checkm_table_is_11_where_the_densities_are_close(self):
        path = self.genome_dir('GCF_1.1', ncbi_table=11)
        summary = os.path.join(self.dir, 'close.tsv')
        with open(summary, 'w') as handle:
            handle.write('user_genome\tbest_tln_table\tcoding_density_4\t'
                         'coding_density_11\n')
            handle.write('GCF_1.1\t11\t86.24689\t86.64953\n')
        rows, _ = G.comparison_rows(
            G.read_translation_table_summary(summary), {'GCF_1.1': path}, {})
        self.assertEqual(rows[0][G.COMPARISON_HEADER.index('checkm_tt')], '11')


class AggregationTests(TempDirCase):
    """The release file is the whole release or absent, never a part of it."""

    def comparison(self, name, *rows):
        path = os.path.join(self.dir, name)
        G.write_comparison(rows, path)
        return path

    def test_headers_are_not_repeated(self):
        first = self.comparison('a.tsv', ('GCF_1.1', '11', '11', 'agree', '', '', 'na'))
        second = self.comparison('b.tsv', ('GCF_2.1', '4', '11', 'conflict', '', '', 'na'))
        out = os.path.join(self.dir, 'all.tsv')
        self.assertEqual(G.concatenate([first, second], out), 2)
        with open(out) as handle:
            lines = handle.read().splitlines()
        self.assertEqual(len(lines), 3)
        self.assertEqual(lines[0].split('\t')[0], 'genome_id')
