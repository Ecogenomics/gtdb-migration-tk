#!/usr/bin/env python3
"""Offline unit tests for the checkm and checkm2 commands in checkm_manager.py.

Both assess the genomes --report says are new or whose genomic FASTA changed
-- and no others -- from the proteins the prodigal command called, taken as
correct. Both cut those genomes into batches in --out_dir, by batching.py.

The programs are stood in for by shell scripts first on PATH, which answer for
their versions as CheckM 1.2.5 and CheckM2 1.1.0 do and write the tables each
step would, one row per genome they were handed. CheckM2's names each genome
for its file's basename less '.faa.gz', as CheckM2 does. Either fails where its output directory lies
under the batch named in FAIL_BATCH. So what is tested is what the commands do
around the programs, not any one build of them.
"""

import logging
import os
import shutil
import stat
import tempfile
import unittest

from gtdb_migration_tk import checkm_manager as M
from gtdb_migration_tk.batching import (FAILED_CANARY, SUCCESS_CANARY,
                                        batchfile_path, read_batchfile)
from gtdb_migration_tk.update_genomes import (STATUS_FASTA_CHANGED,
                                              STATUS_FASTA_UNCHANGED, STATUS_NEW,
                                              STATUS_REMOVED,
                                              STATUS_SEQUENCES_UNCHANGED,
                                              STATUS_TO_CURATE)

FAKE_CHECKM2 = r'''#!/bin/sh
if [ "$1" = "--version" ]; then echo 1.1.0; exit 0; fi
echo "$@" >> "$CALLS"
out=; files=; seen_input=
shift
while [ $# -gt 0 ]; do
    if [ -n "$seen_input" ]; then files="$files $1"
    elif [ "$1" = "--output-directory" ]; then shift; out=$1
    elif [ "$1" = "--input" ]; then seen_input=1
    fi
    shift
done
if [ -n "$FAIL_BATCH" ]; then case "$out" in *"/$FAIL_BATCH/"*) exit 1;; esac; fi
mkdir -p "$out"
printf 'Name\tCompleteness\tContamination\n' > "$out/quality_report.tsv"
for f in $files; do
    printf '%s\t100.0\t0.0\n' "$(basename "$f" .faa.gz)" >> "$out/quality_report.tsv"
done
'''

FAKE_CHECKM = r'''#!/bin/sh
if [ "$1" = "-h" ]; then printf '\n                ...::: CheckM v1.2.5 :::...\n\n'; exit 0; fi
echo "$@" >> "$CALLS"
step=$1; shift
prev=; last=; file=; aln=
for a in "$@"; do
    [ "$prev" = "-f" ] && file=$a
    [ "$prev" = "-a" ] && aln=$a
    prev=$a
done
for a in "$@"; do before=$last; last=$a; done
if [ -n "$FAIL_BATCH" ]; then case "$last" in *"/$FAIL_BATCH/"*) exit 1;; esac; fi
case $step in
    lineage_wf)
        mkdir -p "$last"; : > "$last/lineage.ms"
        for f in "$before"/*.faa.gz; do basename "$f" .faa.gz; done > "$last/bins.txt";;
    tree_qa)
        printf 'Bin Id\tTaxonomy\n' > "$file"
        while read b; do printf '%s\tk__Bacteria\n' "$b" >> "$file"; done < "$last/bins.txt";;
    qa)
        printf 'Bin Id\tCompleteness\n' > "$file"
        while read b; do
            printf '%s\t99.0\n' "$b" >> "$file"
            [ -n "$aln" ] && printf '%s\talignment\n' "$b" >> "$aln"
        done < "$last/bins.txt";;
esac
exit 0
'''


class TempDirCase(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='checkm_manager_test.')
        self.addCleanup(shutil.rmtree, self.dir, True)

        bin_dir = os.path.join(self.dir, 'bin')
        os.makedirs(bin_dir)
        for name, script in (('checkm', FAKE_CHECKM), ('checkm2', FAKE_CHECKM2)):
            path = os.path.join(bin_dir, name)
            with open(path, 'w') as handle:
                handle.write(script)
            os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR)

        self.calls = os.path.join(self.dir, 'calls.txt')
        for name, value in (('PATH', bin_dir + os.pathsep + os.environ['PATH']),
                            ('CALLS', self.calls),
                            ('FAIL_BATCH', '')):
            self.addCleanup(self.restore_env, name, os.environ.get(name))
            os.environ[name] = value

        self.out = os.path.join(self.dir, 'out')
        self.tmp = os.path.join(self.dir, 'tmp')
        logging.getLogger('timestamp').addHandler(logging.NullHandler())

    @staticmethod
    def restore_env(name, value):
        if value is None:
            os.environ.pop(name, None)
        else:
            os.environ[name] = value

    def release(self, genomes, missing=()):
        """A genome_dirs file and report.log.

        genomes is {accession: report outcome}; missing the genomes with no
        proteins.
        """

        dirs = os.path.join(self.dir, 'genome_dirs.tsv')
        report = os.path.join(self.dir, 'report.log')
        with open(dirs, 'w') as fdirs, open(report, 'w') as freport:
            for gid, outcome in genomes.items():
                gdir = os.path.join(self.dir, 'genomes', gid + '_ASM1v1')
                freport.write('{}\t{}\n'.format(gid, outcome))
                if outcome == STATUS_REMOVED:
                    continue
                fdirs.write('{}\t{}\tG{}\n'.format(gid, gdir, gid[4:13]))
                if gid in missing:
                    continue
                os.makedirs(os.path.join(gdir, 'prodigal'))
                with open(os.path.join(gdir, gid + '_ASM1v1_genomic.fna.gz'), 'w') as handle:
                    handle.write('>contig\nACGT\n')
                with open(os.path.join(gdir, 'prodigal', gid + '_protein.faa.gz'), 'w') as handle:
                    handle.write('>gene\nMK\n')
        return dirs, report

    def planned(self, layout):
        """The accessions of every batch planned in --out_dir, in batch order."""
        return [accession for name in sorted(os.listdir(self.out)) if name.startswith('batch_')
                for _, accession in read_batchfile(batchfile_path(
                    os.path.join(self.out, name), layout))]

    def command(self, program, dirs, report, cpus=1, batch_size=10, all_genomes=False):
        manager = program(cpus=cpus, tmp_dir=self.tmp, batch_size=batch_size,
                          heartbeat=3600)
        return manager.run(dirs, report, self.out, all_genomes)

    def table(self, name):
        with open(os.path.join(self.out, name)) as handle:
            header = handle.readline().rstrip('\n').split('\t')
            rows = [line.rstrip('\n').split('\t') for line in handle]
        return header, rows

    def first_column(self, name):
        return [row[0] for row in self.table(name)[1]]

    def not_assessed(self, name):
        header, rows = self.table(name)
        self.assertEqual(header, list(M.NOT_ASSESSED_HEADER))
        return {row[0]: row[1] for row in rows}

    def calls_of(self, step):
        if not os.path.exists(self.calls):
            return []
        with open(self.calls) as handle:
            return [line.split() for line in handle if line.split()[0] == step]

    def batch(self, number):
        return os.path.join(self.out, 'batch_{:06d}'.format(number))


class ChoosingTheGenomes(TempDirCase):
    """Only the genomes --report says are new or changed are planned and assessed."""

    OUTCOMES = {'GCA_000000001.1': STATUS_NEW,
                'GCA_000000002.1': STATUS_FASTA_UNCHANGED,
                'GCA_000000003.1': STATUS_FASTA_CHANGED,
                'GCA_000000004.1': STATUS_SEQUENCES_UNCHANGED,
                'GCA_000000005.1': STATUS_REMOVED,
                'GCA_000000006.1': STATUS_TO_CURATE + ';ValueError: bad'}

    def test_checkm_plans_and_assesses_only_the_new_and_changed_genomes(self):
        dirs, report = self.release(self.OUTCOMES)
        self.command(M.CheckM, dirs, report)

        self.assertEqual(self.planned(M.CHECKM_LAYOUT), ['GCA_000000001.1', 'GCA_000000003.1'])
        [args] = self.calls_of('lineage_wf')
        with open(os.path.join(self.batch(1), 'checkm', 'bins.txt')) as handle:
            self.assertEqual(handle.read().split(),
                             ['GCA_000000001.1_protein', 'GCA_000000003.1_protein'])
        self.assertEqual(self.first_column('checkm.profiles.tsv'),
                         ['GCA_000000001.1_protein', 'GCA_000000003.1_protein'])

    def test_checkm2_plans_and_assesses_only_the_new_and_changed_genomes(self):
        dirs, report = self.release(self.OUTCOMES)
        self.command(M.CheckM2, dirs, report)

        self.assertEqual(self.planned(M.CHECKM2_LAYOUT), ['GCA_000000001.1', 'GCA_000000003.1'])
        [args] = self.calls_of('predict')
        self.assertEqual([os.path.basename(arg) for arg in args[args.index('--input') + 1:]],
                         ['GCA_000000001.1.faa.gz', 'GCA_000000003.1.faa.gz'])
        self.assertEqual(self.first_column(M.CHECKM2_RELEASE_REPORT),
                         ['GCA_000000001.1', 'GCA_000000003.1'])

    def test_all_plans_every_genome_the_release_holds(self):
        dirs, report = self.release(self.OUTCOMES)
        self.command(M.CheckM2, dirs, report, all_genomes=True)

        self.assertEqual(self.planned(M.CHECKM2_LAYOUT),
                         ['GCA_000000001.1', 'GCA_000000002.1', 'GCA_000000003.1',
                          'GCA_000000004.1'])

    def test_a_report_of_none_plans_every_genome_of_the_genome_dirs_file(self):
        dirs, _ = self.release(self.OUTCOMES)
        self.command(M.CheckM2, dirs, 'none')

        self.assertEqual(self.planned(M.CHECKM2_LAYOUT),
                         ['GCA_000000001.1', 'GCA_000000002.1', 'GCA_000000003.1',
                          'GCA_000000004.1', 'GCA_000000006.1'])


class RunningCheckM2(TempDirCase):
    def run_checkm2(self, *args, **kwargs):
        return self.command(M.CheckM2, *args, **kwargs)

    def test_only_genomes_to_regenerate_are_assessed_and_named_by_accession(self):
        dirs, report = self.release({'GCA_000000001.1': STATUS_NEW,
                                     'GCA_000000002.1': STATUS_FASTA_UNCHANGED,
                                     'GCF_000000003.1': STATUS_NEW})
        self.assertTrue(self.run_checkm2(dirs, report))

        self.assertEqual(self.first_column(M.CHECKM2_RELEASE_REPORT),
                         ['GCA_000000001.1', 'GCF_000000003.1'])

    def test_all_genomes_assesses_every_genome_the_release_holds(self):
        dirs, report = self.release({'GCA_000000001.1': STATUS_NEW,
                                     'GCA_000000002.1': STATUS_FASTA_UNCHANGED,
                                     'GCA_000000004.1': STATUS_REMOVED})
        self.run_checkm2(dirs, report, all_genomes=True)

        self.assertEqual(self.first_column(M.CHECKM2_RELEASE_REPORT),
                         ['GCA_000000001.1', 'GCA_000000002.1'])

    def test_checkm2_is_handed_the_prodigal_proteins_the_threads_and_the_tmp_dir(self):
        dirs, report = self.release({'GCA_000000001.1': STATUS_NEW})
        self.run_checkm2(dirs, report, cpus=7)

        [args] = self.calls_of('predict')
        self.assertIn('--genes', args)
        self.assertNotIn('--ttable', args)
        self.assertEqual(args[args.index('--threads') + 1], '7')
        self.assertEqual(args[args.index('--tmpdir') + 1], self.tmp)
        self.assertTrue(args[-1].endswith('GCA_000000001.1.faa.gz'))
        [(planned, _)] = read_batchfile(batchfile_path(
            self.batch(1), M.CHECKM2_LAYOUT))
        self.assertEqual(planned, os.path.join(
            self.dir, 'genomes', 'GCA_000000001.1_ASM1v1', 'prodigal',
            'GCA_000000001.1_protein.faa.gz'))

    def test_a_genome_without_proteins_is_named_rather_than_failing_its_batch(self):
        dirs, report = self.release({'GCA_000000001.1': STATUS_NEW,
                                     'GCA_000000002.1': STATUS_NEW},
                                    missing={'GCA_000000001.1'})
        self.run_checkm2(dirs, report)

        self.assertEqual(self.first_column(M.CHECKM2_RELEASE_REPORT), ['GCA_000000002.1'])
        self.assertEqual(self.not_assessed(M.CHECKM2_NOT_ASSESSED),
                         {'GCA_000000001.1': M.REASON_NO_PROTEINS})

    def test_the_version_is_written_beside_what_it_made(self):
        dirs, report = self.release({'GCA_000000001.1': STATUS_NEW})
        self.run_checkm2(dirs, report)

        with open(os.path.join(self.batch(1), 'checkm2.version')) as handle:
            self.assertEqual(handle.read().strip(), '1.1.0')

    def test_the_genomes_are_cut_into_batches_and_merged_in_order(self):
        genomes = {'GCA_%09d.1' % i: STATUS_NEW for i in range(5)}
        dirs, report = self.release(genomes)
        self.run_checkm2(dirs, report, batch_size=2)

        for number in (1, 2, 3):
            self.assertTrue(os.path.exists(
                os.path.join(self.batch(number), SUCCESS_CANARY)))
        self.assertEqual(len(self.calls_of('predict')), 3)
        self.assertEqual(self.first_column(M.CHECKM2_RELEASE_REPORT), sorted(genomes))


class RestartingCheckM2(TempDirCase):
    def run_checkm2(self, *args, **kwargs):
        return self.command(M.CheckM2, *args, **kwargs)

    def test_a_failed_batch_writes_no_release_file_and_is_repeated_next_run(self):
        genomes = {'GCA_%09d.1' % i: STATUS_NEW for i in range(4)}
        dirs, report = self.release(genomes)

        os.environ['FAIL_BATCH'] = 'batch_000002'
        self.assertFalse(self.run_checkm2(dirs, report, batch_size=2))
        self.assertTrue(os.path.exists(
            os.path.join(self.batch(2), FAILED_CANARY)))
        self.assertFalse(os.path.exists(os.path.join(self.out, M.CHECKM2_RELEASE_REPORT)))

        os.environ['FAIL_BATCH'] = ''
        self.assertTrue(self.run_checkm2(dirs, report, batch_size=2))

        # batch 1 again would be four calls: only the failed batch is repeated
        self.assertEqual(len(self.calls_of('predict')), 3)
        self.assertEqual(self.first_column(M.CHECKM2_RELEASE_REPORT), sorted(genomes))

    def test_a_finished_run_repeated_runs_checkm2_on_nothing(self):
        dirs, report = self.release({'GCA_000000001.1': STATUS_NEW})
        self.run_checkm2(dirs, report)
        self.run_checkm2(dirs, report)

        self.assertEqual(len(self.calls_of('predict')), 1)
        self.assertEqual(self.first_column(M.CHECKM2_RELEASE_REPORT), ['GCA_000000001.1'])


class RunningCheckM(TempDirCase):
    def run_checkm(self, *args, **kwargs):
        return self.command(M.CheckM, *args, **kwargs)

    def test_checkm_is_handed_the_proteins_and_not_left_to_call_genes(self):
        dirs, report = self.release({'GCA_000000001.1': STATUS_NEW})
        self.assertTrue(self.run_checkm(dirs, report, cpus=3))

        [args] = self.calls_of('lineage_wf')
        self.assertIn('--genes', args)
        self.assertEqual(args[args.index('-x') + 1], 'faa.gz')
        self.assertEqual(args[args.index('-t') + 1], '3')
        self.assertEqual(args[args.index('--tmpdir') + 1], self.tmp)

    def test_the_profile_joins_qa_and_tree_qa_for_each_genome(self):
        dirs, report = self.release({'GCA_000000001.1': STATUS_NEW,
                                     'GCA_000000002.1': STATUS_FASTA_UNCHANGED,
                                     'GCA_000000003.1': STATUS_NEW})
        self.run_checkm(dirs, report)

        header, rows = self.table('checkm.profiles.tsv')
        self.assertEqual(header, ['Bin Id', 'Completeness', 'Taxonomy'])
        self.assertEqual(rows, [['GCA_000000001.1_protein', '99.0', 'k__Bacteria'],
                                ['GCA_000000003.1_protein', '99.0', 'k__Bacteria']])
        self.assertEqual(self.first_column('checkm.qa_sh100.tsv'),
                         ['GCA_000000001.1_protein', 'GCA_000000003.1_protein'])
        with open(os.path.join(self.out, M.CHECKM_RELEASE_ALIGNMENT)) as handle:
            self.assertEqual(len(handle.readlines()), 2)

    def test_the_version_is_written_beside_each_batch_it_made(self):
        dirs, report = self.release({'GCA_000000001.1': STATUS_NEW})
        self.run_checkm(dirs, report)

        with open(os.path.join(self.batch(1), 'checkm.version')) as handle:
            self.assertEqual(handle.read().strip(), 'CheckM v1.2.5')

    def test_a_checkm_that_fails_fails_its_batch_and_writes_no_release_file(self):
        dirs, report = self.release({'GCA_000000001.1': STATUS_NEW})
        os.environ['FAIL_BATCH'] = 'batch_000001'
        self.assertFalse(self.run_checkm(dirs, report))

        self.assertTrue(os.path.exists(
            os.path.join(self.batch(1), FAILED_CANARY)))
        self.assertFalse(os.path.exists(os.path.join(self.out, 'checkm.profiles.tsv')))

    def test_the_release_tables_are_the_batches_concatenated_under_one_header(self):
        genomes = {'GCA_%09d.1' % i: STATUS_NEW for i in range(5)}
        dirs, report = self.release(genomes)
        self.run_checkm(dirs, report, batch_size=2)

        proteins = [gid + '_protein' for gid in sorted(genomes)]
        header, rows = self.table('checkm.profiles.tsv')
        self.assertEqual(header, ['Bin Id', 'Completeness', 'Taxonomy'])
        self.assertEqual([row[0] for row in rows], proteins)
        self.assertEqual(self.first_column('checkm.qa_sh100.tsv'), proteins)
        with open(os.path.join(self.out, M.CHECKM_RELEASE_ALIGNMENT)) as handle:
            self.assertEqual([line.split('\t')[0] for line in handle], proteins)

    def test_a_genome_without_proteins_is_named_rather_than_failing_its_batch(self):
        dirs, report = self.release({'GCA_000000001.1': STATUS_NEW,
                                     'GCA_000000002.1': STATUS_NEW},
                                    missing={'GCA_000000001.1'})
        self.run_checkm(dirs, report)

        self.assertEqual(self.first_column('checkm.profiles.tsv'), ['GCA_000000002.1_protein'])
        self.assertEqual(self.not_assessed(M.CHECKM_NOT_ASSESSED),
                         {'GCA_000000001.1': M.REASON_NO_PROTEINS})

    def test_a_batch_with_nothing_to_assess_adds_nothing_to_the_release(self):
        # the second batch holds only a genome with no proteins
        dirs, report = self.release({'GCA_000000001.1': STATUS_NEW,
                                     'GCA_000000002.1': STATUS_NEW,
                                     'GCA_000000003.1': STATUS_NEW},
                                    missing={'GCA_000000003.1'})
        self.assertTrue(self.run_checkm(dirs, report, batch_size=2))

        self.assertEqual(len(self.calls_of('lineage_wf')), 1)
        self.assertEqual(self.first_column('checkm.profiles.tsv'),
                         ['GCA_000000001.1_protein', 'GCA_000000002.1_protein'])
        self.assertEqual(self.not_assessed(M.CHECKM_NOT_ASSESSED),
                         {'GCA_000000003.1': M.REASON_NO_PROTEINS})

    def test_the_release_names_no_genome_left_out_with_a_header_and_no_rows(self):
        dirs, report = self.release({'GCA_000000001.1': STATUS_NEW})
        self.run_checkm(dirs, report)

        self.assertEqual(self.not_assessed(M.CHECKM_NOT_ASSESSED), {})

    def test_the_batches_are_planned_in_the_out_dir_itself(self):
        dirs, report = self.release({'GCA_000000001.1': STATUS_NEW})
        self.run_checkm(dirs, report)

        self.assertTrue(os.path.exists(os.path.join(self.batch(1), SUCCESS_CANARY)))
        self.assertTrue(os.path.exists(os.path.join(self.out, 'checkm.profiles.tsv')))

    def test_an_out_dir_holding_another_commands_batches_is_refused(self):
        dirs, report = self.release({'GCA_000000001.1': STATUS_NEW})
        self.run_checkm(dirs, report)

        # checkm's SUCCESS would otherwise tell checkm2 its batch was done
        with self.assertRaises(RuntimeError):
            self.command(M.CheckM2, dirs, report)
        self.assertEqual(self.calls_of('predict'), [])


if __name__ == '__main__':
    unittest.main()
