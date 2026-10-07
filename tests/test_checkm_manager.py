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
import json
import subprocess
import tempfile
import time
import unittest
from unittest import mock

from gtdb_migration_tk import checkm_manager as M
from gtdb_migration_tk.batching import (FAILED_CANARY, SUCCESS_CANARY,
                                        batchfile_path, read_batchfile)
from gtdb_migration_tk.ncbi_utils import assembly_stats
from gtdb_migration_tk.update_genomes import (STATUS_FASTA_CHANGED,
                                              STATUS_FASTA_UNCHANGED, STATUS_NEW,
                                              STATUS_REMOVED,
                                              STATUS_SEQUENCES_UNCHANGED,
                                              STATUS_TO_CURATE)
from gtdb_migration_tk.utils.common import open_text

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
        # CheckM runs pplacer by bare name from its own PATH, here PPLACER
        [ -n "$PPLACER" ] && "$PPLACER" 0.2
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

        # pplacer as bioconda installs it: bin/pplacer a link to bin/pplacer.exe,
        # and the package named in conda-meta. It is a copy of sleep, since the
        # executable of a running shell script is the shell
        self.pplacer_env = os.path.join(self.dir, 'envs', 'pplacer-1.1.alpha22')
        os.makedirs(os.path.join(self.pplacer_env, 'bin'))
        os.makedirs(os.path.join(self.pplacer_env, 'conda-meta'))
        self.pplacer_exe = os.path.join(self.pplacer_env, 'bin', 'pplacer.exe')
        shutil.copy(shutil.which('sleep'), self.pplacer_exe)
        os.symlink('pplacer.exe', os.path.join(self.pplacer_env, 'bin', 'pplacer'))
        with open(os.path.join(self.pplacer_env, 'conda-meta',
                               'pplacer-1.1.alpha22-hd563303_0.json'), 'w') as handle:
            json.dump({'name': 'pplacer', 'version': '1.1.alpha22',
                       'build': 'hd563303_0'}, handle)

        # pplacer runs for 0.2 seconds, looked for every 20 ms
        watch = mock.patch.object(M, 'WATCH_SECONDS', 0.02)
        watch.start()
        self.addCleanup(watch.stop)

        self.calls = os.path.join(self.dir, 'calls.txt')
        for name, value in (('PATH', bin_dir + os.pathsep + os.environ['PATH']),
                            ('CALLS', self.calls),
                            ('FAIL_BATCH', ''),
                            ('PPLACER', os.path.join(self.pplacer_env, 'bin', 'pplacer'))):
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

    def command(self, program, dirs, report, cpus=1, batch_size=10, all_genomes=False,
                **kwargs):
        manager = program(cpus=cpus, tmp_dir=self.tmp, batch_size=batch_size,
                          heartbeat=3600, **kwargs)
        return manager.run(dirs, report, self.out, all_genomes)

    def table(self, name):
        # a release table is gzipped, a batch's is not and not_assessed is not
        with open_text(os.path.join(self.out, name)) as handle:
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

    def genome_dir(self, gid):
        return os.path.join(self.dir, 'genomes', gid + '_ASM1v1')


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
        self.assertEqual(self.first_column(M.CHECKM_RELEASE_PROFILE),
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

    def test_the_release_report_is_gzipped_and_an_uncompressed_one_left_before_is_removed(self):
        os.makedirs(self.out, exist_ok=True)
        stale = os.path.join(self.out, M.CHECKM2_SUPERSEDED_FILES[0])
        with open(stale, 'w') as handle:
            handle.write('stale\n')
        dirs, report = self.release({'GCA_000000001.1': STATUS_NEW})
        self.assertTrue(self.run_checkm2(dirs, report))

        with open(os.path.join(self.out, M.CHECKM2_RELEASE_REPORT), 'rb') as handle:
            self.assertEqual(handle.read(2), b'\x1f\x8b')
        self.assertFalse(os.path.exists(stale))
        # the batch keeps its own report uncompressed, as CheckM2 wrote it
        with open(os.path.join(self.batch(1), M.CHECKM2_BATCH_REPORT)) as handle:
            self.assertEqual(handle.readline().split('\t')[0], 'Name')

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

    def test_pplacer_is_given_the_cpus_up_to_64_and_the_rest_of_checkm_all_of_them(self):
        dirs, report = self.release({'GCA_000000001.1': STATUS_NEW})
        for cpus, pplacer in ((8, '8'), (64, '64'), (96, '64')):
            # each a run of its own, from nothing
            shutil.rmtree(self.out, True)
            if os.path.exists(self.calls):
                os.remove(self.calls)
            self.assertTrue(self.run_checkm(dirs, report, cpus=cpus))

            [args] = self.calls_of('lineage_wf')
            self.assertEqual(args[args.index('--pplacer_threads') + 1], pplacer)
            self.assertEqual(args[args.index('-t') + 1], str(cpus))
            for qa in self.calls_of('qa'):
                self.assertEqual(qa[qa.index('-t') + 1], str(cpus))

    def test_the_profile_joins_qa_and_tree_qa_for_each_genome(self):
        dirs, report = self.release({'GCA_000000001.1': STATUS_NEW,
                                     'GCA_000000002.1': STATUS_FASTA_UNCHANGED,
                                     'GCA_000000003.1': STATUS_NEW})
        self.run_checkm(dirs, report)

        header, rows = self.table(M.CHECKM_RELEASE_PROFILE)
        self.assertEqual(header, ['Bin Id', 'Completeness', 'Taxonomy'])
        self.assertEqual(rows, [['GCA_000000001.1_protein', '99.0', 'k__Bacteria'],
                                ['GCA_000000003.1_protein', '99.0', 'k__Bacteria']])
        self.assertEqual(self.first_column(M.CHECKM_RELEASE_QA_SH100),
                         ['GCA_000000001.1_protein', 'GCA_000000003.1_protein'])

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
        self.assertFalse(os.path.exists(os.path.join(self.out, M.CHECKM_RELEASE_PROFILE)))

    def test_the_release_tables_are_the_batches_concatenated_under_one_header(self):
        genomes = {'GCA_%09d.1' % i: STATUS_NEW for i in range(5)}
        dirs, report = self.release(genomes)
        self.run_checkm(dirs, report, batch_size=2)

        proteins = [gid + '_protein' for gid in sorted(genomes)]
        header, rows = self.table(M.CHECKM_RELEASE_PROFILE)
        self.assertEqual(header, ['Bin Id', 'Completeness', 'Taxonomy'])
        self.assertEqual([row[0] for row in rows], proteins)
        self.assertEqual(self.first_column(M.CHECKM_RELEASE_QA_SH100), proteins)

    def test_the_release_tables_are_gzipped_and_the_genomes_left_out_are_not(self):
        dirs, report = self.release({'GCA_000000001.1': STATUS_NEW})
        self.run_checkm(dirs, report)

        for name in (M.CHECKM_RELEASE_PROFILE, M.CHECKM_RELEASE_QA_SH100):
            with open(os.path.join(self.out, name), 'rb') as handle:
                self.assertEqual(handle.read(2), b'\x1f\x8b', name)
        with open(os.path.join(self.out, M.CHECKM_NOT_ASSESSED)) as handle:
            self.assertEqual(handle.readline().rstrip('\n').split('\t'), list(M.NOT_ASSESSED_HEADER))

    def test_the_alignments_stay_in_their_batches_and_are_not_gathered_for_the_release(self):
        # nothing reads them, and for r237 the release's copy was 142 GB
        genomes = {'GCA_%09d.1' % i: STATUS_NEW for i in range(3)}
        dirs, report = self.release(genomes)
        self.run_checkm(dirs, report, batch_size=2)

        with open(os.path.join(self.batch(1), M.CHECKM_ALIGNMENT)) as handle:
            self.assertEqual([line.split('\t')[0] for line in handle],
                             [gid + '_protein' for gid in sorted(genomes)[:2]])
        self.assertTrue(os.path.exists(os.path.join(self.batch(2), M.CHECKM_ALIGNMENT)))
        self.assertEqual([name for name in os.listdir(self.out) if 'alignment' in name], [])

    def test_what_an_earlier_run_wrote_for_the_release_is_removed_once_replaced(self):
        # an uncompressed table beside the gzipped one, a few genomes short, is
        # how the stale one gets read; the joined alignments are no longer made
        os.makedirs(self.out, exist_ok=True)
        for name in M.CHECKM_SUPERSEDED_FILES:
            with open(os.path.join(self.out, name), 'w') as handle:
                handle.write('stale\n')
        dirs, report = self.release({'GCA_000000001.1': STATUS_NEW})
        self.run_checkm(dirs, report)

        for name in M.CHECKM_SUPERSEDED_FILES:
            self.assertFalse(os.path.exists(os.path.join(self.out, name)), name)
        self.assertEqual(self.first_column(M.CHECKM_RELEASE_PROFILE), ['GCA_000000001.1_protein'])

    def test_an_earlier_release_file_is_kept_while_a_batch_is_unfinished(self):
        os.makedirs(self.out, exist_ok=True)
        stale = os.path.join(self.out, M.CHECKM_SUPERSEDED_FILES[0])
        with open(stale, 'w') as handle:
            handle.write('stale\n')
        dirs, report = self.release({'GCA_000000001.1': STATUS_NEW})
        os.environ['FAIL_BATCH'] = 'batch_000001'
        self.assertFalse(self.run_checkm(dirs, report))

        self.assertTrue(os.path.exists(stale))

    def test_a_genome_without_proteins_is_named_rather_than_failing_its_batch(self):
        dirs, report = self.release({'GCA_000000001.1': STATUS_NEW,
                                     'GCA_000000002.1': STATUS_NEW},
                                    missing={'GCA_000000001.1'})
        self.run_checkm(dirs, report)

        self.assertEqual(self.first_column(M.CHECKM_RELEASE_PROFILE), ['GCA_000000002.1_protein'])
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
        self.assertEqual(self.first_column(M.CHECKM_RELEASE_PROFILE),
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
        self.assertTrue(os.path.exists(os.path.join(self.out, M.CHECKM_RELEASE_PROFILE)))

    def test_an_out_dir_holding_another_commands_batches_is_refused(self):
        dirs, report = self.release({'GCA_000000001.1': STATUS_NEW})
        self.run_checkm(dirs, report)

        # checkm's SUCCESS would otherwise tell checkm2 its batch was done
        with self.assertRaises(RuntimeError):
            self.command(M.CheckM2, dirs, report)
        self.assertEqual(self.calls_of('predict'), [])


# ------------------------------------------------------------ which version

class WhichVersionMadeTheRelease(TempDirCase):
    """The version beside the release files is the batches', not the last machine's."""

    GENOMES = {'GCA_000000001.1': STATUS_NEW, 'GCA_000000002.1': STATUS_NEW}

    def release_version(self, program):
        with open(os.path.join(self.out, program + '.version')) as handle:
            return handle.read().splitlines()

    def say_batch_was_made_by(self, number, program, version):
        with open(os.path.join(self.batch(number), program + '.version'), 'w') as handle:
            handle.write(version + '\n')

    def test_checkm_writes_its_version_beside_the_release_files(self):
        dirs, report = self.release(self.GENOMES)
        self.command(M.CheckM, dirs, report)

        self.assertEqual(self.release_version('checkm'), ['CheckM v1.2.5'])

    def test_checkm2_writes_its_version_beside_the_release_files(self):
        dirs, report = self.release(self.GENOMES)
        self.command(M.CheckM2, dirs, report)

        self.assertEqual(self.release_version('checkm2'), ['1.1.0'])

    def test_the_version_is_logged_where_the_run_starts(self):
        with self.assertLogs('timestamp', level='INFO') as logs:
            M.CheckM(tmp_dir=self.tmp)
            M.CheckM2(tmp_dir=self.tmp)

        self.assertIn('Using checkm: CheckM v1.2.5.', ' '.join(logs.output))
        self.assertIn('Using checkm2: 1.1.0.', ' '.join(logs.output))

    def test_each_batch_log_says_which_version_assessed_it(self):
        dirs, report = self.release(self.GENOMES)
        # at INFO, as a run's log is
        with self.assertLogs('timestamp', level='INFO'):
            self.command(M.CheckM, dirs, report, batch_size=1)

        for number in (1, 2):
            with open(os.path.join(self.batch(number), M.CHECKM_LAYOUT.log)) as handle:
                self.assertIn('starting with checkm (CheckM v1.2.5).', handle.read())

    def test_the_checkm_release_log_names_checkm_and_pplacer(self):
        dirs, report = self.release(self.GENOMES)
        with self.assertLogs('timestamp', level='INFO') as logs:
            self.command(M.CheckM, dirs, report)

        self.assertIn('Release: 2 genome(s) assessed with checkm (CheckM v1.2.5), '
                      'pplacer (1.1.alpha22).', ' '.join(logs.output))

    def test_the_release_log_names_the_version(self):
        dirs, report = self.release(self.GENOMES)
        with self.assertLogs('timestamp', level='INFO') as logs:
            self.command(M.CheckM2, dirs, report)

        self.assertIn('Release: 2 genome(s) assessed with checkm2 (1.1.0).',
                      ' '.join(logs.output))

    def test_it_is_the_version_the_batches_were_made_with_not_this_machines(self):
        dirs, report = self.release(self.GENOMES)
        self.command(M.CheckM2, dirs, report, batch_size=1)
        for number in (1, 2):
            self.say_batch_was_made_by(number, 'checkm2', '1.0.2')

        self.command(M.CheckM2, dirs, report, batch_size=1)

        self.assertEqual(self.release_version('checkm2'), ['1.0.2'])

    def test_batches_made_by_two_versions_name_both_and_warn_which_made_which(self):
        dirs, report = self.release(self.GENOMES)
        self.command(M.CheckM, dirs, report, batch_size=1)
        self.say_batch_was_made_by(1, 'checkm', 'CheckM v1.2.4')

        with self.assertLogs('timestamp', level='WARNING') as logs:
            self.command(M.CheckM, dirs, report, batch_size=1)

        self.assertEqual(self.release_version('checkm'), ['CheckM v1.2.4', 'CheckM v1.2.5'])
        warning = ' '.join(logs.output)
        self.assertIn('CheckM v1.2.4 by 1 batch(es): batch_000001', warning)
        self.assertIn('CheckM v1.2.5 by 1 batch(es): batch_000002', warning)

    def test_a_release_no_batch_ran_the_program_over_has_no_version(self):
        dirs, report = self.release(self.GENOMES, missing=set(self.GENOMES))
        self.command(M.CheckM, dirs, report)

        self.assertTrue(os.path.exists(os.path.join(self.out, M.CHECKM_NOT_ASSESSED)))
        self.assertFalse(os.path.exists(os.path.join(self.out, 'checkm.version')))


# ------------------------------------------------------------ pplacer

class RecordingPplacer(TempDirCase):
    """The pplacer CheckM ran, which no --version can name, is learned from the run."""

    GENOMES = {'GCA_000000001.1': STATUS_NEW, 'GCA_000000002.1': STATUS_NEW}

    def pplacer_version(self, directory):
        with open(os.path.join(directory, 'pplacer.version')) as handle:
            return handle.read().splitlines()

    def test_the_pplacer_checkm_ran_is_recorded_beside_the_batch(self):
        dirs, report = self.release(self.GENOMES)
        self.command(M.CheckM, dirs, report, batch_size=1)

        for number in (1, 2):
            self.assertEqual(self.pplacer_version(self.batch(number)), ['1.1.alpha22'])

    def test_the_batch_log_names_the_pplacer_and_where_it_ran_from(self):
        dirs, report = self.release(self.GENOMES)
        with self.assertLogs('timestamp', level='INFO'):
            self.command(M.CheckM, dirs, report)

        with open(os.path.join(self.batch(1), M.CHECKM_LAYOUT.log)) as handle:
            self.assertIn('Using pplacer: 1.1.alpha22, run from {}.'.format(
                os.path.realpath(self.pplacer_exe)), handle.read())

    def test_it_is_the_pplacer_checkm_ran_not_one_on_the_toolkits_path(self):
        # CheckM's PATH is set by a wrapper the toolkit never sees; a pplacer the
        # toolkit could find is not the one asked about
        self.assertIsNone(shutil.which('pplacer'))
        dirs, report = self.release(self.GENOMES)
        self.command(M.CheckM, dirs, report)

        self.assertEqual(self.pplacer_version(self.batch(1)), ['1.1.alpha22'])

    def test_the_release_pplacer_version_is_gathered_from_the_batches(self):
        dirs, report = self.release(self.GENOMES)
        self.command(M.CheckM, dirs, report, batch_size=1)

        self.assertEqual(self.pplacer_version(self.out), ['1.1.alpha22'])

    def test_a_pplacer_outside_conda_is_named_and_given_no_version(self):
        shutil.rmtree(os.path.join(self.pplacer_env, 'conda-meta'))
        dirs, report = self.release(self.GENOMES)
        with self.assertLogs('timestamp', level='WARNING') as logs:
            self.assertTrue(self.command(M.CheckM, dirs, report))

        self.assertFalse(os.path.exists(os.path.join(self.batch(1), 'pplacer.version')))
        self.assertIn('pplacer ran from {}'.format(os.path.realpath(self.pplacer_exe)),
                      ' '.join(logs.output))

    def test_a_pplacer_not_seen_to_run_is_said_and_given_no_version(self):
        os.environ['PPLACER'] = ''
        dirs, report = self.release(self.GENOMES)
        with self.assertLogs('timestamp', level='WARNING') as logs:
            self.assertTrue(self.command(M.CheckM, dirs, report))

        self.assertFalse(os.path.exists(os.path.join(self.batch(1), 'pplacer.version')))
        self.assertIn('pplacer was not seen to run', ' '.join(logs.output))

    def test_a_batch_that_recorded_no_pplacer_leaves_the_release_without_one(self):
        dirs, report = self.release(self.GENOMES)
        self.command(M.CheckM, dirs, report, batch_size=1)
        self.assertTrue(os.path.exists(os.path.join(self.out, 'pplacer.version')))
        # as every r237 batch made before pplacer was recorded
        os.remove(os.path.join(self.batch(1), 'pplacer.version'))

        with self.assertLogs('timestamp', level='WARNING') as logs:
            self.command(M.CheckM, dirs, report, batch_size=1)

        self.assertFalse(os.path.exists(os.path.join(self.out, 'pplacer.version')))
        self.assertIn('1 batch(es) that ran checkm recorded no pplacer version, so none '
                      'is written for the release: batch_000001.', ' '.join(logs.output))
        with open(os.path.join(self.out, 'checkm.version')) as handle:
            self.assertEqual(handle.read().strip(), 'CheckM v1.2.5')

    def test_a_batch_made_again_does_not_keep_the_pplacer_of_the_last_attempt(self):
        dirs, report = self.release(self.GENOMES)
        self.command(M.CheckM, dirs, report)
        os.remove(os.path.join(self.batch(1), SUCCESS_CANARY))
        os.environ['PPLACER'] = ''

        self.command(M.CheckM, dirs, report)

        self.assertFalse(os.path.exists(os.path.join(self.batch(1), 'pplacer.version')))

    def test_checkm2_records_no_pplacer(self):
        dirs, report = self.release(self.GENOMES)
        self.command(M.CheckM2, dirs, report)

        self.assertFalse(os.path.exists(os.path.join(self.batch(1), 'pplacer.version')))
        self.assertFalse(os.path.exists(os.path.join(self.out, 'pplacer.version')))


class WatchingAProcess(TempDirCase):
    def test_the_executables_of_every_descendant_are_found(self):
        proc = subprocess.Popen(['sh', '-c', '"$0" 5; true', self.pplacer_exe])
        self.addCleanup(proc.wait)
        self.addCleanup(proc.kill)
        found = set()
        for _ in range(200):
            found = set(M.descendant_executables(proc.pid).values())
            if os.path.realpath(self.pplacer_exe) in found:
                break
            time.sleep(0.01)

        self.assertIn(os.path.realpath(self.pplacer_exe), found)

    def test_a_process_with_no_descendants_has_none(self):
        proc = subprocess.Popen([self.pplacer_exe, '5'])
        self.addCleanup(proc.wait)
        self.addCleanup(proc.kill)

        self.assertEqual(M.descendant_executables(proc.pid), {})

    def test_pplacer_is_known_by_the_name_bioconda_runs_it_under(self):
        self.assertTrue(M.is_program('/env/bin/pplacer.exe', 'pplacer'))
        self.assertTrue(M.is_program('/env/bin/pplacer', 'pplacer'))
        self.assertFalse(M.is_program('/env/bin/guppy', 'pplacer'))


# ------------------------------------------------------------ a genome too large

# GCA_964261755.1, a faecal metagenome deposited as one genome
METAGENOME_BASES = 9528631298


class AGenomeTooLargeToAssess(TempDirCase):
    """Named and left, rather than failing its batch on every machine that takes it."""

    SMALL, LARGE = 'GCA_000000001.1', 'GCA_000000002.1'

    def release_with_a_metagenome(self):
        dirs, report = self.release({self.SMALL: STATUS_NEW, self.LARGE: STATUS_NEW})
        with open(assembly_stats(self.genome_dir(self.LARGE)), 'w') as handle:
            handle.write('all\tall\tall\tall\ttotal-length\t{}\n'.format(METAGENOME_BASES))
        return dirs, report

    def test_checkm_is_not_handed_it_and_the_release_names_it(self):
        dirs, report = self.release_with_a_metagenome()
        self.assertTrue(self.command(M.CheckM, dirs, report))

        self.assertEqual(self.first_column(M.CHECKM_RELEASE_PROFILE), [self.SMALL + '_protein'])
        self.assertEqual(self.not_assessed(M.CHECKM_NOT_ASSESSED),
                         {self.LARGE: M.REASON_GENOME_TOO_LARGE})

    def test_checkm2_is_not_handed_it_and_the_release_names_it(self):
        dirs, report = self.release_with_a_metagenome()
        self.assertTrue(self.command(M.CheckM2, dirs, report))

        self.assertEqual(self.first_column(M.CHECKM2_RELEASE_REPORT), [self.SMALL])
        self.assertEqual(self.not_assessed(M.CHECKM2_NOT_ASSESSED),
                         {self.LARGE: M.REASON_GENOME_TOO_LARGE})

    def test_its_size_in_bases_is_the_detail(self):
        dirs, report = self.release_with_a_metagenome()
        self.command(M.CheckM, dirs, report)

        _, rows = self.table(M.CHECKM_NOT_ASSESSED)
        self.assertEqual(rows, [[self.LARGE, M.REASON_GENOME_TOO_LARGE, str(METAGENOME_BASES)]])

    def test_the_limit_is_the_one_given(self):
        dirs, report = self.release_with_a_metagenome()
        self.command(M.CheckM, dirs, report, max_genome_size=10000)

        self.assertEqual(self.first_column(M.CHECKM_RELEASE_PROFILE),
                         [self.SMALL + '_protein', self.LARGE + '_protein'])
        self.assertEqual(self.not_assessed(M.CHECKM_NOT_ASSESSED), {})

    def test_it_is_planned_so_a_limit_applies_to_batches_already_planned(self):
        """r237's batches were planned before there was a limit."""
        dirs, report = self.release_with_a_metagenome()
        self.command(M.CheckM, dirs, report)

        self.assertEqual(self.planned(M.CHECKM_LAYOUT), [self.SMALL, self.LARGE])

    def test_a_batch_of_nothing_but_it_runs_no_checkm_and_succeeds(self):
        dirs, report = self.release_with_a_metagenome()
        self.assertTrue(self.command(M.CheckM, dirs, report, batch_size=1))

        self.assertEqual(len(self.calls_of('lineage_wf')), 1)
        self.assertTrue(os.path.exists(os.path.join(self.batch(2), SUCCESS_CANARY)))

    def test_the_default_limit_is_the_one_select_genomes_has(self):
        for program in (M.CheckM, M.CheckM2):
            manager = program(tmp_dir=self.tmp)
            self.assertEqual(manager.max_genome_bases, 100 * 1000 * 1000)


if __name__ == '__main__':
    unittest.main()
