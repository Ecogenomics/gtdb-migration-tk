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

"""Offline unit tests for marker_alignment_manager.py -- align_marker_genes, which
was gtdb's `power realign_updated_genomes`.

That ran hmmalign once for each marker of each genome, wrote each genome on its
own connection, passed over a worker that died, and failed with a TypeError where
hmmalign gave no alignment. hmmalign is stood in for by a script that aligns each
sequence to the HMM's length with one insert column after it, and the database
by a fake cursor.
"""

import gzip
import logging
import os
import shutil
import stat
import sys
import tempfile
import unittest
from unittest import mock

from gtdb_migration_tk import __main__ as main_module
from gtdb_migration_tk import main as main_py
from gtdb_migration_tk import marker_alignment_manager as A

# a stand-in hmmalign: the HMM's LENG match states, then one insert column the
# alignment must not keep
FAKE_HMMALIGN = r'''#!{python}
import sys
hmm, fasta = sys.argv[-2], sys.argv[-1]
length = int([l.split()[1] for l in open(hmm) if l.startswith('LENG')][0])
names, seqs = [], {{}}
for line in open(fasta):
    if line.startswith('>'):
        names.append(line[1:].strip()); seqs[names[-1]] = ''
    else:
        seqs[names[-1]] += line.strip()
print('# STOCKHOLM 1.0')
for name in names:
    print('{{}} {{}}'.format(name, seqs[name][:length].ljust(length, '-') + 'q'))
print('#=GC RF ' + 'x' * length + '.')
print('//')
'''

PF = A.Marker(db_id=101, accession='PF00410.20', hmm='', size=6, database='PFAM')
TIGR = A.Marker(db_id=202, accession='TIGR00006', hmm='', size=4, database='TIGR')
TOPHIT_HEADER = 'Gene Id\tTop hits (Family id,e-value,bitscore)\n'


class FakeCursor(object):
    """A database of marker sets, markers and genomes, recording each statement."""

    def __init__(self, markers=(), genomes=(), sets=((1, 'bac120'), (2, 'ar122'))):
        self.sets = list(sets)
        self.markers = list(markers)
        self.genomes = list(genomes)
        self.statements = []
        self.copied = []
        self.result = []
        self.rowcount = 0

    def execute(self, sql, params=None):
        self.statements.append(sql)
        if sql.startswith('SELECT id, name FROM marker_sets'):
            self.result = [s for s in self.sets if s[0] in params[0]]
        elif 'FROM markers m' in sql:
            self.result = [(m.db_id, m.accession, m.hmm, m.size, m.database) for m in self.markers]
        elif 'FROM genomes g' in sql:
            self.result = list(self.genomes)
        elif sql.startswith('INSERT INTO aligned_markers'):
            self.rowcount = len(self.copied[-1])

    def copy_expert(self, sql, handle):
        self.statements.append(sql)
        self.copied.append([line.split('\t') for line in handle.read().splitlines()])

    def fetchall(self):
        return list(self.result)


class FakeConnection(object):
    def __init__(self):
        self.commits = 0
        self.rollbacks = 0

    def commit(self):
        self.commits += 1

    def rollback(self):
        self.rollbacks += 1


class TempDirCase(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='marker_alignment_test.')
        self.addCleanup(shutil.rmtree, self.dir, True)
        logging.getLogger('timestamp').addHandler(logging.NullHandler())

        bin_dir = os.path.join(self.dir, 'bin')
        os.makedirs(bin_dir)
        self.hmmalign = os.path.join(bin_dir, 'hmmalign')
        with open(self.hmmalign, 'w') as handle:
            handle.write(FAKE_HMMALIGN.format(python=sys.executable))
        os.chmod(self.hmmalign, os.stat(self.hmmalign).st_mode | stat.S_IEXEC)
        for name, value in (('HMMALIGN', self.hmmalign),):
            patcher = mock.patch.object(A, name, value)
            patcher.start()
            self.addCleanup(patcher.stop)

        self.pf = PF._replace(hmm=self.hmm('pfam/33.1/PF00410.20.hmm', PF.size))
        self.tigr = TIGR._replace(hmm=self.hmm('tigrfam/15.0/TIGR00006.HMM', TIGR.size))

    def hmm(self, name, length):
        path = os.path.join(self.dir, 'hmms', name)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, 'w') as handle:
            handle.write('HMMER3/f [3.3 | Nov 2019]\nNAME  x\nLENG  {}\n//\n'.format(length))
        return path

    def genome(self, accession, proteins=None, pfam=None, tigrfam=None):
        """A genome directory: proteins {gene: seq}, and top-hit tables of rows (gene, hits)."""
        gdir = os.path.join(self.dir, 'genomes', accession + '_ASM1v1')
        prodigal = os.path.join(gdir, 'prodigal')
        os.makedirs(prodigal)
        if proteins is not None:
            with gzip.open(os.path.join(prodigal, accession + '_protein.faa.gz'), 'wt') as handle:
                for gene, seq in proteins.items():
                    handle.write('>{} # 1 # 2\n{}*\n'.format(gene, seq))
        for database, rows in (('PFAM', pfam), ('TIGR', tigrfam)):
            if rows is None:
                continue
            path = A.tophit_file(prodigal, accession, database)
            os.makedirs(os.path.dirname(path), exist_ok=True)
            with gzip.open(path, 'wt') as handle:
                handle.write(TOPHIT_HEADER)
                for gene, hits in rows:
                    handle.write('{}\t{}\n'.format(gene, hits))
        return gdir

    def manager(self, cursor):
        manager = A.MarkerAlignmentManager.__new__(A.MarkerAlignmentManager)
        manager.logger = logging.getLogger('timestamp')
        manager.cpus, manager.batch_size, manager.tmp_dir = 2, 2, os.path.join(self.dir, 'tmp')
        manager.reclaim, manager.lease, manager.heartbeat = False, 7200, 300
        manager.hmmalign_version = ''
        manager.temp_cur, manager.temp_con = cursor, FakeConnection()
        return manager


class ChoosingTheGenes(TempDirCase):
    def test_the_highest_bitscore_wins_the_first_of_equal_ones_is_kept_and_the_hits_are_counted(self):
        hits = {'PF00410.20': [('g1', '1e-5', 20.0), ('g2', '1e-9', 30.0), ('g3', '1e-9', 30.0), ('g4', '1', 5.0)]}
        chosen = A.choose_genes(hits, {'g1': 'MKA', 'g2': 'MKB', 'g3': 'MKB', 'g4': 'MKC'})
        self.assertEqual(chosen['PF00410.20'], A.ChosenGene(sequence='MKB', evalue='1e-9', bitscore='30.0',
                                                           multiple_hits=True, hit_number=4, unique_genes=3))

    def test_one_hit_is_not_a_multiple_hit(self):
        chosen = A.choose_genes({'TIGR00006': [('g1', '2e-10', 41.5)]}, {'g1': 'MKA'})
        self.assertEqual(chosen['TIGR00006'], A.ChosenGene('MKA', '2e-10', '41.5', False, 1, 1))

    def test_a_tophit_table_is_read_for_the_wanted_markers_only(self):
        gdir = self.genome('GCA_000000001.1', pfam=[('g1', 'PF00410.20,1e-5,20.0;PF99999.1,1e-3,9.0'),
                                                    ('g2', 'PF00410.20,1e-9,30.0')])
        path = A.tophit_file(os.path.join(gdir, 'prodigal'), 'GCA_000000001.1', 'PFAM')
        self.assertTrue(path.endswith('prodigal/pfam_33.1_lite/GCA_000000001.1_pfam_33.1_lite_tophit.tsv.gz'))
        self.assertEqual(A.read_tophits(path, {'PF00410.20'}),
                         {'PF00410.20': [('g1', '1e-5', 20.0), ('g2', '1e-9', 30.0)]})


class ReadingAGenome(TempDirCase):
    JOB_MARKERS = (('PFAM', ('PF00410.20',)), ('TIGR', ('TIGR00006',)))

    def read(self, accession, gdir):
        return A.read_genome((accession, os.path.join(gdir, 'prodigal', accession + '_protein.faa.gz'),
                              self.JOB_MARKERS))

    def test_prodigals_asterisk_is_taken_off_and_each_database_read(self):
        gdir = self.genome('GCA_000000001.1', {'g1': 'MKAL', 'g2': 'MPQ'},
                           pfam=[('g1', 'PF00410.20,1e-5,20.0')], tigrfam=[('g2', 'TIGR00006,1e-9,30.0')])
        accession, genes = self.read('GCA_000000001.1', gdir)
        self.assertEqual((genes['PF00410.20'].sequence, genes['TIGR00006'].sequence), ('MKAL', 'MPQ'))

    def test_a_genome_without_a_protein_file_is_none(self):
        gdir = self.genome('GCA_000000001.1', pfam=[], tigrfam=[])
        self.assertEqual(self.read('GCA_000000001.1', gdir), ('GCA_000000001.1', None))

    def test_proteins_without_a_tophit_table_is_an_error_naming_the_genome(self):
        gdir = self.genome('GCA_000000001.1', {'g1': 'MKAL'}, pfam=[])
        with self.assertRaisesRegex(A.AlignmentError, 'GCA_000000001.1 has called proteins and no top-hit table'):
            self.read('GCA_000000001.1', gdir)

    def test_a_tophit_naming_a_gene_the_proteins_do_not_hold_is_an_error(self):
        gdir = self.genome('GCA_000000001.1', {'g1': 'MKAL'}, pfam=[('g9', 'PF00410.20,1e-5,20.0')], tigrfam=[])
        with self.assertRaisesRegex(A.AlignmentError, 'names gene g9'):
            self.read('GCA_000000001.1', gdir)


class Aligning(TempDirCase):
    def test_only_the_match_states_are_kept(self):
        lines = ['# STOCKHOLM 1.0', 'a MK-.L', 'b MKAqL', '#=GC RF xxx.x', '//']
        self.assertEqual(A.aligned_match_states(lines, {'a', 'b'}), {'a': 'MK-L', 'b': 'MKAL'})

    def test_a_sequence_hmmalign_left_out_or_no_reference_line_is_an_error(self):
        # gtdb's _runHmmAlign then failed on len(None)
        with self.assertRaisesRegex(A.AlignmentError, 'no alignment of 1 sequence'):
            A.aligned_match_states(['a MK', '#=GC RF xx'], {'a', 'b'})
        with self.assertRaisesRegex(A.AlignmentError, 'no #=GC RF line'):
            A.aligned_match_states(['a MK'], {'a'})

    def test_the_genes_of_a_batch_are_aligned_to_a_marker_by_one_hmmalign(self):
        work = os.path.join(self.dir, 'work')
        os.makedirs(work)
        with mock.patch.object(A.subprocess, 'run', wraps=A.subprocess.run) as run:
            aligned = A.align_marker(self.pf, {'0': 'MKALIE', '1': 'MKA'}, work)
        self.assertEqual(run.call_count, 1)
        self.assertEqual(aligned, {'0': 'MKALIE', '1': 'MKA---'})
        self.assertEqual(os.listdir(work), [])

    def test_a_failed_hmmalign_is_an_error_with_its_message(self):
        with open(self.hmmalign, 'w') as handle:
            handle.write('#!/bin/sh\necho "Error: bad HMM" >&2\nexit 1\n')
        with self.assertRaisesRegex(A.AlignmentError, 'exit status 1: Error: bad HMM'):
            A.align_marker(self.pf, {'0': 'MK'}, self.dir)

    def test_an_alignment_of_another_length_than_the_marker_is_an_error(self):
        with self.assertRaisesRegex(A.AlignmentError, 'other than its 9 match states'):
            A.align_marker(self.pf._replace(size=9), {'0': 'MK'}, self.dir)


class TheMarkers(TempDirCase):
    def test_a_set_the_database_does_not_hold_is_refused(self):
        with self.assertRaisesRegex(A.AlignmentError, 'no marker set 7'):
            self.manager(FakeCursor([self.pf])).read_markers([1, 7])

    def test_an_hmm_of_another_length_or_version_or_missing_is_refused_naming_each(self):
        markers = [self.pf._replace(size=99), self.tigr._replace(hmm=self.tigr.hmm.replace('15.0', '14.0')),
                   self.pf._replace(accession='PF00001.1', hmm=os.path.join(self.dir, 'hmms/pfam/33.1/none.hmm'))]
        with self.assertRaises(A.AlignmentError) as raised:
            self.manager(FakeCursor(markers)).read_markers([1])
        message = str(raised.exception)
        self.assertIn('3 marker(s) cannot be aligned to', message)
        self.assertIn('PF00410.20 has 6 match states', message)
        self.assertIn('not of tigrfam 15.0', message)
        self.assertIn('PF00001.1', message)

    def test_the_markers_of_the_sets_are_read(self):
        markers = self.manager(FakeCursor([self.tigr, self.pf])).read_markers([1, 2])
        self.assertEqual([m.accession for m in markers], ['PF00410.20', 'TIGR00006'])


class AligningARelease(TempDirCase):
    def release(self):
        """Three genomes: one with both markers, one with neither hit, one with no proteins."""
        g1 = self.genome('GCA_000000001.1', {'g1': 'MKALIE', 'g2': 'MPQS'},
                         pfam=[('g1', 'PF00410.20,1e-5,20.0')], tigrfam=[('g2', 'TIGR00006,1e-9,30.0'),
                                                                         ('g1', 'TIGR00006,1e-2,10.0')])
        g2 = self.genome('GCA_000000002.1', {'g1': 'MKA'}, pfam=[], tigrfam=[])
        g3 = self.genome('GCA_000000003.1', None, pfam=[], tigrfam=[])
        genome_dirs = os.path.join(self.dir, 'genome_dirs.tsv')
        with open(genome_dirs, 'w') as handle:
            for accession, gdir in (('GCA_000000001.1', g1), ('GCA_000000002.1', g2), ('GCA_000000003.1', g3)):
                handle.write('{}\t{}\tG\n'.format(accession, gdir))
        cursor = FakeCursor([self.pf, self.tigr],
                            genomes=[('GCA_000000001.1', 11), ('GCA_000000002.1', 12), ('GCA_000000003.1', 13)])
        return genome_dirs, cursor

    def run_release(self, genome_dirs, cursor, **kwargs):
        manager = self.manager(cursor)
        out_dir = os.path.join(self.dir, 'out')
        with mock.patch.object(A, 'check_dependencies'), \
                mock.patch.object(A, 'record_program_version', return_value='HMMER 3.3 (Nov 2019)'), \
                self.assertLogs('timestamp', level='INFO') as logged:
            finished = manager.run([1, 2], False, genome_dirs, out_dir, **kwargs)
        return manager, out_dir, finished, [r.getMessage() for r in logged.records]

    def test_every_marker_of_every_genome_has_a_row_and_each_batch_is_its_own_commit(self):
        genome_dirs, cursor = self.release()
        manager, out_dir, finished, _ = self.run_release(genome_dirs, cursor)

        self.assertTrue(finished)
        rows = sorted(row for copied in cursor.copied for row in copied)
        self.assertEqual(rows, [
            ['11', '101', 'MKALIE', 'f', '1e-5', '20.0', '1', '1'],
            ['11', '202', 'MPQS', 't', '1e-9', '30.0', '2', '2'],
            # a marker no gene names is a row of gaps as long as the HMM
            ['12', '101', '------', 'f', '\\N', '\\N', '\\N', '\\N'],
            ['12', '202', '----', 'f', '\\N', '\\N', '\\N', '\\N']])
        # two batches of two genomes, each written and committed on its own
        self.assertEqual(manager.temp_con.commits, 2)
        inserts = [s for s in cursor.statements if s.startswith('INSERT INTO aligned_markers')]
        self.assertEqual(len(inserts), 2)
        self.assertIn('ON CONFLICT (genome_id, marker_id) DO UPDATE SET sequence = EXCLUDED.sequence', inserts[0])

    def test_a_genome_without_a_protein_file_is_counted_and_listed_in_the_output_directory(self):
        genome_dirs, cursor = self.release()
        _, out_dir, finished, messages = self.run_release(genome_dirs, cursor)

        self.assertTrue(finished)
        self.assertNotIn('13', [row[0] for copied in cursor.copied for row in copied])
        with open(os.path.join(out_dir, 'missing_protein_file.tsv')) as handle:
            self.assertEqual(handle.read().splitlines(), ['genome_id', 'GCA_000000003.1'])
        self.assertTrue(any(m.startswith('1 genome(s) have no protein file and could not be aligned')
                            for m in messages))

    def test_a_finished_batch_is_not_done_again(self):
        genome_dirs, cursor = self.release()
        self.run_release(genome_dirs, cursor)
        again = FakeCursor([self.pf, self.tigr], genomes=cursor.genomes)
        manager, _, finished, _ = self.run_release(genome_dirs, again)
        self.assertTrue(finished)
        self.assertEqual((again.copied, manager.temp_con.commits), ([], 0))

    def test_a_batch_that_fails_is_rolled_back_and_the_run_says_so(self):
        # a worker that died was passed over, its genomes never aligned and nothing said
        genome_dirs, cursor = self.release()
        os.remove(A.tophit_file(os.path.join(self.dir, 'genomes', 'GCA_000000001.1_ASM1v1', 'prodigal'),
                                'GCA_000000001.1', 'TIGR'))
        manager, out_dir, finished, messages = self.run_release(genome_dirs, cursor)

        self.assertFalse(finished)
        # one rollback ends reading the genomes, one undoes the failed batch
        self.assertEqual((manager.temp_con.commits, manager.temp_con.rollbacks), (1, 2))
        self.assertTrue(any('failed, nothing of it written' in m and 'GCA_000000001.1 has called proteins' in m
                            for m in messages))
        failed = [d for d in os.listdir(os.path.join(out_dir, 'marker_sets_1_2_new'))
                  if os.path.exists(os.path.join(out_dir, 'marker_sets_1_2_new', d, 'FAILED'))]
        self.assertEqual(failed, ['batch_000001'])

    def test_a_genome_the_genome_dirs_file_does_not_locate_refuses_the_run(self):
        genome_dirs, cursor = self.release()
        cursor.genomes.append(('GCA_000000099.1', 99))
        with self.assertRaisesRegex(A.AlignmentError, '1 genome.* not in .*GCA_000000099.1'):
            self.run_release(genome_dirs, cursor)
        self.assertEqual(cursor.copied, [])

    def test_new_genomes_are_those_with_no_row_for_any_of_the_markers(self):
        cursor = FakeCursor()
        manager = self.manager(cursor)
        manager.select_genomes([self.pf, self.tigr], all_genomes=False)
        self.assertIn('NOT EXISTS (SELECT 1 FROM aligned_markers am WHERE am.genome_id = g.id AND '
                      'am.marker_id = ANY(%s))', cursor.statements[-1])
        manager.select_genomes([self.pf], all_genomes=True)
        self.assertNotIn('aligned_markers', cursor.statements[-1])


class TheCommandLine(TempDirCase):
    def argv(self, *extra):
        return (['align_marker_genes', '--db_service', 'gtdb_r237', '--marker_set_ids', '1', '2',
                 '-g', 'genome_dirs.tsv', '-o', os.path.join(self.dir, 'out'), '-l', 'run.log'] + list(extra))

    def test_one_of_new_or_all_genomes_is_required_and_not_both(self):
        parser = main_module.get_main_parser()
        for extra in ((), ('--new_genomes', '--all_genomes')):
            with self.subTest(extra=extra), mock.patch('sys.stderr'), self.assertRaises(SystemExit) as ended:
                parser.parse_args(self.argv(*extra))
            self.assertEqual(ended.exception.code, 2)

    def test_the_command_is_handed_its_options(self):
        options = main_module.get_main_parser().parse_args(self.argv('--all_genomes', '-c', '8'))
        with mock.patch.object(main_py, 'MarkerAlignmentManager') as manager, \
                mock.patch.object(main_py, 'check_file_exists'):
            manager.return_value.run.return_value = True
            main_py.OptionsParser().parse_options(options)
        manager.assert_called_once_with({'service': 'gtdb_r237'}, 8, 1000, '/tmp', False, 2.0 * 60 * 60)
        manager.return_value.run.assert_called_once_with([1, 2], True, 'genome_dirs.tsv',
                                                         os.path.join(self.dir, 'out'))

    def test_a_failed_batch_or_an_alignment_error_exits_1(self):
        options = main_module.get_main_parser().parse_args(self.argv('--new_genomes'))
        for outcome in (dict(return_value=False), dict(side_effect=A.AlignmentError('no'))):
            with self.subTest(outcome=outcome), \
                    mock.patch.object(main_py, 'MarkerAlignmentManager') as manager, \
                    mock.patch.object(main_py, 'check_file_exists'), \
                    self.assertRaises(SystemExit) as ended:
                manager.return_value.run.configure_mock(**outcome)
                main_py.OptionsParser().parse_options(options)
            self.assertEqual(ended.exception.code, 1)


if __name__ == '__main__':
    unittest.main()
