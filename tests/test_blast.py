#!/usr/bin/env python3
"""Offline unit tests for the BLAST wrapper, biolib_lite/external/blast.py.

Run with the interpreter that has tqdm:

    /opt/centos7/sw/miniconda3/envs/gtdb_migration_tk-r237/bin/python -m unittest discover -s tests -p test_blast.py

blastn, blastp and makeblastdb are each a few lines of shell put first on PATH,
which write what they are told to and exit with the status they are told to, so
what is tested is what the wrapper does with an exit status: until 0.1.36 it ran
them through os.system() and did nothing with it.
"""

import os
import shutil
import stat
import tempfile
import unittest

from gtdb_migration_tk.biolib_lite.external import blast as BL
from gtdb_migration_tk.genometk_lite.rna import RNA


# What each stand-in program does: writes a line to the file after -out, if it
# was given one, then says something on stderr and exits with the status given.
PROGRAM = '''#!/bin/sh
out=""
while [ $# -gt 0 ]; do
    if [ "$1" = "-out" ]; then out="$2"; fi
    shift
done
if [ -n "$out" ]; then printf 'contig_1\\t1500\\tLTP_1\\tx\\t1500\\t1500\\t99.9\\t0\\t2700\\n' > "$out"; fi
echo "{stderr}" >&2
exit {status}
'''


class BlastCase(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='blast_test.')
        self.bin = os.path.join(self.dir, 'bin')
        os.mkdir(self.bin)
        self.addCleanup(os.environ.__setitem__, 'PATH', os.environ['PATH'])
        os.environ['PATH'] = self.bin + os.pathsep + os.environ['PATH']

        for name in ('blastn', 'blastp', 'makeblastdb'):
            self.program(name)

        self.query = os.path.join(self.dir, 'ssu.fna')
        with open(self.query, 'w') as handle:
            handle.write('>contig_1\nACGT\n')
        self.out = os.path.join(self.dir, 'hits.tsv')
        self.db = os.path.join(self.dir, 'db')

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def program(self, name, status=0, stderr=''):
        path = os.path.join(self.bin, name)
        with open(path, 'w') as handle:
            handle.write(PROGRAM.format(status=status, stderr=stderr))
        os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR)


class WhenBlastnSucceeds(BlastCase):

    def test_its_table_is_left_to_be_read(self):
        BL.Blast(1).blastn(self.query, 'db.fna', self.out, output_fmt='custom')

        with open(self.out) as handle:
            self.assertTrue(handle.read().startswith('contig_1\t'))


class WhenBlastnFails(BlastCase):
    """os.system() returned the exit status to nobody, so a blastn that could not
    open its database was read back as a genome with no hits."""

    def setUp(self):
        super().setUp()
        self.program('blastn', status=2,
                     stderr='BLAST Database error: No alias or index file found')

    def test_it_raises_with_what_blastn_said(self):
        with self.assertRaises(BL.BlastError) as caught:
            BL.Blast(1).blastn(self.query, 'db.fna', self.out)

        self.assertIn('status 2', str(caught.exception))
        self.assertIn('No alias or index file found', str(caught.exception))

    def test_what_it_wrote_before_failing_is_removed(self):
        """A partial table read as results is a genome classified wrongly."""
        with self.assertRaises(BL.BlastError):
            BL.Blast(1).blastn(self.query, 'db.fna', self.out)

        self.assertFalse(os.path.exists(self.out))

    def test_it_is_a_runtime_error_the_rrna_workers_already_catch(self):
        self.assertTrue(issubclass(BL.BlastError, RuntimeError))

    def test_an_rrna_classification_is_not_written_from_a_failed_search(self):
        """RNA.classify() went on to write a taxonomy table with no rows, and
        the rRNA commands a canary saying the genome was classified."""
        out_dir = os.path.join(self.dir, 'rna')
        os.mkdir(out_dir)
        taxonomy = os.path.join(self.dir, 'taxonomy.tsv')
        with open(taxonomy, 'w') as handle:
            handle.write('LTP_1\td__Bacteria;p__;c__;o__;f__;g__;s__\n')

        with self.assertRaises(BL.BlastError):
            RNA('ssu', 'bac', 1).classify(self.query, 'db.fna', taxonomy, out_dir)

        self.assertFalse(os.path.exists(os.path.join(out_dir, 'ssu.taxonomy.tsv')))


class TheOtherPrograms(BlastCase):
    """Nothing calls these today; they failed silently in the same way."""

    def test_blastp_raises_where_it_fails(self):
        self.program('blastp', status=1)

        with self.assertRaises(BL.BlastError):
            BL.Blast(1).blastp(self.query, 'db.faa', self.out)

    def test_makeblastdb_raises_where_it_fails(self):
        self.program('makeblastdb', status=1, stderr='FASTA-Reader: Ignoring invalid residues')

        for create in ('create_blastn_db', 'create_blastp_db'):
            with self.assertRaises(BL.BlastError, msg=create):
                getattr(BL.Blast(1, silent=True), create)(self.query, self.db)

    def test_makeblastdb_succeeding_raises_nothing(self):
        BL.Blast(1, silent=True).create_blastn_db(self.query, self.db)

    def test_a_program_that_cannot_be_run_raises_rather_than_exiting(self):
        with self.assertRaises(BL.BlastError):
            BL.run_blast_command([os.path.join(self.dir, 'no_such_program')])


if __name__ == '__main__':
    unittest.main()
