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

__author__ = "Donovan Parks"
__copyright__ = "Copyright 2015"
__credits__ = ["Donovan Parks"]
__license__ = "GPL3"
__maintainer__ = "Donovan Parks"
__email__ = "donovan.parks@gmail.com"
__status__ = "Development"


import os
import gzip
import logging
import re
import subprocess
from collections import namedtuple

from gtdb_migration_tk.biolib_lite.external.execute import check_dependencies


class BlastError(RuntimeError):
    """A BLAST program exited with an error.

    Every program here was run through os.system(), which returns the exit
    status to nobody, so a blastn that could not open its database left an empty
    or partial table behind and the caller went on to read it as a genome with
    no hits -- rna_silva and rna_ltp then wrote a canary saying the genome was
    classified. A RuntimeError, so the callers that already name a genome whose
    search failed and carry on catch it without knowing about BLAST.
    """


def run_blast_command(cmd, output_file=None, quiet=False):
    """Run one BLAST program, raising where it fails.

    Parameters
    ----------
    cmd : list of str
        The program and its arguments, run without a shell, so no path or
        -outfmt string needs quoting.
    output_file : str
        File the program writes its results to. Removed where the program
        fails, since what it wrote before failing is not a result.
    quiet : bool
        Discard what the program prints on stdout.

    @return: None

    Raises
    ------
    BlastError
        The program exited non-zero or was killed, with what it said on stderr.
    """

    try:
        proc = subprocess.run(cmd,
                              stdout=subprocess.DEVNULL if quiet else None,
                              stderr=subprocess.PIPE)
    except OSError as error:
        raise BlastError('{} could not be run: {}'.format(cmd[0], error))

    if proc.returncode != 0:
        if output_file is not None and os.path.exists(output_file):
            os.remove(output_file)
        stderr = proc.stderr.decode('utf-8', 'replace').strip()
        raise BlastError('{} failed with status {}: {}'.format(
            cmd[0], proc.returncode, stderr[:500] or 'no message'))


def get_blastn_version():
    """Returns the version of blastn on the system path.

    Returns
    -------
    str
        The string containing the blastn version.
    """
    try:
        proc = subprocess.Popen(['blastn', '-version'],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                encoding='utf-8')
        stdout, _stderr = proc.communicate()
        version = re.search(r'blastn: (\S+)', stdout)
        if version:
            return version.group(1)
        else:
            return 'unknown'
    except Exception as e:
        print(e)
        return 'unknown'


class Blast():
    """Wrapper for running blast."""

    def __init__(self, cpus, silent=False):
        """Initialization.

        Parameters
        ----------
        cpus : int
            Number of cpus to use.
        """

        self.logger = logging.getLogger('timestamp')

        check_dependencies(['blastn', 'blastp', 'makeblastdb'])

        self.cpus = cpus
        self.silent = silent

        self.output_fmt = {'standard': '6',
                           'custom': '6 qseqid qlen sseqid stitle slen length pident evalue bitscore'}
        self.blastp_tasks = {'blastp', 'blastp-fast', 'blastp-short'}
        self.blastn_tasks = {'blastn', 'blastn-short', 'dc-megablast', 'megablast', 'rmblastn'}

        self.BlastHit = namedtuple('BlastHit', """query_id
                                                subject_id
                                                perc_identity
                                                aln_length
                                                mismatch_count
                                                gap_open_count
                                                query_start
                                                query_end
                                                subject_start
                                                subject_end
                                                evalue
                                                bitscore""")

        self.BlastHitCustom = namedtuple('BlastHitCustom', """query_id
                                                                query_len
                                                                subject_id
                                                                subject_annotation
                                                                subject_len
                                                                alignment_len
                                                                perc_identity
                                                                evalue
                                                                bitscore""")

        self.BlastHitHomologs = namedtuple('BlastHitHomologs', """query_id
                                                                subject_id
                                                                subject_annotation
                                                                perc_identity
                                                                query_perc_aln_len
                                                                subject_perc_aln_len
                                                                evalue
                                                                bitscore""")

    def blastp_cmd(self, query_seqs, prot_db, output_file, evalue=1e-3, max_matches=500, output_fmt='standard', task='blastp'):
        """Get BLASTp command.

        Finds homologs to query sequences using blastp homology search
        against a protein database. Hit can be reported using  either
        the 'standard' table 6 format or the following 'custom' format:
            qseqid qlen sseqid slen length pident evalue bitscore


        Parameters
        ----------
        query_seqs : str
            File containing query sequences.
        prot_db : str
            File containing blastp formatted database.
        output_file : str
            Output file containing blastp results.
        evalue : float
            E-value threshold used to identify homologs.
        max_matches : int
            Maximum hits per query sequence.
        output_fmt : str
            Specified output format of blast table: standard or custom.

        @return: the command as a list of arguments, for run_blast_command().
        """

        assert (output_fmt in self.output_fmt.keys())
        assert (task in self.blastp_tasks)

        return ['blastp', '-num_threads', str(self.cpus),
                '-query', query_seqs, '-db', prot_db, '-out', output_file,
                '-evalue', '%g' % evalue,
                '-max_target_seqs', str(max_matches),
                '-task', task,
                '-outfmt', self.output_fmt[output_fmt]]

    def blastp(self, query_seqs, prot_db, output_file, evalue=1e-3, max_matches=500, output_fmt='standard', task='blastp'):
        """Run BLASTp command, raising BlastError where it fails."""

        run_blast_command(self.blastp_cmd(query_seqs,
                                          prot_db,
                                          output_file,
                                          evalue,
                                          max_matches,
                                          output_fmt,
                                          task),
                          output_file)

    def blastn(self, query_seqs, nucl_db, output_file, evalue=1e-3, max_matches=500, output_fmt='standard', task='megablast'):
        """Apply blastn to query file.

        Finds homologs to query sequences using blastn homology search
        against a nucleotide database. Hit can be reported using  either
        the 'standard' table 6 format or the following 'custom' format:
            qseqid qlen sseqid slen length pident evalue bitscore


        Parameters
        ----------
        query_seqs : str
            File containing query sequences.
        nucl_db : str
            File containing blastn formatted database.
        output_file : str
            Output file containing blastn results.
        evalue : float
            E-value threshold used to identify homologs.
        max_matches : int
            Maximum hits per query sequence.
        output_fmt : str
            Specified output format of blast table: standard or custom.

        @return: None

        Raises
        ------
        BlastError
            blastn failed; output_file is removed rather than left to be read
            as a query with no hits.
        """

        assert (output_fmt in self.output_fmt.keys())
        assert (task in self.blastn_tasks)

        cmd = ['blastn', '-num_threads', str(self.cpus),
               '-query', query_seqs, '-db', nucl_db, '-out', output_file,
               '-evalue', '%g' % evalue,
               '-max_target_seqs', str(max_matches),
               '-task', task,
               '-outfmt', self.output_fmt[output_fmt]]
        run_blast_command(cmd, output_file)

    def create_blastn_db_cmd(self, prot_file, db_file):
        """Get command to create nucleotide database, as a list of arguments."""

        return ['makeblastdb', '-dbtype', 'nucl', '-in', prot_file, '-out', db_file]

    def create_blastn_db(self, prot_file, db_file):
        """Create nucleotide database, raising BlastError where it fails."""

        run_blast_command(self.create_blastn_db_cmd(prot_file, db_file),
                          quiet=self.silent)

    def create_blastp_db_cmd(self, prot_file, db_file):
        """Get command to create protein database, as a list of arguments."""

        return ['makeblastdb', '-dbtype', 'prot', '-in', prot_file, '-out', db_file]

    def create_blastp_db(self, prot_file, db_file):
        """Create protein database, raising BlastError where it fails."""

        run_blast_command(self.create_blastp_db_cmd(prot_file, db_file),
                          quiet=self.silent)

    def read_hit(self, table, table_fmt):
        """Generator function to read hits from a blast output table.

        Parameters
        ----------
        table : str
            Name of table to read.
        table_fmt : str
            Specified output format of blast table: standard or custom.

        Yields
        ------
        namedtuple
            Information about blast hit.
        """

        assert (table_fmt in self.output_fmt)

        if table.endswith('.gz'):
            open_file = gzip.open
        else:
            open_file = open

        if table_fmt == 'standard':
            with open_file(table, 'rt') as f:
                for line in f:
                    line_split = line.split('\t')
                    hit = self.BlastHit(query_id=line_split[0],
                                        subject_id=line_split[1],
                                        perc_identity=float(line_split[2]),
                                        aln_length=int(line_split[3]),
                                        mismatch_count=int(line_split[4]),
                                        gap_open_count=int(line_split[5]),
                                        query_start=int(line_split[6]),
                                        query_end=int(line_split[7]),
                                        subject_start=int(line_split[8]),
                                        subject_end=int(line_split[9]),
                                        evalue=float(line_split[10]),
                                        bitscore=float(line_split[11]))

                    yield hit
        else:
            with open_file(table, 'rt') as f:
                for line in f:
                    line_split = line.split('\t')
                    hit = self.BlastHitCustom(query_id=line_split[0],
                                              query_len=int(line_split[1]),
                                              subject_id=line_split[2],
                                              subject_annotation=line_split[3],
                                              subject_len=int(line_split[4]),
                                              alignment_len=int(line_split[5]),
                                              perc_identity=float(line_split[6]),
                                              evalue=float(line_split[7]),
                                              bitscore=float(line_split[8]))

                    yield hit

    def identify_homologs(self,
                          custom_blast_table,
                          evalue_threshold,
                          perc_identity_threshold,
                          perc_aln_len_threshold):
        """Identify homologs among blast hits based on specified criteria.

        Parameters
        ----------
        custom_blast_table : str
            File containing blast hits in the custom tabular format.
        evalue_threshold : float
            E-value threshold used to define homologs.
        perc_identity_threshold : float
            Percent identity threshold used to define a homologs.
        perc_aln_len_threshold : float
            Alignment length threshold used to define a homologs.

        Returns
        -------
        dict : d[subject_id] -> BlastHitCustom named tuple
            Dictionary with information about blast hits to homologs.
        """

        homologs = {}
        for line in open(custom_blast_table):
            line_split = line.split('\t')

            query_id = line_split[0]
            query_len = int(line_split[1])
            subject_id = line_split[2]
            subject_title = line_split[3]
            subject_len = int(line_split[4])
            align_len = int(line_split[5])
            perc_identity = float(line_split[6])
            evalue = float(line_split[7])
            bitscore = float(line_split[8])

            if evalue <= evalue_threshold and perc_identity >= perc_identity_threshold:
                query_perc_aln_len = align_len * 100.0 / query_len
                subject_perc_aln_len = align_len * 100.0 / subject_len

                if query_perc_aln_len >= perc_aln_len_threshold and subject_perc_aln_len >= perc_aln_len_threshold:
                    prev_hit = homologs.get(subject_id, None)
                    if not prev_hit or bitscore > prev_hit.bitscore:
                        homologs[subject_id] = self.BlastHitHomologs(query_id=query_id,
                                                                     subject_id=subject_id,
                                                                     subject_annotation=subject_title,
                                                                     perc_identity=perc_identity,
                                                                     query_perc_aln_len=query_perc_aln_len,
                                                                     subject_perc_aln_len=subject_perc_aln_len,
                                                                     evalue=evalue,
                                                                     bitscore=bitscore)

        return homologs
