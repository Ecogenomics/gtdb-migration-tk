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

import os
import sys
import logging
import subprocess
from typing import Dict, List, Optional
from collections import defaultdict
from dataclasses import dataclass

import gtdb_migration_tk.biolib_lite.seq_io as seq_io
from gtdb_migration_tk.biolib_lite.taxonomy import Taxonomy
from gtdb_migration_tk.biolib_lite.external.blast import Blast
from gtdb_migration_tk.biolib_lite.common import make_sure_path_exists


@dataclass
class Hit:
    score: float
    evalue: float
    hmm_from: int
    hmm_to: int
    ali_from: int
    ali_to: int
    align_len: int
    rev_comp: bool


DNA_COMPLEMENTS = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}


def rev_comp(seq: str) -> str:
    """Reverse complement a sequence."""

    complement = [DNA_COMPLEMENTS.get(base, base) for base in reversed(seq)]
    return "".join(complement)


class RNA(object):
    """Identify, extract, and classify rRNA genes."""

    def __init__(self, rna_name: str, domain: str, cpus: int):
        """Initialization."""

        self.logger = logging.getLogger('timestamp')

        self.rna_name = rna_name
        self.domain = domain
        self.cpus = cpus

        # E-value threshold for defining valid hits (taken from barrnap 0.9)
        self.evalue_threshold = 1e-6

        # Minmum length of rRNA gene based on minimum proportion of gene length.
        # (minimum proportion of 0.25 and gene lengths taken from barrnap 0.9)
        len_16S_rRNA = 1585
        len_18S_rRNA = 1869
        len_23S_rRNA = 3232
        len_5s_rRNA = 119

        min_len_by_domain = {}
        if rna_name == 'ssu':
            min_len_by_domain['ar'] = 0.25 * len_16S_rRNA
            min_len_by_domain['bac'] = 0.25 * len_16S_rRNA
            min_len_by_domain['euk'] = 0.25 * len_18S_rRNA
        elif rna_name == 'lsu_23S':
            min_len_by_domain['ar'] = 0.25 * len_23S_rRNA
            min_len_by_domain['bac'] = 0.25 * len_23S_rRNA
            min_len_by_domain['euk'] = 0.25 * len_23S_rRNA
        elif rna_name == 'lsu_5S':
            min_len_by_domain['ar'] = 0.25 * len_5s_rRNA
            min_len_by_domain['bac'] = 0.25 * len_5s_rRNA
            min_len_by_domain['euk'] = 0.25 * len_5s_rRNA
        else:
            self.logger.error(f'Unrecognized RNA gene name: {rna_name}')
            sys.exit(1)

        self.min_len = min_len_by_domain[domain]
        self.max_window_len = int(1.2 * len_23S_rRNA)

    def _hmm_search(self, seq_file: str, output_dir: str) -> None:
        """Identify rRNA genes using nhmmer."""

        output_prefix = os.path.join(output_dir, self.rna_name)
        tblout_file = f"{output_prefix}.{self.domain}.txt"

        cmd = [
            "nhmmer",
            "--noali",
            "--cpu", str(self.cpus),
            "-E", str(self.evalue_threshold),
            "--w_length", str(self.max_window_len),
            "-o", os.devnull,
            "--tblout", tblout_file,
            self.hmm_model_file,
            "-"
        ]

        # handle file decompression using zcat if required
        cat_cmd = ["zcat" if seq_file.endswith(".gz") else "cat", seq_file]

        try:
            # start the cat/zcat process
            with subprocess.Popen(cat_cmd, stdout=subprocess.PIPE) as cat_proc:
                # pipe output into nhmmer
                with subprocess.Popen(cmd, stdin=cat_proc.stdout, stderr=subprocess.PIPE) as nhmmer_proc:
                    # allow cat_proc to receive a SIGPIPE if nhmmer_proc exits early
                    if cat_proc.stdout is not None:
                        cat_proc.stdout.close()

                    _, stderr = nhmmer_proc.communicate()

                    if nhmmer_proc.returncode != 0:
                        error_msg = stderr.decode().strip()
                        self.logger.error(f"nhmmer failed with return code {nhmmer_proc.returncode}: {error_msg}")
                        raise RuntimeError(f"nhmmer error: {error_msg}")

        except FileNotFoundError:
            raise FileNotFoundError(f"Required executable (cat/zcat or nhmmer) not found in PATH.")
        except Exception as e:
            self.logger.error(f"An unexpected error occurred during HMM search: {e}")
            raise

    def _read_hits(self, results_file: str) -> Dict[str, List[Hit]]:
        """Parse hits from nhmmer output.

        Parameters
        ----------
        results_file : str
            Output file from nhmmer to parse.

        Returns
        -------
        dict : d[seq_id] -> information about hit
            Information about hits for individual sequences.
        """

        seq_info = defaultdict(list)

        seq_id = None
        with open(results_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue

                tokens = line.strip().split()

                seq_id = tokens[0]

                hmm_from = int(tokens[4])
                hmm_to = int(tokens[5])
                assert hmm_from < hmm_to

                ali_from = int(tokens[6])
                ali_to = int(tokens[7])

                rev_comp = False
                if ali_from > ali_to:
                    rev_comp = True
                    ali_from, ali_to = ali_to, ali_from

                assert ali_from < ali_to

                align_len = int(ali_to) - int(ali_from) + 1

                evalue = float(tokens[12])
                score = float(tokens[13])

                if float(evalue) <= self.evalue_threshold:
                    seq_info[seq_id].append(Hit(
                        score, evalue, hmm_from, hmm_to, ali_from, ali_to, align_len, rev_comp
                    ))

        return seq_info

    def _add_hit(self, seq_hits: Dict[str, Hit], seq_id: str, cur_hit: Hit) -> None:
        """Add hits from individual HMMs generating unique sequence names.

        Parameters
        ----------
        seq_hits : d[seq_id] -> information about hit
            Information about hits for individual sequences.
        seq_id : str
            Sequence identifier with hit to add.
        cur_hit : Hit
            Information about hit to sequence.
        """

        # check if this is the first hit to this sequence
        if seq_id not in seq_hits:
            seq_hits[seq_id] = cur_hit
            return

        # add hit generating a unique sequence name
        base_seq_id = seq_id
        index = 1
        while (True):
            new_seq_id = base_seq_id + '-#' + str(index)
            index += 1

            # hits are not close enough to concatenate
            if new_seq_id not in seq_hits:
                seq_hits[new_seq_id] = cur_hit
                break

    def _identify(self, genome_file: str, output_dir: str) -> Dict[str, Hit]:
        """Identify rRNA genes.

        Parameters
        ----------
        genome_file : str
            Name of fasta file containing nucleotide sequences.
        output_dir : str
            Output directory.

        Returns
        -------
        dict : d[seq_id] -> information about best hit
            Information about best hits.
        """

        # identify rRNAs on contigs/scaffolds
        self._hmm_search(genome_file, output_dir)

        # read HMM hits
        hmm_results_file = os.path.join(output_dir, f'{self.rna_name}.{self.domain}.txt')
        if not os.path.exists(hmm_results_file):
            return {}

        seq_info = self._read_hits(hmm_results_file)

        hits = {}
        for seq_id, seq_hits in seq_info.items():
            for hit in seq_hits:
                self._add_hit(hits, seq_id, hit)

        # filter for min size
        filtered_best_hits = {}
        for seq_id, seq_info in hits.items():
            if seq_info.align_len >= self.min_len:
                filtered_best_hits[seq_id] = seq_info

        return filtered_best_hits

    def _extract(self, genome_file: str, best_hits: Dict[str, Hit], output_dir: str) -> str:
        """Extract rRNA genes.

        Parameters
        ----------
        genome_file : str
            Name of FASTA file containing nucleotide sequences.
        best_hits : d[seq_id] -> information about best hit
            Information about best hits.
        output_dir : str
            Output directory.

        Returns
        -------
        str
            Name of FASTA file containing extractracted sequences.
        """

        # write summary file and putative SSU rRNAs to file
        summary_file = os.path.join(
            output_dir, '%s.hmm_summary.tsv' % self.rna_name)
        summary_out = open(summary_file, 'w')
        summary_out.write('Sequence Id\tHMM\tE-value\tHMM start\tHMM end')
        summary_out.write('\tSequence start\tSequence end\tGene length\tReverse Complement\tSequence length\n')

        ssu_seq_file = os.path.join(output_dir, '%s.fna' % self.rna_name)
        seq_out = open(ssu_seq_file, 'w')

        seqs = seq_io.read(genome_file)

        for seq_id in best_hits:
            orig_seq_id = seq_id
            if '-#' in seq_id:
                seq_id = seq_id[0:seq_id.rfind('-#')]

            hit = best_hits[orig_seq_id]

            seq = seqs[seq_id]
            summary_out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                orig_seq_id,
                self.domain,
                hit.evalue,
                hit.hmm_from,
                hit.hmm_to,
                hit.ali_from,
                hit.ali_to,
                hit.align_len,
                str(hit.rev_comp),
                len(seq)))

            seq_out.write(f'>{orig_seq_id}\n')

            ssu_seq = seq[hit.ali_from - 1:hit.ali_to]
            if hit.rev_comp == True:
                ssu_seq = rev_comp(ssu_seq)

            seq_out.write(ssu_seq + '\n')

        summary_out.close()
        seq_out.close()

        return ssu_seq_file

    def classify(self, seq_file: str, db: str, taxonomy_file: str, output_dir: str) -> None:
        """Classify rRNA genes.

        Parameters
        ----------
        seq_file : str
            Name of fasta file containing rRNA sequences.
        ssu_db : str
            BLAST database of rRNA genes.
        ssu_taxonomy_file : str
            Taxonomy file for genes in the rRNA database.
        output_dir : str
            Output directory.
        """

        # blast sequences against rRNA database
        blast = Blast(self.cpus)
        blast_file = os.path.join(output_dir, f'{self.rna_name}.blastn.tsv')
        blast.blastn(seq_file,
                     db,
                     blast_file,
                     evalue=self.evalue_threshold,
                     max_matches=5,
                     output_fmt='custom')

        # read taxonomy file
        taxonomy = Taxonomy().read(taxonomy_file)

        # write out classification file
        classification_file = os.path.join(
            output_dir, f'{self.rna_name}.taxonomy.tsv')
        fout = open(classification_file, 'w')
        fout.write('query_id\ttaxonomy\tlength\tblast_subject_id')
        fout.write('\tblast_evalue\tblast_bitscore\tblast_align_len\tblast_perc_identity\n')

        processed_query_ids = set()
        for line in open(blast_file):
            line_split = [x.strip() for x in line.split('\t')]
            query_id = line_split[0]

            if query_id in processed_query_ids:
                # A query may have multiple hits to different genes or sections
                # of a gene. BLAST results are organized by bitscore so
                # only the first hit is considered.
                continue

            processed_query_ids.add(query_id)
            query_len = int(line_split[1])
            subject_id = line_split[2]
            align_len = line_split[5]
            perc_identity = line_split[6]
            evalue = line_split[7]
            bitscore = line_split[8]

            taxonomy_str = ';'.join(taxonomy[subject_id])

            fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                query_id,
                taxonomy_str,
                query_len,
                subject_id,
                evalue,
                bitscore,
                align_len,
                perc_identity))

        fout.close()

    def run(self,
            genome_file: str,
            hmm_model_file: str,
            db: Optional[str],
            taxonomy_file: Optional[str],
            output_dir: str) -> None:
        """Identify and extract rRNA gene.

        Parameters
        ----------
        genome_file : str
            Name of fasta file containing nucleotide sequences.
        rna_name : str
            Name of RNA gene for reporting results and prefixing files.
        hmm_model_file : str
            File with HMM model for RNA gene.
        db : str
            BLAST database for classifying rRNA genes.
        taxonomy_file : str
            Taxonomy file for classifying rRNA genes.
        output_dir : str
            Output directory.
        """

        self.hmm_model_file = hmm_model_file

        make_sure_path_exists(output_dir)

        best_hits = self._identify(
            genome_file,
            output_dir)

        if len(best_hits) > 0:
            seq_file = self._extract(genome_file, best_hits, output_dir)

            if db is not None and taxonomy_file is not None and seq_file:
                self.classify(seq_file,
                              db,
                              taxonomy_file,
                              output_dir)

        canary_file = os.path.join(output_dir, self.rna_name + '.canary.txt')
        with open(canary_file, 'w') as filehandle:
            filehandle.write('done.\n')
