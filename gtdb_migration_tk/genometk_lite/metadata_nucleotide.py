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

import logging

import biolib.seq_io as seq_io
import biolib.seq_tk as seq_tk
import biolib.genome_tk as genome_tk


class MetadataNucleotide():
    """Calculate metadata derived from nucleotide sequences."""

    def __init__(self):
        """Initialization"""
        self.logger = logging.getLogger('timestamp')

    def generate(self, genome_file, contig_break):
        """Derive metdata across nucleotide sequences.

        Parameters
        ----------
        genome_file : str
            Name of fasta file containing nucleotide sequences.
        contig_break : int
            Minimum number of ambiguous bases for defining contigs.

        Returns
        -------
        dict : d[metadata_field] -> value
            Map of metadata fields to their respective values.
        dict : d[metadata_field -> description
            Description of each metadata field.
        """

        # calculate nucleotide statistics
        scaffolds = seq_io.read(genome_file)

        nuc_stats = {}
        nuc_desc = {}

        nuc_stats['scaffold_count'] = len(scaffolds)
        nuc_desc['scaffold_count'] = "Number of scaffolds in genome."
        nuc_stats['gc_count'] = genome_tk.gc_count(scaffolds)
        nuc_desc['gc_count'] = "Number of G or C bases in genome."
        nuc_stats['gc_percentage'] = genome_tk.gc(scaffolds) * 100.0
        nuc_desc['gc_percentage'] = "GC content of genome."
        nuc_stats['genome_size'] = sum([len(x) for x in scaffolds.values()])
        nuc_desc['genome_size'] = "Total base pairs in genome including nucleotide bases, ambiguous bases, and gaps."
        nuc_stats['n50_scaffolds'] = seq_tk.N50(scaffolds)
        nuc_desc['n50_scaffolds'] = "Scaffold length at which 50% of total bases in assembly are in scaffolds of that length or greater."
        nuc_stats['l50_scaffolds'] = seq_tk.L50(scaffolds, nuc_stats['n50_scaffolds'])
        nuc_desc['l50_scaffolds'] = "Number of scaffolds longer than, or equal to, the scaffold N50 length."
        nuc_stats['mean_scaffold_length'] = int(seq_tk.mean_length(scaffolds))
        nuc_desc['mean_scaffold_length'] = "Mean length of scaffolds in base pairs."
        nuc_stats['longest_scaffold'] = seq_tk.max_length(scaffolds)
        nuc_desc['longest_scaffold'] = "Number of bases in longest scaffold."

        contigs = seq_tk.identify_contigs(scaffolds, 'N' * contig_break)
        nuc_stats['contig_count'] = len(contigs)
        nuc_desc['contig_count'] = "Number of contigs in genome."
        nuc_stats['ambiguous_bases'] = genome_tk.ambiguous_nucleotides(contigs)
        nuc_desc['ambiguous_bases'] = "Number of ambiguous bases in contigs."
        nuc_stats['total_gap_length'] = genome_tk.ambiguous_nucleotides(scaffolds) - nuc_stats['ambiguous_bases']
        nuc_desc['total_gap_length'] = "Number of ambiguous bases comprising gaps in scaffolds."
        nuc_stats['n50_contigs'] = seq_tk.N50(contigs)
        nuc_desc['n50_contigs'] = "Contig length at which 50% of total bases in assembly are in contigs of that length or greater."
        nuc_stats['l50_contigs'] = seq_tk.L50(contigs, nuc_stats['n50_contigs'])
        nuc_desc['l50_contigs'] = "Number of contigs longer than, or equal to, the contig N50 length."
        nuc_stats['mean_contig_length'] = int(seq_tk.mean_length(contigs))
        nuc_desc['mean_contig_length'] = "Mean length of contigs in base pairs."
        nuc_stats['longest_contig'] = seq_tk.max_length(contigs)
        nuc_desc['longest_contig'] = "Number of bases in longest contig."

        return nuc_stats, nuc_desc
