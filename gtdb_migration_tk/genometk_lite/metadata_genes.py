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

from numpy import (zeros as np_zeros,
                   sum as np_sum)


class GenericFeatureParser():
    """Parses generic feature file (GFF)."""
    def __init__(self, filename):
        self.cds_count = 0
        self.tRNA_count = 0
        self.rRNA_count = 0
        self.rRNA_16S_count = 0
        self.ncRNA_count = 0

        self.genes = {}
        self.last_coding_base = {}

        self._parse(filename)

        self.coding_mask = {}
        for seq_id in self.genes:
            self.coding_mask[seq_id] = self._coding_mask(seq_id)

    def _parse(self, gff_file):
        """Parse GFF file.

        Parameters
        ----------
        gff_file : str
            Generic feature file to parse.
        """

        for line in open(gff_file):
            if line[0] == '#':
                continue

            line_split = line.split('\t')
            if line_split[2] == 'tRNA':
                self.tRNA_count += 1
            elif line_split[2] == 'rRNA':
                self.rRNA_count += 1

                if 'product=16S ribosomal RNA' in line_split[8]:
                    self.rRNA_16S_count += 1
            elif line_split[2] == 'ncRNA':
                self.ncRNA_count += 1
            elif line_split[2] == 'CDS':
                self.cds_count += 1

                seq_id = line_split[0]
                if seq_id not in self.genes:
                    gene_count = 0
                    self.genes[seq_id] = {}
                    self.last_coding_base[seq_id] = 0

                gene_id = seq_id + '_' + str(gene_count)
                gene_count += 1

                start = int(line_split[3])
                end = int(line_split[4])

                self.genes[seq_id][gene_id] = [start, end]
                self.last_coding_base[seq_id] = max(self.last_coding_base[seq_id], end)

    def _coding_mask(self, seq_id):
        """Build mask indicating which bases in a sequences are coding."""

        # safe way to calculate coding bases as it accounts
        # for the potential of overlapping genes
        coding_mask = np_zeros(self.last_coding_base[seq_id])
        for pos in self.genes[seq_id].values():
            coding_mask[pos[0]:pos[1] + 1] = 1

        return coding_mask

    def coding_bases(self, seq_id):
        """Calculate number of coding bases in sequence."""

        # check if sequence has any genes
        if seq_id not in self.genes:
            return 0

        return np_sum(self.coding_mask[seq_id])

    def total_coding_bases(self):
        """Calculate total number of coding bases.

        Returns
        -------
        int
            Number of coding bases.
        """

        coding_bases = 0
        for seq_id in self.genes:
            coding_bases += self.coding_bases(seq_id)

        return int(coding_bases)


class MetadataGenes():
    """Calculate metadata derived from called genes."""

    def __init__(self):
        """Initialization"""
        self.logger = logging.getLogger('timestamp')

    def generate(self, genome_file, gff_file):
        """Derive metdata from gene sequences.

        Parameters
        ----------
        genome_file : str
            Name of fasta file containing nucleotide sequences.
        gff_file : str
            Name of generic feature file describing genes.

        Returns
        -------
        dict : d[metadata_field] -> value
            Map of metadata fields to their respective values.
        dict : d[metadata_field -> description
            Description of each metadata field.
        """

        gff_parser = GenericFeatureParser(gff_file)
        coding_bases = gff_parser.total_coding_bases()

        # calculate nucleotide statistics
        scaffolds = seq_io.read(genome_file)
        genome_size = sum([len(x) for x in scaffolds.values()])

        gene_stats = {}
        gene_desc = {}

        gene_stats['protein_count'] = gff_parser.cds_count
        gene_desc['protein_count'] = "Number of protein coding genes."

        gene_stats['coding_bases'] = coding_bases
        gene_desc['coding_bases'] = "Number of coding bases in genome."

        gene_stats['coding_density'] = float(coding_bases) * 100.0 / genome_size
        gene_desc['coding_density'] = "Percentage of coding bases in genome."

        # This quantities are not being reported here as they are only specified
        # in the NCBI GFF files. A general GFF file (e.g., for prodigal) only
        # specifies identified CDS.
        # gene_stats['tRNA_count'] = gff_parser.tRNA_count
        # gene_desc['tRNA_count'] = "Number of tRNA genes identified in genome."
        # gene_stats['ncRNA_count'] = gff_parser.ncRNA_count
        # gene_desc['ncRNA_count'] = "Number of ncRNA genes identified in genome."
        # gene_stats['rRNA_count'] = gff_parser.rRNA_count
        # gene_desc['rRNA_count'] = "Number of rRNA genes identified in genome."
        # gene_stats['16S_count'] = gff_parser.rRNA_16S_count
        # gene_desc['16S_count'] = "Number of 16S rRNA genes identified in genome."

        return gene_stats, gene_desc
