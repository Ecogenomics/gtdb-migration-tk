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

from gtdb_migration_tk.biolib_lite.external.pypfam.Scan.PfamScan import PfamScan


class PfamSearch(object):
    """Runs PFAM search (Python equivalent of pfam_search.pl) over a set of genomes."""

    def __init__(self, pfam_hmm_dir: str) -> None:
        """Initialization.

        Parameters
        ----------
        pfam_hmm_dir : str
            Directory holding the Pfam HMMs to search against.
        """

        self.cpus_per_genome: int = 1
        self.pfam_hmm_dir: str = pfam_hmm_dir

    def run(self, gene_file: str, output_hit_file: str) -> None:
        """Search one genome's proteins against the Pfam HMMs.

        One genome per call, on one CPU: the caller runs a pool of these rather
        than giving any one search more than a core.

        Parameters
        ----------
        gene_file : str
            Amino acid FASTA of the genome's called genes, uncompressed.
        output_hit_file : str
            Where the hit table is written.

        @return: nothing; the hits are written to output_hit_file.
        """

        pfam_scan = PfamScan(cpu=self.cpus_per_genome, fasta=gene_file, dir=self.pfam_hmm_dir)
        pfam_scan.search()
        pfam_scan.write_results(output_hit_file, None, None, None, None)



