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
import multiprocessing as mp
import os

from gtdb_migration_tk.biolib_lite.external.pypfam.Scan.PfamScan import PfamScan


class PfamSearch(object):
    """Runs pfam_search.pl over a set of genomes.
    Copied from GTDB-Tk 2.1.1 but modified a lot!!"""

    def __init__(self,
                 pfam_hmm_dir):
        """Initialization."""

        self.cpus_per_genome = 1
        self.pfam_hmm_dir = pfam_hmm_dir

    def run(self, gene_file,output_hit_file):
        """Process each data item in parallel."""
        try:
            pfam_scan = PfamScan(cpu=self.cpus_per_genome, fasta=gene_file, dir=self.pfam_hmm_dir)
            pfam_scan.search()
            pfam_scan.write_results(output_hit_file, None, None, None, None)
        except Exception as error:
            raise error



