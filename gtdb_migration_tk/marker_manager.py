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
import argparse
import logging
import ntpath


from biolib.common import remove_extension
from biolib.external.execute import check_dependencies
from biolib.external.prodigal import Prodigal
from biolib.checksum import sha256

from gtdb_migration_tk.prettytable import PrettyTable

from gtdb_migration_tk.tools import Tools


class MarkerManager(object):
    """Create file indicating directory of each genome."""

    def __init__(self, tmp_dir='/tmp/', cpus=1):
        """Initialization."""

        self.tmp_dir = tmp_dir
        self.cpus = cpus

        check_dependencies(['hmmsearch'])

        self.pfam_hmm_dir = '/srv/db/pfam/27/'
        self.protein_file_ext = '_protein.faa'

        self.logger = logging.getLogger('timestamp')

    def run_hmmsearch(self, gtdb_genome_path_file, report, db):
        if db == 'pfam':
            self.pfam_search()

    def pfam_search(self, gtdb_genome_path_file, genome_report):
