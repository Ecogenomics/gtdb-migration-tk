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
import csv
import re
import logging
from collections import defaultdict

from biolib.common import check_file_exists, make_sure_path_exists, check_dir_exists

from gtdb_migration_tk.lpsn import LPSN
from gtdb_migration_tk.bacdive import BacDive
from gtdb_migration_tk.strains import Strains
from gtdb_migration_tk.tools import Tools


class OptionsParser():
    def __init__(self):
        """Initialization"""

        self.logger = logging.getLogger('timestamp')

    def full_lpsn_wf(self, options):
        """Full workflow to parse LPSN."""
        make_sure_path_exists(options.output_dir)
        p = LPSN(options.output_dir)
        p.full_lpsn_wf()

    def pull_html(self, options):
        """Pull all genus.html files."""
        make_sure_path_exists(options.output_dir)
        p = LPSN(options.output_dir)
        p.pull_html()

    def parse_html(self, options):
        """Parse all html files."""
        make_sure_path_exists(options.output_dir)
        p = LPSN(options.output_dir)
        p.parse_html(options.input_dir)

    def download_strains(self, options):
        make_sure_path_exists(options.output_dir)
        p = BacDive(options.output_dir, options.username, options.pwd)
        p.download_strains()

    def generate_date_table(self, options):
        p = Strains()
        p.generate_date_table(options.lpsn_species_info,
                              options.dsmz_species_info,
                              options.straininfo_species_info,
                              options.output_file)

    def generate_type_table(self, options):
        p = Strains(options.output_dir, options.cpus)
        p.generate_type_strain_table(options.metadata_file,
                                     options.ncbi_names,
                                     options.ncbi_nodes,
                                     options.lpsn_dir,
                                     options.dsmz_dir,
                                     options.straininfo_dir,
                                     options.year_table,
                                     options.source_strain)

    def compare_metadata(self, options):
        p = Tools()
        p.compare_metadata(options.previous_metadata_file,
                           options.new_metadata_file)

    def parse_options(self, options):
        """Parse user options and call the correct pipeline(s)"""
        if options.subparser_name == 'lpsn':
            if options.lpsn_subparser_name == 'lpsn_wf':
                self.full_lpsn_wf(options)
            elif options.lpsn_subparser_name == 'parse_html':
                self.parse_html(options)
            elif options.lpsn_subparser_name == 'pull_html':
                self.pull_html(options)
            else:
                self.logger.error('Unknown command: ' +
                                  options.lpsn_subparser_name + '\n')
        elif options.subparser_name == 'bacdive':
            if options.bacdive_subparser_name == 'download_strains':
                self.download_strains(options)
        elif options.subparser_name == 'strains':
            if options.strains_subparser_name == 'date_table':
                self.generate_date_table(options)
            if options.strains_subparser_name == 'type_table':
                self.generate_type_table(options)
        elif options.subparser_name == 'overview':
            self.compare_metadata(options)
        else:
            self.logger.error('Unknown command: ' +
                              options.subparser_name + '\n')
            sys.exit()

        return 0