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

__author__ = 'Pierre Chaumeil'
__copyright__ = 'Copyright 2019'
__credits__ = ['Pierre Chaumeil']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Pierre Chaumeil'
__email__ = 'p.chaumeil@uq.edu.au'
__status__ = 'Development'

import os
import glob
import sys
import argparse
import operator
import logging
import re

import csv
csv.field_size_limit(sys.maxsize)

import pandas as pd
import numpy as np

from biolib.common import check_file_exists, make_sure_path_exists, check_dir_exists

from gtdb_migration_tk.prettytable import PrettyTable


class Tools(object):
    def __init__(self):
        """Initialization."""
        self.logger = logging.getLogger()

    def compare_metadata(self, old_meta_file, new_meta_file):

        old_delimiter = None
        # Parse TSV or CSV file
        for line in open(old_meta_file):
            if len(line.split('\t')) >= len(line.split(',')):
                old_delimiter = '\t'
                break
            else:
                old_delimiter = ','
                break

        old_nested_dict = {}
        with open(old_meta_file, 'r') as omf:
            old_headers = omf.readline().split(old_delimiter)
            if old_delimiter == ',':
                for line in csv.reader(omf):
                    old_nested_dict[line[0]] = {}
                    for i, j in enumerate(line):
                        old_nested_dict[line[0]][old_headers[i]] = str(j)
            else:
                for raw_line in omf:
                    line = raw_line.strip('\n').split('\t')
                    old_nested_dict[line[0]] = {}
                    for i, j in enumerate(line):
                        old_nested_dict[line[0]][old_headers[i]] = str(j)

        self.logger.info('{} parsed'.format(old_meta_file))

        new_delimiter = None
        # Parse TSV or CSV file
        for line in open(old_meta_file):
            if len(line.split('\t')) >= len(line.split(',')):
                new_delimiter = '\t'
                break
            else:
                new_delimiter = ','
                break

        header_summary = {}
        new_nested_dict = {}
        # in the new metadata file
        # we check if the genome id exists, and the columns names exist
        # for each common column name we compare the value for each common
        # genomes and add 1 if they are different
        with open(new_meta_file, 'r') as nmf:
            new_headers = nmf.readline().split(new_delimiter)
            if new_delimiter == ',':
                for line in csv.reader(nmf):
                    if line[0] in old_nested_dict:
                        for i, j in enumerate(line):
                            if new_headers[i] in old_headers:
                                if str(j) != old_nested_dict.get(line[0]).get(new_headers[i]):
                                    header_summary.setdefault(
                                        new_headers[i], []).append(1)
                                else:
                                    header_summary.setdefault(
                                        new_headers[i], []).append(0)

        for k, v in header_summary.items():
            header_summary[k] = round(100 * float(sum(v)) / len(v), 2)

        sorted_d = sorted(header_summary.items(), key=operator.itemgetter(1))

        # We display all common columns from less changes to more changes
        prettyt = PrettyTable()
        prettyt.field_names = ["Column name", "Difference (%)", "Note"]
        for coln in sorted_d:
            prettyt.add_row([coln[0], coln[1], ''])

        print(prettyt)

        removed_columns = set(old_headers) - \
            set(new_headers)
        new_columns = set(new_headers) - set(old_headers)

        print("Deprecated columns:")
        for removed_column in removed_columns:
            print("\t- {}".format(removed_column))

        print("New columns:")
        for new_column in new_columns:
            print("\t- {}".format(new_column))
