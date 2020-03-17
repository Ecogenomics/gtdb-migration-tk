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

    def select_delimiter(self, metafile):
        # Parse TSV or CSV file
        for line in open(metafile):
            if len(line.split('\t')) >= len(line.split(',')):
                return '\t'
            else:
                return ','

    def file_len(self, fname):
        with open(fname) as f:
            for i, l in enumerate(f):
                pass
        return i + 1

    def compare_metadata(self, old_meta_file, new_meta_file, only_ncbi=False):

        old_delimiter = self.select_delimiter(old_meta_file)

        old_nested_dict = {}
        with open(old_meta_file, 'r') as omf:
            old_headers = omf.readline().split(old_delimiter)
            if old_delimiter == ',':
                for line in csv.reader(omf):
                    if (only_ncbi and not line[0].startswith('U_')) or not only_ncbi:
                        old_nested_dict[line[0]] = {}
                        for i, j in enumerate(line):
                            old_nested_dict[line[0]][old_headers[i]] = str(j)
            else:
                for raw_line in omf:
                    line = raw_line.strip('\n').split('\t')
                    if (only_ncbi and not line[0].startswith('U_')) or not only_ncbi:
                        old_nested_dict[line[0]] = {}
                        for i, j in enumerate(line):
                            old_nested_dict[line[0]][old_headers[i]] = str(j)

        self.logger.info('{} parsed'.format(old_meta_file))

        new_delimiter = self.select_delimiter(new_meta_file)

        header_summary = {}
        new_nested_dict = {}
        # in the new metadata file
        # we check if the genome id exists, and the columns names exist
        # for each common column name we compare the value for each common
        # genomes and add 1 if they are different
        number_of_genomes = 0
        with open(new_meta_file, 'r') as nmf:
            new_headers = nmf.readline().split(new_delimiter)
            if new_delimiter == ',':
                for line in csv.reader(nmf):

                    if line[0] in old_nested_dict:
                        number_of_genomes += 1
                        for i, j in enumerate(line):
                            if new_headers[i] in old_headers:
                                if str(j) != old_nested_dict.get(line[0]).get(new_headers[i]):
                                    header_summary.setdefault(
                                        new_headers[i], []).append(1)
                                else:
                                    header_summary.setdefault(
                                        new_headers[i], []).append(0)
                for raw_line in nmf:
                    line = raw_line.strip('\n').split('\t')
                    if line[0] in old_nested_dict:
                        number_of_genomes += 1
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

        print("Based on {} common genomes.".format(number_of_genomes))

        print("Deprecated columns:")
        for removed_column in removed_columns:
            print("\t- {}".format(removed_column))

        print("New columns:")
        for new_column in new_columns:
            print("\t- {}".format(new_column))

    def compare_selected_data(self, old_meta_file, new_meta_file, metafield, output_file, only_ncbi=False):
        old_delimiter = self.select_delimiter(old_meta_file)
        old_nested_dict = {}
        with open(old_meta_file, 'r') as omf:
            old_headers = omf.readline().split(old_delimiter)
            if metafield not in old_headers:
                self.logger.error(f'{metafield} is not in {old_meta_file}')
                sys.exit()

            if old_delimiter == ',':
                for line in csv.reader(omf):
                    if (only_ncbi and not line[0].startswith('U_')) or not only_ncbi:
                        old_nested_dict[line[0]] = str(
                            line[old_headers.index(metafield)])
            else:
                for raw_line in omf:
                    line = raw_line.strip('\n').split('\t')
                    if (only_ncbi and not line[0].startswith('U_')) or not only_ncbi:
                        old_nested_dict[line[0]] = str(
                            line[old_headers.index(metafield)])

        new_delimiter = self.select_delimiter(new_meta_file)
        new_nested_dict = {}
        with open(new_meta_file, 'r') as nmf:
            new_headers = nmf.readline().split(new_delimiter)
            if metafield not in new_headers:
                self.logger.error(f'{metafield} is not in {old_meta_file}')
                sys.exit()
            if new_delimiter == ',':
                for line in csv.reader(nmf):
                    if line[0] in old_nested_dict:
                        new_nested_dict[line[0]] = str(
                            line[new_headers.index(metafield)])
            else:
                for raw_line in nmf:
                    line = raw_line.strip('\n').split('\t')
                    if line[0] in old_nested_dict:
                        new_nested_dict[line[0]] = str(
                            line[new_headers.index(metafield)])

        results = []
        outf = open(output_file, 'w')
        outf.write('genome_id\told_value\tnew_value\tsimilarity\n')
        for k, v in new_nested_dict.items():
            similarity = 'Identical'
            if v != old_nested_dict.get(k):
                similarity = "Different"
            outf.write('{}\n'.format(
                '\t'.join([k, str(old_nested_dict.get(k)), str(v), similarity])))

        self.logger.info('{} parsed'.format(old_meta_file))
