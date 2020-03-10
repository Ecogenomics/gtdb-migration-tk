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


class DirectoryManager(object):
    """Create file indicating directory of each genome."""

    def __init__(self):
        pass

    def normaliseID(self, accession):
        normaccession = "G" + accession[4:accession.find('.', 0)]
        return normaccession

    def run(self, database_dir, output_file):
        """Create file indicating directory of each genome."""

        fout = open(output_file, 'w')

        for code in ['GCA', 'GCF']:
            code_dir = os.path.join(database_dir, code)
            if not os.path.exists(code_dir) or not os.path.isdir(code_dir):
                print('we skip {}'.format(code))
                continue
            for first_three in os.listdir(code_dir):
                onethird_species_dir = os.path.join(code_dir, first_three)
                # print onethird_species_dir
                if os.path.isfile(onethird_species_dir):
                    continue
                for second_three in os.listdir(onethird_species_dir):
                    twothird_species_dir = os.path.join(
                        onethird_species_dir, second_three)
                    # print twothird_species_dir
                    if os.path.isfile(twothird_species_dir):
                        continue
                    for third_three in os.listdir(twothird_species_dir):
                        threethird_species_dir = os.path.join(
                            twothird_species_dir, third_three)
                        # print threethird_species_dir
                        if os.path.isfile(threethird_species_dir):
                            continue
                        for complete_name in os.listdir(threethird_species_dir):
                            full_path = os.path.join(
                                threethird_species_dir, complete_name)
                            if os.path.isfile(full_path):
                                continue

                            accession = complete_name[0:complete_name.find(
                                '_', 4)]

                            fout.write('{}\t{}\t{}\n'.format(accession, os.path.abspath(
                                full_path), self.normaliseID(accession)))

        fout.close()
