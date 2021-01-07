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
import shutil
import logging
from tqdm import tqdm

from gtdb_migration_tk.biolib_lite.common import make_sure_path_exists,canonical_gid
from gtdb_migration_tk.biolib_lite.taxonomy import Taxonomy


class DirectoryManager(object):
    """Create file indicating directory of each genome."""

    def __init__(self):
        self.logger = logging.getLogger('timestamp')

    def generate_genome_dir_file(self, database_dir, output_file):
        """Create file indicating directory of each genome."""

        fout = open(output_file, 'w')

        # initialising progress bar objects
        high_loop = tqdm(range(1))
        mid_loop = tqdm(range(1))
        #low_loop = tqdm(range(1))

        for code in [('GCA','Genbank'), ('GCF','Refseq')]:
            code_dir = os.path.join(database_dir, code[0])
            if not os.path.exists(code_dir) or not os.path.isdir(code_dir):
                print('we skip {}'.format(code[1]))
                continue

            high_loop.refresh()
            high_loop.reset(total=len(os.listdir(code_dir)))

            for first_three in os.listdir(code_dir):
                high_loop.update()  # update outer tqdm
                onethird_species_dir = os.path.join(code_dir, first_three)
                if os.path.isfile(onethird_species_dir):
                    continue
                mid_loop.refresh()  # force print final state
                mid_loop.reset(total=len(os.listdir(onethird_species_dir)))

                for ii,second_three in enumerate(os.listdir(onethird_species_dir)):
                    mid_loop.update()
                    twothird_species_dir = os.path.join(
                        onethird_species_dir, second_three)
                    if os.path.isfile(twothird_species_dir):
                        continue
                    for iii,third_three in enumerate(os.listdir(twothird_species_dir)):
                        threethird_species_dir = os.path.join(
                            twothird_species_dir, third_three)
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
                                full_path), canonical_gid(accession)))

        fout.close()

    def delete_empty_directory(self, genome_path):
        """
        Delete a specific path.

        @param genome_path: path to delete
        @return: True
        """
        if len(os.listdir(genome_path)) == 0:
            os.rmdir(genome_path)
            self.delete_empty_directory(os.path.dirname(genome_path))
        return True

    def clean_ftp(self, new_list_genomes, ftp_genome_dir_file, report_dir, taxonomy_file=None):
        """Clean the FTP directory (remove deprecated genomes not appearing in the FTP folder anymore).

        Parameters
        ----------
        new_list_genomes : str
            Files indicating the Gid present in the new release (comma separated).
        ftp_genome_dir_file : str
            Genome directory file for the FTP server..
        report_dir : str
            Output directory to list reports.
        taxonomy_file : str
            Standardised taxonomy file from NCBI.
        """

        list_of_files = new_list_genomes.split(',')
        genome_in_new_rel = []
        make_sure_path_exists(report_dir)
        for new_genome_file in list_of_files:
            with open(new_genome_file, 'r') as ngf:
                for line in ngf:
                    genome_in_new_rel.append(line.strip().split('\t')[0])

        # read taxonomy file
        taxonomy = {}
        if taxonomy_file is not None:
            taxonomy = Taxonomy().read(taxonomy_file)

        current_ftp_genomes = {}
        with open(ftp_genome_dir_file) as fgdf:
            for line in fgdf:
                infos = line.strip().split('\t')
                current_ftp_genomes[infos[0]] = infos[1]

        deleted_genomes = list(
            set(current_ftp_genomes.keys()) - set(genome_in_new_rel))
        added_genomes = list(set(genome_in_new_rel) -
                             set(current_ftp_genomes.keys()))

        deleted_genome_file = open(os.path.join(
            report_dir, 'deleted_genomes.tsv'), 'w')
        added_genome_file = open(os.path.join(
            report_dir, 'added_genomes.tsv'), 'w')

        self.logger.info('{} genomes have been deleted in the release'.format(
            len(deleted_genomes)))
        self.logger.info('{} genomes have been added in the release'.format(
            len(added_genomes)))

        for idx, deleted_genome in enumerate(deleted_genomes):
            print("{}/{} genomes deleted".format(idx,
                                                 len(deleted_genomes)), end="\r")
            deleted_genome_file.write('{}\n'.format(deleted_genome))
            shutil.rmtree(current_ftp_genomes.get(deleted_genome))
            self.delete_empty_directory(os.path.dirname(
                current_ftp_genomes.get(deleted_genome)))

        for added_genome in added_genomes:
            added_genome_file.write('{}\t{}\n'.format(
                added_genome, taxonomy.get(added_genome, ['N/A'] * 7)[6]))
