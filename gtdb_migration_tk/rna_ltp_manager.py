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
import glob
import os
import shutil
import sys
import logging
import ntpath

import datetime

from gtdb_migration_tk.biolib_lite.parallel import Parallel
from gtdb_migration_tk.genometk_lite.rna import RNA
from gtdb_migration_tk.biolib_lite.common import make_sure_path_exists




class RnaLTPManager(object):
    """Identify, extract, and taxonomically classify 16S rRNA genes against the LTP DB."""

    def __init__(self,cpus,rna_ltp_version,rna_ssu_version,rna_path):
        self.cpus = cpus
        self.rna_ltp_version = rna_ltp_version
        self.rna_ssu_version = rna_ssu_version
        self.rna_path = rna_path
        self.logger = logging.getLogger('timestamp')

        #E-value threshold for defining valid hits.
        self.evalue = 1e-6


        # needed to pick up previously identified and extracted 16S rRNA genes
        self.silva_output_dir = 'rna_silva_{}'.format(self.rna_ssu_version)


        self.ltp_output_dir = 'rna_ltp_{}'.format(self.rna_ltp_version)

        root_path = os.path.join(rna_path,'ltp',str(self.rna_ltp_version))

        self.ltp_ssu_file = os.path.join(root_path,'ltp_{}.fna'.format(str(self.rna_ltp_version)))
        self.ltp_taxonomy_file = os.path.join(root_path,'ltp_{}_taxonomy.tsv'.format(str(self.rna_ltp_version)))

        for item in [self.ltp_ssu_file,self.ltp_taxonomy_file]:
            if not os.path.exists(item):
                print('{} does not exist'.format(item))
                sys.exit(-1)

    def _producer(self, input_files):
        """Process each genome."""

        genome_file, ssu_file = input_files

        full_genome_dir, _ = ntpath.split(genome_file)

        output_dir = os.path.join(full_genome_dir, self.ltp_output_dir)

        rna = RNA(self.cpus,'ssu')
        make_sure_path_exists(output_dir)
        # we clean the directory
        filelist = glob.glob(os.path.join(output_dir, "*"))
        for f in filelist:
            os.remove(f)

        rna.classify(ssu_file,self.ltp_ssu_file,self.ltp_taxonomy_file,self.evalue, output_dir)
        canary_file = os.path.join(full_genome_dir, self.ltp_output_dir, 'ltp.canary.txt')
        # print(canary_file)
        with open(canary_file, 'w') as filehandle:
            filehandle.write(f'Silva version:{self.rna_ssu_version}.\n')
            filehandle.write(f'LTP version:{self.rna_ltp_version}.\n')
            filehandle.write('done.\n')
        return output_dir

    def _progress(self, processed_items, total_items):
        current_time_utc = datetime.datetime.utcnow().replace(microsecond=0)
        if processed_items > 0:
            time_left= (current_time_utc - self.starttime) * (total_items-processed_items)/ processed_items
            return 'Processed {} of {} ({}%) genomes. (ETA {})           '.format(processed_items,
                                                               total_items,
                                                               round(processed_items * 100.0 / total_items,2),time_left)


    def generate_rna_ltp(self, gtdb_genome_path_file,all_genomes=False):
        """Create metadata by parsing assembly stats files."""

        input_files = []

        self.starttime = datetime.datetime.utcnow().replace(microsecond=0)
        input_files = []
        countr = 0
        for line in open(gtdb_genome_path_file):
            countr += 1
            status_str = '{} lines read.'.format(countr)
            sys.stdout.write('%s\r' % status_str)
            sys.stdout.flush()

            line_split = line.strip().split('\t')

            gid = line_split[0]
            gpath = line_split[1]
            assembly_id = os.path.basename(os.path.normpath(gpath))

            ssu_file = os.path.join(
                                    gpath, self.silva_output_dir, 'ssu.fna')
            if os.path.exists(ssu_file):
                canary_file = os.path.join(gpath, self.ltp_output_dir, 'ltp.canary.txt')
                if not all_genomes and os.path.exists(canary_file):
                    continue
                genome_file = os.path.join(gpath, assembly_id + '_genomic.fna')
                input_files.append((genome_file, ssu_file))

        self.logger.info(f'{len(input_files)} ssu files to analyse.')
        # process each genome
        print('Generating metadata for each genome:')
        parallel = Parallel(cpus=self.cpus)
        parallel.run(self._producer,
                     None,
                     input_files,
                     self._progress)
