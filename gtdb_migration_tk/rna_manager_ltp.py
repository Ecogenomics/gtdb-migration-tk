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
import logging
import ntpath
import datetime
from typing import Tuple

from gtdb_migration_tk.genometk_lite.rna import RNA
from gtdb_migration_tk.biolib_lite.parallel import Parallel
from gtdb_migration_tk.biolib_lite.common import make_sure_path_exists, remove_files_in_directory


class RnaManagerLTP(object):
    """Identify, extract, and taxonomically classify 16S rRNA genes against the LTP DB."""

    def __init__(self, ltp_version: str, silva_version: str, rna_path: str, cpus: int) -> None:
        """Initialization."""

        self.cpus = cpus
        self.rna_ltp_version = ltp_version
        self.rna_ssu_version = silva_version
        self.rna_path = rna_path
        self.logger = logging.getLogger('timestamp')

        # needed to pick up previously identified and extracted 16S rRNA genes
        self.silva_output_dir = 'rna_silva_{}'.format(self.rna_ssu_version)

        self.ltp_output_dir = 'rna_ltp_{}'.format(self.rna_ltp_version)

        ltp_root_path = os.path.join(rna_path, str(self.rna_ltp_version))

        self.ltp_ssu_file = os.path.join(
            ltp_root_path, 'ltp_{}.fna'.format(str(self.rna_ltp_version)))
        self.ltp_taxonomy_file = os.path.join(
            ltp_root_path, 'ltp_{}_taxonomy.tsv'.format(str(self.rna_ltp_version)))

        for item in [self.ltp_ssu_file, self.ltp_taxonomy_file]:
            if not os.path.exists(item):
                print('{} does not exist'.format(item))
                sys.exit(-1)

    def _producer(self, input_data: Tuple[str, str, str]) -> str:
        """Process each genome."""

        genome_file, ssu_file, domain = input_data

        full_genome_dir, _ = ntpath.split(genome_file)

        output_dir = os.path.join(full_genome_dir, self.ltp_output_dir)

        rna = RNA('ssu', domain, self.cpus)
        make_sure_path_exists(output_dir)

        rna.classify(ssu_file,
                     self.ltp_ssu_file,
                     self.ltp_taxonomy_file,
                     output_dir)

        canary_file = os.path.join(full_genome_dir, self.ltp_output_dir, 'ltp.canary.txt')
        with open(canary_file, 'w') as filehandle:
            filehandle.write(f'Silva version:{self.rna_ssu_version}.\n')
            filehandle.write(f'LTP version:{self.rna_ltp_version}.\n')
            filehandle.write('done.\n')

        return output_dir

    def _progress(self, processed_items: int, total_items: int) -> str:
        current_time_utc = datetime.datetime.utcnow().replace(microsecond=0)
        if processed_items > 0:
            time_left = (current_time_utc - self.starttime) * \
                (total_items-processed_items) / processed_items

            progress_str = ' - processed {} of {} ({}%) genomes (ETA {})'.format(
                processed_items,
                total_items,
                round(processed_items * 100.0 / total_items, 2), time_left)

            progress_str = progress_str.ljust(72)

            return progress_str

        return ''

    def generate_rna_ltp(self,
                         gtdb_genome_path_file: str,
                         gtdb_domain_file: str,
                         all_genomes: bool = False,
                         remove_prior_results: bool = False) -> None:
        """Create metadata by parsing assembly stats files."""

        # determine domain of each genome
        gid_to_domain = {}
        with open(gtdb_domain_file) as f:
            header = f.readline().strip().split('\t')

            gid_idx = header.index('Genome Id')
            domain_idx = header.index('Predicted domain')
            ncbi_taxonomy_idx = header.index('NCBI taxonomy')

            for line in f:
                tokens = line.strip().split('\t')

                gid = tokens[gid_idx]
                if gid.startswith('GB_') or gid.startswith('RS_'):
                    gid = gid[3:]

                domain = tokens[domain_idx]
                if domain == 'None':
                    # default to domain indicated in NCBI taxonomy
                    domain = tokens[ncbi_taxonomy_idx].split(';')[0]

                if domain == 'd__Bacteria':
                    gid_to_domain[gid] = 'bac'
                elif domain == 'd__Archaea':
                    gid_to_domain[gid] = 'ar'
                else:
                    self.logger.error(f'Unrecognized domain for {gid}: {domain}')
                    sys.exit(1)

        # determine genomes to process
        self.logger.info('Determining genomes to process:')
        input_data = []
        missing_domain = []
        self.starttime = datetime.datetime.utcnow().replace(microsecond=0)
        with open(gtdb_genome_path_file) as f:
            for idx, line in enumerate(f):
                status_str = ' - {} lines read'.format(idx+1)
                sys.stdout.write('%s\r' % status_str)
                sys.stdout.flush()

                line_split = line.strip().split('\t')

                gid = line_split[0]
                gpath = line_split[1]
                assembly_id = os.path.basename(os.path.normpath(gpath))

                ssu_file = os.path.join(gpath, self.silva_output_dir, 'ssu.fna')
                if os.path.exists(ssu_file):
                    cur_ltp_output_dir = os.path.join(gpath, self.ltp_output_dir)
                    canary_file = os.path.join(cur_ltp_output_dir, 'ltp.canary.txt')

                    if not all_genomes and os.path.exists(canary_file):
                        continue

                    if remove_prior_results:
                        # remove all prior results in LTP output directory
                        remove_files_in_directory(cur_ltp_output_dir)

                    if gid not in gid_to_domain:
                        missing_domain.append(gid)
                        domain = 'bac'
                    else:
                        domain = gid_to_domain[gid]

                    genome_file = os.path.join(gpath, assembly_id + '_genomic.fna.gz')
                    input_data.append((genome_file, ssu_file, domain))

        if len(missing_domain) > 0:
            self.logger.warning(
                f'Identified {len(missing_domain):} genomes without a specified domain; assumed to be bacterial.')
            self.logger.warning(f'  - {",".join(missing_domain)}')

        self.logger.info(f'Identified {len(input_data)} 16S rRNA files to analyse.')

        # process each genome
        print('Generating metadata for each genome:')
        parallel = Parallel(cpus=self.cpus)
        parallel.run(self._producer,
                     None,
                     input_data,
                     self._progress)
