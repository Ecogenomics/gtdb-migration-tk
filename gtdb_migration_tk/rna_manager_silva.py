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
import gzip
import shutil
import ntpath
import datetime
import tempfile
import subprocess
import re
from typing import Tuple

from gtdb_migration_tk.genometk_lite.rna import RNA
from gtdb_migration_tk.biolib_lite.external.execute import check_dependencies
from gtdb_migration_tk.biolib_lite.parallel import Parallel
from gtdb_migration_tk.biolib_lite.common import remove_files_in_directory
from gtdb_migration_tk.biolib_lite.seq_io import read_seq
from gtdb_migration_tk.biolib_lite.external.blast import get_blastn_version


class RnaManagerSILVA(object):
    """Identify, extract, and taxonomically classify rRNA genes in genomes."""

    def __init__(self, silva_version: str, rna_path: str, rna_gene: str, cpus: int) -> None:
        """Initialization."""

        self.logger = logging.getLogger('timestamp')

        check_dependencies(['nhmmer', 'blastn'])
        self.logger.info('Using nhmmer v{}.'.format(self.get_nhmmer_version()))
        self.logger.info('Using blastn v{}.'.format(get_blastn_version()))

        self.silva_version = silva_version
        self.rna_path = rna_path
        self.rna_gene = rna_gene
        self.cpus = cpus

        self.silva_output_dir = 'rna_silva_{}-testing'.format(self.silva_version)

        silva_root_path = os.path.join(rna_path, str(self.silva_version))
        self.silva_ssu_taxonomy_file = os.path.join(silva_root_path, 'silva_taxonomy.ssu.tsv')
        self.silva_ssu_ref_file = os.path.join(silva_root_path, f'SILVA_{silva_version}_SSURef_NR99_tax_silva.fasta')

        file_to_check = [self.silva_ssu_taxonomy_file, self.silva_ssu_ref_file]

        if self.rna_gene != 'ssu':
            self.silva_lsu_ref_file = os.path.join(silva_root_path, f'SILVA_{silva_version}_LSURef_tax_silva.fasta')
            self.silva_lsu_taxonomy_file = os.path.join(silva_root_path, 'silva_taxonomy.lsu.tsv')
            file_to_check.append(self.silva_lsu_ref_file)
            file_to_check.append(self.silva_lsu_taxonomy_file)

        for item in file_to_check:
            if not os.path.exists(item):
                self.logger.error(f'{item} does not exist')
                sys.exit(-1)

        if self.rna_gene == 'ssu':
            self.db = self.silva_ssu_ref_file
            self.taxonomy = self.silva_ssu_taxonomy_file
        elif self.rna_gene == 'lsu_23S':
            self.db = self.silva_lsu_ref_file
            self.taxonomy = self.silva_lsu_taxonomy_file
        elif self.rna_gene == 'lsu_5S':
            self.db = None
            self.taxonomy = None
            self.logger.info(
                'We currently do not curate against a 5S database, but do identify these sequences for quality assessment purposes.')

    def get_nhmmer_version(self):
        """Returns the version of nhmmer on the system path.

        Returns
        -------
        str
            The string containing the nhmmer version.
        """
        try:
            proc = subprocess.Popen(['nhmmer', '-h'],
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    encoding='utf-8')
            stdout, _stderr = proc.communicate()
            version = re.search(r'HMMER (\S+)', stdout)
            if version:
                return version.group(1)
            else:
                return 'unknown'
        except Exception as e:
            print(e)
            return 'unknown'

    def _producer(self, input_data: Tuple[str, str]) -> str:
        """Process each genome."""

        genome_file, domain = input_data

        full_genome_dir, _ = ntpath.split(genome_file)
        output_dir = os.path.join(full_genome_dir, self.silva_output_dir)
        self.rna(genome_file, domain, output_dir)

        return output_dir

    def _progress(self, processed_items: int, total_items: int) -> str:
        current_time = datetime.datetime.utcnow().replace(microsecond=0)
        if processed_items > 0:
            time_left = (current_time - self.starttime) * \
                (total_items - processed_items) / processed_items

            progress_str = ' - processed {} of {} ({:.2f}%) genomes (ETA {})'.format(
                processed_items,
                total_items,
                processed_items * 100.0 / total_items,
                time_left)

            progress_str = progress_str.ljust(72)

            return progress_str

        return ''

    def generate_rna_silva(self,
                           gtdb_genome_path_file: str,
                           gtdb_domain_file: str,
                           all_genomes: bool,
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
        self.starttime = datetime.datetime.utcnow().replace(microsecond=0)
        missing_domain = []
        with open(gtdb_genome_path_file) as f:
            for idx, line in enumerate(f):
                statusStr = ' - {} lines read'.format(idx+1)
                sys.stdout.write('%s\r' % statusStr)
                sys.stdout.flush()

                tokens = line.strip().split('\t')

                gid = tokens[0]
                gpath = tokens[1]
                assembly_id = os.path.basename(os.path.normpath(gpath))
                genome_file = os.path.join(gpath, assembly_id + '_genomic.fna.gz')

                cur_silva_output_dir = os.path.join(gpath, self.silva_output_dir)
                canary_file = os.path.join(cur_silva_output_dir, self.rna_gene + '.canary.txt')
                if os.path.exists(canary_file) and not all_genomes:
                    continue

                if remove_prior_results:
                    # remove all prior results in SILVA output directory
                    remove_files_in_directory(cur_silva_output_dir)

                if gid not in gid_to_domain:
                    missing_domain.append(gid)
                    domain = 'bac'
                else:
                    domain = gid_to_domain[gid]

                input_data.append((genome_file, domain))

        if len(missing_domain) > 0:
            self.logger.warning(
                f'Identified {len(missing_domain):} genomes without a specified domain; assumed to be bacterial.')
            self.logger.warning(f'  - {",".join(missing_domain)}')

        # process each genome
        print(f'Generating metadata for {len(input_data):} genomes:')
        parallel = Parallel(cpus=self.cpus)
        parallel.run(self._producer,
                     None,
                     input_data,
                     self._progress)

    def rna(self, genome_file: str, domain: str, output_dir: str) -> None:
        """Identify, extract, and taxonomically classify rRNA genes in genomes."""

        # get HMM directory and HMM models
        file_dir = os.path.dirname(os.path.realpath(__file__))
        hmm_dir = os.path.join(file_dir, 'data_files', 'barrnap')

        rna_models = {}
        rna_models['ssu'] = {'ar': 'ar_16S', 'bac': 'bac_16S', 'euk': 'euk_18S'}
        rna_models['lsu_23S'] = {'ar': 'ar_23S', 'bac': 'bac_23S', 'euk': 'euk_28S'}
        rna_models['lsu_5S'] = {'ar': 'ar_5S', 'bac': 'bac_5S', 'euk': 'euk_5S'}

        hmm_model_prefix = rna_models[self.rna_gene][domain]

        # run each of the rRNA models
        rna = RNA(self.rna_gene, domain, 1)

        # we move the input file to a temp file
        temp_dir = tempfile.mkdtemp(suffix=os.path.basename(genome_file))

        tmp_hmm_model_file = os.path.join(temp_dir, hmm_model_prefix + '.hmm')
        shutil.copy(os.path.join(hmm_dir, hmm_model_prefix + '.hmm'), tmp_hmm_model_file)

        temp_taxonomy = None
        if self.taxonomy is not None:
            temp_taxonomy = os.path.join(temp_dir, os.path.basename(self.taxonomy))
            shutil.copy(self.taxonomy, temp_taxonomy)

        # HACK: create temporary genome file with dummy sequences at start to
        # ensure nhmmer can recognize that file uses DNA. This hack is currently
        # necessary as the --dna flag does not appear to work and nhmmer can fail
        # with an unknown base type error if the first sequence does not contain
        # all four nucleotides. This was an issue as of nhmmer 3.4 and impacts
        # a number of genomes including GCF_903931875.1.
        temp_genome_file = os.path.join(temp_dir, os.path.basename(genome_file))
        fout = gzip.open(temp_genome_file, 'wt')
        fout.write('>dummy\n')
        fout.write('ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n')
        for seq_id, seq in read_seq(genome_file):
            fout.write(f'>{seq_id}\n{seq}\n')
        fout.close()

        temp_output_dir = os.path.join(temp_dir, os.path.basename(output_dir))
        os.mkdir(temp_output_dir)

        try:
            rna.run(temp_genome_file,
                    tmp_hmm_model_file,
                    self.db,
                    temp_taxonomy,
                    temp_output_dir)

            # we copy temp_output_dir to output_dir
            if os.path.exists(output_dir):
                # we copy the files from temp_output_dir to output_dir
                for file in os.listdir(temp_output_dir):
                    shutil.copy(os.path.join(temp_output_dir, file), os.path.join(output_dir, file))
            else:
                shutil.copytree(temp_output_dir, output_dir)

        finally:
            shutil.rmtree(temp_dir)

    def update_silva(self, ssu_ref_file: str, lsu_ref_file: str, output_dir: str) -> None:
        """Update SILVA reference files."""

        fout = open(os.path.join(output_dir, 'silva_taxonomy.ssu.tsv'), 'w')
        for line in open(ssu_ref_file):
            if line[0] == '>':
                line_split = line[1:].strip().split(' ', 1)

                seq_id = line_split[0]
                taxonomy = line_split[1]
                fout.write(f'{seq_id}\t{taxonomy}\n')
        fout.close()

        fout = open(os.path.join(output_dir, 'silva_taxonomy.lsu.tsv'), 'w')
        for line in open(lsu_ref_file):
            if line[0] == '>':
                line_split = line[1:].strip().split(' ', 1)

                seq_id = line_split[0]
                taxonomy = line_split[1]
                fout.write(f'{seq_id}\t{taxonomy}\n')
        fout.close()
