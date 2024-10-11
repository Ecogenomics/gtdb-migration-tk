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
import sys
import logging
import ntpath
import datetime
import tempfile

from collections import defaultdict
import multiprocessing as mp


from tqdm import tqdm

from gtdb_migration_tk.biolib_lite.parallel import Parallel
from gtdb_migration_tk.genometk_lite.rna import RNA


class RnaManager(object):
    """Identify, extract, and taxonomically classify rRNA genes in genomes."""

    def __init__(self, cpus, rna_version, rna_path, rna_gene):
        self.cpus = cpus
        self.rna_version = rna_version
        self.rna_path = rna_path
        self.logger = logging.getLogger('timestamp')

        self.rna_gene = rna_gene
        if rna_gene == 'lsu_5S':
            self.min_len = 80
        else:
            self.min_len = 200

        self.silva_output_dir = 'rna_silva_{}'.format(self.rna_version)

        if self.rna_version:

            root_path = os.path.join(rna_path, str(self.rna_version))

            self.silva_ssu_taxonomy_file = os.path.join(
                root_path, 'silva_taxonomy.ssu.tsv')
            self.silva_ssu_ref_file = os.path.join(root_path, f'SILVA_{rna_version}_SSURef_NR99_tax_silva.fasta')

            file_to_check = [self.silva_ssu_taxonomy_file, self.silva_ssu_ref_file]

            if rna_gene != 'ssu':
                self.silva_lsu_ref_file = os.path.join(root_path, f'SILVA_{rna_version}_LSURef_tax_silva.fasta')
                self.silva_lsu_taxonomy_file = os.path.join(
                    root_path, 'silva_taxonomy.lsu.tsv')
                file_to_check.append(self.silva_lsu_ref_file)
                file_to_check.append(self.silva_lsu_taxonomy_file)

            for item in file_to_check:
                if not os.path.exists(item):
                    print('{} does not exist'.format(item))
                    sys.exit(-1)

    def _producer(self, job):
        """Process each genome."""
        genome_file = job
        full_genome_dir, _ = ntpath.split(genome_file)
        output_dir = os.path.join(full_genome_dir, self.silva_output_dir)

        rna_runner = self.rna(genome_file, output_dir)
        return (genome_file, 'done')

    def _progress(self, processed_items, total_items):
        current_time_utc = datetime.datetime.utcnow().replace(microsecond=0)
        if processed_items > 0:
            time_left = (current_time_utc - self.starttime) * \
                (total_items - processed_items) / processed_items
            return 'Processed {} of {} ({}%) genomes. (ETA {})           '.format(processed_items,
                                                                                  total_items,
                                                                                  round(processed_items * 100.0 / total_items, 2), time_left)

    def generate_rna_silva(self, gtdb_genome_path_file,rerun=False):
        """Create metadata by parsing assembly stats files."""
        # Silva info
        if self.rna_gene == 'ssu':
            self.db = self.silva_ssu_ref_file
            self.taxonomy = self.silva_ssu_taxonomy_file
        elif self.rna_gene == 'lsu_23S':
            self.db = self.silva_lsu_ref_file
            self.taxonomy = self.silva_lsu_taxonomy_file
        elif self.rna_gene == 'lsu_5S':
            self.db = None
            self.taxonomy = None
            print('We currently do not curate against a 5S database, but do identify these sequences for quality assessment purposes.')
        self.output_dir = self.silva_output_dir

        input_files = []

        self.starttime = datetime.datetime.utcnow().replace(microsecond=0)
        input_files = []
        list_genomes_tuples = []
        for line in open(gtdb_genome_path_file):
            list_genomes_tuples.append((line.strip(),rerun))

        with mp.Pool(processes=self.cpus) as pool:
            genome_paths = list(tqdm(pool.imap_unordered(self.rna_parser, list_genomes_tuples),
                                     total=len(list_genomes_tuples), unit='genome'))

        input_files = [x for x in genome_paths if x != 'null']

        # process each genome
        print(f'Generating metadata for each genome: {len(input_files)} genomes')
        processed_items = []
        with mp.Pool(processes=self.cpus) as pool:
            processed_items = list(tqdm(pool.imap_unordered(self._producer, input_files),
                                     total=len(input_files), unit='genome',ncols=100, smoothing=50/len(input_files)))

        # parallel = Parallel(cpus=self.cpus)
        # parallel.run(self._producer,
        #              None,
        #              input_files,
        #              self._progress)

    def rna_parser(self, job):
        line,rerun, = job
        line_split = line.strip().split('\t')

        gid = line_split[0]
        gpath = line_split[1]
        assembly_id = os.path.basename(os.path.normpath(gpath))

        genome_file = os.path.join(gpath, assembly_id + '_genomic.fna.gz')
        if rerun:
            return genome_file

        canary_file = os.path.join(
            gpath, self.output_dir, self.rna_gene + '.canary.txt')
        if not os.path.exists(canary_file):
            return genome_file
        else:
            return 'null'


    def rna(self, genome_file, output_dir):
        #self.logger.info('Identifying, extracting, and classifying rRNA genes.')

        # sanity check length
        if self.rna_gene == 'lsu_5S' and self.min_len > 120:
            self.logger.error(
                'Minimum length was set to %d, but LSU 5S genes are ~120 bp.' % self.min_len)
            sys.exit(-1)

        # get HMM directory and HMM models
        file_dir = os.path.dirname(os.path.realpath(__file__))
        hmm_dir = os.path.join(file_dir, 'data_files', 'barrnap')

        rna_models = defaultdict(list)
        rna_models['ssu'] = ('ar_16S', 'bac_16S', 'euk_18S')
        rna_models['lsu_23S'] = ('ar_23S', 'bac_23S', 'euk_28S')
        rna_models['lsu_5S'] = ('ar_5S', 'bac_5S', 'euk_5S')

        ar_model, bac_model, euk_model = rna_models[self.rna_gene]

        # run each of the rRNA models
        rna = RNA(1, self.rna_gene, self.min_len)

        #we move the input file to a temp file
        temp_dir = tempfile.mkdtemp(suffix=os.path.basename(genome_file))
        temp_genome_file = os.path.join(temp_dir,os.path.basename(genome_file))
        shutil.copy(genome_file,temp_genome_file)

        temp_ar_model = os.path.join(temp_dir,os.path.basename(ar_model)+'.hmm')
        shutil.copy(os.path.join(hmm_dir, ar_model + '.hmm'),temp_ar_model)
        temp_bac_model = os.path.join(temp_dir,os.path.basename(bac_model)+'.hmm')
        shutil.copy(os.path.join(hmm_dir, bac_model + '.hmm'),temp_bac_model)
        temp_euk_model = os.path.join(temp_dir,os.path.basename(euk_model)+'.hmm')
        shutil.copy(os.path.join(hmm_dir, euk_model + '.hmm'),temp_euk_model)



        temp_taxonomy = None
        if self.taxonomy is not None:
            temp_taxonomy = os.path.join(temp_dir,os.path.basename(self.taxonomy))
            shutil.copy(self.taxonomy,temp_taxonomy)


        temp_output_dir = os.path.join(temp_dir,os.path.basename(output_dir))
        os.mkdir(temp_output_dir)

        try:
            rna.run(temp_genome_file,
                    temp_ar_model,
                    temp_bac_model,
                    temp_euk_model,
                    self.db,
                    temp_taxonomy,
                    temp_output_dir)

            #we copy temp_output_dir to output_dir
            if os.path.exists(output_dir):
                # we copy the files from temp_output_dir to output_dir
                for file in os.listdir(temp_output_dir):
                    shutil.copy(os.path.join(temp_output_dir,file),os.path.join(output_dir,file))
            else:
                shutil.copytree(temp_output_dir,output_dir)

        finally:
            shutil.rmtree(temp_dir)

    def update_silva(self,ssu_ref_file,lsu_ref_file,output_dir):
        fout = open(os.path.join(output_dir,'silva_taxonomy.ssu.tsv'), 'w')
        for line in open(ssu_ref_file):
            if line[0] == '>':
                line_split = line[1:].strip().split(' ', 1)

                seq_id = line_split[0]
                taxonomy = line_split[1]
                fout.write(f'{seq_id}\t{taxonomy}\n')
        fout.close()

        fout = open(os.path.join(output_dir,'silva_taxonomy.lsu.tsv'), 'w')
        for line in open(lsu_ref_file):
            if line[0] == '>':
                line_split = line[1:].strip().split(' ', 1)

                seq_id = line_split[0]
                taxonomy = line_split[1]
                fout.write(f'{seq_id}\t{taxonomy}\n')
        fout.close()
