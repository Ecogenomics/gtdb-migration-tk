import gzip
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
import ntpath
import argparse
import subprocess
import multiprocessing as mp
import logging
import tempfile

from biolib.common import remove_extension
from biolib.external.execute import check_dependencies
from biolib.checksum import sha256
from tqdm import tqdm


class tRNAScan(object):
    """Runs Run_tRNAScan-SE over a set of genomes."""

    def __init__(self, gbk_arc_assembly_file, gbk_bac_assembly_file, rfq_arc_assembly_file, rfq_bac_assembly_file,cpus):
        check_dependencies(['tRNAscan-SE'])

        self.genome_file_ext = '_genomic.fna.gz'

        self.logger = logging.getLogger('timestamp')


        self.cpus = cpus

        self.domain_dict = self.parseAssemblySummary(
            gbk_arc_assembly_file, gbk_bac_assembly_file, rfq_arc_assembly_file, rfq_bac_assembly_file)

    def parseAssemblySummary(self, gbk_arc_assembly_file, gbk_bac_assembly_file, rfq_arc_assembly_file, rfq_bac_assembly_file):
        results = {}
        for arcfile in (gbk_arc_assembly_file, rfq_arc_assembly_file):
            with open(arcfile) as arcf:
                arcf.readline()
                for line in arcf:
                    gid = line.split('\t')[0]
                    results[gid] = 'Archaea'
        for bacfile in (gbk_bac_assembly_file, rfq_bac_assembly_file):
            with open(bacfile) as bacf:
                bacf.readline()
                for line in bacf:
                    gid = line.split('\t')[0]
                    results[gid] = 'Bacteria'
        return results

    def __workerThread(self, queueIn, queueOut):
        """Process each data item in parallel."""
        while True:
            genome_file = queueIn.get(block=True, timeout=None)
            if genome_file == None:
                break

            assembly_dir, filename = os.path.split(genome_file)
            trna_dir = os.path.join(assembly_dir, 'trna')
            genome_id = filename[0:filename.find('_', 4)]

            if not os.path.exists(trna_dir):
                os.makedirs(trna_dir)

            output_file = os.path.join(trna_dir, genome_id + '_trna.tsv')
            log_file = os.path.join(trna_dir, genome_id + '_trna.log')
            stats_file = os.path.join(trna_dir, genome_id + '_trna_stats.tsv')


            #because the genome_file file is a zipped file, we need to unzip it in a temporary directory
            temp_dir = tempfile.mkdtemp()

            try:
                temp_gene_file = os.path.join(temp_dir, filename[0:-3])
                with gzip.open(genome_file, 'rb') as f_in:
                    with open(temp_gene_file, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)

                domain_flag = '-B'
                if self.domain_dict.get(genome_id) == 'Archaea':
                    domain_flag = '-A'

                #cmd = 'tRNAscan-SE %s -q -Q -o %s -m %s -l %s %s' % (domain_flag, output_file, stats_file, log_file, genome_file)
                # os.system(cmd)

                cmd_to_run = ['tRNAscan-SE', domain_flag, '-q', '-Q', '-o',
                              output_file, '-m', stats_file, '-l', log_file, temp_gene_file]
                proc = subprocess.Popen(
                    cmd_to_run, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdout, stderr = proc.communicate()
                # print proc.returncode
                if proc.returncode != 0:
                    raise RuntimeError("%r failed, status code %s stdout %r stderr %r" % (
                        cmd_to_run, proc.returncode, stdout, stderr))
                checksum_file = open(output_file + '.sha256', 'w')
                checksum_file.write('{}\n'.format(sha256(output_file)))
                checksum_file.close()
            # we can now delete the temporary directory
            finally:
                shutil.rmtree(temp_dir)

            queueOut.put(genome_file)

    def __writerThread(self, numDataItems, writerQueue):
        """Store or write results of worker threads in a single thread."""

        processedItems = 0
        while True:
            a = writerQueue.get(block=True, timeout=None)
            if a == None:
                break

            processedItems += 1
            statusStr = 'Finished processing %d of %d (%.2f%%) items.' % (processedItems,
                                                                          numDataItems,
                                                                          float(processedItems) * 100 / numDataItems)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

        sys.stdout.write('\n')

    def run(self, gtdb_genome_path_file,all_genomes=False):

        genomes_to_consider = None

        # get path to all genome files
        self.logger.info('Reading genomes.')
        genome_files = []
        list_genomes_to_parse = []
        for line in open(gtdb_genome_path_file):
            line_split = line.strip().split('\t')

            list_genomes_to_parse.append((line_split,all_genomes))

        print(f"number of cpus used:{self.cpus}")

        #populate worker queue with data to process
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()
        manager = mp.Manager()
        return_list = manager.list()


        for f in list_genomes_to_parse:
            workerQueue.put(f)

        for _ in range(self.cpus):
            workerQueue.put(None)

        try:
            workerProc = [mp.Process(target=self.trnascan_parser,
                                     args=(workerQueue, writerQueue,return_list))
                          for _ in range(self.cpus)]
            writeProc = mp.Process(target=self.__writerThread,
                                   args=(len(list_genomes_to_parse), writerQueue))

            writeProc.start()

            for p in workerProc:
                p.start()

            for p in workerProc:
                p.join()

            writerQueue.put(None)
            writeProc.join()

        except:
            for p in workerProc:
                p.terminate()

            writeProc.terminate()


        genome_files =[x for x in return_list if x != 'null']

        #genome_files = [x for x in writerQueue if x != 'null' and x is not None]

        self.logger.info(f' Number of unprocessed genomes: {len(genome_files)}')


        with mp.Pool(processes=self.cpus) as pool:
            genome_paths = list(tqdm(pool.imap_unordered(self.trnascan_worker, genome_files),
                                     total=len(genome_files), unit='genome',ncols=100))


            # workerProc = [mp.Process(target=self.__workerThread,
            #                          args=(workerQueue, writerQueue))
            #               for _ in range(self.cpus)]
            # writeProc = mp.Process(target=self.__writerThread,
            #                        args=(len(genome_files), writerQueue))
            #
            # writeProc.start()
            #
            # for p in workerProc:
            #     p.start()
            #
            # for p in workerProc:
            #     p.join()

        #     writerQueue.put(None)
        #     writeProc.join()
        # except:
        #     for p in workerProc:
        #         p.terminate()
        #
        #     writeProc.terminate()


    def trnascan_parser(self,queueIn, queueOut,return_list):
        while True:
            tuple_infos = queueIn.get(block=True, timeout=None)

            if tuple_infos == None:
                break
            value = None
            line_split, all_genomes = tuple_infos

            gid = line_split[0]
            gpath = line_split[1]
            assembly_id = os.path.basename(os.path.normpath(gpath))

            trna_dir = os.path.join(gpath, 'trna')

            trna_file = os.path.join(trna_dir, gid + '_trna.tsv')
            if all_genomes:
                genome_file = os.path.join(
                    gpath, assembly_id + self.genome_file_ext)
                value=genome_file

            if  value is None and os.path.exists(trna_file):
                # verify checksum
                checksum_file = trna_file + '.sha256'
                if os.path.exists(checksum_file):
                    checksum = sha256(trna_file)
                    cur_checksum = open(checksum_file).readline().strip()
                    if checksum == cur_checksum:
                        # if genomes_to_consider and gid in genomes_to_consider:
                        #     self.logger.warning(f'Genome {gid} is marked as new or modified, but already has tRNAs called.')
                        #     self.logger.warning('Genome is being skipped!')
                        value='null'

                if value is None:
                    self.logger.warning(
                            f'Genome {gid} has tRNAs called, but an invalid checksum ({checksum} for {trna_file} and {cur_checksum} in {checksum_file} and was not flagged.Genome will be reannotated.')

                    # elif genomes_to_consider and (gid not in genomes_to_consider):
                    #     self.logger.warning(f'Genome {gid} has no Pfam annotations, but is also not marked for processing?')
                    #     self.logger.warning('Genome will be reannotated!')

            genome_file = os.path.join(gpath, assembly_id + self.genome_file_ext)
            if os.path.exists(genome_file) and value != 'null':
                if os.stat(genome_file).st_size == 0:
                    self.logger.warning(f'Genome file appears to be empty: {gid}')
                else:
                    value=genome_file
            if value is not None:
                queueOut.put(value)
                return_list.append(value)
            else:
                queueOut.put(value)

    def trnascan_worker(self, job):
        """Process each data item in parallel."""

        genome_file = job

        assembly_dir, filename = os.path.split(genome_file)
        trna_dir = os.path.join(assembly_dir, 'trna')
        genome_id = filename[0:filename.find('_', 4)]

        if not os.path.exists(trna_dir):
            os.makedirs(trna_dir)

        output_file = os.path.join(trna_dir, genome_id + '_trna.tsv')
        log_file = os.path.join(trna_dir, genome_id + '_trna.log')
        stats_file = os.path.join(trna_dir, genome_id + '_trna_stats.tsv')


        #because the genome_file file is a zipped file, we need to unzip it in a temporary directory
        temp_dir = tempfile.mkdtemp()

        try:
            temp_gene_file = os.path.join(temp_dir, filename[0:-3])
            with gzip.open(genome_file, 'rb') as f_in:
                with open(temp_gene_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

            domain_flag = '-B'
            if self.domain_dict.get(genome_id) == 'Archaea':
                domain_flag = '-A'

            #cmd = 'tRNAscan-SE %s -q -Q -o %s -m %s -l %s %s' % (domain_flag, output_file, stats_file, log_file, genome_file)
            # os.system(cmd)

            cmd_to_run = ['tRNAscan-SE', domain_flag, '-q', '-Q', '-o',
                          output_file, '-m', stats_file, '-l', log_file, temp_gene_file]
            proc = subprocess.Popen(
                cmd_to_run, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = proc.communicate()
            # print proc.returncode
            if proc.returncode != 0:
                raise RuntimeError("%r failed, status code %s stdout %r stderr %r" % (
                    cmd_to_run, proc.returncode, stdout, stderr))
            checksum_file = open(output_file + '.sha256', 'w')
            checksum_file.write('{}\n'.format(sha256(output_file)))
            checksum_file.close()
        # we can now delete the temporary directory
        finally:
            shutil.rmtree(temp_dir)
