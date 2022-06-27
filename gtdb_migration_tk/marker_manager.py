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
import logging
import ntpath
import multiprocessing as mp

from collections import defaultdict

from biolib.common import remove_extension, make_sure_path_exists
from biolib.external.execute import check_dependencies
from biolib.checksum import sha256

from gtdb_migration_tk.prettytable import PrettyTable

from gtdb_migration_tk.tools import Tools


class MarkerManager(object):
    """Create file indicating directory of each genome."""

    def __init__(self, tmp_dir='/tmp/', cpus=1):
        """Initialization."""

        self.tmp_dir = tmp_dir
        self.cpus = cpus

        check_dependencies(['prodigal', 'hmmsearch', 'pfam_search.pl'])

        self.tigrfam_hmms = '/srv/whitlam/bio/db/tigrfam/15.0/TIGRFAMs_15.0_HMM/tigrfam.hmm'

        self.pfam_hmm_dir = '/srv/db/pfam/33.1/'
        self.protein_file_ext = '_protein.faa'

        self.logger = logging.getLogger('timestamp')

    def run_hmmsearch(self, gtdb_genome_path_file, report, db):
        extension = ""
        name = ""
        worker = None
        if db == 'pfam':
            marker_folder = 'pfam_33.1'
            full_extension = '_pfam_33.1.tsv'
            symlink_extension = '_pfam.tsv'
            name = 'Pfam'
            worker = self.__pfam_worker
        elif db == 'tigrfam':
            marker_folder = 'tigrfam_15.0'
            full_extension = '_tigrfam_15.0.tsv'
            symlink_extension = '_tigrfam.tsv'
            #extension = '_tigrfam_15.0.tsv'
            name = 'Tigrfam'
            worker = self.__tigrfam_worker
        genomes_to_consider = set()
        for line in open(report):
            line_split = line.strip().split('\t')
            genome_id = line_split[1]

            attributes = line_split[2].split(';')
            for attribute in attributes:
                if attribute == 'new' or attribute == 'modified':
                    genomes_to_consider.add(genome_id)

        self.logger.info(
            f'Identified {len(genomes_to_consider)} genomes as new or modified.')

        # get path to all unprocessed genome gene files
        self.logger.info('Checking genomes.')
        genome_files = []
        countr = 0
        for line in open(gtdb_genome_path_file):
            countr += 1
            statusStr = '{} lines read.'.format(countr)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

            line_split = line.strip().split('\t')

            gid = line_split[0]
            gpath = line_split[1]
            assembly_id = os.path.basename(os.path.normpath(gpath))

            prodigal_dir = os.path.join(gpath, 'prodigal')
            marker_file = os.path.join(
                prodigal_dir, marker_folder, gid + full_extension)
            if os.path.exists(marker_file):
                #print("File exists: {}".format(marker_file))
                # verify checksum
                checksum_file = marker_file + '.sha256'
                if os.path.exists(checksum_file):
                    checksum = sha256(marker_file)
                    cur_checksum = open(
                        checksum_file).readline().strip()
                    if checksum == cur_checksum:
                        if gid in genomes_to_consider:
                            self.logger.warning(
                                f'Genome {gid} is marked as new or modified, but already has {name} annotations.')
                            self.logger.warning('Genome is being skipped!')
                        continue

                self.logger.warning(
                    f'Genome {gid} has {name} annotations, but an invalid checksum and was not marked for reannotation.')
                self.logger.warning(f'Genome will be reannotated.')

            elif gid not in genomes_to_consider:
                self.logger.warning(
                    f'Genome {gid} has no {name} annotations, but is also not marked for processing?')
                self.logger.warning(f'Genome will be reannotated!')

            gene_file = os.path.join(
                prodigal_dir, gid + self.protein_file_ext)
            if os.path.exists(gene_file):
                if os.stat(gene_file).st_size == 0:
                    self.logger.warning(
                        f' Protein file appears to be empty: {gene_file}')
                else:
                    genome_files.append(gene_file)

        self.logger.info(f'Number of unprocessed genomes: {len(genome_files)}')

        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        for f in genome_files:
            workerQueue.put(f)

        for _ in range(self.cpus):
            workerQueue.put(None)

        try:
            workerProc = [mp.Process(target=worker, args=(
                workerQueue, writerQueue)) for _ in range(self.cpus)]
            writeProc = mp.Process(target=self.__progress, args=(
                len(genome_files), writerQueue))

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

            writeProc.terminate

    def run_tophit(self, gtdb_genome_path_file, db):

        extension = ""
        name = ""
        if db == 'pfam':
            marker_version = 'pfam_33.1'
            extension = f'_{marker_version}.tsv'
            tophit_out = f'_{marker_version}_tophit.tsv'
        elif db == 'tigrfam':
            marker_version = 'tigrfam_15.0'
            extension = f'_{marker_version}.tsv'
            tophit_out = f'_{marker_version}_tophit.tsv'

        countr = 0
        for line in open(gtdb_genome_path_file):
            countr += 1
            statusStr = '{} lines read.'.format(countr)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

            line_split = line.strip().split('\t')

            gid = line_split[0]
            gpath = line_split[1]
            assembly_id = os.path.basename(os.path.normpath(gpath))

            prodigal_dir = os.path.join(gpath, 'prodigal')
            output_file = os.path.join(
                prodigal_dir, marker_version, gid + extension)

            gene_file = os.path.join(
                prodigal_dir, gid + self.protein_file_ext)
            if os.path.exists(gene_file):
                if os.stat(gene_file).st_size == 0:
                    self.logger.warning(
                        f' Protein file appears to be empty: {gene_file}')
                else:
                    assembly_dir, filename = os.path.split(gene_file)

                    output_hit_file = os.path.join(assembly_dir, marker_version, filename.replace(
                        self.protein_file_ext, extension))
                    # determine top hits
                    tophit_file = os.path.join(assembly_dir, marker_version, filename.replace(
                        self.protein_file_ext, tophit_out))
                    if db == 'pfam':
                        self._pfam_top_hit(output_hit_file, tophit_file)
                    elif db == 'tigrfam':
                        self._tigr_top_hit(output_hit_file, tophit_file)

    def __progress(self, num_items, queue_out):
        """Store or write results of worker threads in a single thread."""
        processed_items = 0
        while True:
            a = queue_out.get(block=True, timeout=None)
            if a == None:
                break

            processed_items += 1
            statusStr = 'Finished processing %d of %d (%.2f%%) items.' % (
                processed_items, num_items, float(processed_items) * 100 / num_items)
            sys.stdout.write('%s\r' % statusStr)
        sys.stdout.flush()

        sys.stdout.write('\n')

    def __pfam_worker(self, queue_in, queue_out):
        """Process each data item in parallel."""

        pfam_version = 'pfam_33.1'
        pfam_extension = f'_{pfam_version}.tsv'
        pfam_tophit_extension = f'_{pfam_version}_tophit.tsv'

        symlink_pfam_extension = '_pfam.tsv'
        symlink_pfam_tophit_extension = '_pfam_tophit.tsv'

        while True:
            gene_file = queue_in.get(block=True, timeout=None)
            if gene_file == None:
                break

            assembly_dir, filename = os.path.split(gene_file)
            make_sure_path_exists(os.path.join(assembly_dir, pfam_version))

            output_hit_file = os.path.join(
                assembly_dir, pfam_version, filename.replace(self.protein_file_ext, pfam_extension))
            cmd = 'pfam_search.pl -outfile %s -cpu 1 -fasta %s -dir %s' % (
                output_hit_file, gene_file, self.pfam_hmm_dir)
            os.system(cmd)
            # print(cmd)

            # calculate checksum
            checksum = sha256(output_hit_file)
            fout = open(output_hit_file + '.sha256', 'w')
            fout.write(checksum)
            fout.close()

            # determine top hits
            pfam_tophit_file = os.path.join(assembly_dir, pfam_version, filename.replace(
                self.protein_file_ext, pfam_tophit_extension))
            self._pfam_top_hit(output_hit_file, pfam_tophit_file)

            # create symlink in prodigal_folder
            new_hit_link = os.path.join(assembly_dir, filename.replace(
                self.protein_file_ext, symlink_pfam_extension))
            new_tophit_link = os.path.join(assembly_dir, filename.replace(
                self.protein_file_ext, symlink_pfam_tophit_extension))

            # ==================================================================
            # print(f'{new_hit_link} will point to {output_hit_file}')
            # print(f'{new_tophit_link} will point to {pfam_tophit_file}')
            # ==================================================================

            os.symlink(output_hit_file, new_hit_link)
            os.symlink(pfam_tophit_file, new_tophit_link)

            # allow results to be processed or written to file
            queue_out.put(gene_file)

    def _pfam_top_hit(self, pfam_file, pfam_tophit_file):
        """Identify top Pfam hits."""

        tophits = defaultdict(dict)
        for line in open(pfam_file):
            if line[0] == '#' or not line.strip():
                continue

            line_split = line.split()
            gene_id = line_split[0]
            hmm_id = line_split[5]
            evalue = float(line_split[12])
            bitscore = float(line_split[11])
            if gene_id in tophits:
                if hmm_id in tophits[gene_id]:
                    if bitscore > tophits[gene_id][hmm_id][1]:
                        tophits[gene_id][hmm_id] = (evalue, bitscore)
                else:
                    tophits[gene_id][hmm_id] = (evalue, bitscore)
            else:
                tophits[gene_id][hmm_id] = (evalue, bitscore)

        fout = open(pfam_tophit_file, 'w')
        fout.write('Gene Id\tTop hits (Family id,e-value,bitscore)\n')
        for gene_id, hits in tophits.items():
            hit_str = []
            for hmm_id, stats in hits.items():
                hit_str.append(hmm_id + ',' + ','.join(map(str, stats)))
            fout.write('%s\t%s\n' % (gene_id, ';'.join(hit_str)))
        fout.close()

        # calculate checksum
        checksum = sha256(pfam_tophit_file)
        fout = open(pfam_tophit_file + '.sha256', 'w')
        fout.write(checksum)
        fout.close()

    def _tigr_top_hit(self, tigrfam_file, tigrfam_tophit_file):
        """Identify top TIGRfam hits."""

        tophits = {}
        for line in open(tigrfam_file):
            if line[0] == '#' or line[0] == '[':
                continue

            line_split = line.split()
            gene_id = line_split[0]
            hmm_id = line_split[3]
            evalue = float(line_split[4])
            bitscore = float(line_split[5])
            if gene_id in tophits:
                if bitscore > tophits[gene_id][2]:
                    tophits[gene_id] = (hmm_id, evalue, bitscore)
            else:
                tophits[gene_id] = (hmm_id, evalue, bitscore)

        fout = open(tigrfam_tophit_file, 'w')
        fout.write('Gene Id\tTop hits (Family id,e-value,bitscore)\n')
        for gene_id, stats in tophits.items():
            hit_str = ','.join(map(str, stats))
            fout.write('%s\t%s\n' % (gene_id, hit_str))
        fout.close()

        # calculate checksum
        checksum = sha256(tigrfam_tophit_file)
        fout = open(tigrfam_tophit_file + '.sha256', 'w')
        fout.write(checksum)
        fout.close()

    def __tigrfam_worker(self, queue_in, queue_out):
        """Process each data item in parallel."""
        tigrfam_version = 'tigrfam_15.0'
        tigrfam_extension = f'_{tigrfam_version}.tsv'
        tigrfam_tophit_extension = f'_{tigrfam_version}_tophit.tsv'

        symlink_tigrfam_extension = '_tigrfam.tsv'
        symlink_tigrfam_tophit_extension = '_tigrfam_tophit.tsv'

        while True:
            gene_file = queue_in.get(block=True, timeout=None)
            if gene_file == None:
                break

            assembly_dir, filename = os.path.split(gene_file)
            make_sure_path_exists(os.path.join(assembly_dir, tigrfam_version))

            output_hit_file = os.path.join(assembly_dir, tigrfam_version, filename.replace(
                self.protein_file_ext, tigrfam_extension))
            hmmsearch_out = os.path.join(assembly_dir, tigrfam_version, filename.replace(
                self.protein_file_ext, f'_{tigrfam_version}.out'))
            cmd = 'hmmsearch -o %s --tblout %s --noali --notextw --cut_nc --cpu 1 %s %s' % (
                hmmsearch_out, output_hit_file, self.tigrfam_hmms, gene_file)
            os.system(cmd)
            # ==================================================================
            # print(cmd)
            # ==================================================================

            # calculate checksum
            checksum = sha256(output_hit_file)
            fout = open(output_hit_file + '.sha256', 'w')
            fout.write(checksum)
            fout.close()

            # determine top hits
            tigrfam_tophit_file = os.path.join(assembly_dir, tigrfam_version, filename.replace(
                self.protein_file_ext, tigrfam_tophit_extension))
            self._tigr_top_hit(output_hit_file, tigrfam_tophit_file)

            # create symlink in prodigal_folder
            new_hit_link = os.path.join(assembly_dir, filename.replace(
                self.protein_file_ext, symlink_tigrfam_extension))
            new_tophit_link = os.path.join(assembly_dir, filename.replace(
                self.protein_file_ext, symlink_tigrfam_tophit_extension))

            # ==================================================================
            # print(f'{new_hit_link} will point to {output_hit_file}')
            # print(f'{new_tophit_link} will point to {tigrfam_tophit_file}')
            # ==================================================================

            os.symlink(output_hit_file, new_hit_link)
            os.symlink(tigrfam_tophit_file, new_tophit_link)

            # allow results to be processed or written to file
            queue_out.put(gene_file)
