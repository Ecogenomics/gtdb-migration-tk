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
import gzip
import os
import shutil
import sys
import logging
import multiprocessing as mp
import tempfile

from collections import defaultdict
from pathlib import Path

from gtdb_migration_tk.biolib_lite.common import make_sure_path_exists
from gtdb_migration_tk.biolib_lite.external.execute import check_dependencies
from gtdb_migration_tk.biolib_lite.checksum import sha256, sha256_rb
from tqdm import tqdm

from gtdb_migration_tk.biolib_lite.external.pfam_search import PfamSearch
from gtdb_migration_tk.utils.tools import symlink


class MarkerManager(object):
    """Identify marker genes using Pfam and tigrfam HMMs."""

    def __init__(self, tmp_dir='/tmp/', cpus=1):
        """Initialization."""

        self.tmp_dir = tmp_dir
        self.cpus = cpus

        check_dependencies(['prodigal', 'hmmsearch'])

        # identify TIGRfam and Pfam marker genes comprising the bac120, ar122, ar53, or
        # rp2 marker sets using a carefully selected subset of HMMs
        self.tigrfam_hmms = ''
        self.pfam_hmm_dir = ''

        self.protein_file_ext = '_protein.faa.gz'

        self.logger = logging.getLogger('timestamp')

    def run_hmmsearch(self, gtdb_genome_path_file, report, db, folder_suffix,hmm_db_path):
        """Identify marker genes using Pfam and TIGRfam HMMs."""

        name = ""
        worker = None
        if db == 'pfam':
            marker_folder = 'pfam_{}'.format(folder_suffix)
            full_extension = '_pfam_{}.tsv'.format(folder_suffix)
            name = 'Pfam'
            self.pfam_hmm_dir = hmm_db_path
            worker = self.__pfam_worker
        elif db == 'tigrfam':
            marker_folder = 'tigrfam_{}'.format(folder_suffix)
            full_extension = '_tigrfam_{}.tsv'.format(folder_suffix)
            name = 'Tigrfam'
            self.tigrfam_hmms = hmm_db_path
            worker = self.__tigrfam_worker
        full_gz_extension = full_extension + '.gz'

        # get genomes marker as new or modified, and limit
        # marker gene finding to this subset of genomes
        genomes_to_consider = set()
        for line in open(report):
            line_split = line.strip().split('\t')
            genome_id = line_split[1]
            attributes = line_split[2].split(';')

            for attribute in attributes:
                if attribute == 'new' or attribute == 'modified':
                    genomes_to_consider.add(genome_id)

        self.logger.info(f'Identified {len(genomes_to_consider)} genomes as new or modified.')

        # get path to all unprocessed genome gene files
        self.logger.info('Checking genomes.')
        genome_files = []

        with open(gtdb_genome_path_file,'r') as  ggpf:
            for idx,line in enumerate(tqdm(ggpf)):
                gid,gpath,*_ = line.strip().split('\t')

                prodigal_dir = os.path.join(gpath, 'prodigal')
                marker_file = os.path.join(prodigal_dir, marker_folder, gid + full_extension)
                marker_zipped_file = marker_file + '.gz'
                if os.path.exists(marker_zipped_file):
                    # verify checksum
                    checksum_file = marker_file + '.sha256'
                    if os.path.exists(checksum_file):
                        checksum = sha256_rb(gzip.GzipFile(fileobj=open(marker_zipped_file, 'rb')))
                        cur_checksum = open(checksum_file).readline().strip()
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
                    print('Already processed', marker_zipped_file)
                    self.logger.warning(
                        f'Genome {gid} has no {name} annotations, but is also not marked for processing?')
                    self.logger.warning(f'Genome will be reannotated!')

                gene_file = os.path.join(
                    prodigal_dir, gid + self.protein_file_ext)
                if os.path.exists(gene_file):
                    if os.stat(gene_file).st_size == 0:
                        self.logger.warning(f' Protein file appears to be empty: {gene_file}')
                    else:
                        genome_files.append(gene_file)

        self.logger.info(f'Number of unprocessed genomes: {len(genome_files)}')

        # identify marker genes in parallel using HMMs and the HMMER package
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()


        for f in genome_files:
            workerQueue.put(f)

        for _ in range(self.cpus):
            workerQueue.put(None)

        try:
            workerProc = [mp.Process(target=worker, args=(
                workerQueue, writerQueue, folder_suffix)) for _ in range(self.cpus)]
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

    def run_tophit(self, gtdb_genome_path_file, db, folder_name):

        extension = ""
        if db == 'pfam':
            marker_version = 'pfam_{}'.format(folder_name)
            extension = f'_{marker_version}.tsv'
            tophit_out = f'_{marker_version}_tophit.tsv'
        elif db == 'tigrfam':
            marker_version = 'tigrfam_{}'.format(folder_name)
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

            prodigal_dir = os.path.join(gpath, 'prodigal')

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
                    self._parse_top_hit(output_hit_file, tophit_file,db)

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

    def __pfam_worker(self, queue_in, queue_out, folder_name):
        """Process each data item in parallel."""

        prefix = "pfam"
        pfam_version = '{}_{}'.format(prefix,folder_name)
        pfam_extension = f'_{pfam_version}.tsv'
        pfam_extension_gz = f'_{pfam_version}.tsv.gz'
        pfam_tophit_extension = f'_{pfam_version}_tophit.tsv'
        pfam_tophit_extension_gz = f'_{pfam_version}_tophit.tsv.gz'

        if '_lite' in pfam_extension:
            symlink_pfam_extension_gz = f'_{prefix}_lite.tsv.gz'
            symlink_pfam_tophit_extension_gz = f'_{prefix}_lite_tophit.tsv.gz'
        else:
            symlink_pfam_extension_gz = f'_{prefix}.tsv.gz'
            symlink_pfam_tophit_extension_gz = f'_{prefix}_tophit.tsv.gz'

        while True:
            gene_file = queue_in.get(block=True, timeout=None)
            if gene_file == None:
                break

            assembly_dir, filename = os.path.split(gene_file)
            make_sure_path_exists(os.path.join(assembly_dir, pfam_version))

            output_hit_file = os.path.join(
                assembly_dir, pfam_version, filename.replace(self.protein_file_ext, pfam_extension))
            #because the gene file is a zipped file, we need to unzip it in a temporary directory
            temp_dir = tempfile.mkdtemp()
            try:
                temp_gene_file = os.path.join(temp_dir, filename[0:-3])
                with gzip.open(gene_file, 'rb') as f_in:
                    with open(temp_gene_file, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)

                pfam_search = PfamSearch(self.pfam_hmm_dir)
                pfam_search.run(temp_gene_file, output_hit_file)

                # determine top hits
                pfam_tophit_file = os.path.join(assembly_dir, pfam_version, filename.replace(
                    self.protein_file_ext, pfam_tophit_extension))
                self._parse_top_hit(output_hit_file, pfam_tophit_file,prefix)



                # calculate checksum
                checksum = sha256(output_hit_file)
                fout = open(output_hit_file + '.sha256', 'w')
                fout.write(checksum)
                fout.close()

                # archive the pfam file and the tophit file
                with open(output_hit_file, 'rb') as f_in, gzip.open(output_hit_file + '.gz', 'wb') as f_out:
                    f_out.writelines(f_in)
                os.remove(output_hit_file)
                with open(pfam_tophit_file, 'rb') as f_in, gzip.open(pfam_tophit_file + '.gz', 'wb') as f_out:
                    f_out.writelines(f_in)
                os.remove(pfam_tophit_file)


                # create symlink in prodigal_folder
                new_hit_link = os.path.join(assembly_dir, filename.replace(
                    self.protein_file_ext, symlink_pfam_extension_gz))
                new_tophit_link = os.path.join(assembly_dir, filename.replace(
                    self.protein_file_ext, symlink_pfam_tophit_extension_gz))

                output_hit_file_relative = os.path.join(
                    '.', pfam_version, filename.replace(self.protein_file_ext, pfam_extension_gz))
                pfam_tophit_file_relative = os.path.join('.', pfam_version, filename.replace(
                    self.protein_file_ext, pfam_tophit_extension_gz))

                symlink(output_hit_file_relative, new_hit_link, overwrite=True)
                symlink(pfam_tophit_file_relative, new_tophit_link, overwrite=True)


            #we can now delete the temporary directory
            finally:
                shutil.rmtree(temp_dir)

            # allow results to be processed or written to file
            queue_out.put(gene_file)

    def _parse_top_hit(self,input_file,tophit_file,hmmdb):
        """Identify top Pfam and TIGRfam hits."""

        tophits = defaultdict(dict)
        gene_id = None
        for line in open(input_file):
            if hmmdb == 'tigrfam':
                if line[0] == '#' or line[0] == '[':
                    continue
                line_split = line.split()
                gene_id = line_split[0]
                hmm_id = line_split[3]
                evalue = float(line_split[4])
                bitscore = float(line_split[5])

            elif hmmdb == 'pfam':
                if line[0] == '#' or not line.strip():
                    continue
                line_split = line.split()
                gene_id = line_split[0]
                hmm_id = line_split[5]
                evalue = float(line_split[12])
                bitscore = float(line_split[11])


        if gene_id is None:
            self.logger.warning(
                f' No gene id found in {input_file} for hmmdb {hmmdb}')
        elif gene_id in tophits:
            if hmm_id in tophits[gene_id]:
                if bitscore > tophits[gene_id][hmm_id][1]:
                    tophits[gene_id][hmm_id] = (evalue, bitscore)
            else:
                tophits[gene_id][hmm_id] = (evalue, bitscore)
        else:
            tophits[gene_id][hmm_id] = (evalue, bitscore)

        fout = open(tophit_file, 'w')
        fout.write('Gene Id\tTop hits (Family id,e-value,bitscore)\n')
        for gene_id, hits in tophits.items():
            hit_str = []
            for hmm_id, stats in hits.items():
                hit_str.append(hmm_id + ',' + ','.join(map(str, stats)))
            fout.write('%s\t%s\n' % (gene_id, ';'.join(hit_str)))
        fout.close()

        # calculate checksum
        checksum = sha256(tophit_file)
        fout = open(tophit_file + '.sha256', 'w')
        fout.write(checksum)
        fout.close()

    # def _pfam_top_hit(self, pfam_file, pfam_tophit_file):
    #     """Identify top Pfam hits."""
    #
    #     tophits = defaultdict(dict)
    #     for line in open(pfam_file):
    #         if line[0] == '#' or not line.strip():
    #             continue
    #
    #         line_split = line.split()
    #         gene_id = line_split[0]
    #         hmm_id = line_split[5]
    #         evalue = float(line_split[12])
    #         bitscore = float(line_split[11])
    #         if gene_id in tophits:
    #             if hmm_id in tophits[gene_id]:
    #                 if bitscore > tophits[gene_id][hmm_id][1]:
    #                     tophits[gene_id][hmm_id] = (evalue, bitscore)
    #             else:
    #                 tophits[gene_id][hmm_id] = (evalue, bitscore)
    #         else:
    #             tophits[gene_id][hmm_id] = (evalue, bitscore)
    #
    #     fout = open(pfam_tophit_file, 'w')
    #     fout.write('Gene Id\tTop hits (Family id,e-value,bitscore)\n')
    #     for gene_id, hits in tophits.items():
    #         hit_str = []
    #         for hmm_id, stats in hits.items():
    #             hit_str.append(hmm_id + ',' + ','.join(map(str, stats)))
    #         fout.write('%s\t%s\n' % (gene_id, ';'.join(hit_str)))
    #     fout.close()
    #
    #     # calculate checksum
    #     checksum = sha256(pfam_tophit_file)
    #     fout = open(pfam_tophit_file + '.sha256', 'w')
    #     fout.write(checksum)
    #     fout.close()

    # def _tigr_top_hit(self, tigrfam_file, tigrfam_tophit_file):
    #     """Identify top TIGRfam hits."""
    #
    #     tophits = defaultdict(dict)
    #     for line in open(tigrfam_file):
    #         if line[0] == '#' or line[0] == '[':
    #             continue
    #
    #         line_split = line.split()
    #         gene_id = line_split[0]
    #         hmm_id = line_split[3]
    #         evalue = float(line_split[4])
    #         bitscore = float(line_split[5])
    #         if gene_id in tophits:
    #             if bitscore > tophits[gene_id][2]:
    #                 tophits[gene_id] = (hmm_id, evalue, bitscore)
    #         else:
    #             tophits[gene_id] = (hmm_id, evalue, bitscore)
    #
    #     fout = open(tigrfam_tophit_file, 'w')
    #     fout.write('Gene Id\tTop hits (Family id,e-value,bitscore)\n')
    #     for gene_id, stats in tophits.items():
    #         hit_str = ','.join(map(str, stats))
    #         fout.write('%s\t%s\n' % (gene_id, hit_str))
    #     fout.close()
    #
    #     # calculate checksum
    #     checksum = sha256(tigrfam_tophit_file)
    #     fout = open(tigrfam_tophit_file + '.sha256', 'w')
    #     fout.write(checksum)
    #     fout.close()

    def __tigrfam_worker(self, queue_in, queue_out,folder_name):
        """Process each data item in parallel."""

        prefix = "tigrfam"
        tigrfam_version = f'{prefix}_{folder_name}'
        tigrfam_extension = f'_{tigrfam_version}.tsv'
        tigrfam_extension_gz = f'_{tigrfam_version}.tsv.gz'
        tigrfam_out = f'_{tigrfam_version}.out'
        tigrfam_out_gz = f'_{tigrfam_version}.out.gz'
        tigrfam_tophit_extension = f'_{tigrfam_version}_tophit.tsv'
        tigrfam_tophit_extension_gz = f'_{tigrfam_version}_tophit.tsv.gz'


        if '_lite' in tigrfam_extension:
            symlink_tigrfam_extension_gz = f'_{prefix}_lite.tsv.gz'
            symlink_tigrfam_tophit_extension_gz = f'_{prefix}_lite_tophit.tsv.gz'
            symlink_tigrfam_out_gz = f'_{prefix}_lite.out.gz'
        else:
            symlink_tigrfam_extension_gz = f'_{prefix}.tsv.gz'
            symlink_tigrfam_tophit_extension_gz = f'_{prefix}_tophit.tsv.gz'
            symlink_tigrfam_out_gz = f'_{prefix}.out.gz'

        while True:
            gene_file = queue_in.get(block=True, timeout=None)
            if gene_file == None:
                break

            assembly_dir, filename = os.path.split(gene_file)
            make_sure_path_exists(os.path.join(assembly_dir, tigrfam_version))

            output_hit_file = os.path.join(assembly_dir, tigrfam_version, filename.replace(
                self.protein_file_ext, tigrfam_extension))
            out_file = os.path.join(assembly_dir, tigrfam_version, filename.replace(
                self.protein_file_ext, tigrfam_out))

            #because the gene file is a zipped file, we need to unzip it in a temporary directory
            temp_dir = tempfile.mkdtemp()
            try:
                temp_gene_file = os.path.join(temp_dir, filename[0:-3])
                with gzip.open(gene_file, 'rb') as f_in:
                    with open(temp_gene_file, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)

                hmmsearch_out = os.path.join(assembly_dir, tigrfam_version, filename.replace(
                    self.protein_file_ext, f'_{tigrfam_version}.out'))
                cmd = 'hmmsearch -o {} --tblout {} --noali --notextw --cut_nc --cpu 1 {} {}'.format(
                    hmmsearch_out,
                    output_hit_file,
                    self.tigrfam_hmms,
                    temp_gene_file)
                os.system(cmd)


                # determine top hits
                tigrfam_tophit_file = os.path.join(assembly_dir, tigrfam_version, filename.replace(
                    self.protein_file_ext, tigrfam_tophit_extension))
                self._parse_top_hit(output_hit_file, tigrfam_tophit_file,prefix)

                # calculate checksum
                checksum = sha256(output_hit_file)
                fout = open(output_hit_file + '.sha256', 'w')
                fout.write(checksum)
                fout.close()

                # archive the pfam file and the tophit file
                with open(output_hit_file, 'rb') as f_in, gzip.open(output_hit_file + '.gz', 'wb') as f_out:
                    f_out.writelines(f_in)
                os.remove(output_hit_file)
                with open(out_file, 'rb') as f_in, gzip.open(out_file + '.gz', 'wb') as f_out:
                    f_out.writelines(f_in)
                os.remove(out_file)
                with open(tigrfam_tophit_file, 'rb') as f_in, gzip.open(tigrfam_tophit_file + '.gz', 'wb') as f_out:
                    f_out.writelines(f_in)
                os.remove(tigrfam_tophit_file)

                # create symlink in prodigal_folder
                new_hit_link = os.path.join(assembly_dir, filename.replace(
                    self.protein_file_ext, symlink_tigrfam_extension_gz))
                new_tophit_link = os.path.join(assembly_dir, filename.replace(
                    self.protein_file_ext, symlink_tigrfam_tophit_extension_gz))
                new_out_link = os.path.join(assembly_dir, filename.replace(
                    self.protein_file_ext, symlink_tigrfam_out_gz))

                # Symlink needs to be relative to avoid pointing to previous version of Tigrfam when we copy folder
                output_hit_file_relative = os.path.join('.', tigrfam_version, filename.replace(
                    self.protein_file_ext, tigrfam_extension_gz))
                tigrfam_tophit_file_relative = os.path.join('.', tigrfam_version, filename.replace(
                    self.protein_file_ext, tigrfam_tophit_extension_gz))
                out_file_relative = os.path.join('.', tigrfam_version, filename.replace(
                    self.protein_file_ext, tigrfam_out_gz))

                symlink(output_hit_file_relative, new_hit_link,True)
                symlink(tigrfam_tophit_file_relative, new_tophit_link,True)
                symlink(out_file_relative, new_out_link,True)

            #we can now delete the temporary directory
            finally:
                shutil.rmtree(temp_dir)

            # allow results to be processed or written to file
            queue_out.put(gene_file)
