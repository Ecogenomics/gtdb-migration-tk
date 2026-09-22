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
import gzip
import logging
import multiprocessing as mp
import shutil
import tempfile
from collections import defaultdict
from multiprocessing.queues import Queue
from typing import Dict, Optional, Set, Tuple

from tqdm import tqdm

from gtdb_migration_tk.biolib_lite.checksum import sha256, sha256_rb
from gtdb_migration_tk.biolib_lite.common import make_sure_path_exists
from gtdb_migration_tk.biolib_lite.external.execute import check_dependencies
from gtdb_migration_tk.biolib_lite.external.pfam_search import PfamSearch
from gtdb_migration_tk.update_genomes import genomes_to_regenerate
from gtdb_migration_tk.utils.tools import symlink, openfile


# One genome as marker_parser() receives it, the tuple being what mp.Pool.imap_unordered
# can carry: accession, its genome directory, the marker directory within prodigal/, the
# extension the marker table carries, the genomes the release says need annotating, and
# the database's name for the log.
MarkerJob = Tuple[str, str, str, str, Set[str], str]

# Gene ID to the hits kept for it: HMM ID to (e-value, bitscore) for Pfam, where a gene
# keeps its best hit per family, and a single (HMM ID, e-value, bitscore) for TIGRFAM,
# where it keeps one hit overall.
PfamTopHits = Dict[str, Dict[str, Tuple[float, float]]]
TigrTopHits = Dict[str, Tuple[str, float, float]]


class MarkerManager(object):
    """Identify marker genes using Pfam and tigrfam HMMs."""

    def __init__(self, tmp_dir: str = '/tmp/', cpus: int = 1) -> None:
        """Initialization.

        Parameters
        ----------
        tmp_dir : str
            Directory for scratch files; no results are written here.
        cpus : int
            How many genomes are annotated at once.
        """

        self.tmp_dir: str = tmp_dir
        self.cpus: int = cpus

        check_dependencies(['prodigal', 'hmmsearch'])

        # identify TIGRfam and Pfam marker genes comprising the bac120, ar122, ar53, or
        # rp2 marker sets using a carefully selected subset of HMMs. Which of the two is
        # searched is decided per run, so only one of these is ever set.
        self.tigrfam_hmms: str = ''
        self.pfam_hmm_dir: str = ''

        self.protein_file_ext: str = '_protein.faa.gz'

        self.logger: logging.Logger = logging.getLogger('timestamp')

    def run_hmmsearch(self,
                      gtdb_genome_path_file: str,
                      report: str,
                      db: str,
                      dir_suffix: str,
                      hmm_db_path: str) -> None:
        """Identify marker genes using Pfam and TIGRfam HMMs.

        Parameters
        ----------
        gtdb_genome_path_file : str
            genome_dirs file of the release: accession, path, canonical accession.
        report : str
            report.log of the release, read for the genomes needing annotation.
        db : str
            'pfam' or 'tigrfam'.
        dir_suffix : str
            Suffix of the marker directory and files, e.g. 33.1_lite.
        hmm_db_path : str
            The HMMs to search against.

        @return: nothing; results are written into each genome's prodigal/ directory.
        """

        name = ""
        worker = None
        if db == 'pfam':
            marker_dir = 'pfam_{}'.format(dir_suffix)
            full_extension = '_pfam_{}.tsv'.format(dir_suffix)
            name = 'Pfam'
            self.pfam_hmm_dir = hmm_db_path
            worker = self.__pfam_worker
        elif db == 'tigrfam':
            marker_dir = 'tigrfam_{}'.format(dir_suffix)
            full_extension = '_tigrfam_{}.tsv'.format(dir_suffix)
            name = 'Tigrfam'
            self.tigrfam_hmms = hmm_db_path
            worker = self.__tigrfam_worker
        full_gz_extension = full_extension + '.gz'

        # limit marker gene finding to the genomes the release did not bring
        # their derived data with; update_genomes owns what its report means
        genomes_to_consider = genomes_to_regenerate(report)

        self.logger.info(
            f'Identified {len(genomes_to_consider)} genomes whose markers must be searched again.')

        # get path to all unprocessed genome gene files
        self.logger.info('Checking genomes.')
        genome_files = []

        list_genomes_tuples = []
        with open(gtdb_genome_path_file,'r') as  ggpf:
            for idx,line in enumerate(tqdm(ggpf)):
                gid,gpath,*_ = line.strip().split('\t')
                list_genomes_tuples.append((gid,gpath,marker_dir,full_extension,genomes_to_consider,name))

            with mp.Pool(processes=self.cpus) as pool:
                genome_paths = list(tqdm(pool.imap_unordered(self.marker_parser, list_genomes_tuples),
                                         total=len(list_genomes_tuples), unit='genome'))

            # a skipped genome is None, and every None has to go: the queue below ends
            # with one per worker as the signal to stop, and a worker cannot tell a
            # genome that was skipped from the end of the work
            genome_files = [x for x in genome_paths if x is not None]



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
                workerQueue, writerQueue, dir_suffix)) for _ in range(self.cpus)]
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

    def marker_parser(self, job: MarkerJob) -> Optional[str]:
        """Decide whether one genome's markers still have to be searched for.

        Parameters
        ----------
        job : MarkerJob
            One genome, as run_hmmsearch() packed it.

        @return: the protein file to search, or None to skip the genome -- it is
                 already annotated, or it has no protein file to search. One value
                 for both, because run_hmmsearch() does the same thing with them
                 and a second one has only ever been a way of missing one of them.
        """

        gid, gpath,marker_dir,full_extension,genomes_to_consider,name = job
        prodigal_dir = os.path.join(gpath, 'prodigal')
        marker_file = os.path.join(prodigal_dir, marker_dir, gid + full_extension)
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
                    return None

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
                return gene_file

    def run_tophit(self, gtdb_genome_path_file: str, db: str, folder_name: str) -> None:
        """Reduce each genome's marker table to its top hits.

        Parameters
        ----------
        gtdb_genome_path_file : str
            genome_dirs file of the release: accession, path, canonical accession.
        db : str
            'pfam' or 'tigrfam'.
        folder_name : str
            Suffix of the marker directory and files, e.g. 33.1_lite.

        @return: nothing; a tophit file is written beside each marker table.
        """

        extension = ""
        if db == 'pfam':
            marker_version = 'pfam_{}'.format(folder_name)
            extension = f'_{marker_version}.tsv.gz'
            tophit_out = f'_{marker_version}_tophit.tsv'
        elif db == 'tigrfam':
            marker_version = 'tigrfam_{}'.format(folder_name)
            extension = f'_{marker_version}.tsv.gz'
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
                    if not os.path.exists(output_hit_file):
                        self.logger.warning(
                            f'Output file does not exist: {output_hit_file}')
                        continue
                    if db == 'pfam':
                        self._pfam_top_hit(output_hit_file, tophit_file)
                    elif db == 'tigrfam':
                        self._tigr_top_hit(output_hit_file, tophit_file)

                    with open(tophit_file, 'rb') as f_in, gzip.open(tophit_file + '.gz', 'wb') as f_out:
                        f_out.writelines(f_in)
                    os.remove(tophit_file)

    def __progress(self, num_items: int, queue_out: Queue) -> None:
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

    def __pfam_worker(self, queue_in: Queue, queue_out: Queue, folder_name: str) -> None:
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
                print(temp_gene_file)

                # if size of temp_gene_file is 0, then skip hmmsearch
                if os.stat(gene_file).st_size == 0:
                    self.logger.warning('Skipping %s because it is empty' % temp_gene_file)
                    continue

                with gzip.open(gene_file, 'rb') as f_in:
                    with open(temp_gene_file, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)

                pfam_search = PfamSearch(self.pfam_hmm_dir)
                pfam_search.run(temp_gene_file, output_hit_file)

                # determine top hits
                pfam_tophit_file = os.path.join(assembly_dir, pfam_version, filename.replace(
                    self.protein_file_ext, pfam_tophit_extension))
                self._pfam_top_hit(output_hit_file, pfam_tophit_file)



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

    # def _parse_top_hit(self,input_file,tophit_file,hmmdb):
    #     """Identify top Pfam and TIGRfam hits."""
    #
    #     tophits = defaultdict(dict)
    #
    #
    #     for line in openfile(input_file):
    #         gene_id = None
    #         if hmmdb == 'tigrfam':
    #             if line[0] == '#' or line[0] == '[':
    #                 continue
    #             line_split = line.split()
    #             gene_id = line_split[0]
    #             hmm_id = line_split[3]
    #             evalue = float(line_split[4])
    #             bitscore = float(line_split[5])
    #
    #         elif hmmdb == 'pfam':
    #             if line[0] == '#' or not line.strip():
    #                 continue
    #             line_split = line.split()
    #             gene_id = line_split[0]
    #             hmm_id = line_split[5]
    #             evalue = float(line_split[12])
    #             bitscore = float(line_split[11])
    #
    #
    #         if gene_id is None:
    #             self.logger.warning(
    #                 f' No gene id found in {input_file} for hmmdb {hmmdb}')
    #         elif gene_id in tophits:
    #             if hmm_id in tophits[gene_id]:
    #                 if bitscore > tophits[gene_id][hmm_id][1]:
    #                     tophits[gene_id][hmm_id] = (evalue, bitscore)
    #             else:
    #                 tophits[gene_id][hmm_id] = (evalue, bitscore)
    #         else:
    #             tophits[gene_id][hmm_id] = (evalue, bitscore)
    #
    #     fout = open(tophit_file, 'w')
    #     fout.write('Gene Id\tTop hits (Family id,e-value,bitscore)\n')
    #     for gene_id, hits in tophits.items():
    #         hit_str = []
    #         for hmm_id, stats in hits.items():
    #             hit_str.append(hmm_id + ',' + ','.join(map(str, stats)))
    #         fout.write('%s\t%s\n' % (gene_id, ';'.join(hit_str)))
    #     fout.close()
    #
    #     # calculate checksum
    #     checksum = sha256(tophit_file)
    #     fout = open(tophit_file + '.sha256', 'w')
    #     fout.write(checksum)
    #     fout.close()

    def _pfam_top_hit(self, pfam_file: str, pfam_tophit_file: str) -> None:
        """Identify top Pfam hits.

        Parameters
        ----------
        pfam_file : str
            Marker table written by the Pfam search.
        pfam_tophit_file : str
            Where the top hits are written, with a .sha256 beside them.

        @return: nothing; a gene keeps its best hit for each family it matched.
        """

        tophits: PfamTopHits = defaultdict(dict)
        for line in openfile(pfam_file):
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

    def _tigr_top_hit(self, tigrfam_file: str, tigrfam_tophit_file: str) -> None:
        """Identify top TIGRfam hits.

        Parameters
        ----------
        tigrfam_file : str
            Marker table written by the TIGRFAM search.
        tigrfam_tophit_file : str
            Where the top hits are written, with a .sha256 beside them.

        @return: nothing; a gene keeps one hit, the highest scoring of any family.
        """

        tophits: TigrTopHits = {}
        for line in openfile(tigrfam_file):
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

    def __tigrfam_worker(self, queue_in: Queue, queue_out: Queue, folder_name: str) -> None:
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

                # if size of temp_gene_file is 0, then skip hmmsearch
                if os.stat(temp_gene_file).st_size == 0:
                    self.logger.warning('Skipping %s because it is empty' % temp_gene_file)
                    continue

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
                self._tigr_top_hit(output_hit_file, tigrfam_tophit_file)

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
