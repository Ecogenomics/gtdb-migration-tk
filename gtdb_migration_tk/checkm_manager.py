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
import tempfile
import ntpath
import shutil
import logging

from collections import defaultdict

from checkm.util.seqUtils import readFasta
from tqdm import tqdm

from gtdb_migration_tk.biolib_lite.common import get_num_lines


class CheckMManager(object):
    """Apply CheckM to a large set of genomes.

    This script assumes that genomes are stored in individual
    directories in the following format:

    <domain>/<genome_id>/<assembly_id>/<assembly_id>_protein.faa

    where <domain> is either 'archaea' or 'bacteria', and there
      may be multiple assembly_id for a given genome_id. These
      typically represent different strains from a species.

    This is the directory structure which results from extract_ncbi.py.
    """

    def __init__(self,cpus=1):
        """Initialization."""
        self.cpus = cpus
        self.logger = logging.getLogger('timestamp')


    def run_checkm(self, gtdb_genome_path_file, genome_report, output_dir, all_genomes=False):
        """Applying CheckM to genomes."""

        tmp_dir = os.path.join(output_dir, 'genome_chunks')

        if all_genomes:
            self.logger.info('Processing all genomes.')

        if not os.path.exists(tmp_dir):
            # get list of genomes to consider
            if genome_report.lower() != 'none':
                genomes_to_consider = set()
                for line in open(genome_report):
                    line_split = line.strip().split('\t')
                    genome_id = line_split[1]
                    attributes = line_split[2].split(';')
                    if 'removed' in attributes:
                        continue
                    else:
                        for attribute in attributes:
                            if all_genomes or attribute == 'new' or attribute == 'modified':
                                genomes_to_consider.add(genome_id)

                self.logger.info('Identified {} genomes as new or modified.'.format(len(genomes_to_consider)))

            # determine gene files
            gene_files = []

            genome_paths = {}
            for line in open(gtdb_genome_path_file):
                line_split = line.strip().split('\t')
                genome_paths[line_split[0]] = line_split[1]

            for gid in genomes_to_consider:
                gpath = genome_paths[gid]
                gene_file = os.path.join(gpath, 'prodigal', gid + '_protein.faa.gz')
                gene_files.append(gene_file)

            print('  Identified %d gene files.' % len(gene_files))

            # copy genomes in batches of 1000
            print('Partitioning genomes into chunks of 1000.')
            num_chunks = 0
            for i, gene_file in enumerate(gene_files):
                if i % 1000 == 0:
                    chunk_dir = os.path.join(tmp_dir, 'chunk%d' % num_chunks)
                    print(chunk_dir)
                    os.makedirs(chunk_dir)
                    num_chunks += 1

                #shutil.copy(gene_file, os.path.join(chunk_dir, ntpath.basename(gene_file)))
                os.system('ln -s %s %s' % (os.path.abspath(gene_file),
                                           os.path.join(chunk_dir, ntpath.basename(gene_file))))
        else:
            # just determine number of "chunk" directories
            self.logger.info(f'{output_dir}Output directory does exist.')
            num_chunks = 0
            for d in os.listdir(tmp_dir):
                if 'chunk' in d:
                    num_chunks += 1

        # apply CheckM to each set of 1000 genomes
        print('Running CheckM on chunks:')
        for i in range(0, num_chunks):
            print('  Processing chunk %d of %d.' % (i + 1, num_chunks))

            bin_dir = os.path.join(tmp_dir, 'chunk%d' % i)
            checkm_output_dir = os.path.join(output_dir, 'chunk%d' % i)
            if os.path.exists(checkm_output_dir):
                continue

            os.makedirs(checkm_output_dir)
            # check if all genomes in bin_dir are nucleotide or amino acid
            binFiles = self.binFiles(
                bin_dir,'gz')

            #self.checkProteinSeqs(binFiles)


            os.system('checkm lineage_wf --pplacer_threads %d --genes -x faa.gz -t %d %s %s' %
                      (self.cpus, self.cpus, bin_dir, checkm_output_dir))

            tree_qa_file = os.path.join(
                checkm_output_dir, 'tree_qa.o2.chunk%d.tsv' % i)
            os.system('checkm tree_qa -o 2 --tab_table -f %s %s' %
                      (tree_qa_file, checkm_output_dir))

            qa_file = os.path.join(checkm_output_dir, 'qa.chunk%d.tsv' % i)
            os.system('checkm qa -t %d --tab_table -f %s %s %s' % (self.cpus, qa_file,
                                                                   os.path.join(checkm_output_dir, 'lineage.ms'), checkm_output_dir))

            profile_file = os.path.join(
                checkm_output_dir, 'profile.chunk%d.tsv' % i)
            # os.system('checkm join_tables -f %s %s %s' %
            #           (profile_file,qa_file, tree_qa_file)
            self.joinTables([qa_file, tree_qa_file],profile_file)

            qa_file_sh100 = os.path.join(
                checkm_output_dir, 'qa_sh100.chunk%d.tsv' % i)
            alignment_file = os.path.join(
                checkm_output_dir, 'alignment_file.chunk%d.tsv' % i)
            os.system('checkm qa --aai_strain 0.9999 -t %d -a %s --tab_table -f %s %s %s' % (self.cpus,
                                                                                             alignment_file, qa_file_sh100, os.path.join(checkm_output_dir, 'lineage.ms'), checkm_output_dir))

        # create single file with CheckM results
        print('Creating single file with CheckM results.')
        checkm_output = os.path.join(output_dir, 'checkm.profiles.tsv')
        fout = open(checkm_output, 'w')
        for i in range(0, num_chunks):
            profile_file = os.path.join(
                output_dir, 'chunk%d' % i, 'profile.chunk%d.tsv' % i)
            with open(profile_file) as f:
                if i != 0:
                    f.readline()

                for line in f:
                    fout.write(line)
        fout.close()

        # create single file with CheckM strain heterogeneity results at 100%
        print('Creating single file with CheckM results.')
        checkm_output = os.path.join(output_dir, 'checkm.qa_sh100.tsv')
        fout = open(checkm_output, 'w')
        for i in range(0, num_chunks):
            qa_file = os.path.join(output_dir, 'chunk%d' %
                                   i, 'qa_sh100.chunk%d.tsv' % i)
            with open(qa_file) as f:
                if i != 0:
                    f.readline()

                for line in f:
                    fout.write(line)
        fout.close()

        # create single file with CheckM alignments for multi-copy genes
        print('Creating single file with CheckM results.')
        checkm_output = os.path.join(output_dir, 'checkm.alignment_file.tsv')
        fout = open(checkm_output, 'w')
        for i in range(0, num_chunks):
            align_file = os.path.join(
                output_dir, 'chunk%d' % i, 'alignment_file.chunk%d.tsv' % i)
            with open(align_file) as f:
                for line in f:
                    fout.write(line)
        fout.close()

        print('CheckM results written to: %s' % checkm_output)


    def checkProteinSeqs(self,seq_files):
        """Check if files contain sequences in amino acid space.

        Parameters
        ----------
        seq_files : iterable
            Sequence files to check.

        Returns
        -------
        boolean
            True if files can be treated as containing amino acid sequences.
        """

        for seq_file in seq_files:
            if os.stat(seq_file).st_size == 0:
                continue

            if self.isNucleotide(seq_file):
                logger = logging.getLogger('timestamp')
                logger.warning(
                    'File %s appears to contain nucleotide sequences. We delete the file from the input directory' % seq_file)
                os.remove(seq_file)

        return True

    def isNucleotide(self,seq_file, req_perc=0.9, max_seqs_to_read=10):
        """Check if a file contains sequences in nucleotide space.

        The check is performed by looking for the characters in
        {a,c,g,t,n,.,-} and confirming that these comprise the
        majority of a sequences. A set number of sequences are
        read and the file assumed to be not be in nucleotide space
        if none of these sequences are comprised primarily of the
        defined nucleotide set.

        Parameters
        ----------
        seq_file : str
            Name of fasta/q file to read.
        req_perc : float
            Percentage of bases in {a,c,g,t,n,.,-} before
            declaring the sequences as being in nucleotide
            space.
        max_seqs_to_read : int
            Maximum sequences to read before declaring
            sequence file to not be in nucleotide space.

        Returns
        -------
        boolean
            True is sequences are in nucleotide space, or file
            contains no sequences.

        """

        nucleotide_bases = {'a', 'c', 'g', 't'}
        insertion_bases = {'-', '.'}

        seqs = readFasta(seq_file)
        if len(seqs) == 0:
            return True

        seq_count = 0
        for _seq_id, seq in seqs.items():
            seq = seq.lower()

            nt_bases = 0
            for c in (nucleotide_bases | {'n'} | insertion_bases):
                nt_bases += seq.count(c)

            if float(nt_bases) / len(seq) >= req_perc:
                return True

            seq_count += 1
            if seq_count == max_seqs_to_read:
                break

        return False

    def binFiles(self, binInput, binExtension):
        binFiles = []
        binIDs = set()
        isInputDir = True
        if binInput is not None:
            if os.path.isdir(binInput):
                if binExtension[0] != '.':
                    binExtension = '.' + binExtension

                all_files = os.listdir(binInput)
                for f in all_files:
                    if f.endswith(binExtension):
                        binFile = os.path.join(binInput, f)
                        if os.stat(binFile).st_size == 0:
                            self.logger.warning(
                                "Skipping bin %s as it has a size of 0 bytes." % f)
                        else:
                            binFiles.append(binFile)
                            binIDs.add(os.path.basename(binFile))
            else:
                with open(binInput, "r") as oh:
                    for line in oh:
                        files = line.strip().split("\t")
                        binFile = files[1]
                        if not os.path.exists(binFile):
                            self.logger.warning(
                                "Skipping bin %s as it doesn't exists." % binFile)
                        elif os.stat(binFile).st_size == 0:
                            self.logger.warning(
                                "Skipping bin %s as it has a size of 0 bytes." % binFile)
                        else:
                            binFiles.append(binFile)
                            binIDs.add(os.path.basename(binFile))

        if not binFiles:
            if isInputDir:
                self.logger.error(
                    "No bins found. Check the extension (-x) used to identify bins.")
            else:
                self.logger.error(
                    "No bins found. Check the bins input table to verify bins exists.")
            sys.exit(1)

        if len(binIDs) != len(binFiles):
            self.logger.error(
                "There are redundant bin IDs, please check and update.")
            sys.exit(1)

        return sorted(binFiles)

    def joinTables(self, check_tables,check_file):
        self.logger.info('Joining tables containing bin information.')

        # read all tables
        headers = {}
        rows = defaultdict(dict)
        binIds = set()
        for f in check_tables:
            with open(f) as fin:
                headers[f] = [x.strip() for x in fin.readline().split('\t')][1:]

                for line in fin:
                    lineSplit = [x.strip() for x in line.split('\t')]

                    binId = lineSplit[0]
                    binIds.add(binId)

                    for i, header in enumerate(headers[f]):
                        rows[binId][header] = lineSplit[i + 1]

        # write merge table
        oldStdOut = self.reassignStdOut(check_file)

        row = 'Bin Id'
        for f in check_tables:
            row += '\t' + '\t'.join(headers[f])
        print(row)

        for binId in binIds:
            row = binId
            for f in check_tables:
                for header in headers[f]:
                    row += '\t' + rows[binId].get(header, '')
            print(row)

        self.restoreStdOut(check_file, oldStdOut)

        if check_file:
            self.logger.info('\n  Joined table written to: ' + check_file)

    def reassignStdOut(self,outFile):
        """Redirect standard out to a file."""
        oldStdOut = sys.stdout
        if(outFile != ''):
            try:
                # redirect stdout to a file
                sys.stdout = open(outFile, 'w')
            except:
                logger = logging.getLogger()
                logger.error("   [Error] Error diverting stdout to file: " + outFile)
                sys.exit(1)

        return oldStdOut


    def restoreStdOut(self,outFile, oldStdOut):
        """Redirect standard out back to system standard out."""
        if(outFile != ''):
            try:
                # redirect stdout to a file
                sys.stdout.close()
                sys.stdout = oldStdOut
            except:
                logger = logging.getLogger()
                logger.error("   [Error] Error restoring stdout ", outFile)
                sys.exit(1)

    def uniq(self,seq):
        seen = set()
        seen_add = seen.add

        return [x for x in seq if not (x in seen or seen_add(x))]

    def join_checkm_files_releases(self,releasefiles,output_file):
        outf = open(output_file,'w')
        lines = {}
        order_keys = []
        for idx,releasefile in enumerate(releasefiles):
            with open(releasefile,'r') as rf:
                if idx > 0 :
                    rf.readline()
                for line in rf:
                    infos = line.split("\t")
                    lines[infos[0]] = line
                    order_keys.append(infos[0])
        order_keys.reverse()
        reverse_list = self.uniq(order_keys)
        reverse_list.reverse()
        for key in reverse_list:
            outf.write(lines.get(key))
        outf.close()


    def prepare_checkm2_batch(self,checkm_summary_genbank, checkm_summary_refseq,
                              gtdb_metadata_file,gtdb_genome_path_file,output_dir):
        genome_paths = {}
        for line in open(gtdb_genome_path_file):
            line_split = line.strip().split('\t')
            genome_paths[line_split[0]] = line_split[1]

        # read genbank summary
        list_genomes_to_consider = []
        for file in (checkm_summary_genbank, checkm_summary_refseq):
            if file is not None:
                with open(file) as f:
                    headers = f.readline()
                    for line in f:
                        #TEMP use comma instead or tab
                        genome_id = line.strip().split('\t')[0]
                        # TEMP remove GB_GCA_ from genome_id and RS_GCF_ from genome_id
                        genome_id = genome_id.replace('GB_GCA_', 'GCA_').replace('RS_GCF_', 'GCF_')
                        genome_id = genome_id.split('_')[0]+'_'+genome_id.split('_')[1]
                        list_genomes_to_consider.append(genome_id)
        print(list_genomes_to_consider[:10])
        list_genomes_to_consider = set(list_genomes_to_consider)


        self.logger.info('Identified %d genomes to consider.' % len(list_genomes_to_consider))

        # get path to genomic FASTA files
        genomic_files = {}
        with open(gtdb_genome_path_file) as f:
            for line in tqdm(f,total=get_num_lines(gtdb_genome_path_file),ncols=100):
                tokens = line.strip().split('\t')

                gid = tokens[0]
                if gid in list_genomes_to_consider:
                    gp = tokens[1]
                    accn = tokens[1].split('/')[-1]
                    genomic_files[gid] = os.path.join(gp, f'{accn}_genomic.fna.gz')

        # determine coding table for each genome, seperating
        # out genomes using table 11 versus tables 4 or 25
        gid_tt11 = []
        gid_tt4 = []
        gid_tt25 = []
        gid_tt_none = []
        with open(gtdb_metadata_file, 'rt') as f:
            header = f.readline().strip().split('\t')

            trans_table_idx = header.index('ncbi_translation_table')

            for line in tqdm(f,total=get_num_lines(gtdb_metadata_file),ncols=100):
                tokens = [t.strip() for t in line.split('\t')]

                gid = tokens[0]
                if gid not in genomic_files:
                    continue
                trans_table = tokens[trans_table_idx]

                if trans_table == '11':
                    gid_tt11.append(gid)
                elif trans_table == '4':
                    gid_tt4.append(gid)
                elif trans_table == '25':
                    gid_tt25.append(gid)
                elif not trans_table or trans_table == 'none' or trans_table =='':
                    gid_tt_none.append(gid)
                else:
                    print('Error: unrecognized translation table', gid, trans_table)

        print('tt11', len(gid_tt11))
        print('gid_tt4', len(gid_tt4))
        print('gid_tt25', len(gid_tt25))
        print('tt_none', len(gid_tt_none))

        checkm_cmds = []

        #for tt, gid_tt in [('11', gid_tt11), ('4', gid_tt4), ('25', gid_tt25), ('none', gid_tt_none)]:
        for tt, gid_tt in [('4', gid_tt4), ('25', gid_tt25)]:

            print(tt)

            batch_num = 0
            num_genomes = 0
            batch_dir = os.path.join(output_dir,f'checkm2_batch/batch-{batch_num}_tt-{tt}')
            os.makedirs(batch_dir, exist_ok=True)
            out_dir = os.path.join(output_dir,f'checkm2/batch-{batch_num}_tt-{tt}')
            os.makedirs(out_dir, exist_ok=True)
            for gid in gid_tt:
                if num_genomes == 5000:
                    if tt != 'none':
                        cmd = f'checkm2 predict --input {batch_dir} --output-directory {out_dir} --ttable {tt} -x .gz -t {self.cpus}'
                    else:
                        cmd = f'checkm2 predict --input {batch_dir} --output-directory {out_dir} -x .gz -t {self.cpus}'

                    checkm_cmds.append(cmd)

                    num_genomes = 0
                    batch_num += 1
                    batch_dir = os.path.join(output_dir,f'checkm2_batch/batch-{batch_num}_tt-{tt}')
                    os.makedirs(batch_dir, exist_ok=True)
                    out_dir = os.path.join(output_dir,f'checkm2/batch-{batch_num}_tt-{tt}')

                out_file = os.path.join(batch_dir, f'{gid}.fna.gz')
                os.system(f'ln -s {genomic_files[gid]} {out_file}')
                num_genomes += 1

            if tt != 'none':
                cmd = f'checkm2 predict --input {batch_dir} --output-directory {out_dir} --ttable {tt} -x .gz -t {self.cpus}'
            else:
                cmd = f'checkm2 predict --input {batch_dir} --output-directory {out_dir} -x .gz -t {self.cpus}'

            checkm_cmds.append(cmd)

        # run CheckM v2
        fout = open(os.path.join(output_dir,'checkm_cmds.lst'), 'w')
        for cmd in checkm_cmds:
            fout.write(f'{cmd}\n')
        fout.close()


    def join_checkm2_results(self,checkm2_dir,output_dir):
        checkm2_files = []
        # list all directories starting with 'batch-'
        for d in os.listdir(checkm2_dir):
            if d.startswith('batch-'):
                # we get the 'quality_report.tsv' file in the directory
                checkm2_files.append(os.path.join(checkm2_dir, d, 'quality_report.tsv'))

        # join all 'quality_report.tsv' files
        checkm2_output = os.path.join(output_dir, 'checkm2.quality_report_merged.tsv')
        for i, file in enumerate(checkm2_files):
            with open(file) as f:
                if i != 0:
                    f.readline()
                for line in f:
                    with open(checkm2_output, 'a') as out:
                        out.write(line)
        self.logger.info('CheckM2 results written to: %s' % checkm2_output)

        # join only the first 4 columns (Name,Completeness,Contamination,Completeness_Model_Used) of all 'quality_report.tsv' files
        # and rename them genome_id,checkm2_completeness,checkm2_contamination,checkm2_model
        checkm2_output = os.path.join(output_dir, 'checkm2.quality_report_for_database.tsv')
        with open(checkm2_output, 'w') as out:
            out.write('genome_id\tcheckm2_completeness\tcheckm2_contamination\tcheckm2_model\n')
            for i, file in enumerate(checkm2_files):
                with open(file) as f:
                    f.readline()
                    for line in f:
                        infos = line.strip().split('\t')
                        gid = infos[0].replace('.fna','')
                        if 'General Model' in infos[3]:
                            infos[3] = 'General'
                        elif 'Specific Model' in infos[3]:
                            infos[3] = 'Specific'
                        out.write(f'{gid}\t{infos[1]}\t{infos[2]}\t{infos[3]}\n')
        self.logger.info('CheckM2 results for database written to: %s' % checkm2_output)




