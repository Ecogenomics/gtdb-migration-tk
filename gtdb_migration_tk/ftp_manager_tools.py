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
import hashlib
import glob
import gzip
import sys
import argparse
import tempfile
from datetime import datetime
import multiprocessing as mp
import logging


class FTPTools(object):
    def __init__(self, report, genomes_to_review, genome_domain_dict):
        self.reports = ("_assembly_report.txt",
                        "_assembly_stats.txt", "_hashes.txt")

        self.genomic_ext = "_genomic.fna.gz"
        self.protein_ext = "_protein.faa.gz"
        self.cds_ext = "_cds_from_genomic.fna.gz"
        self.rna_ext = "_rna_from_genomic.fna.gz"

        self.extensions = ("_feature_table.txt.gz", "_genomic.gbff.gz",
                           "_genomic.gff.gz", "_protein.gpff.gz", "_wgsmaster.gbff.gz")

        self.fastaExts = (self.genomic_ext, self.protein_ext)
        self.extrafastaExts = (self.cds_ext, self.rna_ext)
        self.allExts = self.fastaExts + self.extensions + self.reports
        self.allbutFasta = self.extensions + self.reports
        self.report = report
        self.genomes_to_review = genomes_to_review
        self.genome_domain_dict = genome_domain_dict
        self.dry_run = True

    def rreplace(self, s, old, new, occurrence):
        '''
         Instead of using the normal replace function, we need to implement our own.
         Some folder are named with a .gz in the middle so we only need to replace the last .gz in the string name
        :param s:
        :param old:
        :param new:
        :param occurrence:
        '''
        li = s.rsplit(old, occurrence)
        return new.join(li)

    def compareGenomes(self, intersect_list, old_dict, new_dict, ftp_directory, new_directory, threads):
        '''
        compare the genomes existing in both folders ( FTP folder and previous gtdb update).

        :param intersect_list:
        :param old_dict:
        :param new_dict:
        '''

        # populate worker queue with data to process
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        for gca_record in intersect_list:
            gtdb_dir = old_dict.get(gca_record)
            ftp_dir = new_dict.get(gca_record)
            target_dir = os.path.join(
                new_directory, ftp_dir.replace(ftp_directory, ''))

            workerQueue.put((gtdb_dir, ftp_dir, target_dir,
                             gca_record))


        for _ in range(threads):
            workerQueue.put((None, None, None, None))

        try:
            workerProc = [mp.Process(target=self.__workerThread, args=(
                workerQueue, writerQueue)) for _ in range(threads)]
            writeProc = mp.Process(target=self.__writerThread, args=(
                len(intersect_list), writerQueue))

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

    def __workerThread(self, queueIn, queueOut):
        """Process each data item in parallel."""
        while True:
            gtdb_dir, ftp_dir, target_dir, gca_record = queueIn.get(
                block=True, timeout=None)
            if gca_record is None:
                break

            status_gca = self.readmd5Checksum(
                gtdb_dir, ftp_dir, target_dir, gca_record)
            queueOut.put(status_gca)

    def __writerThread(self, numgenometoprocess, writerQueue):
        """Store or write results of worker threads in a single thread."""
        processedItems = 0
        while True:
            a = writerQueue.get(block=True, timeout=None)
            if a is None:
                break
            self.report.write(a)

            processedItems += 1
            statusStr = 'Finished processing %d of %d (%.2f%%) genome pairs.' % (
                processedItems, numgenometoprocess, float(processedItems) * 100 / numgenometoprocess)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

        sys.stdout.write('\n')

    def addGenomes(self, added_dict, ftp_dir, new_directory, genome_domain_dict):
        '''
        addGenomes function insert new genomes in the GTDB database. New genomes are present in the FTP folder
        but not in the previous version of GTDB.

        :TODO: Check if the new genome is a new version of an existing genome. in that case we overwrite the previous one
        and keep the same database id
        This will cause a conflict with the removeGenomes function.


        :param added_dict: dictionary of genomes to be added (genome_id:path to genome)
        :param ftp_dir: base directory leading the the FTP repository for refseq
        :param new_directory:base directory leading the new repository for refseq
        '''

        processedItems = 0
        for gcf_record in added_dict:
            processedItems += 1
            statusStr = 'Finished processing %d of %d (%.2f%%) genome to add.' % (
                processedItems, len(added_dict), float(processedItems) * 100 / len(added_dict))
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

            target_dir = os.path.join(
                new_directory, added_dict[gcf_record].replace(ftp_dir, ''))
            shutil.copytree(added_dict[
                gcf_record], target_dir, ignore=shutil.ignore_patterns("*_assembly_structure"))
            self.report.write(
                "{0}\t{1}\tnew\n".format(genome_domain_dict.get(gcf_record).upper(), gcf_record))
            for compressed_file in glob.glob(target_dir + "/*.gz"):
                if not os.path.isdir(compressed_file):
                    inF = gzip.open(compressed_file, 'rb')
                    outF = open(
                        self.rreplace(compressed_file, ".gz", "", 1), 'wb')
                    outF.write(inF.read())
                    inF.close()
                    outF.close()
                    os.remove(compressed_file)
        sys.stdout.write('\n')

    def removeGenomes(self, removed_dict, old_dict_domain):
        '''
        removeGenomes function removes all outdated genomes from the gtdb database
        In addition it tracks the lists(name and owner) that have been modified while deleting those genomes

        :param removed_dict: dictionary of genomes to delete
        '''

        for gca_record in removed_dict:
            self.report.write(
                "{0}\t{1}\tremoved\n".format(old_dict_domain.get(gca_record).upper(), gca_record))

    def readmd5Checksum(self, gtdb_dir, ftp_dir, target_dir, genome_record):
        '''
        Compare the checksum of the file listed in the checksums.txt
        '''
        pathftpmd5 = os.path.join(ftp_dir, "md5checksums.txt")
        pathgtdbmd5 = os.path.join(gtdb_dir, "md5checksums.txt")
        target_pathnewmd5 = os.path.join(target_dir, "md5checksums.txt")
        status = []

        tmp_ftp_dir = tempfile.mkdtemp()
        tmp_target = os.path.join(tmp_ftp_dir, os.path.basename(target_dir))
        shutil.copytree(ftp_dir, tmp_target, symlinks=True,
                        ignore=shutil.ignore_patterns("*_assembly_structure"))
        for compressed_file in glob.glob(tmp_target + "/*.gz"):
            if os.path.isdir(compressed_file) is False:
                inF = gzip.open(compressed_file, 'rb')
                try:
                    outF = open(
                        self.rreplace(compressed_file, ".gz", "", 1), 'wb')
                except IOError:
                    os.chmod(
                        self.rreplace(compressed_file, ".gz", "", 1), 0o775)
                    outF = open(
                        self.rreplace(compressed_file, ".gz", "", 1), 'wb')
                outF.write(inF.read())
                inF.close()
                outF.close()
                os.remove(compressed_file)

        ftpdict, ftpdict_fasta, ftpdict_faa, ftpdict_extra_fasta = self.parse_checksum(
            tmp_target)
        gtdbdict, gtdbdict_fasta, gtdbdict_faa, gtdbdict_extra_fasta = self.parse_checksum(
            gtdb_dir)

        # if the genomic.fna.gz or the protein.faa.gz are missing, we set this
        # record as incomplete
        if len(list(set(ftpdict_fasta.keys()).symmetric_difference(set(gtdbdict_fasta.keys())))) > 0:
            self.genomes_to_review.write("tmp_target:{}\nftpdict_fasta.keys():{}\ngtdb_dir:{}\ngtdbdict_fasta.keys():{}\n\n".format(tmp_target,
                                                                                                                                    ftpdict_fasta.keys(),
                                                                                                                                    gtdb_dir,
                                                                                                                                    gtdbdict_fasta.keys()))
            status.append("incomplete")
            shutil.copytree(tmp_target, target_dir, symlinks=True,
                            ignore=shutil.ignore_patterns("*_assembly_structure"))
            # we unzip of gz file

        else:
            ftp_folder = False
            # check if genomic.fna.gz and protein.faa.gz are similar between
            # previous gtdb and ftp
            for key, value in ftpdict_fasta.items():
                if value != gtdbdict_fasta.get(key):
                    ftp_folder = True

            # if one of the 2 files is different than the previous version , we
            # use the ftp record over the previous gtdb one , we then need to
            # re run the metadata generation
            if ftp_folder:
                if os.path.exists(target_dir):
                    shutil.rmtree(target_dir)
                shutil.copytree(
                    tmp_target, target_dir, symlinks=True,
                    ignore=shutil.ignore_patterns("*_assembly_structure"))
                status.append("modified")

            else:
                # The 2 main fasta files haven't changed so we can copy the old
                # gtdb folder over
                if os.path.exists(target_dir):
                    shutil.rmtree(target_dir)
                shutil.copytree(
                    gtdb_dir, target_dir, symlinks=True,
                    ignore=shutil.ignore_patterns("*_assembly_structure"))
                if not os.path.isdir(os.path.join(target_dir,'prodigal')):
                    status.append("unmodified")



                # We check if all other file of this folder are the same.
                checksum_changed = False

                for key, value in ftpdict_faa.items():
                    if value != gtdbdict_faa.get(key):
                        checksum_changed = True
                        shutil.copy2(
                            os.path.join(tmp_target, key), os.path.join(target_dir, key))
                        status.append("new_protein")

                for key, value in ftpdict.items():
                    if value != gtdbdict.get(key):
                        checksum_changed = True
                        shutil.copy2(
                            os.path.join(tmp_target, key), os.path.join(target_dir, key))
                        status.append("new_metadata")

                for key, value in ftpdict_extra_fasta.items():
                    if value != gtdbdict_extra_fasta.get(key):
                        checksum_changed = True
                        shutil.copy2(
                            os.path.join(tmp_target, key), os.path.join(target_dir, key))
                        status.append("new_cds_rna")
                # we copy the new checksum
                if checksum_changed:
                    try:
                        shutil.copy2(pathgtdbmd5, target_pathnewmd5)
                    except IOError:
                        os.chmod(target_pathnewmd5, 0o664)
                        shutil.copy2(pathftpmd5, target_pathnewmd5)
                for report in self.reports:
                    target_files = glob.glob(
                        os.path.join(target_dir, "*" + report))
                    ftp_files = glob.glob(os.path.join(ftp_dir, "*" + report))
                    if len(target_files) == 1 and len(ftp_files) == 1:
                        status = self.comparesha256(
                            ftp_files[0], target_files[0], status)
                    elif len(target_files) == 0 and len(ftp_files) == 0 and report == '_hashes.txt':
                        status.append("old_folder_dir")
                    elif len(target_files) == 0 and len(ftp_files) == 1 and report == '_hashes.txt':
                        shutil.copy2(ftp_dir[0], ftp_dir[0].replace(
                            ftp_dir.target_dir))
                        status.append("new_hashes")
                    else:
                        print("########")
                        print(target_dir)
                        print(target_files)
                        print(ftp_dir)
                        print(ftp_files)
                        print(f"IT SHOULDN'T HAPPEN ({target_dir},{ftp_dir}) ")
                        print("########")
                        status.append("to_curate")
        status_record = "{0}\t{1}\t{2}\n".format(self.genome_domain_dict.get(
            genome_record).upper(), genome_record, ';'.join([x for x in set(status)]))
        shutil.rmtree(tmp_ftp_dir)
        return status_record

    def comparesha256(self, ftp_file, target_file, status):
        '''
        comparesha256 compares the report file
        :param ftp_file:
        :param target_file:
        :param status:

        '''
        original_checksum = hashlib.md5(
            open(ftp_file, 'rb').read()).hexdigest()
        gtdb_checksum = hashlib.md5(open(target_file, 'rb').read()).hexdigest()
        if original_checksum != gtdb_checksum:
            try:
                shutil.copy2(ftp_file, target_file)
            except IOError:
                os.chmod(target_file, 0o664)
                shutil.copy2(ftp_file, target_file)
            status.append("new_metadata")
        return status

    def parse_checksum(self, pathtodir):
        '''
        parse_checksum function parses the md5 checksum file.
        It returns 2 dictionaries {file:size} : one for the fna and faa files, one for the genbank files
        :param md5File:
        '''

        out_dict, out_dict_fasta, out_dict_faa, out_dict_extra_fasta = {}, {}, {}, {}

        for name in glob.glob(os.path.join(pathtodir, '*')):
            if name.endswith(self.extrafastaExts):
                out_dict_extra_fasta[os.path.basename(
                    name)] = self.sha256Calculator(name)
                os.chmod(name, 0o664)
            elif name.endswith(self.genomic_ext):
                out_dict_fasta[os.path.basename(
                    name)] = self.sha256Calculator(name)
                os.chmod(name, 0o664)
            elif name.endswith(self.protein_ext):
                out_dict_faa[os.path.basename(
                    name)] = self.sha256Calculator(name)
                os.chmod(name, 0o664)
            elif name.endswith(self.allbutFasta):
                out_dict[os.path.basename(name)] = self.sha256Calculator(name)
                os.chmod(name, 0o664)
        return (out_dict, out_dict_fasta, out_dict_faa, out_dict_extra_fasta)

    def sha256Calculator(self, file_path):

        try:
            filereader = open(file_path, "rb")
        except:
            raise Exception("Cannot open Fasta file: " + file_path)
        m = hashlib.sha256()
        for line in filereader:
            m.update(line)
        sha256_checksum = m.hexdigest()
        filereader.close()
        return sha256_checksum
