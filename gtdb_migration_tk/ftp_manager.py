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
from gtdb_migration_tk.ftp_manager_tools import FTPTools


class RefSeqManager(object):

    def __init__(self, new_refseq_folder, dry_run=False, cpus=1):

        self.domains = ["archaea", "bacteria"]
        self.genome_domain_dict = {}

        self.threads = cpus

        self.genomic_ext = "_genomic.fna"
        self.protein_ext = "_protein.faa"
        self.cds_ext = "_cds_from_genomic.fna"
        self.rna_ext = "_rna_from_genomic.fna"

        self.fastaExts = (self.genomic_ext, self.protein_ext)
        self.extrafastaExts = (self.cds_ext, self.rna_ext)

        self.extensions = ("_feature_table.txt", "_genomic.gbff",
                           "_genomic.gff", "_protein.gpff", "_wgsmaster.gbff")
        self.reports = ("_assembly_report.txt",
                        "_assembly_stats.txt", "_hashes.txt")
        self.allExts = self.fastaExts + self.extensions + self.reports
        self.allbutFasta = self.extensions + self.reports

        self.dry_run = dry_run

        self.report_gcf = open(os.path.join(
            new_refseq_folder, "report_gcf.log"), "w", 1)
        self.genomes_to_review = open(os.path.join(
            new_refseq_folder, "gid_to_review.log"), "w", 1)

    def runComparison(self, ftp_refseq, new_refseq, ftp_genome_dirs, old_genome_dirs, archaea_assembly_summary, bacteria_assembly_summary):
        '''
        runComparison function is walking across all directories recursively
        only folder containing latest_assembly_versions but not containing _assembly_structure
        are of interest
        '''

        # old_dict lists all records from the previous gtdb update
        with open(old_genome_dirs, 'r') as old_file:
            old_dict = {old_line.split("\t")[0]: old_line.split("\t")[1].strip()
                        for old_line in old_file}
        print("{}:old_dict loaded".format(
            str(datetime.now())))

        # This is a one of for the transtion from the old directory structure
        # to the new one
        with open(old_genome_dirs, 'r') as old_file:
            old_dict_domain = {old_line.split("\t")[0]: old_line.split("\t")[1].strip().split("/")[8]
                               for old_line in old_file}
        print("{}:old_dict_domain loaded ({} records)".format(
            str(datetime.now()), len(old_dict_domain)))

        # new list list all records from the ftp folder and considered as
        # latest
        new_list = []
        with open(archaea_assembly_summary, 'r') as ftp_assembly_summary_file:
            ftp_assembly_summary_file.readline()
            new_list = [new_line.split(
                "\t")[0] for new_line in ftp_assembly_summary_file if new_line.split("\t")[10] == "latest"]
        self.genome_domain_dict = {arcid: "Archaea" for arcid in new_list}
        with open(bacteria_assembly_summary, 'r') as ftp_assembly_summary_file:
            ftp_assembly_summary_file.readline()
            bacterial_new_list = [new_line.split(
                "\t")[0] for new_line in ftp_assembly_summary_file if new_line.split("\t")[10] == "latest"]
        for bacid in bacterial_new_list:
            self.genome_domain_dict[bacid] = "Bacteria"
        new_list.extend(bacterial_new_list)
        print("{}:new_list loaded ({} records).".format(
            str(datetime.now()), len(new_list)))

        # new dict lists all records from FTP which are in new_list
        new_dict = {}
        cter = 0
        print("loading new_dict.....")
        with open(ftp_genome_dirs, 'r') as new_genome_dirs_file:
            for new_line in new_genome_dirs_file:
                cter += 1
                # print('{}'.format(cter), end='\r')
                new_line_split = new_line.split("\t")
                if new_line_split[0].startswith("GCF") and new_line_split[0] in new_list:
                    new_dict[new_line_split[0]] = new_line_split[1].strip()

        print("{}:new_dict loaded ({} records).".format(
            str(datetime.now()), len(new_dict)))

        ftptools = FTPTools(
            self.report_gcf, self.genomes_to_review, self.genome_domain_dict)

        # delete genomes from the Database
        removed_dict = {removed_key: old_dict[removed_key] for removed_key in list(
            set(old_dict.keys()) - set(new_dict.keys()))}
        ftptools.removeGenomes(removed_dict, old_dict_domain)

        #new genomes in FTP
        added_dict = {added_key: new_dict[added_key] for added_key in list(
            set(new_dict.keys()) - set(old_dict.keys()))}
        ftptools.addGenomes(added_dict, ftp_refseq,
                             new_refseq, self.genome_domain_dict)



        print("{}:Generating intersection list.....".format(str(datetime.now())))
        intersect_list = list(
            set(old_dict.keys()).intersection(set(new_dict.keys())))
        print("{}:Intersection list:{} genomes".format(
            str(datetime.now()), len(intersect_list)))
        ftptools.compareGenomes(
             intersect_list, old_dict, new_dict, ftp_refseq, new_refseq, self.threads)
        self.report_gcf.close()
        self.genomes_to_review.close()


class GenBankManager(object):

    def __init__(self, new_genbank_folder, dry_run=False, cpus=1):
        self.domains = ["archaea", "bacteria"]
        self.genome_domain_dict = {}

        self.genomic_ext = "_genomic.fna.gz"
        self.protein_ext = "_protein.faa.gz"
        self.cds_ext = "_cds_from_genomic.fna.gz"
        self.rna_ext = "_rna_from_genomic.fna.gz"

        self.fastaExts = (self.genomic_ext, self.protein_ext)
        self.extrafastaExts = (self.cds_ext, self.rna_ext)

        self.extensions = ("_feature_table.txt.gz", "_genomic.gbff.gz",
                           "_genomic.gff.gz", "_protein.gpff.gz", "_wgsmaster.gbff.gz")
        self.reports = ("_assembly_report.txt",
                        "_assembly_stats.txt", "_hashes.txt")
        self.allExts = self.fastaExts + self.extensions + self.reports
        self.allbutFasta = self.extensions + self.reports

        self.threads = cpus

        self.dry_run = dry_run

        self.logger = logging.getLogger('timestamp')

        self.report = open(os.path.join(new_genbank_folder,
                                        "extra_gbk_report_gcf.log"), "w")
        self.genomes_to_review = open(os.path.join(
            new_genbank_folder, "gcaid_to_review.log"), "w", 1)
        self.select_gca = open(os.path.join(
            new_genbank_folder, "gca_selection.log"), "w")

    def runComparison(self, ftp_genbank, new_genbank, ftp_genbank_genome_dirs, old_genbank_genome_dirs, new_refseq_genome_dirs, gbk_arc_assembly, gbk_bac_assembly):
        '''
        runComparison function is walking across all directories recursively
        only folder containing latest_assembly_versions but not containing _assembly_structure
        are of interest
        '''

        # old_dict lists all records from the previous gtdb update
        with open(old_genbank_genome_dirs, 'r') as old_file:
            old_dict = {old_line.split("\t")[0]: old_line.split("\t")[1].strip()
                        for old_line in old_file}
        self.logger.info("{}:old_dict loaded".format(str(datetime.now())))

        # This is a one of for the transition from the old directory structure
        # to the new one
        with open(old_genbank_genome_dirs, 'r') as old_file:
            old_dict_domain = {old_line.split("\t")[0]: old_line.split("\t")[1].strip().split("/")[8]
                               for old_line in old_file}
        self.logger.info(
            "{}:old_dict_domain loaded".format(str(datetime.now())))

        listGCA = self.parseAssemblySummary(
            gbk_arc_assembly, gbk_bac_assembly, new_refseq_genome_dirs)
        # new dict lists all records from FTP which are in new_list
        new_dict = {}
        num_lines = sum(1 for line in open(ftp_genbank_genome_dirs))
        processedItems = 0
        with open(ftp_genbank_genome_dirs, 'r') as new_genome_dirs_file:
            for new_line in new_genome_dirs_file:
                processedItems += 1
                statusStr = 'New Dictionary = Finished processing %d of %d (%.2f%%) gca records.' % (
                    processedItems, num_lines, float(processedItems) * 100 / num_lines)
                sys.stdout.write('%s\r' % statusStr)
                sys.stdout.flush()
                new_line_split = new_line.split("\t")
                if new_line_split[0].startswith("GCA") and new_line_split[0] in listGCA:
                    new_dict[new_line_split[0]] = new_line_split[1].strip()
        sys.stdout.write('\n')
        self.logger.info("{}:new_dict loaded".format(str(datetime.now())))

        ftptools = FTPTools(
            self.report, self.genomes_to_review, self.genome_domain_dict)

        # delete genomes from the Database
        removed_dict = {removed_key: old_dict[removed_key] for removed_key in list(
            set(old_dict.keys()) - set(new_dict.keys()))}
        self.logger.info("{0} genomes to remove".format(len(removed_dict)))
        ftptools.removeGenomes(removed_dict, old_dict_domain)

        # new genomes in FTP
        added_dict = {added_key: new_dict[added_key] for added_key in list(
            set(new_dict.keys()) - set(old_dict.keys()))}
        self.logger.info("{0} genomes to add".format(len(added_dict)))
        ftptools.addGenomes(added_dict, ftp_genbank,
                           new_genbank, self.genome_domain_dict)

        self.logger.info(
            "{}:Generating intersection list.....".format(str(datetime.now())))
        intersect_list = list(
            set(old_dict.keys()).intersection(set(new_dict.keys())))
        self.logger.info("{}:Intersection list:{} genomes".format(
            str(datetime.now()), len(intersect_list)))
        ftptools.compareGenomes(
             intersect_list, old_dict, new_dict, ftp_genbank, new_genbank, self.threads)
        self.select_gca.close()
        self.report.close()
        self.genomes_to_review.close()


# Tools

    def parseAssemblySummary(self, gbk_arc_assembly, gbk_bac_assembly, new_refseq_genome_dirs):
        listGCA = []
        dictGCF = self._populateGenomesDict(new_refseq_genome_dirs)
        self.logger.info("parsing of dictionary is done.")

        for domain, assemblyfile in [('archaea', gbk_arc_assembly), ('bacteria', gbk_bac_assembly)]:
            num_lines = sum(1 for line in open(assemblyfile))
            processedItems = 0
            with open(assemblyfile, "r") as sumf:

                # we discard the first line
                sumf.readline()
                for line in sumf:
                    processedItems += 1
                    statusStr = 'Finished processing %d of %d (%.2f%%) gca records.' % (
                        processedItems, num_lines, float(processedItems) * 100 / num_lines)
                    sys.stdout.write('%s\r' % statusStr)
                    sys.stdout.flush()
                    split_line = line.split("\t")
                    gcf_access = 'G' + split_line[17][4:13]
                    full_gca_access = split_line[0]
                    latest = split_line[10]
                    surveillance_info = split_line[20]
                    if 'surveillance' in surveillance_info:
                        continue
                    elif latest == "latest":
                        if not gcf_access.startswith("GCF"):
                            formatted_gcaid = 'G' + full_gca_access[4:13]
                            if formatted_gcaid in dictGCF:
                                self.select_gca.write("{0} skipped because {1} in Refseq (although {0} has no GCF)\n".format(
                                    full_gca_access, gcf_access))
                                continue
                            else:
                                listGCA.append(full_gca_access)
                        else:
                            # if the Refseq folder is empty, we copy the genbank
                            # folder
                            if gcf_access in dictGCF:
                                protein_files = glob.glob(
                                    os.path.join(dictGCF.get(gcf_access), "*_protein.faa"))
                                if len(protein_files) == 0:
                                    self.select_gca.write(
                                        "{0} associated with {1} : {1} missed files in FTP folder\n".format(full_gca_access, gcf_access))
                                    listGCA.append(full_gca_access)
                            else:
                                self.select_gca.write(
                                    "{0} associated with {1} : {1} not present in FTP folder\n".format(full_gca_access, gcf_access))
                                listGCA.append(full_gca_access)
                    self.genome_domain_dict[full_gca_access] = domain
        sys.stdout.write('\n')
        return listGCA

    def _populateGenomesDict(self, genome_dirs_file):
        temp_dict = {}
        with open(genome_dirs_file, "r") as list_dirs:
            for line in list_dirs:
                temp_dict['G' + line.split("\t")[0][4:13]] = line.split(
                    "\t")[1].rstrip()
        return temp_dict
