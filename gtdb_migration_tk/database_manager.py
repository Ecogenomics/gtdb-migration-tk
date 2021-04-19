#!/usr/bin/env python

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

__prog_name__ = 'update_database_from_ftp.py'
__prog_desc__ = ('Update the GTDB with the latest genome downloaded from FTP.' +
                 'Before this update, make sure all metadata have been generated and CheckM did run on all new genomes')

__author__ = 'Pierre Chaumeil'
__copyright__ = 'Copyright 2016'
__credits__ = ['Pierre Chaumeil']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Pierre Chaumeil'
__email__ = 'p.chaumeil@qfab.org'
__status__ = 'Development'

import os
import hashlib
import re
import glob
import datetime
import ntpath
import multiprocessing as mp
import logging

from atpbar import atpbar
from atpbar import flush
import threading


from biolib.common import remove_extension
from dateutil.parser import parse
from gtdb_migration_tk.database_configuration import GenomeDatabaseConnectionFTPUpdate

from atpbar import register_reporter, find_reporter, flush


class DatabaseManager(object):

    def __init__(self, user, hostname, db, password, ftp_download_date, repository, path_to_log, cpus):

        self.repository = repository
        # By default we set the id to genbank (it is either 2 or 3 )
        self.id_database = 3
        if repository == "refseq":
            self.id_database = 2
        self.domains = ["archaea", "bacteria"]
        self.report_database_update = open(os.path.join(path_to_log,
                                                        "report_{0}_{1}_update_db.log".format(repository, ftp_download_date)), "w")

        self.password = password
        self.hostname = hostname
        self.user = user
        self.db = db

        self.temp_con = GenomeDatabaseConnectionFTPUpdate.GenomeDatabaseConnectionFTPUpdate(
            hostname, user, password, db)
        self.temp_con.MakePostgresConnection()
        self.temp_cur = self.temp_con.cursor()
        self.cpus = cpus

        self.logger = logging.getLogger('timestamp')

    def runUpdate(self, checkm, genome_dirs_file, dl_date):

        self.update_date = self.parse_date(dl_date)
        self.dict_existing_records = self._populateExistingRecords()
        self.list_checkm_records = self._populateNewRecords(checkm)
        self.genome_dirs_dict = self._populateGenomeDirs(genome_dirs_file)
        #======================================================================
        # Check if the genome is an existing genome
        # short_checkm_records list the records that are either to add or
        # version
        #======================================================================
        self.logger.info('Update existing genomes.')
        self.short_checkm_records = self._updateExistingGenomes()
        self.logger.info('Time to commit!                     ')
        self.temp_con.commit()
        self.logger.info('Update existing genomes:done.')
        self.logger.info('Add or Version new genomes.')
        self._addOrVersionNewGenomes()
        self.temp_con.commit()
        self.logger.info('Add or Version new genomes:done.')

        # Because we have added and updated script we repopulate
        # dict_existing_records
        dict_existing_records = self._populateExistingRecords()
        self.logger.info('Check Path or Remove records.')

        self._checkPathorRemoveRecord()
        self.logger.info('Check Path or Remove records: Done.')
        self.temp_con.commit()
        self.report_database_update.close()

    def _updateExistingGenomes(self):

        workerQueue = mp.Queue()
        writerQueue = mp.JoinableQueue()
        reportlist = mp.Manager().list()
        shortcheckmlist = mp.Manager().list()
        tasklist = mp.Manager().list()

        for record in self.list_checkm_records:
            workerQueue.put(record)

        for _ in range(self.cpus):
            workerQueue.put(None)

        try:
            workerProc = [mp.Process(target=self.worker_updateExistingGenomes, args=(
                workerQueue, i, tasklist, reportlist, shortcheckmlist, writerQueue)) for i in range(self.cpus)]
            writeProc = mp.Process(target=self.__progress, args=(
                len(self.list_checkm_records), writerQueue))

            writeProc.start()
            for p in workerProc:
                #print('starting', p.name, '=', p.is_alive())
                p.start()

            for p in workerProc:
                #print('stopping', p.name, '=', p.is_alive())
                p.join()

            writerQueue.put(None)
            writeProc.join()

            taskProc = []
            for i, list_sql in enumerate(tasklist):
                subproc = threading.Thread(
                    target=self.task_sql_command, args=(list_sql, i))
                subproc.start()
                taskProc.append(subproc)
            for tp in taskProc:
                tp.join()
            flush()
        except:
            for p in workerProc:
                p.terminate()

            writeProc.terminate

        self.logger.info('We write a report')
        for reportlist_item in reportlist:
            self.report_database_update.write(reportlist_item)
        result = []
        for it in shortcheckmlist:
            result.append(it)
        return result

    def task_sql_command(self, list_sql, process_idx):
        thread_con = GenomeDatabaseConnectionFTPUpdate.GenomeDatabaseConnectionFTPUpdate(
            self.hostname, self.user, self.password, self.db)
        thread_con.MakePostgresConnection()
        thread_cur = thread_con.cursor()

        list_subcommands = list(self.chunks(list_sql, 10))

        for subsql in atpbar(list_subcommands, name="Process-{}".format(process_idx)):
            big_sql_command = ';'.join(subsql)
            thread_cur.execute(big_sql_command)
        thread_con.commit()

    def worker_updateExistingGenomes(self, queue_in, process_idx, tasklist, list_report, list_shortcheckm, queue_out):
        list_sql = []
        while True:
            checkm_record = queue_in.get(block=True, timeout=None)
            if checkm_record == None:
                queue_out.task_done
                break
            if (checkm_record in self.dict_existing_records) and (checkm_record in self.genome_dirs_dict):
                list_report.append("{0}\t{1}\tupdate protein file\n".format(
                    self.repository, checkm_record))
                path_gtdb = re.sub(
                    r"(^.+\/)(GCA\/|GCF\/)", r"\g<2>", self.genome_dirs_dict[checkm_record])
                path_gtdb += "/" + os.path.basename(path_gtdb)
                path_database = re.sub(
                    r"(.+)(_genomic.fna)", r"\g<1>", self.dict_existing_records[checkm_record])
                path_protein_database = re.sub(
                    r"(.+)/(GC._[^_]+)(.*)", r"\g<1>/prodigal/\g<2>_protein.faa", path_database)
                # If the record is in a different folder , we need to change
                # it's path in the database
                if path_database not in path_gtdb:
                    query = "update genomes set fasta_file_location = replace(fasta_file_location, '{0}', '{1}') where name like '{2}'".format(
                            path_database, path_gtdb, checkm_record)
                    list_sql.append(query)
                    query = "update genomes set genes_file_location = '{0}' where name like '{1}'".format(
                            path_protein_database, checkm_record)
                    list_sql.append(query)

                # if the records is in the Checkm folder that means genomics
                # and protein files have changed. We need to re write their
                # sha256 values
                genomic_files = glob.glob(
                    self.genome_dirs_dict[checkm_record] + "/*_genomic.fna")
                if len(genomic_files) == 1:
                    genomic_file = genomic_files[0]
                    new_md5_genomic = self.sha256Calculator(genomic_file)
                    query = "update genomes set fasta_file_sha256 = '{0}' where name like '{1}'".format(
                        new_md5_genomic, checkm_record)
                    list_sql.append(query)
                _genome_path, genome_id = ntpath.split(
                    self.genome_dirs_dict[checkm_record])
                genome_id = genome_id[0:genome_id.find('_', 4)]
                gene_file_path = os.path.join(
                    self.genome_dirs_dict[checkm_record], "prodigal")
                gene_files = glob.glob(gene_file_path + "/*_protein.faa")
                if len(gene_files) == 1:
                    gene_file = gene_files[0]
                    new_md5_gene = self.sha256Calculator(gene_file)
                    query = "update genomes set genes_file_sha256 = '{0}' where name like '{1}'".format(
                        new_md5_gene, checkm_record)
                    list_sql.append(query)
            else:
                list_shortcheckm.append(checkm_record)
            queue_out.put(checkm_record)

        tasklist.append(list_sql)

    def _checkPathorRemoveRecord(self):

        workerQueue = mp.Queue()
        writerQueue = mp.Queue()
        reportlist = mp.Manager().list()
        tasklist = mp.Manager().list()

        for record in self.dict_existing_records:
            workerQueue.put(record)

        for _ in range(self.cpus):
            workerQueue.put(None)

        # try:
        workerProc = [mp.Process(target=self.worker_checkPathorRemoveRecord, args=(
            workerQueue, i, tasklist, reportlist, writerQueue)) for i in range(self.cpus)]
        writeProc = mp.Process(target=self.__progress, args=(
            len(self.dict_existing_records), writerQueue))

        writeProc.start()

        for p in workerProc:
            p.start()

        for p in workerProc:
            p.join()

        writerQueue.put(None)
        writeProc.join()

        taskProc = []
        for i, list_sql in enumerate(tasklist):
            subproc = threading.Thread(
                target=self.task_sql_command, args=(list_sql, i))
            subproc.start()
            taskProc.append(subproc)
        for tp in taskProc:
            tp.join()
        flush()

#=========================================================================
#         except:
#             for p in workerProc:
#                 p.terminate()
#
#             writeProc.terminate
#=========================================================================

        self.logger.info('We write a report')
        for reportlist_item in reportlist:
            self.report_database_update.write(reportlist_item)

    def worker_checkPathorRemoveRecord(self, queue_in, process_idx, tasklist, list_report, queue_out):
        """Process each data item in parallel."""
        list_sql = []
        while True:
            record = queue_in.get(block=True, timeout=None)
            if record == None:
                break
            # if the record was part of the checkm file, it has already been
            # updated
            if record not in self.list_checkm_records:
                if record not in self.genome_dirs_dict:
                    list_sql = self._removeRecord(
                        record, list_report, list_sql)
                else:
                    list_sql = self._checkPathRecord(
                        record, list_sql, list_report, self.dict_existing_records[record], self.genome_dirs_dict[record])
            # allow results to be processed or written to file
            queue_out.put(record)
        tasklist.append(list_sql)

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
            print(statusStr, end='\r')

    def _removeRecord(self, record, list_report, tasklist):
        list_report.append(
            "{0}\t{1}\tremoved\n".format(self.repository, record))
        query_delete = (
            "DELETE FROM genomes WHERE name LIKE '{0}'".format(record))
        tasklist.append(query_delete)
        return tasklist

    def _checkPathRecord(self, record, list_sql, list_report, path_in_db, path_in_folder):
        if path_in_db not in path_in_folder:
            path_in_folder = re.sub(
                r"(^.+\/)(GCA\/|GCF\/)", r"\g<2>", path_in_folder)
            path_in_folder += "/" + os.path.basename(path_in_folder)
            path_in_db = re.sub(r"(.+)(_genomic.fna)", r"\g<1>", path_in_db)
            query = "update genomes set fasta_file_location = replace(fasta_file_location, '{0}', '{1}') where id_at_source like '{2}'".format(
                    path_in_db, path_in_folder, record)
            list_report.append("{0}\t{1}\tupdate path\t{2}\t{3}\n".format(
                self.repository, record, path_in_db, path_in_folder))
            list_sql.append(query)
            pathinfo = path_in_folder.rsplit('/', 1)
            genes_path = os.path.join(
                pathinfo[0], 'prodigal', record + "_protein.faa").replace("\\", "/")
            query = "update genomes set genes_file_location = '{0}' where id_at_source like '{1}'".format(
                genes_path, record)
            list_sql.append(query)
        return list_sql

    def _addOrVersionNewGenomes(self):
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()
        reportlist = mp.Manager().list()
        tasklist = mp.Manager().list()

        for record in self.list_checkm_records:
            workerQueue.put(record)

        for _ in range(self.cpus):
            workerQueue.put(None)

        try:
            workerProc = [mp.Process(target=self.worker_addOrVersionNewGenomes, args=(
                workerQueue, i, tasklist, reportlist, writerQueue)) for i in range(self.cpus)]
            writeProc = mp.Process(target=self.__progress, args=(
                len(self.list_checkm_records), writerQueue))

            writeProc.start()

            for p in workerProc:
                p.start()

            for p in workerProc:
                p.join()

            writerQueue.put(None)
            writeProc.join()

            taskProc = []
            for i, list_sql in enumerate(tasklist):
                subproc = threading.Thread(
                    target=self.task_sql_command, args=(list_sql, i))
                subproc.start()
                taskProc.append(subproc)
            for tp in taskProc:
                tp.join()
            flush()

        except:
            for p in workerProc:
                p.terminate()

            writeProc.terminate

        self.logger.info('We write a report')
        for reportlist_item in reportlist:
            self.report_database_update.write(reportlist_item)

    def worker_addOrVersionNewGenomes(self, queue_in, process_idx, tasklist, list_report, queue_out):
        list_sql = []
        thread_con = GenomeDatabaseConnectionFTPUpdate.GenomeDatabaseConnectionFTPUpdate(
            self.hostname, self.user, self.password, self.db)
        thread_con.MakePostgresConnection()
        thread_cur = thread_con.cursor()

        while True:
            checkm_record = queue_in.get(block=True, timeout=None)
            if checkm_record == None:
                break
            if (checkm_record not in self.dict_existing_records) and (checkm_record in self.genome_dirs_dict):
                check_record_base = checkm_record.rsplit(".", 1)[0]
                id_record = self._checkPreviousVersion(
                    thread_cur, check_record_base)
                if id_record < 0:  # -1
                    # we add the genome to the database
                    list_sql = self._addNewGenomes(
                        checkm_record, list_report, list_sql)
                else:
                    list_sql = self._addNewGenomes(
                        checkm_record, list_report, list_sql, id_record)
            queue_out.put(checkm_record)
        tasklist.append(list_sql)

    def _addNewGenomes(self, checkm_record, list_report, list_sql, id_record=None):
        list_genome_details = [checkm_record]
        list_genome_details.append('')  # description
        list_genome_details.append(True)
        list_genome_details.append(None)
        fasta_file_path = os.path.join(self.genome_dirs_dict[checkm_record],
                                       os.path.basename(self.genome_dirs_dict[checkm_record]) + "_genomic.fna")
        fasta_file_path_shorten = re.sub(
            r"(.+/)(GCA\/|GCF\/)", r"\g<2>", fasta_file_path)
        list_genome_details.append(fasta_file_path_shorten)
        list_genome_details.append(self.sha256Calculator(fasta_file_path))
        list_genome_details.append(self.id_database)
        list_genome_details.append(checkm_record)
        list_genome_details.append(self.update_date)
        list_genome_details.append(True)
        list_genome_details.append(self.update_date)
        _genome_path, genome_id = ntpath.split(
            self.genome_dirs_dict[checkm_record])
        genome_id = genome_id[0:genome_id.find('_', 4)]
        gene_file_path = os.path.join(
            self.genome_dirs_dict[checkm_record], "prodigal", genome_id + "_protein.faa")
        gene_file_path_shorten = re.sub(
            r"(.+/)(GCA\/|GCF\/)", r"\g<2>", gene_file_path)
        list_genome_details.append(gene_file_path_shorten)
        list_genome_details.append(self.sha256Calculator(gene_file_path))
        list_genome_details.append(self.normaliseID(checkm_record))
        if id_record is None:
            list_report.append("{0}\t{1}\tadd\n".format(
                self.repository, checkm_record))
            list_sql.append(self.temp_cur.mogrify("INSERT INTO genomes " +
                                                  "(name,description,owned_by_root,owner_id,fasta_file_location,fasta_file_sha256,genome_source_id,id_at_source,date_added,has_changed,last_update,genes_file_location,genes_file_sha256,formatted_source_id) " +
                                                  "VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)", list_genome_details).decode("utf-8"))
        else:
            list_report.append(
                "{0}\t{1}\tversion\n".format(self.repository, checkm_record))
            list_sql.append(self.temp_cur.mogrify("UPDATE genomes " +
                                                  "SET name = %s,description = %s, " +
                                                  "owned_by_root = %s,owner_id = %s,  " +
                                                  "fasta_file_location = %s,fasta_file_sha256 = %s, " +
                                                  "genome_source_id = %s,id_at_source = %s,  " +
                                                  "date_added = %s,has_changed = %s,  " +
                                                  "last_update = %s,genes_file_location = %s,  " +
                                                  "genes_file_sha256 = %s, " +
                                                  "formatted_source_id = %s WHERE id = {0}; ".format(id_record), list_genome_details).decode("utf-8"))

            list_sql.append(
                "DELETE FROM aligned_markers where genome_id = {0};".format(id_record))
        return list_sql

    def _populateGenomeDirs(self, genome_dirs_file):
        with open(genome_dirs_file, 'r') as gen_file:
            gen_dict = {gen_line.split("\t")[0]: gen_line.split("\t")[1].strip()
                        for gen_line in gen_file}
        return gen_dict

    def _populateNewRecords(self, checkm):
        list_result = []
        checkm_fh = open(checkm, "r")
        checkm_fh.readline()
        for line in checkm_fh:
            full_name = line.split("\t")[0]
            name = full_name.split("_")[0] + "_" + full_name.split("_")[1]
            list_result.append(name)
        return list_result

    def _populateExistingRecords(self):
        self.temp_cur.execute("SELECT  gen.name,gen.fasta_file_location " +
                              "FROM genomes as gen " +
                              "WHERE gen.genome_source_id = {0} ;".format(self.id_database))

        dict_records = {key: os.path.dirname(
            value) for (key, value) in self.temp_cur}
        return dict_records


# Tools

    def _checkPreviousVersion(self, thread_cur, checkm_record):
        thread_cur.execute("SELECT  gen.id " +
                           "FROM genomes as gen " +
                           "WHERE gen.id_at_source like '{0}.%' ;".format(checkm_record))
        list_result = [record for (record,) in thread_cur]
        if len(list_result) > 0:
            return list_result[0]
        else:
            return -1

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

    def parse_date(self, date_text):
        try:
            datetime.datetime.strptime(date_text, '%Y-%m-%d')
        except ValueError:
            raise ValueError("Incorrect data format, should be YYYY-MM-DD")
        return parse(date_text)

    def chunks(self, lst, n):
        """Yield successive n-sized chunks from lst."""
        for i in range(0, len(lst), n):
            yield lst[i:i + n]

    def normaliseID(self, accession):
        normaccession = "G" + accession[4:accession.find('.', 0)]
        return normaccession
