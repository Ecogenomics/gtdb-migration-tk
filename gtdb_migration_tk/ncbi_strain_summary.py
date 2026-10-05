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

__prog_name__ = 'ncbi_strain_summary.py'
__prog_desc__ = 'Parse the strain identifiers and type material status of each genome'

__author__ = 'Pierre Chaumeil'
__copyright__ = 'Copyright 2017'
__credits__ = ['Pierre Chaumeil']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Pierre Chaumeil'
__email__ = 'p.chaumeil@uq.edu.au'
__status__ = 'Development'

"""The strain identifiers and type material status NCBI gives each genome of a release.

The strain identifiers are read from the assembly report's '# Infraspecific
name:' and '# Isolate:' lines, which hold what the assembly summary's
infraspecific_name and isolate columns do, and are split by
ncbi_utils.strain_identifiers(), which strains type_table reads the summaries
with. The GenBank files are not read: the /strain= and /isolate= qualifiers of
<assembly>_genomic.gbff and <assembly>_wgsmaster.gbff were looked for under
those names, which NCBI serves gzipped, so they were never found. Over r237 the
summaries give every strain ID the database holds for all but 51 of 1,047,596
genomes.
"""

import os
import sys
import logging
import multiprocessing as mp
import re

from gtdb_migration_tk.biolib_lite.common import get_num_lines
from gtdb_migration_tk.ncbi_utils import NCBI_NA, open_summary, strain_identifiers

class NCBIStrainParser(object):
    """Extract genes in nucleotide space."""

    def __init__(self, assembly_summary_bacteria_genbank, assembly_summary_archaea_genbank,
                 assembly_summary_bacteria_refseq, assembly_summary_archaea_refseq,cpus):
        self.genbank_dictionary = self.parse_summary(
            assembly_summary_bacteria_genbank, assembly_summary_archaea_genbank)
        self.refseq_dictionary = self.parse_summary(
            assembly_summary_bacteria_refseq, assembly_summary_archaea_refseq)
        self.cpus = cpus

        self.logger = logging.getLogger('timestamp')


    def parse_summary(self, assembly_bacteria_summary, assembly_archaea_summary):
        assembly_summary_dict = {}
        for assembly_summary in [assembly_bacteria_summary, assembly_archaea_summary]:
            with open_summary(assembly_summary) as as_file:
                as_file.readline()
                for line in as_file:
                    if line.startswith('#'):

                        line = line.replace('# ', '')
                        line = line.replace('#', '')
                        headers = line.strip('\n').split('\t')
                        index_genome_id = headers.index('assembly_accession')
                        index_excluded_from_refseq = headers.index('excluded_from_refseq')
                        relation_to_type_material_index = headers.index('relation_to_type_material')
                    else:
                        line_infos = line.strip('\n').split('\t')
                        assembly_summary_dict[line_infos[index_genome_id]
                                              ] = line_infos[relation_to_type_material_index]
        return assembly_summary_dict

    def generate_ncbi_strains_summary(self, genome_dir_file, out_dir):
        outf = open(os.path.join(out_dir, "strain_summary_file.tsv"), "w")
        outf.write(
            "genome_id\tOrganism name\tncbi_strain_identifiers\tncbi_type_material_designation\n")
        number_of_genomes = get_num_lines(genome_dir_file)
        count = 1
        lines_to_process = []
        with open(genome_dir_file, 'r') as genomelistfile:
            for line in genomelistfile:
                sys.stdout.write("{}% complete\r".format(
                    round((float(count) * 100 / number_of_genomes), 3)))
                # sys.stdout.write("{}/{}\r".format(count, num_lines))
                sys.stdout.flush()
                count += 1
                lines_to_process.append(line)


        print(f"number of cpus used:{self.cpus}")

        #populate worker queue with data to process
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()
        manager = mp.Manager()
        return_list = manager.list()
        substr_return_list = manager.list()


        for f in lines_to_process:
            workerQueue.put(f)

        for _ in range(self.cpus):
            workerQueue.put(None)

        try:
            workerProc = [mp.Process(target=self.ncbi_strain_worker,
                                     args=(workerQueue, writerQueue,return_list,substr_return_list))
                          for _ in range(self.cpus)]
            writeProc = mp.Process(target=self.__writerThread,
                                   args=(len(lines_to_process), writerQueue))

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


        list_lines_to_write =[x for x in return_list if x != 'null']

        #Print the substrains of interest:
        self.logger.info("Substrains of interest:")
        for substr in substr_return_list:
            self.logger.info('- '+substr)

        for line_to_write in list_lines_to_write:
            outf.write(line_to_write)
        outf.close()

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


    def ncbi_strain_worker(self, queueIn, queueOut,return_list,substr_return_list):
        while True:
            line = queueIn.get(block=True, timeout=None)

            if line == None:
                break

            line_split = line.strip().split('\t')

            genome_id = line_split[0]
            genome_path = line_split[1]
            genome_dir_id = os.path.basename(os.path.normpath(genome_path))

            species = ''
            infraspecific_names = []
            isolate = ''
            report_file = os.path.join(
                genome_path, genome_dir_id + '_assembly_report.txt')
            if os.path.exists(report_file):
                with open(report_file, 'r') as report:
                    for report_line in report:
                        if report_line.startswith('# Organism name: '):
                            species = re.sub(
                                r'\([^)]*\)', '', report_line.replace("# Organism name:", "")).strip()
                        elif report_line.startswith('# Infraspecific name:'):
                            infraspecific_names.append(
                                report_line.replace("# Infraspecific name:", "").strip())
                        elif report_line.startswith('# Isolate:'):
                            isolate = report_line.replace("# Isolate:", "").strip()

            infraspecific_name = '; '.join(infraspecific_names) or NCBI_NA
            if 'substr.' in infraspecific_name:
                substr_return_list.append("{1} strain {0}".format(infraspecific_name, genome_id))

            strain_ids = strain_identifiers(infraspecific_name, isolate or NCBI_NA)

            if genome_id.startswith('GCA'):
                typemat = self.genbank_dictionary.get(genome_id)
            else:
                typemat = self.refseq_dictionary.get(genome_id)

            line_to_write = "{0}\t{1}\t{2}\t{3}\n".format(
                genome_id, species, ';'.join(strain_ids), typemat)

            queueOut.put(line_to_write)
            return_list.append(line_to_write)
