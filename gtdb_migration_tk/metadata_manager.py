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
import logging
import ntpath
import datetime

from biolib.common import check_file_exists, make_sure_path_exists
from biolib.parallel import Parallel

from gtdb_migration_tk.genometk_lite.metadata_nucleotide import MetadataNucleotide
from gtdb_migration_tk.genometk_lite.metadata_genes import MetadataGenes
from gtdb_migration_tk.genometk_lite.rna import RNA

class MetadataManager(object):
    """Create file indicating directory of each genome."""

    def __init__(self, cpus=1):
        self.cpus = cpus
        self.contig_break = 10
        self.logger = logging.getLogger('timestamp')
        self.starttime = None

    def generate_metadata(self, gtdb_genome_path_file):
        self.starttime = datetime.datetime.utcnow().replace(microsecond=0)
        input_files = []
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

            genome_file = os.path.join(gpath, assembly_id + '_genomic.fna')
            gff_file = os.path.join(gpath, 'prodigal', gid + '_protein.gff')

            input_files.append([genome_file, gff_file])

        # process each genome
        print('Generating metadata for each genome:')
        parallel = Parallel(cpus=self.cpus)
        parallel.run(self._producer,
                     None,
                     input_files,
                     self._progress)

    def _producer(self, input_files):
        """Process each genome."""

        genome_file, gff_file = input_files
        full_genome_dir, _ = ntpath.split(genome_file)

        # clean up old log files
        log_file = os.path.join(full_genome_dir, 'genometk.log')
        if os.path.exists(log_file):
            os.remove(log_file)

        # calculate metadata
        self.nucleotide(genome_file,full_genome_dir)
        self.gene(genome_file,gff_file,full_genome_dir)

        return full_genome_dir

    def nucleotide(self, genome_file,output_dir):
        #self.logger.info('Calculating nucleotide properties of genome.')

        check_file_exists(genome_file)
        make_sure_path_exists(output_dir)

        meta_nuc = MetadataNucleotide()
        metadata_values, metadata_desc = meta_nuc.generate(genome_file,
                                                           self.contig_break)

        # write statistics to file
        output_file = os.path.join(output_dir, 'metadata.genome_nt.tsv')
        fout = open(output_file, 'w')
        for field in sorted(metadata_values.keys()):
            fout.write('%s\t%s\n' % (field, str(metadata_values[field])))
        fout.close()

        # write description to file
        output_file = os.path.join(output_dir, 'metadata.genome_nt.desc.tsv')
        fout = open(output_file, 'w')
        for field in sorted(metadata_desc.keys()):
            fout.write('%s\t%s\t%s\n' % (field,
                                         metadata_desc[field],
                                         type(metadata_values[field]).__name__.upper()))
        fout.close()

    def gene(self, genome_file,gff_file,output_dir):
        #self.logger.info('Calculating gene properties of genome.')

        check_file_exists(genome_file)
        check_file_exists(gff_file)
        make_sure_path_exists(output_dir)

        meta_genes = MetadataGenes()
        metadata_values, metadata_desc = meta_genes.generate(genome_file,
                                                                gff_file)

        # write statistics to file
        output_file = os.path.join(output_dir, 'metadata.genome_gene.tsv')
        fout = open(output_file, 'w')
        for field in sorted(metadata_values.keys()):
            fout.write('%s\t%s\n' % (field, str(metadata_values[field])))
        fout.close()

        # write description to file
        output_file = os.path.join(output_dir, 'metadata.genome_gene.desc.tsv')
        fout = open(output_file, 'w')
        for field in sorted(metadata_desc.keys()):
            fout.write('%s\t%s\t%s\n' % (field,
                                         metadata_desc[field],
                                         type(metadata_values[field]).__name__.upper()))
        fout.close()



    def _progress(self, processed_items, total_items):
        current_time_utc = datetime.datetime.utcnow().replace(microsecond=0)
        if processed_items > 0:
            time_left= (current_time_utc - self.starttime) * (total_items-processed_items)/ processed_items
            return '  Processed {} of {} ({}%) genomes. (ETA {})           '.format(processed_items,
                                                               total_items,
                                                               round(processed_items * 100.0 / total_items,2),time_left)
