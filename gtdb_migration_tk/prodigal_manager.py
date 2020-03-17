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


from biolib.common import remove_extension
from biolib.external.execute import check_dependencies
from biolib.external.prodigal import Prodigal
from biolib.checksum import sha256

from gtdb_migration_tk.prettytable import PrettyTable

from gtdb_migration_tk.tools import Tools


class ProdigalManager(object):
    """Create file indicating directory of each genome."""

    def __init__(self, tmp_dir='/tmp/', cpus=1):
        """Initialization."""

        self.tmp_dir = tmp_dir
        self.cpus = cpus

        check_dependencies(['prodigal'])

        self.logger = logging.getLogger('timestamp')

    def run_prodigal_check(self, gtdb_genome_path_file):
        pretty = PrettyTable()
        pretty.field_names = ['Assembly Id', 'Prodigal Table', 'NCBI Table']
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

            prodigal_translation_table = os.path.join(
                gpath, 'prodigal', 'prodigal_translation_table.tsv')
            gff_file = os.path.join(gpath, assembly_id + '_genomic.gff')

            if os.path.isfile(prodigal_translation_table) and os.path.isfile(gff_file):
                if os.stat(gff_file).st_size == 0:
                    self.logger.warning(
                        'GFF appears to be empty: %s' % gff_file)
                    continue
                prodigal_table = self._parse_prodigal_translation_table(
                    prodigal_translation_table)
                ncbi_table = self._parse_ncbi_translation_table(gff_file)
                if prodigal_table is not None and ncbi_table is not None and prodigal_table != ncbi_table:
                    pretty.add_row([assembly_id, prodigal_table, ncbi_table])
        print(pretty)

    def run(self, gtdb_genome_path_file, all_genomes=False):
        # get path of genomes to process
        genome_paths = {}
        countr = 0

        tools = Tools()

        file_lgth = tools.file_len(gtdb_genome_path_file)

        for line in open(gtdb_genome_path_file):
            countr += 1
            statusStr = '{}/{} ({}%)lines read.'.format(countr,
                                                        file_lgth, round(float(countr) * 100 / file_lgth, 2))
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

            line_split = line.strip().split('\t')
            gid = line_split[0]
            gpath = line_split[1]
            aa_gene_file = os.path.join(
                gpath, 'prodigal', gid + '_protein.faa')

            assembly_id = os.path.basename(os.path.normpath(gpath))
            genome_file = os.path.join(
                gpath, gpath + '_genomic.fna')

            if all_genomes:
                prodigal_dir = os.path.join(gpath, 'prodigal')
                if os.path.exists(prodigal_dir):
                    shutil.rmtree(prodigal_dir)
                genome_paths.append(genome_file)
            else:
                if os.path.exists(aa_gene_file):
                    # verify checksum
                    checksum_file = aa_gene_file + '.sha256'
                    if os.path.exists(checksum_file):
                        checksum = sha256(aa_gene_file)
                        cur_checksum = open(
                            checksum_file).readline().strip()
                        if checksum == cur_checksum:
                            continue
                genome_paths[gid] = gpath
        sys.stdout.write('\n')

        # run Prodigal
        self._run_prodigal(genome_paths)

    def _run_prodigal(self, genome_paths):
        """Run Prodigal on genomes."""

        # get genome path and translation table for each file
        self.logger.info(
            'Determining genomic file and translation table for each of the %d genomes.' % len(genome_paths))
        genome_files = []
        translation_table = {}
        for gid, gpath in genome_paths.items():
            assembly_id = os.path.basename(os.path.normpath(gpath))
            canonical_gid = assembly_id[0:assembly_id.find('_', 4)]

            genome_file = os.path.join(gpath, assembly_id + '_genomic.fna')
            if os.path.exists(genome_file):
                if os.stat(genome_file).st_size == 0:
                    self.logger.warning(
                        'Genomic file appears to be empty: %s' % genome_file)
                    continue

                genome_files.append(genome_file)
            else:
                self.logger.warning(
                    'Genomic file appears to be missing: %s' % genome_file)

            gff_file = os.path.join(gpath, assembly_id + '_genomic.gff')
            if os.path.exists(gff_file):
                if os.stat(gff_file).st_size == 0:
                    self.logger.warning(
                        'GFF appears to be empty: %s' % gff_file)
                    continue

                tt = self._parse_ncbi_translation_table(gff_file)
                if tt:
                    translation_table[canonical_gid] = tt
                else:
                    translation_table[canonical_gid] = None
                    self.logger.warning(
                        'Unable to determine translation table for: %s' % gff_file)

            else:
                self.logger.warning('GFF appears to be missing: %s' % gff_file)

        # run Prodigal on each genome
        self.logger.info('Running Prodigal on %d genomes.' % len(genome_paths))
        prodigal = Prodigal(cpus=self.cpus)
        summary_stats = prodigal.run(genome_files,
                                     translation_table=translation_table,
                                     output_dir=self.tmp_dir)

        # move results into individual genome directories
        self.logger.info('Moving files and calculating checksums.')
        for genome_file in genome_files:
            genome_path, genome_id = ntpath.split(genome_file)
            genome_id = remove_extension(genome_id)
            canonical_gid = genome_id[0:genome_id.find('_', 4)]

            aa_gene_file = os.path.join(self.tmp_dir, genome_id + '_genes.faa')
            nt_gene_file = os.path.join(self.tmp_dir, genome_id + '_genes.fna')
            gff_file = os.path.join(self.tmp_dir, genome_id + '.gff')

            genome_root = genome_id[0:genome_id.find('_', 4)]
            prodigal_path = os.path.join(genome_path, 'prodigal')
            if not os.path.exists(prodigal_path):
                os.makedirs(prodigal_path)
            new_aa_gene_file = os.path.join(
                prodigal_path, genome_root + '_protein.faa')
            new_nt_gene_file = os.path.join(
                prodigal_path, genome_root + '_protein.fna')
            new_gff_file = os.path.join(
                prodigal_path, genome_root + '_protein.gff')

            os.system('mv %s %s' % (aa_gene_file, new_aa_gene_file))
            os.system('mv %s %s' % (nt_gene_file, new_nt_gene_file))
            os.system('mv %s %s' % (gff_file, new_gff_file))

            # save translation table information
            translation_table_file = os.path.join(
                prodigal_path, 'prodigal_translation_table.tsv')
            fout = open(translation_table_file, 'w')
            if canonical_gid in translation_table:
                fout.write('%s\t%d\t%s\n' % ('best_translation_table',
                                             summary_stats[genome_id].best_translation_table,
                                             'used table specified by NCBI'))
            else:
                fout.write('%s\t%d\n' % ('best_translation_table',
                                         summary_stats[genome_id].best_translation_table))
                fout.write('%s\t%.2f\n' % ('coding_density_4',
                                           summary_stats[genome_id].coding_density_4 * 100))
                fout.write('%s\t%.2f\n' % ('coding_density_11',
                                           summary_stats[genome_id].coding_density_11 * 100))
            fout.close()

            checksum = sha256(new_aa_gene_file)
            fout = open(new_aa_gene_file + '.sha256', 'w')
            fout.write(checksum)
            fout.close()

    def _parse_ncbi_translation_table(self, gene_feature_file):
        """Parse translation table from NCBI GFF file."""

        ncbi_transl_table = None
        if os.path.exists(gene_feature_file):
            for line in open(gene_feature_file):
                if line[0] == '#':
                    continue

                if 'transl_table=' in line:
                    trans_num = line[line.rfind('=') + 1:].strip()
                    ncbi_transl_table = int(trans_num)
                    break

        return ncbi_transl_table

    def _parse_prodigal_translation_table(self, prodigal_file):
        """Parse translation table summary generated by run_prodigal function."""
        with open(prodigal_file) as f:
            line = f.readline()
            line_split = line.strip().split()
            return int(line_split[1])
        return None
