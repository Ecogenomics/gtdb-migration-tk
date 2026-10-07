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

__prog_name__ = 'ncbi_assembly_file_metadata.py'
__prog_desc__ = 'Produce filtered metadata file from NCBI assembly metadata file.'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2015'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.0.3'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import logging
import os
import multiprocessing as mp
import random
import re
import string
from collections import defaultdict

from numpy import (zeros as np_zeros,sum as np_sum)
from tqdm import tqdm

from gtdb_migration_tk.biolib_lite.common import get_num_lines
from gtdb_migration_tk.ncbi_utils import read_assembly_summary
from gtdb_migration_tk.utils.common import GZIP_SUFFIX, open_gzip_text, remove_uncompressed
from gtdb_migration_tk.utils.tools import openfile


class GenericFeatureParser():
    """Parses generic feature file (GFF)."""

    def __init__(self, filename):
        self.cds_count = 0
        self.tRNA_count = 0
        self.rRNA_count = 0
        self.rRNA_16S_count = 0
        self.ncRNA_count = 0

        self.genes = {}
        self.last_coding_base = {}

        self._parse(filename)

        self.coding_mask = {}
        for seq_id in self.genes:
            self.coding_mask[seq_id] = self._coding_mask(seq_id)

    def _parse(self, gff_file):
        """Parse GFF file.

        Parameters
        ----------
        gff_file : str
          Generic feature file to parse.
        """

        with openfile(gff_file) as handle:
            for line in handle:
                if line[0] == '#':
                    continue

                line_split = line.split('\t')
                if line_split[2] == 'tRNA':
                    self.tRNA_count += 1
                elif line_split[2] == 'rRNA':
                    self.rRNA_count += 1

                    if 'product=16S ribosomal RNA' in line_split[8]:
                        self.rRNA_16S_count += 1
                elif line_split[2] == 'ncRNA':
                    self.ncRNA_count += 1
                elif line_split[2] == 'CDS':
                    self.cds_count += 1

                    seq_id = line_split[0]
                    if seq_id not in self.genes:
                        self.genes[seq_id] = []
                        self.last_coding_base[seq_id] = 0

                    start = int(line_split[3])
                    end = int(line_split[4])

                    # a list of intervals, not gene IDs of this module's making: the
                    # counter those came from was reset only when a contig was first
                    # met, so a GFF returning to an earlier contig reused IDs that
                    # contig already had and overwrote its own genes. Nothing read them
                    self.genes[seq_id].append([start, end])
                    self.last_coding_base[seq_id] = max(
                        self.last_coding_base[seq_id], end)

    def _coding_mask(self, seq_id):
        """Build mask indicating which bases in a sequences are coding."""

        # safe way to calculate coding bases as it accounts
        # for the potential of overlapping genes
        # last_coding_base + 1 because a GFF counts bases from 1: without the
        # extra entry the final base of the rightmost gene fell off the end
        coding_mask = np_zeros(self.last_coding_base[seq_id] + 1)
        for pos in self.genes[seq_id]:
            coding_mask[pos[0]:pos[1] + 1] = 1

        return coding_mask

    def coding_bases(self, seq_id):
        """Calculate number of coding bases in sequence."""

        # check if sequence has any genes
        if seq_id not in self.genes:
            return 0

        return np_sum(self.coding_mask[seq_id])

    def total_coding_bases(self):
        """Calculate total number of coding bases.

        Returns
        -------
        int
          Number of coding bases.
        """

        coding_bases = 0
        for seq_id in self.genes:
            coding_bases += self.coding_bases(seq_id)

        return int(coding_bases)


# The table parse_ncbi_dir writes in its --out_dir. update_metadata_db loads a
# table of this name against metadata_ncbi_assembly.desc.tsv; under any other
# name --input_folder refuses it. It is not the table parse_ncbi_assemblies
# writes (NCBI_ASSEMBLY_TABLE, ncbi_assembly_summary.tsv).
NCBI_DIR_TABLE = 'ncbi_assembly_metadata.tsv'

# The NCBI files of a genome directory parse_ncbi_dir reads, by the suffix of
# the assembly's name each is found under. A genome missing one is given a row
# all the same, empty in the fields read from it, which update_metadata_db loads
# as NULL; the run ends by saying how many genomes were missing each.
ASSEMBLY_STATS_SUFFIX = '_assembly_stats.txt'
GFF_SUFFIX = '_genomic.gff.gz'
GBFF_SUFFIX = '_genomic.gbff.gz'
NCBI_DIR_FILES = (ASSEMBLY_STATS_SUFFIX, GFF_SUFFIX, GBFF_SUFFIX)

# NCBI publishes no annotation, and so no GFF, for many GenBank assemblies (121
# of 300 genomes of r237 drawn at random): a genome without one is told of, not
# warned of. Every assembly has its statistics and its GenBank flat file.
NCBI_DIR_UNANNOTATED = frozenset({GFF_SUFFIX})

# genomes handed to a worker process at once
NCBI_DIR_CHUNK = 16

# how many genomes a warning names
EXAMPLES = 10


class NCBIMetaDir(object):
    """Create metadata file from the assembly stats file of each NCBI assembly."""

    def __init__(self,cpus=1):
        self.fields = ['Assembly name', 'Organism name',
                       'Taxid', 'Submitter', 'Date']
        self.fields.extend(['BioSample', 'Assembly type',
                            'Release type', 'Assembly level'])
        self.fields.extend(['Genome representation', 'GenBank assembly accession',
                            'RefSeq assembly and GenBank assemblies identical'])

        self.stats = ['molecule-count', 'contig-count',
                      'scaffold-count', 'region-count', 'top-level-count']
        self.stats.extend(['spanned-gaps', 'total-gap-length',
                           'total-length', 'ungapped-length', 'unspanned-gaps'])
        self.stats.extend(['contig-N50', 'scaffold-N50', 'contig-L50', 'component-count',
                           'scaffold-L50', 'scaffold-N75', 'scaffold-N90'])

        self.gff_fields = ['cds_count', 'tRNA_count',
                           'ncRNA_count', 'rRNA_count', 'ssu_count']

        self.gbff_fields = ['translation_table', 'isolation_source', 'geo_loc_name', 'lat_lon','metagenome_source']
        self.stats_info = {}

        self.cpus = cpus
        self.logger = logging.getLogger('timestamp')

    def _randomword(self, length):
        """Generate a random string of lowercase letters to mask internal slashes."""
        return ''.join(random.choice(string.ascii_lowercase) for i in range(length))

    def _parse_assembly_stats(self, assembly_stat_file):
        """Parse data from assembly stats file.

        Parameters
        ----------
        assembly_stat_file : str
          NCBI assembly stat file to parse.

        Returns
        -------
        list
          Parsed metadata in canonical order.
        """

        metadata_fields = [''] * len(self.fields)
        metadata_stats = [''] * len(self.stats)

        file_section = 'Assembly info'
        with open(assembly_stat_file) as handle:
            for line in handle:
                if 'Assembly Statistics Report' in line:
                    file_section = 'ASR'
                elif 'Statistic Types' in line:
                    file_section = 'ST'
                elif 'Sequence-type Description' in line:
                    file_section = 'SD'

                if file_section == 'ASR' and ':' in line:
                    field = line[2:line.find(':')]
                    value = line[line.find(':') + 1:].strip()
                    if field in self.fields:
                        if field == 'Organism name' and '(' in value:
                            metadata_index = self.fields.index(field)
                            metadata_fields[metadata_index] = value[0:value.find(
                                '(')].strip()
                        else:
                            metadata_index = self.fields.index(field)
                            metadata_fields[metadata_index] = value
                elif file_section == 'ST':
                    line_split = line.split('\t')
                    if len(line_split) == 2:
                        field = line_split[0][2:]
                        desc = line_split[1].strip()
                        if field in self.stats:
                            self.stats_info[field] = desc
                elif file_section == 'SD':
                    line_split = line.split('\t')
                    if len(line_split) == 6 and (line_split[0] in ['all', 'Primary Assembly']) and (line_split[1] == 'all') and (line_split[3] == 'all'):
                        field = line_split[4]
                        value = line_split[5].strip()
                        if field in self.stats:
                            metadata_index = self.stats.index(field)
                            metadata_stats[metadata_index] = value

        return metadata_fields, metadata_stats

    def _parse_gff(self, gff_file):
        """Parse statistics from generic feature file (GFF)."""

        metadata_gff = [''] * len(self.gff_fields)

        gff_parser = GenericFeatureParser(gff_file)
        metadata_gff[self.gff_fields.index('cds_count')] = gff_parser.cds_count
        metadata_gff[self.gff_fields.index(
            'tRNA_count')] = gff_parser.tRNA_count
        metadata_gff[self.gff_fields.index(
            'ncRNA_count')] = gff_parser.ncRNA_count
        metadata_gff[self.gff_fields.index(
            'rRNA_count')] = gff_parser.rRNA_count
        metadata_gff[self.gff_fields.index(
            'ssu_count')] = gff_parser.rRNA_16S_count

        return metadata_gff

    def _parse_gbff(self, genbank_file):
        """Parse statistics and metadata from GenBank file."""
        metadata_gbff = [''] * len(self.gbff_fields)

        pattern_gene = re.compile(r"^\s{0,20}\w")
        pattern_source = re.compile(r"^\s{5}source\s{10}")
        source_info_bool = False
        randomstring = self._randomword(10)
        source_info = []

        # We read the file line by line using openfile to handle gzipped files
        with openfile(genbank_file) as handle:
            for line in handle:

                # 1. Extract Translation Table
                if '/transl_table=' in line:
                    translation_table = line[line.rfind('=') + 1:].strip()
                    metadata_gbff[self.gbff_fields.index('translation_table')] = translation_table

                # 2. Extract Source Metadata
                if pattern_source.match(line):
                    source_info_bool = True
                elif pattern_gene.match(line) and source_info_bool:
                    # Stop appending to source_info once we hit the next main feature
                    source_info_bool = False
                elif source_info_bool:
                    # Replace all '/' characters by a random string except the first one
                    # '/' will be used to separate those metadata later on
                    line = re.sub(r"(?!^\/)\/", randomstring, ' '.join(line.split()))
                    source_info.append("{0} ".format(' '.join(line.split())))

        # Process the source_info array if any metadata was found
        if source_info:
            source_info_string = ''.join(x for x in source_info)
            source_info_array = source_info_string.split("/")
            source_info_dict = {}

            for info in source_info_array:
                if "=" in info:
                    try:
                        k, v = info.split("=", 1)
                        source_info_dict[k] = v
                    except Exception as e:
                        print(info)

            # Map the extracted dictionary back to the tracked gbff_fields
            for field in ['isolation_source', 'geo_loc_name', 'lat_lon','metagenome_source']:
                if field in source_info_dict:
                    # Replace special characters and restore the masked '/' strings
                    clean_val = source_info_dict[field].replace('"', '').replace(',', ';').replace("'", " ").replace(
                        randomstring, "/").rstrip()
                    metadata_gbff[self.gbff_fields.index(field)] = clean_val

        return metadata_gbff

    def header(self):
        """The header of NCBI_DIR_TABLE, its columns in the order the worker writes them.

        @return: the header line, newline ended.
        """

        return '\t'.join(['Assembly accession']
                         + ['ncbi_' + x.lower().replace(' ', '_') for x in self.fields]
                         + ['ncbi_' + x.lower().replace('-', '_') for x in self.stats]
                         + ['ncbi_' + x.lower() for x in self.gff_fields]
                         + ['ncbi_' + x.lower() for x in self.gbff_fields]) + '\n'

    def parse_ncbi_dir(self, gtdb_genome_path_file, output_dir):
        """Write the metadata NCBI's own files in each genome directory hold.

        Every genome of the genome_dirs file is given a row. One missing a file
        is empty in the fields read from it, which update_metadata_db loads as
        NULL, and the run ends by warning how many genomes were missing each
        file. Rows are written as the workers finish them, in no particular
        order, rather than held until every genome is done.

        Parameters
        ----------
        gtdb_genome_path_file : str
            genome_dirs file: accession, directory, canonical accession.
        output_dir : str
            Directory NCBI_DIR_TABLE is written to, gzipped (with GZIP_SUFFIX).

        @return: the path of the table written.
        """

        output_file = os.path.join(output_dir, NCBI_DIR_TABLE + GZIP_SUFFIX)
        genome_count = get_num_lines(gtdb_genome_path_file)
        self.logger.info('Reading the NCBI files of {:,} genomes on {:,} processes into {}.'.format(
            genome_count, self.cpus, output_file))

        missing = defaultdict(list)
        written = 0
        with open_gzip_text(output_file) as fout, open(gtdb_genome_path_file) as genomes:
            fout.write(self.header())
            with mp.Pool(processes=self.cpus) as pool:
                for gid, line_to_write, absent in tqdm(
                        pool.imap_unordered(self.ncbi_parser_worker, genomes, chunksize=NCBI_DIR_CHUNK),
                        total=genome_count, ncols=100, smoothing=50 / max(genome_count, 1), unit='genome'):
                    fout.write(line_to_write + '\n')
                    written += 1
                    for suffix in absent:
                        missing[suffix].append(gid)

        for suffix in NCBI_DIR_FILES:
            if missing[suffix]:
                gids = sorted(missing[suffix])
                message = 'Identified {:,} genomes with a missing {} file, e.g.: {}; their fields read ' \
                          'from it are empty.'.format(len(gids), suffix, ', '.join(gids[:EXAMPLES]))
                if suffix in NCBI_DIR_UNANNOTATED:
                    self.logger.info(message + ' NCBI does not annotate every assembly.')
                else:
                    self.logger.warning(message)
        self.logger.info('Wrote the NCBI metadata of {:,} genomes to {}.'.format(written, output_file))

        remove_uncompressed(output_file, self.logger)
        return output_file

    def ncbi_parser_worker(self, line):
        """Read the NCBI files of one genome; run on a worker process.

        Parameters
        ----------
        line : str
            The genome's line of the genome_dirs file.

        @return: the genome's accession, its row (without a newline), and the
                 suffixes of NCBI_DIR_FILES it has no file for.
        """

        line_split = line.strip().split('\t')
        gid = line_split[0]
        gpath = line_split[1]
        assembly_id = os.path.basename(os.path.normpath(gpath))

        absent = []

        assembly_stat_file = os.path.join(gpath, assembly_id + ASSEMBLY_STATS_SUFFIX)
        if os.path.exists(assembly_stat_file):
            metadata_fields, metadata_stats = self._parse_assembly_stats(assembly_stat_file)
        else:
            metadata_fields, metadata_stats = [''] * len(self.fields), [''] * len(self.stats)
            absent.append(ASSEMBLY_STATS_SUFFIX)

        gff_file = os.path.join(gpath, assembly_id + GFF_SUFFIX)
        if os.path.exists(gff_file):
            gff_stats = self._parse_gff(gff_file)
        else:
            gff_stats = [''] * len(self.gff_fields)
            absent.append(GFF_SUFFIX)

        genbank_file = os.path.join(gpath, assembly_id + GBFF_SUFFIX)
        if os.path.exists(genbank_file):
            gbff_stats = self._parse_gbff(genbank_file)
        else:
            gbff_stats = [''] * len(self.gbff_fields)
            absent.append(GBFF_SUFFIX)

        line_to_write = '\t'.join([gid] + metadata_fields + metadata_stats
                                  + [str(v) for v in gff_stats] + [str(v) for v in gbff_stats])

        return gid, line_to_write, absent



# The table parse_ncbi_assemblies writes in its --out_dir. update_metadata_db loads a
# table of this name against metadata_ncbi_assembly_file.desc.tsv, whose fields
# are the ones written here; under any other name --input_folder refuses it, and
# under ncbi_assembly_metadata.tsv it would be read against the description of
# parse_ncbi_dir's fields.
NCBI_ASSEMBLY_TABLE = 'ncbi_assembly_summary.tsv'

class NCBIMeta(object):
    """Create metadata file from the assembly stats file of each NCBI assembly."""

    # the columns of an assembly summary read, in the order their fields are written
    COLUMNS = ('bioproject', 'wgs_master', 'refseq_category', 'species_taxid', 'isolate',
               'version_status', 'seq_rel_date', 'asm_name', 'gbrs_paired_asm',
               'paired_asm_comp', 'excluded_from_refseq', 'relation_to_type_material')

    def __init__(self):
        self.logger = logging.getLogger('timestamp')

        self.fields = {'bioproject': ['ncbi_bioproject'],
                       'wgs_master': ['ncbi_wgs_master', 'ncbi_wgs_formatted'],
                       'refseq_category': ['ncbi_refseq_category'],
                       'species_taxid': ['ncbi_species_taxid'],
                       'isolate': ['ncbi_isolate'],
                       'version_status': ['ncbi_version_status'],
                       'seq_rel_date': ['ncbi_seq_rel_date'],
                       'asm_name': ['ncbi_asm_name'],
                       'gbrs_paired_asm': ['ncbi_gbrs_paired_asm'],
                       'paired_asm_comp': ['ncbi_paired_asm_comp'],
                       'relation_to_type_material': ['ncbi_type_material_designation'],
                       'excluded_from_refseq': ['ncbi_excluded_from_refseq','ncbi_not_used_as_type']}

    def field_values(self, column, value):
        """The values of the fields written from one column of an assembly summary.

        Parameters
        ----------
        column : str
            Column of the assembly summary, e.g. 'wgs_master'.
        value : str
            The genome's value of the column.

        @return: one value for each of self.fields[column].
        """

        if column == 'wgs_master':
            return [value, self.format_wgs(value)]
        if column == 'excluded_from_refseq':
            return ['' if value == 'na' else value, str('not used as type' in value)]
        return [value]

    def parse_assemblies(self, assembly_summary_files, output_dir):
        """Create metadata by parsing NCBI assembly metadata files.

        Every genome of the summaries is written: which of them are loaded into
        the database is update_metadata_db's to decide. A genome is written
        once, from the first summary holding it: a genome in two summaries, or
        a summary given twice, would otherwise have two rows.

        Parameters
        ----------
        assembly_summary_files : sequence of str
            The NCBI assembly summaries the release was selected from, gzipped or
            not, read by column name (ncbi_utils.read_assembly_summary()).
        output_dir : str
            Directory NCBI_ASSEMBLY_TABLE is written to, gzipped (with GZIP_SUFFIX), a row
            per genome.

        @return: the path of the table written.
        """

        output_file = os.path.join(output_dir, NCBI_ASSEMBLY_TABLE + GZIP_SUFFIX)
        self.logger.info('Writing the NCBI metadata of every genome of {:,} assembly '
                         'summaries to {}.'.format(len(assembly_summary_files), output_file))

        # write out metadata
        written = set()
        again = []
        with open_gzip_text(output_file) as fout:
            fout.write('\t'.join(['genome_id'] + [field for column in self.COLUMNS
                                                  for field in self.fields[column]]) + '\n')
            for assembly_file in assembly_summary_files:
                self.logger.info('Reading {}.'.format(assembly_file))
                rows = 0
                before = len(written)
                for row in read_assembly_summary(assembly_file, 'assembly_accession', *self.COLUMNS):
                    rows += 1
                    genome_id = row[0]
                    if genome_id.startswith('GCA_'):
                        genome_id = 'GB_' + genome_id
                    elif genome_id.startswith('GCF_'):
                        genome_id = 'RS_' + genome_id

                    if genome_id in written:
                        again.append((genome_id, assembly_file))
                        continue
                    values = [genome_id]
                    for column, value in zip(self.COLUMNS, row[1:]):
                        values.extend(self.field_values(column, value))
                    fout.write('\t'.join(values) + '\n')
                    written.add(genome_id)
                self.logger.info('  Wrote {:,} of its {:,} genome(s).'.format(len(written) - before, rows))

        if again:
            self.logger.warning('{:,} genome(s) are in more than one assembly summary, e.g. {}; '
                                'each was written once, from the first summary holding it. A '
                                'summary was given twice, or the summaries overlap.'.format(
                                    len(again), ', '.join('{} again in {}'.format(gid, path)
                                                          for gid, path in again[:EXAMPLES])))
        self.logger.info('Wrote the NCBI metadata of {:,} genome(s) from {:,} assembly summaries '
                         'to {}.'.format(len(written), len(assembly_summary_files), output_file))

        remove_uncompressed(output_file, self.logger)
        return output_file

    def format_wgs(self, wgs_accession):
        if not wgs_accession or wgs_accession == "na":
            return ""
        wgs_acc, version = wgs_accession.split('.')
        idx = [ch.isdigit() for ch in wgs_acc].index(True)
        wgs_id = wgs_acc[0:idx] + str(version).zfill(2)
        return wgs_id
