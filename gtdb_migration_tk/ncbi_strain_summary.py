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
with. The NCBI type material status of each genome is the assembly summaries'
relation_to_type_material, read by column name from the summaries the release
was selected from (-n), as select_genomes and strains type_table take them; the
command took them as four files, RefSeq and GenBank of bacteria and archaea
(--rb, --ra, --gb, --ga), until 0.1.55. The GenBank files are not read: the
/strain= and /isolate= qualifiers of
<assembly>_genomic.gbff and <assembly>_wgsmaster.gbff were looked for under
those names, which NCBI serves gzipped, so they were never found. Over r237 the
summaries give every strain ID the database holds for all but 51 of 1,047,596
genomes.
"""

import os
import logging
import multiprocessing as mp
import re

from tqdm import tqdm

from gtdb_migration_tk.biolib_lite.common import get_num_lines
from gtdb_migration_tk.ncbi_utils import NCBI_NA, read_assembly_summary, strain_identifiers
from gtdb_migration_tk.utils.common import GZIP_SUFFIX, open_gzip_text, remove_uncompressed

# the table written in --out_dir, which update_metadata_db loads under this name
STRAIN_SUMMARY_NAME = 'strain_summary_file.tsv'
STRAIN_SUMMARY_HEADER = ('genome_id', 'Organism name', 'ncbi_strain_identifiers',
                         'ncbi_type_material_designation')

# the NCBI file of a genome directory the strain IDs are read from
ASSEMBLY_REPORT_SUFFIX = '_assembly_report.txt'

# genomes handed to a worker process at once
CHUNK = 16

# how many genomes a warning names
EXAMPLES = 10


def read_assembly_report(line):
    """The organism name and strain IDs of one genome, from its assembly report.

    Run on a worker process. The type material status is not looked up here: the
    summaries' table of every NCBI genome stays in the parent rather than being
    handed to each worker.

    Parameters
    ----------
    line : str
        The genome's line of the genome_dirs file.

    @return: (accession, organism name, strain IDs, its infraspecific name where
             it names a substrain or else None, whether it has an assembly report).
    """

    line_split = line.strip().split('\t')
    genome_id = line_split[0]
    genome_path = line_split[1]
    genome_dir_id = os.path.basename(os.path.normpath(genome_path))

    species = ''
    infraspecific_names = []
    isolate = ''
    report_file = os.path.join(genome_path, genome_dir_id + ASSEMBLY_REPORT_SUFFIX)
    has_report = os.path.exists(report_file)
    if has_report:
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
    substrain = infraspecific_name if 'substr.' in infraspecific_name else None
    strain_ids = strain_identifiers(infraspecific_name, isolate or NCBI_NA)

    return genome_id, species, strain_ids, substrain, has_report

class NCBIStrainParser(object):
    """Extract genes in nucleotide space."""

    def __init__(self, assembly_summary_files, cpus):
        """Initialization.

        Parameters
        ----------
        assembly_summary_files : sequence of str
            The NCBI assembly summaries the release was selected from, gzipped
            or not, RefSeq and GenBank alike: an accession names its database.
        cpus : int
            Worker processes reading the assembly reports.

        @return: None
        """

        self.logger = logging.getLogger('timestamp')
        self.type_material = self.parse_summaries(assembly_summary_files)
        self.cpus = cpus

    def parse_summaries(self, assembly_summary_files):
        """The NCBI type material status of every genome of the summaries.

        Parameters
        ----------
        assembly_summary_files : sequence of str
            NCBI assembly summaries, read by column name
            (ncbi_utils.read_assembly_summary()).

        @return: accession -> relation_to_type_material, e.g.
                 {'GCF_000005845.2': 'assembly from type material', ...}.
        """

        type_material = {}
        for assembly_summary in assembly_summary_files:
            before = len(type_material)
            for accession, relation in read_assembly_summary(assembly_summary, 'assembly_accession',
                                                             'relation_to_type_material'):
                type_material[accession] = relation
            self.logger.info('Read the NCBI type material status of {:,} genomes from {}.'.format(
                len(type_material) - before, assembly_summary))
        return type_material

    def generate_ncbi_strains_summary(self, genome_dir_file, out_dir):
        """Write the strain IDs and NCBI type material status of every genome.

        The assembly reports are read on self.cpus processes and each row is
        written as it comes back, in no particular order, rather than held until
        every genome is done. A worker that fails stops the run: the table is
        then incomplete, and a run that went on would report success. The table
        is written under a temporary name and given its own only once every
        genome is in it, so a failed run leaves no table, or the one an earlier
        run wrote, rather than part of one for update_metadata_db to load. A genome
        in none of the summaries has no type material status, written empty,
        which update_metadata_db loads as NULL; it was written as the text None.
        One with no assembly report has no organism name or strain IDs. Each is
        counted in a closing warning.

        Parameters
        ----------
        genome_dir_file : str
            genome_dirs file: accession, directory, canonical accession.
        out_dir : str
            Directory STRAIN_SUMMARY_NAME is written to, gzipped (with GZIP_SUFFIX).

        @return: the path of the table written.
        """

        output_file = os.path.join(out_dir, STRAIN_SUMMARY_NAME + GZIP_SUFFIX)
        genome_count = get_num_lines(genome_dir_file)
        self.logger.info('Reading the assembly reports of {:,} genomes on {:,} processes into {}.'.format(
            genome_count, self.cpus, output_file))

        no_report = []
        no_status = []
        substrains = []
        written = 0
        partial = output_file + '.partial'
        try:
            with open_gzip_text(partial) as outf, open(genome_dir_file) as genomes:
                outf.write('\t'.join(STRAIN_SUMMARY_HEADER) + '\n')
                with mp.Pool(processes=max(1, self.cpus)) as pool:
                    for genome_id, species, strain_ids, substrain, has_report in tqdm(
                            pool.imap_unordered(read_assembly_report, genomes, chunksize=CHUNK),
                            total=genome_count, ncols=100, smoothing=50 / max(genome_count, 1), unit='genome'):
                        typemat = self.type_material.get(genome_id)
                        if typemat is None:
                            no_status.append(genome_id)
                        if not has_report:
                            no_report.append(genome_id)
                        if substrain is not None:
                            substrains.append('{} strain {}'.format(genome_id, substrain))
                        outf.write('{}\t{}\t{}\t{}\n'.format(
                            genome_id, species, ';'.join(strain_ids), typemat if typemat is not None else ''))
                        written += 1
        except BaseException:
            if os.path.exists(partial):
                os.remove(partial)
            raise
        os.replace(partial, output_file)
        remove_uncompressed(output_file, self.logger)

        self.logger.info("Substrains of interest:")
        for substr in sorted(substrains):
            self.logger.info('- ' + substr)
        if no_report:
            self.logger.warning('Identified {:,} genomes with a missing {} file, e.g.: {}; their organism '
                                'name and strain IDs are empty.'.format(
                                    len(no_report), ASSEMBLY_REPORT_SUFFIX,
                                    ', '.join(sorted(no_report)[:EXAMPLES])))
        if no_status:
            self.logger.warning('Identified {:,} genomes in none of the assembly summaries, e.g.: {}; their '
                                'ncbi_type_material_designation is empty. The summaries are not those '
                                'the release was selected from, or one is missing.'.format(
                                    len(no_status), ', '.join(sorted(no_status)[:EXAMPLES])))
        self.logger.info('Wrote the strain IDs and NCBI type material status of {:,} genomes to {}.'.format(
            written, output_file))

        return output_file
