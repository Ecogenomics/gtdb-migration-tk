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

"""Which genomes of a release NCBI marks as a MAG, a SAG or an environmental genome.

The command is parse_ncbi_genome_category, named ncbi_genome_category until
0.1.62; the table and its column keep that name, which is the database field's.
As parse_ncbi_dir does, it writes a table of a fixed name into -o/--output_dir,
gzipped, which update_metadata_db --input_folder knows and loads against
metadata_ncbi_genome_category.desc.tsv; it was written to whatever file -o named,
and loaded by hand with --metadata_table.

The assembly summaries say so for most of them, in excluded_from_refseq
('derived from metagenome', 'derived from single cell', 'derived from
environmental source'): 814,000 of r237's 1,346,119 genomes. That is decided in
the parent, from the summaries the release was selected from (-n, as
select_genomes and ncbi_strains take them, by column name). The command took one
GenBank and one RefSeq summary until 0.1.62, where r237's are four, split by
domain, so two of them left the other domain's genomes to the GenBank files.

Every other genome is looked for in its GenBank file, on worker processes, and
only in the first record's header and source feature
(ncbi_utils.gbff_first_record_head()): in r226 that search added 104 SAGs to the
470,000 genomes the summaries categorised, every one found on an
/isolation_source or /note of the source feature or a reference's TITLE, and
over 5,103 r237 genomes, those 104 among them, the table is the one the whole
file gave. The whole file was read, 3.3 MB of gzip a genome over those r237
genomes: on 16
processes the first record's head is read at 157 genomes a second, the file
server's pace for one machine, and the whole file at 88, about an hour for
r237's 532,119 genomes where it was nearer two.

The hours went elsewhere. Each genome's work was sent with the summaries'
category of every genome of the release, a 16 MB dictionary over r237, pickled
by the queue's one feeder thread (0.13 s) and unpickled by a worker (0.19 s)
for each of 1.35M genomes: some 49 hours of pickling that no number of --cpus
shortened, where r226's smaller dictionary took 8 to 13 hours a part. A worker is now handed a genome_dirs line alone, and only for the
genomes the summaries leave.

The table is written under a temporary name and given its own only once every
genome is decided, so a failed run leaves no table, or the one an earlier run
wrote. A worker that fails stops the run, where it was caught by a bare except,
or ended that worker alone, and the table was written without the genome and
the run exited 0. An excluded_from_refseq value naming a 'derived from' this
command does not know is refused before any GenBank file is read, every such
value named: a new kind of genome is for a person to place.
"""

import gzip
import logging
import multiprocessing as mp
import os
import re
from collections import Counter

from tqdm import tqdm

from gtdb_migration_tk.ncbi_utils import gbff_first_record_head, genomic_gbff, read_assembly_summary
from gtdb_migration_tk.utils.common import GZIP_SUFFIX, open_gzip_text, remove_uncompressed

# the table written in --output_dir, gzipped (with GZIP_SUFFIX), one row a
# genome with a category; update_metadata_db --input_folder loads it under this
# name against metadata_ncbi_genome_category.desc.tsv, which names
# ncbi_genome_category alone, so source is not loaded
CATEGORY_TABLE_NAME = 'ncbi_genome_category.tsv'
CATEGORY_HEADER = ('genome_id', 'ncbi_genome_category', 'source')

# beside the table: the GenBank file lines a category was found on. It is not a
# table, and update_metadata_db passes over it (NOT_TABLES)
EVIDENCE_NAME = 'ncbi_genome_category_evidence.tsv'
EVIDENCE_HEADER = ('genome_id', 'evidence')

MAG, SAG, ENV = 'MAG', 'SAG', 'ENV'

# what the table says of each, in order of precedence: a genome found to be
# both a MAG and a SAG is a MAG, and warned of
CATEGORY = ((MAG, 'derived from metagenome'),
            (SAG, 'derived from single cell'),
            (ENV, 'derived from environmental sample'))

# excluded_from_refseq: the phrase giving each category, in the order looked
# for, and the 'derived from' that is none
SUMMARY_PHRASES = (('derived from single cell', SAG),
                   ('derived from metagenome', MAG),
                   ('derived from environmental source', ENV))
NOT_A_CATEGORY = ('derived from surveillance project',)

# the source column: where a genome's category was found
SOURCE_SUMMARY = 'assembly report'
SOURCE_GBFF = 'GBFF file'
SOURCE_USER = 'user genome'

# a genome named so is a user's, taken to be a MAG
USER_GENOME_PREFIX = 'U_'

METAGENOME_RE = re.compile(r'derived from(\w|\s)*metagenome', re.I)

# genomes handed to a worker process at once
CHUNK = 16

# how many genomes or values a warning or error names
EXAMPLES = 10


class GenomeCategoryError(ValueError):
    """The summaries hold an excluded_from_refseq value this command cannot place."""


def summary_category(excluded):
    """The category an excluded_from_refseq value gives a genome.

    Parameters
    ----------
    excluded : str
        The genome's excluded_from_refseq, e.g. 'derived from metagenome; partial'.

    @return: MAG, SAG or ENV; None where it gives none.

    Raises
    ------
    GenomeCategoryError
        It names a 'derived from' this command does not know.
    """

    for phrase, category in SUMMARY_PHRASES:
        if phrase in excluded:
            return category
    if 'derived' in excluded and not any(phrase in excluded for phrase in NOT_A_CATEGORY):
        raise GenomeCategoryError(excluded)
    return None


def gbff_categories(lines):
    """The categories the head of a GenBank file's first record gives a genome.

    Parameters
    ----------
    lines : iterable of str
        The lines of the file, as an open file gives them; only the first
        record's head is read (ncbi_utils.gbff_first_record_head()).

    @return: (set of MAG, SAG and ENV found, [each line one was found on, stripped]).
    """

    categories = set()
    evidence = []
    for line in gbff_first_record_head(lines):
        if 'metagenome' in line:
            if '/metagenome_source' in line or METAGENOME_RE.search(line):
                categories.add(MAG)
                evidence.append(line.strip())
        elif 'single cell' in line:
            categories.add(SAG)
            evidence.append(line.strip())
        elif '/environmental_sample' in line or 'derived from environmental source' in line:
            categories.add(ENV)
            evidence.append(line.strip())

    return categories, evidence


def read_gbff(line):
    """The categories a genome's GenBank file gives it; run on a worker process.

    The file is read with any byte that is not UTF-8 replaced, as nothing looked
    for is anything but ASCII: one stray byte in a header stops the run otherwise.

    Parameters
    ----------
    line : str
        The genome's line of the genome_dirs file.

    @return: (accession, set of categories, [evidence lines], whether it has a GenBank file).
    """

    genome_id, genome_dir = line.rstrip('\n').split('\t')[:2]
    try:
        with open_gbff(genomic_gbff(genome_dir)) as handle:
            categories, evidence = gbff_categories(handle)
    except FileNotFoundError:
        return genome_id, set(), [], False

    return genome_id, categories, evidence, True


def open_gbff(path):
    """A gzipped GenBank file opened as text, any byte that is not UTF-8 replaced.

    @return: the open file.
    """

    return gzip.open(path, 'rt', encoding='utf-8', errors='replace')


class GenomeType(object):
    """Identify genomes marked by NCBI as being a MAG, a SAG or an environmental genome."""

    def __init__(self, cpus=1):
        """Initialization.

        Parameters
        ----------
        cpus : int
            Worker processes reading the GenBank files.

        @return: None
        """

        self.logger = logging.getLogger('timestamp')
        self.cpus = cpus

    def summary_categories(self, assembly_summary_files, release):
        """The category the summaries give each genome of the release, and which they hold.

        Parameters
        ----------
        assembly_summary_files : sequence of str
            NCBI assembly summaries, RefSeq and GenBank alike, gzipped or not.
        release : set of str
            The accessions of the release.

        @return: (accession -> MAG, SAG or ENV, for the genomes given one;
                 set of the release's accessions any summary holds).

        Raises
        ------
        GenomeCategoryError
            A value names a 'derived from' this command does not know: every one
            is named, with how many of the release's genomes hold it.
        """

        categories = {}
        in_summaries = set()
        unknown = Counter()
        for assembly_summary in assembly_summary_files:
            read = 0
            for accession, excluded in read_assembly_summary(assembly_summary, 'assembly_accession',
                                                             'excluded_from_refseq'):
                if accession not in release:
                    continue
                read += 1
                in_summaries.add(accession)
                try:
                    category = summary_category(excluded)
                except GenomeCategoryError:
                    unknown[excluded] += 1
                    continue
                if category is not None:
                    categories[accession] = category
            self.logger.info('Read excluded_from_refseq for {:,} genomes of the release from {}.'.format(
                read, assembly_summary))

        if unknown:
            raise GenomeCategoryError(
                'The assembly summaries hold {:,} excluded_from_refseq value(s) naming a "derived from" this '
                'command does not know, which it cannot place as a MAG, SAG or environmental genome: {}. '
                'Add each to SUMMARY_PHRASES or NOT_A_CATEGORY in ncbi_genome_category.py.'.format(
                    len(unknown), '; '.join('"{}" ({:,} genomes)'.format(value, count)
                                            for value, count in unknown.most_common())))

        return categories, in_summaries

    def run(self, assembly_summary_files, genome_file, output_dir):
        """Write the NCBI genome category of each genome of a release that has one.

        Parameters
        ----------
        assembly_summary_files : sequence of str
            The NCBI assembly summaries the release was selected from (-n).
        genome_file : str
            genome_dirs file: accession, directory, canonical accession.
        output_dir : str
            Directory CATEGORY_TABLE_NAME is written to, gzipped (with
            GZIP_SUFFIX): CATEGORY_HEADER, a row for each genome with a
            category, in accession order. The GenBank file lines each category
            was found on are written beside it, EVIDENCE_NAME.

        @return: the path of the table written.

        Raises
        ------
        GenomeCategoryError
            As summary_categories(); nothing is written.
        """

        with open(genome_file) as genomes:
            lines = [line for line in genomes if line.strip()]
        release = {line.split('\t', 1)[0] for line in lines}
        self.logger.info('Deciding the NCBI genome category of {:,} genomes.'.format(len(release)))

        categories, in_summaries = self.summary_categories(assembly_summary_files, release)
        for category, described in CATEGORY:
            self.logger.info(' - the assembly summaries mark {:,} genomes {}.'.format(
                sum(1 for c in categories.values() if c == category), described))

        rows = {}
        for genome_id, category in categories.items():
            rows[genome_id] = (category, SOURCE_SUMMARY)

        to_read = []
        for line in lines:
            genome_id = line.split('\t', 1)[0]
            if genome_id.startswith(USER_GENOME_PREFIX):
                rows[genome_id] = (MAG, SOURCE_USER)
            elif genome_id not in categories:
                to_read.append(line)

        no_summary = sorted(release - in_summaries - {g for g in release if g.startswith(USER_GENOME_PREFIX)})
        self.logger.info('Reading the GenBank files of the {:,} genomes the summaries do not categorise '
                         'on {:,} processes.'.format(len(to_read), self.cpus))

        evidence_rows = []
        no_gbff = []
        both = []
        if to_read:
            with mp.Pool(processes=max(1, self.cpus)) as pool:
                for genome_id, found, evidence, has_gbff in tqdm(
                        pool.imap_unordered(read_gbff, to_read, chunksize=CHUNK),
                        total=len(to_read), ncols=100, smoothing=50 / len(to_read), unit='genome'):
                    if not has_gbff:
                        no_gbff.append(genome_id)
                    if MAG in found and SAG in found:
                        both.append(genome_id)
                    for category, _ in CATEGORY:
                        if category in found:
                            rows[genome_id] = (category, SOURCE_GBFF)
                            break
                    evidence_rows.extend((genome_id, line) for line in evidence)

        from_gbff = Counter(category for category, source in rows.values() if source == SOURCE_GBFF)
        for category, described in CATEGORY:
            self.logger.info(' - the GenBank files mark {:,} more {}.'.format(from_gbff[category], described))

        described = dict(CATEGORY)
        output_file = os.path.join(output_dir, CATEGORY_TABLE_NAME + GZIP_SUFFIX)
        self.write(output_file, CATEGORY_HEADER,
                   ((genome_id, described[category], source)
                    for genome_id, (category, source) in sorted(rows.items())))
        remove_uncompressed(output_file, self.logger)
        self.write(os.path.join(output_dir, EVIDENCE_NAME), EVIDENCE_HEADER, sorted(evidence_rows))

        if both:
            self.logger.warning('Identified {:,} genomes whose GenBank file marks them both a MAG and a SAG, '
                                'e.g.: {}; each is written as a MAG.'.format(
                                    len(both), ', '.join(sorted(both)[:EXAMPLES])))
        if no_gbff:
            self.logger.warning('Identified {:,} genomes with a missing _genomic.gbff.gz file, e.g.: {}; '
                                'they are written with no category.'.format(
                                    len(no_gbff), ', '.join(sorted(no_gbff)[:EXAMPLES])))
        if no_summary:
            self.logger.warning('Identified {:,} genomes in none of the assembly summaries, e.g.: {}; only '
                                'their GenBank files decided them. The summaries are not those the release '
                                'was selected from, or one is missing.'.format(
                                    len(no_summary), ', '.join(no_summary[:EXAMPLES])))
        self.logger.info('Wrote the NCBI genome category of {:,} genomes to {}.'.format(len(rows), output_file))

        return output_file

    @staticmethod
    def write(path, header, rows):
        """Write a table under a temporary name, then give it its own.

        Parameters
        ----------
        path : str
            The table, gzipped where it ends with GZIP_SUFFIX.
        header : tuple of str
            Its columns.
        rows : iterable of tuple of str
            Its rows.

        @return: None
        """

        partial = path + '.partial'
        try:
            with (open_gzip_text(partial) if path.endswith(GZIP_SUFFIX) else open(partial, 'w')) as handle:
                handle.write('\t'.join(header) + '\n')
                for row in rows:
                    handle.write('\t'.join(row) + '\n')
        except BaseException:
            if os.path.exists(partial):
                os.remove(partial)
            raise
        os.replace(partial, path)
