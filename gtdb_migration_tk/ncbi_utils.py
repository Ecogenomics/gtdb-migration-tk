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

"""
ncbi_utils.py -- read the NCBI assembly summary files, and hold what the NCBI
commands know in common about NCBI's files.

An assembly summary file (assembly_summary.txt) is the table NCBI publishes
describing every assembly it holds: a block of '#'-prefixed comments, a
'#assembly_accession ...' header, then one tab-separated row per genome.

Two parts of this package read those tables. ncbi_genome_sync.py reads ftp_path to
mirror the genomes from the NCBI FTP site, and select_genomes.py reads
version_status, gbrs_paired_asm, and excluded_from_refseq to decide which
genomes belong in a GTDB release. Both need the same thing from the file, so
the reading lives here and the two callers differ only in what they do with a
row.

Beside the reader sit the few facts about NCBI's files that more than one
command depends on: the accession prefixes of the two databases, the shape of
a line of md5checksums.txt, and the test for an assembly NCBI lists but does
not serve. They live here rather than in whichever command first needed them
so that no command module imports another. ncbi_genome_sync.py imports from
this module and from nothing else in the package, and this module imports
nothing from the package at all; keep both true, or the sync stops being
runnable as the standalone script it started as.

A summary file may be gzipped or not, and open_summary() takes either, so no
caller has to know which it was handed.

Columns are located BY NAME from the header, never by position: NCBI has grown
assembly_summary.txt from 23 fields to 38, and a positional reader silently
starts reading the wrong column the next time the table is revised. Downstream
that mirrors the wrong files, or builds a release from the wrong set of
genomes, with nothing to indicate anything went wrong. For the same reason a
header is required rather than assumed.
"""

import collections
import gzip
import re
from typing import Dict, Iterator, List, Optional, Sequence, Tuple


class BadInput(ValueError):
    """Rejected input: a URL outside the mirror, or a table without the required header."""


# Column names that identify the header row of an assembly summary file.
#
# A summary file is read for either the accession or the FTP path, and the
# outputs written by ncbi_genome_sync carry a subset of the NCBI columns, so a row is
# recognised as the header when it names either. Which columns a caller in fact
# requires is its own business, and is stated through `required`.
HEADER_FIELDS = ('assembly_accession', 'ftp_path')


def open_summary(assembly_summary: str, encoding: str = 'utf-8'):
    """Open an assembly summary file for reading, gzipped or not.

    GTDB stores these files gzipped -- the GenBank bacteria summary alone is over
    a gigabyte of highly repetitive text -- while NCBI serves them, and older
    releases hold them, uncompressed. Every reader of an assembly summary goes
    through here so that neither form has to be thought about anywhere else.

    The encoding is stated rather than left to the locale: organism names carry
    non-ASCII characters, and under a C locale the default would be ASCII and the
    read would fail partway through a perfectly good file.

    Parameters
    ----------
    assembly_summary : str
        NCBI assembly summary file, optionally gzipped.
    encoding : str
        Text encoding of the file.

    @return: open text handle.
    """

    if assembly_summary.endswith('.gz'):
        return gzip.open(assembly_summary, 'rt', encoding=encoding)

    return open(assembly_summary, encoding=encoding)


def summary_columns(line: str,
                    required: Sequence[str] = ()) -> Optional[Dict[str, int]]:
    """Map column name to index for an assembly summary header row.

    The comment block above the header ("##  See ftp://...README_assembly_summary.txt")
    is '#'-prefixed too, so the header is identified by content rather than by
    position: it is the '#' line naming one of the HEADER_FIELDS columns. A
    genome row names neither and is likewise not mistaken for a header.

    The '#' is stripped before the surrounding whitespace, not after, because
    NCBI has written the marker both ways: '#assembly_accession' in the summary
    files it publishes today and '# assembly_accession' in those it published
    until 2020. Stripping in the other order leaves the older spelling as
    ' assembly_accession', which matches no column name, so the header was
    recognised through its ftp_path column and then rejected as missing the
    accession column it in fact carries.

    Parameters
    ----------
    line : str
        Line from an NCBI assembly summary file.
    required : sequence
        Column names the caller cannot proceed without.

    @return: dict of column name to index, or None if this is not a header row.

    Raises
    ------
    BadInput
        The line is a header, but is missing a required column. Returning None
        would leave the caller reading the rest of the file as headerless.
    """

    fields = [field.strip('#').strip() for field in line.split('\t')]
    if not any(name in fields for name in HEADER_FIELDS):
        return None

    missing = [name for name in required if name not in fields]
    if missing:
        raise BadInput('header is missing the {} column(s): {}'.format(
            ', '.join(missing), line))

    return {name: index for index, name in enumerate(fields)}


def summary_field(fields: List[str], columns: Dict[str, int], name: str) -> str:
    """Read a single named field from a row of an assembly summary file.

    Parameters
    ----------
    fields : list
        Tab-separated values of one row.
    columns : dict
        Column name to index, as returned by summary_columns().
    name : str
        Name of the column to read.

    @return: value of the column, or an empty string if the column is absent
        from the header or the row stops short of it.
    """

    index = columns.get(name, -1)

    return fields[index].strip() if 0 <= index < len(fields) else ''


def read_summary_rows(assembly_summary: str,
                      required: Sequence[str] = ()) -> Iterator[Tuple[int, List[str], Dict[str, int]]]:
    """Read the genome rows of an NCBI assembly summary file.

    Comment lines and blank lines are passed over, and the column map in force
    is yielded with each row so fields can be read by name. Summary files carry
    a single header, but one is re-read whenever it appears so that
    concatenated tables are handled.

    Parameters
    ----------
    assembly_summary : str
        NCBI assembly summary file.
    required : sequence
        Column names the caller cannot proceed without.

    @return: iterator over (line number, fields, column map), one per genome.

    Raises
    ------
    BadInput
        Genomes appear before any header, leaving no way to tell which field is
        which. Guessing at column positions is how a summary revision corrupts
        a mirror, or a release, silently.
    """

    columns = None
    with open_summary(assembly_summary) as summary_file:
        for line_number, line in enumerate(summary_file, 1):
            line = line.rstrip('\n').rstrip('\r')

            if not line.strip():
                continue

            if line.startswith('#'):
                columns = summary_columns(line, required) or columns
                continue

            if columns is None:
                raise BadInput(
                    '{}:{}: genomes appear before any "#assembly_accession ..." header; '
                    'this is not an NCBI assembly summary file'.format(
                        assembly_summary, line_number))

            yield line_number, line.split('\t'), columns


def count_summary_rows(assembly_summary: str) -> int:
    """Count the genome rows of an assembly summary file.

    This counts what read_summary_rows() yields, and so must skip exactly what it skips:
    blank lines, and the '#' block of comments and the header. Counting LINES instead --
    which is what a general purpose line counter does -- overstates the total by the size
    of that block, so a progress bar sized from it stops two short of its own total and
    never reaches 100%. The genomes were all read; only the total was wrong.

    It costs one pass over the file, as counting lines did, and the caller is about to
    make a second, far more expensive one.

    Parameters
    ----------
    assembly_summary : str
        NCBI assembly summary file, optionally gzipped.

    @return: number of genome rows in the file.
    """

    with open_summary(assembly_summary) as summary_file:
        return sum(1 for line in summary_file
                   if line.strip() and not line.startswith('#'))


def read_assembly_summary(assembly_summary: str,
                          *field_names: str) -> Iterator[Tuple[str, ...]]:
    """Read the requested fields of each genome in an NCBI assembly summary file.

    Columns absent from the header yield an empty string, as older summary
    files do not carry every field in use today.

    Parameters
    ----------
    assembly_summary : str
        NCBI assembly summary file.
    field_names : str
        Names of the columns to read, e.g. 'assembly_accession'.

    @return: iterator over the requested field values, one tuple per genome.
    """

    for _, fields, columns in read_summary_rows(assembly_summary,
                                                required=('assembly_accession',)):
        yield tuple(summary_field(fields, columns, name) for name in field_names)


# Where NCBI serves everything this toolkit fetches: the assembly summaries and
# the taxonomy (ncbi_metadata_sync) and the genomes themselves (ncbi_genome_sync).
NCBI_HOST = 'ftp.ncbi.nlm.nih.gov'
NCBI_URL = 'https://' + NCBI_HOST

# The two NCBI databases GTDB draws on, and the three names each goes by: the
# directory NCBI serves it from (genomes/<name>/...), the label it is written
# with in prose and logs, and the prefix of its accessions. RefSeq comes first,
# and every consumer keeps that order: select_genomes must know which genomes
# RefSeq covers before it can judge a GenBank assembly.
NCBIDatabase = collections.namedtuple('NCBIDatabase', 'name label prefix')
REFSEQ = NCBIDatabase('refseq', 'RefSeq', 'GCF')
GENBANK = NCBIDatabase('genbank', 'GenBank', 'GCA')
NCBI_DATABASES = (REFSEQ, GENBANK)

# The prefixes on their own, as most callers want them: select_genomes reads a
# genome's database from its accession, update_genomes runs once per prefix.
REFSEQ_PREFIX = REFSEQ.prefix
GENBANK_PREFIX = GENBANK.prefix

# NCBI calls every assembly summary file assembly_summary.txt and tells them
# apart by directory, so GTDB puts the domain and database back into the name as
# it saves them. ncbi_metadata_sync writes these names and select_genomes reads
# the database back out of them, so both go through the two functions below
# rather than each spelling the convention for itself.
ASSEMBLY_SUMMARY_NAME = 'assembly_summary_{domain}_{database}.txt'

# One line of NCBI's md5checksums.txt: the MD5, whitespace, the file name. The
# sync reads the manifest to verify what it fetched, and update_genomes reads
# the mirror's and the previous release's copies to tell whether a genomic FASTA
# changed, so the two must parse the same lines.
MD5_LINE_RE = re.compile(r"^([0-9a-f]{32})\s+(.+)$")


# NCBI's null: what an empty field holds in an assembly summary, and what GTDB
# writes in the same position of the tables it derives from one.
NCBI_NA = 'na'

# The columns of an assembly summary that identify a genome and say whether and
# where NCBI serves it, in this order. They are what ncbi_genome_sync needs to
# mirror a genome, so they open every table GTDB hands it: the selection
# select_genomes writes, and the .fail and .bad files the sync writes for
# itself. Each of those builds its header from here, so the three tables cannot
# drift apart.
GENOME_COLUMNS = ('assembly_accession', 'ftp_path', 'version_status', 'excluded_from_refseq')


def table_header(*columns: str) -> str:
    """The header line of a GTDB table that the assembly summary readers read.

    It is '#'-prefixed so those readers, which skip comment lines, find the
    column names on it the way they find them on an NCBI summary, while the
    line still names its columns for anyone opening the file.

    Parameters
    ----------
    columns : str
        Column names, in order.

    @return: the header line, without a newline.
    """

    return '#' + '\t'.join(columns)


def has_ftp_path(ftp_path: str) -> bool:
    """Report whether NCBI serves a directory for an assembly.

    NCBI writes 'na' in ftp_path for an assembly it lists but does not serve, and the
    column can be empty in an older file. Such a genome cannot be mirrored, so it is not
    selected: the table select_genomes writes is the list ncbi_genome_sync fetches from,
    and a row with nothing to fetch would be reported as skipped by every run of it
    forever.

    select_genomes applies this when choosing a genome, and
    ncbi_genome_sync.read_assembly_summary() when reading the selection back to decide
    a row has "no usable ftp_path". The two must agree, or the selection would promise
    genomes the sync then refuses; sharing the one function makes them agree by
    construction rather than by keeping two copies of the expression in step.

    Parameters
    ----------
    ftp_path : str
        Value of the ftp_path column of an assembly summary file.

    @return: True if the assembly has a directory at NCBI.
    """

    return bool(ftp_path) and ftp_path.lower() != NCBI_NA


def assembly_summary_filename(domain: str, database: NCBIDatabase) -> str:
    """The name GTDB saves one NCBI assembly summary under.

    The name ends in .gz because ncbi_metadata_sync compresses the files as they
    are downloaded; every reader goes through open_summary(), which takes either
    form, and assembly_summary_database() reads the name with or without it.

    Parameters
    ----------
    domain : str
        NCBI directory the file describes: archaea, bacteria or fungi.
    database : NCBIDatabase
        Database the file describes.

    @return: file name, e.g. assembly_summary_bacteria_refseq.txt.gz.
    """

    return ASSEMBLY_SUMMARY_NAME.format(domain=domain, database=database.name) + '.gz'


def assembly_summary_database(filename: str) -> Optional[NCBIDatabase]:
    """The NCBI database an assembly summary file describes, read from its name.

    Only the suffix is read, so a file from an older release, held uncompressed
    or under a longer name, is placed the same way as one this toolkit wrote.

    Parameters
    ----------
    filename : str
        Path or name of an assembly summary file, gzipped or not.

    @return: the database, or None if the name says neither.
    """

    name = filename.rsplit('/', 1)[-1]
    if name.endswith('.gz'):
        name = name[:-len('.gz')]

    # the part of the template after the domain, e.g. _refseq.txt
    suffix = ASSEMBLY_SUMMARY_NAME.split('{domain}')[1]
    for database in NCBI_DATABASES:
        if name.endswith(suffix.format(database=database.name)):
            return database

    return None
