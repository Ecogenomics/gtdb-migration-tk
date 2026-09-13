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
ncbi_utils.py -- read the NCBI assembly summary files.

An assembly summary file (assembly_summary.txt) is the table NCBI publishes
describing every assembly it holds: a block of '#'-prefixed comments, a
'#assembly_accession ...' header, then one tab-separated row per genome.

Two parts of this package read those tables. ncbi_sync.py reads ftp_path to
mirror the genomes from the NCBI FTP site, and ftp_manager.py reads
version_status, gbrs_paired_asm, and excluded_from_refseq to decide which of
the mirrored genomes belong in a GTDB release. Both need the same thing from
the file, so the reading lives here and the two callers differ only in what
they do with a row.

Columns are located BY NAME from the header, never by position: NCBI has grown
assembly_summary.txt from 23 fields to 38, and a positional reader silently
starts reading the wrong column the next time the table is revised. Downstream
that mirrors the wrong files, or builds a release from the wrong set of
genomes, with nothing to indicate anything went wrong. For the same reason a
header is required rather than assumed.
"""

import os
from typing import Dict, Iterator, List, Optional, Sequence, Tuple


class BadInput(ValueError):
    """Rejected input: a URL outside the mirror, or a table without the required header."""


# Column names that identify the header row of an assembly summary file.
#
# A summary file is read for either the accession or the FTP path, and the
# outputs written by ncbi_sync carry a subset of the NCBI columns, so a row is
# recognised as the header when it names either. Which columns a caller in fact
# requires is its own business, and is stated through `required`.
HEADER_FIELDS = ('assembly_accession', 'ftp_path')


def summary_columns(line: str,
                    required: Sequence[str] = ()) -> Optional[Dict[str, int]]:
    """Map column name to index for an assembly summary header row.

    The comment block above the header ("##  See ftp://...README_assembly_summary.txt")
    is '#'-prefixed too, so the header is identified by content rather than by
    position: it is the '#' line naming one of the HEADER_FIELDS columns. A
    genome row names neither and is likewise not mistaken for a header.

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

    fields = [field.strip().lstrip('#') for field in line.split('\t')]
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
    with open(assembly_summary) as summary_file:
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


def genome_assembly_file(genome_dir: str) -> str:
    """Path of the genomic FASTA file expected within a genome directory.

    NCBI names every file of an assembly after the directory holding it, e.g.
    GCF_002287175.1_ASM228717v1/GCF_002287175.1_ASM228717v1_genomic.fna.gz, and
    the GTDB directories mirror the NCBI layout. The name of the genomic FASTA
    file therefore follows from the directory alone, without consulting the
    assembly summary file. Files are gzipped once they are part of GTDB.

    Parameters
    ----------
    genome_dir : str
        Genome directory, named for the assembly it holds.

    @return: path of the genomic FASTA file the directory should contain.
    """

    return os.path.join(genome_dir, os.path.basename(genome_dir) + '_genomic.fna.gz')
