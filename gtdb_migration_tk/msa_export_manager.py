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

"""Export the concatenated marker alignment of a set of genomes: export_msa and list_marker_sets.

export_msa writes, for every genome of --genome_batchfile, its aligned markers of
the sets of --marker_set_ids, concatenated, with the genome's ARB metadata: the
concatenation and the files of 'gtdb tree create --no_tree --no_alignment',
under the same names, and the columns of each marker as partitions. list_marker_sets is 'gtdb marker_sets view --all', the sets
--marker_set_ids is chosen from.

THE FULL MSA OF EVERY GENOME ASKED FOR

gtdb filtered the genomes (CheckM quality, exclusion, the share of the MSA they
have residues in, with genomes guaranteed against all three) and trimmed the
columns (a random sample of those passing gap and consensus filters). Neither is
done here: the MSA holds every genome of the batchfile and every column of every
marker. Which genomes and which columns a tree is inferred from is decided
downstream, from the full MSA, where it can be decided again without asking the
database again; and so there is no mask, no list of filtered genomes, and no
seed.

IT ONLY READS

'gtdb tree create' first ran assignToRepresentative(), giving every genome
without a gtdb_representative one by AAI, and aligned every marker of every
genome given that aligned_markers had no row for, both committed, so that making
a tree changed the database. gtdb_representative is the toolkit's to write
(update_propagated_tax, update_reps_db), as aligned_markers is (align_marker_genes).
The export is one read-only REPEATABLE READ transaction, so the genomes and their
alignments are of one snapshot.

A GENOME NOT ALIGNED IS REFUSED

A genome without a row (with a sequence) in aligned_markers for each of the
markers is refused, before anything is written: they are listed beside the log
(UNALIGNED_NAME) and the run ends naming align_marker_genes. gtdb aligned them
itself; a genome with some of the rows ended 'gtdb tree create' with a KeyError,
and one with none was dropped from the files with a WARNING. A genome the
database does not hold, and a marker set it does not hold, are refused the same
way (UNKNOWN_GENOMES_NAME).

THE ARB METADATA, AND ITS TABLE

The ARB metadata is written without the alignment (aligned_seq= is left blank),
as --no_alignment wrote it, but msa_aa_count and msa_length are those of the
genome's sequence in <prefix>_concatenated.faa.gz; gtdb counted them from the
blanked sequence, giving every genome 0. <prefix>_metadata.tsv.gz holds the same
fields, in the same order, one row a genome under a header, so that the metadata
can be read without an ARB parser: db_name, every column of metadata_view but id
with ncbi_genus and ncbi_species after ncbi_taxonomy, gtdb_cluster_ncbi_species,
then the five msa_ fields. Its values are the ARB record's, NCBI's organism name
and the canonical clustered genomes included, except that a float is written to
three decimal places (the ARB record: '%.4g'), '=' is kept (the ARB record: '/',
which ARB's parser needs), and NULL is 'none', as export_metadata writes it.

THE PARTITIONS OF THE MSA

gene_partitions.txt names the columns of each marker in <prefix>_concatenated.faa.gz,
in RAxML's partition format ('AA, PFAM_PF00380.20 = 1-121'), for locusmasker
(--partitions gene_partitions.txt --datatype aa) or a partitioned tree. It can be
written because the MSA is untrimmed: each marker is a run of columns as long as
its HMM. Each partition is marked AA, amino acids, rather than given a
substitution model, which is chosen where a tree is inferred; locusmasker reads
the name and the range and ignores what precedes them. A marker is named as
<prefix>_multi_hits.tsv names it, which locusmasker requires to be unique.

THE ORDER IS FIXED

Every file lists the genomes ordered by accession, where gtdb gave them in the
order of a set. The alignments and both metadata files are gzipped, as gtdb's
were not, their headers recording no time, so that two exports of an unchanged
database are the same files.

READING aligned_markers IN PARALLEL

A genome's rows are read through aligned_markers' primary key, ALIGNED_CHUNK
genomes a query. Over r237 a chunk of 1,000 bacterial representatives took 13 s
from disk and 1 s from the database's cache: about 13,500 random reads of 8 kB
blocks at a millisecond each, the rows of a genome's markers being spread over
the 99 GB table by when each set was aligned. Neither COPY nor ordering the
genomes by id made a difference, and Python's share was 0.2 s. So the time is
the storage's latency, and is cut by asking it more than one thing at a time:
-c/--cpus connections each read a chunk. Over 5,000 to 20,000 genomes not read
before, 1 connection read 80 genomes a second, 8 478, 16 629 and 32 1,247, the
storage not yet saturated. A sequential scan of the table was not taken: it
reads at about 99 MB/s, 18 minutes for any export.
Each connection joins the main transaction's snapshot (pg_export_snapshot(),
SET TRANSACTION SNAPSHOT), so the export is still of one snapshot, and runs
without JIT compilation, which cost 0.1 s a chunk and saved nothing.

READING metadata_view

metadata_view joins representative_stats_view by accession, which is made for
every genome of the database whatever the query asks for: 100 genomes asked for
by id took 96 s over r237, 350 by accession 230 s, and leaving out the view's
two cluster columns still 64 s. So the view is read once, in full, through a
server-side cursor, and the rows of the genomes exported kept: about six minutes
over r237, however few the genomes.
"""

import contextlib
import logging
import os
import queue
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from decimal import Decimal
from typing import Dict, Iterable, List, NamedTuple, Optional, Sequence, Set, Tuple

from psycopg2.extensions import ISOLATION_LEVEL_REPEATABLE_READ

from gtdb_migration_tk.biolib_lite.common import canonical_gid
from gtdb_migration_tk.biolib_lite.logger import log_directory
from gtdb_migration_tk.database_configuration.GenomeDatabaseConnectionFTPUpdate import \
    GenomeDatabaseConnectionFTPUpdate
from gtdb_migration_tk.utils.common import GZIP_SUFFIX, open_gzip_text
from gtdb_migration_tk.utils.prettytable import PrettyTable

# the alignment's gap
GAP = '-'

# what became of a marker in a genome, as <prefix>_multi_hits.tsv says it
SINGLE = 'Single'
MULTI_IDENTICAL = 'Multi_Identical'
MULTIPLE = 'Multiple'
MISSING = 'Missing'

# the fields of the ARB record after the metadata, and what is added to the metadata
MSA_FIELDS = ('msa_gene_count', 'msa_num_marker_genes', 'msa_aa_count', 'msa_length', 'multiple_homologs')
CLUSTER_SPECIES_FIELD = 'gtdb_cluster_ncbi_species'
NCBI_RANK_FIELDS = ('ncbi_genus', 'ncbi_species')
DB_NAME = 'db_name'
NULL = 'none'

# the files written, each <out_dir>/<prefix><suffix>
DEFAULT_PREFIX = 'gtdb'
MULTI_HITS_SUFFIX = '_multi_hits.tsv'
MSA_SUFFIX = '_concatenated.faa.gz'
ARB_METADATA_SUFFIX = '_arb_metadata.txt.gz'
ARB_FILTER_SUFFIX = '_arb_filter.ift'
TSV_METADATA_SUFFIX = '_metadata.tsv.gz'
# with --individual, <out_dir>/<prefix>_<marker's id in its database><suffix>
INDIVIDUAL_SUFFIX = '.faa.gz'
MARKERS_INFO_SUFFIX = '_markers_info.tsv'
# the columns of each marker, in <out_dir>, whatever the prefix
PARTITIONS_NAME = 'gene_partitions.txt'
# what a partition is marked: amino acids, its model left to the tree's inference
PARTITION_MODEL = 'AA'

# beside the log, where a run is refused
UNKNOWN_GENOMES_NAME = 'export_msa_unknown_genomes.tsv'
UNALIGNED_NAME = 'export_msa_unaligned_genomes.tsv'

# what an output is written to until it is complete
PARTIAL_SUFFIX = '.partial'

# genomes a query of aligned_markers asks for, through its primary key
ALIGNED_CHUNK = 1000
# how a connection reading aligned_markers joins the main transaction
EXPORT_SNAPSHOT = 'SELECT pg_export_snapshot()'
SET_SNAPSHOT = 'SET TRANSACTION SNAPSHOT %s'
NO_JIT = 'SET jit = off'
# rows of metadata_view fetched at a time
VIEW_ITERSIZE = 10000
# how many genomes a line of the log names
EXAMPLES = 10

SOURCE_PREFIXES = 'SELECT external_id_prefix FROM genome_sources'
GENOMES_BY_ID_AT_SOURCE = (
    'SELECT g.id, gs.external_id_prefix, g.id_at_source, g.formatted_source_id '
    'FROM genomes g JOIN genome_sources gs ON gs.id = g.genome_source_id '
    'WHERE g.id_at_source = ANY(%s)')
MARKER_SETS_BY_ID = 'SELECT id FROM marker_sets WHERE id = ANY(%s)'
# gtdb's order: by database, then the marker's ID in it, as PostgreSQL sorts them
MARKERS_OF_SETS = (
    'SELECT DISTINCT m.id, m.name, m.description, m.id_in_database, m.size, md.external_id_prefix '
    'FROM markers m JOIN marker_databases md ON md.id = m.marker_database_id '
    'JOIN marker_set_contents c ON c.marker_id = m.id '
    'WHERE c.set_id = ANY(%s) '
    'ORDER BY md.external_id_prefix ASC, m.id_in_database ASC')
ALIGNED_ROWS = (
    'SELECT genome_id, marker_id, sequence, multiple_hits, evalue, unique_genes '
    'FROM aligned_markers WHERE genome_id = ANY(%s) AND marker_id = ANY(%s) AND sequence IS NOT NULL')
METADATA_VIEW = 'SELECT * FROM metadata_view'
NCBI_TAXONOMY_BY_ID_AT_SOURCE = (
    "SELECT gs.external_id_prefix || '_' || g.id_at_source, m.ncbi_taxonomy "
    'FROM genomes g JOIN genome_sources gs ON gs.id = g.genome_source_id '
    'LEFT JOIN metadata_taxonomy m ON m.id = g.id '
    'WHERE g.id_at_source = ANY(%s)')
LIST_MARKER_SETS = (
    'SELECT s.id, s.name, s.description, s.owned_by_root, u.username, count(c.marker_id) '
    'FROM marker_sets s LEFT JOIN users u ON u.id = s.owner_id '
    'LEFT JOIN marker_set_contents c ON c.set_id = s.id '
    'GROUP BY s.id, u.username ORDER BY s.id')
MARKER_SET_HEADER = ('set_id', 'name', 'description', 'owner', 'marker_count')


class MSAExportError(ValueError):
    """What the run was given cannot be exported; nothing has been written."""


class Marker(NamedTuple):
    """A marker of the MSA, as gtdb named it (external_id: e.g. PFAM_PF00380.20)."""

    db_id: int
    external_id: str
    id_in_database: str
    name: str
    description: str
    size: int


class AlignedRow(NamedTuple):
    """A genome's row of aligned_markers for one marker."""

    sequence: str
    multiple_hits: bool
    evalue: Optional[str]
    unique_genes: Optional[int]


class GenomeName(NamedTuple):
    """A genome as the files name it: accession (GB_GCA_...) and db_name (G...)."""

    accession: str
    db_name: str


class AlignedGenome(NamedTuple):
    """A genome's markers concatenated, and what became of each marker."""

    sequence: str
    statuses: List[str]


def read_genome_file(path: str) -> List[str]:
    """The genomes of a batchfile, as gtdb read one: the first column, '#' lines skipped.

    Parameters
    ----------
    path : str
        One genome a line; further tab-separated columns are ignored.

    @return: the genomes, in the file's order, each once.
    """

    genomes = []
    with open(path) as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            genome = line.strip().split('\t')[0].strip()
            if genome:
                genomes.append(genome)
    return list(dict.fromkeys(genomes))


def flatten(values: Optional[Iterable]) -> List:
    """The items of an option given as comma-separated lists, perhaps more than once.

    Parameters
    ----------
    values : iterable or None
        Each a list (as the argparse types of __main__ give them) or a string to split on ','.

    @return: the items, in order, each once.
    """

    items = []
    for value in values or []:
        parts = value.split(',') if isinstance(value, str) else value
        for part in parts:
            if isinstance(part, str):
                part = part.strip()
                if not part:
                    continue
            if part not in items:
                items.append(part)
    return items


def source_id(name: str, prefixes: Iterable[str]) -> Tuple[Optional[str], str]:
    """Split a genome as named into its source prefix, where it has one, and genomes.id_at_source.

    Parameters
    ----------
    name : str
        e.g. GB_GCA_000003645.1, or GCA_000003645.1.
    prefixes : iterable of str
        genome_sources.external_id_prefix, e.g. GB, RS.

    @return: e.g. ('GB', 'GCA_000003645.1'), or (None, 'GCA_000003645.1').
    """

    head, _, rest = name.partition('_')
    if rest and head in set(prefixes):
        return head, rest
    return None, name


def resolve_genomes(names: Iterable[str], held: Dict[Tuple[str, str], int],
                    prefixes: Iterable[str]) -> Tuple[Dict[str, int], List[str]]:
    """The database's genome for each name, and the names it holds none for.

    A name without a prefix is the genome of whichever source holds it, and
    unknown where two do.

    Parameters
    ----------
    names : iterable of str
        Genomes as named.
    held : dict
        (external_id_prefix, id_at_source) -> genomes.id.
    prefixes : iterable of str
        genome_sources.external_id_prefix.

    @return: ({name: genomes.id}, [names not held]).
    """

    prefixes = list(prefixes)
    found, unknown = {}, []
    for name in names:
        prefix, at_source = source_id(name, prefixes)
        if prefix is not None:
            ids = [held[(prefix, at_source)]] if (prefix, at_source) in held else []
        else:
            ids = [held[(p, at_source)] for p in prefixes if (p, at_source) in held]
        if len(ids) == 1:
            found[name] = ids[0]
        elif name not in unknown:
            unknown.append(name)
    return found, unknown


def concatenate(rows: Dict[int, AlignedRow], markers: Sequence[Marker]) -> AlignedGenome:
    """A genome's markers concatenated, and what became of each, as gtdb concatenated them.

    A marker no gene was found for (no evalue), or that several different genes
    were, is a run of gaps as long as the marker.

    Parameters
    ----------
    rows : dict
        marker db_id -> the genome's row; one for every marker.
    markers : sequence of Marker
        In the order of the MSA.

    @return: AlignedGenome.

    Raises
    ------
    MSAExportError
        A sequence used is not as long as its marker, and so would shift every
        column after it.
    """

    parts, statuses = [], []
    for marker in markers:
        row = rows[marker.db_id]
        unique_genes = row.unique_genes if row.unique_genes is not None else 0
        present = bool(row.evalue)
        if not present:
            statuses.append(MISSING)
        elif not row.multiple_hits:
            statuses.append(SINGLE)
        elif unique_genes == 1:
            statuses.append(MULTI_IDENTICAL)
        else:
            statuses.append(MULTIPLE)

        if present and (not row.multiple_hits or unique_genes == 1):
            if len(row.sequence) != marker.size:
                raise MSAExportError('The aligned sequence of {} is {:,} long, where the marker is {:,}.'.format(
                    marker.external_id, len(row.sequence), marker.size))
            parts.append(row.sequence)
        else:
            parts.append(GAP * marker.size)
    return AlignedGenome(''.join(parts), statuses)


def gene_partitions(markers: Sequence[Marker], model: str = PARTITION_MODEL) -> str:
    """The columns of each marker of the MSA, in RAxML's partition format.

    Parameters
    ----------
    markers : sequence of Marker
        In the order of the MSA.
    model : str
        Given to every partition.

    @return: e.g. 'AA, PFAM_PF00380.20 = 1-121\nAA, PFAM_PF00410.20 = 122-247\n';
             1-based and inclusive, together every column.

    Raises
    ------
    MSAExportError
        Two markers of one name, or a name RAxML's format cannot hold.
    """

    lines, start, seen = [], 1, set()
    for marker in markers:
        name = marker.external_id
        if name in seen or any(c in name for c in ',=#') or name != name.strip():
            raise MSAExportError('Marker {!r} cannot be named in a partition file: it repeats another, or holds a '
                                 "',', '=' or '#'.".format(name))
        seen.add(name)
        end = start + marker.size - 1
        lines.append('{}, {} = {}-{}\n'.format(model, name, start, end))
        start = end + 1
    return ''.join(lines)


def metadata_fields(view_columns: Sequence[str]) -> List[str]:
    """The fields of the ARB record before the msa_ fields, in order.

    @return: the columns of metadata_view but id, with ncbi_genus and
             ncbi_species after ncbi_taxonomy, then gtdb_cluster_ncbi_species.
    """

    fields = [c for c in view_columns if c != 'id']
    index = fields.index('ncbi_taxonomy') + 1
    fields[index:index] = list(NCBI_RANK_FIELDS)
    fields.append(CLUSTER_SPECIES_FIELD)
    return fields


def ncbi_species_of(ncbi_taxonomy: Optional[str]) -> Optional[str]:
    """The species of an NCBI taxonomy string, its last rank.

    @return: e.g. 's__Escherichia coli', or None where there is no taxonomy.
    """

    if not ncbi_taxonomy:
        return None
    return ncbi_taxonomy.split(';')[-1].strip()


def cluster_ncbi_species(clustered_genomes: Optional[str], species: Dict[str, Optional[str]]) -> str:
    """gtdb_cluster_ncbi_species: the NCBI species of a cluster's genomes, most common first.

    Parameters
    ----------
    clustered_genomes : str or None
        gtdb_clustered_genomes: accessions separated by ';'.
    species : dict
        accession -> NCBI species; a genome not in it is not counted.

    @return: e.g. 's__Escherichia coli: 12; s__Escherichia sp.: 1', or ''.
    """

    counts: Dict[str, int] = {}
    for genome in clustered_genomes.split(';') if clustered_genomes else []:
        sp = species.get(genome)
        if sp is not None:
            counts[sp] = counts.get(sp, 0) + 1
    return '; '.join('{}: {}'.format(sp, n) for sp, n in sorted(counts.items(), key=lambda kv: kv[1], reverse=True))


def genome_values(metadata: Dict[str, object], db_name: str, cluster_species: str) -> Dict[str, object]:
    """A genome's values for metadata_fields(), as the ARB record gives them.

    For an NCBI genome (db_name G...), organism_name is NCBI's, ncbi_genus and
    ncbi_species are taken from ncbi_taxonomy, and gtdb_clustered_genomes is
    written in canonical form.

    @return: {field: value}.
    """

    values = dict(metadata)
    values.pop('id', None)
    genus = species = None
    if db_name.startswith('G'):
        ranks = (values.get('ncbi_taxonomy') or '').split(';')
        if len(ranks) >= 7:
            genus, species = ranks[5].strip(), ranks[6].strip()
        values['organism_name'] = values.get('ncbi_organism_name')
        clustered = values.get('gtdb_clustered_genomes')
        if clustered is not None:
            values['gtdb_clustered_genomes'] = ';'.join(canonical_gid(g) for g in clustered.split(';'))
    values['ncbi_genus'], values['ncbi_species'] = genus, species
    values[CLUSTER_SPECIES_FIELD] = cluster_species
    return values


def msa_values(sequence: str, statuses: Sequence[str]) -> Dict[str, int]:
    """The msa_ fields of a genome, counted from the sequence written.

    @return: {field: value}.
    """

    return {'msa_gene_count': sum(1 for s in statuses if s == SINGLE),
            'msa_num_marker_genes': len(statuses),
            'msa_aa_count': len(sequence) - sequence.count(GAP),
            'msa_length': len(sequence),
            'multiple_homologs': sum(1 for s in statuses if s == MULTIPLE)}


def arb_value(value) -> str:
    """A value as the ARB record writes it: a float '%.4g', '=' as '/', NULL 'None'.

    @return: the text.
    """

    if isinstance(value, float):
        value = '%.4g' % value
    value = str(value)
    if value:
        value = value.replace('=', '/')
    return value


def tsv_value(value) -> str:
    """A value as the TSV writes it: a float to three decimal places, NULL 'none'.

    A tab or line break inside a value is written as a space, so that a row is a line.

    @return: the text.
    """

    if value is None:
        return NULL
    if isinstance(value, (float, Decimal)):
        return '%.3f' % value
    return str(value).replace('\t', ' ').replace('\r', ' ').replace('\n', ' ')


def arb_record(db_name: str, fields: Sequence[str], values: Dict[str, object], msa: Dict[str, int]) -> str:
    """A genome's ARB record, its alignment left blank.

    @return: the record, BEGIN to END and a blank line.
    """

    lines = ['BEGIN', '{}={}'.format(DB_NAME, db_name)]
    lines.extend('{}={}'.format(f, arb_value(values.get(f))) for f in fields)
    lines.extend('{}={}'.format(f, msa[f]) for f in MSA_FIELDS)
    lines.extend(['aligned_seq=', 'END', '', ''])
    return '\n'.join(lines)


def tsv_header(fields: Sequence[str]) -> str:
    """The header of <prefix>_metadata.tsv.gz.

    @return: the line.
    """

    return '\t'.join([DB_NAME] + list(fields) + list(MSA_FIELDS)) + '\n'


def tsv_row(db_name: str, fields: Sequence[str], values: Dict[str, object], msa: Dict[str, int]) -> str:
    """A genome's row of <prefix>_metadata.tsv.gz.

    @return: the line.
    """

    return '\t'.join([db_name] + [tsv_value(values.get(f)) for f in fields]
                     + [str(msa[f]) for f in MSA_FIELDS]) + '\n'


def arb_filter(fields: Sequence[str]) -> str:
    """The ARB import filter of the metadata file, as gtdb wrote it.

    @return: the filter's text.
    """

    out = ['AUTODETECT\t"BEGIN"\n\n', 'BEGIN\t"BEGIN*"\n\n']

    def match(field_name, written):
        out.append('MATCH\t"%s\\=*"\n' % field_name)
        out.append('\tSRT "*\\=="\n')
        out.append('\tWRITE "%s"\n\n' % written)

    match(DB_NAME, 'name')
    # organism name near the top, for convenience
    match('organism_name', 'organism_name')
    for f in list(fields) + list(MSA_FIELDS):
        if f != 'organism_name':
            match(f, f)
    out.append('SEQUENCEAFTER\t"multiple_homologs*"\n')
    out.append('SEQUENCESRT\t"*\\=="\n')
    out.append('SEQUENCEEND\t"END"\n\n')
    out.append('END\t"END"\n')
    return ''.join(out)


@contextlib.contextmanager
def output_file(path: str):
    """A file written beside its name and renamed onto it once complete; gzipped where it ends in .gz.

    A gzipped file's header records no time, so the same text gives the same bytes
    (utils.common.open_gzip_text()).

    @return: a text stream, as a context manager.
    """

    partial = path + PARTIAL_SUFFIX
    try:
        with (open_gzip_text(partial) if path.endswith(GZIP_SUFFIX) else open(partial, 'w')) as handle:
            yield handle
        os.replace(partial, path)
    except BaseException:
        if os.path.exists(partial):
            os.remove(partial)
        raise


def write_beside_log(name: str, header: str, lines: Iterable[str]) -> str:
    """A list of what a run was refused over, beside the log.

    @return: the file's path.
    """

    path = os.path.join(log_directory(), name)
    with open(path, 'w') as handle:
        handle.write(header)
        handle.writelines(lines)
    return path


class MSAExportManager(object):
    """export_msa and list_marker_sets, in one read-only transaction."""

    def __init__(self, database, cpus: int = 1):
        """Connect, in a read-only REPEATABLE READ transaction.

        Parameters
        ----------
        database : dict
            How to reach the database, as utils.common.database_keywords() gives it.
        cpus : int
            Connections reading aligned_markers at once.

        @return: None
        """

        self.logger = logging.getLogger('timestamp')
        self.database = database
        self.cpus = max(1, cpus)
        self.readers: List[GenomeDatabaseConnectionFTPUpdate] = []
        self.reader_cursors: 'queue.Queue' = queue.Queue()
        self.temp_con = GenomeDatabaseConnectionFTPUpdate(database)
        self.temp_con.MakePostgresConnection()
        # set before the first statement opens the transaction
        self.temp_con.conn.set_session(isolation_level=ISOLATION_LEVEL_REPEATABLE_READ, readonly=True)
        self.temp_cur = self.temp_con.cursor()

    # ---- what is read from the database

    def read_marker_sets(self, marker_set_ids: Sequence[int]) -> Set[int]:
        """The marker sets of those given that the database holds.

        @return: {marker_sets.id}.
        """

        self.temp_cur.execute(MARKER_SETS_BY_ID, (list(marker_set_ids),))
        return {row[0] for row in self.temp_cur.fetchall()}

    def read_markers(self, marker_set_ids: Sequence[int]) -> List[Marker]:
        """The markers of the sets, each once, in gtdb's order.

        @return: [Marker].
        """

        self.temp_cur.execute(MARKERS_OF_SETS, (list(marker_set_ids),))
        return [Marker(db_id, '{}_{}'.format(prefix, id_in_database), id_in_database, name, description, size)
                for db_id, name, description, id_in_database, size, prefix in self.temp_cur.fetchall()]

    def read_source_prefixes(self) -> List[str]:
        """genome_sources.external_id_prefix, e.g. ['GB', 'RS'].

        @return: the prefixes.
        """

        self.temp_cur.execute(SOURCE_PREFIXES)
        return [row[0] for row in self.temp_cur.fetchall()]

    def read_genomes(self, at_source: Sequence[str]) -> List[Tuple[int, str, str, str]]:
        """The genomes whose id_at_source is one of those given.

        @return: [(genomes.id, external_id_prefix, id_at_source, formatted_source_id)].
        """

        self.temp_cur.execute(GENOMES_BY_ID_AT_SOURCE, (list(at_source),))
        return list(self.temp_cur.fetchall())

    def open_readers(self) -> None:
        """Open self.cpus connections for aligned_markers, each in the main transaction's snapshot.

        @return: None
        """

        self.temp_cur.execute(EXPORT_SNAPSHOT)
        snapshot = self.temp_cur.fetchone()[0]
        for _ in range(self.cpus):
            reader = GenomeDatabaseConnectionFTPUpdate(self.database)
            reader.MakePostgresConnection()
            self.readers.append(reader)
            reader.conn.set_session(isolation_level=ISOLATION_LEVEL_REPEATABLE_READ, readonly=True)
            cursor = reader.cursor()
            # the first statement of the transaction, as it must be
            cursor.execute(SET_SNAPSHOT, (snapshot,))
            cursor.execute(NO_JIT)
            self.reader_cursors.put(cursor)

    def close_readers(self) -> None:
        """Close the connections open_readers() opened.

        @return: None
        """

        while self.readers:
            reader = self.readers.pop()
            try:
                reader.rollback()
            finally:
                reader.ClosePostgresConnection()

    def read_aligned_rows(self, genome_ids: Sequence[int],
                          marker_ids: Sequence[int]) -> Dict[int, Dict[int, AlignedRow]]:
        """The rows of aligned_markers, with a sequence, of the genomes for the markers.

        Read on whichever of open_readers()' connections is free.

        @return: {genomes.id: {markers.id: AlignedRow}}.
        """

        cursor = self.reader_cursors.get()
        try:
            cursor.execute(ALIGNED_ROWS, (list(genome_ids), list(marker_ids)))
            fetched = cursor.fetchall()
        finally:
            self.reader_cursors.put(cursor)
        rows: Dict[int, Dict[int, AlignedRow]] = {}
        for genome_id, marker_id, sequence, multiple_hits, evalue, unique_genes in fetched:
            rows.setdefault(genome_id, {})[marker_id] = AlignedRow(sequence, multiple_hits, evalue, unique_genes)
        return rows

    def read_metadata(self, genome_ids: Set[int]) -> Tuple[List[str], Dict[int, Dict[str, object]]]:
        """The rows of metadata_view of the genomes, read through the whole view (see READING metadata_view).

        @return: (the view's columns, {genomes.id: {column: value}}).
        """

        cursor = self.temp_con.conn.cursor(name='export_msa_metadata_view')
        cursor.itersize = VIEW_ITERSIZE
        try:
            cursor.execute(METADATA_VIEW)
            columns: Optional[List[str]] = None
            metadata = {}
            for row in cursor:
                if columns is None:
                    columns = [d[0] for d in cursor.description]
                    id_index = columns.index('id')
                if row[id_index] in genome_ids:
                    metadata[row[id_index]] = dict(zip(columns, row))
            if columns is None:
                columns = [d[0] for d in cursor.description]
        finally:
            cursor.close()
        return columns, metadata

    def read_ncbi_taxonomy(self, accessions: Iterable[str], prefixes: Sequence[str]) -> Dict[str, Optional[str]]:
        """The NCBI taxonomy of the genomes named.

        @return: {accession: ncbi_taxonomy}, only for genomes the database holds.
        """

        accessions = set(accessions)
        if not accessions:
            return {}
        self.temp_cur.execute(NCBI_TAXONOMY_BY_ID_AT_SOURCE,
                              (sorted({source_id(a, prefixes)[1] for a in accessions}),))
        return {accession: taxonomy for accession, taxonomy in self.temp_cur.fetchall() if accession in accessions}

    # ---- list_marker_sets

    def list_marker_sets(self) -> List[tuple]:
        """Print every marker set of the database: its ID, name, description, owner and size.

        @return: the rows printed.
        """

        try:
            self.temp_cur.execute(LIST_MARKER_SETS)
            rows = [(set_id, name, description, 'root' if owned_by_root else username, count)
                    for set_id, name, description, owned_by_root, username, count in self.temp_cur.fetchall()]
        finally:
            self.temp_con.rollback()

        table = PrettyTable()
        table.field_names = list(MARKER_SET_HEADER)
        for row in rows:
            table.add_row([NULL if v is None else v for v in row])
        print(table)
        self.logger.info('The database holds {:,} marker sets.'.format(len(rows)))
        return rows

    # ---- export_msa

    def resolve(self, genomes: Sequence[str]) -> Tuple[List[int], Dict[int, GenomeName], List[str]]:
        """Every genome named, as the database holds it, ordered by accession.

        @return: ([genomes.id], {genomes.id: GenomeName}, prefixes).

        Raises
        ------
        MSAExportError
            A genome named that the database does not hold.
        """

        prefixes = self.read_source_prefixes()
        held, names = {}, {}
        for genome_id, prefix, at_source, formatted in self.read_genomes(
                sorted({source_id(n, prefixes)[1] for n in genomes})):
            held[(prefix, at_source)] = genome_id
            names[genome_id] = GenomeName('{}_{}'.format(prefix, at_source), formatted)
        found, unknown = resolve_genomes(genomes, held, prefixes)
        if unknown:
            report = write_beside_log(UNKNOWN_GENOMES_NAME, 'genome_id\n', (n + '\n' for n in unknown))
            raise MSAExportError('{:,} genome(s) named are not in the database, e.g. {}; every one is in {}.'.format(
                len(unknown), ', '.join(unknown[:EXAMPLES]), report))
        # a genome named twice, with and without its prefix, is exported once
        order = sorted(set(found.values()), key=lambda g: names[g].accession)
        return order, names, prefixes

    def align_chunk(self, chunk: Sequence[int], markers: Sequence[Marker],
                    names: Dict[int, GenomeName]) -> Tuple[Dict[int, AlignedGenome], List[Tuple[int, List[Marker]]]]:
        """Read and concatenate the markers of a chunk of genomes; run on a reader's thread.

        @return: ({genomes.id: AlignedGenome}, [(genomes.id, [markers it has no row for])]).

        Raises
        ------
        MSAExportError
            A sequence of the wrong length.
        """

        rows = self.read_aligned_rows(chunk, [m.db_id for m in markers])
        aligned, unaligned = {}, []
        for genome_id in chunk:
            genome_rows = rows.get(genome_id, {})
            lacking = [m for m in markers if m.db_id not in genome_rows]
            if lacking:
                unaligned.append((genome_id, lacking))
                continue
            try:
                aligned[genome_id] = concatenate(genome_rows, markers)
            except MSAExportError as exc:
                raise MSAExportError('{}: {}'.format(names[genome_id].accession, exc)) from None
        return aligned, unaligned

    def align(self, genome_ids: List[int], markers: Sequence[Marker],
              names: Dict[int, GenomeName]) -> Dict[int, AlignedGenome]:
        """Concatenate the markers of each genome, chunks read on self.cpus connections at once.

        Raises
        ------
        MSAExportError
            A genome without a row for each marker (UNALIGNED_NAME), or a sequence of the wrong length.

        @return: {genomes.id: AlignedGenome}.
        """

        aligned: Dict[int, AlignedGenome] = {}
        unaligned: List[Tuple[int, List[Marker]]] = []
        chunks = [genome_ids[i:i + ALIGNED_CHUNK] for i in range(0, len(genome_ids), ALIGNED_CHUNK)]

        started = time.time()
        self.open_readers()
        try:
            with ThreadPoolExecutor(max_workers=self.cpus) as pool:
                futures = {pool.submit(self.align_chunk, chunk, markers, names): len(chunk) for chunk in chunks}
                done = done_chunks = 0
                try:
                    for future in as_completed(futures):
                        chunk_aligned, chunk_unaligned = future.result()
                        aligned.update(chunk_aligned)
                        unaligned.extend(chunk_unaligned)
                        done += futures[future]
                        done_chunks += 1
                        if done_chunks % 20 == 0 or done == len(genome_ids):
                            elapsed = time.time() - started
                            self.logger.info('Read the aligned markers of {:,} of {:,} genomes ({:.0f} s, {:.0f} '
                                             'genomes/s).'.format(done, len(genome_ids), elapsed,
                                                                  done / max(elapsed, 1e-9)))
                except BaseException:
                    # the chunks not yet started are not read
                    for future in futures:
                        future.cancel()
                    raise
        finally:
            self.close_readers()

        if unaligned:
            unaligned.sort(key=lambda u: names[u[0]].accession)
            report = write_beside_log(
                UNALIGNED_NAME, 'genome_id\tmissing_markers\tmarkers\n',
                ('{}\t{}\t{}\n'.format(names[g].accession, len(m), ','.join(x.external_id for x in m))
                 for g, m in unaligned))
            raise MSAExportError(
                '{:,} genome(s) have no row in aligned_markers for some of the {:,} markers, e.g. {}: run '
                'align_marker_genes --new_genomes for these marker sets first. Every one is in {}.'.format(
                    len(unaligned), len(markers), ', '.join(names[g].accession for g, _ in unaligned[:EXAMPLES]),
                    report))
        return aligned

    def export_msa(self, genomes: Sequence[str], marker_set_ids: Sequence[int], out_dir: str,
                   prefix: str = DEFAULT_PREFIX, individual: bool = False) -> int:
        """Write the full MSA of the genomes over the markers of the sets, and their metadata, to out_dir.

        Parameters
        ----------
        genomes : sequence of str
            As named: GB_GCA_... or GCA_....
        marker_set_ids : sequence of int
        out_dir : str
            Made if it is not there.
        prefix : str
            Of each file written.
        individual : bool
            Also write the alignment of each marker to a file of its own.

        @return: the number of genomes in the MSA.

        Raises
        ------
        MSAExportError
            What was given cannot be exported; nothing has been written but a
            list of the genomes it was refused over, beside the log.
        """

        try:
            return self._export_msa(genomes, marker_set_ids, out_dir, prefix, individual)
        finally:
            # nothing to commit: the transaction only read
            self.temp_con.rollback()

    def _export_msa(self, genomes, marker_set_ids, out_dir, prefix, individual):
        if not prefix or os.sep in prefix:
            raise MSAExportError('--prefix must name a file, not {!r}.'.format(prefix))
        if not genomes:
            raise MSAExportError('--genome_batchfile names no genomes.')
        path_prefix = os.path.join(out_dir, prefix)

        # the markers
        unknown_sets = sorted(set(marker_set_ids) - self.read_marker_sets(marker_set_ids))
        if unknown_sets:
            raise MSAExportError('--marker_set_ids names marker set(s) not in the database: {}; '
                                 'list_marker_sets lists those it holds.'.format(', '.join(map(str, unknown_sets))))
        markers = self.read_markers(marker_set_ids)
        if not markers:
            raise MSAExportError('Marker set(s) {} hold no markers.'.format(', '.join(map(str, marker_set_ids))))
        # made before anything is read, as a marker it cannot name refuses the run
        partitions = gene_partitions(markers)
        self.logger.info('Exporting the MSA of {:,} markers ({:,} columns) of marker set(s) {}.'.format(
            len(markers), sum(m.size for m in markers), ', '.join(map(str, marker_set_ids))))

        # the genomes and their alignments
        order, names, prefixes = self.resolve(genomes)
        self.logger.info('Reading the aligned markers of {:,} genomes on {} connection(s).'.format(
            len(order), self.cpus))
        aligned = self.align(order, markers, names)

        # the metadata, read before anything is written
        self.logger.info('Reading metadata_view.')
        started = time.time()
        view_columns, metadata = self.read_metadata(set(order))
        missing = [g for g in order if g not in metadata]
        if missing:
            raise MSAExportError('metadata_view has no row for {:,} genome(s), e.g. {}.'.format(
                len(missing), ', '.join(names[g].accession for g in missing[:EXAMPLES])))
        self.logger.info('Read the metadata of {:,} genomes in {:.0f} s.'.format(len(metadata), time.time() - started))
        self.logger.info('Determining NCBI species assignments for genomes in GTDB clusters.')
        clustered = {a for g in order for a in (metadata[g].get('gtdb_clustered_genomes') or '').split(';') if a}
        species = {a: ncbi_species_of(t) for a, t in self.read_ncbi_taxonomy(clustered, prefixes).items()}
        not_held = clustered - set(species)
        if not_held:
            self.logger.warning('{:,} genome(s) of gtdb_clustered_genomes are not in the database, e.g. {}, and '
                                'are not counted in {}.'.format(len(not_held), ', '.join(sorted(not_held)[:EXAMPLES]),
                                                                CLUSTER_SPECIES_FIELD))

        # the files
        os.makedirs(out_dir, exist_ok=True)
        self.logger.info('Writing the MSA, multiple hits, and ARB metadata of {:,} genomes.'.format(len(order)))
        fields = metadata_fields(view_columns)
        with output_file(path_prefix + MSA_SUFFIX) as msa_out, \
                output_file(path_prefix + MULTI_HITS_SUFFIX) as hits_out, \
                output_file(path_prefix + ARB_METADATA_SUFFIX) as arb_out, \
                output_file(path_prefix + TSV_METADATA_SUFFIX) as tsv_out:
            hits_out.write('\t'.join(['Genome_ID'] + [m.external_id for m in markers]) + '\n')
            tsv_out.write(tsv_header(fields))
            for genome_id in order:
                db_name = names[genome_id].db_name
                genome = aligned[genome_id]
                msa_out.write('>{}\n{}\n'.format(db_name, genome.sequence))
                hits_out.write('{}\t{}\n'.format(db_name, '\t'.join(genome.statuses)))
                values = genome_values(metadata[genome_id], db_name,
                                       cluster_ncbi_species(metadata[genome_id].get('gtdb_clustered_genomes'),
                                                            species))
                counts = msa_values(genome.sequence, genome.statuses)
                arb_out.write(arb_record(db_name, fields, values, counts))
                tsv_out.write(tsv_row(db_name, fields, values, counts))
        with output_file(path_prefix + ARB_FILTER_SUFFIX) as handle:
            handle.write(arb_filter(fields))
        with output_file(os.path.join(out_dir, PARTITIONS_NAME)) as handle:
            handle.write(partitions)

        with output_file(path_prefix + MARKERS_INFO_SUFFIX) as handle:
            handle.write('Marker Id\tName\tDescription\tLength (bp)\tSingle copy (%)\tUbiquity (%)\n')
            for i, marker in enumerate(markers):
                single = sum(1 for g in order if aligned[g].statuses[i] == SINGLE)
                present = sum(1 for g in order if aligned[g].statuses[i] != MISSING)
                handle.write('\t'.join([marker.external_id, marker.name, marker.description, str(marker.size),
                                        '%.2f' % (single * 100.0 / len(order)),
                                        '%.2f' % (present * 100.0 / len(order))]) + '\n')

        if individual:
            self.logger.info('Writing individual alignments.')
            start = 0
            for marker in markers:
                end = start + marker.size
                with output_file('{}_{}{}'.format(path_prefix, marker.id_in_database, INDIVIDUAL_SUFFIX)) as handle:
                    for genome_id in order:
                        handle.write('>{}\n{}\n'.format(names[genome_id].db_name,
                                                        aligned[genome_id].sequence[start:end]))
                start = end

        self.logger.info('Wrote the MSA of {:,} genomes, {:,} columns, to {}.'.format(
            len(order), sum(m.size for m in markers), path_prefix + MSA_SUFFIX))
        return len(order)
