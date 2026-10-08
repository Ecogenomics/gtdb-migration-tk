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
import re
import sys
import glob
import logging
from collections import Counter, defaultdict
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

from psycopg2.extras import execute_values

from gtdb_migration_tk.biolib_lite.common import canonical_gid
from gtdb_migration_tk.biolib_lite.taxonomy import Taxonomy
from gtdb_migration_tk.database_configuration import GenomeDatabaseConnectionFTPUpdate
# reset_unwritten() lives beside one_transaction, update_propagated_tax using it too
from gtdb_migration_tk.database_configuration.GenomeDatabaseConnectionFTPUpdate import (
    RESET_UNWRITTEN, WRITTEN_TABLE, one_transaction, reset_unwritten)
from gtdb_migration_tk.gtdb_lite.gtdb_importer import (SKIP, UNKNOWN_EXAMPLES, GTDBImporter,
                                                       UnknownGenomesError, id_at_source)
from gtdb_migration_tk.utils.common import GZIP_SUFFIX, open_text

# rows a statement of update_type_designation hands the server at a time
PAGE_SIZE = 10000

# How many genomes of a metadata table update_metadata_db holds at once. A table
# is read as it is written, each field of this many genomes handed to the
# importer before the next are read, so what is held is bounded by this rather
# than by the release: r237's NCBI metadata is 1.35M genomes of 37 fields, and
# releases only grow. Every chunk is written in the command's one transaction,
# so a run stopped part way still leaves the database as it was.
CHUNK_GENOMES = 100000

# the types a description gives a field that only a whole number can be written to
INTEGER_TYPES = ('INT', 'INTEGER')
WHOLE_NUMBER = re.compile(r'^[+-]?[0-9]+$')

# how many lines or genomes an error names
EXAMPLES = 10

# Files a command writes beside its tables that are reports, not tables, so are
# passed over in --input_folder rather than refused as tables this command does
# not know: ncbi_strains writes substrains.tsv beside strain_summary_file.tsv.gz,
# which r237 wrote into the folder update_metadata_db loads, and
# parse_ncbi_genome_category the GenBank file lines it found each category on
# beside ncbi_genome_category.tsv.gz.
NOT_TABLES = ('substrains.tsv', 'ncbi_genome_category_evidence.tsv')


# update_ncbi_tax_db: what each of its three files is written to, in the order
# they are loaded and the error file names them
NCBI_TAX_FIELDS = (('metadata_ncbi', 'ncbi_organism_name'),
                   ('metadata_taxonomy', 'ncbi_taxonomy'),
                   ('metadata_taxonomy', 'ncbi_taxonomy_unfiltered'))

# update_ncbi_tax_db: in --out_dir, each genome without one of the three, and
# what a missing value is marked with there
NCBI_TAX_MISSING_NAME = 'ncbi_tax_missing.tsv'
MISSING = 'missing'


def taxonomy_string(value: str) -> str:
    """A taxonomy as update_ncbi_tax_db writes it: no trailing ';', no space around a rank.

    The same string biolib_lite's Taxonomy().read() gave, joined again with ';',
    without its failing on an empty one.

    Parameters
    ----------
    value : str
        e.g. 'd__Bacteria; p__Pseudomonadota;'.

    @return: e.g. 'd__Bacteria;p__Pseudomonadota', or '' for an empty value.
    """

    value = value.strip()
    if value.endswith(';'):
        value = value[:-1]
    return ';'.join(rank.strip() for rank in value.split(';')) if value else ''


class MetadataTableError(ValueError):
    """A metadata table, or its description, that update_metadata_db will not load."""


def read_descriptions(paths: Sequence[str]) -> Dict[str, Tuple[str, str]]:
    """The type and database table of each field the description files name.

    A table update_metadata_db knows may be described by more than one file
    (strain_summary_file.tsv by two); their fields are merged so that the table
    is read once.

    Parameters
    ----------
    paths : sequence of str
        Description files: field, description, type, table, tab-separated.

    @return: field -> (type, table), e.g. {'ncbi_taxid': ('INT', 'metadata_ncbi')}.

    Raises
    ------
    MetadataTableError
        A line has fewer than four columns, or two lines give one field a
        different type or table.
    """

    described: Dict[str, Tuple[str, str]] = {}
    for path in paths:
        with open(path) as handle:
            for line_number, line in enumerate(handle, start=1):
                if not line.strip():
                    continue
                columns = line.rstrip('\n').split('\t')
                if len(columns) < 4:
                    raise MetadataTableError('{} line {:,} has {} column(s); a description gives a '
                                             'field, its description, its type and its table.'.format(
                                                 path, line_number, len(columns)))
                field, data_type, table = columns[0].strip(), columns[2].strip(), columns[3].strip()
                if field in described and described[field] != (data_type, table):
                    raise MetadataTableError('{} describes {} as {} in {}, and an earlier description '
                                             'as {} in {}.'.format(path, field, data_type, table,
                                                                   *described[field]))
                described[field] = (data_type, table)
    return described


def whole_number(value: str) -> Optional[str]:
    """The value of an INT field as a whole number, or None if it is not one.

    Parameters
    ----------
    value : str
        e.g. '12', '12.0' or '0.0'; '12.5', 'na' and '' are not whole numbers.

    @return: e.g. '12', or None.
    """

    value = value.strip()
    if WHOLE_NUMBER.match(value):
        return value
    try:
        number = float(value)
    except ValueError:
        return None
    return str(int(number)) if number.is_integer() else None


def read_genome_list(path: str) -> Set[str]:
    """The genomes of a --genome_list file, as genomes.id_at_source names them.

    Parameters
    ----------
    path : str
        A table whose first column, tab or comma separated, names the genomes,
        e.g. a metadata file exported from GTDB; GB_GCA_000003645.1 and
        GCA_000003645.1 are the same genome.

    @return: e.g. {'GCA_000003645.1', ...}.
    """

    genomes = set()
    with open(path) as handle:
        for line in handle:
            if not line.strip():
                continue
            separator = '\t' if len(line.split('\t')) >= len(line.split(',')) else ','
            genomes.add(id_at_source(line.rstrip().split(separator)[0].strip()))
    genomes.discard(None)
    return genomes


def confirm_partial_load(genome_list: str, command: str = 'update_metadata_db') -> None:
    """Ask before loading a genome list's genomes after setting the fields to NULL for every genome.

    Asked on the terminal whatever runs it: with no one to answer (nohup, a
    script) the answer is end of file, and the run ends there, having changed
    nothing.

    Parameters
    ----------
    genome_list : str
        The --genome_list file.
    command : str
        The command asking, as its messages name it.

    @return: None, where the answer is yes.

    Raises
    ------
    SystemExit
        The answer is anything but yes, or there is none; the run exits 1.
    """

    question = ('--genome_list is given without --do_not_null_field: every field loaded will be '
                'set to NULL for every genome in the database, and written again only for the '
                'genomes in {}. The metadata of every other genome will be removed. '
                'Proceed? [y/n] '.format(genome_list))
    try:
        answer = input(question)
    except EOFError:
        sys.exit('{}: no answer to whether to proceed (is there no terminal?); '
                 'nothing was changed. Give --do_not_null_field to keep the other genomes\' '
                 'metadata.'.format(command))
    if answer.strip().lower() not in ('y', 'yes'):
        sys.exit('{}: not proceeding; nothing was changed.'.format(command))

# what update_type_designation writes
TYPE_STRAIN_OF_SPECIES = 'type strain of species'
NOT_USED_AS_TYPE = 'not used as type'
SEQCODE_SOURCE = 'Seqcode'


def type_designation_changes(rows: Iterable[Tuple]) -> Tuple[Dict[int, str], Dict[int, str]]:
    """What update_type_designation changes, decided without the database.

    A genome whose species is valid under the SeqCode is the type strain of its
    species, and 'Seqcode' is added to the sources of that designation. A genome
    NCBI excludes from RefSeq as derived from a metagenome and not used as type
    is not used as type, whatever the SeqCode says.

    Parameters
    ----------
    rows : iterable of tuple
        (id, seqcode_species_status, ncbi_excluded_from_refseq,
        gtdb_type_designation_ncbi_taxa, gtdb_type_designation_ncbi_taxa_sources)
        of each genome, as update_type_designation selects them.

    @return: the new gtdb_type_designation_ncbi_taxa of each genome it changes,
             and the new gtdb_type_designation_ncbi_taxa_sources, each by id. The
             sources keep their order, 'Seqcode' added last where it is not
             already there: they went through a set, and came out in an order
             that differed between runs.
    """

    designations: Dict[int, str] = {}
    sources: Dict[int, str] = {}
    for genome_id, seqcode_status, excluded, _designation, typed_sources in rows:
        if seqcode_status is not None and 'Valid' in seqcode_status:
            designations[genome_id] = TYPE_STRAIN_OF_SPECIES
            parts = [part for part in (typed_sources or '').split(';') if part]
            sources[genome_id] = ';'.join(dict.fromkeys(parts + [SEQCODE_SOURCE]))

        if (excluded is not None and 'not used as type' in excluded
                and 'derived from metagenome' in excluded):
            designations[genome_id] = NOT_USED_AS_TYPE

    return designations, sources

class MetadataDatabaseManager(object):

    def __init__(self, database):
        """Initialization."""
        self.logger = logging.getLogger('timestamp')
        self.description_table = {'metadata_gene.tsv':['metadata_gene.desc.tsv'],
                                  'metadata_nt.tsv':['metadata_nt.desc.tsv'],
                                  'metadata_ssu_silva.tsv':['metadata_rna.table.desc.tsv','metadata_sequence.desc.tsv'],
                                  'metadata_lsu_silva_23s.tsv':['metadata_rna.table.desc.tsv','metadata_sequence.desc.tsv'],
                                  'metadata_lsu_5S.tsv':['metadata_rna.table.desc.tsv','metadata_sequence.desc.tsv'],
                                  'metadata_ssu_silva_count.tsv':['metadata_ssu_count.desc.tsv'],
                                  'metadata_lsu_silva_23s_count.tsv':['metadata_ssu_count.desc.tsv'],
                                  'metadata_lsu_5S_count.tsv':['metadata_ssu_count.desc.tsv'],
                                  'metadata_trna_count.tsv':['metadata_trna.desc.tsv'],
                                  'ncbi_assembly_summary.tsv':['metadata_ncbi_assembly_file.desc.tsv'],
                                  'strain_summary_file.tsv':['metadata_ncbi_assembly.desc.tsv','metadata_ncbi_assembly_file.desc.tsv'],
                                  'ncbi_assembly_metadata.tsv':['metadata_ncbi_assembly.desc.tsv'],
                                  'ncbi_genome_category.tsv':['metadata_ncbi_genome_category.desc.tsv']
                                  }

        self.temp_con = GenomeDatabaseConnectionFTPUpdate.GenomeDatabaseConnectionFTPUpdate(database)
        self.temp_con.MakePostgresConnection()
        self.temp_cur = self.temp_con.cursor()

    def tables_to_load(self, table_folder: Optional[str], table_file: Optional[str],
                       table_file_desc: Optional[str]) -> List[Tuple[str, List[str]]]:
        """The tables a run loads, each with its description files; nothing is read but names.

        Parameters
        ----------
        table_folder : str or None
            --input_folder: every .tsv in it, each a table this command knows.
        table_file, table_file_desc : str or None
            --metadata_table and --metadata_table_desc.

        @return: [(table, [description file, ...]), ...], the folder's in name order.

        Raises
        ------
        MetadataTableError
            Neither or both of a folder and a table are given, a table without
            its description, an empty folder, or a folder holding a table this
            command does not know -- every one named, before anything is written.
        """

        if (table_folder is None) == (table_file is None):
            raise MetadataTableError('Give --input_folder, or --metadata_table and --metadata_table_desc; '
                                     'not both, and not neither.')
        if table_file is not None:
            if table_file_desc is None:
                raise MetadataTableError('--metadata_table {} is given without --metadata_table_desc, '
                                         'which says what each of its fields is.'.format(table_file))
            return [(table_file, [table_file_desc])]

        desc_table_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                      'data_files', 'table_description')
        # a table is known by its name without GZIP_SUFFIX: create_tables,
        # parse_ncbi_assemblies, parse_ncbi_dir, parse_ncbi_genome_category and
        # ncbi_strains write them gzipped
        tables = sorted(glob.glob(os.path.join(table_folder, '*.tsv'))
                        + glob.glob(os.path.join(table_folder, '*.tsv' + GZIP_SUFFIX)))
        reports = [table for table in tables if os.path.basename(table) in NOT_TABLES]
        if reports:
            self.logger.info('Passing over {}, not a table.'.format(
                ', '.join(os.path.basename(report) for report in reports)))
            tables = [table for table in tables if table not in reports]
        if not tables:
            raise MetadataTableError('{} holds no .tsv or .tsv.gz table.'.format(table_folder))

        def name(table):
            base = os.path.basename(table)
            return base[:-len(GZIP_SUFFIX)] if base.endswith(GZIP_SUFFIX) else base

        unknown = [os.path.basename(table) for table in tables if name(table) not in self.description_table]
        if unknown:
            raise MetadataTableError('{} holds {:,} table(s) this command does not know: {}. It knows {}, '
                                     'gzipped or not.'.format(table_folder, len(unknown), ', '.join(unknown),
                                                             ', '.join(sorted(self.description_table))))
        both = sorted(name(table) for table in tables if table.endswith(GZIP_SUFFIX)
                      and table[:-len(GZIP_SUFFIX)] in tables)
        if both:
            raise MetadataTableError('{} holds {:,} table(s) both gzipped and not: {}. Each would be loaded '
                                     'twice; remove the one that is not of this release.'.format(
                                         table_folder, len(both), ', '.join(both)))
        return [(table, [os.path.join(desc_table_dir, desc) for desc in self.description_table[name(table)]])
                for table in tables]

    @one_transaction
    def process_metadata_files(self, genome_list_file, do_not_null_field=False, table_folder=None,
                               table_file=None, table_file_desc=None):
        """Load metadata tables into the database, in one transaction.

        A genome is loaded where it is in --genome_list, or, without one, in the
        genomes table: a table may cover more genomes than the database holds
        (parse_ncbi_assemblies' covers every genome of NCBI's summaries), and
        the importer refuses a field naming a genome it does not hold. A genome
        of --genome_list the database does not hold is still refused.

        Parameters
        ----------
        genome_list_file : str or None
            --genome_list.
        do_not_null_field : bool
            Keep what the fields hold for genomes the tables do not write.
        table_folder, table_file, table_file_desc : str or None
            As tables_to_load() takes them.

        @return: None

        Raises
        ------
        MetadataTableError
            A table or description that cannot be loaded as it is; nothing is written.
        UnknownGenomesError
            A genome of --genome_list is not in the database; nothing is written.
        """

        tables = self.tables_to_load(table_folder, table_file, table_file_desc)
        importer = GTDBImporter(self.temp_cur)
        if genome_list_file:
            keep = read_genome_list(genome_list_file)
            keep_source = genome_list_file
        else:
            keep = importer.genomes()
            keep_source = 'the genomes table'
        self.logger.info('Loading {:,} table(s) for the {:,} genomes of {}.'.format(
            len(tables), len(keep), keep_source))

        for table_idx, (metadata_file, description_files) in enumerate(tables, start=1):
            self.logger.info('Processing file {}/{}: {}'.format(table_idx, len(tables), metadata_file))
            self.update_metadata_db(metadata_file, description_files, keep, keep_source,
                                    importer, do_not_null_field)

    def update_metadata_db(self, metadata_file: str, description_files: Sequence[str], keep: Set[str],
                           keep_source: str, importer: GTDBImporter, do_not_null_field: bool) -> None:
        """Load one metadata table, inside the caller's transaction.

        The table is read once, as it is written, CHUNK_GENOMES genomes at a
        time, each field of a chunk handed to the importer before the next
        chunk is read. A table is refused -- the transaction rolled back, so
        nothing of the run is written -- where a row has more or fewer columns
        than the header, a genome is named twice, or an INT field holds a value
        that is not a whole number: each was passed over, the row's fields left
        unset, the last row kept, or 12.5 written as 12.

        Parameters
        ----------
        metadata_file : str
            Tab-separated, gzipped or not: the genome, then a column per field.
        description_files : sequence of str
            The descriptions of its fields (read_descriptions()).
        keep : set of str
            The genomes to load, as genomes.id_at_source names them.
        keep_source : str
            Where keep came from, for the log.
        importer : GTDBImporter
            The importer of the caller's transaction.
        do_not_null_field : bool
            Keep what the fields hold for genomes the table does not write.

        @return: None

        Raises
        ------
        MetadataTableError
            The table cannot be loaded as it is.
        """

        descriptions = read_descriptions(description_files)
        with open_text(metadata_file) as f:
            header = [column.strip() for column in f.readline().rstrip('\n').split('\t')]
            loaded = [(index, field) for index, field in enumerate(header)
                      if index > 0 and field in descriptions]
            undescribed = [field for field in header[1:] if field not in descriptions]
            self.logger.info('{} holds {:,} field(s) to load: {}.'.format(
                metadata_file, len(loaded), ', '.join(field for _, field in loaded)))
            if undescribed:
                self.logger.info('{:,} column(s) of {} are in no description and are not loaded: {}.'.format(
                    len(undescribed), metadata_file, ', '.join(undescribed)))

            written = Counter()
            # the genomes of a field left without a value, which with seen
            # say whom the field was written for (reset_unwritten())
            empty: Dict[str, Set[str]] = {field: set() for _, field in loaded}
            chunk: Dict[str, List[Tuple[str, str]]] = {field: [] for _, field in loaded}
            in_chunk = 0
            seen: Set[str] = set()
            rows = 0
            skipped = 0

            def flush():
                for _, field in loaded:
                    if chunk[field]:
                        data_type, table = descriptions[field]
                        importer.import_metadata_to_db(table, field, data_type, chunk[field])
                        written[field] += len(chunk[field])
                        chunk[field] = []

            for line_number, line in enumerate(f, start=2):
                if not line.strip():
                    continue
                row = line.rstrip('\n').split('\t')
                if len(row) != len(header):
                    raise MetadataTableError('{} line {:,} ({}) has {} column(s) where its header has {}; '
                                             'nothing was written.'.format(
                                                 metadata_file, line_number, row[0], len(row), len(header)))
                rows += 1
                genome_id = id_at_source(row[0].strip())
                if genome_id not in keep:
                    skipped += 1
                    continue
                if genome_id in seen:
                    raise MetadataTableError('{} names {} more than once, again on line {:,}; the table holds '
                                             'more than one run. Nothing was written.'.format(
                                                 metadata_file, genome_id, line_number))
                seen.add(genome_id)

                for index, field in loaded:
                    value = row[index]
                    if not value.strip():
                        empty[field].add(genome_id)
                        continue
                    if descriptions[field][0].upper() in INTEGER_TYPES:
                        number = whole_number(value)
                        if number is None:
                            raise MetadataTableError('{} line {:,} gives {} {} the value {!r}, which is not a '
                                                     'whole number, and {} is {}; nothing was written.'.format(
                                                         metadata_file, line_number, genome_id, field, value,
                                                         field, descriptions[field][0]))
                        value = number
                    chunk[field].append((genome_id, value))

                in_chunk += 1
                if in_chunk >= CHUNK_GENOMES:
                    flush()
                    in_chunk = 0
            flush()

        self.logger.info('Read {:,} genomes of {}: {:,} loaded, {:,} not in {} and skipped.'.format(
            rows, metadata_file, rows - skipped, skipped, keep_source))
        for _, field in loaded:
            self.logger.info('Wrote {}.{} for {:,} genomes.'.format(descriptions[field][1], field, written[field]))

        # the field set to NULL for every genome holding a value this table did
        # not write, unless asked not to. It is in the transaction the new
        # values are written in, so a failure leaves the values the fields held
        if not do_not_null_field:
            for _, field in loaded:
                cleared = reset_unwritten(self.temp_cur, descriptions[field][1], field,
                                          (gid for gid in seen if gid not in empty[field]))
                self.logger.info('Set {}.{} to NULL for {:,} genome(s) holding a value this table did '
                                 'not give them.'.format(descriptions[field][1], field, cleared))

    @one_transaction
    def update_reps(self, final_cluster_file):
        """Update representative genomes of species clusters.

        The cluster file is read whole, and every genome it names found in the
        database, before the fields are set to NULL. The NULLs were committed
        first, and a genome of the file the database did not hold then raised a
        KeyError, which left the database with no representatives at all.
        """

        # mark all genomes as not being representatives and get translation
        # between canonical genome IDs and NCBI accessions
        q = ("SELECT accession FROM metadata_view")
        self.temp_cur.execute(q)

        gid_to_ncbi_accn = {}
        is_rep = {}
        for r in self.temp_cur:
            ncbi_accn = r[0]
            gid_to_ncbi_accn[canonical_gid(ncbi_accn)] = ncbi_accn
            is_rep[ncbi_accn] = False

        # determine representative assignment of genomes
        genome_rep_data = []
        not_held = []
        num_sp_reps = 0
        with open(final_cluster_file) as f:
            headers = f.readline().strip().split('\t')

            rep_index = headers.index('Representative')
            clustered_genomes_index = headers.index('Clustered genomes')

            for line in f:
                line_split = line.strip().split('\t')

                rep_gid = line_split[rep_index]
                if rep_gid not in gid_to_ncbi_accn:
                    not_held.append(rep_gid)
                    continue
                rep_accn = gid_to_ncbi_accn[rep_gid]
                if len(line_split) > clustered_genomes_index:
                    gids = [gid.strip() for gid in line_split[clustered_genomes_index].split(',')]
                    for gid in gids:
                        if gid not in gid_to_ncbi_accn:
                            not_held.append(gid)
                            continue
                        ncbi_accn = gid_to_ncbi_accn[gid]
                        genome_rep_data.append((ncbi_accn, rep_accn))

                genome_rep_data.append((rep_accn, rep_accn))

                is_rep[rep_accn] = True
                num_sp_reps += 1

        if not_held:
            raise UnknownGenomesError(
                '{:,} genome(s) of {} are not in the database, e.g. {}; the database '
                'holds another release, or update_db has not been run.'.format(
                    len(not_held), final_cluster_file, ', '.join(not_held[:UNKNOWN_EXAMPLES])))

        # clear representative fields
        self.logger.info('Setting GTDB representative fields to NULL.')
        q = ("UPDATE metadata_taxonomy SET gtdb_representative = NULL, gtdb_genome_representative = NULL")
        self.temp_cur.execute(q)

        print(f'Identified {num_sp_reps:,} species clusters.')
        print('Identified {:,} genomes marked as representatives.'.format(sum([1 for rid in is_rep if is_rep[rid]])))
        gtdbimporter = GTDBImporter(self.temp_cur)
        gtdbimporter.import_metadata_to_db('metadata_taxonomy', 'gtdb_genome_representative', 'TEXT', genome_rep_data)

        # mark representative genomes
        is_rep_data = []
        for rep_accn, rep_status in is_rep.items():
            is_rep_data.append((rep_accn, str(rep_status)))

        gtdbimporter.import_metadata_to_db('metadata_taxonomy', 'gtdb_representative', 'BOOLEAN', is_rep_data)

    @one_transaction
    def add_surveillance_genomes(self, genome_list):
        """Replace the surveillance genomes with those of a list.

        The table is emptied and filled in one transaction, so a list that cannot
        be inserted leaves the genomes it held. It was emptied and committed
        first, after a question that was answered 'y' in the code rather than
        asked. Blank lines are not genomes, and a genome listed twice is inserted
        once.
        """

        genomes = []
        with open(genome_list) as glf:
            for line in glf:
                gid = line.strip()
                if gid:
                    genomes.append(gid)
        unique = list(dict.fromkeys(genomes))

        self.logger.info('Replacing the surveillance genomes with the {:,} of {}{}.'.format(
            len(unique), genome_list,
            ' ({:,} listed more than once)'.format(len(genomes) - len(unique))
            if len(unique) < len(genomes) else ''))
        self.temp_cur.execute('TRUNCATE survey_genomes')

        q_add = "INSERT INTO survey_genomes(canonical_gid) VALUES (%s) "
        self.temp_cur.executemany(q_add, [(gid,) for gid in unique])

    @one_transaction
    def update_type_designation(self):
        """Set the type designation of genomes from SeqCode and NCBI's exclusions.

        What changes is decided first (type_designation_changes()) and written in
        two statements. Each genome was updated, and committed, on its own, with
        its values written into the SQL; a sources value holding a quote broke the
        statement.
        """

        self.temp_cur.execute("SELECT mn.id,seq.seqcode_species_status,mn.ncbi_excluded_from_refseq,"
                              "mtm.gtdb_type_designation_ncbi_taxa,mtm.gtdb_type_designation_ncbi_taxa_sources "
                              "from metadata_ncbi mn "
                              "LEFT JOIN metadata_seqcode seq USING (id) "
                              "LEFT JOIN metadata_type_material mtm USING (id)")
        rows = self.temp_cur.fetchall()
        self.logger.info('Loaded {:,} genomes.'.format(len(rows)))

        designations, sources = type_designation_changes(rows)

        execute_values(
            self.temp_cur,
            'UPDATE metadata_type_material AS m SET gtdb_type_designation_ncbi_taxa = v.designation '
            'FROM (VALUES %s) AS v(id, designation) WHERE m.id = v.id',
            list(designations.items()), template='(%s::integer, %s::text)', page_size=PAGE_SIZE)
        execute_values(
            self.temp_cur,
            'UPDATE metadata_type_material AS m SET gtdb_type_designation_ncbi_taxa_sources = v.sources '
            'FROM (VALUES %s) AS v(id, sources) WHERE m.id = v.id',
            list(sources.items()), template='(%s::integer, %s::text)', page_size=PAGE_SIZE)

        for designation, count in sorted(Counter(designations.values()).items()):
            self.logger.info("Set gtdb_type_designation_ncbi_taxa to '{}' for {:,} genomes.".format(
                designation, count))
        self.logger.info("Added '{}' to gtdb_type_designation_ncbi_taxa_sources of {:,} genomes "
                         "valid under the SeqCode.".format(SEQCODE_SOURCE, len(sources)))


class NCBITaxDatabaseManager(object):
    """Add organism name to GTDB."""

    def __init__(self, database):
        """Initialization."""
        self.logger = logging.getLogger('timestamp')

        self.temp_con = GenomeDatabaseConnectionFTPUpdate.GenomeDatabaseConnectionFTPUpdate(database)
        self.temp_con.MakePostgresConnection()
        self.temp_cur = self.temp_con.cursor()


    def load_two_column_file(self, path: str, table: str, field: str, keep: Set[str],
                             importer: GTDBImporter, normalise=str.strip) -> Set[str]:
        """Write one field from a file of genome and value, in chunks, inside the caller's transaction.

        Parameters
        ----------
        path : str
            Tab-separated, no header: the genome (GB_GCA_... or GCA_...) and its value.
        table, field : str
            Where the value is written, e.g. metadata_ncbi, ncbi_organism_name.
        keep : set of str
            The genomes to write, as genomes.id_at_source names them.
        importer : GTDBImporter
            The importer of the caller's transaction.
        normalise : callable
            What is made of a value before it is written, e.g. taxonomy_string.

        @return: the genomes of keep given a value; one with no line, or an
                 empty value, is not.

        Raises
        ------
        MetadataTableError
            The file names a genome of keep twice; nothing is written.
        """

        written: Set[str] = set()
        chunk: List[Tuple[str, str]] = []
        rows = 0
        with open_text(path) as handle:
            for line_number, line in enumerate(handle, start=1):
                if not line.strip():
                    continue
                rows += 1
                genome, _, value = line.rstrip('\n').partition('\t')
                genome_id = id_at_source(genome.strip())
                if genome_id not in keep:
                    continue
                if genome_id in written:
                    raise MetadataTableError('{} names {} more than once, again on line {:,}; nothing was '
                                             'written.'.format(path, genome_id, line_number))
                value = normalise(value)
                if not value:
                    continue
                written.add(genome_id)
                chunk.append((genome_id, value))
                if len(chunk) >= CHUNK_GENOMES:
                    importer.import_metadata_to_db(table, field, 'TEXT', chunk, unknown=SKIP)
                    chunk = []
        if chunk:
            importer.import_metadata_to_db(table, field, 'TEXT', chunk, unknown=SKIP)

        self.logger.info('Wrote {}.{} for {:,} of the {:,} genomes of {}.'.format(
            table, field, len(written), rows, path))
        return written

    @one_transaction
    def update_ncbi_tax_db(self, organism_name_file: str, filtered_file: str, unfiltered_file: str,
                           genome_list_file: Optional[str], out_dir: str, do_not_null_field: bool = False) -> None:
        """Write the NCBI organism name, taxonomy and unfiltered taxonomy of each genome, in one transaction.

        A genome is written where it is in --genome_list, or, without one, in the
        genomes table: the files cover every assembly NCBI holds. Each genome so
        chosen without one of the three -- no line in its file, or an empty
        value -- is written to NCBI_TAX_MISSING_NAME in out_dir, which is written
        whether or not there are any, and counted in a WARNING for each field.

        Parameters
        ----------
        organism_name_file : str
            Genome and organism name, as parse_ncbi_taxonomy writes them.
        filtered_file, unfiltered_file : str
            Genome and standardised taxonomy, and genome and taxonomy as NCBI gives it.
        genome_list_file : str or None
            --genome_list: a table whose first column names the genomes.
        out_dir : str
            Directory the error file is written to.
        do_not_null_field : bool
            Keep what the fields hold for genomes the files do not write.

        @return: None

        Raises
        ------
        MetadataTableError
            A file names a genome twice; nothing is written.
        """

        importer = GTDBImporter(self.temp_cur)
        if genome_list_file:
            keep = read_genome_list(genome_list_file)
            keep_source = genome_list_file
        else:
            keep = importer.genomes()
            keep_source = 'the genomes table'
        self.logger.info('Writing the NCBI organism names and taxonomies of the {:,} genomes of {}.'.format(
            len(keep), keep_source))

        present = []
        for (table, field), path, normalise in zip(NCBI_TAX_FIELDS,
                                                   (organism_name_file, filtered_file, unfiltered_file),
                                                   (str.strip, taxonomy_string, taxonomy_string)):
            given = self.load_two_column_file(path, table, field, keep, importer, normalise)
            if not do_not_null_field:
                cleared = reset_unwritten(self.temp_cur, table, field, given)
                self.logger.info('Set {}.{} to NULL for {:,} genome(s) holding a value {} did not give '
                                 'them.'.format(table, field, cleared, path))
            present.append(given)

        self.write_missing(keep, present, keep_source, out_dir)

    def write_missing(self, keep: Set[str], present: Sequence[Set[str]], keep_source: str, out_dir: str) -> str:
        """Write each genome without an organism name, taxonomy or unfiltered taxonomy.

        Parameters
        ----------
        keep : set of str
            The genomes that should have all three.
        present : sequence of set of str
            For each of NCBI_TAX_FIELDS, the genomes given a value.
        keep_source : str
            Where keep came from, for the log.
        out_dir : str
            Directory NCBI_TAX_MISSING_NAME is written to.

        @return: the path of the file written.
        """

        path = os.path.join(out_dir, NCBI_TAX_MISSING_NAME)
        missing = sorted(gid for gid in keep if any(gid not in given for given in present))
        with open(path, 'w') as handle:
            handle.write('\t'.join(['genome_id'] + [field for _table, field in NCBI_TAX_FIELDS]) + '\n')
            for gid in missing:
                handle.write('\t'.join([gid] + [MISSING if gid not in given else '' for given in present]) + '\n')

        for (_table, field), given in zip(NCBI_TAX_FIELDS, present):
            without = sorted(keep - given)
            if without:
                self.logger.warning('{:,} genome(s) of {} have no {}, e.g. {}; each is in {}.'.format(
                    len(without), keep_source, field, ', '.join(without[:EXAMPLES]), path))
        if not missing:
            self.logger.info('Every genome of {} has an organism name, a taxonomy and an unfiltered '
                             'taxonomy.'.format(keep_source))
        return path
