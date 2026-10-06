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
import glob
import logging
from collections import Counter, defaultdict
from typing import Dict, Iterable, Tuple

from psycopg2.extras import execute_values

from gtdb_migration_tk.biolib_lite.common import canonical_gid
from gtdb_migration_tk.biolib_lite.taxonomy import Taxonomy
from gtdb_migration_tk.database_configuration import GenomeDatabaseConnectionFTPUpdate
from gtdb_migration_tk.database_configuration.GenomeDatabaseConnectionFTPUpdate import one_transaction
from gtdb_migration_tk.gtdb_lite.gtdb_importer import (SKIP, UNKNOWN_EXAMPLES, GTDBImporter,
                                                       UnknownGenomesError)
from gtdb_migration_tk.utils.common import open_text

# rows a statement of update_type_designation hands the server at a time
PAGE_SIZE = 10000

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
                                  'metadata_ssu_gg.tsv':['metadata_rna.table.desc.tsv'],
                                  'metadata_ssu_silva.tsv':['metadata_rna.table.desc.tsv','metadata_sequence.desc.tsv'],
                                  'metadata_lsu_silva_23s.tsv':['metadata_rna.table.desc.tsv','metadata_sequence.desc.tsv'],
                                  'metadata_lsu_5S.tsv':['metadata_rna.table.desc.tsv','metadata_sequence.desc.tsv'],
                                  'metadata_ssu_silva_count.tsv':['metadata_ssu_count.desc.tsv'],
                                  'metadata_lsu_silva_23s_count.tsv':['metadata_ssu_count.desc.tsv'],
                                  'metadata_lsu_5S_count.tsv':['metadata_ssu_count.desc.tsv'],
                                  'metadata_trna_count.tsv':['metadata_trna.desc.tsv'],
                                  'ncbi_assembly_summary.tsv':['metadata_ncbi_assembly_file.desc.tsv'],
                                  'strain_summary_file.tsv':['metadata_ncbi_assembly.desc.tsv','metadata_ncbi_assembly_file.desc.tsv'],
                                  'ncbi_assembly_metadata.tsv':['metadata_ncbi_assembly.desc.tsv']
                                  }

        self.temp_con = GenomeDatabaseConnectionFTPUpdate.GenomeDatabaseConnectionFTPUpdate(database)
        self.temp_con.MakePostgresConnection()
        self.temp_cur = self.temp_con.cursor()

    @one_transaction
    def process_metadata_files(self,genome_list_file,do_not_null_field=False,table_folder=None,table_file=None,table_file_desc=None):
        file_dir = os.path.dirname(os.path.realpath(__file__))
        desc_table_dir = os.path.join(file_dir, 'data_files', 'table_description')
        if table_folder is not None:
            list_tsv_files = glob.glob(os.path.join(table_folder, '*.tsv'))
            #list_tsv_files = ['metadata_trna_count.tsv']
            for tsv_file in list_tsv_files:
                if os.path.basename(tsv_file) not in self.description_table:
                    print(f'{os.path.basename(tsv_file)} is not a standard table')
                    sys.exit(-1)
            for tsv_idx,tsv_file in enumerate(list_tsv_files):
                self.logger.info(f'Processing file {tsv_idx+1}/{len(list_tsv_files)}: {os.path.basename(tsv_file)}')
                for desc_file in self.description_table.get(os.path.basename(tsv_file)):
                    metadata_desc_file = os.path.join(desc_table_dir,desc_file)
                    self.update_metadata_db(tsv_file,metadata_desc_file,genome_list_file,do_not_null_field)
        elif table_file is not None:
            self.logger.info("Processing specific file")
            self.update_metadata_db(table_file, table_file_desc, genome_list_file, do_not_null_field)

    def update_metadata_db(self,metadata_file,metadata_desc_file,genome_list_file,do_not_null_field):
        # get fields in metadata file
        gtdbimporter = GTDBImporter(self.temp_cur)
        self.logger.info('Parsing metadata file: %s' % metadata_file)
        with open_text(metadata_file) as f:
            metadata_fields = f.readline().strip().split('\t')[1:]
        self.logger.info(
            'Metadata file contains {} fields.'.format(len(metadata_fields)))
        self.logger.info('Fields: %s' % ', '.join(metadata_fields))

        # get database table and data type of each metadata field
        metadata_type = {}
        metadata_table = {}
        with open(metadata_desc_file) as f:
            for line in f:
                line_split = line.strip('\n').split('\t')
                field = line_split[0]
                if field in metadata_fields:
                    metadata_type[field] = line_split[2]
                    metadata_table[field] = line_split[3]
        self.logger.info('Identified {} matching fields in metadata description file.'.format(
            len(metadata_table)))
        self.logger.info('Fields: %s' % ', '.join(metadata_table))

        # set fields to NULL for every genome, unless asked not to. This asked
        # [y/n] first, which a run under nohup or with no terminal could not
        # answer; it is in the transaction the new values are written in, so a
        # failure leaves the values the fields held rather than NULL
        if not do_not_null_field:
            for field in metadata_table:
                self.logger.info('Setting {}.{} to NULL for every genome.'.format(
                    metadata_table[field], field))
                q = ("UPDATE {} SET {} = NULL".format(
                    metadata_table[field], field))
                self.temp_cur.execute(q)

        # get genomes to process
        genome_list = set()
        if genome_list_file:
            for line in open(genome_list_file):
                if len(line.split('\t')) >= len(line.split(',')):
                    genome_list.add(line.rstrip().split('\t')[0])
                else:
                    genome_list.add(line.rstrip().split(',')[0])
        self.logger.info('Processing {} genomes.'.format(len(genome_list)))

        # read metadata file
        metadata = defaultdict(lambda: defaultdict(str))
        with open_text(metadata_file) as f:
            fields = [x.strip() for x in f.readline().split('\t')]

            for line in f:
                line_split = line.rstrip('\n').split('\t')

                genome_id = line_split[0]
                # print line_split
                for i, value in enumerate(line_split[1:]):
                    metadata[fields[i + 1]][genome_id] = value

        # add each field to the database
        for field in metadata:
            data_to_commit = []

            if field not in metadata_type:
                continue

            data_type = metadata_type[field]
            table = metadata_table[field]

            records_to_update = 0
            for orig_genome_id, value in metadata[field].items():

                try:
                    if float(value) and data_type in ['INT', 'INTEGER']:
                        # assume specified data type is correct and that we may need
                        # to cast floats to integers
                        value = str(int(float(value)))
                except:
                    pass

                if value.strip():
                    genome_id = str(orig_genome_id)
                    if genome_id.startswith('GCA_'):
                        genome_id = 'GB_' + genome_id
                    elif genome_id.startswith('GCF_'):
                        genome_id = 'RS_' + genome_id

                    if (not genome_list
                            or genome_id in genome_list
                            or orig_genome_id in genome_list):
                        data_to_commit.append((genome_id, value))
                        records_to_update += 1

            self.logger.info('Updating {} for {} genomes.'.format(
                field, records_to_update))

            # print(f'Committing {len(data_to_commit)} records to database for field {field} in table {table}.')
            # print(f'Example record: {data_to_commit[0]}')
            # print(f'Data type: {data_type}')

            gtdbimporter.import_metadata_to_db(table, field, data_type, data_to_commit)
            self.logger.info(f'Finished updating {field} for {records_to_update} genomes.')

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


    # set a field to NULL for every genome, in the transaction its new values are
    # written in; this asked [y/n] first, which a run with no terminal could not answer
    def set_field_to_null(self,metadata_table,field):
        self.logger.info('Setting {}.{} to NULL for every genome.'.format(metadata_table, field))
        q = ("UPDATE {} SET {} = NULL".format(
            metadata_table, field))
        self.temp_cur.execute(q)

    @one_transaction
    def update_ncbitax_db(self, organism_name_file,filtered_file,unfiltered_file, genome_list_file,do_not_null_field=False):
        """Add organism name to database."""
        gtdbimporter = GTDBImporter(self.temp_cur)

        genome_list = set()
        data_to_commit = []
        if genome_list_file:
            for line in open(genome_list_file):
                if '\t' in line:
                    genome_list.add(line.rstrip().split('\t')[0])
                else:
                    genome_list.add(line.rstrip().split(',')[0])

        # add full taxonomy string to database
        records_to_update =0
        for line in open(organism_name_file):
            line_split = line.strip().split('\t')

            gid = line_split[0]
            org_name = line_split[1]
            if genome_list_file and gid not in genome_list:
                continue

            data_to_commit.append((gid, org_name))
            records_to_update += 1

        if not do_not_null_field:
            self.set_field_to_null('metadata_ncbi', 'ncbi_organism_name')
        self.logger.info('Updating {} for {} genomes.'.format(
            'ncbi_organism_name', records_to_update))
        gtdbimporter.import_metadata_to_db('metadata_ncbi', 'ncbi_organism_name', 'TEXT', data_to_commit,
                                           unknown=SKIP)

        taxonomy = Taxonomy().read(filtered_file)
        data_filtered_to_commit = []
        records_to_update =0
        # add full taxonomy string to database
        for genome_id, taxa in taxonomy.items():
            if genome_id.startswith('GCA_'):
                genome_id = 'GB_' + genome_id
            elif genome_id.startswith('GCF_'):
                genome_id = 'RS_' + genome_id

            if genome_list_file and genome_id not in genome_list:
                continue
            taxa_str = ';'.join(taxa)
            data_filtered_to_commit.append((genome_id, taxa_str))
            records_to_update += 1
        if not do_not_null_field:
            self.set_field_to_null('metadata_taxonomy', 'ncbi_taxonomy')
        self.logger.info('Updating {} for {} genomes.'.format(
            'ncbi_taxonomy', records_to_update))
        gtdbimporter.import_metadata_to_db('metadata_taxonomy', 'ncbi_taxonomy', 'TEXT', data_filtered_to_commit,
                                           unknown=SKIP)

        # read taxonomy file
        unfiltered_taxonomy = Taxonomy().read(unfiltered_file)
        data_unfiltered_to_commit = []

        # add full taxonomy string to database
        records_to_update =0
        for genome_id, taxa in unfiltered_taxonomy.items():
            if genome_id.startswith('GCA_'):
                genome_id = 'GB_' + genome_id
            elif genome_id.startswith('GCF_'):
                genome_id = 'RS_' + genome_id

            if genome_list_file and genome_id not in genome_list:
                continue

            taxa_str = ';'.join(taxa)
            data_unfiltered_to_commit.append((genome_id, taxa_str))
            records_to_update += 1
        if not do_not_null_field:
            self.set_field_to_null('metadata_taxonomy', 'ncbi_taxonomy_unfiltered')
        self.logger.info('Updating {} for {} genomes.'.format(
            'ncbi_taxonomy_unfiltered', records_to_update))
        gtdbimporter.import_metadata_to_db('metadata_taxonomy', 'ncbi_taxonomy_unfiltered', 'TEXT', data_unfiltered_to_commit,
                                           unknown=SKIP)



