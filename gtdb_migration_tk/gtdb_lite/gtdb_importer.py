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

"""Writing one metadata field of many genomes, through the database's upsert().

Every metadata command writes through import_metadata_to_db(): update_metadata_db,
update_checkm_db, update_reps_db, update_ncbi_tax_db, update_propagated_tax and
add_taxonomy_to_database.

WHAT A FAILURE DOES

Until 0.1.47 any error here was printed and passed over. PostgreSQL then ignores
every statement of the transaction until it ends, and commit() on such a
transaction rolls it back without raising -- so a command whose write failed
finished as though it had written, and wrote nothing. An error is now raised, and
the command's transaction is rolled back (database_configuration.one_transaction).

GENOMES THE DATABASE DOES NOT HOLD

upsert() finds a genome by genomes.id_at_source, and gives a genome it does not
find a NULL id, which every metadata table refuses: it raises, failing the whole
field, so one genome missing from the database lost a field's values for every
genome. The genomes are therefore checked here, before upsert() is called, and
what happens to those it does not hold is the caller's to say. Input meant to
be the release's own genomes -- a metadata table, CheckM's results, a taxonomy
-- is refused (REFUSE): a genome the database does not hold says the input is
of another release, or that update_db has not been run. Input that is a superset
by design, such as NCBI's organism names, which cover every assembly NCBI holds,
has those genomes skipped and counted (SKIP).

upsert() is handed a genome as the database holds it in genomes.id_at_source,
GCA_000003645.1, and a genome is taken named that way or as GTDB names it, with a
leading GB_ or RS_ (GB_GCA_000003645.1), which is taken off. Everything before
the first '_' was taken off until 0.1.47, whatever it was, so a genome named
without its prefix was cut to 000003645.1, which is no genome, and its field
failed; add_taxonomy_to_database is documented as taking such names.
"""

import csv
import logging
import os
from typing import Iterable, List, Optional, Set, Tuple

from gtdb_migration_tk.biolib_lite.logger import log_directory

# what import_metadata_to_db() does with a genome the database does not hold
REFUSE = 'refuse'
SKIP = 'skip'

# the prefixes GTDB names NCBI's genomes with, which genomes.id_at_source is without
GTDB_PREFIXES = ('GB_', 'RS_')

# how many of the genomes the database does not hold an error names; all of them
# are written to a file beside the log
UNKNOWN_EXAMPLES = 10


class UnknownGenomesError(ValueError):
    """Genomes to be written are not in the database."""


def id_at_source(genome_id: Optional[str]) -> Optional[str]:
    """What upsert() is handed for a genome: the ID without a GTDB prefix.

    Parameters
    ----------
    genome_id : str or None
        The genome as GTDB names it, e.g. GB_GCA_000003645.1, or as NCBI does,
        e.g. GCA_000003645.1.

    @return: e.g. GCA_000003645.1, the genome's genomes.id_at_source, a leading
             GB_ or RS_ taken off and anything else left as it is; None where
             the ID is None or empty.
    """

    if not genome_id:
        return None
    for prefix in GTDB_PREFIXES:
        if genome_id.startswith(prefix):
            return genome_id[len(prefix):]
    return genome_id


class GTDBImporter(object):

    def __init__(self, temp_cur):
        """Initialization.

        Parameters
        ----------
        temp_cur : psycopg2.extensions.cursor
            Cursor of the command's transaction.

        @return: None
        """

        self.temp_cur = temp_cur
        self.logger = logging.getLogger('timestamp')
        self._genomes: Optional[Set[str]] = None

    def genomes(self) -> Set[str]:
        """Every genome the database holds, as genomes.id_at_source names it.

        Read once for each importer, the first time it is asked.

        @return: e.g. {'GCA_000003645.1', 'GCF_000003135.1', ...}.
        """

        if self._genomes is None:
            self.temp_cur.execute('SELECT id_at_source FROM genomes')
            self._genomes = {gid for (gid,) in self.temp_cur.fetchall()}
        return self._genomes

    def import_metadata_to_db(self, table: str, field: str, typemeta: str,
                              data_list: Iterable[Tuple[str, str]],
                              unknown: str = REFUSE) -> int:
        """Write one field of metadata for a list of genomes.

        Parameters
        ----------
        table : str
            Table holding the field, e.g. metadata_taxonomy.
        field : str
            The field, e.g. gtdb_domain.
        typemeta : str
            The field's type in the database, e.g. TEXT or FLOAT.
        data_list : iterable of (str, str)
            Each genome, as GTDB names it (GB_GCA_000003645.1) or as NCBI does
            (GCA_000003645.1), and its value.
        unknown : str
            REFUSE or SKIP: what is done with genomes the database does not hold.

        @return: the number of genomes skipped as not in the database.

        Raises
        ------
        UnknownGenomesError
            unknown is REFUSE and a genome is not in the database; nothing is
            written.
        """

        if unknown not in (REFUSE, SKIP):
            raise ValueError('unknown must be {} or {}, not {}.'.format(REFUSE, SKIP, unknown))

        genomes = self.genomes()
        source_ids: List[str] = []
        values: List[str] = []
        not_held: List[str] = []
        for genome_id, value in data_list:
            source_id = id_at_source(genome_id)
            if source_id in genomes:
                source_ids.append(source_id)
                values.append(value)
            else:
                not_held.append(str(genome_id))

        if not_held:
            written = self.write_unknown_genomes(table, field, not_held)
            if unknown == REFUSE:
                raise UnknownGenomesError(
                    '{:,} genome(s) to be written to {}.{} are not in the database, '
                    'e.g. {}; the database holds another release, or update_db has '
                    'not been run. Every one is in {}.'.format(
                        len(not_held), table, field,
                        ', '.join(not_held[:UNKNOWN_EXAMPLES]), written))
            self.logger.info('Skipped {:,} genome(s) not in the database for {}.{} '
                             '(listed in {}).'.format(len(not_held), table, field, written))

        if source_ids:
            self.temp_cur.execute('SELECT upsert(%s::regclass, %s, %s, %s, %s)',
                                  (table, field, typemeta, source_ids, values))

        return len(not_held)

    @staticmethod
    def write_unknown_genomes(table: str, field: str, genome_ids: List[str]) -> str:
        """Write the genomes of a field the database does not hold beside the log.

        Parameters
        ----------
        table, field : str
            What was being written.
        genome_ids : list of str
            The genomes, as they were given.

        @return: the file written, e.g. <log dir>/unknown_genomes.metadata_ncbi.ncbi_organism_name.tsv.
        """

        path = os.path.join(log_directory(),
                            'unknown_genomes.{}.{}.tsv'.format(table, field))
        with open(path, 'w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
            writer.writerow(['genome_id'])
            writer.writerows([gid] for gid in genome_ids)
        return path
