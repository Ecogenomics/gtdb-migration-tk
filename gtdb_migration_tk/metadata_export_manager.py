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

"""Export the metadata of every genome of the database: export_metadata.

This is what 'gtdb metadata export' did, the export the release's metadata files
and gtdb_species_clusters' input were taken by: every column of metadata_view
but its internal id, one row a genome, as CSV or as tab-separated text writing
NULL as 'none'. It differs from it in three ways.

IT ONLY READS

'gtdb metadata export' first ran assignToRepresentative(), which aligned the
markers of every genome whose gtdb_representative was NULL, gave it a
representative by AAI and committed, so that an export changed the database.
gtdb_representative is the toolkit's to write: update_propagated_tax writes it
for every genome from the previous release, and update_reps_db for every genome
from the release's species clusters (r237: none NULL). The export is one
read-only transaction, REPEATABLE READ so that what it says of the genomes and
what it writes are of the same snapshot. A genome without a gtdb_representative
is counted in a WARNING and listed beside the log (NO_REPRESENTATIVE_NAME), and
exported as it is.

THE ROWS ARE IN ORDER

The rows were in whatever order metadata_view gave, which may change between two
exports of the same database. They are ordered by accession, and a gzipped
export is written without a time in its header (utils.common.open_gzip_text()),
so two exports of an unchanged database are the same file.

A FAILED EXPORT LEAVES NO FILE

The export is written beside the output and renamed onto it once complete; one
that fails part way is removed, leaving any earlier export of the same name as it
was. A truncated table that reads as a complete one is worse than none.
"""

import logging
import os
import time
from typing import List, Tuple

from psycopg2.extensions import ISOLATION_LEVEL_REPEATABLE_READ

from gtdb_migration_tk.biolib_lite.common import make_sure_path_exists
from gtdb_migration_tk.biolib_lite.logger import log_directory
from gtdb_migration_tk.database_configuration.GenomeDatabaseConnectionFTPUpdate import \
    GenomeDatabaseConnectionFTPUpdate
from gtdb_migration_tk.utils.common import GZIP_SUFFIX, open_gzip_text

# the view exported, the columns of it that are not, and the column the rows are ordered by
METADATA_VIEW = 'metadata_view'
NOT_EXPORTED = ('id',)
ORDER_BY = 'accession'

# --format: COPY's options for each, as 'gtdb metadata export' gave them
FORMATS = {
    'csv': "FORMAT csv, HEADER",
    'tab': "FORMAT csv, HEADER, DELIMITER E'\\t', NULL 'none'",
}
DEFAULT_FORMAT = 'tab'

# each genome of the database without a gtdb_representative, as GTDB names it; a
# genome with no metadata_taxonomy row has none
NO_REPRESENTATIVE = ("SELECT gs.external_id_prefix || '_' || g.id_at_source "
                     'FROM genomes g JOIN genome_sources gs ON gs.id = g.genome_source_id '
                     'LEFT JOIN metadata_taxonomy m ON m.id = g.id '
                     'WHERE m.gtdb_representative IS NULL ORDER BY 1')
NUM_GENOMES = 'SELECT count(*) FROM genomes'

# beside the log: the genomes without a gtdb_representative
NO_REPRESENTATIVE_NAME = 'export_no_gtdb_representative.tsv'

# what an export is written to until it is complete
PARTIAL_SUFFIX = '.partial'

# how many genomes a line of the log names
EXAMPLES = 10


def quoted(name: str) -> str:
    """A column name quoted for SQL.

    @return: e.g. '"gtdb_representative"'.
    """

    return '"' + name.replace('"', '""') + '"'


def export_query(columns: List[str], out_format: str) -> str:
    """The COPY statement exporting the columns of metadata_view, ordered by accession.

    Parameters
    ----------
    columns : list of str
        The columns of metadata_view, in its order; NOT_EXPORTED are left out.
    out_format : str
        A key of FORMATS.

    @return: the statement.

    Raises
    ------
    ValueError
        An unknown format, or a view without the column the rows are ordered by.
    """

    if out_format not in FORMATS:
        raise ValueError('The format must be one of {}, not {!r}.'.format(', '.join(FORMATS), out_format))
    if ORDER_BY not in columns:
        raise ValueError('{} has no column {!r} to order the rows by.'.format(METADATA_VIEW, ORDER_BY))
    exported = ', '.join(quoted(column) for column in columns if column not in NOT_EXPORTED)
    return 'COPY (SELECT {} FROM {} ORDER BY {}) TO STDOUT WITH ({})'.format(
        exported, METADATA_VIEW, quoted(ORDER_BY), FORMATS[out_format])


def open_output(path: str, output_file: str):
    """Open the file an export is written to, gzipped where the output's name ends in .gz.

    Parameters
    ----------
    path : str
        The file written: output_file + PARTIAL_SUFFIX until it is complete, so
        it is the output's name that says whether to gzip.
    output_file : str
        The export's name.

    @return: a text stream, to be used as a context manager.
    """

    if output_file.endswith(GZIP_SUFFIX):
        return open_gzip_text(path)
    return open(path, 'w', encoding='utf-8', newline='')


class MetadataExportManager(object):
    """export_metadata: the metadata of every genome of the database, written to a file."""

    def __init__(self, database):
        """Connect, in a read-only REPEATABLE READ transaction.

        Parameters
        ----------
        database : dict
            How to reach the database, as utils.common.database_keywords() gives it.

        @return: None
        """

        self.logger = logging.getLogger('timestamp')
        self.temp_con = GenomeDatabaseConnectionFTPUpdate(database)
        self.temp_con.MakePostgresConnection()
        # set before the first statement opens the transaction
        self.temp_con.conn.set_session(isolation_level=ISOLATION_LEVEL_REPEATABLE_READ, readonly=True)
        self.temp_cur = self.temp_con.cursor()

    def columns(self) -> List[str]:
        """The columns of metadata_view, in its order.

        @return: e.g. ['id', 'accession', ...].
        """

        self.temp_cur.execute('SELECT * FROM {} LIMIT 0'.format(METADATA_VIEW))
        return [column[0] for column in self.temp_cur.description]

    def check_representatives(self) -> Tuple[int, int]:
        """Count, in a WARNING, the genomes without a gtdb_representative, each listed beside the log.

        @return: (genomes without one, genomes of the database).
        """

        self.temp_cur.execute(NUM_GENOMES)
        num_genomes = self.temp_cur.fetchone()[0]
        self.temp_cur.execute(NO_REPRESENTATIVE)
        missing = [row[0] for row in self.temp_cur.fetchall()]
        if missing:
            report = os.path.join(log_directory(), NO_REPRESENTATIVE_NAME)
            with open(report, 'w') as handle:
                handle.write('genome_id\n')
                for genome_id in missing:
                    handle.write(genome_id + '\n')
            self.logger.warning(
                '{:,} of the {:,} genomes of the database have no gtdb_representative, e.g. {}: '
                'update_reps_db has not been run since they were added. They are exported '
                'without one; gtdb metadata export gave each a representative by AAI first. '
                'Every one is in {}.'.format(len(missing), num_genomes, ', '.join(missing[:EXAMPLES]), report))
        else:
            self.logger.info('Every one of the {:,} genomes of the database has a gtdb_representative.'.format(
                num_genomes))
        return len(missing), num_genomes

    def export_metadata(self, output_file: str, out_format: str = DEFAULT_FORMAT) -> int:
        """Write the metadata of every genome of the database to output_file.

        Parameters
        ----------
        output_file : str
            The file written, gzipped where it ends in .gz; its directory is made.
        out_format : str
            'tab' (tab-separated, NULL as 'none') or 'csv' (NULL empty).

        @return: the number of genomes exported.
        """

        try:
            columns = self.columns()
            query = export_query(columns, out_format)
            self.check_representatives()

            make_sure_path_exists(os.path.dirname(output_file))
            partial = output_file + PARTIAL_SUFFIX
            started = time.time()
            try:
                with open_output(partial, output_file) as handle:
                    self.temp_cur.copy_expert(query, handle)
                os.replace(partial, output_file)
            except BaseException:
                if os.path.exists(partial):
                    os.remove(partial)
                raise
            exported = self.temp_cur.rowcount
            self.logger.info('Exported the metadata of {:,} genomes, {} columns of {}, as {}, to {} in {:.0f} s.'.format(
                exported, len([c for c in columns if c not in NOT_EXPORTED]), METADATA_VIEW, out_format,
                output_file, time.time() - started))
            return exported
        finally:
            # nothing to commit: the transaction only read
            self.temp_con.rollback()
