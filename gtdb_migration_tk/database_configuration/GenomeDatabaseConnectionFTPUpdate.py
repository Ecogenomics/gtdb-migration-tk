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

import functools
import io
from typing import Iterable


# The temporary table reset_unwritten() puts the genomes just written in.
WRITTEN_TABLE = 'gtdb_written_genomes'

# What reset_unwritten() sets to NULL: a genome holding a value that this run gave
# none. A genome given a value is rewritten by upsert() whatever it held, so
# setting it to NULL first only rewrote the row once more.
RESET_UNWRITTEN = ('UPDATE {table} AS m SET {field} = NULL FROM genomes g '
                   'WHERE g.id = m.id AND m.{field} IS NOT NULL '
                   'AND NOT EXISTS (SELECT 1 FROM ' + WRITTEN_TABLE + ' w WHERE w.id_at_source = g.id_at_source)')


def reset_unwritten(cur, table: str, field: str, written: Iterable[str]) -> int:
    """Set a field to NULL for every genome holding a value this run did not write.

    Done once the field is written, in the same transaction, which leaves what
    setting it to NULL for every genome first did -- the new values, and NULL
    for every other genome -- without rewriting the rows of the genomes given a
    value. PostgreSQL writes a new version of every row an UPDATE touches: the
    NULLs alone rewrote metadata_ncbi, 1.3 GB, for each field, and upsert()
    then rewrote it again, each version kept to the end of the transaction.

    Parameters
    ----------
    cur : cursor
        The cursor of the caller's transaction.
    table, field : str
        e.g. metadata_ncbi, ncbi_organism_name.
    written : iterable of str
        The genomes given a value, as genomes.id_at_source names them.

    @return: the number of genomes set to NULL.
    """

    cur.execute('CREATE TEMPORARY TABLE IF NOT EXISTS {} (id_at_source TEXT) ON COMMIT DROP'.format(WRITTEN_TABLE))
    cur.execute('TRUNCATE {}'.format(WRITTEN_TABLE))
    cur.copy_expert('COPY {} (id_at_source) FROM STDIN'.format(WRITTEN_TABLE),
                    io.StringIO(''.join(genome + '\n' for genome in written)))
    cur.execute('ANALYZE {}'.format(WRITTEN_TABLE))
    cur.execute(RESET_UNWRITTEN.format(table=table, field=field))
    return cur.rowcount


def one_transaction(method):
    """A manager's command, whose writes are committed together or not at all.

    The metadata commands committed as they went: the fields set to NULL, then
    each field or taxonomic rank written. A run that failed part way left the
    database with some fields of the new release and some of the old, or a field
    NULL for every genome, its new values never written. The method's writes are
    now committed once, when it returns, and rolled back when it raises -- a
    sys.exit() among them.

    The manager holds its connection as temp_con, a GenomeDatabaseConnectionFTPUpdate,
    and the method makes no commit of its own.

    Parameters
    ----------
    method : function
        A method of the manager.

    @return: the method, wrapped.
    """

    @functools.wraps(method)
    def in_one_transaction(self, *args, **kwargs):
        try:
            result = method(self, *args, **kwargs)
        except BaseException:
            self.temp_con.rollback()
            raise
        self.temp_con.commit()
        return result

    return in_one_transaction


class GenomeDatabaseConnectionFTPUpdate(object):

    def __init__(self, database):
        """Initialization.

        Parameters
        ----------
        database : dict
            How to reach the database, as utils.common.database_keywords() gives it.

        @return: None
        """

        self.conn = None
        self.database = database

    # Opens a connection to the PostgreSQL database. The libpq keywords are
    # handed over as they are rather than written into a connection string,
    # where a password holding a space or a quote was mangled.
    #
    # Returns:
    #   No return value.
    def MakePostgresConnection(self):
        import psycopg2 as pg

        self.conn = pg.connect(**self.database)

    # Function: ClosePostgresConnection
    # Closes an open connection to the PostgreSQL database.
    #
    # Returns:
    #   No return value.
    def ClosePostgresConnection(self):
        if self.IsPostgresConnectionActive():
            self.conn.close()
            self.conn = None

    # Function: IsPostgresConnectionActive
    # Check if the connection to the PostgreSQL database is active.
    #
    # Returns:
    #   True if connection is active, False otherwise
    def IsPostgresConnectionActive(self):
        if self.conn is not None:
            cur = self.conn.cursor()
            try:
                cur.execute("SELECT count(*) from users")
            except:
                return False
            cur.close()
            return True
        else:
            return False

    # Convenience methods to the pg connection
    def commit(self):
        return self.conn.commit()

    def rollback(self):
        return self.conn.rollback()

    def cursor(self):
        return self.conn.cursor()
