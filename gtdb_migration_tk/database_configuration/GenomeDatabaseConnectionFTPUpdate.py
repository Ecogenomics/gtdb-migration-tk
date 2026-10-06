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
