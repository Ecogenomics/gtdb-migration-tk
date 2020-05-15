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

import psycopg2 as pg


class GenomeDatabaseConnectionFTPUpdate(object):

    def __init__(self, hostname, user, password, db_name):
        self.conn = None
        self.hostname = hostname
        self.user = user
        self.password = password
        self.db_name = db_name

    # Opens a connection to the PostgreSQL database
    #
    # Returns:
    #   No return value.
    def MakePostgresConnection(self):
        conn_string = "dbname=%s user=%s host=%s password=%s" % (
            self.db_name, self.user,
            self.hostname, self.password
        )
        self.conn = pg.connect(conn_string)

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
