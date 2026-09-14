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
config.py -- settings that change from one GTDB release to the next.

The Pfam and TIGRFAM releases GTDB annotates against, and the SILVA and LTP
releases it classifies rRNA genes against, change every couple of years. Each
version number names the directory inside a genome directory that the results
of that database are written to, so a version left behind in one place produces
a release whose commands read and write different directories. The version
constants below are the one place each is declared.

The marker versions are the defaults of the --folder_suffix option of hmmsearch
and top_hit, through MARKER_FOLDER_SUFFIX, so those commands write the
directories this file declares unless told otherwise. The rRNA versions name the
directories in GTDB_DERIVED_DIRS_TO_COPY, the derived data carried across from
the previous release when a genome's sequence is unchanged; the rna_silva and
rna_ltp commands still take their version on the command line
(--silva_version, --ltp_version), and the value passed there must agree with
this file for the carried-over directories to be the ones later steps read.
"""

# Version of each marker database GTDB currently annotates against.
PFAM_VERSION = '33.1'
TIGRFAM_VERSION = '15.0'

# Version of each rRNA database GTDB currently classifies against. These name the
# directory holding the results, as several releases of a database may sit side
# by side in a genome directory.
SILVA_VERSION = '138.2'
LTP_VERSION = '10_2024'

# Default --folder_suffix of hmmsearch and top_hit for each marker database: the
# directory they write is pfam_<suffix>/ or tigrfam_<suffix>/ inside prodigal/,
# and the files in it carry the same suffix. GTDB searches the reduced ("lite")
# marker sets, hence _lite. Keyed by the value of --db.
MARKER_FOLDER_SUFFIX = {
    'pfam': '{}_lite'.format(PFAM_VERSION),
    'tigrfam': '{}_lite'.format(TIGRFAM_VERSION),
}

# Directories of derived data within a genome directory, carried across when a
# genome is taken from the previous release rather than from NCBI. A genome is
# only carried across when its FASTA files are unchanged, which is what makes the
# derived data still valid; copying it is what saves the release from calling
# genes and searching rRNA databases again for every genome GTDB already holds.
GTDB_DERIVED_DIRS_TO_COPY = (
    'prodigal',
    'rna_silva_{}'.format(SILVA_VERSION),
    'trna',
    'rna_ltp_{}'.format(LTP_VERSION),
)
