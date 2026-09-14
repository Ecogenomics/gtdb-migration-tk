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

The Pfam and TIGRFAM releases GTDB annotates against change every couple of
years. Their version numbers appear in the name of the directory each genome's
search results are written to, and in the name of every file within it, so a
version number left behind in one place produces genome directories whose
symlinks point at annotations that are not there. Setting the two constants
below is the whole of that change.

Note that the search results themselves are produced by the hmmsearch and
top_hit commands, which take the marker version on the command line rather than
from here (see --folder_suffix). The values below describe what the genome
directories are expected to contain once those commands have run, and they must
agree with what was passed to them.

The SILVA and LTP releases the rRNA genes of a genome are classified against
change on the same footing, and for the same reason: their version numbers name
the directory the results are written to inside each genome directory. Those
directories are taken on the command line too, by rna_silva and rna_ltp
(--silva_version, --ltp_version).
"""

# Version of each marker database GTDB currently annotates against.
PFAM_VERSION = '33.1'
TIGRFAM_VERSION = '15.0'

# Version of each rRNA database GTDB currently classifies against. These name the
# directory holding the results, as several releases of a database may sit side
# by side in a genome directory.
SILVA_VERSION = '138.2'
LTP_VERSION = '10_2024'

# Directory within a genome directory holding its search results. GTDB searches
# the reduced ("lite") marker sets, hence the suffix.
PFAM_MARKER_DIR = 'pfam_{}_lite'.format(PFAM_VERSION)
TIGRFAM_MARKER_DIR = 'tigrfam_{}_lite'.format(TIGRFAM_VERSION)

# Suffixes of the search result files held in those directories.
PFAM_EXT = '_{}.tsv.gz'.format(PFAM_MARKER_DIR)
PFAM_TOPHIT_EXT = '_{}_tophit.tsv.gz'.format(PFAM_MARKER_DIR)
TIGRFAM_EXT = '_{}.tsv.gz'.format(TIGRFAM_MARKER_DIR)
TIGRFAM_TOPHIT_EXT = '_{}_tophit.tsv.gz'.format(TIGRFAM_MARKER_DIR)
TIGRFAM_OUT_EXT = '_{}.out.gz'.format(TIGRFAM_MARKER_DIR)

# Names of the symlinks pointing at those files. These carry no version number
# by design: downstream code opens a genome's Pfam hits without knowing which
# Pfam release produced them, so only the symlink target changes at a release.
PFAM_SYMLINK_EXT = '_pfam_lite.tsv.gz'
PFAM_TOPHIT_SYMLINK_EXT = '_pfam_lite_tophit.tsv.gz'
TIGRFAM_SYMLINK_EXT = '_tigrfam_lite.tsv.gz'
TIGRFAM_TOPHIT_SYMLINK_EXT = '_tigrfam_lite_tophit.tsv.gz'
TIGRFAM_OUT_SYMLINK_EXT = '_tigrfam_lite.out.gz'

# Uncompressed HMMER output to gzip when a genome is taken from the FTP site.
# These name the full marker sets rather than the lite ones, as that is what a
# genome directory carries before the lite search results are generated.
HMMER_EXTS_TO_GZIP = (
    '_pfam_{}.tsv'.format(PFAM_VERSION),
    '_pfam_{}_tophit.tsv'.format(PFAM_VERSION),
    '_tigrfam_{}.out'.format(TIGRFAM_VERSION),
    '_tigrfam_{}.tsv'.format(TIGRFAM_VERSION),
    '_tigrfam_{}_tophit.tsv'.format(TIGRFAM_VERSION),
)

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
