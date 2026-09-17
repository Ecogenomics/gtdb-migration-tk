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
directory_manager.py -- index the NCBI mirror.

list_genomes works on the local NCBI mirror rather than on the GTDB release
directories, and exists because the mirror is a deep tree that is expensive to
walk. NCBI nests a genome by the nine digits of its accession,
three at a time, so GCA_000001405.28 lives at

    <mirror>/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/

and finding a genome means four levels of os.listdir(). Walking that once and
writing the result to a file is what makes the rest of the toolkit tractable.

That file -- the genome_dirs file written by generate_genome_dir_file() and
consumed by nearly every other command -- is the lingua franca between steps:
one line per genome, accession, absolute path, canonical accession. Downstream
readers split on tabs and ignore any further columns (see
update_genomes.UpdateGenomes.load_genome_dirs), so columns may be appended but
never reordered.

The file describes one release, so the walk is filtered by the table
select_genomes writes: a directory holding a genome that table does not list is
passed over, the tree being in general a GTDB release accumulated over several
cycles. Whether the tree holds what it should is a separate question, answered
by ncbi_genome_sync --verify against the same table -- it has the manifests, so
it can check the files and not merely the directories.

The mirror is held on NFS, where the cost of the walk is one network round trip
per directory read rather than any computation, so the walk is spread over
threads (--cpus): the round trips then overlap instead of being paid one after
another. Measured on release220 that is worth roughly 8x, and it saturates by
about four threads, so raising --cpus far beyond its default buys nothing.

Pruning the mirror of genomes a release has dropped used to be a second
command here (clean_ftp). It is now the first step of ncbi_genome_sync, which
holds the selection that defines what the mirror should contain and so can
remove what it does not list before fetching anything; see REMOVAL in
ncbi_genome_sync.py.
"""

import os
import logging
from concurrent.futures import ThreadPoolExecutor
from typing import List, Set, Tuple

from tqdm import tqdm

from gtdb_migration_tk.ncbi_utils import assembly_accession
from gtdb_migration_tk.biolib_lite.common import canonical_gid
from gtdb_migration_tk.ncbi_utils import read_assembly_summary


class DirectoryManager(object):
    """Index the genomes held on the NCBI mirror.

    Backs 'list_genomes', which writes the genome_dirs file the rest of the
    toolkit reads: where on disk each genome of the new release is held.
    """

    def __init__(self) -> None:
        """Instantiate the manager.

        @return: None
        """

        self.logger = logging.getLogger('timestamp')

    def read_selected_genomes(self, gtdb_selected_genomes: str) -> Set[str]:
        """Read the accessions making up the new release.

        The release is defined by the table select_genomes writes, which is read
        through ncbi_utils like every other table of this shape -- so it is taken
        gzipped or plain, and its columns are located by name.

        Parameters
        ----------
        gtdb_selected_genomes : str
            Table of the genomes selected for the new release.

        @return: set of accessions present in the new release.
        """

        return {accession for accession, in
                read_assembly_summary(gtdb_selected_genomes, 'assembly_accession')}

    def _list_subdirs(self, path: str) -> List[str]:
        """List the subdirectories of a directory.

        Every level of the mirror holds files alongside its directories, and
        only the directories are part of the accession layout. This is also the
        unit the walk is built from: one directory read, which on the NFS
        mounted mirror is one network round trip.

        Parameters
        ----------
        path : str
            Directory to list.

        @return: sorted paths of the subdirectories of path.
        """

        subdirs = []
        for name in os.listdir(path):
            full_path = os.path.join(path, name)
            if os.path.isdir(full_path):
                subdirs.append(full_path)

        return sorted(subdirs)

    def _genomes_under_triplet(self, second_triplet_dir: str) -> List[Tuple[str, str]]:
        """Find every genome directory beneath a second digit triplet.

        This is the unit of work handed to a thread. It covers the two deepest
        levels of the layout, which is where all but a few thousand of the
        mirror's directories are, so splitting the walk here is what puts the
        round trips in flight concurrently.

        Parameters
        ----------
        second_triplet_dir : str
            Directory of a second digit triplet, e.g. <mirror>/GCF/005/435.

        @return: list of (accession, genome directory) pairs beneath it.
        """

        genomes = []
        for third_triplet_dir in self._list_subdirs(second_triplet_dir):
            for genome_dir in self._list_subdirs(third_triplet_dir):
                # the leaf directory is named <accession>_<assembly name>,
                # e.g. GCA_000001405.28_GRCh38.p13; cutting at the first
                # underscore after the GCA_/GCF_ prefix leaves the versioned
                # accession and drops the assembly name, which may itself
                # contain underscores
                complete_name = os.path.basename(genome_dir)
                accession = assembly_accession(complete_name)
                genomes.append((accession, genome_dir))

        return genomes

    def generate_genome_dir_file(self,
                                 database_dir: str,
                                 output_file: str,
                                 gtdb_selected_genomes: str,
                                 cpus: int = 8) -> None:
        """Create file indicating directory of each genome.

        Walks the four levels of the layout described in the module docstring and
        writes one line per genome of the new release, giving where it is held.

        Only genomes the selection lists are written. A directory holding anything
        else is passed over: the tree may be a GTDB release built over several
        cycles, and the genome_dirs file describes THIS release. Accessions carry
        their version, so a genome revised since the tree was built does not match
        and is not written.

        Whether the tree holds what it should is not decided here. That is
        ncbi_genome_sync's --verify, which compares a mirror against the same
        selection and has the manifests to check the files as well as the
        directories. What this reports is a count, so a short file is noticed.

        Parameters
        ----------
        database_dir : str
            Root of the tree, holding GCA/ and GCF/ subdirectories.
        output_file : str
            Genome directory file to write (accession, absolute path, canonical accession).
        gtdb_selected_genomes : str
            Table of the genomes selected for the new release.
        cpus : int
            Number of threads walking the tree concurrently.

        @return: None
        """

        selected = self.read_selected_genomes(gtdb_selected_genomes)
        self.logger.info('Release: {:,} genomes.'.format(len(selected)))

        # resolved once, so every path built below is already absolute and the
        # walk need not call os.path.abspath() for each of the millions of genomes
        database_dir = os.path.abspath(database_dir)

        first_triplet_dirs = []
        for code, archive in [('GCA', 'GenBank'), ('GCF', 'RefSeq')]:
            code_dir = os.path.join(database_dir, code)

            # a mirror may hold only one of the two archives, so a missing
            # GCA/ or GCF/ is normal rather than an error
            if not os.path.isdir(code_dir):
                self.logger.info('Skipping {}, no {} directory.'.format(archive, code_dir))
                continue

            first_triplet_dirs.extend(self._list_subdirs(code_dir))

        # accession -> genome directory, for the two reports below
        genomes_on_disk = {}

        with ThreadPoolExecutor(max_workers=cpus) as pool:
            # the second triplets are all enumerated before the walk starts so
            # that the progress bar has a real total to count against. There are
            # tens of thousands of them holding a comparable number of genomes
            # each, whereas the 90 first triplets differ in size by more than an
            # order of magnitude, which made both the bar and its estimate of
            # the time remaining next to useless
            second_triplets_by_first = list(pool.map(self._list_subdirs, first_triplet_dirs))

            progress = tqdm(total=sum(len(s) for s in second_triplets_by_first),
                            bar_format='{desc:<12.12}{percentage:3.0f}%|{bar:20}{r_bar}')

            with open(output_file, 'w') as fout:
                # one first triplet is submitted at a time, so at most a thousand
                # or so units of work are queued however large the mirror is
                for second_triplet_dirs in second_triplets_by_first:
                    # map() yields in submission order, so each result can be
                    # paired back up with its directory, and the genome_dirs file
                    # does not depend on the order the threads happen to finish in
                    for second_triplet_dir, genomes in zip(
                            second_triplet_dirs,
                            pool.map(self._genomes_under_triplet, second_triplet_dirs)):
                        progress.set_description_str(
                            '/'.join(second_triplet_dir.split(os.sep)[-3:]))
                        progress.update()

                        for accession, genome_dir in genomes:
                            if accession not in selected:
                                continue
                            genomes_on_disk[accession] = genome_dir
                            fout.write('{}\t{}\t{}\n'.format(
                                accession, genome_dir, canonical_gid(accession)))

            progress.close()

        self.logger.info('Indexed {:,} of {:,} genomes in {}'.format(
            len(genomes_on_disk), len(selected), output_file))

        # not a verification -- that is ncbi_genome_sync --verify -- but a file
        # quietly short of the release it claims to describe is worth a line
        missing = len(selected) - len(genomes_on_disk)
        if missing:
            self.logger.warning('{:,} genomes of the release have no directory under {} and '
                                'are not in {}; ncbi_genome_sync --verify says which.'.format(
                                    missing, database_dir, output_file))
