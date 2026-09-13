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
directory_manager.py -- index the NCBI mirror, and prune genomes it has dropped.

Both commands here work on the local NCBI mirror rather than on the GTDB
release directories, and both exist because the mirror is a deep tree that is
expensive to walk. NCBI nests a genome by the nine digits of its accession,
three at a time, so GCA_000001405.28 lives at

    <mirror>/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/

and finding a genome means four levels of os.listdir(). Walking that once and
writing the result to a file is what makes the rest of the toolkit tractable.

That file -- the genome_dirs file written by generate_genome_dir_file() and
consumed by nearly every other command -- is the lingua franca between steps:
one line per genome, accession, absolute path, canonical accession. Downstream
readers split on tabs and ignore any further columns (see
ftp_manager._populate_genomes_dict), so columns may be appended but never
reordered.

Walking the mirror is also the only way to find out what it does not hold, so
generate_genome_dir_file() checks the directories it finds against the genome
lists for the new release. A genome on those lists with no directory was never
mirrored, and a directory for a genome absent from them is left over from an
earlier release; either one silently distorts the release built from the file,
so both are reported rather than left for a later command to trip over, in
files named after the genome_dirs file with '-missing' and '-extra' appended.

The mirror is held on NFS, where the cost of the walk is one network round trip
per directory read rather than any computation, so the walk is spread over
threads (--cpus): the round trips then overlap instead of being paid one after
another. Measured on release220 that is worth roughly 8x, and it saturates by
about four threads, so raising --cpus far beyond its default buys nothing.

clean_ftp() is the other direction: NCBI suppresses assemblies between
releases, and the mirror keeps serving them until something deletes them. It
compares the genome_dirs file describing what the mirror currently holds
against the genome lists for the new release, removes what has been dropped,
and writes both differences to a report so the deletion can be reviewed after
the fact. It deletes real data, so it reports before and while it acts, never
afterwards.
"""

import os
import logging
import shutil
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import List, Optional, Set, Tuple, Union

from tqdm import tqdm

from gtdb_migration_tk.biolib_lite.common import make_sure_path_exists, canonical_gid
from gtdb_migration_tk.biolib_lite.taxonomy import Taxonomy


class DirectoryManager(object):
    """Index the genomes held on the NCBI mirror, and delete the ones NCBI has dropped.

    Backs two commands: 'list_genomes', which writes the genome_dirs file the
    rest of the toolkit reads, and 'clean_ftp', which prunes the mirror of
    genomes missing from the new release. Both are told what the new release
    contains, and both compare that against what the mirror actually holds.
    """

    def __init__(self) -> None:
        """Instantiate the manager.

        @return: None
        """

        self.logger = logging.getLogger('timestamp')

    def read_genome_lists(self, new_list_genomes: List[str]) -> Set[str]:
        """Read the accessions making up the new release.

        The release may be described by more than one file; the accessions from
        all of them form a single set. Both commands in this module take the
        same lists and read them the same way, so they read them here.

        Parameters
        ----------
        new_list_genomes : list of str
            Files indicating the Gid present in the new release.

        @return: set of accessions present in the new release.
        """

        accessions = set()
        for new_genome_file in new_list_genomes:
            with open(new_genome_file, 'r') as ngf:
                for line in ngf:
                    if line.startswith('#'):
                        continue
                    accessions.add(line.strip().split('\t')[0])

        return accessions

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
                accession = complete_name[0:complete_name.find('_', 4)]
                genomes.append((accession, genome_dir))

        return genomes

    def generate_genome_dir_file(self,
                                 database_dir: str,
                                 output_file: str,
                                 new_list_genomes: List[str],
                                 cpus: int = 8) -> None:
        """Create file indicating directory of each genome.

        Walks the four levels of the NCBI mirror layout described in the module
        docstring, writes one line per genome directory found, and reports any
        genome in the new release with no directory on the mirror, or directory
        on the mirror holding a genome absent from the new release. Those two
        sets are written beside the genome_dirs file, suffixed '-missing' and
        '-extra'.

        Accessions are compared including their version, so a genome NCBI has
        revised since the mirror was last synced is reported twice: once as
        missing at its new version, once as unexpected at its old one.

        All three files are written whatever the comparison finds, as the
        discrepancies are for the operator to judge.

        Parameters
        ----------
        database_dir : str
            Root of the local NCBI mirror, holding GCA/ and GCF/ subdirectories.
        output_file : str
            Genome directory file to write (accession, absolute path, canonical accession).
        new_list_genomes : list of str
            Files indicating the Gid present in the new release.
        cpus : int
            Number of threads walking the mirror concurrently.

        @return: None
        """

        genomes_in_new_rel = self.read_genome_lists(new_list_genomes)

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
                            genomes_on_disk[accession] = genome_dir
                            fout.write('{}\t{}\t{}\n'.format(
                                accession, genome_dir, canonical_gid(accession)))

            progress.close()

        self.logger.info('Indexed {:,} genome directories.'.format(len(genomes_on_disk)))

        # what the release expects but the mirror does not have, and what the
        # mirror has but the release does not expect
        missing = sorted(genomes_in_new_rel - genomes_on_disk.keys())
        unexpected = sorted(genomes_on_disk.keys() - genomes_in_new_rel)

        # both are written even when empty, so a caller can count on them
        missing_file = output_file + '-missing'
        with open(missing_file, 'w') as fout:
            for accession in missing:
                fout.write('{}\n'.format(accession))

        unexpected_file = output_file + '-extra'
        with open(unexpected_file, 'w') as fout:
            for accession in unexpected:
                fout.write('{}\t{}\n'.format(accession, genomes_on_disk[accession]))

        if missing:
            self.logger.warning('{:,} genomes in the new release have no directory: {}'.format(
                len(missing), missing_file))
        if unexpected:
            self.logger.warning('{:,} directories hold a genome not in the new release: {}'.format(
                len(unexpected), unexpected_file))
        if not missing and not unexpected:
            self.logger.info('All {:,} genomes in the new release have a directory.'.format(
                len(genomes_in_new_rel)))

    def delete_empty_directory(self, genome_path: Union[str, Path]) -> bool:
        """
        Delete a specific path.

        Removing a genome leaves the digit-triplet directories that held it
        behind, so this walks back up the tree deleting each ancestor that the
        removal has emptied. It stops at the first ancestor that still has
        contents, which keeps it inside the mirror: every level above the
        triplets holds the other archives or the other genomes.

        @param genome_path: path to delete
        @return: True
        """
        if Path(genome_path).exists() and len(os.listdir(genome_path)) == 0:
            os.rmdir(genome_path)
            self.delete_empty_directory(os.path.dirname(genome_path))
        return True

    def clean_ftp(self,
                  new_list_genomes: List[str],
                  ftp_genome_dir_file: str,
                  report_dir: str,
                  taxonomy_file: Optional[str] = None) -> None:
        """Clean the FTP directory (remove deprecated genomes not appearing in the FTP folder anymore).

        Deletes from the mirror every genome the previous release held that the
        new release does not, and reports both that set and the genomes newly
        added.

        Parameters
        ----------
        new_list_genomes : list of str
            Files indicating the Gid present in the new release.
        ftp_genome_dir_file : str
            Genome directory file for the FTP server..
        report_dir : str
            Output directory to list reports.
        taxonomy_file : str, optional
            Standardised taxonomy file from NCBI.

        @return: None
        """

        make_sure_path_exists(report_dir)
        genome_in_new_rel = self.read_genome_lists(new_list_genomes)

        # read taxonomy file
        # only used to name the added genomes in the report; without it they
        # are still reported, as 'N/A'
        taxonomy = {}
        if taxonomy_file is not None:
            taxonomy = Taxonomy().read(taxonomy_file)

        # what the mirror holds now: accession -> genome directory
        current_ftp_genomes = {}
        with open(ftp_genome_dir_file) as fgdf:
            for line in fgdf:
                infos = line.strip().split('\t')
                current_ftp_genomes[infos[0]] = infos[1]

        # sorted() so the two reports have a stable line order and can be diffed between runs
        deleted_genomes = sorted(current_ftp_genomes.keys() - genome_in_new_rel)
        added_genomes = sorted(genome_in_new_rel - current_ftp_genomes.keys())

        self.logger.info('{:,} genomes have been deleted in the release'.format(len(deleted_genomes)))
        self.logger.info('{:,} genomes have been added in the release'.format(len(added_genomes)))

        # each genome is recorded before it is removed, so an interrupted run
        # leaves a report of what it had already deleted
        with open(os.path.join(report_dir, 'deleted_genomes.tsv'), 'w') as deleted_genome_file:
            for idx, deleted_genome in enumerate(deleted_genomes, 1):
                print("{:,}/{:,} genomes deleted".format(idx,
                                                     len(deleted_genomes)), end="\r")
                deleted_genome_file.write('{}\n'.format(deleted_genome))

                # a genome listed in the genome_dirs file may already be gone
                # from disk, which is not an error
                genome_dir = Path(current_ftp_genomes[deleted_genome])
                if genome_dir.is_dir():
                    shutil.rmtree(genome_dir)
                self.delete_empty_directory(genome_dir.parent)

        # report the species name of each added genome, index 6 being the
        # species rank of the seven Taxonomy() returns
        with open(os.path.join(report_dir, 'added_genomes.tsv'), 'w') as added_genome_file:
            for added_genome in added_genomes:
                added_genome_file.write('{}\t{}\n'.format(added_genome, taxonomy.get(added_genome, ['N/A'] * 7)[6]))
