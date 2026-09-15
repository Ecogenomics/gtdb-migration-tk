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
select_genomes.py -- select the NCBI genomes belonging in a GTDB release.

SelectedGenomesManager answers which genomes a release wants, and answers it from
the summary files alone: no previous release, no FTP mirror, no genome
directories. It is what the select_genomes command runs, so the genomes wanted in
a release can be settled and reviewed before any of them is copied. The mirror is
then built from this selection by ncbi_genome_sync.py, and the release directories
updated from the mirror by update_genomes.py.

Which genomes NCBI considers latest, which belong to large multi-isolate projects,
and which are paired with a RefSeq assembly are all read from the NCBI assembly
summary files by ncbi_utils.py, which ncbi_genome_sync.py reads them with too.
"""

import os
import sys
import gzip
import logging
from typing import Dict, List, Set, Tuple

from tqdm import tqdm

from gtdb_migration_tk.biolib_lite.common import canonical_gid
from gtdb_migration_tk.ncbi_utils import (
    GENBANK_PREFIX, REFSEQ_PREFIX, count_summary_rows, has_ftp_path,
    read_assembly_summary)


# Suffixes of the assembly summary file names of the two NCBI databases. NCBI
# names each file for the database it describes, and the per-domain copies GTDB
# works from keep that suffix: assembly_summary_archaea_refseq.txt,
# assembly_summary_bacteria_genbank.txt. Reading the suffix is what lets each
# file be parsed once, RefSeq files before GenBank files, rather than every file
# being parsed once per database.
REFSEQ_SUFFIX = '_refseq.txt'
GENBANK_SUFFIX = '_genbank.txt'

# Values of excluded_from_refseq marking an assembly as one of the thousands of
# near-identical isolates NCBI sequences for outbreak and pathogen surveillance.
#
# NCBI renamed this annotation: summary files written before the change carry
# 'derived from surveillance project', and those written after it carry 'large
# multi-isolate project', for the same genomes. Both are matched, so a selection
# made from an archived set of summary files excludes what a selection made from
# a current set excludes. The column holds several such annotations separated by
# semicolons, hence the substring test.
MULTI_ISOLATE_TAGS = ('large multi-isolate project', 'surveillance')

# Placeholder for the notes column of the selected genome table, matching the
# convention NCBI uses for an empty field in the summary files themselves.
NO_NOTE = 'na'

# Table of the genomes selected for a new GTDB release, written by the
# select_genomes command into its output directory. It is gzipped, as the files
# of a release are once they are part of GTDB, and because the table runs to one
# row per genome in NCBI: a few million lines of highly repetitive accessions and
# FTP paths, which compress to roughly a tenth of their size.
SELECTED_GENOMES_FILE = 'gtdb_selected_genomes.tsv.gz'

# One row of that table: accession, ftp_path, version_status,
# excluded_from_refseq, gbrs_paired_asm, notes.
SelectedRow = Tuple[str, str, str, str, str, str]

# Header of the table. It is '#'-prefixed so the file is read by the same
# readers as an NCBI assembly summary file, which skip comment lines, while
# still naming its columns for anyone opening it.
#
# The first four columns are exactly the four ncbi_genome_sync reads, in the
# order it writes its own .fail and .bad files, so this table can be handed
# straight to it as the list of genomes to mirror. It needs assembly_accession
# and ftp_path, and renders assembly_status.txt from version_status and
# excluded_from_refseq; the last two columns it simply ignores, as every reader
# of these tables locates columns by name.
#
# The notes column carries what would otherwise be a per genome warning in the
# log. It exists because NCBI's summary files routinely pair a GenBank assembly
# with a RefSeq assembly they do not themselves list: too common to report one
# line at a time, and too material to drop, since it is the reason a genome that
# appears to be covered by RefSeq was taken from GenBank instead.
SELECTED_GENOMES_HEADER = ('#assembly_accession\tftp_path\tversion_status'
                           '\texcluded_from_refseq\tgbrs_paired_asm\tnotes')


def is_multi_isolate(excluded_from_refseq: str) -> bool:
    """Report whether an assembly belongs to a large multi-isolate project.

    Parameters
    ----------
    excluded_from_refseq : str
        Value of the excluded_from_refseq column of an assembly summary file.

    @return: True if the assembly is annotated as part of such a project.
    """

    return any(tag in excluded_from_refseq for tag in MULTI_ISOLATE_TAGS)


def write_selected_genomes(selected: List[SelectedRow], output_file: str) -> None:
    """Write the table of genomes selected for a new GTDB release.

    Rows are sorted by accession rather than left in the order the summary files
    were read, so that the tables of two releases can be compared directly to
    see what the new release gained and lost.

    Parameters
    ----------
    selected : list
        Selected genomes as (accession, ftp_path, version_status,
        excluded_from_refseq, gbrs_paired_asm, notes) rows.
    output_file : str
        Gzipped table to write.
    """

    with gzip.open(output_file, 'wt') as table:
        table.write(SELECTED_GENOMES_HEADER + '\n')
        for row in sorted(selected):
            table.write('\t'.join(row) + '\n')


class SelectedGenomesManager:
    """Select the NCBI genomes belonging in a GTDB release, from the summary files alone.

    This is the one decision in a release about which genomes are wanted, and
    it is made without reference to a previous release or to the FTP mirror: the
    NCBI assembly summary files are the only input, and the answer is a table of
    the genomes wanted. Nothing is copied, compared, or deleted, so the selection
    can be produced and reviewed before any part of a release is built; the sync
    then mirrors exactly this table, and GenomeManager updates the release from
    that mirror without deciding anything further.

    Every RefSeq assembly is wanted. A GenBank assembly is wanted only where
    RefSeq does not already cover the genome and the assembly is not part of a
    large multi-isolate project, those projects being the thousands of
    near-identical isolates NCBI excludes from RefSeq and which GTDB does not
    want either.

    A genome NCBI lists but does not serve -- ftp_path 'na' -- is not selected
    from either database: there is nothing to mirror, and the table this writes
    is what ncbi_genome_sync fetches from. That cuts both ways, and the second
    way is the point: RefSeq covers a genome only if the RefSeq assembly was
    itself selected, so a GenBank assembly whose RefSeq counterpart has no
    ftp_path is NOT covered and is selected in its place. The alternative would
    drop the genome from the release entirely -- NCBI holds it, under its
    GenBank accession, and nothing else would carry it.

    Two assumptions the rule rests on are verified rather than trusted, since a
    summary file that breaks either yields a plausible looking but wrong
    selection: that NCBI flags no RefSeq assembly as a large multi-isolate
    project, and that every assembly listed is the latest version of itself. A
    genome breaking either is reported and left out of the selection; neither is
    treated as fatal, so a single odd row cannot cost a whole run.

    Which database a summary file describes is read from its name, so each file
    is parsed once: the RefSeq files first, since a GenBank assembly cannot be
    judged until the genomes RefSeq covers are known, then the GenBank files.

    The table it writes carries the four columns ncbi_genome_sync reads, so the
    selection is directly the list of genomes to mirror: nothing has to be joined
    back against the NCBI summary files to sync what was chosen.

A third discrepancy is recorded but not acted on. NCBI regularly pairs a
    GenBank assembly with a RefSeq assembly that its own summary files do not
    list, which is an error in those files rather than anything about the genome.
    The genome is selected, since nothing else covers it, and the pairing is
    written to the notes column of its row: reporting each one in the log would
    bury the failures that do need attention. A GenBank assembly taken in place
    of a RefSeq assembly that WAS listed but not selected carries a note too, and
    a different one, naming the reason the RefSeq assembly was passed over --
    those two look identical in the table otherwise, and only one of them is
    NCBI's mistake.
    """

    def __init__(self, output_dir: str) -> None:
        """Record where the selection and its log are to be written.

        Parameters
        ----------
        output_dir : str
            Output directory for the selected genome table.
        """

        self.output_dir = output_dir
        self.logger = logging.getLogger('timestamp')

    def _read_summary(self, assembly_summary: str):
        """Read the fields of an assembly summary file needed to select genomes.

        Parameters
        ----------
        assembly_summary : str
            NCBI assembly summary file.

        @return: iterator over (accession, version_status, gbrs_paired_asm,
            excluded_from_refseq, ftp_path), one tuple per genome.
        """

        records = read_assembly_summary(assembly_summary,
                                        'assembly_accession',
                                        'version_status',
                                        'gbrs_paired_asm',
                                        'excluded_from_refseq',
                                        'ftp_path')

        return tqdm(records,
                    total=count_summary_rows(assembly_summary),
                    desc='Selecting genomes from {}'.format(
                        os.path.basename(assembly_summary)))

    def _is_latest(self, accession: str, version_status: str) -> bool:
        """Report whether NCBI considers an assembly the latest version of itself.

        A superseded assembly in the input is a sign the summary files describe
        something other than the genomes of the new release, so it is reported
        rather than quietly passed over.

        Parameters
        ----------
        accession : str
            Assembly accession of the genome.
        version_status : str
            Value of the version_status column.

        @return: True if the assembly is flagged as the latest version.
        """

        if version_status == 'latest':
            return True

        self.logger.warning(
            '{} has version_status "{}" rather than "latest" and was not selected.'.format(
                accession, version_status or 'na'))

        return False

    def group_by_database(self, new_list_genomes: List[str]) -> Tuple[List[str], List[str]]:
        """Split the assembly summary files into the RefSeq files and the GenBank files.

        A file whose name says neither is fatal rather than ignored: skipping it
        would leave every genome it lists out of the release, which looks exactly
        like a successful run of a smaller release.

        Parameters
        ----------
        new_list_genomes : list of str
            NCBI assembly summary files describing the new release.

        @return: tuple of the RefSeq files and the GenBank files.
        """

        refseq_files, genbank_files = [], []

        for assembly_summary in new_list_genomes:
            # the suffix is what names the database, and sits behind the .gz on a
            # compressed file; ncbi_metadata_sync writes these gzipped, older
            # releases hold them plain, and both are read the same way
            name = os.path.basename(assembly_summary)
            if name.endswith('.gz'):
                name = name[:-len('.gz')]

            if name.endswith(REFSEQ_SUFFIX):
                refseq_files.append(assembly_summary)
            elif name.endswith(GENBANK_SUFFIX):
                genbank_files.append(assembly_summary)
            else:
                self.logger.error(
                    'Cannot tell which NCBI database {} describes: the name of an assembly '
                    'summary file must end with {} or {}, optionally gzipped.'.format(
                        os.path.basename(assembly_summary), REFSEQ_SUFFIX, GENBANK_SUFFIX))
                sys.exit()

        # either list being empty is legal, but it is far more often a mistake in
        # the file list than a release genuinely drawn from one database
        for database, files in (('RefSeq', refseq_files), ('GenBank', genbank_files)):
            if not files:
                self.logger.warning(
                    'No {} assembly summary files were given; the selection will '
                    'contain no {} genomes.'.format(database, database))

        return refseq_files, genbank_files

    def select_refseq(self,
                      refseq_files: List[str]
                      ) -> Tuple[List[SelectedRow], Set[str], Dict[str, str]]:
        """Select the RefSeq genomes of the new release.

        Every RefSeq assembly is wanted, so this is a matter of reading the
        summary files. NCBI is not expected to flag any RefSeq assembly as a
        large multi-isolate project; one that is flagged is reported and
        dropped, as it would otherwise enter the release by the very route the
        GenBank rule exists to close.

        This runs before select_genbank, which needs the genomes RefSeq covers.

        Parameters
        ----------
        refseq_files : list of str
            NCBI RefSeq assembly summary files describing the new release.

        @return: tuple of the selected rows (accession, ftp_path, version_status,
            excluded_from_refseq, gbrs_paired_asm, notes), the canonical accessions
            they cover, and the canonical accession of every RefSeq genome that was
            read but NOT selected, mapped to why. select_genbank needs both: the
            first says which genomes need no GenBank copy, the second explains, in
            the notes column, why a genome that appears to have a RefSeq assembly
            was taken from GenBank anyway.
        """

        selected = []
        covered = set()
        dropped = {}
        read = 0
        misfiled = 0
        no_ftp_path = 0

        for assembly_summary in refseq_files:
            for accession, version_status, paired_asm, excluded, ftp_path in self._read_summary(
                    assembly_summary):

                if not accession.startswith(REFSEQ_PREFIX):
                    # the file is named for RefSeq but this row is not a RefSeq
                    # assembly, so the file list says something the table does not
                    misfiled += 1
                    continue

                read += 1

                gid = canonical_gid(accession)

                if not self._is_latest(accession, version_status):
                    dropped[gid] = 'version_status is "{}"'.format(version_status or 'na')
                    continue

                if is_multi_isolate(excluded):
                    self.logger.warning(
                        '[Unexpected] RefSeq genome {} is annotated as "{}" and was not '
                        'selected.'.format(accession, excluded))
                    dropped[gid] = 'annotated as "{}"'.format(excluded)
                    continue

                if not has_ftp_path(ftp_path):
                    # NCBI lists it but serves no directory, so there is nothing to
                    # mirror. Recorded in `dropped` so that its GenBank counterpart,
                    # which may well have a directory, is selected in its place.
                    no_ftp_path += 1
                    dropped[gid] = 'NCBI serves no directory for it'
                    continue

                selected.append((accession, ftp_path, version_status, excluded,
                                 paired_asm, NO_NOTE))
                covered.add(gid)

        self.logger.info('RefSeq: selected {:,} of {:,} genomes.'.format(len(selected), read))
        if no_ftp_path:
            self.logger.warning('RefSeq: {:,} genomes ignored, NCBI serves no directory for '
                                'them (ftp_path is na); where a GenBank assembly exists it is '
                                'selected instead.'.format(no_ftp_path))
        if misfiled:
            self.logger.warning('RefSeq: {:,} rows in the RefSeq assembly summary files are '
                                'not RefSeq assemblies and were ignored.'.format(misfiled))

        return selected, covered, dropped

    def select_genbank(self,
                       genbank_files: List[str],
                       covered: Set[str],
                       dropped: Dict[str, str]) -> List[SelectedRow]:
        """Select the GenBank genomes needed to supplement the RefSeq genomes.

        A GenBank assembly is wanted when RefSeq does not already cover the
        genome and the assembly is not part of a large multi-isolate project.
        Coverage is decided by canonical accession against the RefSeq genomes
        actually selected above, not by the gbrs_paired_asm field, so a GenBank
        assembly paired with a RefSeq assembly that is absent from the new
        release is still selected -- including when that assembly was read and
        passed over because NCBI serves no directory for it, which is the case
        this rule exists to catch. Where the two disagree NCBI describes a RefSeq
        assembly its own summary files do not list; that is common enough to be
        recorded in the notes column of the affected row rather than reported
        genome by genome.

        Parameters
        ----------
        genbank_files : list of str
            NCBI GenBank assembly summary files describing the new release.
        covered : set
            Canonical accessions of the selected RefSeq genomes.
        dropped : dict
            Canonical accession to the reason its RefSeq assembly was not selected.

        @return: list of selected rows (accession, ftp_path, version_status,
            excluded_from_refseq, gbrs_paired_asm, notes).
        """

        selected = []
        read = 0
        multi_isolate = 0
        unpaired = 0
        misfiled = 0
        no_ftp_path = 0
        rescued = 0

        for assembly_summary in genbank_files:
            for accession, version_status, paired_asm, excluded, ftp_path in self._read_summary(
                    assembly_summary):

                if not accession.startswith(GENBANK_PREFIX):
                    misfiled += 1
                    continue

                read += 1

                gid = canonical_gid(accession)

                if not self._is_latest(accession, version_status):
                    continue

                if gid in covered:
                    continue

                if not has_ftp_path(ftp_path):
                    # nothing to mirror. Unlike the RefSeq case this rescues
                    # nothing: RefSeq was already asked first and did not cover it
                    no_ftp_path += 1
                    continue

                if is_multi_isolate(excluded):
                    multi_isolate += 1
                    continue

                note = NO_NOTE
                if paired_asm.startswith(REFSEQ_PREFIX):
                    # The row appears to have a RefSeq counterpart, yet here it is
                    # being taken from GenBank; the note says why, and the two
                    # reasons are not the same kind of thing.
                    reason = dropped.get(gid)
                    if reason is None:
                        # NCBI names a RefSeq assembly that its own summary files
                        # do not list: an error in those files, nothing about the
                        # genome
                        unpaired += 1
                        note = ('paired RefSeq assembly {} absent from the assembly '
                                'summary files'.format(paired_asm))
                    else:
                        # the RefSeq assembly was listed and passed over, so the
                        # GenBank copy is the only one this release can carry
                        rescued += 1
                        note = ('paired RefSeq assembly {} was not selected: {}'
                                .format(paired_asm, reason))

                selected.append((accession, ftp_path, version_status, excluded,
                                 paired_asm, note))

        self.logger.info('GenBank: selected {:,} of {:,} genomes.'.format(len(selected), read))
        self.logger.info('GenBank: {:,} genomes skipped as part of a large multi-isolate '
                         'project.'.format(multi_isolate))
        if no_ftp_path:
            self.logger.warning('GenBank: {:,} genomes ignored, NCBI serves no directory for '
                                'them (ftp_path is na).'.format(no_ftp_path))
        if rescued:
            self.logger.info('GenBank: {:,} genomes selected in place of a RefSeq assembly '
                             'that was read but not selected; the reason is in the notes '
                             'column of each.'.format(rescued))
        if misfiled:
            self.logger.warning('GenBank: {:,} rows in the GenBank assembly summary files are '
                                'not GenBank assemblies and were ignored.'.format(misfiled))
        if unpaired:
            # common enough that it is reported as a count rather than per genome:
            # the affected rows carry a note, and the fault is NCBI's
            self.logger.info('GenBank: {:,} selected genomes name a paired RefSeq assembly '
                             'absent from the assembly summary files; each is noted in the '
                             'table.'.format(unpaired))

        return selected

    def run(self, new_list_genomes: List[str]) -> None:
        """Write the table of genomes selected for the new release.

        Parameters
        ----------
        new_list_genomes : list of str
            NCBI assembly summary files describing the new release.
        """

        refseq_files, genbank_files = self.group_by_database(new_list_genomes)
        self.logger.info('Reading {:,} RefSeq and {:,} GenBank assembly summary file(s).'.format(
            len(refseq_files), len(genbank_files)))

        refseq, covered, dropped = self.select_refseq(refseq_files)
        genbank = self.select_genbank(genbank_files, covered, dropped)

        output_file = os.path.join(self.output_dir, SELECTED_GENOMES_FILE)
        write_selected_genomes(refseq + genbank, output_file)

        self.logger.info('Selected {:,} genomes for the new release: {}'.format(
            len(refseq) + len(genbank), output_file))
