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
update_genomes.py -- update the GTDB genome directories to match the NCBI mirror.

Each GTDB release is built by comparing the genomes held by the previous release
with the genomes NCBI currently offers on its FTP site. Genomes NCBI no longer
offers are removed, genomes new to NCBI are copied across, and genomes held by
both are checked for a changed genomic FASTA since the previous release. The
bookkeeping of which genomes fall into which of those three groups is done by
UpdateGenomes; the copying, comparing and reporting that follows, by FTPTools.

UpdateGenomes decides nothing beyond what the genome directory files of the
mirror and of the previous release already say. The mirror is built from the
selection, and the selection is where a genome is accepted or passed over, so
every genome the mirror holds is wanted.

A release can also be built with no previous release behind it at all. A FRESH
run takes every genome the mirror holds, copies it into the new release and
reports it as new: nothing is compared, no derived data is carried across, and
nothing is removed, there being nothing there to be removed. It is the run that
starts a release from the NCBI data alone -- because there is no previous release
to inherit from, or because everything derived from the genomes is to be
regenerated regardless of what did not change. It reads no previous genome
directory file, so the previous release is never named or opened, and the mirror
is the only thing it consults.

RefSeq and GenBank are updated in ONE pass over each genome directory file. They
were once two runs of this code, each filtered to an accession prefix, because
the GenBank run read the RefSeq run's output: a GenBank assembly was dropped when
its paired RefSeq assembly had been taken into the new release. That decision now
belongs to select_genomes.py, which makes it against the RefSeq genomes it selects
in the same pass, and nothing here consults one genome while handling another --
every genome is added, removed or compared on its own accession alone. Two runs
bought nothing after that but a second comparison pool idling until the first had
finished.

What the split did give was a RefSeq figure and a GenBank figure to compare, the
two churning very differently between releases, and that is kept: every count
this module logs is broken down by database. The reports are not, describing the
release as one thing.

The release tree does keep them apart. NCBI nests a genome under the archive its
accession names and the nine digits of that accession, three at a time, all of it
below one all/ directory:

    <mirror>/all/GCA/047/639/395/GCA_047639395.1_ASM4763939v1/

A GTDB release keeps that nesting but splits the two databases at the top, so a
genome's release path is its mirror path with all/ replaced by the database it
belongs to:

    <release>/genbank/GCA/047/639/395/GCA_047639395.1_ASM4763939v1/

which makes each database a tree that can be indexed, moved or handed on by
itself, as it could not be while the two were interleaved under all/.

The run writes the release's genome_dirs file itself, as genome_dirs.tsv in the
output directory, in the format list_genomes writes and every later command reads
(accession, absolute path, canonical accession). Walking the tree back with
list_genomes afterwards would only recover what this run already knew: it placed
every directory. Only genomes whose directory was actually written are in it --
see FTPTools.record_genome_dir() -- and a dry run, having written none, writes no
genome_dirs.tsv either.

An interrupted run is resumed with --resume, which reads that genome_dirs.tsv back
and does again only what it does not name. The file is what makes this possible
without a bookkeeping file of its own: a genome is written to it once its
directory is there, so a run stopped part way -- the disk filled, the job hit a
wall clock -- leaves behind an exact account of what survived it, and the genomes
it failed on or never reached are simply the ones missing from it. A directory it
was in the middle of writing is not named there either, so the genome is placed
again over whatever was got through.

What a resumed run writes is what one uninterrupted run would have written. The
reports are rewritten rather than appended to, because the rows the earlier run
wrote for the genomes it did NOT finish must not survive into the release's
account of itself -- a genome that failed for want of disk is to be described by
what happens to it this time, not by that. The rows for the genomes it did finish
are carried across before the rest of the release is handled, so report.log
describes every genome once through rather than the fragment the second run
happened to do.
"""

import os
import gzip
import shutil
import hashlib
import logging
import collections
import multiprocessing as mp
from concurrent.futures import FIRST_COMPLETED, Future, ThreadPoolExecutor, wait
from contextlib import contextmanager, ExitStack
from multiprocessing.queues import Queue
from typing import (Collection, Dict, Iterable, Iterator, List, Optional,
                    TextIO, Tuple)

from tqdm import tqdm

from gtdb_migration_tk import config
from gtdb_migration_tk.biolib_lite.common import canonical_gid
from gtdb_migration_tk.ncbi_utils import (GENOMIC_FASTA_EXT, MD5_MANIFEST,
                                          NCBI_DATABASES, NCBIDatabase,
                                          file_md5, read_md5_manifest)
from gtdb_migration_tk.utils.common import count_lines


# Outcomes a genome held by both the previous release and NCBI can have, as
# written to the report. Constants rather than literals so the tests, and any
# reader of the report, name the same strings this module writes. A dry run
# reports the same outcomes as a real one; it only refrains from copying.
STATUS_FASTA_UNCHANGED = 'genomic FASTA file unchanged'
STATUS_FASTA_CHANGED = 'genomic FASTA file changed'

# The genomic FASTA is a different FILE but the same GENOME: NCBI rewrote the
# deflines -- a new organism name, a corrected strain, a new assembly name --
# without touching a base, and the MD5 it publishes changed with them. Believing
# that MD5 alone throws away gene calls, rRNA classifications and tRNA scans that
# are still exactly right, and has every release recompute thousands of genomes
# for a change in somebody's free text. Between STATUS_FASTA_UNCHANGED and
# STATUS_FASTA_CHANGED because that is what it is: the file changed, the genome
# did not, and the derived data is carried across as though it had not.
STATUS_SEQUENCES_UNCHANGED = 'genomic FASTA sequences unchanged'

# A genome the comparison itself could not make, written as 'to_curate;<exception>'.
# The reason varies -- a manifest with no genomic FASTA entry, a directory that will
# not read -- so the counts fold every one of them into this single outcome.
STATUS_TO_CURATE = 'to_curate'

# How many distinct to_curate reasons the run names in the log. A systematic failure
# has one cause and a handful of variants; the report holds every row regardless.
CURATE_REASONS_LOGGED = 5

# How many copies are queued per thread while genomes are being added. Submitting
# every genome at once would build one future per genome before the first copy
# finished -- some 850 MB of them for a release of half a million genomes, which
# is exactly the release --fresh is for -- and a few per thread is all it takes
# for none of them to sit idle between one copy and the next.
COPY_QUEUE_DEPTH = 4

# The comparison bar carries four running counts as well as the usual rate and ETA,
# so it is given more room than the plain bars elsewhere in this module.
COMPARE_BAR_WIDTH = 165

# The database an accession beginning with neither prefix is counted under. It
# cannot come from a mirror built by ncbi_genome_sync, but a genome directory file
# is written by walking a tree, and while this command ran once per prefix such a
# row was dropped by both runs without a word. It is now updated like any other
# genome, and named here so that an unexpected one is seen rather than assumed.
UNKNOWN_DATABASE = 'unrecognised'

# What the listener process hands back once every comparison has been reported:
# the outcome counts the summary is built from, those same counts per database,
# and what the genomes that could not be compared at all blamed.
ComparisonTally = collections.namedtuple(
    'ComparisonTally', 'counts by_database reasons examples')


def tsv_safe(text: str) -> str:
    """Collapse whitespace so a value cannot break the TSV it is written into.

    An exception message carries whatever the filesystem put in it, newlines
    included, and the report is a table someone greps and cuts.

    Parameters
    ----------
    text : str
        Any value bound for a report column.

    @return: the value with every run of whitespace reduced to one space.
    """

    return ' '.join(str(text).split())


def curate_status(exc: Exception) -> str:
    """The outcome written for a genome that could not be compared at all.

    The exception TYPE alone says almost nothing -- a ValueError is a manifest
    with no genomic FASTA entry, an OSError is a directory that would not read --
    so the message goes in too. It is what turns a report full of identical
    to_curate rows into one that says what to fix.

    Parameters
    ----------
    exc : Exception
        What the comparison raised.

    @return: '<type>: <message>' under the to_curate outcome.
    """

    return '{};{}: {}'.format(STATUS_TO_CURATE, type(exc).__name__, tsv_safe(exc))


def report_outcome(row: str) -> str:
    """The outcome a report row records, as the counts group them.

    Parameters
    ----------
    row : str
        A report row, '<accession>\t<outcome>\n'.

    @return: the outcome, with a to_curate reason dropped so that every genome
             that could not be compared counts as one outcome however it failed.
    """

    outcome = row.rstrip('\n').split('\t')[-1]
    if outcome.startswith(STATUS_TO_CURATE + ';') or outcome == STATUS_TO_CURATE:
        return STATUS_TO_CURATE
    return outcome


def curate_reason(row: str) -> Tuple[str, str]:
    """What a to_curate row blames, split into what to group by and what to show.

    Grouping is by exception TYPE, not by message: a message names the genome's own
    path, so grouping by it gives one group per genome and a summary as long as the
    report. A run that fails systematically -- a previous release whose directories
    have moved, say -- fails the same way 16,000 times, and the type says that in
    one line while the message shows what it looked like.

    Parameters
    ----------
    row : str
        A to_curate report row.

    @return: (exception type, full reason) as written, either possibly ''.
    """

    outcome = row.rstrip('\n').split('\t')[-1]
    _, _, reason = outcome.partition(';')
    kind, _, _ = reason.partition(':')
    return kind.strip(), reason


def report_accession(row: str) -> str:
    """The genome a report row is about.

    Parameters
    ----------
    row : str
        A report row, '<accession>\t<outcome>\n'.

    @return: the accession, the first column of the row.
    """

    return row.split('\t')[0]


def database_of(gid: str) -> Optional[NCBIDatabase]:
    """Which NCBI database an accession belongs to.

    Read from the prefix rather than by slicing, so that the two prefixes are
    named in ncbi_utils.py alone and an accession matching neither is visible as
    such instead of falling in with whichever database was tested for first.

    Parameters
    ----------
    gid : str
        Genome accession, e.g. GCF_000000001.1.

    @return: the database, or None for an accession beginning with neither prefix.
    """

    for database in NCBI_DATABASES:
        if gid.startswith(database.prefix):
            return database

    return None


def database_label(gid: str) -> str:
    """The database an accession belongs to, as the counts name it.

    Parameters
    ----------
    gid : str
        Genome accession, e.g. GCF_000000001.1.

    @return: the database's label, or UNKNOWN_DATABASE for an accession beginning
             with neither prefix.
    """

    database = database_of(gid)

    return database.label if database is not None else UNKNOWN_DATABASE


def tally_by_database(gids: Iterable[str]) -> Dict[str, int]:
    """Count genomes by the database of their accession.

    BOTH databases are always named, one with no genomes included: a breakdown
    that drops an empty database reads as though it were never looked at, and the
    counts of a release are worth being seen to add up. UNKNOWN_DATABASE is named
    only when something fell under it, as it never should.

    Parameters
    ----------
    gids : iterable of str
        Genome accessions.

    @return: dict of database label to number of genomes, in the order
             ncbi_utils.NCBI_DATABASES gives, anything unrecognised last.
    """

    counted = collections.Counter(database_label(gid) for gid in gids)

    tally = {database.label: counted.pop(database.label, 0)
             for database in NCBI_DATABASES}
    tally.update(sorted(counted.items()))

    return tally


def count_by_database(gids: Iterable[str]) -> str:
    """The breakdown of a count of genomes, as the log writes it.

    Parameters
    ----------
    gids : iterable of str
        Genome accessions.

    @return: '12,345 RefSeq, 6,789 GenBank'.
    """

    return ', '.join('{:,} {}'.format(count, label)
                     for label, count in tally_by_database(gids).items())


def genomic_fasta(genome_dir: str) -> str:
    """The genomic FASTA NCBI serves for a genome.

    NCBI names the file for the assembly, <accession>_<asm_name>_genomic.fna.gz,
    and names the genome directory <accession>_<asm_name>, so the file wanted is
    known exactly and is named rather than searched for -- _cds_from_genomic.fna.gz
    and _rna_from_genomic.fna.gz end the same way and are different files.

    Parameters
    ----------
    genome_dir : str
        Genome directory, of the mirror or of a release.

    @return: path of the genomic FASTA in that directory, which may not exist.
    """

    assembly = os.path.basename(os.path.normpath(genome_dir))

    return os.path.join(genome_dir, assembly + GENOMIC_FASTA_EXT)


def sequences_md5(fasta_file: str) -> str:
    """Hash what a genomic FASTA says the GENOME is, rather than what the file is.

    NCBI reissues a genomic FASTA whose sequences have not changed at all: the
    deflines are rewritten as an organism is renamed or an assembly relabelled,
    and the MD5 NCBI publishes changes with them. That MD5 is the whole file, so
    it cannot tell that apart from a genome that was actually resequenced. This
    can, and it is only ever reached once the published MD5s have already
    disagreed -- a genome whose file is unchanged is never read at all.

    Three differences are deliberately invisible here, being the file's and not
    the genome's:

      * everything after the first whitespace of a defline, which is NCBI's free
        text -- the organism, the strain, the assembly, the plasmid note;
      * line wrapping, the sequence being hashed with its newlines removed;
      * base case, soft masking not changing which base is at a position.

    The contig ID itself, the first token, is NOT invisible, and neither is the
    division between one sequence and the next: the derived data carried across
    on the strength of this -- Prodigal's gene calls above all -- names the
    contigs it was called on and gives coordinates within them. A genome whose
    contig was renamed or whose contigs were merged or reordered has derived data
    that no longer lines up with it, however unchanged the bases, so it must
    hash differently and be regenerated.

    The file is read and hashed as it decompresses, nothing being written and no
    more than a line held: this runs on every genome whose MD5 disagreed, which
    on a large update is tens of thousands of them across --cpus workers.

    Parameters
    ----------
    fasta_file : str
        Gzipped genomic FASTA.

    @return: hex MD5 of the canonical form of its sequences.
    """

    digest = hashlib.md5()
    with gzip.open(fasta_file, 'rb') as handle:
        for line in handle:
            if line.startswith(b'>'):
                # the ID alone, on a line of its own, so that the boundary between
                # two sequences is part of what is hashed and concatenating them
                # could never produce the same digest
                fields = line[1:].split()
                digest.update(b'>' + (fields[0] if fields else b'') + b'\n')
            else:
                digest.update(line.strip().upper())

    return digest.hexdigest()


def release_genome_dir(new_genome_dir: str, mirror_relpath: str, gid: str) -> str:
    """Where a genome held by the mirror belongs in the new release.

    The mirror keeps NCBI's own layout, every genome of both databases nested by
    archive and accession digits below one all/ directory. A release keeps the
    nesting but splits the databases at the top, so all/ is replaced by the
    database the accession names and the rest of the path is carried over
    unchanged:

        all/GCA/047/639/395/GCA_047639395.1_ASM4763939v1
        -> genbank/GCA/047/639/395/GCA_047639395.1_ASM4763939v1

    The tail is taken from the archive directory rather than by dropping the
    first component, so the same path comes out whether --ftp_directory named the
    mirror root or the all/ directory inside it. The nesting itself is never
    rebuilt from the accession: it is NCBI's to define, the sync laid it down
    from the URLs NCBI published, and a release that reshaped it would no longer
    be the thing --verify checks.

    Parameters
    ----------
    new_genome_dir : str
        Output directory for the new release.
    mirror_relpath : str
        Genome directory of the mirror, relative to the mirror root.
    gid : str
        Accession of the genome, which names the database.

    @return: directory this genome is to be written to.

    Raises
    ------
    ValueError
        The accession belongs to neither database, or the path does not run
        through that database's archive directory. Either means the mirror is
        not laid out as the sync lays it out, and every genome placed from here
        on would be placed somewhere unintended.
    """

    database = database_of(gid)
    if database is None:
        raise ValueError('{} belongs to neither {}'.format(
            gid, ' nor '.join(d.label for d in NCBI_DATABASES)))

    parts = mirror_relpath.split(os.sep)
    if database.prefix not in parts:
        raise ValueError('{} is not held under a {} directory of the mirror: {}'.format(
            gid, database.prefix, mirror_relpath))

    archive = parts.index(database.prefix)

    return os.path.join(new_genome_dir, database.name, *parts[archive:])


def genome_dirs_row(gid: str, genome_dir: str) -> str:
    """One line of the genome_dirs file describing the new release.

    The format is list_genomes\'s, which every later command reads: see
    directory_manager.py. Columns may be appended but never reordered.

    Parameters
    ----------
    gid : str
        Accession of the genome, e.g. GCA_003004725.1.
    genome_dir : str
        Absolute path of the genome\'s directory in the new release.

    @return: \'<accession>\\t<path>\\t<canonical accession>\\n\'.
    """

    return '{}\t{}\t{}\n'.format(gid, genome_dir, canonical_gid(gid))


class UpdateGenomes:
    """Update the GTDB copy of the NCBI genomes from the mirror.

    Every genome held by the FTP mirror is of interest, the mirror being a copy
    of the genomes selected for the release, so genome selection is simply a
    matter of reading its genome directory file. Genomes are then added,
    removed, or compared relative to the previous GTDB release. Genomes are
    tracked as dictionaries mapping an accession to its genome directory.

    run_comparison() builds a release from the previous one; run_fresh() builds
    one from the mirror alone, adding every genome the mirror holds and
    comparing nothing. Both write the same reports through release_reports().

    One instance updates the whole release, RefSeq and GenBank together: what is
    done with a genome follows from its own accession and from nothing else, so
    the database it belongs to changes nothing here. It is still worth seeing
    apart, and every count logged carries the breakdown.
    """

    def __init__(self,
                 new_genome_dir: str,
                 dry_run: bool = False,
                 cpus: int = 1,
                 resume: bool = False) -> None:
        """Record where the release is written and how it is to be built.

        Parameters
        ----------
        new_genome_dir : str
            Output directory for the new release, where reports are written.
        dry_run : bool
            Report the changes that would be made without modifying any files.
        cpus : int
            Number of genomes compared, and copied, at once.
        resume : bool
            Continue an interrupted run, handling again only the genomes its
            genome_dirs file does not name.
        """

        # resolved here, so the genome_dirs file holds absolute paths however the
        # output directory was given on the command line, as list_genomes does
        self.new_genome_dir = os.path.abspath(new_genome_dir)
        self.dry_run = dry_run
        self.cpus = cpus
        self.resume = resume
        self.logger = logging.getLogger('timestamp')

    def report_file(self) -> str:
        """Report recording the fate of every genome of the release.

        @return: path of the report.
        """

        return os.path.join(self.new_genome_dir, 'report.log')

    def review_file(self) -> str:
        """Report recording the genomes needing manual attention.

        @return: path of the report.
        """

        return os.path.join(self.new_genome_dir, 'to_review.log')

    def genome_dirs_file(self) -> str:
        """Genome directory file describing the new release.

        @return: path of the genome_dirs file.
        """

        return os.path.join(self.new_genome_dir, 'genome_dirs.tsv')

    def completed_genomes(self) -> Dict[str, str]:
        """The genomes an interrupted run has already placed in the release.

        Read from the genome_dirs file of the OUTPUT directory, which names a
        genome only once its directory has been written -- see
        FTPTools.record_genome_dir() -- and is therefore the account of what an
        interrupted run left behind rather than of what it set out to do. A genome
        it failed on, or never reached, is absent from it and so is done again.

        The rows are taken as they stand rather than checked against the mirror.
        Verifying 1.2M placed genomes costs hours of stat calls on a mirror held
        over NFS, against seconds to read this file, and it would be answering a
        question the file has already answered: the row exists because the copy
        returned. Whether a tree holds what it should is ncbi_genome_sync --verify,
        which is where that check belongs and where it is already written.

        A row that never finished being written is the one thing a killed run can
        leave behind here, and it is dropped -- the genome it half-names is then
        placed again like any other the file does not name.

        @return: accession to genome directory for the genomes already placed;
                 empty when not resuming, or when there is no file to resume from.
        """

        if not self.resume:
            return {}

        genome_dirs_file = self.genome_dirs_file()
        if not os.path.exists(genome_dirs_file):
            self.logger.info(
                'Nothing to resume from at {}: the release is built from the start.'.format(
                    genome_dirs_file))
            return {}

        placed = {}
        with open(genome_dirs_file) as handle:
            for line in handle:
                fields = line.rstrip('\n').split('\t')
                if not line.endswith('\n') or len(fields) < 2:
                    self.logger.warning(
                        'Ignoring an unfinished row of {}: {!r}.'.format(
                            genome_dirs_file, line))
                    continue
                placed[fields[0]] = fields[1]

        self.logger.info('Resuming from {}: {:,} genomes are already in the release: {}.'.format(
            genome_dirs_file, len(placed), count_by_database(placed)))

        return placed

    def carried_rows(self, report_file: str, placed: Dict[str, str]) -> List[str]:
        """The rows of an interrupted run's report that still describe the release.

        What the earlier run said about a genome it FINISHED still stands, and is
        carried into the report the resumed run writes. What it said about one it
        did not finish does not: that genome is handled again, and the row it gets
        this time is the one the release is described by.

        Parameters
        ----------
        report_file : str
            Report of the interrupted run, read before it is opened for writing.
        placed : dict
            Accessions the interrupted run placed, from completed_genomes().

        @return: rows to write at the head of the new report, in their original
                 order; empty when not resuming or when there is no report to read.
        """

        if not placed or not os.path.exists(report_file):
            return []

        with open(report_file) as handle:
            return [line for line in handle
                    if line.endswith('\n') and line.split('\t')[0] in placed]

    def skip_completed(self,
                       genomes: Collection[str],
                       placed: Dict[str, str],
                       step: str) -> List[str]:
        """The genomes of a step that an interrupted run has not already placed.

        Parameters
        ----------
        genomes : collection
            Accessions the step would handle were the run not being resumed.
        placed : dict
            Accessions the interrupted run placed, from completed_genomes().
        step : str
            What the step does with them, for the line logged: 'add' or 'compare'.

        @return: accessions left to handle, in the order given.
        """

        remaining = [gid for gid in genomes if gid not in placed]
        if placed:
            self.logger.info(
                'Resuming: {:,} of the {:,} genomes to {} are already in the release, '
                'leaving {:,} to {}: {}.'.format(
                    len(genomes) - len(remaining), len(genomes), step,
                    len(remaining), step, count_by_database(remaining)))

        return remaining

    def load_genome_dirs(self, genome_dirs_file: str) -> Dict[str, str]:
        """Read the genomes of a release from a genome directory file.

        The file describes a whole release, RefSeq and GenBank together, and the
        whole of it is read: both databases are updated in the one pass.

        Parameters
        ----------
        genome_dirs_file : str
            Genome directory file (accession, path) of the mirror or of a release.

        @return: dict of accession to genome directory.
        """

        self.logger.info('Reading genomes from {}:'.format(genome_dirs_file))

        genome_paths = {}
        with open(genome_dirs_file, 'r') as f:
            for line in tqdm(f, total=count_lines(genome_dirs_file)):
                gid, path, *_ = line.split('\t')
                genome_paths[gid] = path.strip()

        self.logger.info(' - identified {:,} genomes: {}'.format(
            len(genome_paths), count_by_database(genome_paths)))

        return genome_paths

    def generate_genomes_to_remove(self,
                                   new_genomes: Dict[str, str],
                                   old_genomes: Dict[str, str]) -> Dict[str, str]:
        """Identify genomes present in the previous release, but no longer on the NCBI FTP site.

        Parameters
        ----------
        new_genomes : dict
            Accession to genome directory for genomes currently on the FTP site.
        old_genomes : dict
            Accession to genome directory for genomes in the previous release.

        @return: dict of accession to genome directory for genomes to remove.
        """

        removed_genomes = {gid: path for gid, path in old_genomes.items()
                           if gid not in new_genomes}
        self.logger.info('Identified {:,} genomes to remove: {}.'.format(
            len(removed_genomes), count_by_database(removed_genomes)))

        return removed_genomes

    def generate_genomes_to_add(self,
                                new_genomes: Dict[str, str],
                                old_genomes: Dict[str, str]) -> Dict[str, str]:
        """Identify genomes new to the NCBI FTP site since the previous release.

        Parameters
        ----------
        new_genomes : dict
            Accession to genome directory for genomes currently on the FTP site.
        old_genomes : dict
            Accession to genome directory for genomes in the previous release.

        @return: dict of accession to genome directory for genomes to add.
        """

        added_genomes = {gid: path for gid, path in new_genomes.items()
                         if gid not in old_genomes}
        self.logger.info('Identified {:,} genomes to add: {}.'.format(
            len(added_genomes), count_by_database(added_genomes)))

        return added_genomes

    def generate_genomes_to_compare(self,
                                    new_genomes: Dict[str, str],
                                    old_genomes: Dict[str, str]) -> List[str]:
        """Identify genomes common to the previous release and the NCBI FTP site.

        These genomes are candidates for an update as their files may have
        changed since the previous release.

        Parameters
        ----------
        new_genomes : dict
            Accession to genome directory for genomes currently on the FTP site.
        old_genomes : dict
            Accession to genome directory for genomes in the previous release.

        @return: list of accessions to compare between the two releases.
        """

        shared_genomes = list(old_genomes.keys() & new_genomes.keys())
        self.logger.info('Identified {:,} genomes to compare: {}.'.format(
            len(shared_genomes), count_by_database(shared_genomes)))

        return shared_genomes

    @contextmanager
    def release_reports(self, placed: Optional[Dict[str, str]] = None) -> Iterator['FTPTools']:
        """Open the files a run writes, and hand back what writes them.

        The reports are opened for the duration of the run and line buffered, so
        a run that fails part way through leaves behind the account of what it
        had done rather than an empty file. That is what a resumed run reads, and
        what it then writes: the rows an interrupted run left for the genomes it
        finished are carried into the new reports here, BEFORE the files are
        opened for writing, which truncates them.

        A dry run never opens the genome_dirs file, so a resume can be run dry to
        report what is left without putting the record it resumes from at risk.

        Shared by both ways of building a release -- from the previous release
        and afresh from the mirror -- because what a run writes does not depend
        on how it decided what to write: the same two reports, and the same
        genome_dirs file naming what was placed.

        Parameters
        ----------
        placed : dict
            Accessions an interrupted run already placed, from completed_genomes();
            None or empty for a run that is not being resumed.

        @return: context manager yielding the FTPTools the run is to use.
        """

        placed = placed or {}
        # read while the files still hold the interrupted run's account of itself
        carried = [self.carried_rows(report, placed)
                   for report in (self.report_file(), self.review_file(),
                                  self.genome_dirs_file())]

        with ExitStack() as reports:
            report = reports.enter_context(open(self.report_file(), 'w', 1))
            genomes_to_review = reports.enter_context(open(self.review_file(), 'w', 1))

            # a dry run writes no genome directory, so it has none to describe and
            # writes no genome_dirs file; the reports say what it would have done
            genome_dirs = None
            if not self.dry_run:
                genome_dirs = reports.enter_context(
                    open(self.genome_dirs_file(), 'w', 1))

            for rows, handle in zip(carried, (report, genomes_to_review, genome_dirs)):
                if handle is not None:
                    handle.writelines(rows)

            yield FTPTools(report, genomes_to_review, self.dry_run, genome_dirs,
                           self.resume)

            # inside the stack, so the handle is still open and line buffered: every
            # row written is on disk, and reading them back is how the count is of
            # what the file HOLDS rather than of what was meant to go into it
            if not self.dry_run:
                self.log_genome_dirs()

    def run_comparison(self,
                       ftp_dir: str,
                       ftp_genome_dirs: str,
                       old_genome_dirs: str) -> None:
        """Update the GTDB genome directories to match the mirror.

        The genomes of the mirror are compared to those of the previous GTDB
        release. Genomes no longer at NCBI are recorded as removed, new genomes
        are copied across, and genomes common to both are checked for a changed
        genomic FASTA and carried over with their derived data when it is not.

        Parameters
        ----------
        ftp_dir : str
            Root of the NCBI FTP mirror, replaced by the new release directory
            to place each genome.
        ftp_genome_dirs : str
            Genome directory file (accession, path) for the mirror.
        old_genome_dirs : str
            Genome directory file (accession, path) for the previous release.
        """

        self.logger.info('Updating the GTDB genome directories from the NCBI mirror.')

        placed = self.completed_genomes()
        with self.release_reports(placed) as ftptools:
            old_genomes = self.load_genome_dirs(old_genome_dirs)
            new_genomes = self.load_genome_dirs(ftp_genome_dirs)

            # a genome that has gone from NCBI is not in the release, so it is not
            # in what was placed and is reported again rather than carried across
            removed_genomes = self.generate_genomes_to_remove(new_genomes, old_genomes)
            ftptools.remove_genomes(removed_genomes)

            added_genomes = self.generate_genomes_to_add(new_genomes, old_genomes)
            added_genomes = {gid: added_genomes[gid] for gid
                             in self.skip_completed(added_genomes, placed, 'add')}
            ftptools.add_genomes(added_genomes, ftp_dir, self.new_genome_dir, self.cpus)

            shared_genomes = self.skip_completed(
                self.generate_genomes_to_compare(new_genomes, old_genomes),
                placed, 'compare')
            ftptools.compare_genomes(shared_genomes, old_genomes, new_genomes,
                                     ftp_dir, self.new_genome_dir, self.cpus)

    def run_fresh(self, ftp_dir: str, ftp_genome_dirs: str) -> None:
        """Build the GTDB genome directories afresh from the mirror alone.

        Every genome the mirror holds is copied into the new release and reported
        as new. No previous release is read, nothing is compared, and no derived
        data is carried across: the release starts from the NCBI data, and
        everything derived from it is to be produced again. Nothing is recorded
        as removed either, a run with no previous release behind it having
        nothing it could be dropping.

        Every genome of the release is copied here, so --cpus is what says how
        long the run takes: it is the number of genomes copied at once, and a
        copy is round trips to the file server rather than work for a CPU.

        What the run writes is what a run made against a previous release writes:
        every genome has the outcome 'new' in report.log, to_review.log is
        written and stays empty -- nothing having been looked for in a previous
        release to be found missing -- and a dry run reports the release it would
        build without copying a file.

        Parameters
        ----------
        ftp_dir : str
            Root of the NCBI FTP mirror, replaced by the new release directory
            to place each genome.
        ftp_genome_dirs : str
            Genome directory file (accession, path) for the mirror.
        """

        self.logger.info(
            'Building the GTDB genome directories afresh from the NCBI mirror, '
            'with no comparison to a previous release.')

        placed = self.completed_genomes()
        with self.release_reports(placed) as ftptools:
            new_genomes = self.load_genome_dirs(ftp_genome_dirs)

            # against an empty previous release: every genome the mirror holds is
            # new to a release that has nothing behind it, and reaching that
            # through the one method that decides it either way keeps the line
            # logged here the line the other run logs
            added_genomes = self.generate_genomes_to_add(new_genomes, {})
            added_genomes = {gid: added_genomes[gid] for gid
                             in self.skip_completed(added_genomes, placed, 'add')}
            ftptools.add_genomes(added_genomes, ftp_dir, self.new_genome_dir, self.cpus)

    def log_genome_dirs(self) -> None:
        """Report the genome directory file the run wrote.

        The release is what this file says it is from here on: it is what the
        call genes, marker and metadata steps are pointed at, and what the next
        release's update reads as its previous release. A count short of the
        genomes that went into the release is the first sign something was not
        placed, so it is stated rather than left to be noticed later.

        @return: None
        """

        with open(self.genome_dirs_file()) as handle:
            written = [line.split('\t')[0] for line in handle]

        self.logger.info('Wrote {:,} genomes to {}: {}.'.format(
            len(written), self.genome_dirs_file(), count_by_database(written)))


class FTPTools():
    """Carry out the genome directory changes required by a GTDB release.

    The decision about which genomes belong in a release is made by
    select_genomes.py, and the bookkeeping of what changed since the previous
    release by UpdateGenomes above; this class performs the resulting work. Genomes new to
    NCBI are copied into the new release, genomes that have gone are recorded,
    and for genomes held by both the previous release and NCBI the mirror's copy
    is taken and, if the genomic FASTA is unchanged, the derived data of the
    previous release (config.GTDB_DERIVED_DIRS_TO_COPY) is carried across with
    it, as that data took the most effort to produce and is still valid.

    Every genome handled is described in the report file. A genome held by both
    sides has one of the outcomes named by the STATUS_* constants above, or
    'to_curate;<exception>' if it could not be compared at all.

    Set dry_run to describe the changes in the reports without touching any file.
    Only the copying is withheld: every genome is still compared, so the reports
    and the counts logged at the end are the ones a real run would produce, which
    is what lets a release be sized before any of it is assembled.
    """

    def __init__(self,
                 report: TextIO,
                 genomes_to_review: TextIO,
                 dry_run: bool,
                 genome_dirs: Optional[TextIO] = None,
                 resume: bool = False) -> None:
        """Record the files to be written.

        Parameters
        ----------
        report : file
            Open file recording the fate of every genome in the release.
        genomes_to_review : file
            Open file recording genomes needing manual attention.
        dry_run : bool
            Report the changes that would be made without modifying any files.
        genome_dirs : file
            Open genome directory file describing the new release, or None to
            write none, as a dry run does.
        resume : bool
            Continue an interrupted run, whose part-written genome directories
            are to be replaced rather than refused.
        """

        self.report = report
        self.genomes_to_review = genomes_to_review
        self.dry_run = dry_run
        self.genome_dirs = genome_dirs
        self.resume = resume
        self.logger = logging.getLogger('timestamp')

    def record_genome_dir(self, gid: str, genome_dir: str) -> None:
        """Record where the new release holds a genome, for the genome_dirs file.

        Called for a genome whose directory has been WRITTEN, and for no other, so
        the file never names a path that is not there. That rules out three: a
        removed genome, which is not in the release at all; a genome that could
        not be compared, which was never copied, both MD5s being read before the
        copytree so that a failure leaves nothing behind; and every genome of a
        dry run, which is handed no file to write to.

        Parameters
        ----------
        gid : str
            Accession of the genome.
        genome_dir : str
            Directory the genome was written to.

        @return: None
        """

        if self.genome_dirs is not None:
            self.genome_dirs.write(genome_dirs_row(gid, genome_dir))

    def copy_genome(self, source: str, target: str) -> None:
        """Copy a genome's directory out of the mirror and into the release.

        copytree refuses a directory that is already there, which is the right
        answer for a run that is not being resumed: a release being written over
        the top of another is a mistake, and one caught before it has taken the
        wrong genome's files for its own. A RESUMED run meets that directory
        legitimately -- the genomes the interrupted run was copying when it
        stopped are not in its genome_dirs file, so they are placed again, and
        what is on disk is however far each copy had got. So the directory is
        removed first rather than the refusal being turned off, which would leave
        a half-copied genome with the files the new copy does not overwrite.

        Parameters
        ----------
        source : str
            Genome directory on the NCBI FTP mirror.
        target : str
            Genome directory to write for the new release.

        @return: None
        """

        if self.resume and os.path.exists(target):
            shutil.rmtree(target)
        shutil.copytree(source, target, symlinks=True)

    def add_genomes(self,
                    added_genomes: Dict[str, str],
                    ftp_dir: str,
                    new_directory: str,
                    cpus: int = 1) -> None:
        """Copy genomes new to NCBI into the new release.

        These genomes are on the FTP site but were not in the previous release,
        so there is nothing to compare against and the mirror's directory is
        taken whole. The mirror holds only the files GTDB keeps, the sync having
        fetched nothing else, so nothing is filtered here -- as nothing is in
        compare_genome_directories, which copies a shared genome the same way.

        The copies are made by cpus THREADS, where the comparison uses processes.
        Copying a genome directory is a few dozen round trips to an NFS server
        and no arithmetic at all, so the thread that issues one spends its time
        waiting in a system call with the GIL released, and threads are what
        overlap that waiting: 200 genome directories of the size NCBI serves --
        15 files and some 6 MB apiece -- take 17.1s copied one at a time over
        NFS, 5.2s copied eight at a time and 4.2s copied 32 at a time, the link
        itself being the ceiling from there. Threads also keep the reports where
        they belong -- this process holds the only handle on them and writes
        every row itself, so there is nothing for the listener process of
        compare_genomes to do here. An update adds a few thousand genomes and
        would hardly notice; a --fresh run adds the entire release, and copying
        one genome at a time is the whole of its work.

        Only COPY_QUEUE_DEPTH copies per thread are queued at a time, the next
        genome being submitted as one finishes. Submitting the lot up front is
        one line shorter and builds a future per genome before the first copy has
        finished, which for the release --fresh is meant for is hundreds of
        megabytes of them.

        The report row of a genome is written before its copy is attempted, as it
        is in a run made one genome at a time: the report is the account of every
        genome the run HANDLED. The genome_dirs row is written once the copy has
        returned, so it names only genomes that are there. A copy that fails
        stops the run, whatever is queued behind it dropped rather than left to
        finish: a release quietly short of a genome is worse than one that did
        not finish being built.

        Parameters
        ----------
        added_genomes : dict
            Accession to genome directory for the genomes to add.
        ftp_dir : str
            Base directory of the FTP mirror, replaced to form the target path.
        new_directory : str
            Base directory of the new release.
        cpus : int
            Number of genomes copied at once.
        """

        targets = {}
        for gid, path_record in added_genomes.items():
            targets[gid] = release_genome_dir(
                new_directory, os.path.relpath(path_record, ftp_dir), gid)
            self.report.write("{0}\tnew\n".format(gid))

        if self.dry_run:
            return

        threads = max(1, cpus)
        copying = {}
        pbar = tqdm(total=len(added_genomes), desc='Adding new genomes', ncols=100)
        pool = ThreadPoolExecutor(max_workers=threads)
        try:
            for gid, source in added_genomes.items():
                self.record_copied(copying, targets, pbar,
                                   threads * COPY_QUEUE_DEPTH)
                copying[pool.submit(self.copy_genome, source, targets[gid])] = gid

            # and the copies still running when the last genome was submitted
            self.record_copied(copying, targets, pbar, 1)
        except Exception:
            # shutting the pool down on its own would first run every copy already
            # queued behind the one that failed; shutdown(cancel_futures=True) would
            # say this in one line, but it arrived in 3.9 and the package supports 3.8
            for future in copying:
                future.cancel()
            raise
        finally:
            pool.shutdown()
            pbar.close()

    def record_copied(self,
                      copying: Dict[Future, str],
                      targets: Dict[str, str],
                      pbar: tqdm,
                      limit: int) -> None:
        """Wait until fewer than limit copies are outstanding, recording those done.

        Each copy that has finished is taken out of copying and written to the
        genome_dirs file, so that file names a genome once its directory is
        there. Called with the depth of the queue to make room for the next
        genome, and with 1 to wait for the last of them.

        Parameters
        ----------
        copying : dict
            Future of each copy still outstanding, to the accession it copies;
            every future that has finished is removed.
        targets : dict
            Accession to the directory the genome is copied to.
        pbar : tqdm
            Progress bar, advanced once per genome copied.
        limit : int
            Return once fewer than this many copies are outstanding.

        Raises
        ------
        Exception
            Whatever a copy raised, so that a genome that would not copy stops
            the run rather than being missing from the release.
        """

        while len(copying) >= limit:
            done, _ = wait(copying, return_when=FIRST_COMPLETED)
            for future in done:
                gid = copying.pop(future)
                future.result()
                self.record_genome_dir(gid, targets[gid])
                pbar.update()

    def remove_genomes(self, removed_genomes: Dict[str, str]) -> None:
        """Record the genomes NCBI no longer offers.

        These genomes are in the previous release but have gone from the FTP
        site. Nothing is deleted here: the new release is assembled in a fresh
        directory, so a genome is dropped by not being copied into it, and this
        report is the record of which genomes that applies to.

        Parameters
        ----------
        removed_genomes : dict
            Accession to genome directory for the genomes to drop.
        """

        for gid in removed_genomes:
            self.report.write("{0}\tremoved\n".format(gid))

    def compare_genomes(self,
                        shared_genomes: List[str],
                        old_genome_dirs: Dict[str, str],
                        new_genome_dirs: Dict[str, str],
                        ftp_directory: str,
                        new_directory: str,
                        threads: int) -> None:
        """Compare genomes held by both the previous release and the NCBI FTP site.

        Each genome is examined by a worker process which decides whether the
        derived data (e.g. Prodigal results) from the previous GTDB release should
        be retained. All NCBI data files are always copied from the NCBI FTP mirror 
        as these are the latest versions of these files.

        Parameters
        ----------
        shared_genomes : list
            Accessions held by both the previous release and the FTP site.
        old_genome_dirs : dict
            Accession to genome directory for the previous release.
        new_genome_dirs : dict
            Accession to genome directory for the FTP mirror.
        ftp_directory : str
            Base directory of the FTP mirror, replaced to form the target path.
        new_directory : str
            Base directory of the new release.
        threads : int
            Number of worker processes.
        """

        # populate worker queue with data to process
        worker_queue = mp.Queue()
        writer_queue = mp.Queue()
        tally_queue = mp.Queue()                 # the listener's counts, on its way back

        # the listener writes the genome_dirs row of each genome it reports as
        # compared, so it is given the same paths the workers are handed
        targets = {}

        for gca_record in shared_genomes:
            gtdb_dir = old_genome_dirs.get(gca_record)
            ftp_dir = new_genome_dirs.get(gca_record)
            target_dir = release_genome_dir(
                new_directory, os.path.relpath(ftp_dir, ftp_directory), gca_record)
            targets[gca_record] = target_dir

            worker_queue.put((gtdb_dir, ftp_dir, target_dir, gca_record))

        for _ in range(threads):
            worker_queue.put((None, None, None, None))

        # bound before the try, so the handler cannot fail with NameError when
        # creating the processes is itself what raised
        worker_proc, write_proc = [], None

        try:
            worker_proc = [mp.Process(target=self.__worker_thread, args=(
                worker_queue, writer_queue)) for _ in range(threads)]
            write_proc = mp.Process(target=self.__listener,
                                    args=(len(shared_genomes), writer_queue,
                                          tally_queue, targets))
            write_proc.start()

            for p in worker_proc:
                p.start()

            for p in worker_proc:
                p.join()

            writer_queue.put(None)
            # taken BEFORE the join: a process that has put an item on a queue does not
            # exit until it has been drained, so joining first can deadlock
            tally = tally_queue.get()
            write_proc.join()
        except Exception:
            for p in worker_proc:
                p.terminate()

            if write_proc is not None:
                write_proc.terminate()

            # genomes are left uncompared, which must not read as a completed run
            raise

        # a worker that died took its share of the genomes with it, and those
        # genomes are simply absent from the report rather than flagged
        failed = [p.exitcode for p in worker_proc if p.exitcode != 0]
        if failed:
            raise RuntimeError(
                '{} of {} comparison processes failed (exit codes: {}); '
                'the report is incomplete'.format(
                    len(failed), len(worker_proc), ', '.join(str(c) for c in failed)))

        # only once the report is known to be complete: a partial tally read as a whole
        # one would say a release was compared when most of it was not
        self.log_comparison(tally)

    def log_comparison(self, tally: ComparisonTally) -> None:
        """Report how the shared genomes came out, by outcome.

        The point of the comparison is how much derived data the release inherits,
        and that is a number the report file makes you count for yourself. A dry
        run compares in full and reports the same figures as the real run would,
        so this is what says whether a release is a small update or a large one
        before any file is copied.

        The two outcomes that carry derived data across are reported apart rather
        than added together. They are the same result reached by different means,
        and how often the second is reached is worth watching: it is the count of
        genomes NCBI reissued without changing, which is how much work a release
        would have repeated for nothing had it believed the published MD5.

        EVERY outcome is named, including the genomes that could not be compared,
        and the three are stated as adding to the total. An earlier version gave
        the total and then only the two FASTA outcomes, which read as a
        contradiction -- "Compared 16,843 shared genomes: 0 unchanged, 0 changed"
        -- when in truth every one of them had failed.

        The outcomes are then given per database, RefSeq and GenBank being
        updated in one pass since 0.1.7. That is the one thing the two runs this
        command used to make reported for free, and it is worth keeping: the two
        databases differ by an order of magnitude in size and churn differently
        between releases, so a single figure can hide a database that changed
        wholesale or not at all.

        Parameters
        ----------
        tally : ComparisonTally
            Outcome counts, the same counts per database, and the reasons of the
            genomes that could not be compared, as the listener tallied them.
        """

        unchanged = tally.counts.get(STATUS_FASTA_UNCHANGED, 0)
        sequences = tally.counts.get(STATUS_SEQUENCES_UNCHANGED, 0)
        changed = tally.counts.get(STATUS_FASTA_CHANGED, 0)
        curate = tally.counts.get(STATUS_TO_CURATE, 0)
        carried = 'would be carried across' if self.dry_run else 'carried across'
        self.logger.info(
            '{}Compared {:,} shared genomes: {:,} with an unchanged genomic FASTA '
            '(derived data {}), {:,} whose sequences were unchanged despite a '
            'differing MD5 (derived data {}), {:,} changed (derived data left '
            'behind, to be regenerated), {:,} that could not be compared.'.format(
                'DRY RUN: ' if self.dry_run else '',
                unchanged + sequences + changed + curate,
                unchanged, carried, sequences, carried, changed, curate))

        # both databases named, one with nothing shared included, so the lines are
        # always seen to add to the total above; anything unrecognised comes last
        # and only if a genome fell under it
        labels = [database.label for database in NCBI_DATABASES]
        labels += [label for label in sorted(tally.by_database) if label not in labels]
        for label in labels:
            outcomes = tally.by_database.get(label, {})
            self.logger.info(
                ' - {}: {:,} unchanged, {:,} sequences unchanged, {:,} changed, '
                '{:,} not compared.'.format(
                    label,
                    outcomes.get(STATUS_FASTA_UNCHANGED, 0),
                    outcomes.get(STATUS_SEQUENCES_UNCHANGED, 0),
                    outcomes.get(STATUS_FASTA_CHANGED, 0),
                    outcomes.get(STATUS_TO_CURATE, 0)))

        if not curate:
            return

        self.logger.warning(
            'warning: {:,} genome(s) could not be compared at all and carry no data '
            'into the release; they are the {} rows of {}.'.format(
                curate, STATUS_TO_CURATE, self.report.name
                if hasattr(self.report, 'name') else 'the report'))
        # the reasons, commonest first: a failure this systematic has one cause, and
        # naming it here saves reading a report with a row per genome to find it
        for kind, count in sorted(tally.reasons.items(),
                                  key=lambda kv: -kv[1])[:CURATE_REASONS_LOGGED]:
            self.logger.warning('  {:,} x {}'.format(
                count, tally.examples.get(kind) or kind or '(no reason given)'))
        if len(tally.reasons) > CURATE_REASONS_LOGGED:
            self.logger.warning('  ... and {:,} further kind(s) of failure; see {}'.format(
                len(tally.reasons) - CURATE_REASONS_LOGGED,
                self.report.name if hasattr(self.report, 'name') else 'the report'))

    def __worker_thread(self, queue_in: Queue, queue_out: Queue) -> None:
        """Compare one genome at a time until the queue is exhausted.

        Parameters
        ----------
        queue_in : multiprocessing.Queue
            Genomes to compare, as (gtdb_dir, ftp_dir, target_dir, accession);
            an accession of None ends the worker.
        queue_out : multiprocessing.Queue
            Report rows produced by the comparisons.
        """

        while True:
            gtdb_dir, ftp_dir, target_dir, gca_record = queue_in.get(
                block=True, timeout=None)
            if gca_record is None:
                break

            # a genome that cannot be compared is reported for curation; letting
            # it kill the worker silently drops every genome queued behind it
            try:
                status_gca = self.compare_genome_directories(gtdb_dir, ftp_dir, target_dir, gca_record)
            except Exception as e:
                status_gca = "{0}\t{1}\n".format(gca_record, curate_status(e))

            queue_out.put(status_gca)

    def __listener(self,
                   num_genomes: int,
                   writer_queue: Queue,
                   tally_queue: Queue,
                   targets: Dict[str, str]) -> None:
        """Write the outcome of every comparison to the report, and count the outcomes.

        The report is written by this process alone, so rows from the workers
        cannot interleave. Counting here rather than in the parent costs nothing:
        every row passes through this loop already, so the alternative is reading
        the whole report back once it is written.

        The genome_dirs row of a compared genome is written here for the same
        reason, and it is written HERE rather than by the worker that copied the
        genome because that is the one place a row cannot be torn in half by
        another process writing at the same moment.

        Parameters
        ----------
        num_genomes : int
            Number of genomes being compared, used to size the progress bar.
        writer_queue : multiprocessing.Queue
            Report rows from the workers, terminated by None.
        tally_queue : multiprocessing.Queue
            Receives the ComparisonTally once every row has been written; this
            process has its own copy of the counters, so they cannot be returned.
        targets : dict
            Accession to the directory the genome was written to, for the
            genome_dirs file.
        """

        counts = collections.Counter()
        # the same outcomes again, per database: the counts are what the two runs
        # this command used to make gave apart, and RefSeq and GenBank churn
        # differently enough between releases that one figure hides the other
        by_database = collections.defaultdict(collections.Counter)
        reasons = collections.Counter()
        examples = {}
        # wider than the other bars in this module to leave room for the running
        # counts; a comparison of a full release takes long enough that how it is
        # going matters, and a run failing wholesale shows in the first seconds
        pbar = tqdm(total=num_genomes, desc='Comparing shared genomes', ncols=COMPARE_BAR_WIDTH)
        for item in iter(writer_queue.get, None):
            self.report.write(item)
            gid = report_accession(item)
            outcome = report_outcome(item)
            counts[outcome] += 1
            by_database[database_label(gid)][outcome] += 1
            if outcome != STATUS_TO_CURATE:
                # anything but a failure means the mirror's directory was copied,
                # whether or not the derived data of the previous release went with it
                self.record_genome_dir(gid, targets[gid])
            if outcome == STATUS_TO_CURATE:
                # grouped by exception type, so the log can name the cause of a
                # systematic failure instead of leaving it to be found in a report of
                # 16,000 rows; one message per type is kept to show what it looked like
                kind, reason = curate_reason(item)
                reasons[kind] += 1
                examples.setdefault(kind, reason)
            # The same three numbers the summary ends with, as they accumulate, and in
            # the same order: passed as a dict because set_postfix SORTS its keyword
            # arguments, which would put them in alphabetical order instead.
            # refresh=False: the string is stored and drawn on the bar's own refresh,
            # not once per genome.
            pbar.set_postfix({
                'unchanged': '{:,}'.format(counts[STATUS_FASTA_UNCHANGED]),
                'seqs unchanged': '{:,}'.format(counts[STATUS_SEQUENCES_UNCHANGED]),
                'changed': '{:,}'.format(counts[STATUS_FASTA_CHANGED]),
                'failed': '{:,}'.format(counts[STATUS_TO_CURATE]),
            }, refresh=False)
            pbar.update()
        pbar.close()
        # plain dicts, a defaultdict carrying its factory through the queue being
        # more than the parent needs to be handed
        tally_queue.put(ComparisonTally(
            dict(counts),
            {label: dict(outcomes) for label, outcomes in by_database.items()},
            dict(reasons),
            examples))

    def compare_genome_directories(self,
                                   prev_gtdb_dir: str,
                                   ftp_dir: str,
                                   target_dir: str,
                                   genome_record: str) -> str:
        """Build the new release's copy of a genome held by both GTDB and NCBI.

        The mirror's directory is copied whole, md5checksums.txt included: what
        NCBI serves for a genome is what the new release carries, whether or not
        anything changed. The question is only whether the derived data of the
        previous release can come with it, and the genomic FASTA decides that.
        Its MD5 is read from the md5checksums.txt of each directory rather than
        computed, so the file is never opened and nothing has to be decompressed;
        NCBI published both sums, and the sync verified the mirror against its
        copy. If the two agree the sequence is unchanged, so gene calls, rRNA
        classifications and tRNA scans made on it still hold and the directories
        named by config.GTDB_DERIVED_DIRS_TO_COPY are copied across from the
        previous release.

        If they differ, the sequences are compared before anything is thrown
        away. NCBI reissues a genomic FASTA with rewritten deflines and identical
        sequences often enough to matter, and the published MD5, being of the
        whole file, changes with them; taking it at its word had every release
        recompute thousands of genomes that had not changed. sequences_md5() asks
        the narrower question -- the same contigs, under the same IDs, with the
        same bases? -- and when the answer is yes the derived data is carried
        across exactly as for an unchanged file, under its own outcome. Only when
        the sequences themselves differ is the derived data left behind, to be
        regenerated.

        The previous release's MD5 is not always there to be read. NCBI has served
        a md5checksums.txt naming files of its own that have nothing to do with
        the genome, and a release built while it did carries that manifest still.
        previous_genomic_fasta_md5() hashes the FASTA itself in that case, which
        is what the manifest was standing in for, so a genome is not recomputed
        over a file NCBI has since corrected. Only if the previous release cannot
        produce a readable FASTA either is the genome built from the mirror alone,
        as a new genome is. Both are noted in the review report: the previous
        release holds something that wants looking at, whichever way the genome
        went.

        A derived directory the previous release lacks is not an error, as a
        genome may not have had every step run on it; it is noted in the review
        report so that the step can be run this time.

        Under dry_run the comparison is made in full and reported, so the report
        of a dry run is the report the real run would write; only the copying is
        withheld, from the mirror and from the previous release alike.

        Parameters
        ----------
        prev_gtdb_dir : str
            Genome directory in the previous GTDB release.
        ftp_dir : str
            Genome directory on the NCBI FTP mirror.
        target_dir : str
            Genome directory to create for the new release.
        genome_record : str
            Accession of the genome.

        @return: report row of accession and outcome, one of the STATUS_* constants.
        """

        ftp_md5 = self.genomic_fasta_md5(ftp_dir)
        prev_md5 = self.previous_genomic_fasta_md5(prev_gtdb_dir, genome_record)

        if not self.dry_run:
            # a rerun must not inherit derived data from a run made before the
            # FASTA changed, so an existing target is replaced rather than added to
            if os.path.exists(target_dir):
                shutil.rmtree(target_dir)
            shutil.copytree(ftp_dir, target_dir, symlinks=True)

        if prev_md5 is None:
            # nothing the previous release holds can be shown to be this genome, so
            # there is nothing to carry across: it is built from the mirror alone,
            # exactly as a genome new to NCBI is, and reported as changed because
            # that is what the release does with it
            return '{}\t{}\n'.format(genome_record, STATUS_FASTA_CHANGED)

        if ftp_md5 != prev_md5:
            # the published MD5s disagree, which is not yet a reason to throw the
            # derived data away: NCBI rewrites deflines without touching a base, and
            # the MD5 of the whole file cannot tell that from a resequenced genome.
            # Only the sequences themselves can, so they are what is compared.
            if sequences_md5(genomic_fasta(ftp_dir)) != sequences_md5(genomic_fasta(prev_gtdb_dir)):
                return '{}\t{}\n'.format(genome_record, STATUS_FASTA_CHANGED)
            status = STATUS_SEQUENCES_UNCHANGED
        else:
            status = STATUS_FASTA_UNCHANGED

        for derived in config.GTDB_DERIVED_DIRS_TO_COPY:
            source = os.path.join(prev_gtdb_dir, derived)
            if not os.path.isdir(source):
                # noted on a dry run too: it is part of what the run would report
                self.genomes_to_review.write(
                    '{}\tno {} directory in the previous release: {}\n'.format(
                        genome_record, derived, prev_gtdb_dir))
                continue
            if not self.dry_run:
                # symlinks kept as symlinks: prodigal/ holds version-free links to the
                # marker results beside them, and following them would copy each
                # Pfam and TIGRFAM file twice
                shutil.copytree(source, os.path.join(target_dir, derived), symlinks=True)

        return '{}\t{}\n'.format(genome_record, status)

    def previous_genomic_fasta_md5(self,
                                   prev_gtdb_dir: str,
                                   genome_record: str) -> Optional[str]:
        """The MD5 of the previous release's genomic FASTA, hashed if the manifest cannot say.

        The manifest is asked first, as it is for the mirror, because reading it
        costs nothing. What it cannot answer is asked of the file: NCBI has
        published a md5checksums.txt listing its own scratch files and no genome
        file at all, and a release built while that was being served holds it to
        this day. The genome is not in doubt -- the FASTA is there and the mirror
        publishes an MD5 for it -- only the previous release's record of it, so
        the record is recomputed rather than the genome abandoned. Eight genomes
        of release 232 are in exactly that state.

        A previous release directory that is not there at all is a different
        thing: not one genome's bad manifest but an --old_genome_dirs_file naming
        a tree that has moved, which every genome of the run is about to hit. It
        is left to fail so that the run says so once, rather than quietly
        rebuilding the whole release as though nothing had come before.

        Parameters
        ----------
        prev_gtdb_dir : str
            Genome directory in the previous GTDB release.
        genome_record : str
            Accession of the genome, for the review report.

        @return: hex MD5 of the previous release's genomic FASTA, or None if the
            previous release cannot describe the genome at all.
        """

        try:
            return self.genomic_fasta_md5(prev_gtdb_dir)
        except (ValueError, OSError) as manifest_error:
            if not os.path.isdir(prev_gtdb_dir):
                raise

            try:
                checksum = file_md5(genomic_fasta(prev_gtdb_dir))
            except OSError as fasta_error:
                self.genomes_to_review.write(
                    '{}\ttreated as new: the previous release can neither name nor hold '
                    'its genomic FASTA ({}; {})\n'.format(
                        genome_record, tsv_safe(manifest_error), tsv_safe(fasta_error)))
                return None

            self.genomes_to_review.write(
                '{}\tgenomic FASTA of the previous release hashed, its manifest being '
                'unusable ({})\n'.format(genome_record, tsv_safe(manifest_error)))
            return checksum

    def genomic_fasta_md5(self, genome_dir: str) -> str:
        """Read the MD5 NCBI publishes for the genomic FASTA of a genome.

        The entry wanted is the file genomic_fasta() names, looked up by name
        rather than searched for. The manifest is read by
        ncbi_utils.read_md5_manifest(), as the sync reads it to mirror and verify.
        A manifest that cannot name the FASTA is an error here; what the previous
        release does about it is previous_genomic_fasta_md5()'s to decide, and the
        mirror's manifest not naming it is left to fail, being a genome the sync
        has something to answer for.

        Parameters
        ----------
        genome_dir : str
            Genome directory holding an md5checksums.txt.

        @return: hex MD5 of the genomic FASTA, as recorded in the manifest.
        """

        wanted = os.path.basename(genomic_fasta(genome_dir))
        manifest = os.path.join(genome_dir, MD5_MANIFEST)

        with open(manifest) as handle:
            for checksum, name in read_md5_manifest(handle):
                if name == wanted:
                    return checksum

        raise ValueError('{} has no entry for {}'.format(manifest, wanted))
