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

"""Bring the genomes table of the GTDB database into line with a new release.

WHAT DECIDES A GENOME

update_genomes says what became of every genome of the release in report.log,
one outcome per genome, and writes genome_dirs.tsv for the genomes the release
holds. Those two files decide everything here: each genome is decided once, from
its outcome and whether the database holds it, and nothing is inferred.

Until 0.1.41 the genomes to add were read from a CheckM profile, CheckM having
been run on exactly the new and changed genomes; the command was run once per
NCBI database (--repository); and the work was three passes, each committed on
its own, the failure of one caught and passed over. The profile stopped being
that list when CheckM began to leave genomes out (no proteins, too large); it
never named a genome whose deflines alone had changed, so that genome kept a
hash of a file it no longer had; and a genome missing from the database that it
did not name was never added. A release now covers both databases, so each
genome's source is read from its accession's prefix and looked up by name in
genome_sources.

    outcome                      in the database           not in the database
    new                          (cannot be)               added, or takes over its
                                                           predecessor's row
    genomic FASTA file changed   sequences changed         added
    sequences unchanged          updated                   added
    genomic FASTA file unchanged updated, or unchanged     added
    removed, to_curate           deleted                   nothing to do

A genome the database does not hold whose directory has no genomic FASTA is
not added (no genomic FASTA): NCBI publishes some assemblies as their reports
alone (r237: GCA_056491145.1 and GCF_056491165.1), the database cannot record
a genome without a file, and a release from 0.1.41 on reports such a genome for
curation rather than holding it. A new version without one does not take over
its predecessor's row. A genome the database DOES hold whose genomic FASTA is
gone stops the run when its files are hashed: that is a release that lost a
file, not a genome NCBI never published.

A genome the database holds and the report does not name is left as it is and
counted in the log: the report names every genome of the previous release, so
such a row is one this command did not put there.

A NEW VERSION

A new version of an accession (GCA_001317685.3, its predecessor .2 removed)
takes over its predecessor's row rather than being added beside it, so the
genome keeps its id and with it its place in every curated genome list. Its
aligned markers are deleted, being made from sequences it no longer has, as are
those of a genome whose sequences changed under the same accession. A
predecessor that is still in the release is left alone, and the new version is
added as a genome of its own.

THE HASHES

fasta_file_sha256 and genes_file_sha256 hold, despite their names, the SHA-1 of
the DECOMPRESSED file, as prodigal and hmmsearch record theirs (biolib_lite's
sha256_rb() is SHA-1, to match the original GTDB code). The content is hashed
rather than the gzip file because the gzip header records a time and a file name
and its bytes depend on the compressor: the same proteins gzipped twice hash
differently. A genome's files are hashed where they are new or may have changed:
both files of a genome added or versioned or whose sequences changed, the
genomic FASTA of a genome whose deflines NCBI rewrote, and a protein file the
database holds no hash for. --rehash_all hashes every genome of the release,
putting right the rows hashed with SHA-256 before 2022 and the genomic hashes
that no longer match their file.

ONE TRANSACTION

The files are hashed before the database is written to, on --cpus processes,
and the update is then made in one transaction: the release is in the database
whole or not at all, and the hours of hashing are not spent holding locks.
--dry_run makes every change and rolls it back, so the statements are tried
against the real tables and the report says what a run would have done.

has_changed is TRUE for the genomes whose sequences this release brought and
FALSE for every other NCBI genome, so that it says what changed in this release
and nothing older: a genome the report calls new or changed, and one whose
date_added is this update's download date, which is a genome an earlier run of
the same update added.

Deleting a genome deletes its metadata, its aligned markers and its place in
every curated genome list, the foreign keys cascading. The lists each deleted
genome was in are appended to genome_lists_affected.tsv, with when the run began
and whether it was a dry run, and are on disk before the transaction commits,
so that no crash loses them.

STOPPED PART WAY

Nothing reaches the database until the one COMMIT, so a run stopped at any
point -- an exception, Ctrl-C, a machine reset -- leaves it as it was or wholly
updated. The connection asks the server to abort a transaction left idle for
IDLE_IN_TRANSACTION_TIMEOUT, which is what a transaction becomes when the
machine running it is reset: without it, the server holds its locks until TCP
gives up on the connection, hours later, and a run started again waits behind
them. TCP keepalives are asked for on the client's side too, so that a run does
not wait for ever on a server that has gone.

Running the update again once it has committed -- when a crash left it unknown
whether the commit landed -- does what it did the first time, and nothing more:
a genome the report calls new or changed that the database already holds with
these very files was put there by the earlier run (already updated), keeps its
aligned markers and is not rewritten, and the genomes deleted are no longer
there to delete.

The hashes are appended as they are made to update_db_hashes.tsv.gz in
--out_dir (HashCache), with each file's size and modification time, and read
back by the next run, which hashes only the files not in it or changed since:
an interrupted run, or a run after a dry run, does not spend the hours again.
"""

import gzip
import logging
import multiprocessing as mp
import os
import time
import zlib
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from contextlib import ExitStack
from datetime import datetime
from typing import (AbstractSet, Dict, Iterable, List, NamedTuple, Optional,
                    Sequence, Set, Tuple)

import psycopg2
from psycopg2.extras import execute_values
from tqdm import tqdm

from gtdb_migration_tk.batching import read_genome_dirs
from gtdb_migration_tk.biolib_lite.checksum import sha256_rb
from gtdb_migration_tk.biolib_lite.common import canonical_gid, make_sure_path_exists
from gtdb_migration_tk.ncbi_utils import NCBI_DATABASES
from gtdb_migration_tk.update_genomes import (STATUS_FASTA_CHANGED,
                                             STATUS_FASTA_UNCHANGED,
                                             STATUS_IN_RELEASE, STATUS_NEW,
                                             STATUS_REGENERATE,
                                             STATUS_REMOVED,
                                             STATUS_SEQUENCES_UNCHANGED,
                                             STATUS_TO_CURATE, count_by_database,
                                             database_of, genomic_fasta,
                                             report_outcomes)
from gtdb_migration_tk.utils.common import protein_fasta

# What became of each genome of the database, as the report names it.
ACTION_ADDED = 'added'
ACTION_VERSIONED = 'versioned'
ACTION_SEQUENCES_CHANGED = 'sequences changed'
ACTION_UPDATED = 'updated'
ACTION_UNCHANGED = 'unchanged'
ACTION_DELETED = 'deleted'
ACTION_REPLACED = 'replaced by new version'
ACTION_NOT_IN_DATABASE = 'not in database'
ACTION_NOT_IN_REPORT = 'not in report'
ACTION_NO_GENOMIC_FASTA = 'no genomic FASTA'
ACTION_ALREADY_UPDATED = 'already updated'

# The actions after which a genome's sequences are new to the database: the
# genomes has_changed is TRUE for, and whose aligned markers have to be made.
ACTIONS_CHANGING_SEQUENCES = frozenset((ACTION_ADDED, ACTION_VERSIONED,
                                        ACTION_SEQUENCES_CHANGED))

# Every outcome a report can record, and so every one decided below. A report
# naming another is refused: an outcome added to update_genomes has to be given
# a meaning here before update_db acts on a release holding it.
KNOWN_OUTCOMES = frozenset((STATUS_NEW, STATUS_FASTA_CHANGED,
                            STATUS_SEQUENCES_UNCHANGED, STATUS_FASTA_UNCHANGED,
                            STATUS_REMOVED, STATUS_TO_CURATE))

# A genome directory of a release, as the database records its files: the path
# below the release, from the accession's prefix down, e.g.
# GCA/024/206/075/GCA_024206075.2_ASM2420607v2. The prefix is this many
# components above the directory, however deep the release itself lies.
DATABASE_PATH_DEPTH = 5

REPORT_NAME = 'update_db_report.tsv'
REPORT_HEADER = ('genome_id', 'report_outcome', 'action', 'detail')
LISTS_AFFECTED_NAME = 'genome_lists_affected.tsv'
LISTS_AFFECTED_HEADER = ('genome_id', 'list_id', 'list_name', 'run_started', 'dry_run')

# The hashes of earlier runs, and how often the ones being made are put on disk.
HASH_CACHE_NAME = 'update_db_hashes.tsv.gz'
HASH_CACHE_FLUSH_SECONDS = 30.0

# How long the server lets this update's transaction sit idle before aborting it
# and releasing its locks: the gaps between its statements are seconds, and a
# transaction idle for longer is one whose client has gone. TCP keepalives ask
# the same of the connection from the client's side.
IDLE_IN_TRANSACTION_TIMEOUT = '10min'
KEEPALIVES = dict(keepalives=1, keepalives_idle=60, keepalives_interval=15, keepalives_count=4)

# Rows sent to the database per statement.
PAGE_SIZE = 1000

# The most accessions an error message names; the rest are counted.
NAMED_GENOMES = 10


class UpdateDbError(RuntimeError):
    """A release update_db refuses to put into the database."""


class DatabaseGenome(NamedTuple):
    """A genome as the genomes table holds it."""

    id: int
    name: str
    source_id: int
    fasta_location: str
    fasta_hash: str
    genes_location: Optional[str]
    genes_hash: Optional[str]
    date_added: Optional[datetime] = None


class Decision(NamedTuple):
    """What update_db does with one genome, before any file is read."""

    accession: str
    outcome: str
    action: str
    row: Optional[DatabaseGenome] = None    # the row acted on: its own, or its predecessor's
    genome_dir: Optional[str] = None
    hash_genomic: bool = False
    hash_genes: bool = False
    detail: str = ''


class GenomeFiles(NamedTuple):
    """Where a genome's files are, as the database records them, and their hashes."""

    fasta_location: str
    fasta_hash: Optional[str]
    genes_location: Optional[str]
    genes_hash: Optional[str]


def name_genomes(accessions: Iterable[str]) -> str:
    """Genomes as an error names them: the first NAMED_GENOMES, the rest counted.

    @return: e.g. 'GCA_000000001.1, GCA_000000002.1 and 3 more'.
    """

    accessions = sorted(accessions)
    named = ', '.join(accessions[:NAMED_GENOMES])
    if len(accessions) > NAMED_GENOMES:
        named += ' and {:,} more'.format(len(accessions) - NAMED_GENOMES)
    return named


def versionless(accession: str) -> str:
    """An accession without its version, e.g. GCA_001317685 for GCA_001317685.3.

    The prefix is kept, so a RefSeq genome is never taken for a version of its
    GenBank counterpart.

    @return: the accession up to its last '.'.
    """

    return accession.rsplit('.', 1)[0]


def database_location(path: str, genome_dir: str, accession: str) -> str:
    """A file of a genome directory, as the database records where it is.

    Parameters
    ----------
    path : str
        The file, inside genome_dir.
    genome_dir : str
        Genome directory of the release, as genome_dirs.tsv gives it.
    accession : str
        The genome's accession.

    @return: the path below the release, e.g.
             GCA/024/206/075/GCA_024206075.2_ASM2420607v2/GCA_024206075.2_ASM2420607v2_genomic.fna.gz.

    Raises
    ------
    UpdateDbError
        The directory is not laid out as NCBI lays out a genome, under its
        accession's prefix and three levels of its digits.
    """

    parts = os.path.normpath(genome_dir).split(os.sep)
    if (len(parts) < DATABASE_PATH_DEPTH
            or parts[-DATABASE_PATH_DEPTH] != accession[:3]
            or not parts[-1].startswith(accession)):
        raise UpdateDbError(
            '{} is not where a release keeps {}: its directory should be '
            '<release>/<database>/{}/<digits>/<digits>/<digits>/{}_<assembly>.'.format(
                genome_dir, accession, accession[:3], accession))

    root = os.sep.join(parts[:-DATABASE_PATH_DEPTH])
    return os.path.relpath(path, root)


def content_hash(path: str) -> Optional[str]:
    """The hash the database records for a gzipped file: SHA-1 of its content.

    Parameters
    ----------
    path : str
        A gzipped file.

    @return: the hex digest, or None where there is no file.
    """

    try:
        raw = open(path, 'rb')
    except FileNotFoundError:
        return None
    # sha256_rb() closes the GzipFile it is handed, which leaves open the file
    # beneath it
    with raw:
        return sha256_rb(gzip.GzipFile(fileobj=raw))


def file_identity(path: str) -> Optional[Tuple[int, int]]:
    """What says a file is the one a hash was made of: its size and modification time.

    @return: (size in bytes, modification time in ns), or None where there is no file.
    """

    try:
        stat = os.stat(path)
    except FileNotFoundError:
        return None
    return stat.st_size, stat.st_mtime_ns


def hash_file(path: str) -> Tuple[str, Optional[int], Optional[int], Optional[str]]:
    """Hash one file; a worker of the pool.

    The file is identified before it is read, so that one rewritten while it is
    hashed is recorded with the older time and hashed again by the next run.

    @return: (path, size, modification time in ns, hash), the last three None
             where there is no file.
    """

    identity = file_identity(path)
    if identity is None:
        return path, None, None, None
    return (path,) + identity + (content_hash(path),)


class HashCache(object):
    """The hashes earlier runs made, so that a run started again need not remake them.

    One gzipped TSV of path, size, modification time (ns) and hash. A hash is
    taken from it only while the file's size and modification time are those it
    was made of. Hashes are appended as they are made, as a gzip member of their
    own, and put on disk every HASH_CACHE_FLUSH_SECONDS, so a run stopped hard
    loses that much at most. The member a crash cut short is read as far as it
    goes, and the file is written again whole when it is read, so that what is
    appended next follows a gzip file that ends where it should.
    """

    def __init__(self, path: str) -> None:
        self.path = path
        self.entries: Dict[str, Tuple[int, int, str]] = {}
        self.raw = None
        self.handle = None
        self.flushed = 0.0

    def load(self) -> bool:
        """Read what earlier runs hashed, and write it back whole.

        @return: True where the file ended part way through, as a crash leaves it.
        """

        if not os.path.exists(self.path):
            return False

        cut_short = False
        try:
            with gzip.open(self.path, 'rt') as handle:
                for line in handle:
                    if not line.endswith('\n'):
                        cut_short = True
                        break
                    fields = line.rstrip('\n').split('\t')
                    try:
                        self.entries[fields[0]] = (int(fields[1]), int(fields[2]), fields[3])
                    except (IndexError, ValueError):
                        cut_short = True
        except (EOFError, OSError, zlib.error):
            cut_short = True

        self.rewrite()
        return cut_short

    def rewrite(self) -> None:
        """Write every entry to a new file, and put it in place of the old one."""

        temporary = self.path + '.tmp'
        with open(temporary, 'wb') as raw:
            with gzip.GzipFile(fileobj=raw, mode='wb') as handle:
                for path, (size, mtime, digest) in self.entries.items():
                    handle.write('{}\t{}\t{}\t{}\n'.format(path, size, mtime, digest).encode())
            raw.flush()
            os.fsync(raw.fileno())
        os.replace(temporary, self.path)

    def lookup(self, path: str, identity: Optional[Tuple[int, int]]) -> Optional[str]:
        """The hash of a file, where it is the file the hash was made of.

        @return: the hash, or None.
        """

        entry = self.entries.get(path)
        if entry is None or identity is None or entry[:2] != identity:
            return None
        return entry[2]

    def __enter__(self) -> 'HashCache':
        self.raw = open(self.path, 'ab')
        self.handle = gzip.GzipFile(fileobj=self.raw, mode='wb')
        self.flushed = time.monotonic()
        return self

    def add(self, path: str, size: int, mtime: int, digest: str) -> None:
        """Record a hash just made, putting it on disk with the others in time."""

        self.entries[path] = (size, mtime, digest)
        self.handle.write('{}\t{}\t{}\t{}\n'.format(path, size, mtime, digest).encode())
        if time.monotonic() - self.flushed >= HASH_CACHE_FLUSH_SECONDS:
            self.flush()

    def flush(self) -> None:
        # Z_SYNC_FLUSH: what is on disk can be decompressed up to here
        self.handle.flush(zlib.Z_SYNC_FLUSH)
        os.fsync(self.raw.fileno())
        self.flushed = time.monotonic()

    def __exit__(self, *exc) -> None:
        self.handle.close()
        self.raw.flush()
        os.fsync(self.raw.fileno())
        self.raw.close()


def missing_genomic_fasta(accessions: Iterable[str], genome_dirs: Dict[str, str],
                          threads: int = 1) -> Set[str]:
    """The genomes whose release directory holds no genomic FASTA.

    One stat each, on threads, being NFS round trips.

    Parameters
    ----------
    accessions : iterable of str
        Genomes to look at.
    genome_dirs : dict
        accession -> genome directory.
    threads : int
        Stats made at once.

    @return: the accessions with no genomic FASTA.
    """

    accessions = list(accessions)
    with ThreadPoolExecutor(max_workers=max(1, threads)) as pool:
        present = pool.map(lambda acc: os.path.isfile(genomic_fasta(genome_dirs[acc])), accessions)
        return {acc for acc, there in zip(accessions, present) if not there}


def plan_update(outcomes: Dict[str, str],
                genome_dirs: Dict[str, str],
                database: Dict[str, DatabaseGenome],
                rehash_all: bool = False,
                no_genomic_fasta: AbstractSet[str] = frozenset()) -> List[Decision]:
    """Decide what becomes of every genome of the release and of the database.

    Nothing is read but the arguments, so the whole of what update_db decides is
    here, and testable without a database.

    Parameters
    ----------
    outcomes : dict
        accession -> outcome, as report_outcomes() reads report.log.
    genome_dirs : dict
        accession -> genome directory, from genome_dirs.tsv.
    database : dict
        name -> the NCBI genomes the database holds.
    rehash_all : bool
        Hash the files of every genome of the release.
    no_genomic_fasta : set
        Genomes the database does not hold whose directory has no genomic FASTA,
        as missing_genomic_fasta() finds them.

    @return: one decision per genome of the report, then one per genome of the
             database the report does not name.

    Raises
    ------
    UpdateDbError
        The report names an outcome not known here, or the report and
        genome_dirs.tsv disagree about which genomes the release holds.
    """

    unknown = Counter(outcome for outcome in outcomes.values() if outcome not in KNOWN_OUTCOMES)
    if unknown:
        raise UpdateDbError(
            'The report records outcomes update_db does not know what to do with: '
            '{}.'.format(', '.join('{} ({:,})'.format(o, n) for o, n in sorted(unknown.items()))))

    neither = [acc for acc in outcomes if database_of(acc) is None]
    if neither:
        raise UpdateDbError('{:,} accession(s) of the report belong to neither RefSeq nor '
                            'GenBank: {}.'.format(len(neither), name_genomes(neither)))

    in_release = {acc for acc, outcome in outcomes.items() if outcome in STATUS_IN_RELEASE}
    no_dir = in_release - genome_dirs.keys()
    not_in_report = genome_dirs.keys() - in_release
    if no_dir or not_in_report:
        raise UpdateDbError(
            'The report and the genome_dirs file are of different releases, or one '
            'is incomplete: {:,} genome(s) the report puts in the release have no '
            'directory ({}), and {:,} with a directory are not in the release by the '
            'report ({}).'.format(len(no_dir), name_genomes(no_dir) or 'none',
                                  len(not_in_report), name_genomes(not_in_report) or 'none'))

    # the database holds one version of an accession: GTDB has always replaced
    # a genome's row with its new version rather than kept both
    predecessors = {}
    for name, row in database.items():
        predecessors.setdefault(versionless(name), []).append(row)

    decisions = []
    replaced = {}
    for accession in sorted(in_release):
        outcome = outcomes[accession]
        genome_dir = genome_dirs[accession]
        row = database.get(accession)

        if row is None and accession in no_genomic_fasta:
            decisions.append(Decision(accession, outcome, ACTION_NO_GENOMIC_FASTA, None,
                                      genome_dir, detail='not added'))
        elif row is None:
            older = [r for r in predecessors.get(versionless(accession), ())
                     if outcomes.get(r.name) not in STATUS_IN_RELEASE and r.name not in replaced]
            if outcome == STATUS_NEW and len(older) == 1:
                replaced[older[0].name] = accession
                decisions.append(Decision(accession, outcome, ACTION_VERSIONED, older[0],
                                          genome_dir, True, True,
                                          'replaces {}'.format(older[0].name)))
            else:
                detail = '' if outcome == STATUS_NEW else 'in the release but not the database'
                decisions.append(Decision(accession, outcome, ACTION_ADDED, None,
                                          genome_dir, True, True, detail))
        elif outcome in (STATUS_NEW, STATUS_FASTA_CHANGED):
            decisions.append(Decision(accession, outcome, ACTION_SEQUENCES_CHANGED, row,
                                      genome_dir, True, True))
        else:
            # deflines rewritten: the genomic FASTA is a different file, the
            # proteins are carried across as they were
            hash_genomic = rehash_all or outcome == STATUS_SEQUENCES_UNCHANGED
            hash_genes = rehash_all or row.genes_hash is None
            decisions.append(Decision(accession, outcome, ACTION_UPDATED, row,
                                      genome_dir, hash_genomic, hash_genes))

    for accession in sorted(set(outcomes) - in_release):
        outcome = outcomes[accession]
        row = database.get(accession)
        if row is None:
            decisions.append(Decision(accession, outcome, ACTION_NOT_IN_DATABASE))
        elif accession in replaced:
            decisions.append(Decision(accession, outcome, ACTION_REPLACED, row,
                                      detail='by {}'.format(replaced[accession])))
        else:
            decisions.append(Decision(accession, outcome, ACTION_DELETED, row))

    for name in sorted(database.keys() - outcomes.keys()):
        decisions.append(Decision(name, '', ACTION_NOT_IN_REPORT, database[name]))

    return decisions


def genome_files(decision: Decision,
                 hashes: Dict[str, Tuple[Optional[str], Optional[str]]]) -> GenomeFiles:
    """Where a genome of the release has its files, as the database records them.

    A file not hashed keeps the hash the database holds for it; a protein file
    that is not there is recorded as none.

    Parameters
    ----------
    decision : Decision
        A genome of the release.
    hashes : dict
        accession -> (genomic hash, protein hash) of the genomes hashed.

    @return: the locations and hashes to record.
    """

    accession, genome_dir, row = decision.accession, decision.genome_dir, decision.row
    genomic_hash, genes_hash = hashes.get(accession, (None, None))

    fasta_location = database_location(genomic_fasta(genome_dir), genome_dir, accession)
    if not decision.hash_genomic:
        genomic_hash = row.fasta_hash

    genes_location = database_location(protein_fasta(accession, genome_dir), genome_dir, accession)
    if not decision.hash_genes:
        genes_hash = row.genes_hash
    if genes_hash is None:
        genes_location = None

    return GenomeFiles(fasta_location, genomic_hash, genes_location, genes_hash)


def brought_by_this_update(decision: Decision, download_date: datetime) -> bool:
    """Whether a genome the database holds is one this release brought.

    One the report calls new or changed, or one whose date_added is this update's
    download date, which an earlier run of the same update added: has_changed is
    TRUE for it however many times the update is run.

    @return: True for such a genome.
    """

    if decision.outcome in STATUS_REGENERATE:
        return True
    added = decision.row.date_added if decision.row is not None else None
    return added is not None and added.date() == download_date.date()


def changed_columns(row: DatabaseGenome, files: GenomeFiles) -> List[str]:
    """The columns of a genome's row that its files no longer agree with.

    @return: the names, in the order of the table, e.g. ['fasta_file_sha256'].
    """

    changed = []
    for column, old, new in (('fasta_file_location', row.fasta_location, files.fasta_location),
                             ('fasta_file_sha256', row.fasta_hash, files.fasta_hash),
                             ('genes_file_location', row.genes_location, files.genes_location),
                             ('genes_file_sha256', row.genes_hash, files.genes_hash)):
        if old != new:
            changed.append(column)
    return changed


class DatabaseManager(object):
    """Bring the NCBI genomes of the GTDB database into line with a release."""

    def __init__(self, hostname: str, user: str, password: str, db: str, cpus: int = 1) -> None:
        """Initialization.

        Parameters
        ----------
        hostname, user, password, db : str
            How to reach the database. They are handed to psycopg2 as they are
            rather than written into a connection string, so a password holding a
            space or a quote is not mangled.
        cpus : int
            Processes the files are hashed on.

        @return: None
        """

        self.connection = dict(host=hostname, user=user, password=password, dbname=db)
        self.cpus = cpus
        self.logger = logging.getLogger('timestamp')

    def connect(self):
        """A connection whose abandoned transaction the server aborts; see STOPPED PART WAY."""
        return psycopg2.connect(
            **self.connection, **KEEPALIVES,
            options='-c idle_in_transaction_session_timeout={}'.format(IDLE_IN_TRANSACTION_TIMEOUT))

    def ncbi_sources(self, cur) -> Dict[str, int]:
        """The genome_sources id of each NCBI database, by accession prefix.

        @return: e.g. {'GCF': 2, 'GCA': 3}.
        """

        cur.execute('SELECT id, name FROM genome_sources')
        by_name = {name: source_id for source_id, name in cur.fetchall()}
        missing = [db.label for db in NCBI_DATABASES if db.label not in by_name]
        if missing:
            raise UpdateDbError('genome_sources names no source {}; it names {}.'.format(
                ' or '.join(missing), ', '.join(sorted(by_name))))
        return {db.prefix: by_name[db.label] for db in NCBI_DATABASES}

    def database_genomes(self, cur, source_ids: Sequence[int]) -> Dict[str, DatabaseGenome]:
        """The NCBI genomes the database holds, by name.

        @return: name -> DatabaseGenome.
        """

        cur.execute('SELECT id, name, genome_source_id, fasta_file_location, '
                    'fasta_file_sha256, genes_file_location, genes_file_sha256, date_added '
                    'FROM genomes WHERE genome_source_id = ANY(%s)', (list(source_ids),))
        return {row[1]: DatabaseGenome(*row) for row in cur.fetchall()}

    def hash_files(self, decisions: Sequence[Decision],
                   cache_file: Optional[str] = None) -> Dict[str, Tuple[Optional[str], Optional[str]]]:
        """Hash the files the decisions ask for, on self.cpus processes.

        Parameters
        ----------
        decisions : sequence of Decision
            What is to be done with each genome.
        cache_file : str
            The HashCache of earlier runs, read and added to; None for none.

        @return: accession -> (genomic hash, protein hash).

        Raises
        ------
        UpdateDbError
            A genome of the release has no genomic FASTA to hash.
        """

        # (accession, 0 for the genomic FASTA and 1 for the proteins, path)
        files = []
        for d in decisions:
            if d.hash_genomic:
                files.append((d.accession, 0, genomic_fasta(d.genome_dir)))
            if d.hash_genes:
                files.append((d.accession, 1, protein_fasta(d.accession, d.genome_dir)))
        self.logger.info('Hashing the files of {:,} genome(s): {:,} genomic and {:,} '
                         'protein FASTA file(s).'.format(
                             len({f[0] for f in files}), sum(1 for f in files if f[1] == 0),
                             sum(1 for f in files if f[1] == 1)))
        if not files:
            return {}

        digests = {}
        todo = [path for _, _, path in files]
        cache = HashCache(cache_file) if cache_file else None
        if cache is not None:
            if cache.load():
                self.logger.warning(
                    'warning: {} ends part way through, as a run stopped hard leaves it; '
                    'the {:,} hash(es) before that are used.'.format(cache_file, len(cache.entries)))
            if cache.entries:
                # one stat a file, on threads, to tell a file hashed before from one
                # rewritten since
                with ThreadPoolExecutor(max_workers=max(1, self.cpus)) as pool:
                    identities = dict(zip(todo, pool.map(file_identity, todo)))
                for path in todo:
                    digest = cache.lookup(path, identities[path])
                    if digest is not None:
                        digests[path] = digest
                todo = [path for path in todo if path not in digests]
                self.logger.info('{:,} of the {:,} file(s) were hashed by an earlier run and '
                                 'are taken from {}, leaving {:,} to hash.'.format(
                                     len(digests), len(files), cache_file, len(todo)))

        if todo:
            with ExitStack() as stack:
                if cache is not None:
                    stack.enter_context(cache)
                pool = stack.enter_context(mp.Pool(self.cpus))
                for path, size, mtime, digest in tqdm(
                        pool.imap_unordered(hash_file, todo, chunksize=64),
                        total=len(todo), ncols=100, leave=False, desc='Hashing files'):
                    digests[path] = digest
                    if cache is not None and digest is not None:
                        cache.add(path, size, mtime, digest)

        found = {}
        for accession, which, path in files:
            found.setdefault(accession, [None, None])[which] = digests.get(path)
        hashes = {accession: tuple(pair) for accession, pair in found.items()}

        missing = [d.accession for d in decisions
                   if d.hash_genomic and hashes[d.accession][0] is None]
        if missing:
            # a genome to be added was looked at before planning, so this is one
            # the database holds: the release has lost its file
            raise UpdateDbError(
                '{:,} genome(s) the database holds have no genomic FASTA in the directory '
                'genome_dirs gives them, so the release has lost the file: {}.'.format(
                    len(missing), name_genomes(missing)))

        return hashes

    def run(self,
            genome_dirs_file: str,
            report_file: str,
            download_date: datetime,
            out_dir: str,
            rehash_all: bool = False,
            dry_run: bool = False) -> None:
        """Update the NCBI genomes of the database to a release.

        Parameters
        ----------
        genome_dirs_file : str
            genome_dirs.tsv of the release, as update_genomes wrote it.
        report_file : str
            report.log of the release, as update_genomes wrote it.
        download_date : datetime
            When the release's genomes were taken from NCBI: date_added of a
            genome added or versioned, and last_update of every row changed.
        out_dir : str
            Directory the report of the update is written to.
        rehash_all : bool
            Hash the files of every genome of the release.
        dry_run : bool
            Make every change and roll it back.

        @return: None
        """

        make_sure_path_exists(out_dir)
        run_started = datetime.now().isoformat(timespec='seconds')
        # refused now rather than after hours of hashing
        self.check_lists_affected(out_dir)
        outcomes = report_outcomes(report_file)
        genome_dirs = dict(read_genome_dirs(genome_dirs_file))
        self.logger.info('The report names {:,} genome(s) ({}); the release holds {:,}.'.format(
            len(outcomes), count_by_database(outcomes), len(genome_dirs)))

        # read, then closed while the files are hashed: hours, which a connection
        # held open is not worth, and the pool's processes are forked without it
        with self.connect() as conn, conn.cursor() as cur:
            sources = self.ncbi_sources(cur)
            database = self.database_genomes(cur, sources.values())
        conn.close()
        self.logger.info('The database holds {:,} NCBI genome(s) ({}).'.format(
            len(database), count_by_database(database)))

        # a genome the database would gain needs a sequence to record; one it
        # holds is checked when its files are hashed
        unpublished = missing_genomic_fasta(
            [acc for acc in genome_dirs if acc not in database], genome_dirs, self.cpus)
        decisions = plan_update(outcomes, genome_dirs, database, rehash_all, unpublished)
        hashes = self.hash_files(decisions, os.path.join(out_dir, HASH_CACHE_NAME))
        decisions = self.finish_decisions(decisions, hashes)

        conn = self.connect()
        try:
            with conn.cursor() as cur:
                self.write_lists_affected(cur, decisions, out_dir, run_started, dry_run)
                self.apply(cur, decisions, hashes, sources, download_date)
            if dry_run:
                conn.rollback()
                self.logger.info('Dry run: every change was made and rolled back.')
            else:
                conn.commit()
                self.logger.info('Committed the update begun {}.'.format(run_started))
        except BaseException:
            conn.rollback()
            raise
        finally:
            conn.close()

        self.write_report(decisions, out_dir)
        self.log_counts(decisions, dry_run)

    def finish_decisions(self, decisions: Sequence[Decision],
                         hashes: Dict[str, Tuple[Optional[str], Optional[str]]]) -> List[Decision]:
        """Tell an updated genome from an unchanged one, now its files are hashed.

        @return: the decisions, an update that changes nothing made ACTION_UNCHANGED
                 and the rest given the columns they change and any missing file
                 as their detail.
        """

        finished = []
        for decision in decisions:
            if decision.genome_dir is None or decision.action == ACTION_NO_GENOMIC_FASTA:
                finished.append(decision)
                continue

            files = genome_files(decision, hashes)
            notes = [decision.detail] if decision.detail else []
            if files.genes_hash is None:
                notes.append('no protein file')

            if decision.action == ACTION_UPDATED:
                changed = changed_columns(decision.row, files)
                if not changed:
                    finished.append(decision._replace(action=ACTION_UNCHANGED,
                                                      detail='; '.join(notes)))
                    continue
                notes.append(', '.join(changed))
            elif (decision.action == ACTION_SEQUENCES_CHANGED
                  and not changed_columns(decision.row, files)):
                # the database already holds these very files: an earlier run of
                # this update put them there
                notes.append('by an earlier run of this update')
                finished.append(decision._replace(action=ACTION_ALREADY_UPDATED,
                                                  detail='; '.join(notes)))
                continue

            finished.append(decision._replace(detail='; '.join(notes)))

        return finished

    def check_lists_affected(self, out_dir: str) -> None:
        """Refuse a genome_lists_affected.tsv whose columns are not this version's.

        Rows are appended to it, so one of another layout would end up holding
        two; checked where the run starts rather than once the files are hashed.

        @return: None
        """

        path = os.path.join(out_dir, LISTS_AFFECTED_NAME)
        if not os.path.exists(path) or os.path.getsize(path) == 0:
            return
        with open(path) as handle:
            header = handle.readline().rstrip('\n').split('\t')
        if tuple(header) != LISTS_AFFECTED_HEADER:
            raise UpdateDbError(
                '{} has the columns {}, where this version appends rows of {}; move it '
                'aside to keep it.'.format(path, ', '.join(header), ', '.join(LISTS_AFFECTED_HEADER)))

    def write_lists_affected(self, cur, decisions: Sequence[Decision], out_dir: str,
                             run_started: str, dry_run: bool) -> None:
        """Name the curated genome lists each genome to be deleted is in.

        Appended, so that what one run records survives the next, and put on disk
        before the transaction is committed, so that a crash at the commit cannot
        lose it: a run that did not commit is told from one that did by its log.

        @return: None
        """

        ids = {d.row.id: d.accession for d in decisions if d.action == ACTION_DELETED}
        rows = []
        if ids:
            cur.execute('SELECT glc.genome_id, gl.id, gl.name FROM genome_list_contents glc '
                        'JOIN genome_lists gl ON gl.id = glc.list_id '
                        'WHERE glc.genome_id = ANY(%s)', (list(ids),))
            rows = sorted((ids[genome_id], str(list_id), name, run_started, str(dry_run))
                          for genome_id, list_id, name in cur.fetchall())

        path = os.path.join(out_dir, LISTS_AFFECTED_NAME)
        new_file = not os.path.exists(path) or os.path.getsize(path) == 0
        with open(path, 'a') as handle:
            if new_file:
                handle.write('\t'.join(LISTS_AFFECTED_HEADER) + '\n')
            for row in rows:
                handle.write('\t'.join(row) + '\n')
            handle.flush()
            os.fsync(handle.fileno())

        if rows:
            self.logger.warning(
                'warning: {:,} genome(s) to be deleted are in {:,} curated genome list(s), '
                'which lose them; each is named in {}.'.format(
                    len({r[0] for r in rows}), len({r[1] for r in rows}), path))

    def apply(self, cur, decisions: Sequence[Decision],
              hashes: Dict[str, Tuple[Optional[str], Optional[str]]],
              sources: Dict[str, int], download_date: datetime) -> None:
        """Make the decisions in the database, inside the caller's transaction.

        @return: None
        """

        by_action = {}
        for decision in decisions:
            by_action.setdefault(decision.action, []).append(decision)

        # has_changed says what this update changed and nothing older
        cur.execute('UPDATE genomes SET has_changed = FALSE '
                    'WHERE genome_source_id = ANY(%s) AND has_changed', (list(sources.values()),))

        deleted = [d.row.id for d in by_action.get(ACTION_DELETED, ())]
        if deleted:
            cur.execute('DELETE FROM genomes WHERE id = ANY(%s)', (deleted,))

        versioned = []
        for d in by_action.get(ACTION_VERSIONED, ()):
            f = genome_files(d, hashes)
            versioned.append((d.row.id, d.accession, sources[d.accession[:3]],
                              canonical_gid(d.accession), f.fasta_location, f.fasta_hash,
                              f.genes_location, f.genes_hash, download_date))
        if versioned:
            execute_values(
                cur,
                'UPDATE genomes AS g SET name = v.name, id_at_source = v.name, '
                'genome_source_id = v.source, formatted_source_id = v.formatted, '
                'fasta_file_location = v.fasta, fasta_file_sha256 = v.fasta_hash, '
                'genes_file_location = v.genes, genes_file_sha256 = v.genes_hash, '
                'date_added = v.added, last_update = v.added, has_changed = TRUE '
                'FROM (VALUES %s) AS v(id, name, source, formatted, fasta, fasta_hash, '
                'genes, genes_hash, added) WHERE g.id = v.id',
                versioned,
                template='(%s::integer, %s::text, %s::integer, %s::text, %s::text, %s::text, '
                         '%s::text, %s::text, %s::timestamp)',
                page_size=PAGE_SIZE)

        updated = []
        for action, has_changed in ((ACTION_SEQUENCES_CHANGED, True), (ACTION_UPDATED, False)):
            for d in by_action.get(action, ()):
                f = genome_files(d, hashes)
                updated.append((d.row.id, f.fasta_location, f.fasta_hash,
                                f.genes_location, f.genes_hash, download_date, has_changed))
        if updated:
            execute_values(
                cur,
                'UPDATE genomes AS g SET fasta_file_location = v.fasta, '
                'fasta_file_sha256 = v.fasta_hash, genes_file_location = v.genes, '
                'genes_file_sha256 = v.genes_hash, last_update = v.updated, '
                'has_changed = v.has_changed '
                'FROM (VALUES %s) AS v(id, fasta, fasta_hash, genes, genes_hash, updated, '
                'has_changed) WHERE g.id = v.id',
                updated,
                template='(%s::integer, %s::text, %s::text, %s::text, %s::text, %s::date, '
                         '%s::boolean)',
                page_size=PAGE_SIZE)

        remade = ([d.row.id for d in by_action.get(ACTION_VERSIONED, ())]
                  + [d.row.id for d in by_action.get(ACTION_SEQUENCES_CHANGED, ())])
        if remade:
            cur.execute('DELETE FROM aligned_markers WHERE genome_id = ANY(%s)', (remade,))

        # a genome this release brought that the database already holds as it
        # should: put there by an earlier run of this same update
        brought = [d.row.id for d in decisions
                   if d.action in (ACTION_ALREADY_UPDATED, ACTION_UPDATED, ACTION_UNCHANGED)
                   and brought_by_this_update(d, download_date)]
        if brought:
            cur.execute('UPDATE genomes SET has_changed = TRUE WHERE id = ANY(%s)', (brought,))

        added = []
        for d in by_action.get(ACTION_ADDED, ()):
            f = genome_files(d, hashes)
            added.append((d.accession, '', True, None, f.fasta_location, f.fasta_hash,
                          sources[d.accession[:3]], d.accession, download_date, True,
                          download_date, f.genes_location, f.genes_hash,
                          canonical_gid(d.accession)))
        if added:
            execute_values(
                cur,
                'INSERT INTO genomes (name, description, owned_by_root, owner_id, '
                'fasta_file_location, fasta_file_sha256, genome_source_id, id_at_source, '
                'date_added, has_changed, last_update, genes_file_location, '
                'genes_file_sha256, formatted_source_id) VALUES %s',
                added, page_size=PAGE_SIZE)

    def write_report(self, decisions: Sequence[Decision], out_dir: str) -> None:
        """Write what became of every genome, one row each.

        @return: None
        """

        path = os.path.join(out_dir, REPORT_NAME)
        with open(path, 'w') as handle:
            handle.write('\t'.join(REPORT_HEADER) + '\n')
            for d in decisions:
                handle.write('\t'.join((d.accession, d.outcome, d.action, d.detail)) + '\n')
        self.logger.info('Wrote what became of each genome to {}.'.format(path))

    def log_counts(self, decisions: Sequence[Decision], dry_run: bool) -> None:
        """Say how many genomes each action took, by database.

        @return: None
        """

        by_action = {}
        for d in decisions:
            by_action.setdefault(d.action, []).append(d.accession)

        self.logger.info('What the update {} to each genome:'.format(
            'would have done' if dry_run else 'did'))
        for action in (ACTION_ADDED, ACTION_VERSIONED, ACTION_SEQUENCES_CHANGED,
                       ACTION_UPDATED, ACTION_UNCHANGED, ACTION_DELETED,
                       ACTION_REPLACED, ACTION_NOT_IN_DATABASE, ACTION_NOT_IN_REPORT,
                       ACTION_NO_GENOMIC_FASTA, ACTION_ALREADY_UPDATED):
            accessions = by_action.get(action, ())
            self.logger.info('  {}: {:,} ({}).'.format(
                action, len(accessions), count_by_database(accessions)))

        unpublished = by_action.get(ACTION_NO_GENOMIC_FASTA, ())
        if unpublished:
            self.logger.warning(
                'warning: {:,} genome(s) of the release have no genomic FASTA and were not '
                'added: {}. NCBI publishes some assemblies without their sequence; such a '
                'genome is added by a later release, once it is published.'.format(
                    len(unpublished), name_genomes(unpublished)))

        missing_genes = [d.accession for d in decisions if 'no protein file' in d.detail]
        if missing_genes:
            self.logger.warning(
                'warning: {:,} genome(s) of the release have no protein file and are '
                'recorded with none: {}.'.format(len(missing_genes), name_genomes(missing_genes)))

        not_in_report = by_action.get(ACTION_NOT_IN_REPORT, ())
        if not_in_report:
            self.logger.warning(
                'warning: {:,} genome(s) of the database are not named by the report and '
                'were left as they are: {}.'.format(len(not_in_report), name_genomes(not_in_report)))
