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
batching.py -- cutting a release into batches and sharing them between machines.

A release is a million-odd genomes and the work over it takes days, so it is cut
into batches and several machines are pointed at one output directory to divide
them. That is the same problem whichever command is doing the work, and this is
the one implementation of it: trans_table predicts a translation table for each
genome of a batch, prodigal calls its genes, and neither has any business owning
a second copy of the rules about who holds what.

What differs between them is what a batch is FOR and what it leaves behind, which
is why nothing here knows about gTranslate or Prodigal. What differs only in
name -- what the batchfile is called, what the batch's own log is called --
travels in a BatchLayout, so that a batch directory of either command is read by
this code and neither command's files are named by it.

BATCHES
The batches are settled before any of them is processed, from the genomes sorted
by accession, so which genomes are in batch N follows from the set of genomes and
not from the order the genome_dirs file happens to be written in. Once the
batchfiles exist they are authoritative: a later run reuses them, since
partitioning a release again that has gained a genome would move genomes between
batches that are already finished.

CANARIES
The state of a batch is held in files in the batch's own directory, so that
another machine can see it and so that it outlives the process that wrote it.
RUNNING says a machine is working on the batch and names it; SUCCESS says the
batch is finished and its results are complete; FAILED says it was attempted and
something went wrong, and what an earlier attempt said is kept beside it.

THE CLAIM IS A LEASE
A claim is made by linking a uniquely named file onto RUNNING rather than by
creating RUNNING directly: this runs against NFS, where O_EXCL has never been the
operation two machines can race on safely and link() is. The machine holding a
batch touches RUNNING every few minutes, and a claim untouched for the lease is
taken by whichever machine next comes to the batch, whatever host made it. A
claim of this host's whose process has gone is taken at once. The lease is
measured against the FILE SERVER's clock, not against each machine's, so the
machines sharing an output directory need not agree about the time.

THE LOG OF A BATCH IS KEPT WITH THE BATCH
Every machine writes what it does to its own --log, and several machines sharing
one log file over NFS do not append to it: they overwrite one another and leave
the file full of holes. So what happens to a batch is also written into the
batch's own directory, which no other machine writes to. Whatever became of the
run that started a batch, the batch says.
"""

import contextlib
import datetime
import gzip
import logging
import os
import shutil
import socket
import tempfile
import threading
import time
import uuid
from concurrent.futures import ThreadPoolExecutor
from typing import (Callable, Dict, Iterator, List, NamedTuple, Optional,
                    Sequence, Tuple, TypeVar)

from tqdm import tqdm

from gtdb_migration_tk.ncbi_utils import assembly_total_length, genomic_fasta
from gtdb_migration_tk.utils.common import MBP, open_text


class BatchLayout(NamedTuple):
    """What a command calls the files of its own batches.

    The mechanism is shared and the names are not: a trans_table batch directory
    and a prodigal one are read by the same code and are told apart by what is in
    them. batchfiles is in preference order, so that a batch planned by an older
    version, under a name since changed, is still found rather than looking
    unplanned -- a release whose batches looked unplanned would be partitioned
    again with its batches already done.
    """

    batchfiles: Tuple[str, ...]
    log: str


# Genomes per batch. A batch is the unit of restart and of sharing, so it is
# small enough that losing one is not a day and large enough that the claiming is
# nothing against the work.
DEFAULT_BATCH_SIZE = 10000

# Batch directories are numbered rather than named for the genomes they hold: the
# accessions of a batch are in its batchfile, and a name is a thing to sort by.
BATCH_DIR_PREFIX = 'batch_'
BATCH_DIR_FORMAT = BATCH_DIR_PREFIX + '{:06d}'

# The accessions of a batch whose genomic FASTA is not where the release says it
# is. Written only when there are some, so the file being there at all says a
# batch had something wrong with it.
MISSING_NAME = 'missing_genomic_fasta.tsv'

# The state of a batch. A command may write canaries of its own beside these --
# trans_table writes PREDICTED to say the hours of a batch are over -- but these
# three are what claiming and restarting are decided by.
RUNNING_CANARY = 'RUNNING'
SUCCESS_CANARY = 'SUCCESS'
FAILED_CANARY = 'FAILED'

# How often the machine holding a batch says it is still there, and how long a
# claim outlives the last thing said. The interval is small against the hours a
# batch takes and the lease is large against the interval, so a claim is freed
# only when a machine has really stopped saying anything.
HEARTBEAT_SECONDS = 300
CLAIM_LEASE_SECONDS = 2 * 60 * 60

# The extension a gzipped file of either command is named with, so that what a
# file IS and what it is called cannot drift apart.
GZIP_EXT = '.gz'

STATE_PENDING = 'pending'
STATE_RUNNING = 'running'
STATE_SUCCESS = 'success'
STATE_FAILED = 'failed'

# How many of the release's genomic FASTA files are asked about at once while the
# batches are planned. The question is one stat per genome and nothing else, so
# what it costs is round trips to the file server and not CPU. Measured against
# r237 over the NFS the genomes are held on, one at a time takes 12-23 ms a
# genome -- 4 to 8 hours for 1.35M genomes -- against 3-7 ms with 32 outstanding
# at once, so a few hours become one or two. Beyond 32 nothing further was
# measurable: the limit is the server and the load it is already under, not the
# number of threads asking. Threads rather than processes because a stat spends
# its time in the kernel waiting and brings back one number.
STAT_THREADS = 32

# How many genomes are handed to the pool at a time. Executor.map() submits every
# item it is given before the first result can be read, one future per genome, so
# the whole release at once builds a million-odd futures before a single answer
# comes back. A chunk is large enough that no thread waits for the next one to be
# cut and small enough to be nothing in memory.
STAT_CHUNK = 50000

# A genome in whatever form the command holding it keeps it.
T = TypeVar('T')


def read_genome_dirs(gtdb_genome_path_file: str) -> List[Tuple[str, str]]:
    """Read the genomes of a release from its genome_dirs file.

    The file is the headerless accession / directory / canonical accession TSV
    list_genomes and update_genomes write. Rows are split on tabs and further
    columns ignored, as every other reader of it does, so a column appended later
    does not reach here.

    Parameters
    ----------
    gtdb_genome_path_file : str
        genome_dirs file of the release.

    @return: (accession, genome directory) for each genome, in the order read.
    """

    genomes = []
    with open(gtdb_genome_path_file) as handle:
        # leave=False: the bar is worth having while a release of half a million
        # genomes is read and worth nothing afterwards, and a bar left behind
        # sits in the middle of the log saying what has already been reported
        for line in tqdm(handle, ncols=100, leave=False, desc='Reading genomes'):
            line = line.strip()
            if not line:
                continue
            tokens = line.split('\t')
            genomes.append((tokens[0], tokens[1]))

    return genomes


def fasta_size(fasta: str) -> int:
    """Size of a genomic FASTA, and 0 where there is no file to have one.

    One stat rather than os.path.exists() and os.path.getsize(), which ask the
    file server the same question twice; over NFS, and once per genome of a
    release, that is half the cost of planning a run. A file that has gone
    between the two calls also raises from the second, where here it is simply a
    genome with nothing to process.

    Parameters
    ----------
    fasta : str
        Path of the genomic FASTA.

    @return: size in bytes, or 0 if the file is absent or cannot be read.
    """

    try:
        return os.stat(fasta).st_size
    except OSError:
        return 0


def split_by_fasta(rows: Sequence[Tuple[str, str]],
                   threads: int = STAT_THREADS) -> Tuple[List[Tuple[str, str]], List[str]]:
    """Sort a batch's genomes into those that can be asked about and those that cannot.

    A genome is asked about only where its genomic FASTA is on disk and is not
    empty. gTranslate checks the paths of a batchfile before it starts and
    refuses the WHOLE batch if one of them is missing, so a single absent file
    would cost the other ten thousand genomes of the batch -- and would cost them
    again on every retry, the batch failing identically each time. Filtering here
    leaves the batch to run and the accession to be named.

    The files are asked about many at a time, the genomes being held over NFS
    where a stat is a round trip to a server rather than a lookup in a cache.
    Nothing is computed here to be divided up: what the pool is for is having many
    round trips outstanding at once rather than one.

    The answer does not depend on how many threads asked: results are read back in
    the order the genomes were given, which is the order gTranslate is handed them.

    Parameters
    ----------
    rows : sequence of tuple
        (FASTA path, accession) for the genomes of a batch.
    threads : int
        Stat calls to keep in flight at once.

    @return: (present, missing), the rows to hand gTranslate and the accessions
             of the genomes left out.
    """

    present, missing = [], []
    # leave=False, as the bar reading the genome_dirs file is: it is worth having
    # while the batch is checked over and worth nothing once it is
    with tqdm(total=len(rows), ncols=100, leave=False,
              desc='Checking genomes') as pbar, \
            ThreadPoolExecutor(max_workers=max(1, threads)) as pool:
        for start in range(0, len(rows), STAT_CHUNK):
            chunk = rows[start:start + STAT_CHUNK]
            for (fasta, accession), size in zip(
                    chunk, pool.map(fasta_size, [fasta for fasta, _ in chunk])):
                if size > 0:
                    present.append((fasta, accession))
                else:
                    missing.append(accession)
                pbar.update()

    return present, missing


def count_bases(fasta: str) -> int:
    """The number of bases in a genomic FASTA, gaps included.

    Counted a line at a time rather than with biolib_lite's reader, which calls
    sys.exit() on a file it cannot read -- from a thread of the pool that would
    be the end of the run rather than of one genome.

    Parameters
    ----------
    fasta : str
        Path of the gzipped genomic FASTA.

    @return: number of bases.
    """

    bases = 0
    with gzip.open(fasta, 'rt') as handle:
        for line in handle:
            if not line.startswith('>'):
                bases += len(line.strip())

    return bases


def genome_size(genome_dir: str) -> Optional[int]:
    """The size of a genome assembly, in bases.

    NCBI's assembly statistics are read where they are there, which is one small
    file; the genomic FASTA is counted only where they are not, since that means
    reading every base of it.

    Parameters
    ----------
    genome_dir : str
        Genome directory of the release.

    @return: the size in bases, or None where neither file can say.
    """

    size = assembly_total_length(genome_dir)
    if size is not None:
        return size

    fasta = genomic_fasta(genome_dir)
    if fasta_size(fasta) == 0:
        return None

    try:
        return count_bases(fasta)
    except (OSError, EOFError, UnicodeDecodeError):
        # a FASTA that cannot be read is not known to be too large; it is left to
        # the command, whose own failure path names it
        return None


def split_by_genome_size(items: Sequence[T],
                         genome_dir_of: Callable[[T], str],
                         max_bases: int,
                         logger: logging.Logger,
                         threads: int = STAT_THREADS) -> Tuple[List[T], List[Tuple[T, int]]]:
    """Leave out the genomes too large for a command to be given.

    A genome assembly of billions of bases is a metagenome deposited as one
    genome, and a command run over it for each of its genes -- tRNAscan-SE,
    nhmmer, a blastn per rRNA gene against SILVA and the LTP -- is not seconds but
    a day or more, which holds its whole batch unfinished for that long. It is
    left out and named rather than worked on.

    A genome whose size cannot be told is kept: what is left out is what is known
    to be too large.

    Parameters
    ----------
    items : sequence
        The genomes to consider, in whatever form the command holds them.
    genome_dir_of : callable
        The genome directory of an item.
    max_bases : int
        Largest assembly processed, in bases.
    logger : logging.Logger
        Where each genome left out is named, with its size.
    threads : int
        Genomes asked about at once, the files being held over NFS.

    @return: (kept, too_large), the items to process in the order given and the
             (item, size in bases) of those left out.
    """

    kept, too_large = [], []
    with tqdm(total=len(items), ncols=100, leave=False,
              desc='Checking genome sizes') as pbar, \
            ThreadPoolExecutor(max_workers=max(1, threads)) as pool:
        for start in range(0, len(items), STAT_CHUNK):
            chunk = items[start:start + STAT_CHUNK]
            for item, size in zip(
                    chunk, pool.map(genome_size, [genome_dir_of(item) for item in chunk])):
                if size is not None and size > max_bases:
                    too_large.append((item, size))
                else:
                    kept.append(item)
                pbar.update()

    for item, size in too_large:
        logger.warning('warning: {} is {:,.1f} Mbp, larger than the {:,.1f} Mbp '
                       'maximum genome size, and is not processed.'.format(
                           os.path.basename(genome_dir_of(item)),
                           size / MBP, max_bases / MBP))

    return kept, too_large


def write_batchfile(rows: Sequence[Tuple[str, str]], batchfile: str,
                    compress: bool = False) -> None:
    """Write the two-column batchfile of a batch.

    The genome ID given is the accession the genome_dirs file names, so every row
    of the prediction table can be matched back to the genome directory it came
    from without canonicalising anything.

    Parameters
    ----------
    rows : sequence of tuple
        (FASTA path, accession) for each genome to ask about.
    batchfile : str
        File to write.
    compress : bool
        Write it gzipped, which the batch's own plan is and the copy handed to
        gTranslate is not: gTranslate reads a batchfile with a plain open().

    @return: None
    """

    with (gzip.open(batchfile, 'wt') if compress else open(batchfile, 'w')) as handle:
        for fasta, accession in rows:
            handle.write('{}\t{}\n'.format(fasta, accession))


def read_batchfile(batchfile: str) -> List[Tuple[str, str]]:
    """Read back the genomes of a batch.

    A batch is read from its own batchfile and not from the genome_dirs file, so
    that a batch is self-contained: the machine that processes it needs to agree
    with the machine that planned it about which genomes it holds, and the
    batchfile is that agreement written down.

    Parameters
    ----------
    batchfile : str
        Batchfile of one batch.

    @return: (FASTA path, accession) for each genome of the batch.
    """

    rows = []
    with open_text(batchfile) as handle:
        for line in handle:
            line = line.rstrip('\n')
            if not line:
                continue
            fasta, _, accession = line.partition('\t')
            rows.append((fasta, accession))

    return rows


def batchfile_path(batch_dir: str, layout: 'BatchLayout') -> str:
    """Where a batch's plan is, whichever version of this command wrote it.

    The plan is gzipped now and was not before, and a finished output directory
    is still read: the command is run again over one to pick up the work that has
    since been added to it. Asked of a batch that has neither, the answer is
    where the plan would be written, so that a caller's error names the file it
    was looking for.

    Parameters
    ----------
    batch_dir : str
        Batch directory.

    @return: path of the batchfile.
    """

    for name in layout.batchfiles:
        path = os.path.join(batch_dir, name)
        if os.path.exists(path):
            return path

    return os.path.join(batch_dir, layout.batchfiles[0])


def read_accessions(path: str) -> List[str]:
    """Read a file of one accession per line.

    Parameters
    ----------
    path : str
        File to read, which may not exist.

    @return: the accessions, and an empty list where there is no file.
    """

    try:
        with open(path) as handle:
            return [line.strip() for line in handle if line.strip()]
    except OSError:
        return []


def batch_dir_names(out_dir: str, layout: 'BatchLayout') -> List[str]:
    """The batch directories already planned under an output directory.

    Parameters
    ----------
    out_dir : str
        Output directory of the run.

    @return: paths of the batch directories holding a batchfile, in batch order.
    """

    if not os.path.isdir(out_dir):
        return []

    found = []
    for name in sorted(os.listdir(out_dir)):
        path = os.path.join(out_dir, name)
        if name.startswith(BATCH_DIR_PREFIX) and os.path.isdir(path):
            if os.path.exists(batchfile_path(path, layout)):
                found.append(path)

    return found


def create_batches(rows: Sequence[Tuple[str, str]],
                   batch_size: int,
                   out_dir: str,
                   layout: 'BatchLayout') -> List[str]:
    """Cut the release into batches and write the batchfile of each.

    Every batchfile is written before any batch is processed, so that the plan is
    complete the moment the first genome is worked on and a second machine
    starting later finds the same batches rather than making its own.

    Parameters
    ----------
    rows : sequence of tuple
        (FASTA path, accession) for every genome of the release, in the order the
        batches are to be cut from.
    batch_size : int
        Genomes per batch.
    out_dir : str
        Output directory of the run.

    @return: paths of the batch directories created, in batch order.
    """

    created = []
    for index, start in enumerate(range(0, len(rows), batch_size), start=1):
        batch_dir = os.path.join(out_dir, BATCH_DIR_FORMAT.format(index))
        os.makedirs(batch_dir, exist_ok=True)
        write_batchfile(rows[start:start + batch_size],
                        os.path.join(batch_dir, layout.batchfiles[0]), compress=True)
        created.append(batch_dir)

    return created


def canary_payload(**extra: object) -> str:
    """What a canary file says, beyond the fact that it exists.

    A claim has to name its owner for another machine to know whose it is, and
    for this machine to recognise a claim of its own left behind by a process
    that is no longer running.

    Parameters
    ----------
    extra : dict
        Further fields to record, one per line.

    @return: the contents of the canary file.
    """

    fields = [('host', socket.gethostname()),
              ('pid', os.getpid()),
              ('time', datetime.datetime.now().isoformat(timespec='seconds'))]
    fields += sorted(extra.items())

    return ''.join('{}\t{}\n'.format(key, value) for key, value in fields)


def read_canary(path: str) -> Dict[str, str]:
    """Read a canary file back.

    Parameters
    ----------
    path : str
        Canary file.

    @return: its fields, empty where the file has gone or cannot be read.
    """

    fields = {}
    try:
        with open(path) as handle:
            for line in handle:
                key, _, value = line.rstrip('\n').partition('\t')
                fields[key] = value
    except OSError:
        pass

    return fields


def process_alive(pid: str) -> bool:
    """Whether a process of this host is still running.

    Only ever asked about a PID this host recorded. A PID from another host says
    nothing here and is not looked up: PIDs are reused, and taking a batch from a
    machine that is still working on it costs more than leaving it.

    Parameters
    ----------
    pid : str
        Process ID, as the canary recorded it.

    @return: True if the process exists or cannot be ruled out, False if it is
             certainly gone.
    """

    try:
        os.kill(int(pid), 0)
    except ProcessLookupError:
        return False
    except (ValueError, TypeError):
        # not a PID at all, so nothing can be concluded from it
        return True
    except PermissionError:
        # running as another user, which means running
        return True

    return True


def server_time(directory: str) -> float:
    """What time it is by the clock of the file server holding a directory.

    A lease is only as good as the clock it is measured against, and the machines
    sharing an --out_dir have a clock each. What they do share is the file server,
    so a file is created in the directory and the time the server gives it is
    taken as now. Skew between the machines then cannot expire a live claim or
    hold a dead one.

    Parameters
    ----------
    directory : str
        Directory to ask about, which is the batch directory holding the claim.

    @return: the server's idea of now, as a POSIX timestamp; this machine's own
             clock where the directory cannot be written to.
    """

    handle, temp = None, None
    try:
        handle, temp = tempfile.mkstemp(prefix='.now.', dir=directory)
        return os.fstat(handle).st_mtime
    except OSError:
        return time.time()
    finally:
        if handle is not None:
            try:
                os.close(handle)
            except OSError:
                pass
        if temp is not None:
            try:
                os.unlink(temp)
            except OSError:
                pass


def claim_age(running_file: str) -> Optional[float]:
    """How long it is since the machine holding a batch last said so.

    Parameters
    ----------
    running_file : str
        The RUNNING canary of a batch.

    @return: seconds since the claim was last touched, or None if it has gone.
    """

    try:
        touched = os.stat(running_file).st_mtime
    except OSError:
        return None

    return max(0.0, server_time(os.path.dirname(running_file)) - touched)


class Heartbeat(object):
    """Touch a claim while its batch is worked on, so that it does not expire.

    The beating stops when the process holding the batch stops, whether it
    returns, is killed or wedges, which is the whole point: a claim outlives the
    process that made it by one lease and no longer.
    """

    def __init__(self, running_file: str, interval: float = HEARTBEAT_SECONDS) -> None:
        """Initialization.

        Parameters
        ----------
        running_file : str
            The RUNNING canary to keep alive.
        interval : float
            Seconds between touches.

        @return: None
        """

        self.running_file = running_file
        self.interval = interval
        self.stop = threading.Event()
        self.thread = None

    def beat(self) -> None:
        """Touch the claim until asked to stop.

        @return: None
        """

        # wait() returns True only when it was set, so the loop ends the moment
        # the batch does rather than after one more interval
        while not self.stop.wait(self.interval):
            try:
                os.utime(self.running_file, None)
            except OSError:
                # the claim has gone, which another machine taking the batch or
                # the batch finishing both look like; there is nothing to keep
                return

    def __enter__(self) -> 'Heartbeat':
        self.thread = threading.Thread(target=self.beat, daemon=True)
        self.thread.start()
        return self

    def __exit__(self, *exc_info: object) -> None:
        self.stop.set()
        if self.thread is not None:
            self.thread.join(timeout=self.interval)




@contextlib.contextmanager
def batch_log(batch_dir: str, logger: logging.Logger,
              layout: 'BatchLayout') -> Iterator[None]:
    """Write what happens to a batch into the batch's own directory as well.

    Parameters
    ----------
    batch_dir : str
        Batch directory, which takes the command's own log.
    logger : logging.Logger
        The logger to tee, which is the 'timestamp' logger of the run.
    layout : BatchLayout
        Names of the command whose batches these are.

    @return: a context in which the logger also writes to the batch.
    """

    handler = logging.FileHandler(os.path.join(batch_dir, layout.log), 'a')
    handler.setFormatter(logging.Formatter(
        fmt='[%(asctime)s] %(levelname)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'))
    logger.addHandler(handler)
    try:
        yield
    finally:
        logger.removeHandler(handler)
        handler.close()


def age_phrase(age: Optional[float]) -> str:
    """How long ago something was, as a log line says it.

    Parameters
    ----------
    age : float
        Seconds ago, or None where there is nothing to say.

    @return: a phrase naming the time, for a log message.
    """

    if age is None:
        return 'never'
    if age < 90:
        return '{:.0f}s ago'.format(age)
    if age < 5400:
        return '{:.0f}m ago'.format(age / 60)

    return '{:.1f}h ago'.format(age / 3600)


def batch_state(batch_dir: str) -> str:
    """What has happened to a batch.

    Parameters
    ----------
    batch_dir : str
        Batch directory.

    @return: one of the STATE_* constants.
    """

    if os.path.exists(os.path.join(batch_dir, SUCCESS_CANARY)):
        return STATE_SUCCESS
    if os.path.exists(os.path.join(batch_dir, RUNNING_CANARY)):
        return STATE_RUNNING
    if os.path.exists(os.path.join(batch_dir, FAILED_CANARY)):
        return STATE_FAILED

    return STATE_PENDING


def stale_claim(running_file: str,
                lease: float = CLAIM_LEASE_SECONDS) -> bool:
    """Whether a claim has been given up by the machine that made it.

    Two things say so. A claim of this host's whose process has gone is dead and
    known to be dead, which is what a reset leaves behind on the machine that
    reset. Any claim that has not been touched for a lease is dead as well: the
    machine holding a batch says so every HEARTBEAT_SECONDS for as long as it
    works, so silence for far longer than that is a machine that stopped, and
    whether it stopped by dying, by being killed or by wedging on a mount is
    neither knowable from here nor worth knowing.

    The second rule is what makes a batch recoverable from ANOTHER machine, and
    it is also the more reliable of the two: PIDs are reused, so after a reset
    the PID a claim names is as likely to belong to something new as to be
    missing, and a liveness check then holds a dead claim forever.

    Parameters
    ----------
    running_file : str
        The RUNNING canary of a batch.
    lease : float
        Seconds a claim survives without being touched.

    @return: True if the claim can be taken over without being asked to.
    """

    fields = read_canary(running_file)
    if (fields.get('host') == socket.gethostname()
            and not process_alive(fields.get('pid', ''))):
        return True

    age = claim_age(running_file)

    return age is not None and age > lease


def keep_failure_record(batch_dir: str) -> None:
    """Move a previous attempt's FAILED aside instead of deleting it.

    A batch is retried by claiming it, and the claim has to clear FAILED or the
    batch would still read as failed while it runs. Deleting it takes with it the
    only record of why the batch failed, which on a batch that fails the same way
    every time is the thing a person needs to read. It is kept under the time it
    was cleared, beside the batch it belongs to.

    Parameters
    ----------
    batch_dir : str
        Batch directory.

    @return: None
    """

    failed = os.path.join(batch_dir, FAILED_CANARY)
    kept = '{}.{}'.format(failed, datetime.datetime.now().strftime('%Y%m%dT%H%M%S'))
    try:
        os.rename(failed, kept)
    except OSError:
        pass


def release_claim(batch_dir: str) -> None:
    """Give up a claim without saying anything about how the batch went.

    What an interrupted run leaves: the batch was neither finished nor tried and
    found wanting, and the machine that held it is about to stop. Releasing it
    has the next run take it up rather than wait out the lease.

    Parameters
    ----------
    batch_dir : str
        Batch directory.

    @return: None
    """

    try:
        os.unlink(os.path.join(batch_dir, RUNNING_CANARY))
    except OSError:
        pass


def claim_batch(batch_dir: str,
                reclaim: bool = False,
                lease: float = CLAIM_LEASE_SECONDS) -> bool:
    """Take a batch for this machine, if no other machine holds it.

    The claim is made by linking a uniquely named file onto RUNNING rather than by
    creating RUNNING directly: this runs against NFS, where O_EXCL has never been
    the operation that two machines can race on safely and link() is. link()
    failing is read back rather than believed, for the same reason -- an NFS
    client can be told a link failed that in fact succeeded, and the link count of
    the file it made is what settles it.

    Parameters
    ----------
    batch_dir : str
        Batch directory to claim.
    reclaim : bool
        Take a batch another machine holds before its claim has expired. Only
        ever right when that machine is known not to be working on it.
    lease : float
        Seconds a claim survives without being touched.

    @return: True if this machine now holds the batch.
    """

    running_file = os.path.join(batch_dir, RUNNING_CANARY)

    if os.path.exists(running_file):
        if not (reclaim or stale_claim(running_file, lease)):
            return False
        try:
            os.unlink(running_file)
        except OSError:
            return False

    tmp_file = os.path.join(batch_dir, '.{}.{}.{}'.format(
        RUNNING_CANARY, os.getpid(), uuid.uuid4().hex))
    with open(tmp_file, 'w') as handle:
        handle.write(canary_payload())

    try:
        os.link(tmp_file, running_file)
        claimed = True
    except OSError:
        # the link may have been made even so, and the link count says whether
        claimed = os.stat(tmp_file).st_nlink == 2
    finally:
        try:
            os.unlink(tmp_file)
        except OSError:
            pass

    # a previous attempt on this batch is no longer what happened to it, though
    # what it had to say about itself is kept
    if claimed:
        keep_failure_record(batch_dir)

    return claimed


def finish_batch(batch_dir: str, **extra: object) -> None:
    """Record that a batch finished, and give up the claim on it.

    SUCCESS is written before RUNNING is removed. The other order leaves a moment
    in which the batch looks unclaimed and unfinished, which is the one state that
    would have a second machine repeat it.

    Parameters
    ----------
    batch_dir : str
        Batch directory.
    extra : dict
        Further fields to record in the canary.

    @return: None
    """

    with open(os.path.join(batch_dir, SUCCESS_CANARY), 'w') as handle:
        handle.write(canary_payload(**extra))

    try:
        os.unlink(os.path.join(batch_dir, RUNNING_CANARY))
    except OSError:
        pass


def fail_batch(batch_dir: str, reason: str) -> None:
    """Record that a batch was attempted and did not finish.

    The claim is given up, so the batch is retried by the next run without
    --reclaim: this machine is known not to be working on it, which is exactly
    what --reclaim exists to assert about a machine that cannot be asked.

    Parameters
    ----------
    batch_dir : str
        Batch directory.
    reason : str
        What went wrong, recorded for whoever reads the directory.

    @return: None
    """

    with open(os.path.join(batch_dir, FAILED_CANARY), 'w') as handle:
        handle.write(canary_payload(reason=' '.join(str(reason).split())))

    try:
        os.unlink(os.path.join(batch_dir, RUNNING_CANARY))
    except OSError:
        pass


def genome_fasta_of(accession: str, genome_dir: str) -> str:
    """The file a batch is planned around by default: the genomic FASTA.

    What a batch is planned around is the file the work reads, which differs by
    command: trans_table and prodigal read the genomic FASTA NCBI served, while
    hmmsearch reads the proteins prodigal called. It is the first column of the
    batchfile and it is what split_by_fasta() stats, so a command that works from
    another file says so with its own function rather than being given a
    batchfile naming files it will not open.

    Parameters
    ----------
    accession : str
        Accession of the genome, unused here and taken for the signature.
    genome_dir : str
        Genome directory of the release.

    @return: path of the genomic FASTA in that directory.
    """

    return genomic_fasta(genome_dir)


def plan_batches(gtdb_genome_path_file: str,
             out_dir: str,
             batch_size: int,
             layout: BatchLayout,
             logger: logging.Logger,
             genome_file: Callable[[str, str], str] = genome_fasta_of) -> List[str]:
    """Settle which genomes are in which batch, once for every machine.

    A plan already under the output directory is used as it stands. It is what
    another machine is working from and what the finished batches were cut
    from, and partitioning a release again that has since gained or lost a
    genome would move genomes between batches that are already done.

    Parameters
    ----------
    gtdb_genome_path_file : str
        genome_dirs file of the release.
    out_dir : str
        Output directory of the run.
    batch_size : int
        Genomes per batch.
    layout : BatchLayout
        What the command calls the files of its own batches.
    logger : logging.Logger
        Where the plan is reported.
    genome_file : callable
        (accession, genome directory) to the file the work reads, which becomes
        the first column of the batchfile. The default is the genomic FASTA.

    @return: paths of the batch directories, in batch order.
    """

    existing = batch_dir_names(out_dir, layout)
    if existing:
        logger.warning(
            'warning: {:,} batch(es) are already planned under {}; using them '
            'and not regenerating the batchfiles. Remove the batch directories '
            'to partition the release again.'.format(len(existing), out_dir))
        return existing

    genomes = read_genome_dirs(gtdb_genome_path_file)
    logger.info('Read {:,} genomes from {}.'.format(
        len(genomes), gtdb_genome_path_file))

    # sorted so that which genomes are in batch N follows from the set of
    # genomes, and not from the order the genome_dirs file was written in
    genomes.sort(key=lambda genome: genome[0])

    # named, not stat-ed: whether the file is there is asked of each batch as
    # it is run, where ten thousand stat calls are nothing against the hours
    # the work then takes, rather than of the whole release here, where on
    # r237 it is hours over NFS before a single batch directory exists and a
    # run stopped in it has no plan to resume from. It also leaves which
    # genomes share a batch following from the genome_dirs file alone.
    rows = [(genome_file(accession, genome_dir), accession)
            for accession, genome_dir in genomes]
    if not rows:
        raise RuntimeError(
            '{} names no genomes.'.format(gtdb_genome_path_file))

    batches = create_batches(rows, batch_size, out_dir, layout)
    logger.info('Planned {:,} genomes as {:,} batch(es) of up to {:,}.'.format(
        len(rows), len(batches), batch_size))

    return batches


def write_table(rows: Sequence[Sequence[str]], path: str,
                header: Sequence[str],
                compress: bool = False) -> None:
    """Write a headered TSV, of a batch or of a release.

    The file is written whether or not there are any rows: a batch that finished
    with nothing to report says so with a header and no rows, and the release
    file is then the concatenation of every batch's, however many they found.

    Parameters
    ----------
    rows : sequence of sequence of str
        Rows, as comparison_rows(), conflicts_from_comparison() or
        annotate_conflicts() returned them.
    path : str
        File to write.
    header : sequence of str
        Column names: COMPARISON_HEADER, CONFLICT_HEADER for a batch, or
        CONFLICT_HEADER_CHECKM2 for the release, which carries the CheckM2
        columns as well.
    compress : bool
        Write it gzipped, which the comparison is and the conflicts are not.

    @return: None
    """

    with (gzip.open(path, 'wt') if compress else open(path, 'w')) as handle:
        handle.write('\t'.join(header) + '\n')
        for row in rows:
            handle.write('\t'.join(row) + '\n')


def concatenate(files: Sequence[str], path: str, compress: bool = False) -> int:
    """Join the tables of every batch into one, keeping a single header.

    Parameters
    ----------
    files : sequence of str
        Files to join, each with the same header, in batch order.
    path : str
        File to write.
    compress : bool
        Write it gzipped, which the release summary is and the conflicts are not:
        the summary is a row per genome of the release and the conflicts are a few
        hundred rows meant to be looked at.

    @return: number of rows written, the header not counted.
    """

    written = 0
    with (gzip.open(path, 'wt') if compress else open(path, 'w')) as out:
        for index, name in enumerate(files):
            # open_text: a batch's comparison is gzipped and its conflicts are
            # not, and gTranslate's summary is whatever gTranslate wrote
            with open_text(name) as handle:
                header = handle.readline()
                if index == 0:
                    out.write(header)
                for line in handle:
                    if line.strip():
                        out.write(line)
                        written += 1

    return written


def tally_reasons(rows: Sequence[Tuple[str, str]]) -> Dict[str, int]:
    """How many genomes there are of each reason.

    Parameters
    ----------
    rows : sequence of tuple
        (accession, reason), as no_prediction_rows() returned them.

    @return: reason to the number of genomes with it.
    """

    counts = {}
    for _, reason in rows:
        counts[reason] = counts.get(reason, 0) + 1

    return counts
