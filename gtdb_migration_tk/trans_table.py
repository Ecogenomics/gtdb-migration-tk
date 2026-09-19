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
trans_table.py -- predict the translation table of each genome with gTranslate.

gTranslate is run as a SUBPROCESS rather than imported. It requires Python >= 3.12
where this toolkit supports 3.8, and it is installed in an environment of its own
holding pinned versions of scikit-learn, xgboost and lightgbm that the predictions
depend on; importing it would make those the toolkit's own dependencies and its
floor the toolkit's floor. What crosses between the two is a file of genomes and a
file of predictions, so a subprocess is all the coupling the work needs.

BATCHES
The release is cut into batches of --batch_size genomes, each with a directory of
its own under --out_dir, because a release takes days to predict and nothing that
takes days should have to be started again from the beginning. A batch is the unit
of work, of restart and of sharing between machines: several machines may be given
the same --out_dir and will divide the release between them without being told
which part to take, and a machine that resets loses the batch it was on rather
than the run.

The batches are settled BEFORE any of them is processed, and are the plan every
machine works from. The genomes are sorted by accession first, so which genomes
are in batch N follows from the set of genomes alone and not from the order a
genome_dirs file happens to be written in. Once the batchfiles exist they are
authoritative: a later run reuses them and says so rather than partitioning the
release again, since a second partition of a release that has gained a genome
would move genomes between batches that are already finished.

The plan is cut from the genome_dirs file and nothing else. Whether a genome's
FASTA is actually on disk is asked of each batch as that batch is run, not of
the release beforehand: gTranslate checks the paths of a batchfile itself and
refuses the whole batch if one is missing, so the check has to happen, but asked
of 1.35M genomes over NFS it is hours in which nothing is written and nothing
can be resumed, where asked of ten thousand it is a minute against a batch that
then runs for hours. A genome left out is named in the batch's own directory.

CANARIES
A batch directory carries its state in files, because a file is what two machines
can both see and what survives the process that wrote it:

    RUNNING   a machine is working on this batch, and names itself, its PID and
              when it started
    SUCCESS   the batch finished and its results are complete
    FAILED    the batch was attempted and gTranslate returned non-zero

RUNNING is created by linking a uniquely named file to it, not by opening it with
O_EXCL: the release lives on NFS, where link() is the operation that is atomic
across machines. A batch already claimed is left to the machine that claimed it,
and reported at the end rather than waited for. The exception is a claim this host
made in a process that no longer exists, which is what a machine reset leaves
behind and is reclaimed automatically; a claim from another host cannot be
checked from here and is taken only when --reclaim says to take it. A batch that
FAILED is retried by the next run without a flag, the machine that failed it
having cleared its own claim.

THE COMPARISON
gTranslate predicts a table; NCBI declares one in the GFF it serves for a genome
that it has annotated. Each batch writes ncbi_tt_comparison.tsv over the genomes
where both are known, saying for each whether the two agree, and carrying the
coding densities the prediction was made from and the genome's NCBI taxonomy.
What the comparison is for is the disagreements: a genome whose genes GTDB would
call under a table NCBI does not agree with, and whether those genomes fall
together in the taxonomy. The table NCBI declares is read by
ncbi_utils.ncbi_translation_table(), which prodigal reads it with too.

Once every batch has succeeded the run also writes the comparison and the
prediction summary for the whole release at the top of --out_dir, so that a
release finished across several machines is one file to read.
"""

import csv
import datetime
import logging
import os
import shutil
import socket
import subprocess
import uuid
from concurrent.futures import ThreadPoolExecutor
from typing import Dict, List, Optional, Sequence, Tuple

from tqdm import tqdm

from gtdb_migration_tk.biolib_lite.external.execute import check_dependencies
from gtdb_migration_tk.biolib_lite.common import canonical_gid
from gtdb_migration_tk.ncbi_utils import (GENOMIC_FASTA_EXT, NCBI_NA,
                                          genomic_gff, ncbi_translation_table)
from gtdb_migration_tk.utils.common import (TT_SUMMARY_DENSITY_4,
                                            TT_SUMMARY_DENSITY_11,
                                            TT_SUMMARY_TABLE,
                                            checkm_translation_table,
                                            read_translation_table_summary)


# The executable, looked up on PATH rather than given a path of its own: it is
# installed as a module here, and the module is what puts it on PATH along with
# prodigal and GTRANSLATE_MODEL_PATH, none of which this command can supply.
GTRANSLATE_BIN = 'gtranslate'

# The subcommand run. gTranslate's other subcommands either train models or plot
# what one predicted, neither of which belongs in a migration.
DETECT_TABLE = 'detect_table'

# Genomes per batch. Large enough that the cost of starting gTranslate and loading
# its classifiers is nothing against the genomes it then processes, small enough
# that a machine lost mid-batch costs hours and not days.
DEFAULT_BATCH_SIZE = 10000

# Batch directories are numbered rather than named for the genomes they hold: the
# accessions of a batch are in its batchfile, and a name is a thing to sort by.
BATCH_DIR_PREFIX = 'batch_'
BATCH_DIR_FORMAT = BATCH_DIR_PREFIX + '{:06d}'

# Written into the batch directory rather than a temporary one: it is the record
# of which genomes the batch is, and it is what a later run and another machine
# read to agree on that without partitioning the release again.
BATCHFILE_NAME = 'gtranslate_batchfile.tsv'

# What gTranslate is handed when some genome of the batch has no FASTA to process,
# and the accessions left out of it. Written only in that case, so the file being
# there at all says a batch had something wrong with it. BATCHFILE_NAME stays the
# record of which genomes the batch IS, which is what the comparison reads.
PRESENT_BATCHFILE_NAME = 'gtranslate_batchfile_present.tsv'
MISSING_NAME = 'missing_genomic_fasta.tsv'

# The state of a batch, held in files so that another machine can see it and so
# that it outlives the process that wrote it.
RUNNING_CANARY = 'RUNNING'
SUCCESS_CANARY = 'SUCCESS'
FAILED_CANARY = 'FAILED'

STATE_PENDING = 'pending'
STATE_RUNNING = 'running'
STATE_SUCCESS = 'success'
STATE_FAILED = 'failed'

# gTranslate names its summary for its --prefix, which this command passes through.
DEFAULT_PREFIX = 'gtranslate'
SUMMARY_SUFFIX = '.translation_table_summary.tsv'

COMPARISON_NAME = 'ncbi_tt_comparison.tsv'

# checkm_tt sits beside the other two tables rather than at the end: the three
# are the answers to one question, and result reports on the two that are
# predictions of what the genome uses. It is the table the coding density rule
# alone would choose, which is what Prodigal and CheckM do unaided, and it cannot
# express table 25 at all -- a genome gTranslate calls 25 is one the old rule was
# never able to get right.
COMPARISON_HEADER = ('genome_id', 'gtranslate_tt', 'ncbi_tt', 'checkm_tt',
                     'result', 'coding_density_4', 'coding_density_11',
                     'ncbi_taxonomy')

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

# The two outcomes a comparison has. A genome NCBI declares no table for is not
# one of them: it is left out of the file, there being nothing to compare.
AGREE = 'agree'
CONFLICT = 'conflict'


def genomic_fasta(genome_dir: str) -> str:
    """The genomic FASTA NCBI serves for a genome.

    NCBI names the file for the assembly and names the genome directory the same,
    so the file is named rather than searched for -- _cds_from_genomic.fna.gz and
    _rna_from_genomic.fna.gz end the same way and are different files.

    Parameters
    ----------
    genome_dir : str
        Genome directory, of a release or of the mirror.

    @return: path of the genomic FASTA in that directory, which may not exist.
    """

    assembly = os.path.basename(os.path.normpath(genome_dir))

    return os.path.join(genome_dir, assembly + GENOMIC_FASTA_EXT)


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


def write_batchfile(rows: Sequence[Tuple[str, str]], batchfile: str) -> None:
    """Write the two-column batchfile gTranslate reads.

    The genome ID given is the accession the genome_dirs file names, so every row
    of the prediction table can be matched back to the genome directory it came
    from without canonicalising anything.

    Parameters
    ----------
    rows : sequence of tuple
        (FASTA path, accession) for each genome to ask about.
    batchfile : str
        File to write.

    @return: None
    """

    with open(batchfile, 'w') as handle:
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
    with open(batchfile) as handle:
        for line in handle:
            line = line.rstrip('\n')
            if not line:
                continue
            fasta, _, accession = line.partition('\t')
            rows.append((fasta, accession))

    return rows


def check_batch_fastas(batch_dir: str,
                       threads: int = STAT_THREADS
                       ) -> Tuple[str, List[Tuple[str, str]], List[str]]:
    """Check over the genomes of one batch, just before gTranslate is run on it.

    The check belongs to the batch rather than to the plan. Asked of the whole
    release up front it is hours of stat calls before a single batch directory
    exists -- on r237, 1.35M genomes over NFS -- and a run stopped in it has
    nothing to resume from, the plan not yet written. Asked of a batch it is ten
    thousand stat calls against a batch that then runs for hours, it happens
    while other machines are already working, and a machine lost during it costs
    one batch. It also leaves the batch boundaries following from the genome_dirs
    file alone, rather than from what stat said on the day the plan was cut.

    Nothing is written where every genome has its FASTA, which is the normal case
    and the one where BATCHFILE_NAME is handed straight to gTranslate.

    Parameters
    ----------
    batch_dir : str
        Batch directory, holding the batchfile the plan cut.
    threads : int
        Stat calls to keep in flight at once.

    @return: (batchfile, present, missing) -- the file to hand gTranslate, the
             rows it names, and the accessions left out of it.
    """

    batchfile = os.path.join(batch_dir, BATCHFILE_NAME)
    present, missing = split_by_fasta(read_batchfile(batchfile), threads)

    if not missing:
        return batchfile, present, missing

    write_batchfile(present, os.path.join(batch_dir, PRESENT_BATCHFILE_NAME))
    with open(os.path.join(batch_dir, MISSING_NAME), 'w') as handle:
        for accession in missing:
            handle.write('{}\n'.format(accession))

    return os.path.join(batch_dir, PRESENT_BATCHFILE_NAME), present, missing


def summary_name(prefix: Optional[str] = None) -> str:
    """The name gTranslate gives its translation table summary.

    Parameters
    ----------
    prefix : str
        The --prefix the run passes gTranslate, or None for gTranslate's default.

    @return: filename of the summary within a batch directory.
    """

    return (prefix or DEFAULT_PREFIX) + SUMMARY_SUFFIX


def batch_dir_names(out_dir: str) -> List[str]:
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
            if os.path.exists(os.path.join(path, BATCHFILE_NAME)):
                found.append(path)

    return found


def create_batches(rows: Sequence[Tuple[str, str]],
                   batch_size: int,
                   out_dir: str) -> List[str]:
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
                        os.path.join(batch_dir, BATCHFILE_NAME))
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


def stale_claim(running_file: str) -> bool:
    """Whether a claim was left behind by a process of this host that has gone.

    This is what a machine reset leaves: a RUNNING file naming a PID on this host
    that no longer exists. A claim from another host is never stale here, however
    old, because nothing on this machine can tell a reset host from a busy one.

    Parameters
    ----------
    running_file : str
        The RUNNING canary of a batch.

    @return: True if the claim can be taken over without being asked to.
    """

    fields = read_canary(running_file)
    if not fields:
        return False

    return (fields.get('host') == socket.gethostname()
            and not process_alive(fields.get('pid', '')))


def claim_batch(batch_dir: str, reclaim: bool = False) -> bool:
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
        Take a batch another host has claimed. Only ever right when that host is
        known not to be working on it.

    @return: True if this machine now holds the batch.
    """

    running_file = os.path.join(batch_dir, RUNNING_CANARY)

    if os.path.exists(running_file):
        if not (reclaim or stale_claim(running_file)):
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

    # a previous attempt on this batch is no longer what happened to it
    if claimed:
        try:
            os.unlink(os.path.join(batch_dir, FAILED_CANARY))
        except OSError:
            pass

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


def detect_table_command(batchfile: str,
                         out_dir: str,
                         cpus: int = 1,
                         tmp_dir: Optional[str] = None,
                         force: bool = False,
                         keep_called_genes: bool = False,
                         prefix: Optional[str] = None,
                         custom_model_path: Optional[str] = None) -> List[str]:
    """Build the gtranslate detect_table command line.

    An option the run was not given is left off the command line rather than
    passed with a default of this module's choosing, so gTranslate's own defaults
    remain the defaults and do not have to be tracked here as it changes.

    Parameters
    ----------
    batchfile : str
        Batchfile naming the genomes to process.
    out_dir : str
        Directory gTranslate writes its results to, which is the batch directory.
    cpus : int
        Number of genomes processed at once.
    tmp_dir : str
        Directory for gTranslate's intermediate files, or None for its default.
    force : bool
        Carry on when a single genome fails rather than stopping the batch.
    keep_called_genes : bool
        Keep the genes called under the predicted table.
    prefix : str
        Prefix of gTranslate's output files, or None for its default.
    custom_model_path : str
        Classifiers to predict with, or None to use GTRANSLATE_MODEL_PATH.

    @return: the command as a list of arguments, ready for subprocess.
    """

    cmd = [GTRANSLATE_BIN, DETECT_TABLE,
           '--batchfile', batchfile,
           '--out_dir', out_dir,
           '--cpus', str(cpus)]

    if tmp_dir:
        cmd += ['--tmpdir', tmp_dir]
    if prefix:
        cmd += ['--prefix', prefix]
    if custom_model_path:
        cmd += ['--custom_model_path', custom_model_path]
    if force:
        cmd += ['--force']
    if keep_called_genes:
        cmd += ['--keep_called_genes']

    return cmd


def read_taxonomy(taxonomy_file: str) -> Dict[str, str]:
    """Read the standardised NCBI taxonomy of the genomes.

    The file is the two-column accession / semicolon-separated lineage TSV
    ncbi_metadata_sync writes. Each genome is recorded under the accession as
    given AND under its canonical form, so that a GenBank genome of a release
    finds the lineage the taxonomy holds against its RefSeq counterpart; an
    accession given exactly is preferred to a canonical match.

    Parameters
    ----------
    taxonomy_file : str
        Standardised NCBI taxonomy file.

    @return: accession, and canonical accession, to lineage.
    """

    exact, canonical = {}, {}
    with open(taxonomy_file) as handle:
        for line in tqdm(handle, ncols=100, leave=False, desc='Reading taxonomy'):
            line = line.rstrip('\n')
            if not line:
                continue
            accession, _, lineage = line.partition('\t')
            if not lineage:
                continue
            exact[accession] = lineage
            canonical.setdefault(canonical_gid(accession), lineage)

    canonical.update(exact)

    return canonical


def lineage_of(accession: str, taxonomy: Dict[str, str]) -> str:
    """The NCBI lineage of a genome.

    Parameters
    ----------
    accession : str
        Accession of the genome.
    taxonomy : dict
        Taxonomy as read_taxonomy() returned it.

    @return: the lineage, or NCBI's null where the taxonomy does not hold one.
    """

    if accession in taxonomy:
        return taxonomy[accession]

    return taxonomy.get(canonical_gid(accession), NCBI_NA)


def comparison_rows(predictions: Dict[str, Dict[str, str]],
                    genome_dirs: Dict[str, str],
                    taxonomy: Dict[str, str]) -> Tuple[List[Tuple[str, ...]], int]:
    """Compare what gTranslate predicted against what NCBI declares.

    Only a genome NCBI declares a table for is reported: a genome NCBI has not
    annotated has nothing to compare against, and a row saying so would be a row
    per unannotated genome of the release saying nothing.

    Parameters
    ----------
    predictions : dict
        Predictions as read_predictions() returned them.
    genome_dirs : dict
        Accession to genome directory, for the genomes of the batch.
    taxonomy : dict
        Taxonomy as read_taxonomy() returned it.

    @return: (rows, no_ncbi_table), the comparison rows in accession order and
             the number of genomes NCBI declared no table for.
    """

    rows, no_ncbi_table = [], 0
    for accession in sorted(predictions):
        genome_dir = genome_dirs.get(accession)
        if genome_dir is None:
            continue

        ncbi_table = ncbi_translation_table(genomic_gff(genome_dir))
        if ncbi_table is None:
            no_ncbi_table += 1
            continue

        predicted = predictions[accession]
        table = predicted.get(TT_SUMMARY_TABLE, '')
        density_4 = predicted.get(TT_SUMMARY_DENSITY_4, '')
        density_11 = predicted.get(TT_SUMMARY_DENSITY_11, '')
        result = AGREE if table.strip() == str(ncbi_table) else CONFLICT
        checkm_table = checkm_translation_table(density_4, density_11)

        rows.append((accession,
                     table,
                     str(ncbi_table),
                     str(checkm_table) if checkm_table else NCBI_NA,
                     result,
                     density_4,
                     density_11,
                     lineage_of(accession, taxonomy)))

    return rows, no_ncbi_table


def write_comparison(rows: Sequence[Tuple[str, ...]], path: str) -> None:
    """Write the comparison of a batch, or of the release.

    Parameters
    ----------
    rows : sequence of tuple
        Comparison rows, as comparison_rows() returned them.
    path : str
        File to write.

    @return: None
    """

    with open(path, 'w') as handle:
        handle.write('\t'.join(COMPARISON_HEADER) + '\n')
        for row in rows:
            handle.write('\t'.join(row) + '\n')


def concatenate(files: Sequence[str], path: str) -> int:
    """Join the tables of every batch into one, keeping a single header.

    Parameters
    ----------
    files : sequence of str
        Files to join, each with the same header, in batch order.
    path : str
        File to write.

    @return: number of rows written, the header not counted.
    """

    written = 0
    with open(path, 'w') as out:
        for index, name in enumerate(files):
            with open(name) as handle:
                header = handle.readline()
                if index == 0:
                    out.write(header)
                for line in handle:
                    if line.strip():
                        out.write(line)
                        written += 1

    return written


class GTranslate(object):
    """Predict the translation table of each genome of a release, in batches."""

    def __init__(self,
                 cpus: int = 1,
                 batch_size: int = DEFAULT_BATCH_SIZE,
                 tmp_dir: Optional[str] = None,
                 force: bool = False,
                 keep_called_genes: bool = False,
                 prefix: Optional[str] = None,
                 custom_model_path: Optional[str] = None,
                 reclaim: bool = False) -> None:
        """Initialization.

        Parameters
        ----------
        cpus : int
            Number of genomes gTranslate processes at once.
        batch_size : int
            Genomes per batch.
        tmp_dir : str
            Directory for gTranslate's intermediate files, or None for its default.
        force : bool
            Carry on when a single genome fails rather than stopping the batch.
        keep_called_genes : bool
            Keep the genes called under the predicted table.
        prefix : str
            Prefix of gTranslate's output files, or None for its default.
        custom_model_path : str
            Classifiers to predict with, or None to use GTRANSLATE_MODEL_PATH.
        reclaim : bool
            Take over batches another machine claimed and did not finish.

        @return: None
        """

        self.cpus = cpus
        self.batch_size = batch_size
        self.tmp_dir = tmp_dir
        self.force = force
        self.keep_called_genes = keep_called_genes
        self.prefix = prefix
        self.custom_model_path = custom_model_path
        self.reclaim = reclaim

        check_dependencies(['gtranslate', 'prodigal'])

        self.logger = logging.getLogger('timestamp')

    def plan_batches(self, gtdb_genome_path_file: str, out_dir: str) -> List[str]:
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

        @return: paths of the batch directories, in batch order.
        """

        existing = batch_dir_names(out_dir)
        if existing:
            self.logger.warning(
                'warning: {:,} batch(es) are already planned under {}; using them '
                'and not regenerating the batchfiles. Remove the batch directories '
                'to partition the release again.'.format(len(existing), out_dir))
            return existing

        genomes = read_genome_dirs(gtdb_genome_path_file)
        self.logger.info('Read {:,} genomes from {}.'.format(
            len(genomes), gtdb_genome_path_file))

        # sorted so that which genomes are in batch N follows from the set of
        # genomes, and not from the order the genome_dirs file was written in
        genomes.sort(key=lambda genome: genome[0])

        # named, not stat-ed: whether the file is there is asked of each batch as
        # it is run, where ten thousand stat calls are nothing against the hours
        # gTranslate then spends, rather than of the whole release here, where on
        # r237 it is hours over NFS before a single batch directory exists and a
        # run stopped in it has no plan to resume from. It also leaves which
        # genomes share a batch following from the genome_dirs file alone.
        rows = [(genomic_fasta(genome_dir), accession)
                for accession, genome_dir in genomes]
        if not rows:
            raise RuntimeError(
                '{} names no genomes.'.format(gtdb_genome_path_file))

        batches = create_batches(rows, self.batch_size, out_dir)
        self.logger.info('Planned {:,} genomes as {:,} batch(es) of up to {:,}.'.format(
            len(rows), len(batches), self.batch_size))

        return batches

    def run_gtranslate(self, batch_dir: str) -> None:
        """Run gTranslate over the genomes of one batch.

        Parameters
        ----------
        batch_dir : str
            Batch directory, which holds the batchfile and takes the results.

        @return: None
        """

        # gTranslate checks the paths of a batchfile before it starts and refuses
        # the whole batch if one of them is missing, so a genome whose FASTA is
        # not there costs the other ten thousand -- and costs them again on every
        # retry, the batch failing identically each time. It is left out here and
        # named in the batch directory instead.
        batchfile, present, missing = check_batch_fastas(batch_dir)
        if missing:
            self.logger.warning(
                'warning: {:,} genome(s) of {} have no genomic FASTA and were left '
                'out; the first is {}. They are named in {}.'.format(
                    len(missing), os.path.basename(batch_dir), missing[0],
                    MISSING_NAME))
        if not present:
            raise RuntimeError(
                'None of the {:,} genomes of {} has a genomic FASTA to process.'.format(
                    len(missing), batch_dir))

        cmd = detect_table_command(batchfile,
                                   batch_dir,
                                   cpus=self.cpus,
                                   tmp_dir=self.tmp_dir,
                                   force=self.force,
                                   keep_called_genes=self.keep_called_genes,
                                   prefix=self.prefix,
                                   custom_model_path=self.custom_model_path)
        self.logger.info('Command: {}'.format(' '.join(cmd)))

        # gTranslate's output is left to the terminal rather than read back a
        # line at a time and logged again. Its progress bars redraw one line with
        # a carriage return, and a carriage return is a line ending to a reader,
        # so re-logging what it prints turns each bar into one line per redraw and
        # the bar into a wall. Nothing is lost by not capturing it: gTranslate
        # writes its own gtranslate.log into the batch directory, and writes it
        # without the bars, which is the better record of a batch anyway.
        # --silent is honoured by discarding the output rather than by showing it.
        silent = getattr(self.logger, 'is_silent', False)
        proc = subprocess.run(cmd,
                              stdout=subprocess.DEVNULL if silent else None,
                              stderr=subprocess.STDOUT if silent else None)

        if proc.returncode != 0:
            raise RuntimeError('{} returned exit code {}.'.format(
                GTRANSLATE_BIN, proc.returncode))

    def compare_batch(self, batch_dir: str, taxonomy: Dict[str, str]) -> int:
        """Compare a batch's predictions against the tables NCBI declares.

        The genome directories are taken from the batch's own batchfile, so the
        comparison asks about the genomes the batch was run on rather than about
        whatever a genome_dirs file says now.

        Parameters
        ----------
        batch_dir : str
            Batch directory.
        taxonomy : dict
            Taxonomy as read_taxonomy() returned it.

        @return: number of genomes compared.
        """

        genome_dirs = {accession: os.path.dirname(fasta) for fasta, accession
                       in read_batchfile(os.path.join(batch_dir, BATCHFILE_NAME))}

        predictions = read_translation_table_summary(
            os.path.join(batch_dir, summary_name(self.prefix)))

        rows, no_ncbi_table = comparison_rows(predictions, genome_dirs, taxonomy)
        write_comparison(rows, os.path.join(batch_dir, COMPARISON_NAME))

        result = COMPARISON_HEADER.index('result')
        conflicts = sum(1 for row in rows if row[result] == CONFLICT)
        self.logger.info(
            'Compared {:,} genomes: {:,} agree with NCBI, {:,} conflict; '
            '{:,} genome(s) have no table from NCBI to compare.'.format(
                len(rows), len(rows) - conflicts, conflicts, no_ncbi_table))

        return len(rows)

    def aggregate(self, batches: Sequence[str], out_dir: str) -> None:
        """Write the comparison and the summary for the whole release.

        Written only once every batch has succeeded, so that the files at the top
        of the output directory are either the whole release or absent, and never
        a part of it that reads like the whole.

        Parameters
        ----------
        batches : sequence of str
            Every batch directory of the run, in batch order.
        out_dir : str
            Output directory of the run.

        @return: None
        """

        unfinished = [batch for batch in batches
                      if batch_state(batch) != STATE_SUCCESS]
        if unfinished:
            self.logger.info(
                '{:,} of {:,} batch(es) are done; the release comparison is '
                'written once they all are.'.format(
                    len(batches) - len(unfinished), len(batches)))
            return

        for name in (COMPARISON_NAME, summary_name(self.prefix)):
            written = concatenate([os.path.join(batch, name) for batch in batches],
                                  os.path.join(out_dir, name))
            self.logger.info('Wrote {:,} rows to {}.'.format(
                written, os.path.join(out_dir, name)))

    def run(self, gtdb_genome_path_file: str, taxonomy_file: str, out_dir: str) -> bool:
        """Predict the translation table of every genome of a release.

        Parameters
        ----------
        gtdb_genome_path_file : str
            genome_dirs file of the release, accession and genome directory per line.
        taxonomy_file : str
            Standardised NCBI taxonomy file, for the comparison.
        out_dir : str
            Directory the batches and their results are written to.

        @return: True once every batch this machine took has finished.
        """

        binary = shutil.which(GTRANSLATE_BIN)
        self.logger.info('Using {}.'.format(binary))

        taxonomy = read_taxonomy(taxonomy_file)
        self.logger.info('Read the NCBI taxonomy of {:,} genomes from {}.'.format(
            len(taxonomy), taxonomy_file))

        batches = self.plan_batches(gtdb_genome_path_file, out_dir)

        done, held, failed = 0, 0, 0
        for index, batch_dir in enumerate(batches, start=1):
            label = 'Batch {:,} of {:,} ({})'.format(
                index, len(batches), os.path.basename(batch_dir))

            if batch_state(batch_dir) == STATE_SUCCESS:
                self.logger.info('{}: already finished, skipping.'.format(label))
                continue

            if not claim_batch(batch_dir, self.reclaim):
                owner = read_canary(os.path.join(batch_dir, RUNNING_CANARY))
                held += 1
                self.logger.info('{}: held by {} since {}, skipping.'.format(
                    label, owner.get('host', 'another machine'),
                    owner.get('time', 'an unknown time')))
                continue

            self.logger.info('{}: starting.'.format(label))
            try:
                self.run_gtranslate(batch_dir)
                compared = self.compare_batch(batch_dir, taxonomy)
            except Exception as exc:
                failed += 1
                fail_batch(batch_dir, str(exc))
                self.logger.error('{}: failed and will be retried by a later '
                                  'run: {}'.format(label, exc))
                continue

            finish_batch(batch_dir, compared=compared)
            done += 1
            self.logger.info('{}: done.'.format(label))

        self.logger.info(
            '{:,} batch(es) finished here, {:,} held by another machine, '
            '{:,} failed.'.format(done, held, failed))

        self.aggregate(batches, out_dir)

        if failed:
            raise RuntimeError(
                '{:,} batch(es) failed; they are the directories holding a {} '
                'file and are retried by running the command again.'.format(
                    failed, FAILED_CANARY))

        return True
