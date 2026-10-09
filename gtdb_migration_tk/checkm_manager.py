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

"""Estimate the completeness and contamination of the genomes of a release with
CheckM (checkm) and CheckM2 (checkm2).

WHY THE WORK IS CUT INTO BATCHES

CheckM runs for hours over a thousand genomes and CheckM2 over a few thousand,
and a release brings tens of thousands of new ones, so the run is days and
belongs on several machines. The batching is the machinery trans_table,
prodigal, hmmsearch and trnascan share, in batching.py: the genomes are
partitioned once under --out_dir, and a machine claims a batch directory before
working on it, so several machines given the same --out_dir divide the work
between them without being told which genomes to take. Each batch keeps its
own log, and its SUCCESS canary is what says it is done.

Unlike those commands, what CheckM makes lives in the batch and not in the
genome directories, and nothing beside a genome can say it was assessed. So the
genomes to assess are decided once, when the batches are planned, from the
release report -- those update_genomes says are to be regenerated, or with --all
every genome the release holds -- and only those are planned
(batching.plan_batches(accessions=...)), so that a batch's size is the size of
its work. Once planned the batches are authoritative, as they are for every
batched command: a run with another report or --all over the same --out_dir
uses the plan already there. The batches are planned in --out_dir itself, and
the release files are written beside them once every batch has succeeded, and
never a part of them that reads like the whole. checkm and checkm2 are each
given an --out_dir of their own: one that holds batch directories another
command planned is refused, since that command's SUCCESS would tell this one
it had nothing to do.

THE RELEASE FILES

The release tables are gzipped (checkm.profiles.tsv.gz, checkm.qa_sh100.tsv.gz,
checkm2.quality_report.tsv.gz), being a row per genome of the release; what says
which genomes were left out, and the version files, are a few lines meant to be
read, and are not. A table an earlier run left uncompressed is removed once the
gzipped one is written (SUPERSEDED_RELEASE_FILES), as trans_table removes its
own: two tables a few genomes apart, one of them stale, is how the stale one
gets read.

The alignments of multi-copy genes CheckM writes (qa -a) stay in their batches
and are not gathered for the release. Nothing reads them: they were joined into
checkm.alignment_file.tsv until 0.1.48, and for r237 that was 142 GB, 137 GB of
it from one batch, copied again by every run that wrote the release files. They
are the record of the genes behind each contamination estimate, so CheckM is
still asked for them; checkm.alignment_file.tsv is removed as superseded, since
one an earlier run left would not hold the batches made since.

THE PROTEINS

Both commands are handed the proteins the prodigal command called (--genes),
and take them as correct: they are called once for a release, by one program,
under the table trans_table chose, and what CheckM and CheckM2 assess is then
the same proteome hmmsearch and the metadata are made from. Left to call genes
themselves, both would call them under whichever of tables 4 and 11 codes more
of the genome, a rule that cannot express table 25.

Handed genes, CheckM2 reports only its estimates, not the statistics it takes
from calling genes itself (the table used, coding density, genome size, GC,
N50 and the rest); its models read their features from the proteins either
way.

ONE BAD GENOME DOES NOT COST A BATCH

A genome with no protein file, or an empty one, is named in
the batch's not_assessed.tsv, gathered into checkm_not_assessed.tsv or
checkm2_not_assessed.tsv for the release, and left; the program is run over the
rest. The program failing fails the batch, which the next run repeats.

A genome larger than --max_genome_size is named and left too, as trnascan,
rna_silva and rna_ltp leave one: it is a metagenome deposited as one genome,
whose completeness and contamination say nothing about any one organism. It is sized when its batch is assessed rather than when the batches are
planned, so that a limit given to a run over batches already planned applies to
them. A batch that finishes without it is done, and a larger limit does not go
back to it: the genome is assessed later by removing that batch's SUCCESS.

WHICH VERSION MADE THE RELEASE

The program's version is asked once where a run starts and logged, and written
beside each batch it assessed (checkm.version, checkm2.version), as for every
external program. The release files are made from batches that several machines
ran, over days, and the machine that writes them may have run none of them, so
the version written beside them is gathered from the batches' own files, not
asked of the program the last machine runs. Where the batches disagree every
version is written, one to a line, and the batches each made are named in the
log: a release whose estimates come from two versions of CheckM says so rather
than claiming one.

PPLACER

CheckM places each genome in its reference tree with pplacer, and which pplacer
matters: 1.1.alpha20, which checkm-genome 1.2.5 pins, dies part-way through
placing some batches, and r237 was finished with alpha22. Its version is the one
the run cannot ask. CheckM runs pplacer by bare name, from the PATH CheckM was
started with, which a wrapper that puts CheckM's environment on PATH sets inside
itself where the toolkit never sees it; and bioconda's pplacer answers --version
with 'dev' (or the git describe of whatever directory it runs in), alpha20 and
alpha22 alike. So nothing is asked of a pplacer before the run. While each
CheckM step runs, the processes it starts are watched, and the executable of
any pplacer among them is read from /proc: that is the pplacer that placed the
batch, wherever it was found. Its version is the conda package it was installed
from (utils.common.conda_package_version()), written as pplacer.version beside
checkm.version. A pplacer whose version cannot be told is named in the batch's
log and given no version file, which is better than one that guesses.

The release's pplacer.version is gathered from the batches as checkm.version
is, with one more rule: a batch that ran CheckM but recorded no pplacer -- made
before pplacer was recorded, or in which pplacer was never seen -- leaves the
release with no pplacer.version at all, since a version claimed for the whole
release would be a guess about that batch.
"""

import logging
import os
import shutil
import subprocess
import threading
from collections import defaultdict
from typing import Dict, List, NamedTuple, Optional, Sequence, Set, Tuple


from gtdb_migration_tk.batching import (BATCH_DIR_PREFIX, CLAIM_LEASE_SECONDS, HEARTBEAT_SECONDS,
                                        RUNNING_CANARY, STATE_SUCCESS,
                                        SUCCESS_CANARY, STAT_THREADS,
                                        BatchLayout, Heartbeat, age_phrase,
                                        batch_log, batch_state, batchfile_path,
                                        claim_age, claim_batch, concatenate,
                                        fail_batch, finish_batch,
                                        plan_batches, read_batchfile, read_canary,
                                        release_claim, split_by_fasta,
                                        split_by_genome_size, tally_reasons,
                                        write_table)
from gtdb_migration_tk.biolib_lite.common import make_sure_path_exists
from gtdb_migration_tk.biolib_lite.external.execute import check_dependencies
from gtdb_migration_tk.update_genomes import (genomes_in_release,
                                             genomes_to_regenerate)
from gtdb_migration_tk.utils.common import (DEFAULT_MAX_GENOME_SIZE, MBP,
                                            conda_package_version,
                                            protein_fasta,
                                            record_program_version,
                                            version_file, write_version_file)

# The programs, as they are called. Each one's version goes into the log and, as
# checkm.version or checkm2.version, into each batch it was run over and, once
# the release files are written, beside them.
CHECKM = 'checkm'
CHECKM2 = 'checkm2'

# CheckM's placement program, whose version is learned from the process CheckM
# starts rather than asked (see PPLACER), and recorded as pplacer.version. The
# watch looks for it this often: placement takes minutes, so it is not missed.
PPLACER = 'pplacer'
WATCH_SECONDS = 2.0

# The most batches a log line names; the rest are counted.
NAMED_BATCHES = 20

# What each command calls the files of its batches. There is no older name to
# look for: neither command has had batches before.
CHECKM_LAYOUT = BatchLayout(batchfiles=('checkm_batchfile.tsv.gz',), log='checkm.log')
CHECKM2_LAYOUT = BatchLayout(batchfiles=('checkm2_batchfile.tsv.gz',), log='checkm2.log')

# Genomes per batch. CheckM2 loads its models and searches the whole of its
# DIAMOND database once per run whatever the number of genomes, so its batches
# are larger; CheckM's are the chunks of 1,000 it has always been run over.
DEFAULT_CHECKM_BATCH_SIZE = 1000
DEFAULT_CHECKM2_BATCH_SIZE = 5000

# The most threads pplacer is given, whatever --cpus says. pplacer runs badly
# on more than 64, so a machine given -c 96 runs the rest of CheckM on 96 and
# places the genomes in the tree on 64.
PPLACER_MAX_THREADS = 64

# The genomes a batch left out, and why, and the same gathered for the release.
NOT_ASSESSED_NAME = 'not_assessed.tsv'
CHECKM_NOT_ASSESSED = 'checkm_not_assessed.tsv'
CHECKM2_NOT_ASSESSED = 'checkm2_not_assessed.tsv'
NOT_ASSESSED_HEADER = ('genome_id', 'reason', 'detail')
REASON_NO_PROTEINS = 'no_protein_file'
REASON_GENOME_TOO_LARGE = 'genome_too_large'

# Within a checkm batch: the proteins linked in for CheckM, what it writes, and
# the tables made from that, with the release file each is gathered into,
# gzipped. The alignment of multi-copy genes stays in its batch (see THE RELEASE
# FILES).
CHECKM_INPUT_DIR = 'input'
CHECKM_OUTPUT_DIR = 'checkm'
CHECKM_PROTEIN_EXT = 'faa.gz'
CHECKM_TREE_QA = 'tree_qa.o2.tsv'
CHECKM_QA = 'qa.tsv'
CHECKM_PROFILE = 'profile.tsv'
CHECKM_QA_SH100 = 'qa_sh100.tsv'
CHECKM_ALIGNMENT = 'alignment_file.tsv'
CHECKM_RELEASE_PROFILE = 'checkm.profiles.tsv.gz'
CHECKM_RELEASE_QA_SH100 = 'checkm.qa_sh100.tsv.gz'
CHECKM_RELEASE_FILES = ((CHECKM_PROFILE, CHECKM_RELEASE_PROFILE),
                        (CHECKM_QA_SH100, CHECKM_RELEASE_QA_SH100))

# What an earlier run wrote for the release that the files above replace: the
# tables uncompressed, and the batches' alignments joined
CHECKM_SUPERSEDED_FILES = ('checkm.profiles.tsv', 'checkm.qa_sh100.tsv',
                           'checkm.alignment_file.tsv')

# Within a checkm2 batch: the proteins linked in, what CheckM2 wrote, and a
# copy of its report beside the batch's other files. The input is beside the output rather than in
# it, since --force empties the output directory before CheckM2 starts. A
# genome is linked as <accession>.faa.gz: CheckM2 labels a genome with its
# file's basename less '.faa.gz', which gives the report the accession alone
# rather than prodigal's <accession>_protein.
CHECKM2_INPUT_DIR = 'input'
CHECKM2_OUTPUT_DIR = 'checkm2'
CHECKM2_LINK_EXT = '.faa.gz'
CHECKM2_REPORT = 'quality_report.tsv'
CHECKM2_BATCH_REPORT = 'checkm2.quality_report.tsv'
CHECKM2_RELEASE_REPORT = 'checkm2.quality_report.tsv.gz'
CHECKM2_SUPERSEDED_FILES = ('checkm2.quality_report.tsv',)


class BatchCounts(NamedTuple):
    """What assessing one batch came to, recorded in its SUCCESS canary."""

    assessed: int
    not_assessed: int


def accessions_rows(genome_report: str, all_genomes: bool = False) -> Optional[Set[str]]:
    """The genomes whose quality is to be estimated, from the release report.

    Parameters
    ----------
    genome_report : str
        update_genomes report.log, or 'none' for every genome in the release.
    all_genomes : bool
        Every genome the release holds, not only those to regenerate.

    @return: the accessions, or None for every genome of the genome_dirs file.
    """

    if genome_report.lower() == 'none':
        return None

    # --all asks for every genome the release HOLDS, not every row of the
    # report: a removed genome, and one that could not be compared, have no
    # directory to look in
    if all_genomes:
        return set(genomes_in_release(genome_report))

    return set(genomes_to_regenerate(genome_report))


def foreign_batches(out_dir: str, layout: BatchLayout) -> List[str]:
    """Batch directories under an output directory that another command planned.

    Parameters
    ----------
    out_dir : str
        Output directory of the run.
    layout : BatchLayout
        What this command calls the files of its batches.

    @return: names of the batch directories holding none of its batchfiles.
    """

    if not os.path.isdir(out_dir):
        return []

    return sorted(name for name in os.listdir(out_dir)
                  if name.startswith(BATCH_DIR_PREFIX)
                  and os.path.isdir(os.path.join(out_dir, name))
                  and not os.path.exists(batchfile_path(os.path.join(out_dir, name), layout)))


def checkm2_command(proteins: Sequence[str], out_dir: str,
                    threads: int, tmp_dir: str) -> List[str]:
    """The CheckM2 command run over the proteins of one batch.

    No --ttable: CheckM2 ignores it when handed genes, the proteins having been
    called under their table by prodigal already. The files are named on the command line,
    not given as a directory: CheckM2 decompresses gzipped files only when they
    are listed.

    Parameters
    ----------
    proteins : sequence of str
        Staged protein files, named for their accessions.
    out_dir : str
        Directory CheckM2 writes its report and intermediates to.
    threads : int
        Threads CheckM2 is given.
    tmp_dir : str
        Directory CheckM2 decompresses the proteins into.

    @return: the command as a list of arguments, ready for subprocess.
    """

    # --input takes the rest of the command line, so it goes last
    return [CHECKM2, 'predict',
            '--genes',
            '--threads', str(threads),
            '--tmpdir', tmp_dir,
            '--remove_intermediates',
            '--force',
            '--output-directory', out_dir,
            '--input'] + list(proteins)


def join_tables(tables: Sequence[str], path: str) -> None:
    """Join CheckM tables on their first column, the bin, as checkm join_tables does.

    Parameters
    ----------
    tables : sequence of str
        Tab-separated tables with a header, the bin ID first.
    path : str
        File to write.

    @return: None
    """

    headers = []
    rows = defaultdict(dict)
    bins = []
    for table in tables:
        with open(table) as handle:
            header = [field.strip() for field in handle.readline().split('\t')][1:]
            headers.append(header)
            for line in handle:
                fields = [field.strip() for field in line.split('\t')]
                if fields[0] not in rows:
                    bins.append(fields[0])
                for name, value in zip(header, fields[1:]):
                    rows[fields[0]][name] = value

    with open(path, 'w') as handle:
        handle.write('\t'.join(['Bin Id'] + [name for header in headers
                                             for name in header]) + '\n')
        for bin_id in bins:
            handle.write('\t'.join([bin_id] + [rows[bin_id].get(name, '')
                                               for header in headers
                                               for name in header]) + '\n')


def descendant_executables(pid: int) -> Dict[int, str]:
    """The executables of every process descended from one.

    Read from /proc, so Linux alone; elsewhere there are none. A process that
    ends while it is being read, or whose executable cannot be read, is left
    out.

    Parameters
    ----------
    pid : int
        The process whose descendants are wanted.

    @return: each descendant's executable, the symlinks followed, by process ID.
    """

    children = defaultdict(list)
    try:
        entries = os.listdir('/proc')
    except OSError:
        return {}
    for entry in entries:
        if not entry.isdigit():
            continue
        try:
            with open(os.path.join('/proc', entry, 'stat')) as handle:
                stat = handle.read()
        except OSError:
            continue
        # the name between the parentheses may hold spaces and parentheses of its
        # own; the state and the parent's ID follow the last ')'
        try:
            children[int(stat[stat.rindex(')') + 1:].split()[1])].append(int(entry))
        except (ValueError, IndexError):
            continue

    executables = {}
    pending = list(children.get(pid, ()))
    while pending:
        child = pending.pop()
        pending.extend(children.get(child, ()))
        try:
            executables[child] = os.path.realpath(os.readlink(
                os.path.join('/proc', str(child), 'exe')))
        except OSError:
            continue

    return executables


def is_program(executable: str, program: str) -> bool:
    """Whether an executable is a program, as bioconda installs it.

    bioconda's bin/pplacer is a link to bin/pplacer.exe, which is what runs.

    @return: True where the executable is <program> or <program>.<ext>.
    """

    return os.path.splitext(os.path.basename(executable))[0] == program


class ProgramWatch(object):
    """The executables a program ran as, among the descendants of a process.

    A context manager: a thread looks every WATCH_SECONDS while it is open, and
    once more as it closes.
    """

    def __init__(self, pid: int, program: str, interval: Optional[float] = None) -> None:
        """Initialization.

        Parameters
        ----------
        pid : int
            The process whose descendants are watched.
        program : str
            The program looked for, as it is called.
        interval : float
            Seconds between looks; WATCH_SECONDS by default.

        @return: None
        """

        self.pid = pid
        self.program = program
        self.interval = WATCH_SECONDS if interval is None else interval
        self.executables: Set[str] = set()
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._watch, daemon=True)

    def look(self) -> None:
        for executable in descendant_executables(self.pid).values():
            if is_program(executable, self.program):
                self.executables.add(executable)

    def _watch(self) -> None:
        while True:
            self.look()
            if self._stop.wait(self.interval):
                return

    def __enter__(self) -> 'ProgramWatch':
        self._thread.start()
        return self

    def __exit__(self, *exc) -> None:
        self._stop.set()
        self._thread.join()


def name_batches(names: Sequence[str]) -> str:
    """Batches as a log line names them: the first NAMED_BATCHES, the rest counted.

    @return: e.g. 'batch_000001, batch_000002 and 303 more'.
    """

    named = ', '.join(names[:NAMED_BATCHES])
    if len(names) > NAMED_BATCHES:
        named += ' and {:,} more'.format(len(names) - NAMED_BATCHES)
    return named


def batch_versions(batches: Sequence[str], program: str) -> Dict[str, List[str]]:
    """The versions of a program that made a release's batches.

    A batch in which no genome could be assessed never ran the program and has
    no version file, and so is in none of them.

    Parameters
    ----------
    batches : sequence of str
        Batch directories, in batch order.
    program : str
        The program, as it is called.

    @return: each version, as its batches' files state it, with the names of the
             batches it made in batch order.
    """

    versions = defaultdict(list)
    for batch_dir in batches:
        path = version_file(batch_dir, program)
        if os.path.exists(path):
            with open(path) as handle:
                versions[handle.read().strip()].append(os.path.basename(batch_dir))

    return dict(versions)


class BatchedQuality(object):
    """What checkm and checkm2 share: the batches, the claims, and the release files.

    A subclass names its program and layout, the file a batch is planned
    around, how a batch is assessed, the batch tables concatenated into the
    release's (RELEASE_TABLES, as (batch file, release file), the release file
    gzipped), and what an earlier run wrote that those replace
    (SUPERSEDED_RELEASE_FILES).
    """

    PROGRAM = None
    # every program whose version is gathered for the release; PROGRAM first
    RECORDED_PROGRAMS = ()
    LAYOUT = None
    RELEASE_TABLES = ()
    SUPERSEDED_RELEASE_FILES = ()
    NOT_ASSESSED_RELEASE = None

    def __init__(self,
                 cpus: int = 1,
                 tmp_dir: str = '/tmp/',
                 batch_size: int = DEFAULT_CHECKM_BATCH_SIZE,
                 reclaim: bool = False,
                 lease: float = CLAIM_LEASE_SECONDS,
                 heartbeat: float = HEARTBEAT_SECONDS,
                 max_genome_size: float = DEFAULT_MAX_GENOME_SIZE) -> None:
        """Initialization.

        Parameters
        ----------
        cpus : int
            Threads the program is given.
        tmp_dir : str
            Directory for the program's scratch files; no results are written here.
        batch_size : int
            Genomes per batch.
        reclaim : bool
            Take over a batch another machine holds before its claim has expired.
        lease : float
            Seconds a claim survives without the machine holding it saying so.
        heartbeat : float
            Seconds between this machine saying so about a batch of its own.
        max_genome_size : float
            Largest genome assembly assessed, in Mbp.

        @return: None
        """

        check_dependencies([self.PROGRAM])

        self.logger = logging.getLogger('timestamp')
        self.version = record_program_version(self.PROGRAM)

        self.cpus = cpus
        self.tmp_dir = tmp_dir
        self.batch_size = batch_size
        self.reclaim = reclaim
        self.lease = lease
        self.heartbeat = heartbeat
        self.max_genome_bases = int(max_genome_size * MBP)

        # made here rather than by the first batch that wants it, so that a
        # --tmp_dir that cannot be made is met before a batch is claimed
        make_sure_path_exists(self.tmp_dir)

    def genome_file(self, accession: str, genome_dir: str) -> str:
        """The file a batch is planned around: what the program reads."""
        raise NotImplementedError

    def assess_batch(self, batch_dir: str, rows: Sequence[Tuple[str, str]]) -> int:
        """Run the program over the genomes of a batch that have proteins.

        @return: the number of genomes assessed.
        """
        raise NotImplementedError

    def run(self,
            gtdb_genome_path_file: str,
            genome_report: str,
            out_dir: str,
            all_genomes: bool = False) -> bool:
        """Estimate the quality of the genomes of a release, batch by batch.

        Parameters
        ----------
        gtdb_genome_path_file : str
            genome_dirs file of the release.
        genome_report : str
            update_genomes report.log, or 'none' for every genome in the release.
        out_dir : str
            Directory the batches and the release files are written to.
        all_genomes : bool
            Every genome the release holds, not only those to regenerate. It
            decides what is planned, and so has no effect on batches already
            planned.

        @return: True where every batch this machine took finished, False where
                 one failed and is left to a later run.
        """

        foreign = foreign_batches(out_dir, self.LAYOUT)
        if foreign:
            raise RuntimeError(
                '{} holds {:,} batch director(ies) {} did not plan, e.g. {}; give '
                'each command an --out_dir of its own.'.format(
                    out_dir, len(foreign), self.PROGRAM, foreign[0]))

        accessions = accessions_rows(genome_report, all_genomes)
        if accessions is not None:
            self.logger.info('The report names {:,} genome(s) to estimate the '
                             'quality of.'.format(len(accessions)))
        batches = plan_batches(gtdb_genome_path_file, out_dir, self.batch_size,
                               self.LAYOUT, self.logger,
                               genome_file=self.genome_file, accessions=accessions)

        done, held, failed = 0, 0, 0
        for index, batch_dir in enumerate(batches, start=1):
            label = 'Batch {:,} of {:,} ({})'.format(
                index, len(batches), os.path.basename(batch_dir))

            if batch_state(batch_dir) == STATE_SUCCESS:
                self.logger.info('{}: already finished, skipping.'.format(label))
                continue

            if not claim_batch(batch_dir, self.reclaim, self.lease):
                owner = read_canary(os.path.join(batch_dir, RUNNING_CANARY))
                held += 1
                self.logger.info('{}: held by {} since {}, last heard from {}, '
                                 'skipping.'.format(
                                     label, owner.get('host', 'another machine'),
                                     owner.get('time', 'an unknown time'),
                                     age_phrase(claim_age(
                                         os.path.join(batch_dir, RUNNING_CANARY)))))
                continue

            # the batch has its own log from here, since this is where anything
            # happens to it and every machine of a run writes its own --log
            with batch_log(batch_dir, self.logger, self.LAYOUT):
                self.logger.info('{}: starting with {} ({}).'.format(
                    label, self.PROGRAM, self.version))
                try:
                    with Heartbeat(os.path.join(batch_dir, RUNNING_CANARY),
                                   self.heartbeat):
                        counts = self.process_batch(batch_dir)
                except KeyboardInterrupt:
                    # the machine holding it is stopping, so the batch is handed
                    # back rather than left to sit out its lease
                    release_claim(batch_dir)
                    self.logger.error('{}: interrupted; the claim is given up.'.format(label))
                    raise
                except Exception as exc:
                    failed += 1
                    fail_batch(batch_dir, str(exc))
                    self.logger.error('{}: failed and will be retried by a later '
                                      'run: {}'.format(label, exc))
                    continue

                finish_batch(batch_dir, assessed=counts.assessed,
                             not_assessed=counts.not_assessed)
                done += 1
                self.logger.info('{}: done.'.format(label))

        self.logger.info(
            '{:,} batch(es) finished here, {:,} held by another machine, '
            '{:,} failed.'.format(done, held, failed))

        self.aggregate(batches, out_dir)

        # a batch that failed has already said why, in its own log and in its
        # FAILED file
        if failed:
            self.logger.error(
                '{:,} batch(es) failed; they are the directories holding a FAILED '
                'file and are retried by running the command again.'.format(failed))
            return False

        return True

    def process_batch(self, batch_dir: str) -> BatchCounts:
        """Assess one batch, and name the genomes it left out.

        Parameters
        ----------
        batch_dir : str
            Batch directory.

        @return: what the batch came to, which run() records in its canary.
        """

        rows = read_batchfile(batchfile_path(batch_dir, self.LAYOUT))

        # a genome whose file is not there has nothing to assess; it is named and
        # left rather than stopping the rest of the batch
        present, missing = split_by_fasta(rows, STAT_THREADS)
        not_assessed = [(accession, REASON_NO_PROTEINS, '') for accession in missing]

        # a batch's file is <genome dir>/prodigal/<accession>_protein.faa.gz
        present, too_large = split_by_genome_size(
            present, lambda row: os.path.dirname(os.path.dirname(row[0])),
            self.max_genome_bases, self.logger, STAT_THREADS)
        not_assessed.extend((accession, REASON_GENOME_TOO_LARGE, str(size))
                            for (_, accession), size in too_large)

        assessed = self.assess_batch(batch_dir, present)

        not_assessed.sort()
        write_table(not_assessed, os.path.join(batch_dir, NOT_ASSESSED_NAME),
                    header=NOT_ASSESSED_HEADER)
        if not_assessed:
            self.logger.warning(
                '{:,} genome(s) of this batch were not assessed: {}.'.format(
                    len(not_assessed),
                    '; '.join('{:,} {}'.format(count, reason)
                              for reason, count in sorted(tally_reasons(
                                  [(gid, reason) for gid, reason, _ in not_assessed]).items()))))

        return BatchCounts(assessed=assessed, not_assessed=len(not_assessed))

    def aggregate(self, batches: Sequence[str], out_dir: str) -> None:
        """Write the release files from the batches, once every batch has succeeded.

        As trans_table writes its own: written only once every batch has
        succeeded, so that the files at the top of the output directory are
        either the whole release or absent, and never a part of it that reads
        like the whole; each table is the batches' own concatenated in batch
        order under a single header; and whichever machine finishes last writes
        them, having counted the release from every batch's SUCCESS canary.

        Parameters
        ----------
        batches : sequence of str
            Every batch directory of the run, in batch order.
        out_dir : str
            Directory the release files are written to.

        @return: None
        """

        unfinished = [batch for batch in batches if batch_state(batch) != STATE_SUCCESS]
        if unfinished:
            self.logger.info(
                '{:,} of {:,} batch(es) are done; the release files are written '
                'once they all are.'.format(len(batches) - len(unfinished), len(batches)))
            return

        for batch_name, release_name in self.RELEASE_TABLES:
            # a batch in which no genome could be assessed never ran the program
            # and so has no table to add
            path = os.path.join(out_dir, release_name)
            written = concatenate([os.path.join(batch, batch_name) for batch in batches
                                   if os.path.exists(os.path.join(batch, batch_name))],
                                  path, compress=True)
            self.logger.info('Wrote {:,} rows to {}.'.format(written, path))

        self.remove_superseded(out_dir)
        self.report_not_assessed(batches, out_dir)
        versions = [(program, self.write_release_version(batches, out_dir, program))
                    for program in self.RECORDED_PROGRAMS]

        assessed = 0
        for batch_dir in batches:
            try:
                assessed += int(read_canary(os.path.join(batch_dir, SUCCESS_CANARY))['assessed'])
            except (KeyError, ValueError):
                pass
        self.logger.info('Release: {:,} genome(s) assessed with {}.'.format(
            assessed, ', '.join(
                '{} ({})'.format(program, '; '.join(found)) if found else program
                for program, found in versions
                if found or program == self.PROGRAM)))

    def remove_superseded(self, out_dir: str) -> None:
        """Remove what an earlier run wrote for the release that this run's files replace.

        Called once the release tables are written, so that a run that fails
        before then leaves the earlier files as they were. See THE RELEASE FILES.

        Parameters
        ----------
        out_dir : str
            Directory the release files are written to.

        @return: None
        """

        for name in self.SUPERSEDED_RELEASE_FILES:
            path = os.path.join(out_dir, name)
            if os.path.exists(path):
                os.remove(path)
                self.logger.info('Removed {}, which an earlier run wrote and this '
                                 'run\'s release files replace.'.format(path))

    def write_release_version(self, batches: Sequence[str], out_dir: str,
                              program: str) -> List[str]:
        """Record beside the release files the version of a program that made them.

        Gathered from the batches' version files rather than taken from this
        run, since other machines ran most of the batches; see WHICH VERSION
        MADE THE RELEASE. Where no batch ran the program, nothing made the
        release files and no version is written, and one an earlier aggregation
        left is removed. A program the command runs through another, as CheckM
        runs pplacer, is written only where every batch that ran the command
        recorded it (see PPLACER).

        Parameters
        ----------
        batches : sequence of str
            Every batch directory of the run, in batch order.
        out_dir : str
            Directory the release files are written to.
        program : str
            One of RECORDED_PROGRAMS.

        @return: the versions written, in order; none where none was.
        """

        versions = batch_versions(batches, program)
        path = version_file(out_dir, program)

        unrecorded = []
        if program != self.PROGRAM:
            recorded = {name for names in versions.values() for name in names}
            unrecorded = [name for names in batch_versions(batches, self.PROGRAM).values()
                          for name in names if name not in recorded]

        if not versions or unrecorded:
            if os.path.exists(path):
                os.remove(path)
            if unrecorded:
                self.logger.warning(
                    '{:,} batch(es) that ran {} recorded no {} version, so none '
                    'is written for the release: {}. They were made before {} was '
                    'recorded, or {} was not seen to run; their logs say which. A batch '
                    'is made again by removing its SUCCESS.'.format(
                        len(unrecorded), self.PROGRAM, program,
                        name_batches(sorted(unrecorded)), program, program))
            return []

        ordered = sorted(versions)
        write_version_file(out_dir, program, '\n'.join(ordered))
        if len(ordered) == 1:
            self.logger.info('The release was made with {} {}; wrote {}.'.format(
                program, ordered[0], path))
        else:
            self.logger.warning(
                'The batches of the release were made with {:,} versions of '
                '{}, all of them written to {}: {}. A batch is made again, with the '
                'version this machine runs, by removing its SUCCESS.'.format(
                    len(ordered), program, path,
                    '; '.join('{} by {:,} batch(es): {}'.format(
                        version, len(versions[version]), name_batches(versions[version]))
                        for version in ordered)))

        return ordered

    def report_not_assessed(self, batches: Sequence[str], out_dir: str) -> None:
        """Name the genomes of the release that were not assessed, and why.

        Each batch names its own in not_assessed.tsv; this is the one file that
        says it for the release. It is written whether or not there are any, so
        that a release with nothing left out says so rather than leaving the
        question open.

        Parameters
        ----------
        batches : sequence of str
            Every batch directory of the run.
        out_dir : str
            Directory the release files are written to.

        @return: None
        """

        rows = []
        for batch in batches:
            with open(os.path.join(batch, NOT_ASSESSED_NAME)) as handle:
                handle.readline()
                rows.extend(tuple(line.rstrip('\n').split('\t'))
                            for line in handle if line.strip())
        rows.sort()

        path = os.path.join(out_dir, self.NOT_ASSESSED_RELEASE)
        write_table(rows, path, header=NOT_ASSESSED_HEADER)

        if not rows:
            self.logger.info('Every genome was assessed; wrote {} with no rows.'.format(path))
            return

        by_reason = tally_reasons([(row[0], row[1]) for row in rows])
        self.logger.warning(
            '{:,} genome(s) of the release were not assessed and are named '
            'in {}: {}.'.format(
                len(rows), path,
                '; '.join('{:,} {}'.format(count, reason)
                          for reason, count in sorted(by_reason.items()))))

    def run_program(self, cmd: Sequence[str], watch: Optional[str] = None) -> Set[str]:
        """Run one step of the program, failing the batch where it fails.

        Parameters
        ----------
        cmd : sequence of str
            The command.
        watch : str
            A program the step may start, whose executables are to be learned.

        @return: the executables the watched program ran as; none where nothing
                 was watched.
        """

        self.logger.info('Command: {}'.format(' '.join(cmd[:12])
                                              + (' ...' if len(cmd) > 12 else '')))
        silent = getattr(self.logger, 'is_silent', False)
        proc = subprocess.Popen(list(cmd),
                                stdout=subprocess.DEVNULL if silent else None,
                                stderr=subprocess.STDOUT if silent else None)
        executables = set()
        try:
            if watch:
                with ProgramWatch(proc.pid, watch) as watching:
                    proc.wait()
                executables = watching.executables
            else:
                proc.wait()
        finally:
            # a step interrupted is not left running behind the batch it was for
            if proc.poll() is None:
                proc.kill()
                proc.wait()

        if proc.returncode != 0:
            raise RuntimeError('{} {} returned exit code {}.'.format(
                cmd[0], cmd[1], proc.returncode))

        return executables


class CheckM(BatchedQuality):
    """CheckM over the Prodigal proteins of the genomes of a release."""

    PROGRAM = CHECKM
    RECORDED_PROGRAMS = (CHECKM, PPLACER)
    LAYOUT = CHECKM_LAYOUT
    RELEASE_TABLES = CHECKM_RELEASE_FILES
    SUPERSEDED_RELEASE_FILES = CHECKM_SUPERSEDED_FILES
    NOT_ASSESSED_RELEASE = CHECKM_NOT_ASSESSED

    def genome_file(self, accession: str, genome_dir: str) -> str:
        return protein_fasta(accession, genome_dir)

    def assess_batch(self, batch_dir, rows):
        # what an earlier attempt at the batch left is removed, so a retry
        # starts from nothing rather than from a CheckM run stopped part-way
        input_dir = os.path.join(batch_dir, CHECKM_INPUT_DIR)
        output_dir = os.path.join(batch_dir, CHECKM_OUTPUT_DIR)
        for path in [input_dir, output_dir] + [os.path.join(batch_dir, name) for name in (
                CHECKM_TREE_QA, CHECKM_QA, CHECKM_PROFILE, CHECKM_QA_SH100, CHECKM_ALIGNMENT)] + [
                version_file(batch_dir, program) for program in self.RECORDED_PROGRAMS]:
            if os.path.isdir(path):
                shutil.rmtree(path)
            elif os.path.exists(path):
                os.remove(path)

        if not rows:
            return 0

        # linked under the name prodigal gave them, which CheckM takes the bin ID
        # from, so the tables name <accession>_protein as they always have
        os.makedirs(input_dir)
        for gene_file, _ in rows:
            os.symlink(os.path.abspath(gene_file),
                       os.path.join(input_dir, os.path.basename(gene_file)))

        threads = str(self.cpus)
        pplacer_threads = str(min(self.cpus, PPLACER_MAX_THREADS))
        lineage_ms = os.path.join(output_dir, 'lineage.ms')
        tree_qa = os.path.join(batch_dir, CHECKM_TREE_QA)
        qa = os.path.join(batch_dir, CHECKM_QA)
        # pplacer places the genomes in lineage_wf; every step is watched, so
        # that a pplacer run anywhere is not missed
        pplacers = set()
        pplacers |= self.run_program(
            [CHECKM, 'lineage_wf', '--pplacer_threads', pplacer_threads, '--genes',
             '-x', CHECKM_PROTEIN_EXT, '-t', threads, '--tmpdir', self.tmp_dir,
             input_dir, output_dir], watch=PPLACER)
        pplacers |= self.run_program(
            [CHECKM, 'tree_qa', '-o', '2', '--tab_table', '-f', tree_qa,
             '--tmpdir', self.tmp_dir, output_dir], watch=PPLACER)
        pplacers |= self.run_program(
            [CHECKM, 'qa', '-t', threads, '--tab_table', '-f', qa,
             '--tmpdir', self.tmp_dir, lineage_ms, output_dir], watch=PPLACER)
        join_tables([qa, tree_qa], os.path.join(batch_dir, CHECKM_PROFILE))
        pplacers |= self.run_program(
            [CHECKM, 'qa', '--aai_strain', '0.9999', '-t', threads,
             '-a', os.path.join(batch_dir, CHECKM_ALIGNMENT),
             '--tab_table', '-f', os.path.join(batch_dir, CHECKM_QA_SH100),
             '--tmpdir', self.tmp_dir, lineage_ms, output_dir], watch=PPLACER)

        shutil.rmtree(input_dir)
        write_version_file(batch_dir, CHECKM, self.version)
        self.record_pplacer(batch_dir, pplacers)

        return len(rows)

    def record_pplacer(self, batch_dir: str, executables: Set[str]) -> None:
        """Write beside a batch the version of the pplacer that placed its genomes.

        Parameters
        ----------
        batch_dir : str
            Batch directory.
        executables : set of str
            The executables pplacer ran as, as the watch saw them.

        @return: None
        """

        if not executables:
            self.logger.warning(
                'pplacer was not seen to run, so no {} version is recorded '
                'for this batch.'.format(PPLACER))
            return

        versions = set()
        for executable in sorted(executables):
            version = conda_package_version(executable, PPLACER)
            if version is None:
                self.logger.warning(
                    'pplacer ran from {}, which is not in a conda environment '
                    'holding the {} package, so its version cannot be told and none is '
                    'recorded for this batch.'.format(executable, PPLACER))
                return
            self.logger.info('Using {}: {}, run from {}.'.format(PPLACER, version, executable))
            versions.add(version)

        write_version_file(batch_dir, PPLACER, '\n'.join(sorted(versions)))


class CheckM2(BatchedQuality):
    """CheckM2 over the Prodigal proteins of the genomes of a release."""

    PROGRAM = CHECKM2
    RECORDED_PROGRAMS = (CHECKM2,)
    LAYOUT = CHECKM2_LAYOUT
    RELEASE_TABLES = ((CHECKM2_BATCH_REPORT, CHECKM2_RELEASE_REPORT),)
    SUPERSEDED_RELEASE_FILES = CHECKM2_SUPERSEDED_FILES
    NOT_ASSESSED_RELEASE = CHECKM2_NOT_ASSESSED

    def genome_file(self, accession: str, genome_dir: str) -> str:
        return protein_fasta(accession, genome_dir)

    def assess_batch(self, batch_dir, rows):
        # what an earlier attempt at the batch left is removed first
        input_dir = os.path.join(batch_dir, CHECKM2_INPUT_DIR)
        output_dir = os.path.join(batch_dir, CHECKM2_OUTPUT_DIR)
        batch_report = os.path.join(batch_dir, CHECKM2_BATCH_REPORT)
        for path in (input_dir, output_dir):
            if os.path.isdir(path):
                shutil.rmtree(path)
        for path in (batch_report, version_file(batch_dir, CHECKM2)):
            if os.path.exists(path):
                os.remove(path)

        if not rows:
            return 0

        os.makedirs(input_dir)
        links = []
        for gene_file, accession in rows:
            link = os.path.join(input_dir, accession + CHECKM2_LINK_EXT)
            os.symlink(os.path.abspath(gene_file), link)
            links.append(link)

        self.run_program(checkm2_command(links, output_dir, self.cpus, self.tmp_dir))

        # CheckM2 writes its report last, so a run that exited cleanly without
        # one did not finish
        report = os.path.join(output_dir, CHECKM2_REPORT)
        if not os.path.exists(report):
            raise RuntimeError('{} wrote no {}.'.format(CHECKM2, CHECKM2_REPORT))

        shutil.copyfile(report, batch_report)
        shutil.rmtree(input_dir)
        write_version_file(batch_dir, CHECKM2, self.version)

        return len(rows)
