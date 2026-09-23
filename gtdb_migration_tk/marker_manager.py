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
import os
import sys
import gzip
import logging
import multiprocessing as mp
import shutil
import subprocess
import tempfile
from collections import defaultdict
from multiprocessing.queues import Queue
from typing import Callable, Dict, List, NamedTuple, Optional, Sequence, Set, Tuple

from tqdm import tqdm

from gtdb_migration_tk.batching import (CLAIM_LEASE_SECONDS,
                                        DEFAULT_BATCH_SIZE, HEARTBEAT_SECONDS,
                                        RUNNING_CANARY, STATE_SUCCESS,
                                        SUCCESS_CANARY, STAT_THREADS,
                                        BatchLayout, Heartbeat, age_phrase,
                                        batch_log, batch_state, batchfile_path,
                                        claim_age, claim_batch, concatenate,
                                        fail_batch, finish_batch, plan_batches,
                                        read_batchfile, read_canary,
                                        release_claim, split_by_fasta,
                                        tally_reasons, write_table)
from gtdb_migration_tk.biolib_lite.checksum import sha256, sha256_rb
from gtdb_migration_tk.biolib_lite.common import make_sure_path_exists
from gtdb_migration_tk.biolib_lite.external.execute import check_dependencies
from gtdb_migration_tk.biolib_lite.external.pfam_search import PfamSearch
from gtdb_migration_tk.update_genomes import genomes_to_regenerate
from gtdb_migration_tk.utils.common import PROTEIN_FASTA_EXT, protein_fasta
from gtdb_migration_tk.utils.tools import symlink, openfile


# What this command calls the files of its own batches. The batchfile is the
# record of which genomes a batch is, and its first column is each genome's
# protein FASTA -- the file this command reads, where trans_table and prodigal
# read the genomic FASTA. --out_dir holds nothing else of the run's results: the
# marker tables go into the genome directories as they always have, which is why
# two machines on different batches never write to the same place. There is no
# older name to look for: this command has never had batches before.
BATCHFILE_NAME = 'hmmsearch_batchfile.tsv.gz'
BATCH_LOG_NAME = 'hmmsearch.log'
LAYOUT = BatchLayout(batchfiles=(BATCHFILE_NAME,), log=BATCH_LOG_NAME)

# The batches of a run live under <out_dir>/<marker directory>/, e.g.
# <out_dir>/pfam_33.1_lite/batch_000001/. One run of this command searches ONE
# database at ONE version, and its batches are finished when that database has
# been searched; a run of the other database over the same release has different
# work to do for the same genomes. Sharing a directory would have the SUCCESS of
# a Pfam batch tell a TIGRFAM run that batch was done. The suffix is in the name
# for the same reason: annotating against a new Pfam release is new work.

# Genomes of a batch that came out of it with no marker table, and the same
# gathered for the release. A genome with no called proteins cannot be searched,
# and the step after this one should not have to look in 135 batch directories to
# find out which genomes those were.
NOT_SEARCHED_NAME = 'not_searched.tsv'
NOT_SEARCHED_RELEASE_NAME = 'hmmsearch_not_searched.tsv'
NOT_SEARCHED_HEADER = ('genome_id', 'reason')
REASON_NO_PROTEINS = 'no_protein_file'

# One genome as marker_parser() receives it, the tuple being what mp.Pool.imap_unordered
# can carry: accession, its protein FASTA as the batchfile names it, the marker
# directory within prodigal/, the extension the marker table carries, the genomes the
# release says need annotating, the database's name for the log, and whether the run was
# told to search every genome again.
MarkerJob = Tuple[str, str, str, str, Set[str], str, bool]

# Gene ID to the hits kept for it: HMM ID to (e-value, bitscore) for Pfam, where a gene
# keeps its best hit per family, and a single (HMM ID, e-value, bitscore) for TIGRFAM,
# where it keeps one hit overall.
PfamTopHits = Dict[str, Dict[str, Tuple[float, float]]]
TigrTopHits = Dict[str, Tuple[str, float, float]]

# What --hmm_db_path is checked to be before a batch is claimed. The two values of
# --db want different things of it and neither says so when handed the other:
# PfamScan is given a DIRECTORY and looks for Pfam-A.hmm inside it, hmmsearch is
# given the TIGRFAM HMM FILE itself. Given the directory, hmmsearch reads it as a
# file that "appears to be empty", exits non-zero and writes no marker table, and
# the run met that only one call later, in a worker, as a FileNotFoundError on the
# table -- naming the missing file rather than the reason for it. In r237 that cost
# seven hours on five machines and 303,096 genomes, every batch of them marked
# SUCCESS. The first five bytes of an HMM library are its format line, so being
# sure costs a stat and a read.
HMM_MAGIC = b'HMMER'
PFAM_LIBRARY = 'Pfam-A.hmm'


class BadHmmDatabase(ValueError):
    """--hmm_db_path is not the HMMs the database asked for on --db is searched from."""


class BatchCounts(NamedTuple):
    """What searching the markers of one batch came to.

    Recorded in the batch's SUCCESS canary, because the release totals are added
    up from the batches and a machine that ran the last batch has searched the
    markers of none of the others.
    """

    searched: int
    already_searched: int
    not_searched: int


class MarkerSetup(NamedTuple):
    """What searching one marker database needs, settled once for the whole run.

    Which database is searched decides four things at once -- where the results
    go, what they are called, which worker runs and what the log calls it -- and
    they are settled together rather than at each place one of them is wanted.
    """

    marker_dir: str
    extension: str
    name: str
    worker: Callable


def check_hmm_file(path: str) -> None:
    """Refuse a file that is not an HMM library.

    Parameters
    ----------
    path : str
        File hmmsearch or hmmscan would be pointed at.

    @return: None; BadHmmDatabase where the file is not one hmmsearch can read.
    """

    try:
        with open(path, 'rb') as handle:
            magic = handle.read(len(HMM_MAGIC))
    except OSError as exc:
        raise BadHmmDatabase('{} cannot be read: {}'.format(path, exc))

    if magic != HMM_MAGIC:
        raise BadHmmDatabase(
            '{} does not begin with {}, so it is not an HMM library hmmsearch '
            'can read.'.format(path, HMM_MAGIC.decode()))


def check_hmm_db(db: str, hmm_db_path: str) -> None:
    """Refuse --hmm_db_path now where the search would fail on every genome.

    Called before the release is cut into batches, so a path the search cannot
    use costs a second at startup rather than a run that claims batches for hours
    and finishes them with nothing in them.

    Parameters
    ----------
    db : str
        'pfam' or 'tigrfam'.
    hmm_db_path : str
        --hmm_db_path: the directory holding Pfam-A.hmm for 'pfam', the HMM file
        itself for 'tigrfam'.

    @return: None; BadHmmDatabase naming what was wanted where it is not that.
    """

    if db == 'pfam':
        if not os.path.isdir(hmm_db_path):
            raise BadHmmDatabase(
                "--db pfam searches a DIRECTORY holding {}, and --hmm_db_path "
                "{} is not a directory.".format(PFAM_LIBRARY, hmm_db_path))

        library = os.path.join(hmm_db_path, PFAM_LIBRARY)
        if not os.path.isfile(library):
            raise BadHmmDatabase(
                "--db pfam looks for {} inside --hmm_db_path, and there is no "
                "{}.".format(PFAM_LIBRARY, library))

        check_hmm_file(library)
        return

    if os.path.isdir(hmm_db_path):
        # the r237 mistake, and the one the two --db values invite: the directory
        # that IS --hmm_db_path for Pfam holds the file that is --hmm_db_path here
        suggestion = os.path.join(hmm_db_path, 'tigrfam.hmm')
        raise BadHmmDatabase(
            "--db tigrfam searches ONE HMM FILE rather than a directory of them, "
            "and --hmm_db_path {} is a directory.{}".format(
                hmm_db_path,
                ' Did you mean {}?'.format(suggestion)
                if os.path.isfile(suggestion) else ''))

    if not os.path.isfile(hmm_db_path):
        raise BadHmmDatabase(
            "--db tigrfam searches the HMM file --hmm_db_path names, and there is "
            "no {}.".format(hmm_db_path))

    check_hmm_file(hmm_db_path)


class MarkerManager(object):
    """Identify marker genes using Pfam and tigrfam HMMs."""

    def __init__(self,
                 tmp_dir: str = '/tmp/',
                 cpus: int = 1,
                 batch_size: int = DEFAULT_BATCH_SIZE,
                 reclaim: bool = False,
                 lease: float = CLAIM_LEASE_SECONDS,
                 heartbeat: float = HEARTBEAT_SECONDS) -> None:
        """Initialization.

        Parameters
        ----------
        tmp_dir : str
            Directory each genome's proteins are decompressed into before they
            are searched, one directory per genome in flight and removed after.
            No results are written here.
        cpus : int
            How many genomes are annotated at once, and so how many decompressed
            proteomes sit in tmp_dir at once.
        batch_size : int
            Genomes per batch.
        reclaim : bool
            Take over a batch another machine holds before its claim has expired.
        lease : float
            Seconds a claim survives without the machine holding it saying so.
        heartbeat : float
            Seconds between this machine saying so about a batch of its own.

        @return: None
        """

        self.tmp_dir: str = tmp_dir
        self.cpus: int = cpus
        self.batch_size: int = batch_size
        self.reclaim: bool = reclaim
        self.lease: float = lease
        self.heartbeat: float = heartbeat

        check_dependencies(['prodigal', 'hmmsearch'])

        # made here rather than by the first worker that wants it: a --tmp_dir
        # that cannot be made would otherwise be met once per genome, inside a
        # batch already claimed, and would fail every batch this machine took
        make_sure_path_exists(self.tmp_dir)

        # identify TIGRfam and Pfam marker genes comprising the bac120, ar122, ar53, or
        # rp2 marker sets using a carefully selected subset of HMMs. Which of the two is
        # searched is decided per run, so only one of these is ever set.
        self.tigrfam_hmms: str = ''
        self.pfam_hmm_dir: str = ''

        self.protein_file_ext: str = PROTEIN_FASTA_EXT

        self.logger: logging.Logger = logging.getLogger('timestamp')

    def marker_setup(self, db: str, dir_suffix: str, hmm_db_path: str) -> MarkerSetup:
        """Settle everything that follows from which marker database is searched.

        The HMMs are put on the instance rather than passed down, because the
        workers are processes forked from it and read them from there.

        Parameters
        ----------
        db : str
            'pfam' or 'tigrfam'.
        dir_suffix : str
            Suffix of the marker directory and files, e.g. 33.1_lite.
        hmm_db_path : str
            The HMMs to search against.

        @return: where the results go, what they are called, what the log calls
                 the database, and the worker that searches it.
        """

        if db in ('pfam', 'tigrfam'):
            # here rather than in either worker: run_hmmsearch() settles the
            # database before it plans the batches, so HMMs the search cannot use
            # are met once, at the start, by the machine that was mistyped at --
            # not ninety-six times a batch by workers that die one after another
            check_hmm_db(db, hmm_db_path)

        if db == 'pfam':
            self.pfam_hmm_dir = hmm_db_path
            return MarkerSetup(marker_dir='pfam_{}'.format(dir_suffix),
                               extension='_pfam_{}.tsv'.format(dir_suffix),
                               name='Pfam',
                               worker=self.__pfam_worker)

        if db == 'tigrfam':
            self.tigrfam_hmms = hmm_db_path
            return MarkerSetup(marker_dir='tigrfam_{}'.format(dir_suffix),
                               extension='_tigrfam_{}.tsv'.format(dir_suffix),
                               name='Tigrfam',
                               worker=self.__tigrfam_worker)

        # argparse limits --db to the two, so this is a caller that went around it
        # rather than a user; it used to leave marker_dir unbound and fail later
        # with a NameError naming nothing
        raise ValueError("--db is 'pfam' or 'tigrfam', not {!r}.".format(db))

    def run_hmmsearch(self,
                      gtdb_genome_path_file: str,
                      report: str,
                      db: str,
                      dir_suffix: str,
                      hmm_db_path: str,
                      out_dir: str,
                      all_genomes: bool = False) -> bool:
        """Identify marker genes using Pfam and TIGRfam HMMs.

        The release is cut into batches under --out_dir and a batch is claimed
        before it is worked on, so several machines can be pointed at one
        --out_dir and will divide the release between them. The batching is the
        same machinery trans_table and prodigal use, in batching.py. What
        --out_dir holds is the state of the run and nothing else: the marker
        tables go into each genome's own prodigal/ directory, as they always
        have, which is why two machines on different batches never write to the
        same place.

        The batches of a run are held under <out_dir>/<marker directory>/, so one
        --out_dir carries a Pfam run and a TIGRFAM run of the same release
        without either reading the other's canaries. They are different work over
        the same genomes, and a batch finished for one is not finished for the
        other.

        Whether a genome is searched is decided by the marker table it already
        has, not by the report: a genome is skipped where its table is there AND
        its checksum agrees. The report says which genomes the release did not
        carry derived data for, and a genome that disagrees with it -- annotated
        already though the release calls it new, or unannotated though it does
        not -- is searched and said so in the log, since the disagreement is
        worth seeing and neither answer is worth withholding the work over.

        Parameters
        ----------
        gtdb_genome_path_file : str
            genome_dirs file of the release, accession and genome directory per line.
        report : str
            report.log of the release, read for the genomes needing annotation.
        db : str
            'pfam' or 'tigrfam'.
        dir_suffix : str
            Suffix of the marker directory and files, e.g. 33.1_lite.
        hmm_db_path : str
            The HMMs to search against.
        out_dir : str
            Directory the batches and the state of the run are written to.
        all_genomes : bool
            Search every genome again, discarding the marker tables that are there.

        @return: True where every batch this machine took finished, False where
                 one failed and is left to a later run.
        """

        setup = self.marker_setup(db, dir_suffix, hmm_db_path)

        genomes_to_consider = genomes_to_regenerate(report)
        self.logger.info(
            '{:,} genome(s) of the release did not carry their derived data '
            'across.'.format(len(genomes_to_consider)))

        # the batches of THIS database at THIS version; another database's run
        # over the same release has its own, beside these
        state_dir = os.path.join(out_dir, setup.marker_dir)
        batches = plan_batches(gtdb_genome_path_file, state_dir, self.batch_size,
                               LAYOUT, self.logger, genome_file=protein_fasta)

        if all_genomes:
            self.logger.warning(
                'warning: --all discards the marker table of every genome and '
                'searches it again, so batches already finished are done again '
                'too. Without it a finished batch is skipped.')

        done, held, failed = 0, 0, 0
        for index, batch_dir in enumerate(batches, start=1):
            label = 'Batch {:,} of {:,} ({})'.format(
                index, len(batches), os.path.basename(batch_dir))

            if not all_genomes and batch_state(batch_dir) == STATE_SUCCESS:
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
            with batch_log(batch_dir, self.logger, LAYOUT):
                self.logger.info('{}: starting.'.format(label))
                try:
                    with Heartbeat(os.path.join(batch_dir, RUNNING_CANARY),
                                   self.heartbeat):
                        counts = self.search_batch(batch_dir, setup, dir_suffix,
                                                   genomes_to_consider, all_genomes)
                except KeyboardInterrupt:
                    # the machine holding it is stopping, so the batch is handed
                    # back rather than left to sit out its lease
                    release_claim(batch_dir)
                    self.logger.error('{}: interrupted; the claim is given up and '
                                      'the batch carries on where it stopped.'.format(label))
                    raise
                except Exception as exc:
                    failed += 1
                    fail_batch(batch_dir, str(exc))
                    self.logger.error('{}: failed and will be retried by a later '
                                      'run: {}'.format(label, exc))
                    continue

                finish_batch(batch_dir,
                             searched=counts.searched,
                             already_searched=counts.already_searched,
                             not_searched=counts.not_searched)
                done += 1
                self.logger.info('{}: done.'.format(label))

        self.logger.info(
            '{:,} batch(es) finished here, {:,} held by another machine, '
            '{:,} failed.'.format(done, held, failed))

        self.aggregate(batches, state_dir, setup.name)

        # a batch that failed has already said why, in its own log and in its
        # FAILED file, and a run of several machines over days should not end by
        # printing a traceback out of the argparse frames
        if failed:
            self.logger.error(
                '{:,} batch(es) failed; they are the directories holding a FAILED '
                'file and are retried by running the command again.'.format(failed))
            return False

        return True

    def search_batch(self,
                     batch_dir: str,
                     setup: MarkerSetup,
                     dir_suffix: str,
                     genomes_to_consider: Set[str],
                     all_genomes: bool = False) -> BatchCounts:
        """Search the markers of one batch, and record the genomes that got none.

        The genomes are taken from the batch's own batchfile, so the work asks
        about the genomes the batch was cut from rather than about whatever a
        genome_dirs file says now.

        Two passes, because deciding is cheap and searching is not: the first
        sorts the batch into genomes whose marker tables are already there and
        vouched for and genomes that are not, the second searches what is left.
        Doing it per batch rather than over the release means a rerun re-reads the
        marker tables of one batch at a time instead of all 1.35M before it starts.

        Parameters
        ----------
        batch_dir : str
            Batch directory.
        setup : MarkerSetup
            Which database is searched and where its results go.
        dir_suffix : str
            Suffix of the marker directory and files, which the workers rebuild
            their own filenames from.
        genomes_to_consider : set
            Genomes the release did not carry derived data for.
        all_genomes : bool
            Discard the marker table of every genome of the batch and search again.

        @return: what the batch came to, which run_hmmsearch() records in its canary.
        """

        rows = read_batchfile(batchfile_path(batch_dir, LAYOUT))

        # a genome whose proteins are not there has nothing to search; it is named
        # and left rather than stopping the other ten thousand of the batch
        present, missing = split_by_fasta(rows, STAT_THREADS)
        not_searched = [(accession, REASON_NO_PROTEINS) for accession in missing]

        jobs = [(accession, proteins, setup.marker_dir, setup.extension,
                 genomes_to_consider, setup.name, all_genomes)
                for proteins, accession in present]
        with mp.Pool(processes=self.cpus) as pool:
            decided = list(tqdm(pool.imap_unordered(self.marker_parser, jobs),
                                total=len(jobs), unit='genome', ncols=100,
                                leave=False, desc='Checking marker tables'))

        # every skipped genome is None and every None has to go: the queue the
        # workers draw from ends with one per worker as the signal to stop, and a
        # worker cannot tell a genome that was skipped from the end of the work
        genome_files = [proteins for proteins in decided if proteins is not None]
        already_searched = len(jobs) - len(genome_files)
        self.logger.info(
            '{:,} genome(s) require {} annotation; {:,} already have valid '
            'results.'.format(len(genome_files), setup.name, already_searched))

        self.search_markers(genome_files, setup.worker, dir_suffix)

        not_searched.sort()
        write_table(not_searched, os.path.join(batch_dir, NOT_SEARCHED_NAME),
                    header=NOT_SEARCHED_HEADER)

        if not_searched:
            self.logger.warning(
                'warning: {:,} genome(s) of this batch have no marker table: '
                '{}.'.format(len(not_searched),
                             '; '.join('{:,} {}'.format(count, reason)
                                       for reason, count
                                       in sorted(tally_reasons(not_searched).items()))))

        return BatchCounts(searched=len(genome_files),
                           already_searched=already_searched,
                           not_searched=len(not_searched))

    def search_markers(self,
                       genome_files: Sequence[str],
                       worker: Callable,
                       dir_suffix: str) -> None:
        """Run the HMMs over the genomes that need them, on self.cpus processes.

        The queue ends with one None per worker, which is how each is told there
        is no more work; nothing else in it may be None.

        Parameters
        ----------
        genome_files : sequence of str
            Protein FASTA of each genome to search.
        worker : callable
            __pfam_worker or __tigrfam_worker.
        dir_suffix : str
            Suffix the worker rebuilds its own filenames from.

        @return: None
        """

        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        for f in genome_files:
            workerQueue.put(f)

        for _ in range(self.cpus):
            workerQueue.put(None)

        workerProc = [mp.Process(target=worker, args=(
            workerQueue, writerQueue, dir_suffix)) for _ in range(self.cpus)]
        writeProc = mp.Process(target=self.__progress, args=(
            len(genome_files), writerQueue))

        try:
            writeProc.start()

            for p in workerProc:
                p.start()

            for p in workerProc:
                p.join()

            # how a worker ended was never looked at. A worker that raises is
            # gone, and the genomes still in the queue go unsearched with it, but
            # the batch was finished all the same and its canary said
            # searched=10,000 -- the length of the work list, which is what was
            # handed out and not what came back. In r237 every worker of every
            # batch died on the first genome, on HMMs hmmsearch could not read,
            # and the run marked all 135 batches SUCCESS.
            died = [p.exitcode for p in workerProc if p.exitcode]

            writerQueue.put(None)
            writeProc.join()

            if died:
                raise RuntimeError(
                    '{:,} of {:,} search process(es) ended in error, the first '
                    'with exit status {}; the genomes they held have no marker '
                    'table, and the batch is failed rather than finished. What '
                    'the search said is above this in the log.'.format(
                        len(died), len(workerProc), died[0]))
        except BaseException:
            # raised on, rather than swallowed: a batch whose search died has not
            # searched its genomes, and returning quietly here would have
            # run_hmmsearch() write it a SUCCESS canary that no other machine
            # would ever look behind
            for p in workerProc:
                p.terminate()
            writeProc.terminate()
            raise
        finally:
            # a worker that died left the work it never took in the queue, and a
            # queue still holding data holds the process open: the feeder thread
            # is joined at exit and waits forever for a reader that has gone. The
            # r237 run sat there for hours after its last batch, with every
            # worker dead, its log finished and nothing left to do.
            workerQueue.cancel_join_thread()
            writerQueue.cancel_join_thread()

    def aggregate(self, batches: Sequence[str], state_dir: str, name: str) -> None:
        """Report the release, once every batch of this database has succeeded.

        Written only when they all have, so that the file at the top of the
        directory is either the whole release or absent, and never a part of it
        that reads like the whole.

        Parameters
        ----------
        batches : sequence of str
            Every batch directory of the run, in batch order.
        state_dir : str
            Where this database's batches are held, under --out_dir.
        name : str
            What the log calls the database.

        @return: None
        """

        unfinished = [batch for batch in batches
                      if batch_state(batch) != STATE_SUCCESS]
        if unfinished:
            self.logger.info(
                '{:,} of {:,} batch(es) are done; the release files are '
                'written once they all are.'.format(
                    len(batches) - len(unfinished), len(batches)))
            return

        path = os.path.join(state_dir, NOT_SEARCHED_RELEASE_NAME)
        written = concatenate(
            [os.path.join(batch, NOT_SEARCHED_NAME) for batch in batches], path)

        totals = {'searched': 0, 'already_searched': 0}
        for batch_dir in batches:
            canary = read_canary(os.path.join(batch_dir, SUCCESS_CANARY))
            for field in totals:
                try:
                    totals[field] += int(canary[field])
                except (KeyError, ValueError):
                    pass

        # of the release and not of this machine: the counts are added up from
        # every batch's canary, and the batches were shared out
        self.logger.info(
            'Release: {:,} genome(s) were searched against {}, {:,} already had '
            'valid results.'.format(totals['searched'], name,
                                    totals['already_searched']))

        if written:
            self.logger.warning(
                'warning: {:,} genome(s) of the release have no marker table and '
                'are named in {}.'.format(written, path))
        else:
            self.logger.info(
                'Every genome of the release has a marker table; wrote {} with no '
                'rows.'.format(path))

    def marker_parser(self, job: MarkerJob) -> Optional[str]:
        """Decide whether one genome's markers still have to be searched for.

        Parameters
        ----------
        job : MarkerJob
            One genome, as search_batch() packed it.

        @return: the protein file to search, or None to skip the genome, it
                 being annotated already. One value for skipping, because
                 search_batch() does the same thing with every genome it skips
                 and a second one has only ever been a way of missing one of them.
        """

        gid, gene_file, marker_dir, full_extension, genomes_to_consider, name, all_genomes = job

        if all_genomes:
            return gene_file

        prodigal_dir = os.path.dirname(gene_file)
        marker_file = os.path.join(prodigal_dir, marker_dir, gid + full_extension)
        marker_zipped_file = marker_file + '.gz'
        if os.path.exists(marker_zipped_file):
            # verify checksum
            checksum_file = marker_file + '.sha256'
            if os.path.exists(checksum_file):
                with open(marker_zipped_file, 'rb') as raw:
                    checksum = sha256_rb(gzip.GzipFile(fileobj=raw))
                with open(checksum_file) as handle:
                    cur_checksum = handle.readline().strip()
                if checksum == cur_checksum:
                    if gid in genomes_to_consider:
                        self.logger.warning(
                            f'Genome {gid} is marked as new or modified, but already has {name} annotations.')
                        self.logger.warning('Genome is being skipped!')
                    return None

            # reached when the .sha256 is absent as well as when it disagrees:
            # neither says the annotations are wrong, but neither shows them to be
            # right. Whether the release expected them is said truthfully rather
            # than assumed -- a table that cannot be vouched for usually belongs
            # to a genome the release DOES call new, its run having been
            # interrupted partway through writing it
            if gid in genomes_to_consider:
                self.logger.warning(
                    f'Genome {gid} has {name} annotations with no valid checksum, and is marked as new or modified.')
            else:
                self.logger.warning(
                    f'Genome {gid} has {name} annotations with no valid checksum, though it is not marked for reannotation.')
            self.logger.warning(f'Genome will be reannotated.')

        elif gid not in genomes_to_consider:
            self.logger.warning(
                f'Genome {gid} has no {name} annotations, but is also not marked for processing?')
            self.logger.warning(f'Genome will be reannotated!')

        return gene_file

    def run_tophit(self, gtdb_genome_path_file: str, db: str, folder_name: str) -> None:
        """Reduce each genome's marker table to its top hits.

        Parameters
        ----------
        gtdb_genome_path_file : str
            genome_dirs file of the release: accession, path, canonical accession.
        db : str
            'pfam' or 'tigrfam'.
        folder_name : str
            Suffix of the marker directory and files, e.g. 33.1_lite.

        @return: nothing; a tophit file is written beside each marker table.
        """

        extension = ""
        if db == 'pfam':
            marker_version = 'pfam_{}'.format(folder_name)
            extension = f'_{marker_version}.tsv.gz'
            tophit_out = f'_{marker_version}_tophit.tsv'
        elif db == 'tigrfam':
            marker_version = 'tigrfam_{}'.format(folder_name)
            extension = f'_{marker_version}.tsv.gz'
            tophit_out = f'_{marker_version}_tophit.tsv'

        countr = 0
        for line in open(gtdb_genome_path_file):
            countr += 1
            statusStr = '{} lines read.'.format(countr)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

            line_split = line.strip().split('\t')

            gid = line_split[0]
            gpath = line_split[1]

            prodigal_dir = os.path.join(gpath, 'prodigal')

            gene_file = os.path.join(
                prodigal_dir, gid + self.protein_file_ext)
            if os.path.exists(gene_file):
                if os.stat(gene_file).st_size == 0:
                    self.logger.warning(
                        f' Protein file appears to be empty: {gene_file}')
                else:
                    assembly_dir, filename = os.path.split(gene_file)

                    output_hit_file = os.path.join(assembly_dir, marker_version, filename.replace(
                        self.protein_file_ext, extension))
                    # determine top hits
                    tophit_file = os.path.join(assembly_dir, marker_version, filename.replace(
                        self.protein_file_ext, tophit_out))
                    if not os.path.exists(output_hit_file):
                        self.logger.warning(
                            f'Output file does not exist: {output_hit_file}')
                        continue
                    if db == 'pfam':
                        self._pfam_top_hit(output_hit_file, tophit_file)
                    elif db == 'tigrfam':
                        self._tigr_top_hit(output_hit_file, tophit_file)

                    with open(tophit_file, 'rb') as f_in, gzip.open(tophit_file + '.gz', 'wb') as f_out:
                        f_out.writelines(f_in)
                    os.remove(tophit_file)

    def __progress(self, num_items: int, queue_out: Queue) -> None:
        """Store or write results of worker threads in a single thread."""
        processed_items = 0
        while True:
            a = queue_out.get(block=True, timeout=None)
            if a == None:
                break

            processed_items += 1
            statusStr = 'Finished processing %d of %d (%.2f%%) items.' % (
                processed_items, num_items, float(processed_items) * 100 / num_items)
            sys.stdout.write('%s\r' % statusStr)
        sys.stdout.flush()

        sys.stdout.write('\n')

    def __pfam_worker(self, queue_in: Queue, queue_out: Queue, folder_name: str) -> None:
        """Process each data item in parallel."""

        prefix = "pfam"
        pfam_version = '{}_{}'.format(prefix,folder_name)
        pfam_extension = f'_{pfam_version}.tsv'
        pfam_extension_gz = f'_{pfam_version}.tsv.gz'
        pfam_tophit_extension = f'_{pfam_version}_tophit.tsv'
        pfam_tophit_extension_gz = f'_{pfam_version}_tophit.tsv.gz'

        if '_lite' in pfam_extension:
            symlink_pfam_extension_gz = f'_{prefix}_lite.tsv.gz'
            symlink_pfam_tophit_extension_gz = f'_{prefix}_lite_tophit.tsv.gz'
        else:
            symlink_pfam_extension_gz = f'_{prefix}.tsv.gz'
            symlink_pfam_tophit_extension_gz = f'_{prefix}_tophit.tsv.gz'

        while True:
            gene_file = queue_in.get(block=True, timeout=None)
            if gene_file == None:
                break

            assembly_dir, filename = os.path.split(gene_file)
            make_sure_path_exists(os.path.join(assembly_dir, pfam_version))

            output_hit_file = os.path.join(
                assembly_dir, pfam_version, filename.replace(self.protein_file_ext, pfam_extension))
            #because the gene file is a zipped file, we need to unzip it in a temporary directory
            temp_dir = tempfile.mkdtemp(dir=self.tmp_dir)
            try:
                temp_gene_file = os.path.join(temp_dir, filename[0:-3])

                # if size of temp_gene_file is 0, then skip hmmsearch
                if os.stat(gene_file).st_size == 0:
                    self.logger.warning('Skipping %s because it is empty' % temp_gene_file)
                    continue

                with gzip.open(gene_file, 'rb') as f_in:
                    with open(temp_gene_file, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)

                pfam_search = PfamSearch(self.pfam_hmm_dir)
                pfam_search.run(temp_gene_file, output_hit_file)

                # determine top hits
                pfam_tophit_file = os.path.join(assembly_dir, pfam_version, filename.replace(
                    self.protein_file_ext, pfam_tophit_extension))
                self._pfam_top_hit(output_hit_file, pfam_tophit_file)



                # calculate checksum
                checksum = sha256(output_hit_file)
                fout = open(output_hit_file + '.sha256', 'w')
                fout.write(checksum)
                fout.close()

                # archive the pfam file and the tophit file
                with open(output_hit_file, 'rb') as f_in, gzip.open(output_hit_file + '.gz', 'wb') as f_out:
                    f_out.writelines(f_in)
                os.remove(output_hit_file)
                with open(pfam_tophit_file, 'rb') as f_in, gzip.open(pfam_tophit_file + '.gz', 'wb') as f_out:
                    f_out.writelines(f_in)
                os.remove(pfam_tophit_file)


                # create symlink in prodigal_folder
                new_hit_link = os.path.join(assembly_dir, filename.replace(
                    self.protein_file_ext, symlink_pfam_extension_gz))
                new_tophit_link = os.path.join(assembly_dir, filename.replace(
                    self.protein_file_ext, symlink_pfam_tophit_extension_gz))

                output_hit_file_relative = os.path.join(
                    '.', pfam_version, filename.replace(self.protein_file_ext, pfam_extension_gz))
                pfam_tophit_file_relative = os.path.join('.', pfam_version, filename.replace(
                    self.protein_file_ext, pfam_tophit_extension_gz))

                symlink(output_hit_file_relative, new_hit_link, overwrite=True)
                symlink(pfam_tophit_file_relative, new_tophit_link, overwrite=True)


            #we can now delete the temporary directory
            finally:
                shutil.rmtree(temp_dir)

            # allow results to be processed or written to file
            queue_out.put(gene_file)

    def _pfam_top_hit(self, pfam_file: str, pfam_tophit_file: str) -> None:
        """Identify top Pfam hits.

        Parameters
        ----------
        pfam_file : str
            Marker table written by the Pfam search.
        pfam_tophit_file : str
            Where the top hits are written, with a .sha256 beside them.

        @return: nothing; a gene keeps its best hit for each family it matched.
        """

        tophits: PfamTopHits = defaultdict(dict)
        for line in openfile(pfam_file):
            if line[0] == '#' or not line.strip():
                continue

            line_split = line.split()
            gene_id = line_split[0]
            hmm_id = line_split[5]
            evalue = float(line_split[12])
            bitscore = float(line_split[11])
            if gene_id in tophits:
                if hmm_id in tophits[gene_id]:
                    if bitscore > tophits[gene_id][hmm_id][1]:
                        tophits[gene_id][hmm_id] = (evalue, bitscore)
                else:
                    tophits[gene_id][hmm_id] = (evalue, bitscore)
            else:
                tophits[gene_id][hmm_id] = (evalue, bitscore)

        fout = open(pfam_tophit_file, 'w')
        fout.write('Gene Id\tTop hits (Family id,e-value,bitscore)\n')
        for gene_id, hits in tophits.items():
            hit_str = []
            for hmm_id, stats in hits.items():
                hit_str.append(hmm_id + ',' + ','.join(map(str, stats)))
            fout.write('%s\t%s\n' % (gene_id, ';'.join(hit_str)))
        fout.close()

        # calculate checksum
        checksum = sha256(pfam_tophit_file)
        fout = open(pfam_tophit_file + '.sha256', 'w')
        fout.write(checksum)
        fout.close()

    def _tigr_top_hit(self, tigrfam_file: str, tigrfam_tophit_file: str) -> None:
        """Identify top TIGRfam hits.

        Parameters
        ----------
        tigrfam_file : str
            Marker table written by the TIGRFAM search.
        tigrfam_tophit_file : str
            Where the top hits are written, with a .sha256 beside them.

        @return: nothing; a gene keeps one hit, the highest scoring of any family.
        """

        tophits: TigrTopHits = {}
        for line in openfile(tigrfam_file):
            if line[0] == '#' or line[0] == '[':
                continue

            line_split = line.split()
            gene_id = line_split[0]
            hmm_id = line_split[3]
            evalue = float(line_split[4])
            bitscore = float(line_split[5])
            if gene_id in tophits:
                if bitscore > tophits[gene_id][2]:
                    tophits[gene_id] = (hmm_id, evalue, bitscore)
            else:
                tophits[gene_id] = (hmm_id, evalue, bitscore)

        fout = open(tigrfam_tophit_file, 'w')
        fout.write('Gene Id\tTop hits (Family id,e-value,bitscore)\n')
        for gene_id, stats in tophits.items():
            hit_str = ','.join(map(str, stats))
            fout.write('%s\t%s\n' % (gene_id, hit_str))
        fout.close()

        # calculate checksum
        checksum = sha256(tigrfam_tophit_file)
        fout = open(tigrfam_tophit_file + '.sha256', 'w')
        fout.write(checksum)
        fout.close()

    def __tigrfam_worker(self, queue_in: Queue, queue_out: Queue, folder_name: str) -> None:
        """Process each data item in parallel."""

        prefix = "tigrfam"
        tigrfam_version = f'{prefix}_{folder_name}'
        tigrfam_extension = f'_{tigrfam_version}.tsv'
        tigrfam_extension_gz = f'_{tigrfam_version}.tsv.gz'
        tigrfam_out = f'_{tigrfam_version}.out'
        tigrfam_out_gz = f'_{tigrfam_version}.out.gz'
        tigrfam_tophit_extension = f'_{tigrfam_version}_tophit.tsv'
        tigrfam_tophit_extension_gz = f'_{tigrfam_version}_tophit.tsv.gz'


        if '_lite' in tigrfam_extension:
            symlink_tigrfam_extension_gz = f'_{prefix}_lite.tsv.gz'
            symlink_tigrfam_tophit_extension_gz = f'_{prefix}_lite_tophit.tsv.gz'
            symlink_tigrfam_out_gz = f'_{prefix}_lite.out.gz'
        else:
            symlink_tigrfam_extension_gz = f'_{prefix}.tsv.gz'
            symlink_tigrfam_tophit_extension_gz = f'_{prefix}_tophit.tsv.gz'
            symlink_tigrfam_out_gz = f'_{prefix}.out.gz'

        while True:
            gene_file = queue_in.get(block=True, timeout=None)
            if gene_file == None:
                break

            assembly_dir, filename = os.path.split(gene_file)
            make_sure_path_exists(os.path.join(assembly_dir, tigrfam_version))

            output_hit_file = os.path.join(assembly_dir, tigrfam_version, filename.replace(
                self.protein_file_ext, tigrfam_extension))
            out_file = os.path.join(assembly_dir, tigrfam_version, filename.replace(
                self.protein_file_ext, tigrfam_out))

            #because the gene file is a zipped file, we need to unzip it in a temporary directory
            temp_dir = tempfile.mkdtemp(dir=self.tmp_dir)
            try:
                temp_gene_file = os.path.join(temp_dir, filename[0:-3])
                with gzip.open(gene_file, 'rb') as f_in:
                    with open(temp_gene_file, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)

                # if size of temp_gene_file is 0, then skip hmmsearch
                if os.stat(temp_gene_file).st_size == 0:
                    self.logger.warning('Skipping %s because it is empty' % temp_gene_file)
                    continue

                hmmsearch_out = os.path.join(assembly_dir, tigrfam_version, filename.replace(
                    self.protein_file_ext, f'_{tigrfam_version}.out'))
                cmd = ['hmmsearch', '-o', hmmsearch_out, '--tblout', output_hit_file,
                       '--noali', '--notextw', '--cut_nc', '--cpu', '1',
                       self.tigrfam_hmms, temp_gene_file]

                # os.system() threw the exit status away, so a search that failed
                # was met one call later, as a FileNotFoundError on the marker
                # table it never wrote: the traceback named the missing table and
                # not the reason there was none. hmmsearch says the reason on
                # stderr, so the reason is what is raised.
                search = subprocess.run(cmd, stderr=subprocess.PIPE,
                                        universal_newlines=True)
                if search.returncode != 0:
                    raise RuntimeError(
                        'hmmsearch exited {} searching {} against {}: {}'.format(
                            search.returncode, gene_file, self.tigrfam_hmms,
                            ' '.join(search.stderr.split()) or 'it said nothing'))

                # determine top hits
                tigrfam_tophit_file = os.path.join(assembly_dir, tigrfam_version, filename.replace(
                    self.protein_file_ext, tigrfam_tophit_extension))
                self._tigr_top_hit(output_hit_file, tigrfam_tophit_file)

                # calculate checksum
                checksum = sha256(output_hit_file)
                fout = open(output_hit_file + '.sha256', 'w')
                fout.write(checksum)
                fout.close()

                # archive the pfam file and the tophit file
                with open(output_hit_file, 'rb') as f_in, gzip.open(output_hit_file + '.gz', 'wb') as f_out:
                    f_out.writelines(f_in)
                os.remove(output_hit_file)
                with open(out_file, 'rb') as f_in, gzip.open(out_file + '.gz', 'wb') as f_out:
                    f_out.writelines(f_in)
                os.remove(out_file)
                with open(tigrfam_tophit_file, 'rb') as f_in, gzip.open(tigrfam_tophit_file + '.gz', 'wb') as f_out:
                    f_out.writelines(f_in)
                os.remove(tigrfam_tophit_file)

                # create symlink in prodigal_folder
                new_hit_link = os.path.join(assembly_dir, filename.replace(
                    self.protein_file_ext, symlink_tigrfam_extension_gz))
                new_tophit_link = os.path.join(assembly_dir, filename.replace(
                    self.protein_file_ext, symlink_tigrfam_tophit_extension_gz))
                new_out_link = os.path.join(assembly_dir, filename.replace(
                    self.protein_file_ext, symlink_tigrfam_out_gz))

                # Symlink needs to be relative to avoid pointing to previous version of Tigrfam when we copy folder
                output_hit_file_relative = os.path.join('.', tigrfam_version, filename.replace(
                    self.protein_file_ext, tigrfam_extension_gz))
                tigrfam_tophit_file_relative = os.path.join('.', tigrfam_version, filename.replace(
                    self.protein_file_ext, tigrfam_tophit_extension_gz))
                out_file_relative = os.path.join('.', tigrfam_version, filename.replace(
                    self.protein_file_ext, tigrfam_out_gz))

                symlink(output_hit_file_relative, new_hit_link,True)
                symlink(tigrfam_tophit_file_relative, new_tophit_link,True)
                symlink(out_file_relative, new_out_link,True)

            #we can now delete the temporary directory
            finally:
                shutil.rmtree(temp_dir)

            # allow results to be processed or written to file
            queue_out.put(gene_file)
