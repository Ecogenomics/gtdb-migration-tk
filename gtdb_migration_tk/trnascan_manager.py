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

"""Identify the tRNAs of every genome of a release with tRNAscan-SE.

WHY THE RELEASE IS CUT INTO BATCHES

tRNAscan-SE is seconds per genome and a release is a million-odd genomes, so the
run is days and belongs on several machines. The batching is the machinery
trans_table, prodigal and hmmsearch share, in batching.py: the release is
partitioned once under --out_dir, and a machine claims a batch directory before
working on it, so several machines given the same --out_dir divide the release
between them without being told which genomes to take. What --out_dir holds is
the state of the run and nothing else -- the tRNAs go into each genome's own
trna/ directory, as they always have, which is why two machines on different
batches never write to the same place.

WHAT DECIDES THE WORK

The checksum beside a genome's own results, and not a report. Each scan writes
<gid>_trna.tsv and a <gid>_trna.tsv.sha256 next to it, and a genome is skipped
where both are there and agree. This is the same rule hmmsearch follows and it is
the right one here for a reason worth writing down: trna is in
config.GTDB_DERIVED_DIRS_TO_COPY, so a genome whose sequences did not change
carries its tRNAs across from the previous release along with its checksum, and
is skipped without anything having to look up what became of it. The command
therefore takes no --report, and the commented-out report parse that sat in this
module until 0.1.29 was never needed.

WHY A DOMAIN IS LOOKED UP AT ALL, AND WHERE IT COMES FROM

tRNAscan-SE searches with a bacterial or an archaeal model and the two give
different answers, so every genome must be told which it is, and a genome told
wrong does not fail -- it gets a worse answer, which goes on into the tRNA counts
of the metadata tables. The domain therefore comes from the same two files
rna_silva reads it from, by utils.common.read_domains():

  --gtdb_domain_file    GTDB's own Predicted domain, from its marker genes
  --taxonomy_file       the standardised NCBI taxonomy, as a fallback

GTDB's call is preferred because it is made from the genome rather than from
where NCBI filed it, and it is the one that catches a genome under the wrong
domain at NCBI; the NCBI lineage answers for the genomes GTDB has no prediction
for, which is what 'None' in that column means. The taxonomy is matched on the
accession and then on its canonical form, so a GenBank genome finds the lineage
recorded against its RefSeq counterpart.

Until 0.1.29 the domain came from which of four NCBI assembly summary files an
accession appeared in -- nothing in those files says 'bacteria', so the answer
was asserted by which argument each file was passed as, and swapping two of them
on the command line would have scanned every archaeon as a bacterium in silence.
A genome neither file answers for is still scanned as a bacterium, which is what
this command has always done, but the run says how many of those there were.

ONE BAD GENOME DOES NOT COST A BATCH

A genome whose genomic FASTA is missing or empty is named in the batch's
not_scanned.tsv and left, and so is one tRNAscan-SE itself fails on; neither
takes the other ten thousand genomes of the batch down with it, and neither is
retried for ever by a batch that fails identically every time. A genome larger
than --max_genome_size is named and left too: it is a metagenome deposited as
one genome, and would hold its batch for as long as it took to scan. Its size is
asked only of the genomes that are to be scanned, so a rerun over a finished
release reads no statistics, and a genome whose tRNAs are already there is not
named for a scan it no longer needs. A batch in which
every genome was to be scanned and none could be is failed rather than recorded
as a success, because that is not a release of difficult genomes -- it is
tRNAscan-SE not working on this machine.
"""

import functools
import gzip
import logging
import multiprocessing as mp
import os
import shutil
import subprocess
import tempfile
from typing import Dict, List, NamedTuple, Optional, Sequence, Tuple

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
                                        split_by_genome_size, tally_reasons,
                                        write_table)
from gtdb_migration_tk.biolib_lite.checksum import sha256
from gtdb_migration_tk.biolib_lite.common import make_sure_path_exists
from gtdb_migration_tk.biolib_lite.external.execute import check_dependencies
# DOMAIN_BACTERIA and the domain file's names are imported so that they are still
# found here, where they were declared until rna_silva came to read them too.
from gtdb_migration_tk.utils.common import (DOMAIN_ARCHAEA, DOMAIN_BACTERIA,  # noqa: F401
                                            DOMAIN_FILE_DOMAIN,
                                            DEFAULT_MAX_GENOME_SIZE,
                                            DOMAIN_FILE_GENOME, MBP, NO_PREDICTION,
                                            domain_of, read_domains,
                                            record_program_version,
                                            write_version_file)


# What this command calls the files of its own batches. The batchfile's first
# column is each genome's genomic FASTA, which is what tRNAscan-SE reads. There
# is no older name to look for: this command has never had batches before.
BATCHFILE_NAME = 'trnascan_batchfile.tsv.gz'
BATCH_LOG_NAME = 'trnascan.log'
LAYOUT = BatchLayout(batchfiles=(BATCHFILE_NAME,), log=BATCH_LOG_NAME)

# Genomes of a batch that came out of it with no tRNAs, and the same gathered for
# the release. The step after this one should not have to look in 135 batch
# directories to find out which genomes those were.
NOT_SCANNED_NAME = 'not_scanned.tsv'
NOT_SCANNED_RELEASE_NAME = 'trnascan_not_scanned.tsv'
NOT_SCANNED_HEADER = ('genome_id', 'reason')
REASON_NO_GENOMIC_FASTA = 'no_genomic_fasta'
REASON_TRNASCAN_FAILED = 'trnascan_failed'
REASON_GENOME_TOO_LARGE = 'genome_too_large'

# Where a genome's tRNAs are written, within its own directory, and what they are
# called. metadata_manager reads <gid>_trna_stats.tsv from here by the accession
# the genome_dirs file gives, so these names are the agreement between the two
# commands.
TRNA_DIR = 'trna'
TRNA_EXT = '_trna.tsv'
TRNA_LOG_EXT = '_trna.log'
TRNA_STATS_EXT = '_trna_stats.tsv'
CHECKSUM_EXT = '.sha256'

# The program, as it is called. Its version goes into the log and, as
# trnascan-se.version, into the trna/ directory of each genome it scanned.
TRNASCAN = 'tRNAscan-SE'

# The two models tRNAscan-SE searches with. A genome of neither domain is scanned
# with the bacterial one, which is what this command has always done. How the
# domain is read, and the names of the domains and of the domain file's columns,
# are in utils/common.py, which rna_silva reads them from too.
ARCHAEAL_FLAG = '-A'
BACTERIAL_FLAG = '-B'


class TrnaJob(NamedTuple):
    """One genome as the workers receive it.

    The domain flag is settled in the parent, where the table naming it is, so
    that a worker carries what it needs rather than a copy of the release's
    domains.
    """

    accession: str
    genome_file: str
    domain_flag: str


class BatchCounts(NamedTuple):
    """What scanning one batch came to.

    Recorded in the batch's SUCCESS canary, because the release totals are added
    up from the batches and a machine that ran the last batch has scanned none of
    the others.
    """

    scanned: int
    already_scanned: int
    not_scanned: int


# The two things a worker of the pool does, as functions of the module rather than
# methods of tRNAScan. The pool pickles the function it is handed with every
# genome, and a bound method pickles its instance -- the domain table of the whole
# release with it, 2.7M entries and 47 MB, for each of the ten thousand genomes of
# a batch. The parent pickles them one at a time, 0.6 s apiece, so a batch spent
# hours 'Checking tRNAs' that is seconds of work, and the scan after it could feed
# its workers fewer than two genomes a second however many -c it was given. The
# domain is settled in the parent and travels in the TrnaJob, so a worker needs
# nothing of the instance. They log where the instance would have: 'timestamp'.
LOGGER = logging.getLogger('timestamp')


def trnascan_parser(job: TrnaJob) -> Optional[TrnaJob]:
    """Decide whether a genome's tRNAs still need identifying.

    A genome is skipped where its tRNA table is there AND its checksum
    agrees. A table whose checksum does not agree was written by a run that
    was interrupted partway through it, and is scanned again.

    Parameters
    ----------
    job : TrnaJob
        The genome to consider.

    @return: the job where the genome is to be scanned, None where its tRNAs
             are already there and vouched for.
    """

    trna_file = os.path.join(os.path.dirname(job.genome_file), TRNA_DIR,
                             job.accession + TRNA_EXT)
    checksum_file = trna_file + CHECKSUM_EXT

    if not (os.path.exists(trna_file) and os.path.exists(checksum_file)):
        return job

    checksum = str(sha256(trna_file))
    with open(checksum_file) as handle:
        recorded = handle.readline().strip()

    if checksum == recorded:
        return None

    LOGGER.warning(
        'warning: {} has tRNAs called with an invalid checksum ({} against '
        '{} in {}); the genome is scanned again.'.format(
            job.accession, checksum, recorded, checksum_file))

    return job

def trnascan_worker(job: TrnaJob, tmp_dir: str, version: str) -> Optional[str]:
    """Identify the tRNAs of one genome.

    The genome is decompressed into --tmp_dir first, tRNAscan-SE reading a
    plain FASTA, and the copy is removed however the scan ends. The results
    go into the genome's own directory, which is why two machines on
    different batches never write to the same place.

    Parameters
    ----------
    job : TrnaJob
        The genome to scan, and the model to scan it with.
    tmp_dir : str
        Directory the genome is decompressed into, one directory per genome.
    version : str
        tRNAscan-SE's version, recorded beside the tRNAs.

    @return: None where the genome was scanned, its accession where
             tRNAscan-SE failed on it -- which is reported and left rather
             than taking the rest of the batch down.
    """

    trna_dir = os.path.join(os.path.dirname(job.genome_file), TRNA_DIR)
    make_sure_path_exists(trna_dir)

    output_file = os.path.join(trna_dir, job.accession + TRNA_EXT)
    log_file = os.path.join(trna_dir, job.accession + TRNA_LOG_EXT)
    stats_file = os.path.join(trna_dir, job.accession + TRNA_STATS_EXT)

    # tRNAscan-SE reads a plain FASTA, and the genomes are held gzipped
    temp_dir = tempfile.mkdtemp(dir=tmp_dir)
    try:
        genome_copy = os.path.join(
            temp_dir, os.path.basename(job.genome_file)[:-len('.gz')])
        with gzip.open(job.genome_file, 'rb') as f_in:
            with open(genome_copy, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        # -o the tRNAs, -m the statistics metadata reads, -l the log
        command = [TRNASCAN, job.domain_flag, '-q', '-Q',
                   '-o', output_file, '-m', stats_file, '-l', log_file,
                   genome_copy]
        proc = subprocess.Popen(command, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        if proc.returncode != 0:
            LOGGER.warning(
                'warning: tRNAscan-SE failed on {} with status {}: {}'.format(
                    job.accession, proc.returncode,
                    (stderr or stdout).decode('utf-8', 'replace').strip()[:200]))
            return job.accession

        write_version_file(trna_dir, TRNASCAN, version)

        # written last, so a table without one is a scan that was interrupted
        with open(output_file + CHECKSUM_EXT, 'w') as handle:
            handle.write('{}\n'.format(sha256(output_file)))
    except Exception as error:
        LOGGER.warning(
            'warning: {} could not be scanned: {}'.format(job.accession, error))
        return job.accession
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)

    return None


class tRNAScan(object):
    """Runs tRNAscan-SE over the genomes of a release."""

    def __init__(self,
                 gtdb_domain_file: str,
                 taxonomy_file: str,
                 cpus: int = 1,
                 tmp_dir: str = '/tmp/',
                 batch_size: int = DEFAULT_BATCH_SIZE,
                 reclaim: bool = False,
                 lease: float = CLAIM_LEASE_SECONDS,
                 heartbeat: float = HEARTBEAT_SECONDS,
                 max_genome_size: float = DEFAULT_MAX_GENOME_SIZE) -> None:
        """Initialization.

        Parameters
        ----------
        gtdb_domain_file : str
            GTDB domain report, read for the domain predicted from each genome's
            marker genes.
        taxonomy_file : str
            Standardised NCBI taxonomy, read for the domain of the genomes GTDB
            has no prediction for.
        cpus : int
            How many genomes are scanned at once.
        tmp_dir : str
            Directory each genome is decompressed into before it is scanned, one
            directory per genome in flight and removed after. No results are
            written here.
        batch_size : int
            Genomes per batch.
        reclaim : bool
            Take over a batch another machine holds before its claim has expired.
        lease : float
            Seconds a claim survives without the machine holding it saying so.
        heartbeat : float
            Seconds between this machine saying so about a batch of its own.
        max_genome_size : float
            Largest genome assembly scanned, in Mbp.

        @return: None
        """

        check_dependencies([TRNASCAN])

        self.logger: logging.Logger = logging.getLogger('timestamp')
        self.version: str = record_program_version(TRNASCAN)

        self.cpus: int = cpus
        self.tmp_dir: str = tmp_dir
        self.batch_size: int = batch_size
        self.reclaim: bool = reclaim
        self.lease: float = lease
        self.heartbeat: float = heartbeat
        self.max_genome_bases: int = int(max_genome_size * MBP)

        # made here rather than by the first worker that wants it: a --tmp_dir
        # that cannot be made would otherwise be met once per genome, inside a
        # batch already claimed, and would fail every batch this machine took
        make_sure_path_exists(self.tmp_dir)

        self.domains: Dict[str, str] = read_domains(gtdb_domain_file,
                                                    taxonomy_file)

    def domain_flag(self, accession: str) -> str:
        """Which model tRNAscan-SE searches this genome with.

        Parameters
        ----------
        accession : str
            Genome accession, as the genome_dirs file names it.

        @return: the tRNAscan-SE flag, bacterial for a genome of no known domain.
        """

        if domain_of(self.domains, accession) == DOMAIN_ARCHAEA:
            return ARCHAEAL_FLAG

        return BACTERIAL_FLAG

    def run(self,
            gtdb_genome_path_file: str,
            out_dir: str,
            all_genomes: bool = False) -> bool:
        """Identify the tRNAs of every genome of a release.

        The release is cut into batches under --out_dir and a batch is claimed
        before it is worked on, so several machines can be pointed at one
        --out_dir and will divide the release between them.

        Parameters
        ----------
        gtdb_genome_path_file : str
            genome_dirs file of the release, accession and genome directory per line.
        out_dir : str
            Directory the batches and the state of the run are written to.
        all_genomes : bool
            Scan every genome again, discarding the tRNAs that are there.

        @return: True where every batch this machine took finished, False where
                 one failed and is left to a later run.
        """

        batches = plan_batches(gtdb_genome_path_file, out_dir, self.batch_size,
                               LAYOUT, self.logger)

        if all_genomes:
            self.logger.warning(
                'warning: --all discards the tRNAs of every genome and scans it '
                'again, so batches already finished are done again too. Without '
                'it a finished batch is skipped.')

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
                        counts = self.scan_batch(batch_dir, all_genomes)
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
                             scanned=counts.scanned,
                             already_scanned=counts.already_scanned,
                             not_scanned=counts.not_scanned)
                done += 1
                self.logger.info('{}: done.'.format(label))

        self.logger.info(
            '{:,} batch(es) finished here, {:,} held by another machine, '
            '{:,} failed.'.format(done, held, failed))

        self.aggregate(batches, out_dir)

        # a batch that failed has already said why, in its own log and in its
        # FAILED file, and a run of several machines over days should not end by
        # printing a traceback out of the argparse frames
        if failed:
            self.logger.error(
                '{:,} batch(es) failed; they are the directories holding a FAILED '
                'file and are retried by running the command again.'.format(failed))
            return False

        return True

    def scan_batch(self, batch_dir: str, all_genomes: bool = False) -> BatchCounts:
        """Identify the tRNAs of one batch, and record the genomes that got none.

        The genomes are taken from the batch's own batchfile, so the work asks
        about the genomes the batch was cut from rather than about whatever a
        genome_dirs file says now.

        Two passes, because deciding is cheap and scanning is not: the first
        sorts the batch into genomes whose tRNAs are already there and vouched
        for and genomes that are not, the second scans what is left.

        Parameters
        ----------
        batch_dir : str
            Batch directory.
        all_genomes : bool
            Discard the tRNAs of every genome of the batch and scan again.

        @return: what the batch came to, which run() records in its canary.
        """

        rows = read_batchfile(batchfile_path(batch_dir, LAYOUT))

        # a genome whose sequences are not there has nothing to scan; it is named
        # and left rather than stopping the other ten thousand of the batch
        present, missing = split_by_fasta(rows, STAT_THREADS)
        not_scanned = [(accession, REASON_NO_GENOMIC_FASTA) for accession in missing]

        jobs = [TrnaJob(accession, genome_file, self.domain_flag(accession))
                for genome_file, accession in present]

        unknown = [job.accession for job in jobs
                   if domain_of(self.domains, job.accession) is None]
        if unknown:
            self.logger.warning(
                'warning: {:,} genome(s) of this batch have no domain in either '
                'the GTDB domain file or the NCBI taxonomy and are scanned as '
                'bacteria, e.g. {}.'.format(
                    len(unknown), ', '.join(sorted(unknown)[:3])))

        if all_genomes:
            to_scan, already_scanned = jobs, 0
        else:
            with mp.Pool(processes=self.cpus) as pool:
                decided = list(tqdm(pool.imap_unordered(trnascan_parser, jobs),
                                    total=len(jobs), unit='genome', ncols=100,
                                    leave=False, desc='Checking tRNAs'))
            to_scan = [job for job in decided if job is not None]
            already_scanned = len(jobs) - len(to_scan)

        to_scan, too_large = split_by_genome_size(
            to_scan, lambda job: os.path.dirname(job.genome_file),
            self.max_genome_bases, self.logger, STAT_THREADS)
        not_scanned.extend((job.accession, REASON_GENOME_TOO_LARGE)
                           for job, _ in too_large)

        self.logger.info(
            '{:,} genome(s) require tRNA identification; {:,} already have valid '
            'results.'.format(len(to_scan), already_scanned))

        not_scanned.extend(self.scan_genomes(to_scan))

        not_scanned.sort()
        write_table(not_scanned, os.path.join(batch_dir, NOT_SCANNED_NAME),
                    header=NOT_SCANNED_HEADER)

        if not_scanned:
            self.logger.warning(
                'warning: {:,} genome(s) of this batch have no tRNAs: {}.'.format(
                    len(not_scanned),
                    '; '.join('{:,} {}'.format(count, reason)
                              for reason, count
                              in sorted(tally_reasons(not_scanned).items()))))

        scanned = len(to_scan) - sum(1 for _, reason in not_scanned
                                     if reason == REASON_TRNASCAN_FAILED)

        return BatchCounts(scanned=scanned,
                           already_scanned=already_scanned,
                           not_scanned=len(not_scanned))

    def scan_genomes(self, jobs: Sequence[TrnaJob]) -> List[Tuple[str, str]]:
        """Run tRNAscan-SE over the genomes of a batch that need it.

        Parameters
        ----------
        jobs : sequence of TrnaJob
            The genomes to scan.

        @return: (accession, reason) for each genome tRNAscan-SE failed on.

        Raises
        ------
        RuntimeError
            Every one of several genomes failed, which is tRNAscan-SE not working
            on this machine rather than a batch of difficult genomes, and is a
            reason to fail the batch instead of recording it as a success with
            nothing in it. One genome failing on its own is not: a batch whose
            last unscanned genome is one tRNAscan-SE refuses would otherwise fail
            identically on every retry, which is what naming the genome and
            carrying on exists to avoid.
        """

        if not jobs:
            return []

        with mp.Pool(processes=self.cpus) as pool:
            worker = functools.partial(trnascan_worker, tmp_dir=self.tmp_dir,
                                       version=self.version)
            results = list(tqdm(pool.imap_unordered(worker, jobs),
                                total=len(jobs), unit='genome', ncols=100,
                                desc='Identifying tRNAs'))

        failures = [(accession, REASON_TRNASCAN_FAILED)
                    for accession in results if accession is not None]

        if len(jobs) > 1 and len(failures) == len(jobs):
            raise RuntimeError(
                'tRNAscan-SE failed on every one of the {:,} genome(s) it was '
                'given.'.format(len(jobs)))

        return failures

    def aggregate(self, batches: Sequence[str], out_dir: str) -> None:
        """Report the release, once every batch has succeeded.

        Written only when they all have, so that the file at the top of the
        directory is either the whole release or absent, and never a part of it
        that reads like the whole.

        Parameters
        ----------
        batches : sequence of str
            Every batch directory of the run, in batch order.
        out_dir : str
            Directory the batches are held under.

        @return: None
        """

        unfinished = [batch for batch in batches
                      if batch_state(batch) != STATE_SUCCESS]
        if unfinished:
            self.logger.info(
                '{:,} of {:,} batch(es) are done; the release files are written '
                'once they all are.'.format(
                    len(batches) - len(unfinished), len(batches)))
            return

        path = os.path.join(out_dir, NOT_SCANNED_RELEASE_NAME)
        written = concatenate(
            [os.path.join(batch, NOT_SCANNED_NAME) for batch in batches], path)

        totals = {'scanned': 0, 'already_scanned': 0}
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
            'Release: {:,} genome(s) had their tRNAs identified, {:,} already had '
            'valid results.'.format(totals['scanned'], totals['already_scanned']))

        if written:
            self.logger.warning(
                'warning: {:,} genome(s) of the release have no tRNAs and are '
                'named in {}.'.format(written, path))
        else:
            self.logger.info(
                'Every genome of the release has tRNAs; wrote {} with no '
                'rows.'.format(path))
