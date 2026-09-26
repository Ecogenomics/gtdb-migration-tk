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

"""Classify the 16S rRNA genes rna_silva extracted against the LTP.

WHY THE RELEASE IS CUT INTO BATCHES

blastn is seconds per genome and a release is a million-odd genomes, so the run
is days and belongs on several machines. The batching is the machinery
trans_table, prodigal, hmmsearch, trnascan and rna_silva share, in batching.py:
the release is partitioned once under --out_dir and a machine claims a batch
directory before working on it. What --out_dir holds is the state of the run and
nothing else -- the classifications go into each genome's own rna_ltp directory,
as they always have.

The batches live under <out_dir>/rna_ltp_<ltp version>-silva_<ssu version>/, as
both versions decide what a genome's results are: the LTP is what the genes are
classified against, and the SILVA version names the rna_silva directory the genes
are read from.

WHAT A BATCH IS PLANNED AROUND

rna_silva's ssu.fna, which is the file this command reads -- not the genomic
FASTA, which it never opens. A genome with no ssu.fna is one of two things, and
rna_silva's own canary tells them apart: where ssu.canary.txt is there, rna_silva
searched the genome and found no 16S rRNA gene, so there is nothing to classify
and nothing wrong, and it is counted but not named; where it is not, rna_silva has
not searched the genome yet, and it is named in not_classified.tsv as
ssu_not_identified, so that the order the two commands ran in is not mistaken for
a release of genomes without a 16S gene.

A GENOME TOO LARGE TO CLASSIFY

A genome larger than --max_genome_size is named as genome_too_large and left, as
rna_silva names and leaves one larger than its own: it is a metagenome deposited
as one genome. rna_silva writes no canary for a genome it left out, so a genome
with no ssu.fna and no canary is sized too, and one rna_silva left for its size
is named for that rather than as one rna_silva has yet to search. The size is
asked of no other genome, a genome already classified or found to have no 16S
gene having nothing left to be spared.

WHAT DECIDES THE WORK

ltp.canary.txt beside a genome's results, as it always has been: the rna_ltp
directory is in config.GTDB_DERIVED_DIRS_TO_COPY, so a genome whose sequences did
not change carries its classification across from the previous release with the
canary among them. The classification is made in --tmp_dir and copied into place
with the canary last, so that a run stopped partway through the copy leaves a
genome classified again rather than one that looks finished with half its files.

WHY THERE IS NO DOMAIN

Until 0.1.36 the command read the domain of each genome from --gtdb_domain_file,
and handed it to genometk_lite's RNA, whose classify() never reads it: the domain
chooses the HMM and the minimum gene length when a gene is searched for, which is
rna_silva's work, and the LTP classification of a gene already extracted is the
same whatever domain it is said to be of. The file and the lookup are gone.

ONE BAD GENOME DOES NOT COST A BATCH

A genome blastn fails on is named as blastn_failed and left. BLAST's wrapper
raised on nothing until 0.1.36: blastn was run through os.system(), a failed
search left an empty table, and the genome was written a canary saying it had
been classified with no hits. A batch in which several genomes were to be
classified and every one failed is failed, that being blastn or the LTP database
not working on this machine rather than a batch of difficult genomes.
"""

import functools
import logging
import multiprocessing as mp
import os
import shutil
import sys
import tempfile
from typing import List, NamedTuple, Optional, Sequence, Tuple

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
from gtdb_migration_tk.biolib_lite.common import (make_sure_path_exists,
                                                  remove_files_in_directory)
from gtdb_migration_tk.biolib_lite.external.execute import check_dependencies
from gtdb_migration_tk.genometk_lite.rna import RNA
from gtdb_migration_tk.utils.common import (DEFAULT_MAX_GENOME_SIZE, MBP,
                                            record_program_version,
                                            write_version_file)

# The program, as it is called. Its version goes into the log and, as
# blastn.version, into the rna_ltp directory of each genome it classified.
BLASTN = 'blastn'

# What this command calls the files of its own batches. The batchfile's first
# column is each genome's ssu.fna, which is what blastn reads.
BATCHFILE_NAME = 'rna_ltp_batchfile.tsv.gz'
BATCH_LOG_NAME = 'rna_ltp.log'
LAYOUT = BatchLayout(batchfiles=(BATCHFILE_NAME,), log=BATCH_LOG_NAME)

# Genomes of a batch that came out of it unclassified, and the same gathered for
# the run once every batch has SUCCESS. A genome rna_silva found no 16S gene in
# is not among them: it has nothing to classify.
NOT_CLASSIFIED_NAME = 'not_classified.tsv'
NOT_CLASSIFIED_RELEASE_NAME = 'rna_ltp_not_classified.tsv'
NOT_CLASSIFIED_HEADER = ('genome_id', 'reason')
REASON_SSU_NOT_IDENTIFIED = 'ssu_not_identified'
REASON_BLASTN_FAILED = 'blastn_failed'
REASON_GENOME_TOO_LARGE = 'genome_too_large'

# Where rna_silva left the 16S rRNA genes, and what says it searched for them.
# The name is rna_silva's RESULTS_DIR_FORMAT, repeated here because the command
# modules do not import one another; tests/test_rna_manager_ltp.py holds the two
# to agreeing, and to the directory config.py carries across between releases.
SSU_RESULTS_DIR_FORMAT = 'rna_silva_{}'
SSU_FASTA = 'ssu.fna'
SSU_CANARY = 'ssu.canary.txt'

# Where the classification goes inside each genome directory, what it is called,
# and what says it is done.
LTP_RESULTS_DIR_FORMAT = 'rna_ltp_{}'
LTP_RESULTS_PREFIX = 'ssu.'
LTP_CANARY = 'ltp.canary.txt'

# RNA takes a domain it does not read when it classifies, so it is handed this.
UNUSED_DOMAIN = 'bac'


class LtpJob(NamedTuple):
    """One genome as the workers receive it."""

    accession: str
    ssu_file: str


class LtpSettings(NamedTuple):
    """What every genome of the run is classified with, handed to each worker."""

    ltp_ssu_file: str
    ltp_taxonomy_file: str
    ltp_version: str
    ssu_version: str
    results_dir: str
    tmp_dir: str
    blastn_version: str
    remove_prior: bool


class BatchCounts(NamedTuple):
    """What classifying one batch came to, recorded in its SUCCESS canary."""

    classified: int
    already_classified: int
    no_ssu_gene: int
    not_classified: int


LOGGER = logging.getLogger('timestamp')


def ssu_fasta(ssu_version: str, accession: str, genome_dir: str) -> str:
    """The file a batch is planned around: the 16S genes rna_silva extracted.

    Parameters
    ----------
    ssu_version : str
        SILVA version naming rna_silva's results directory.
    accession : str
        Accession of the genome, unused here and taken for the signature.
    genome_dir : str
        Genome directory of the release.

    @return: path of ssu.fna in the genome's rna_silva directory.
    """

    return os.path.join(genome_dir, SSU_RESULTS_DIR_FORMAT.format(ssu_version),
                        SSU_FASTA)


def genome_dir_of(ssu_file: str) -> str:
    """The genome directory holding the rna_silva directory an ssu.fna is in.

    Parameters
    ----------
    ssu_file : str
        The genome's ssu.fna, inside its rna_silva directory.

    @return: path of the genome directory.
    """

    return os.path.dirname(os.path.dirname(ssu_file))


def output_dir_of(ssu_file: str, results_dir: str) -> str:
    """The rna_ltp directory of the genome whose ssu.fna this is.

    Parameters
    ----------
    ssu_file : str
        The genome's ssu.fna, inside its rna_silva directory.
    results_dir : str
        Name of the rna_ltp directory, e.g. rna_ltp_10_2024.

    @return: path of the rna_ltp directory inside the genome directory.
    """

    return os.path.join(genome_dir_of(ssu_file), results_dir)


def ltp_parser(job: LtpJob, results_dir: str) -> Optional[LtpJob]:
    """Decide whether a genome's 16S genes still need classifying.

    Parameters
    ----------
    job : LtpJob
        The genome to consider.
    results_dir : str
        Name of the rna_ltp directory inside the genome directory.

    @return: the job where the genome is to be classified, None where its
             canary says it already was.
    """

    canary = os.path.join(output_dir_of(job.ssu_file, results_dir), LTP_CANARY)

    return None if os.path.exists(canary) else job


def ltp_worker(job: LtpJob, settings: LtpSettings) -> Optional[str]:
    """Classify the 16S genes of one genome against the LTP.

    Parameters
    ----------
    job : LtpJob
        The genome to classify.
    settings : LtpSettings
        What every genome of the run is classified with.

    @return: None where the genome was classified, its accession where blastn
             failed on it -- which is reported and left rather than taking the
             rest of the batch down.
    """

    output_dir = output_dir_of(job.ssu_file, settings.results_dir)
    try:
        if settings.remove_prior and os.path.isdir(output_dir):
            remove_files_in_directory(output_dir)

        temp_dir = tempfile.mkdtemp(dir=settings.tmp_dir, prefix=job.accession)
        try:
            RNA('ssu', UNUSED_DOMAIN, 1).classify(job.ssu_file,
                                                  settings.ltp_ssu_file,
                                                  settings.ltp_taxonomy_file,
                                                  temp_dir)
            with open(os.path.join(temp_dir, LTP_CANARY), 'w') as handle:
                handle.write(f'Silva version:{settings.ssu_version}.\n')
                handle.write(f'LTP version:{settings.ltp_version}.\n')
                handle.write('done.\n')

            install_results(temp_dir, output_dir, settings)
        finally:
            shutil.rmtree(temp_dir, ignore_errors=True)
    except Exception as error:
        LOGGER.warning('warning: the 16S rRNA genes of {} could not be classified '
                       'against the LTP: {}'.format(job.accession, error))
        return job.accession

    return None


def install_results(temp_dir: str, output_dir: str, settings: LtpSettings) -> None:
    """Copy one genome's classification into place, the canary last.

    What an earlier classification left is removed first, so that nothing of it
    outlives the one replacing it.

    Parameters
    ----------
    temp_dir : str
        Where the classification was made.
    output_dir : str
        The genome's rna_ltp directory.
    settings : LtpSettings
        What the genome was classified with.

    @return: None
    """

    make_sure_path_exists(output_dir)

    for name in os.listdir(output_dir):
        if name.startswith(LTP_RESULTS_PREFIX) or name == LTP_CANARY:
            os.remove(os.path.join(output_dir, name))

    for name in sorted(os.listdir(temp_dir)):
        if name != LTP_CANARY:
            shutil.copy(os.path.join(temp_dir, name), os.path.join(output_dir, name))

    write_version_file(output_dir, BLASTN, settings.blastn_version)

    shutil.copy(os.path.join(temp_dir, LTP_CANARY),
                os.path.join(output_dir, LTP_CANARY))


def ssu_searched(ssu_file: str) -> bool:
    """Whether rna_silva searched a genome for its 16S gene.

    Parameters
    ----------
    ssu_file : str
        Where the genome's ssu.fna would be.

    @return: True where rna_silva's canary is beside where ssu.fna would be.
    """

    return os.path.exists(os.path.join(os.path.dirname(ssu_file), SSU_CANARY))


class RnaManagerLTP(object):
    """Classify the 16S rRNA genes of a release against the LTP."""

    def __init__(self,
                 ltp_version: str,
                 silva_version: str,
                 rna_path: str,
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
        ltp_version : str
            LTP release classified against, e.g. 10_2024.
        silva_version : str
            SILVA release whose rna_silva directory the 16S genes are read from.
        rna_path : str
            Directory holding a directory per LTP release.
        cpus : int
            How many genomes are classified at once.
        tmp_dir : str
            Directory each genome is classified in before its results are copied
            into place.
        batch_size : int
            Genomes per batch.
        reclaim : bool
            Take over a batch another machine holds before its claim has expired.
        lease : float
            Seconds a claim survives without the machine holding it saying so.
        heartbeat : float
            Seconds between this machine saying so about a batch of its own.
        max_genome_size : float
            Largest genome assembly classified, in Mbp.

        @return: None
        """

        self.logger: logging.Logger = logging.getLogger('timestamp')

        check_dependencies([BLASTN])
        self.blastn_version: str = record_program_version(BLASTN)

        self.ltp_version: str = ltp_version
        self.ssu_version: str = silva_version
        self.cpus: int = cpus
        self.tmp_dir: str = tmp_dir
        self.batch_size: int = batch_size
        self.reclaim: bool = reclaim
        self.lease: float = lease
        self.heartbeat: float = heartbeat
        self.max_genome_bases: int = int(max_genome_size * MBP)

        self.results_dir: str = LTP_RESULTS_DIR_FORMAT.format(ltp_version)

        ltp_root_path = os.path.join(rna_path, str(ltp_version))
        self.ltp_ssu_file: str = os.path.join(
            ltp_root_path, 'ltp_{}.fna'.format(ltp_version))
        self.ltp_taxonomy_file: str = os.path.join(
            ltp_root_path, 'ltp_{}_taxonomy.tsv'.format(ltp_version))

        for item in [self.ltp_ssu_file, self.ltp_taxonomy_file]:
            if not os.path.exists(item):
                self.logger.error('{} does not exist'.format(item))
                sys.exit(-1)

        # made here rather than by the first worker that wants it, so that a
        # --tmp_dir that cannot be made is met before a batch is claimed
        make_sure_path_exists(self.tmp_dir)

    def run_dir(self, out_dir: str) -> str:
        """Where the batches of this LTP and SILVA version live under --out_dir.

        Parameters
        ----------
        out_dir : str
            Output directory of the run.

        @return: <out_dir>/rna_ltp_<ltp version>-silva_<ssu version>.
        """

        return os.path.join(out_dir, '{}-silva_{}'.format(self.results_dir,
                                                         self.ssu_version))

    def run(self,
            gtdb_genome_path_file: str,
            out_dir: str,
            all_genomes: bool = False,
            remove_prior: bool = False) -> bool:
        """Classify the 16S rRNA genes of every genome of a release.

        Parameters
        ----------
        gtdb_genome_path_file : str
            genome_dirs file of the release, accession and genome directory per line.
        out_dir : str
            Directory the batches and the state of the run are written to.
        all_genomes : bool
            Classify every genome again, whatever its canary says.
        remove_prior : bool
            Empty the rna_ltp directory of each genome before it is classified.

        @return: True where every batch this machine took finished, False where
                 one failed and is left to a later run.
        """

        run_dir = self.run_dir(out_dir)
        batches = plan_batches(gtdb_genome_path_file, run_dir, self.batch_size,
                               LAYOUT, self.logger,
                               genome_file=functools.partial(ssu_fasta,
                                                             self.ssu_version))

        if all_genomes:
            self.logger.warning(
                'warning: --all classifies every genome again, so batches already '
                'finished are done again too. Without it a finished batch is '
                'skipped.')

        settings = LtpSettings(ltp_ssu_file=self.ltp_ssu_file,
                               ltp_taxonomy_file=self.ltp_taxonomy_file,
                               ltp_version=self.ltp_version,
                               ssu_version=self.ssu_version,
                               results_dir=self.results_dir,
                               tmp_dir=self.tmp_dir,
                               blastn_version=self.blastn_version,
                               remove_prior=remove_prior)

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

            with batch_log(batch_dir, self.logger, LAYOUT):
                self.logger.info('{}: starting.'.format(label))
                try:
                    with Heartbeat(os.path.join(batch_dir, RUNNING_CANARY),
                                   self.heartbeat):
                        counts = self.classify_batch(batch_dir, settings, all_genomes)
                except KeyboardInterrupt:
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
                             classified=counts.classified,
                             already_classified=counts.already_classified,
                             no_ssu_gene=counts.no_ssu_gene,
                             not_classified=counts.not_classified)
                done += 1
                self.logger.info('{}: done.'.format(label))

        self.logger.info(
            '{:,} batch(es) finished here, {:,} held by another machine, '
            '{:,} failed.'.format(done, held, failed))

        self.aggregate(batches, run_dir)

        if failed:
            self.logger.error(
                '{:,} batch(es) failed; they are the directories holding a FAILED '
                'file and are retried by running the command again.'.format(failed))
            return False

        return True

    def classify_batch(self, batch_dir: str, settings: LtpSettings,
                       all_genomes: bool = False) -> BatchCounts:
        """Classify one batch, and record the genomes that could not be.

        Parameters
        ----------
        batch_dir : str
            Batch directory.
        settings : LtpSettings
            What every genome of the run is classified with.
        all_genomes : bool
            Classify every genome of the batch again, whatever its canary says.

        @return: what the batch came to, which run() records in its canary.
        """

        rows = read_batchfile(batchfile_path(batch_dir, LAYOUT))
        ssu_of = {accession: ssu_file for ssu_file, accession in rows}

        # a genome with no ssu.fna has nothing to classify; rna_silva's canary
        # says whether that is because it has no 16S gene or because rna_silva
        # has not looked yet, and only the second is worth naming
        present, missing = split_by_fasta(rows, STAT_THREADS)
        unsearched = [accession for accession in missing
                      if not ssu_searched(ssu_of[accession])]
        no_ssu_gene = len(missing) - len(unsearched)

        # rna_silva writes no canary for a genome it left out for its size, and
        # that genome is named for its size rather than as one not yet searched
        unsearched, too_large = split_by_genome_size(
            unsearched, lambda accession: genome_dir_of(ssu_of[accession]),
            self.max_genome_bases, self.logger, STAT_THREADS)
        not_classified = [(accession, REASON_SSU_NOT_IDENTIFIED)
                          for accession in unsearched]
        not_classified.extend((accession, REASON_GENOME_TOO_LARGE)
                              for accession, _ in too_large)

        jobs = [LtpJob(accession, ssu_file) for ssu_file, accession in present]

        if all_genomes:
            to_classify, already_classified = jobs, 0
        else:
            parser = functools.partial(ltp_parser, results_dir=self.results_dir)
            with mp.Pool(processes=self.cpus) as pool:
                decided = list(tqdm(pool.imap_unordered(parser, jobs),
                                    total=len(jobs), unit='genome', ncols=100,
                                    leave=False, desc='Checking LTP results'))
            to_classify = [job for job in decided if job is not None]
            already_classified = len(jobs) - len(to_classify)

        to_classify, too_large = split_by_genome_size(
            to_classify, lambda job: genome_dir_of(job.ssu_file),
            self.max_genome_bases, self.logger, STAT_THREADS)
        not_classified.extend((job.accession, REASON_GENOME_TOO_LARGE)
                              for job, _ in too_large)

        self.logger.info(
            '{:,} genome(s) require classification against the LTP; {:,} already '
            'have results and {:,} have no 16S rRNA gene.'.format(
                len(to_classify), already_classified, no_ssu_gene))

        not_classified.extend(self.classify_genomes(to_classify, settings))

        not_classified.sort()
        write_table(not_classified, os.path.join(batch_dir, NOT_CLASSIFIED_NAME),
                    header=NOT_CLASSIFIED_HEADER)

        if not_classified:
            self.logger.warning(
                'warning: {:,} genome(s) of this batch were not classified: {}.'.format(
                    len(not_classified),
                    '; '.join('{:,} {}'.format(count, reason)
                              for reason, count
                              in sorted(tally_reasons(not_classified).items()))))

        classified = len(to_classify) - sum(1 for _, reason in not_classified
                                            if reason == REASON_BLASTN_FAILED)

        return BatchCounts(classified=classified,
                           already_classified=already_classified,
                           no_ssu_gene=no_ssu_gene,
                           not_classified=len(not_classified))

    def classify_genomes(self, jobs: Sequence[LtpJob],
                         settings: LtpSettings) -> List[Tuple[str, str]]:
        """Classify the genomes of a batch that need it.

        Parameters
        ----------
        jobs : sequence of LtpJob
            The genomes to classify.
        settings : LtpSettings
            What every genome of the run is classified with.

        @return: (accession, reason) for each genome blastn failed on.

        Raises
        ------
        RuntimeError
            Every one of several genomes failed, which is blastn or the LTP
            database not working on this machine rather than a batch of
            difficult genomes.
        """

        if not jobs:
            return []

        worker = functools.partial(ltp_worker, settings=settings)
        with mp.Pool(processes=self.cpus) as pool:
            results = list(tqdm(pool.imap_unordered(worker, jobs),
                                total=len(jobs), unit='genome', ncols=100,
                                desc='Classifying 16S rRNA genes'))

        failures = [(accession, REASON_BLASTN_FAILED)
                    for accession in results if accession is not None]

        if len(jobs) > 1 and len(failures) == len(jobs):
            raise RuntimeError(
                'blastn failed on every one of the {:,} genome(s) it was '
                'given.'.format(len(jobs)))

        return failures

    def aggregate(self, batches: Sequence[str], run_dir: str) -> None:
        """Report the release, once every batch has succeeded.

        Parameters
        ----------
        batches : sequence of str
            Every batch directory of the run, in batch order.
        run_dir : str
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

        path = os.path.join(run_dir, NOT_CLASSIFIED_RELEASE_NAME)
        written = concatenate(
            [os.path.join(batch, NOT_CLASSIFIED_NAME) for batch in batches], path)

        totals = {'classified': 0, 'already_classified': 0, 'no_ssu_gene': 0}
        for batch_dir in batches:
            canary = read_canary(os.path.join(batch_dir, SUCCESS_CANARY))
            for field in totals:
                try:
                    totals[field] += int(canary[field])
                except (KeyError, ValueError):
                    pass

        self.logger.info(
            'Release: {:,} genome(s) were classified against the LTP, {:,} already '
            'had results and {:,} have no 16S rRNA gene.'.format(
                totals['classified'], totals['already_classified'],
                totals['no_ssu_gene']))

        if written:
            self.logger.warning(
                'warning: {:,} genome(s) of the release were not classified and '
                'are named in {}.'.format(written, path))
        else:
            self.logger.info(
                'Every genome of the release with a 16S rRNA gene was classified; '
                'wrote {} with no rows.'.format(path))
