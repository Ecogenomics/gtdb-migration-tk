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

"""Identify, extract and classify the rRNA genes of every genome of a release.

WHY THE RELEASE IS CUT INTO BATCHES

nhmmer and blastn are seconds per genome and a release is a million-odd genomes,
three times over -- once each for --rna_gene ssu, lsu_23S and lsu_5S -- so the
run is days and belongs on several machines. The batching is the machinery
trans_table, prodigal, hmmsearch and trnascan share, in batching.py: the release
is partitioned once under --out_dir, and a machine claims a batch directory
before working on it, so several machines given the same --out_dir divide the
release between them without being told which genomes to take. What --out_dir
holds is the state of the run and nothing else -- the rRNA genes go into each
genome's own rna_silva directory, as they always have, which is why two machines
on different batches never write to the same place.

The batches of a run live under <out_dir>/<rna_silva directory>/<rna_gene>/, for
the reason hmmsearch keeps a directory per database: an ssu batch and an lsu_23S
batch cover the same genomes and are different work, and sharing a directory
would have one's SUCCESS tell the other it had nothing to do. The SILVA version
is in the path for the same reason.

WHAT DECIDES THE WORK

The canary beside a genome's own results, <rna_gene>.canary.txt, which
genometk_lite's RNA writes once the gene has been searched for, extracted and
classified -- whether or not the genome turned out to have one. This is the rule
this command has always followed and it is the right one for the reason trnascan
gives for its checksum: the rna_silva directory is in
config.GTDB_DERIVED_DIRS_TO_COPY, so a genome whose sequences did not change
carries its results across from the previous release with the canary among them,
and is skipped without anything having to look up what became of it.

The results are made in --tmp_dir and copied into the genome directory with the
canary copied LAST, so that a run stopped partway through the copy leaves a
genome that is searched again rather than one that looks finished with half its
files. The files an earlier search of the same gene left are removed first, so a
genome searched again that no longer has a gene does not keep the taxonomy of
the one it had.

WHERE THE DOMAIN COMES FROM

The rRNA HMMs are per domain, bac_16S against ar_16S and so on, so every genome
must be told which it is, and a genome told wrong does not fail -- it gets a
worse answer. The domain is read as trnascan reads it, by
utils.common.read_domains():

  --gtdb_domain_file    GTDB's own Predicted domain, from its marker genes
  --taxonomy_file       the standardised NCBI taxonomy, as a fallback

Until 0.1.35 the fallback was the NCBI taxonomy column of the domain file itself,
and a genome whose domain was neither Bacteria nor Archaea stopped the run. The
standardised NCBI taxonomy is the one trnascan reads and is matched on the
accession and then on its canonical form, so a GenBank genome finds the lineage
recorded against its RefSeq counterpart. A genome neither file answers for is
still searched as a bacterium, which is what this command has always done, and
each batch says how many of those it had.

ONE BAD GENOME DOES NOT COST A BATCH

A genome whose genomic FASTA is missing or empty is named in the batch's
not_searched.tsv and left, and so is one nhmmer or blastn fails on; neither takes
the other ten thousand genomes of the batch down with it. A genome larger than
--max_genome_size is named and left too, before nhmmer is given it: it is a
metagenome deposited as one genome, and GCA_964261755.1, at 9,529 Mbp, held an
r237 ssu batch for a day on one single-threaded blastn classifying the 16S genes
nhmmer found in it. Only the genomes still to be searched are sized, so one
already searched is not named for work it no longer needs. A batch in which every
genome was to be searched and none could be is failed rather than recorded as a
success, because that is the tools not working on this machine rather than a
batch of difficult genomes.

WHY THE WORKERS ARE FUNCTIONS OF THE MODULE

The pool pickles the function it is handed with every genome, and a bound method
pickles its instance -- the domain table of the whole release with it, which is
what cost trnascan hours a batch until 0.1.34. The domain is settled in the
parent and travels in the RnaJob, and what every genome shares travels in one
small RnaSettings, so a task is a few hundred bytes.
"""

import functools
import gzip
import logging
import multiprocessing as mp
import os
import shutil
import sys
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
from gtdb_migration_tk.biolib_lite.common import (make_sure_path_exists,
                                                  remove_files_in_directory)
from gtdb_migration_tk.biolib_lite.external.execute import check_dependencies
from gtdb_migration_tk.biolib_lite.seq_io import read_seq
from gtdb_migration_tk.genometk_lite.rna import RNA
from gtdb_migration_tk.utils.common import (DEFAULT_MAX_GENOME_SIZE,
                                            DOMAIN_ARCHAEA, MBP, domain_of,
                                            read_domains,
                                            record_program_version,
                                            write_version_file)

# The programs, as they are called. Their versions go into the log and, as
# nhmmer.version and blastn.version, into the rna_silva directory of each genome
# -- blastn's only where it ran, which is where nhmmer found a gene to classify
# and there is a database to classify it against.
NHMMER = 'nhmmer'
BLASTN = 'blastn'

# What this command calls the files of its own batches. The batchfile's first
# column is each genome's genomic FASTA, which is what nhmmer reads. There is no
# older name to look for: this command has never had batches before.
BATCHFILE_NAME = 'rna_silva_batchfile.tsv.gz'
BATCH_LOG_NAME = 'rna_silva.log'
LAYOUT = BatchLayout(batchfiles=(BATCHFILE_NAME,), log=BATCH_LOG_NAME)

# Genomes of a batch that came out of it unsearched, and the same gathered for
# the run once every batch has SUCCESS.
NOT_SEARCHED_NAME = 'not_searched.tsv'
NOT_SEARCHED_RELEASE_NAME = 'rna_silva_not_searched.tsv'
NOT_SEARCHED_HEADER = ('genome_id', 'reason')
REASON_NO_GENOMIC_FASTA = 'no_genomic_fasta'
REASON_SEARCH_FAILED = 'rna_search_failed'
REASON_GENOME_TOO_LARGE = 'genome_too_large'

# The directory inside each genome directory the results go into, for a SILVA
# version. It is the agreement with three readers: metadata_manager reads
# rna_silva_<version>/, config.py carries that directory across between releases,
# and rna_ltp reads its ssu.fna from there. Until 0.1.36 this was
# rna_silva_<version>-testing/, left from when the rRNA workflow was rewritten to
# follow Prokka's so that its results sat beside the old workflow's, and none of
# the three read what this command wrote.
RESULTS_DIR_FORMAT = 'rna_silva_{}'

# What RNA writes once a gene has been searched for, and what decides the work.
CANARY_EXT = '.canary.txt'

# The HMM each domain's genome is searched with, for each rRNA gene.
HMM_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                       'data_files', 'barrnap')
RNA_MODELS = {'ssu': {'ar': 'ar_16S', 'bac': 'bac_16S', 'euk': 'euk_18S'},
              'lsu_23S': {'ar': 'ar_23S', 'bac': 'bac_23S', 'euk': 'euk_28S'},
              'lsu_5S': {'ar': 'ar_5S', 'bac': 'bac_5S', 'euk': 'euk_5S'}}

# How RNA spells the two domains read_domains() returns. A genome of neither is
# searched as a bacterium, which is what this command has always done.
RNA_ARCHAEA = 'ar'
RNA_BACTERIA = 'bac'


class RnaJob(NamedTuple):
    """One genome as the workers receive it.

    The domain is settled in the parent, where the table naming it is, so that a
    worker carries what it needs rather than a copy of the release's domains.
    """

    accession: str
    genome_file: str
    domain: str


class RnaSettings(NamedTuple):
    """What every genome of the run is searched with, handed to each worker."""

    rna_gene: str
    db: Optional[str]
    taxonomy: Optional[str]
    results_dir: str
    tmp_dir: str
    nhmmer_version: str
    blastn_version: str
    remove_prior: bool


class BatchCounts(NamedTuple):
    """What searching one batch came to.

    Recorded in the batch's SUCCESS canary, because the release totals are added
    up from the batches and a machine that ran the last batch has searched none
    of the others.
    """

    searched: int
    already_searched: int
    not_searched: int


LOGGER = logging.getLogger('timestamp')


def output_dir_of(genome_file: str, results_dir: str) -> str:
    """The directory a genome's rRNA results go into.

    Parameters
    ----------
    genome_file : str
        The genome's genomic FASTA, in its genome directory.
    results_dir : str
        Name of the results directory, e.g. rna_silva_138.2.

    @return: path of the results directory inside the genome directory.
    """

    return os.path.join(os.path.dirname(genome_file), results_dir)


def rna_parser(job: RnaJob, rna_gene: str, results_dir: str) -> Optional[RnaJob]:
    """Decide whether a genome's rRNA gene still needs searching for.

    Parameters
    ----------
    job : RnaJob
        The genome to consider.
    rna_gene : str
        The gene of this run.
    results_dir : str
        Name of the results directory inside the genome directory.

    @return: the job where the genome is to be searched, None where its canary
             says the gene was already searched for.
    """

    canary = os.path.join(output_dir_of(job.genome_file, results_dir),
                          rna_gene + CANARY_EXT)

    return None if os.path.exists(canary) else job


def rna_worker(job: RnaJob, settings: RnaSettings) -> Optional[str]:
    """Identify, extract and classify the rRNA gene of one genome.

    Parameters
    ----------
    job : RnaJob
        The genome to search, and the domain to search it as.
    settings : RnaSettings
        What every genome of the run is searched with.

    @return: None where the genome was searched, its accession where nhmmer or
             blastn failed on it -- which is reported and left rather than
             taking the rest of the batch down.
    """

    output_dir = output_dir_of(job.genome_file, settings.results_dir)
    try:
        # every file of the directory, the other genes' results included, as
        # --remove has always meant; there is nothing to remove where the genome
        # has never been searched
        if settings.remove_prior and os.path.isdir(output_dir):
            remove_files_in_directory(output_dir)

        search_genome(job, settings, output_dir)
    except Exception as error:
        LOGGER.warning('warning: the {} gene of {} could not be searched for: '
                       '{}'.format(settings.rna_gene, job.accession, error))
        return job.accession

    return None


def search_genome(job: RnaJob, settings: RnaSettings, output_dir: str) -> None:
    """Run RNA over one genome in --tmp_dir and copy its results into place.

    Parameters
    ----------
    job : RnaJob
        The genome to search.
    settings : RnaSettings
        What every genome of the run is searched with.
    output_dir : str
        The genome's results directory.

    @return: None
    """

    hmm_model_prefix = RNA_MODELS[settings.rna_gene][job.domain]
    rna = RNA(settings.rna_gene, job.domain, 1)

    temp_dir = tempfile.mkdtemp(dir=settings.tmp_dir,
                                suffix=os.path.basename(job.genome_file))
    try:
        tmp_hmm_model_file = os.path.join(temp_dir, hmm_model_prefix + '.hmm')
        shutil.copy(os.path.join(HMM_DIR, hmm_model_prefix + '.hmm'),
                    tmp_hmm_model_file)

        temp_taxonomy = None
        if settings.taxonomy is not None:
            temp_taxonomy = os.path.join(temp_dir,
                                         os.path.basename(settings.taxonomy))
            shutil.copy(settings.taxonomy, temp_taxonomy)

        # HACK: create temporary genome file with dummy sequences at start to
        # ensure nhmmer can recognize that file uses DNA. This hack is currently
        # necessary as the --dna flag does not appear to work and nhmmer can fail
        # with an unknown base type error if the first sequence does not contain
        # all four nucleotides. This was an issue as of nhmmer 3.4 and impacts
        # a number of genomes including GCF_903931875.1.
        temp_genome_file = os.path.join(temp_dir,
                                        os.path.basename(job.genome_file))
        with gzip.open(temp_genome_file, 'wt') as fout:
            fout.write('>dummy\n')
            fout.write('ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n')
            for seq_id, seq in read_seq(job.genome_file):
                fout.write(f'>{seq_id}\n{seq}\n')

        temp_output_dir = os.path.join(temp_dir, os.path.basename(output_dir))
        os.mkdir(temp_output_dir)

        rna.run(temp_genome_file,
                tmp_hmm_model_file,
                settings.db,
                temp_taxonomy,
                temp_output_dir)

        install_results(temp_output_dir, output_dir, settings)
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


def install_results(temp_output_dir: str, output_dir: str,
                    settings: RnaSettings) -> None:
    """Copy one gene's results into the genome directory, the canary last.

    What an earlier search of the same gene left is removed first, so that a
    genome searched again keeps nothing its last search found and this one did
    not. The other genes' results are left as they are.

    Parameters
    ----------
    temp_output_dir : str
        Where RNA wrote the results.
    output_dir : str
        The genome's results directory.
    settings : RnaSettings
        What the gene was searched with.

    @return: None
    """

    make_sure_path_exists(output_dir)

    prefix = settings.rna_gene + '.'
    for name in os.listdir(output_dir):
        if name.startswith(prefix):
            os.remove(os.path.join(output_dir, name))

    canary = settings.rna_gene + CANARY_EXT
    for name in sorted(os.listdir(temp_output_dir)):
        if name != canary:
            shutil.copy(os.path.join(temp_output_dir, name),
                        os.path.join(output_dir, name))

    write_version_file(output_dir, NHMMER, settings.nhmmer_version)
    if os.path.exists(os.path.join(output_dir, f'{settings.rna_gene}.blastn.tsv')):
        write_version_file(output_dir, BLASTN, settings.blastn_version)

    shutil.copy(os.path.join(temp_output_dir, canary),
                os.path.join(output_dir, canary))


class RnaManagerSILVA(object):
    """Identify, extract, and taxonomically classify rRNA genes in genomes."""

    def __init__(self,
                 silva_version: str,
                 rna_path: str,
                 rna_gene: str,
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
        silva_version : str
            SILVA release classified against, e.g. 138.2.
        rna_path : str
            Directory holding a directory per SILVA release.
        rna_gene : str
            ssu, lsu_23S or lsu_5S.
        gtdb_domain_file : str
            GTDB domain report, read for the domain predicted from each genome's
            marker genes.
        taxonomy_file : str
            Standardised NCBI taxonomy, read for the domain of the genomes GTDB
            has no prediction for.
        cpus : int
            How many genomes are searched at once.
        tmp_dir : str
            Directory each genome is searched in, one directory per genome in
            flight and removed after.
        batch_size : int
            Genomes per batch.
        reclaim : bool
            Take over a batch another machine holds before its claim has expired.
        lease : float
            Seconds a claim survives without the machine holding it saying so.
        heartbeat : float
            Seconds between this machine saying so about a batch of its own.
        max_genome_size : float
            Largest genome assembly searched, in Mbp.

        @return: None
        """

        self.logger: logging.Logger = logging.getLogger('timestamp')

        check_dependencies([NHMMER, BLASTN])
        self.nhmmer_version: str = record_program_version(NHMMER)
        self.blastn_version: str = record_program_version(BLASTN)

        self.silva_version: str = silva_version
        self.rna_gene: str = rna_gene
        self.cpus: int = cpus
        self.tmp_dir: str = tmp_dir
        self.batch_size: int = batch_size
        self.reclaim: bool = reclaim
        self.lease: float = lease
        self.heartbeat: float = heartbeat
        self.max_genome_bases: int = int(max_genome_size * MBP)

        self.results_dir: str = RESULTS_DIR_FORMAT.format(silva_version)

        silva_root_path = os.path.join(rna_path, str(silva_version))
        silva_ssu_taxonomy_file = os.path.join(silva_root_path, 'silva_taxonomy.ssu.tsv')
        silva_ssu_ref_file = os.path.join(silva_root_path, f'SILVA_{silva_version}_SSURef_NR99_tax_silva.fasta')

        file_to_check = [silva_ssu_taxonomy_file, silva_ssu_ref_file]

        if self.rna_gene != 'ssu':
            silva_lsu_ref_file = os.path.join(silva_root_path, f'SILVA_{silva_version}_LSURef_tax_silva.fasta')
            silva_lsu_taxonomy_file = os.path.join(silva_root_path, 'silva_taxonomy.lsu.tsv')
            file_to_check.append(silva_lsu_ref_file)
            file_to_check.append(silva_lsu_taxonomy_file)

        for item in file_to_check:
            if not os.path.exists(item):
                self.logger.error(f'{item} does not exist')
                sys.exit(-1)

        self.db: Optional[str] = None
        self.taxonomy: Optional[str] = None
        if self.rna_gene == 'ssu':
            self.db = silva_ssu_ref_file
            self.taxonomy = silva_ssu_taxonomy_file
        elif self.rna_gene == 'lsu_23S':
            self.db = silva_lsu_ref_file
            self.taxonomy = silva_lsu_taxonomy_file
        elif self.rna_gene == 'lsu_5S':
            self.logger.info(
                'We currently do not curate against a 5S database, but do identify these sequences for quality assessment purposes.')

        # made here rather than by the first worker that wants it: a --tmp_dir
        # that cannot be made would otherwise be met once per genome, inside a
        # batch already claimed, and would fail every batch this machine took
        make_sure_path_exists(self.tmp_dir)

        self.domains: Dict[str, str] = read_domains(gtdb_domain_file,
                                                    taxonomy_file)

    def domain(self, accession: str) -> str:
        """Which domain's HMM this genome is searched with.

        Parameters
        ----------
        accession : str
            Genome accession, as the genome_dirs file names it.

        @return: 'ar' or 'bac' as RNA spells them, bacterial for a genome of no
                 known domain.
        """

        if domain_of(self.domains, accession) == DOMAIN_ARCHAEA:
            return RNA_ARCHAEA

        return RNA_BACTERIA

    def run_dir(self, out_dir: str) -> str:
        """Where the batches of this gene and SILVA version live under --out_dir.

        Parameters
        ----------
        out_dir : str
            Output directory of the run.

        @return: <out_dir>/<results directory>/<rna_gene>.
        """

        return os.path.join(out_dir, self.results_dir, self.rna_gene)

    def run(self,
            gtdb_genome_path_file: str,
            out_dir: str,
            all_genomes: bool = False,
            remove_prior: bool = False) -> bool:
        """Search every genome of a release for the rRNA gene of this run.

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
            Search every genome again, whatever its canary says.
        remove_prior : bool
            Empty the results directory of each genome before it is searched,
            the other genes' results included.

        @return: True where every batch this machine took finished, False where
                 one failed and is left to a later run.
        """

        run_dir = self.run_dir(out_dir)
        batches = plan_batches(gtdb_genome_path_file, run_dir, self.batch_size,
                               LAYOUT, self.logger)

        if all_genomes:
            self.logger.warning(
                'warning: --all searches every genome again, so batches already '
                'finished are done again too. Without it a finished batch is '
                'skipped.')

        settings = RnaSettings(rna_gene=self.rna_gene,
                               db=self.db,
                               taxonomy=self.taxonomy,
                               results_dir=self.results_dir,
                               tmp_dir=self.tmp_dir,
                               nhmmer_version=self.nhmmer_version,
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

            # the batch has its own log from here, since this is where anything
            # happens to it and every machine of a run writes its own --log
            with batch_log(batch_dir, self.logger, LAYOUT):
                self.logger.info('{}: starting.'.format(label))
                try:
                    with Heartbeat(os.path.join(batch_dir, RUNNING_CANARY),
                                   self.heartbeat):
                        counts = self.search_batch(batch_dir, settings, all_genomes)
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

        self.aggregate(batches, run_dir)

        # a batch that failed has already said why, in its own log and in its
        # FAILED file, and a run of several machines over days should not end by
        # printing a traceback out of the argparse frames
        if failed:
            self.logger.error(
                '{:,} batch(es) failed; they are the directories holding a FAILED '
                'file and are retried by running the command again.'.format(failed))
            return False

        return True

    def search_batch(self, batch_dir: str, settings: RnaSettings,
                     all_genomes: bool = False) -> BatchCounts:
        """Search one batch, and record the genomes that could not be searched.

        The genomes are taken from the batch's own batchfile, so the work asks
        about the genomes the batch was cut from rather than about whatever a
        genome_dirs file says now. Two passes, because deciding is cheap and
        searching is not.

        Parameters
        ----------
        batch_dir : str
            Batch directory.
        settings : RnaSettings
            What every genome of the run is searched with.
        all_genomes : bool
            Search every genome of the batch again, whatever its canary says.

        @return: what the batch came to, which run() records in its canary.
        """

        rows = read_batchfile(batchfile_path(batch_dir, LAYOUT))

        # a genome whose sequences are not there has nothing to search; it is
        # named and left rather than stopping the other ten thousand of the batch
        present, missing = split_by_fasta(rows, STAT_THREADS)
        not_searched = [(accession, REASON_NO_GENOMIC_FASTA) for accession in missing]

        jobs = [RnaJob(accession, genome_file, self.domain(accession))
                for genome_file, accession in present]

        unknown = [job.accession for job in jobs
                   if domain_of(self.domains, job.accession) is None]
        if unknown:
            self.logger.warning(
                'warning: {:,} genome(s) of this batch have no domain in either '
                'the GTDB domain file or the NCBI taxonomy and are searched as '
                'bacteria, e.g. {}.'.format(
                    len(unknown), ', '.join(sorted(unknown)[:3])))

        if all_genomes:
            to_search, already_searched = jobs, 0
        else:
            parser = functools.partial(rna_parser, rna_gene=self.rna_gene,
                                       results_dir=self.results_dir)
            with mp.Pool(processes=self.cpus) as pool:
                decided = list(tqdm(pool.imap_unordered(parser, jobs),
                                    total=len(jobs), unit='genome', ncols=100,
                                    leave=False, desc='Checking rRNA results'))
            to_search = [job for job in decided if job is not None]
            already_searched = len(jobs) - len(to_search)

        to_search, too_large = split_by_genome_size(
            to_search, lambda job: os.path.dirname(job.genome_file),
            self.max_genome_bases, self.logger, STAT_THREADS)
        not_searched.extend((job.accession, REASON_GENOME_TOO_LARGE)
                            for job, _ in too_large)

        self.logger.info(
            '{:,} genome(s) require {} identification; {:,} already have '
            'results.'.format(len(to_search), self.rna_gene, already_searched))

        not_searched.extend(self.search_genomes(to_search, settings))

        not_searched.sort()
        write_table(not_searched, os.path.join(batch_dir, NOT_SEARCHED_NAME),
                    header=NOT_SEARCHED_HEADER)

        if not_searched:
            self.logger.warning(
                'warning: {:,} genome(s) of this batch were not searched: {}.'.format(
                    len(not_searched),
                    '; '.join('{:,} {}'.format(count, reason)
                              for reason, count
                              in sorted(tally_reasons(not_searched).items()))))

        searched = len(to_search) - sum(1 for _, reason in not_searched
                                        if reason == REASON_SEARCH_FAILED)

        return BatchCounts(searched=searched,
                           already_searched=already_searched,
                           not_searched=len(not_searched))

    def search_genomes(self, jobs: Sequence[RnaJob],
                       settings: RnaSettings) -> List[Tuple[str, str]]:
        """Search the genomes of a batch that need it.

        Parameters
        ----------
        jobs : sequence of RnaJob
            The genomes to search.
        settings : RnaSettings
            What every genome of the run is searched with.

        @return: (accession, reason) for each genome the search failed on.

        Raises
        ------
        RuntimeError
            Every one of several genomes failed, which is nhmmer or blastn not
            working on this machine rather than a batch of difficult genomes. One
            genome failing on its own is not: a batch whose last unsearched
            genome is a bad one would otherwise fail identically on every retry.
        """

        if not jobs:
            return []

        worker = functools.partial(rna_worker, settings=settings)
        with mp.Pool(processes=self.cpus) as pool:
            results = list(tqdm(pool.imap_unordered(worker, jobs),
                                total=len(jobs), unit='genome', ncols=100,
                                desc='Identifying {} genes'.format(self.rna_gene)))

        failures = [(accession, REASON_SEARCH_FAILED)
                    for accession in results if accession is not None]

        if len(jobs) > 1 and len(failures) == len(jobs):
            raise RuntimeError(
                'The {} search failed on every one of the {:,} genome(s) it was '
                'given.'.format(self.rna_gene, len(jobs)))

        return failures

    def aggregate(self, batches: Sequence[str], run_dir: str) -> None:
        """Report the release, once every batch has succeeded.

        Written only when they all have, so that the file at the top of the
        directory is either the whole release or absent, and never a part of it
        that reads like the whole.

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

        path = os.path.join(run_dir, NOT_SEARCHED_RELEASE_NAME)
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
            'Release: {:,} genome(s) were searched for {} genes, {:,} already had '
            'results.'.format(totals['searched'], self.rna_gene,
                              totals['already_searched']))

        if written:
            self.logger.warning(
                'warning: {:,} genome(s) of the release were not searched and are '
                'named in {}.'.format(written, path))
        else:
            self.logger.info(
                'Every genome of the release was searched; wrote {} with no '
                'rows.'.format(path))

    @staticmethod
    def update_silva(ssu_ref_file: str, lsu_ref_file: str, output_dir: str) -> None:
        """Update SILVA reference files.

        A static method since it needs none of what a run of rna_silva is set up
        with: update_silva built an RnaManagerSILVA with no SILVA directory to
        call it on, and the constructor raised a TypeError joining None to a path.

        Parameters
        ----------
        ssu_ref_file : str
            SILVA SSURef NR99 FASTA.
        lsu_ref_file : str
            SILVA LSURef FASTA.
        output_dir : str
            Directory the taxonomy files are written to.

        @return: None
        """

        with open(os.path.join(output_dir, 'silva_taxonomy.ssu.tsv'), 'w') as fout, \
                open(ssu_ref_file) as fin:
            for line in fin:
                if line[0] == '>':
                    line_split = line[1:].strip().split(' ', 1)

                    seq_id = line_split[0]
                    taxonomy = line_split[1]
                    fout.write(f'{seq_id}\t{taxonomy}\n')

        with open(os.path.join(output_dir, 'silva_taxonomy.lsu.tsv'), 'w') as fout, \
                open(lsu_ref_file) as fin:
            for line in fin:
                if line[0] == '>':
                    line_split = line[1:].strip().split(' ', 1)

                    seq_id = line_split[0]
                    taxonomy = line_split[1]
                    fout.write(f'{seq_id}\t{taxonomy}\n')
