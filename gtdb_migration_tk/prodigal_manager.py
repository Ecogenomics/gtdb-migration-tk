import os
import csv
import gzip
import logging
import multiprocessing as mp
import shutil
from typing import Dict, List, NamedTuple, Optional, Sequence, Tuple

from tqdm import tqdm

from gtdb_migration_tk.batching import (CLAIM_LEASE_SECONDS,
                                        DEFAULT_BATCH_SIZE, HEARTBEAT_SECONDS,
                                        RUNNING_CANARY, STATE_SUCCESS,
                                        SUCCESS_CANARY,
                                        STAT_THREADS, BatchLayout, Heartbeat,
                                        age_phrase, batch_log, batch_state,
                                        batchfile_path, claim_batch, claim_age,
                                        concatenate, fail_batch, finish_batch,
                                        plan_batches, read_batchfile,
                                        read_canary, release_claim,
                                        split_by_fasta, tally_reasons,
                                        write_table)
from gtdb_migration_tk.biolib_lite.checksum import sha256_rb
from gtdb_migration_tk.biolib_lite.external.execute import check_dependencies
from gtdb_migration_tk.biolib_lite.external.prodigal import Prodigal, ProdigalTask
from gtdb_migration_tk.ncbi_utils import GENOMIC_FASTA_EXT
from gtdb_migration_tk.utils.common import (TT_SUMMARY_TABLE,
                                            read_translation_table_summary)



# Columns of the correction file --tt_override takes. It is written by hand
# against a genome that was predicted wrongly, so it is two columns and no more.
OVERRIDE_GENOME = 'genome_id'
OVERRIDE_TABLE = 'translation_table'

# What a genome's prodigal_translation_table.tsv names as the source of its
# table, so that a corrected genome says on disk that it was corrected.
SOURCE_PREDICTED = 'predicted by gTranslate'
SOURCE_OVERRIDE = 'specified by --tt_override'

# What this command calls the files of its own batches. The batchfile is the
# record of which genomes a batch is; --out_dir holds nothing else of the run's
# results, the called genes going into the genome directories as they always
# have. There is no older name to look for: this command has never had batches
# before.
BATCHFILE_NAME = 'prodigal_batchfile.tsv.gz'
BATCH_LOG_NAME = 'prodigal.log'
LAYOUT = BatchLayout(batchfiles=(BATCHFILE_NAME,), log=BATCH_LOG_NAME)

# The genomes of a batch that came out of it with no genes, and the same gathered
# for the release. Three things can leave a genome uncalled and they are
# different failures: no table was predicted for it, its genomic FASTA is not
# where the release says it is, or Prodigal was run over it and produced nothing.
# The next command needs to know which genomes have no proteins, and it should
# not have to look in 135 directories to find out.
NOT_CALLED_NAME = 'not_called.tsv'
NOT_CALLED_RELEASE_NAME = 'prodigal_not_called.tsv'
NOT_CALLED_HEADER = ('genome_id', 'reason')
REASON_NO_TABLE = 'no_translation_table'
REASON_NO_FASTA = 'no_genomic_fasta'
REASON_FAILED = 'prodigal_failed'

# Genomes Prodigal would not call in single mode and called in meta mode instead,
# and the same gathered for the release. Their genes are called from Prodigal's
# precalculated parameters rather than from a model trained on the genome, which
# is a fact about those proteins that nothing downstream could otherwise recover:
# the proteome looks like any other. The reason each one fell back is carried with
# it, since 'too many regions of N' is a draft assembly full of gaps and says so.
META_FALLBACK_NAME = 'meta_fallback.tsv'
META_FALLBACK_RELEASE_NAME = 'prodigal_meta_fallback.tsv'
META_FALLBACK_HEADER = ('genome_id', 'reason')

# What a genome's own prodigal_translation_table.tsv calls the mode, on the line
# written only for a genome whose genes were called in meta mode after single mode
# refused it. A genome called as asked gets no such line, so the file of every
# other genome of the release is what it has always been.
MODE_FIELD = 'prodigal_mode'
MODE_META = 'meta'


def has_proteins(aa_gene_file: str) -> bool:
    """Whether a genome's protein FASTA holds anything at all.

    Read rather than stat: the file is gzipped, and a gzip of nothing is still
    forty-odd bytes of header and trailer. Only the first block is decompressed.

    Parameters
    ----------
    aa_gene_file : str
        Gzipped protein FASTA of one genome.

    @return: True where it decompresses to at least one byte. A file that is not
             there, or that cannot be read as a gzip, holds no proteins either.
    """

    try:
        with gzip.open(aa_gene_file, 'rb') as handle:
            return bool(handle.read(1))
    except OSError:
        return False


class BatchCounts(NamedTuple):
    """What calling the genes of one batch came to.

    Recorded in the batch's SUCCESS canary, because the release totals are added
    up from the batches and a machine that ran the last batch has called the
    genes of none of the others.
    """

    called: int
    already_called: int
    not_called: int


# How many genomes without a translation table are named in the warning. The
# whole list is trans_table's gtranslate_no_prediction.tsv, which is where to
# read it; this is enough to recognise the kind of thing being left out.
MISSING_TABLES_LOGGED = 5


def read_translation_tables(trans_table_file: str,
                            tt_override_file: Optional[str] = None
                            ) -> Tuple[Dict[str, int], Dict[str, str]]:
    """Read the translation table to call each genome's genes under.

    The tables are what trans_table predicted, read from the summary gTranslate
    writes; the corrections are a two-column file written by hand against the
    genomes that were predicted wrongly, and they replace a prediction rather than
    being weighed against it. A correction for a genome that was never predicted
    is kept: a genome can be corrected into the run, and refusing it would mean a
    genome no classifier could answer for could never be called at all.

    Parameters
    ----------
    trans_table_file : str
        gtranslate.translation_table_summary.tsv, as trans_table writes it.
    tt_override_file : str
        Corrections, with genome_id and translation_table columns, or None.

    @return: (tables, sources), the table to use for each genome and where that
             table came from, for the genome's own record on disk.
    """

    tables, sources = {}, {}
    for genome, row in read_translation_table_summary(trans_table_file).items():
        try:
            tables[genome] = int(row[TT_SUMMARY_TABLE])
        except (KeyError, TypeError, ValueError):
            continue
        sources[genome] = SOURCE_PREDICTED

    if tt_override_file:
        with open(tt_override_file) as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                genome = (row.get(OVERRIDE_GENOME) or '').strip()
                if not genome:
                    continue
                tables[genome] = int(row[OVERRIDE_TABLE])
                sources[genome] = SOURCE_OVERRIDE

    return tables, sources


class ProdigalManager(object):
    """Call genes with Prodigal for the genomes of a release.

    The translation table of each genome is given rather than worked out here:
    trans_table predicts it with gTranslate and compares it against the table NCBI
    declares, so nothing in this module reads a GFF or weighs a coding density any
    more. It is handed the genome_dirs file of the release and the predictions, and
    writes into each genome's own directory.
    """

    def __init__(self,
                 tmp_dir: str = '/tmp/',
                 cpus: int = 1,
                 batch_size: int = DEFAULT_BATCH_SIZE,
                 reclaim: bool = False,
                 lease: float = CLAIM_LEASE_SECONDS,
                 heartbeat: float = HEARTBEAT_SECONDS) -> None:
        """Initialization.

        Prodigal is checked for here rather than when it is first called, so a
        run over half a million genomes fails at once if the binary is absent
        instead of after the pool has walked the whole release.

        Parameters
        ----------
        tmp_dir : str
            Directory Prodigal writes to before the results are moved into each
            genome directory. It holds one set of outputs per genome in flight.
        cpus : int
            Number of genomes gene called at once, and the size of the pool that
            decides which genomes need it.
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

        self.tmp_dir = tmp_dir
        self.cpus = cpus
        self.batch_size = batch_size
        self.reclaim = reclaim
        self.lease = lease
        self.heartbeat = heartbeat

        check_dependencies(['prodigal'])

        self.logger = logging.getLogger('timestamp')

    def run(self,
            gtdb_genome_path_file: str,
            trans_table_file: str,
            out_dir: str,
            tt_override_file: Optional[str] = None,
            all_genomes: bool = False) -> bool:
        """Call genes for every genome of a release that still needs them.

        The translation table of each genome is given rather than chosen: trans_table
        predicted it with gTranslate, and --tt_override corrects the predictions that
        were wrong. Prodigal is handed that table, so it calls the genes once instead
        of calling them under tables 4 and 11 and keeping whichever coded more of the
        genome, which is a rule that cannot express table 25 at all.

        The release is cut into batches under --out_dir and a batch is claimed
        before it is worked on, so several machines can be pointed at one --out_dir
        and will divide the release between them. The batching is the same
        machinery trans_table uses, in batching.py. What --out_dir holds is the
        state of the run and nothing else: the genes go into each genome's own
        prodigal/ directory, as they always have, which is why two machines on
        different batches never write to the same place.

        A genome with no table is not called. gTranslate returns no prediction for a
        few genomes of a release -- eight of r237's 1.35M, some of which have no
        genomic FASTA to predict from at all -- and there is nothing to call their
        genes under: Prodigal choosing a table by coding density is the very thing
        the summary is handed over to prevent. They are recorded and left, and
        --tt_override is how one is given a table and called after all.

        A summary covering NO genome of the release stops the run before any batch
        is claimed. That is the wrong file rather than a few unpredictable genomes,
        and carrying on would call nothing at all and report that the run had
        finished.

        Parameters
        ----------
        gtdb_genome_path_file : str
            genome_dirs file of the release, accession and genome directory per line.
        trans_table_file : str
            gtranslate.translation_table_summary.tsv.gz, as trans_table writes it.
        out_dir : str
            Directory the batches and the state of the run are written to.
        tt_override_file : str
            Corrections to those predictions, or None.
        all_genomes : bool
            Call genes for every genome, discarding the previous results, rather than
            only for those that need it (default = False).

        @return: True where every batch this machine took finished, False where one
                 failed and is left to a later run.
        """

        tables, sources = read_translation_tables(trans_table_file, tt_override_file)
        corrected = sum(1 for source in sources.values() if source == SOURCE_OVERRIDE)
        self.logger.info('Read translation tables for {:,} genomes from {}{}.'.format(
            len(tables), trans_table_file,
            ', {:,} of them corrected by {}'.format(corrected, tt_override_file)
            if tt_override_file else ''))

        batches = self.plan_batches(gtdb_genome_path_file, out_dir, tables)

        if all_genomes:
            self.logger.warning(
                'warning: --all_genomes discards the results of every genome and '
                'calls its genes again, so batches already finished are done again '
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
                        counts = self.call_batch(batch_dir, tables, sources,
                                                 all_genomes)
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
                             called=counts.called,
                             already_called=counts.already_called,
                             not_called=counts.not_called)
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

    def plan_batches(self,
                     gtdb_genome_path_file: str,
                     out_dir: str,
                     tables: Dict[str, int]) -> List[str]:
        """Settle which genomes are in which batch, once for every machine.

        The release is checked against the predictions BEFORE any batch is cut, so
        that a --trans_table for another release is met here rather than by every
        batch in turn reporting that it called nothing.

        Parameters
        ----------
        gtdb_genome_path_file : str
            genome_dirs file of the release.
        out_dir : str
            Output directory of the run.
        tables : dict
            Accession to translation table, as read_translation_tables() gave them.

        @return: paths of the batch directories, in batch order.
        """

        covered = 0
        total = 0
        with open(gtdb_genome_path_file) as handle:
            for line in handle:
                if line.strip():
                    total += 1
                    if line.split('\t')[0] in tables:
                        covered += 1

        if total and not covered:
            raise RuntimeError(
                'None of the {:,} genomes of this release has a translation table. '
                'That file is for another release, or trans_table has not been run '
                'over this one.'.format(total))

        return plan_batches(gtdb_genome_path_file, out_dir, self.batch_size,
                            LAYOUT, self.logger)


    def call_batch(self,
                   batch_dir: str,
                   tables: Dict[str, int],
                   sources: Dict[str, str],
                   all_genomes: bool = False) -> BatchCounts:
        """Call the genes of one batch, and record the genomes that got none.

        The genomes are taken from the batch's own batchfile, so the work asks
        about the genomes the batch was cut from rather than about whatever a
        genome_dirs file says now.

        Two passes, because deciding is cheap and calling is not: the first sorts
        the batch into genomes whose proteins are already there and vouched for
        and genomes that are not, the second calls what is left. A release is
        mostly carried over from the one before, so the first pass is what keeps a
        run proportional to the genomes that are actually new -- and doing it per
        batch rather than over the release means a rerun re-reads the proteins of
        one batch at a time instead of all 1.35M before it starts.

        Parameters
        ----------
        batch_dir : str
            Batch directory.
        tables : dict
            Accession to the translation table its genes are to be called under.
        sources : dict
            Accession to where that table came from.
        all_genomes : bool
            Discard what is there and call every genome of the batch again.

        @return: what the batch came to, which run() records in its canary.
        """

        rows = read_batchfile(batchfile_path(batch_dir, LAYOUT))

        # a genome with no table has nothing to call its genes under, and the
        # coding density rule is the very thing the summary is handed over to
        # prevent being fallen back on
        not_called = [(accession, REASON_NO_TABLE) for _, accession in rows
                      if accession not in tables]
        with_table = [(fasta, accession) for fasta, accession in rows
                      if accession in tables]

        present, missing = split_by_fasta(with_table, STAT_THREADS)
        not_called += [(accession, REASON_NO_FASTA) for accession in missing]

        # named for what it asks, which is not what split_by_fasta() above asks:
        # that one stats the genomic FASTA going in, this one reads the called
        # genes already sitting there. Two bars a batch saying the same thing
        # would leave a reader unable to tell which pass they were watching.
        work = [(accession, os.path.dirname(fasta), all_genomes)
                for fasta, accession in present]
        with mp.Pool(processes=self.cpus) as pool:
            decided = list(tqdm(pool.imap_unordered(self.prodigal_parser, work),
                                total=len(work), unit='genome', ncols=100,
                                leave=False, desc='Checking called genes'))

        genome_paths = [answer for answer in decided if answer != ('null', 'null')]
        already_called = len(work) - len(genome_paths)
        self.logger.info(
            '{:,} genome(s) require gene calling; {:,} already have valid '
            'Prodigal results.'.format(len(genome_paths), already_called))

        called, fell_back = self.run_prodigal(genome_paths, tables, sources)

        # written whether or not anything fell back, so that a batch in which
        # nothing did says so rather than leaving a reader to wonder
        write_table(sorted(fell_back.items()),
                    os.path.join(batch_dir, META_FALLBACK_NAME),
                    header=META_FALLBACK_HEADER)

        # a genome the wrapper was given that has no proteins afterwards was
        # tried and produced nothing; it is named rather than left to be found by
        # the next command
        not_called += [(accession, REASON_FAILED) for accession, _ in genome_paths
                       if accession not in called]

        not_called.sort()
        write_table(not_called, os.path.join(batch_dir, NOT_CALLED_NAME),
                    header=NOT_CALLED_HEADER)

        if not_called:
            self.logger.warning(
                'warning: {:,} genome(s) of this batch have no called genes: '
                '{}.'.format(len(not_called),
                             '; '.join('{:,} {}'.format(count, reason)
                                       for reason, count
                                       in sorted(tally_reasons(not_called).items()))))

        return BatchCounts(called=len(called),
                           already_called=already_called,
                           not_called=len(not_called))

    def aggregate(self, batches: Sequence[str], out_dir: str) -> None:
        """Report the release, once every batch has succeeded.

        Written only when they all have, so that the file at the top of the output
        directory is either the whole release or absent, and never a part of it
        that reads like the whole.

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
                '{:,} of {:,} batch(es) are done; the release files are '
                'written once they all are.'.format(
                    len(batches) - len(unfinished), len(batches)))
            return

        path = os.path.join(out_dir, NOT_CALLED_RELEASE_NAME)
        written = concatenate(
            [os.path.join(batch, NOT_CALLED_NAME) for batch in batches], path)

        totals = {'called': 0, 'already_called': 0}
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
            'Release: {:,} genome(s) had their genes called, {:,} already had '
            'valid Prodigal results.'.format(totals['called'],
                                             totals['already_called']))

        if written:
            self.logger.warning(
                'warning: {:,} genome(s) of the release have no called genes and '
                'are named in {}.'.format(written, path))
        else:
            self.logger.info(
                'Every genome of the release has called genes; wrote {} with no '
                'rows.'.format(path))

        fallback_path = os.path.join(out_dir, META_FALLBACK_RELEASE_NAME)
        fell_back = concatenate(
            [os.path.join(batch, META_FALLBACK_NAME) for batch in batches],
            fallback_path)

        if fell_back:
            self.logger.warning(
                'warning: {:,} genome(s) of the release were called in Prodigal\'s '
                'meta mode, single mode having refused them, and are named in '
                '{}.'.format(fell_back, fallback_path))
        else:
            self.logger.info(
                'Every genome of the release was called in the mode it was asked '
                'for; wrote {} with no rows.'.format(fallback_path))

    def prodigal_parser(self, data: Tuple[str, str, bool]) -> Tuple[str, str]:
        """Decide whether one genome still needs its genes called.

        Called in a worker of the pool run() opens, once per genome, which is why
        the arguments arrive as a single tuple: imap_unordered() passes one item.
        The genome is not processed here; the answer is a work list for
        run_prodigal(), so the decision costs a stat of the protein file and, at
        most, one pass over it.

        Under all_genomes the answer is always yes, and the existing prodigal directory is
        removed first so that the rerun cannot leave a genome holding some files
        from this run and some from the last: Prodigal's output is moved in file by
        file, and the marker results and their version-free symlinks written beside
        it are not rewritten at all.

        Otherwise a genome is skipped only if its proteins are there AND vouched
        for. The sha256 is taken of the DECOMPRESSED stream, as run_prodigal()
        wrote it, so the genome is judged by the amino acids and not by the gzip
        container, which recompressing or copying about can change on its own. A
        missing .sha256 or one that disagrees means the genes are called again:
        neither says the proteins are wrong, but neither shows them to be right,
        and the cost of recalling a genome is small against a release built on
        proteins nothing can account for.

        Parameters
        ----------
        data : tuple
            (gid, gpath, all_genomes) for one genome: its accession, its genome
            directory, and whether the run was told to redo every genome.

        @return: (gid, gpath) if the genes are to be called, else the sentinel
                 ('null', 'null'), which run() drops from the work list.
        """

        gid, gpath, all_genomes = data
        
        if all_genomes:
            prodigal_dir = os.path.join(gpath, 'prodigal')
            if os.path.exists(prodigal_dir):
                shutil.rmtree(prodigal_dir)
            return (gid, gpath)
        else:
            aa_gene_file = os.path.join(gpath, 'prodigal', gid + '_protein.faa.gz')
            checksum_file = aa_gene_file[0:-3] + '.sha256'

            if os.path.exists(aa_gene_file) and os.path.exists(checksum_file):
                # a proteome with nothing in it is not a result, however exactly
                # the digest beside it agrees. Prodigal leaves an empty file where
                # it failed, and da39a3ee... is the digest of nothing, so the two
                # agreed and the genome was called valid and skipped -- by this
                # run and by the five releases before it. GCA_000722275.1 has
                # carried an empty proteome since 2020 that way, never called
                # again and never named as failed
                if not has_proteins(aa_gene_file):
                    return (gid, gpath)

                # verify checksum
                checksum = sha256_rb(gzip.GzipFile(fileobj=open(aa_gene_file, 'rb')))
                cur_checksum = open(checksum_file).readline().strip()
                if checksum == cur_checksum:
                    return ('null', 'null')

            return (gid, gpath)

    def run_prodigal(self,
                     genome_paths: List[Tuple[str, str]],
                     tables: Dict[str, int],
                     sources: Dict[str, str]) -> List[str]:
        """Call genes for the genomes run() decided need them, and file the results.

        The table of each genome is passed to Prodigal, which then calls the genes
        once under it rather than twice and choosing.

        The paths are settled here, so each worker writes its genome's results into
        that genome's own prodigal/ directory under the names the rest of the
        toolkit expects, <accession>_protein.{faa,fna,gff}.gz, and writes the digest
        of the DECOMPRESSED protein file beside it as prodigal_parser() reads it
        back. tmp_dir is where a worker makes its own scratch directory, holding
        one genome's uncompressed gene calls until that genome is done.

        A genome whose genomic FASTA is missing or empty is warned about and left
        out of the work list rather than stopping the run: a release is built from
        whatever the mirror holds, and one unreadable genome should not cost the
        gene calling of every genome behind it.

        Parameters
        ----------
        genome_paths : list of tuple
            (accession, genome directory) for each genome to call genes for, as
            prodigal_parser() returned them.
        tables : dict
            Accession to the translation table its genes are to be called under.
        sources : dict
            Accession to where that table came from, recorded with it.

        @return: (the accessions that have proteins afterwards, the accessions
                 called in meta mode to what single mode said about them). A genome
                 the wrapper was given that has none was tried and produced
                 nothing, and call_batch() names it rather than leaving it to be
                 found by the next command.
        """

        self.logger.info(
            'Determining genomic file for each of the {:,} genomes.'.format(
                len(genome_paths)))

        # the paths are settled here because what the files are called is GTDB's
        # convention, not the wrapper's, and the wrapper writes them where they
        # belong rather than into a scratch directory for this process to move
        tasks = []
        for gid, gpath in tqdm(genome_paths, ncols=100, leave=False,
                               desc='Locating genomes'):
            if gid == 'null':
                continue
            assembly_id = os.path.basename(os.path.normpath(gpath))

            genome_file = os.path.join(gpath, assembly_id + GENOMIC_FASTA_EXT)
            if not os.path.exists(genome_file):
                self.logger.warning(
                    'Genomic file appears to be missing: %s' % genome_file)
                continue
            if os.stat(genome_file).st_size == 0:
                self.logger.warning(
                    'Genomic file appears to be empty: %s' % genome_file)
                continue

            prodigal_path = os.path.join(gpath, 'prodigal')
            aa_gene_file = os.path.join(prodigal_path, gid + '_protein.faa.gz')
            tasks.append(ProdigalTask(
                genome_id=gid,
                genome_file=genome_file,
                translation_table=tables[gid],
                aa_gene_file=aa_gene_file,
                nt_gene_file=os.path.join(prodigal_path, gid + '_protein.fna.gz'),
                gff_file=os.path.join(prodigal_path, gid + '_protein.gff.gz'),
                checksum_file=aa_gene_file[0:-3] + '.sha256',
                tmp_root=self.tmp_dir,
                called_genes=False,
                meta=False,
                closed_ends=False))

        # the genomes whose FASTA was found, not the genomes asked about: a
        # release short of a FASTA reported the full count and then called
        # genes for fewer, with nothing saying so
        self.logger.info('Running Prodigal on {:,} genomes.'.format(len(tasks)))
        summary_stats, refused = Prodigal(cpus=self.cpus).run(tasks)

        for gid, reason in sorted(refused.items()):
            self.logger.warning(
                'Prodigal would not call {}: {}'.format(gid, reason))

        # every genome of a batch refused is Prodigal not working on this machine
        # rather than a batch of difficult genomes, and the batch is failed for a
        # later run to take. One genome refused on its own is not: a batch whose
        # last uncalled genome is a bad one would fail identically on every retry
        if len(tasks) > 1 and len(refused) == len(tasks):
            raise RuntimeError(
                'Prodigal would not call any of the {:,} genome(s) of this batch; '
                'the first said: {}'.format(
                    len(tasks), sorted(refused.items())[0][1]))

        # everything else a genome's directory holds was written by the worker
        # that called its genes; this is the one record the wrapper has no
        # business knowing about
        self.logger.info('Recording the translation table of each genome.')
        called = []
        fell_back = {}
        for task in tasks:
            # the proteins are what the genome is carried into the release by, so
            # they are what says the genome was called and not the exit status of
            # a worker that may have written nothing. A genome Prodigal refused
            # still HAS a protein file -- the empty one the last release gave it,
            # which is what had it skipped here in the first place -- so what is
            # asked is whether there are proteins in it
            if task.genome_id in refused or not has_proteins(task.aa_gene_file):
                continue

            called.append(task.genome_id)
            stats = summary_stats[task.genome_id]
            table_file = os.path.join(os.path.dirname(task.aa_gene_file),
                                      'prodigal_translation_table.tsv')
            with open(table_file, 'w') as handle:
                handle.write('{}\t{}\t{}\n'.format(
                    'best_translation_table',
                    stats.best_translation_table,
                    sources.get(task.genome_id, SOURCE_PREDICTED)))

                # only for a genome that fell back, so that the file of every
                # other genome of the release is byte for byte what it was
                if stats.meta_fallback:
                    fell_back[task.genome_id] = stats.meta_fallback
                    handle.write('{}\t{}\t{}\n'.format(
                        MODE_FIELD, MODE_META,
                        'single mode failed: ' + stats.meta_fallback))

        if fell_back:
            self.logger.warning(
                'warning: {:,} genome(s) were called in Prodigal\'s meta mode, '
                'single mode having refused them.'.format(len(fell_back)))

        return called, fell_back
