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

"""align_marker_genes: the aligned marker genes GTDB's trees are built from.

It was `gtdb -r power realign_updated_genomes` of the gtdb package until 0.1.66
(PowerUserManager.RealignNCBIgenomes() and AlignedMarkerManager), which read the
genomes from the database, the files from the paths the gtdb package's Config
named, and wrote each genome's rows on its own connection. This command takes
the release's genome_dirs file, as the toolkit's other commands do, and writes
each batch of genomes in one transaction.

WHAT A GENOME'S ROWS ARE
For each marker of the sets asked for (--marker_set_ids), the gene of the
genome's top-hit table naming the marker with the highest bitscore, the first
of equal ones, is aligned to the marker's HMM, and only the HMM's match states
are kept (the 'x' columns of hmmalign's #=GC RF line), so the sequence is as
long as the HMM. multiple_hits says more than one gene named the marker,
hit_number how many, unique_genes how many distinct sequences they were. A
marker no gene names is written as a row of gaps as long as the HMM: a genome is
not expected to have every marker, incomplete or not, and the row says the
marker was looked for and is absent where no row would say it was never looked
for.

ONE HMMALIGN A MARKER A BATCH
The gtdb package ran hmmalign once for each marker of each genome, some 70M
processes over r237's 303,097 new genomes. hmmalign aligns each sequence to the
profile on its own, and only the insert columns of its output depend on the
other sequences of the file, which are the columns not kept; so the genes of a
batch's genomes are aligned to a marker by one hmmalign. Over 200 r237 genomes
and 908 genes of five markers, every alignment is the one hmmalign gives the gene
alone.

WHICH GENOMES
--new_genomes are the NCBI genomes of the database with no row in aligned_markers
for any marker of the sets: those update_db added, and those whose sequences
changed, whose rows update_db removes. A genome with rows for some markers of
the sets is not one of them. --all_genomes is every NCBI genome of the database,
its rows of those markers written again. The gtdb package also required
gtdb_representative to be set, which update_propagated_tax sets for every genome.

MISSING FILES
A genome without called proteins (prodigal/<gid>_protein.faa.gz) cannot be
aligned: it is passed over, counted, and listed in MISSING_PROTEINS_NAME in
--out_dir. One with proteins and no top-hit table, or whose top-hit table names a
gene its proteins do not hold, fails its batch: hmmsearch or top_hit has not been
run for it, or was run on other proteins.

BATCHES AND TRANSACTIONS
The genomes are cut into batches under --out_dir as batching.py cuts them for
the other long commands, so a run stopped part way, or several machines, carry on
from the batches finished. Each batch's rows are written in one transaction
(INSERT ... ON CONFLICT DO UPDATE) as the batch ends, and SUCCESS written after
it commits: a genome has all its rows or none, and a batch stopped between the
two is done again, its rows written again as they were.
"""

import io
import logging
import multiprocessing as mp
import os
import shutil
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from typing import Dict, List, NamedTuple, Optional, Sequence, Set, Tuple

from gtdb_migration_tk.batching import (CLAIM_LEASE_SECONDS, HEARTBEAT_SECONDS, RUNNING_CANARY,
                                        STATE_SUCCESS, BatchLayout, Heartbeat, age_phrase, batch_log,
                                        batch_state, claim_age, claim_batch, fail_batch, finish_batch,
                                        plan_batches, read_batchfile, read_canary, read_genome_dirs,
                                        release_claim, write_table)
from gtdb_migration_tk.biolib_lite.common import make_sure_path_exists
from gtdb_migration_tk.biolib_lite.external.execute import check_dependencies
from gtdb_migration_tk.config import MARKER_DIR_SUFFIX, PFAM_VERSION, TIGRFAM_VERSION
from gtdb_migration_tk.database_configuration import GenomeDatabaseConnectionFTPUpdate
from gtdb_migration_tk.utils.common import open_text, protein_fasta, record_program_version, write_version_file

BATCHFILE_NAME = 'align_marker_genes_batchfile.tsv.gz'
BATCH_LOG_NAME = 'align_marker_genes.log'
LAYOUT = BatchLayout(batchfiles=(BATCHFILE_NAME,), log=BATCH_LOG_NAME)

# Genomes a batch, which is also a transaction: about 233,000 rows for bac120
# and ar122, a few hundred MB of sequence held while the batch is aligned.
DEFAULT_BATCH_SIZE = 1000

# the genomes without called proteins, of a batch and, gathered, of the run
MISSING_PROTEINS_NAME = 'missing_protein_file.tsv'
MISSING_PROTEINS_HEADER = ('genome_id',)

# what the database calls the two marker databases (marker_databases.external_id_prefix),
# the directory of a genome's prodigal/ their top hits are in, and the version
# config.py says their HMMs are of
MARKER_DATABASES = {'PFAM': ('pfam', PFAM_VERSION), 'TIGR': ('tigrfam', TIGRFAM_VERSION)}

# genomes.genome_source_id of a user's genome; every other genome is NCBI's
USER_GENOME_SOURCE = 1

HMMALIGN = 'hmmalign'
HMM_MAGIC = b'HMMER'
GAP = '-'
RF_LINE = '#=GC RF'
MATCH_STATE = 'x'

# the batch's rows, copied into this temporary table and written in one statement
NEW_ROWS_TABLE = 'gtdb_new_aligned_markers'
ALIGNED_COLUMNS = ('genome_id', 'marker_id', 'sequence', 'multiple_hits', 'evalue', 'bitscore',
                   'hit_number', 'unique_genes')

# how many genomes or problems an error names
EXAMPLES = 10


class AlignmentError(ValueError):
    """The markers or genomes cannot be aligned as the database and files stand."""


class Marker(NamedTuple):
    """One marker gene of the sets aligned, as the database holds it."""

    db_id: int
    accession: str
    hmm: str
    size: int
    database: str


class ChosenGene(NamedTuple):
    """The gene a genome's top hits give a marker."""

    sequence: str
    evalue: str
    bitscore: str
    multiple_hits: bool
    hit_number: int
    unique_genes: int


def read_fasta_gz(path: str) -> Dict[str, str]:
    """The sequences of a (gzipped) protein FASTA, without Prodigal's trailing '*'.

    @return: gene ID -> sequence.
    """

    sequences: Dict[str, List[str]] = {}
    current: Optional[List[str]] = None
    with open_text(path) as handle:
        for line in handle:
            if line.startswith('>'):
                current = sequences.setdefault(line[1:].split(None, 1)[0], [])
            elif current is not None:
                current.append(line.strip())
    # Prodigal ends each protein with '*', which some tools downstream of the
    # alignment do not take
    return {gene: ''.join(parts).rstrip('*') for gene, parts in sequences.items()}


def tophit_file(prodigal_dir: str, accession: str, database: str) -> str:
    """A genome's top-hit table of a marker database, as top_hit wrote it.

    @return: e.g. prodigal/pfam_33.1_lite/<gid>_pfam_33.1_lite_tophit.tsv.gz.
    """

    name = '{}_{}'.format(MARKER_DATABASES[database][0], MARKER_DIR_SUFFIX[MARKER_DATABASES[database][0]])
    return os.path.join(prodigal_dir, name, '{}_{}_tophit.tsv.gz'.format(accession, name))


def read_tophits(path: str, wanted: Set[str]) -> Dict[str, List[Tuple[str, str, float]]]:
    """Each wanted marker's hits in a top-hit table, in the order the table gives them.

    A table is a header line, then a gene and its hits, 'family,e-value,bitscore',
    separated by ';' where a gene hit more than one family.

    @return: marker accession -> [(gene, e-value as written, bitscore), ...].
    """

    hits: Dict[str, List[Tuple[str, str, float]]] = {}
    with open_text(path) as handle:
        handle.readline()
        for line in handle:
            columns = line.rstrip('\n').split('\t')
            if len(columns) < 2:
                continue
            for hit in columns[1].split(';'):
                marker, evalue, bitscore = hit.split(',')
                if marker in wanted:
                    hits.setdefault(marker, []).append((columns[0], evalue, float(bitscore)))
    return hits


def choose_genes(hits: Dict[str, List[Tuple[str, str, float]]],
                 proteins: Dict[str, str]) -> Dict[str, ChosenGene]:
    """The gene each marker is aligned from, as the gtdb package chose it.

    The highest bitscore wins and the first of equal ones is kept.

    @return: marker accession -> ChosenGene.

    Raises
    ------
    KeyError
        A gene of the hits is not among the proteins.
    """

    chosen = {}
    for marker, marker_hits in hits.items():
        best = marker_hits[0]
        for hit in marker_hits[1:]:
            if hit[2] > best[2]:
                best = hit
        sequences = [proteins[gene] for gene, _, _ in marker_hits]
        chosen[marker] = ChosenGene(sequence=proteins[best[0]], evalue=best[1], bitscore=str(best[2]),
                                    multiple_hits=len(marker_hits) > 1, hit_number=len(marker_hits),
                                    unique_genes=len(set(sequences)))
    return chosen


def read_genome(job: Tuple[str, str, Tuple[Tuple[str, Tuple[str, ...]], ...]]):
    """A genome's chosen gene for each marker; run on a worker process.

    Parameters
    ----------
    job : (accession, protein FASTA, ((database, (marker accession, ...)), ...))

    @return: (accession, marker accession -> ChosenGene), or (accession, None)
             where the genome has no protein file.

    Raises
    ------
    AlignmentError
        A top-hit table is missing, or names a gene the proteins do not hold.
    """

    accession, proteins_path, by_database = job
    try:
        proteins = read_fasta_gz(proteins_path)
    except FileNotFoundError:
        return accession, None

    chosen: Dict[str, ChosenGene] = {}
    prodigal_dir = os.path.dirname(proteins_path)
    for database, markers in by_database:
        path = tophit_file(prodigal_dir, accession, database)
        try:
            hits = read_tophits(path, set(markers))
        except FileNotFoundError:
            raise AlignmentError('{} has called proteins and no top-hit table {}; run hmmsearch and top_hit '
                                 'for it.'.format(accession, path))
        try:
            chosen.update(choose_genes(hits, proteins))
        except KeyError as exc:
            raise AlignmentError('{} names gene {} in {}, which is not among its proteins {}; the top hits '
                                 'are of other proteins.'.format(accession, exc.args[0], path, proteins_path))
    return accession, chosen


def aligned_match_states(lines: Sequence[str], names: Set[str]) -> Dict[str, str]:
    """Each sequence's match-state columns from hmmalign's Pfam-format output.

    @return: sequence name -> its residues and deletions ('-') at the columns
             the #=GC RF line marks 'x'.

    Raises
    ------
    AlignmentError
        The output has no #=GC RF line, or not every sequence.
    """

    rows: Dict[str, str] = {}
    reference = None
    for line in lines:
        if line.startswith(RF_LINE):
            reference = line.rsplit(None, 1)[-1]
        elif line and not line.startswith('#') and not line.startswith('//'):
            name, _, aligned = line.partition(' ')
            if name in names:
                rows[name] = aligned.strip()
    if reference is None:
        raise AlignmentError('hmmalign gave no {} line.'.format(RF_LINE))
    missing = names - set(rows)
    if missing:
        raise AlignmentError('hmmalign gave no alignment of {:,} sequence(s), e.g. {}.'.format(
            len(missing), ', '.join(sorted(missing)[:EXAMPLES])))
    columns = [i for i, state in enumerate(reference) if state == MATCH_STATE]
    return {name: ''.join(aligned[i] for i in columns) for name, aligned in rows.items()}


def align_marker(marker: Marker, sequences: Dict[str, str], work_dir: str) -> Dict[str, str]:
    """Align the genes of a batch to one marker's HMM, by one hmmalign.

    Parameters
    ----------
    marker : Marker
    sequences : dict
        Sequence name, unique in the batch, -> protein sequence.
    work_dir : str
        Directory the FASTA handed to hmmalign is written in.

    @return: sequence name -> its alignment, marker.size long.

    Raises
    ------
    AlignmentError
        hmmalign failed, or gave an alignment of another length.
    """

    fasta = os.path.join(work_dir, '{}.faa'.format(marker.db_id))
    with open(fasta, 'w') as handle:
        for name, sequence in sequences.items():
            handle.write('>{}\n{}\n'.format(name, sequence))
    proc = subprocess.run([HMMALIGN, '--outformat', 'Pfam', marker.hmm, fasta],
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    os.remove(fasta)
    if proc.returncode != 0:
        raise AlignmentError('hmmalign failed on {} ({}), exit status {}: {}'.format(
            marker.accession, marker.hmm, proc.returncode, ' '.join(proc.stderr.split())))
    aligned = aligned_match_states(proc.stdout.splitlines(), set(sequences))
    wrong = [name for name, sequence in aligned.items() if len(sequence) != marker.size]
    if wrong:
        raise AlignmentError('hmmalign aligned {:,} gene(s) to {} with other than its {} match states.'.format(
            len(wrong), marker.accession, marker.size))
    return aligned


def copy_text(value) -> str:
    """A value as COPY's text format reads it.

    @return: str, NULL as \\N.
    """

    if value is None:
        return '\\N'
    if isinstance(value, bool):
        return 't' if value else 'f'
    return str(value).replace('\\', '\\\\').replace('\t', '\\t').replace('\n', '\\n').replace('\r', '\\r')


class MarkerAlignmentManager(object):
    """Align the marker genes of the database's genomes and write them to aligned_markers."""

    def __init__(self, database, cpus: int = 1, batch_size: int = DEFAULT_BATCH_SIZE, tmp_dir: str = '/tmp',
                 reclaim: bool = False, lease: float = CLAIM_LEASE_SECONDS,
                 heartbeat: float = HEARTBEAT_SECONDS) -> None:
        """Initialization.

        Parameters
        ----------
        database : dict
            How to reach the database, as utils.common.database_keywords() gives it.
        cpus : int
            Genomes read, and hmmaligns run, at once.
        batch_size : int
            Genomes per batch, and so per transaction.
        tmp_dir : str
            Directory the FASTA handed to hmmalign is written in.
        reclaim, lease, heartbeat
            As the other batched commands take them.

        @return: None
        """

        self.logger = logging.getLogger('timestamp')
        self.cpus = max(1, cpus)
        self.batch_size = batch_size
        self.tmp_dir = tmp_dir
        self.reclaim = reclaim
        self.lease = lease
        self.heartbeat = heartbeat
        self.hmmalign_version = ''

        self.temp_con = GenomeDatabaseConnectionFTPUpdate.GenomeDatabaseConnectionFTPUpdate(database)
        self.temp_con.MakePostgresConnection()
        self.temp_cur = self.temp_con.cursor()

    def read_markers(self, marker_set_ids: Sequence[int]) -> List[Marker]:
        """The markers of the sets, each checked against its HMM.

        @return: the markers, by accession.

        Raises
        ------
        AlignmentError
            A set the database does not hold or that holds no markers, or a
            marker whose HMM is missing, is not an HMM, is of another length than
            the database says, or of another version than config.py names; every
            one is named.
        """

        self.temp_cur.execute('SELECT id, name FROM marker_sets WHERE id = ANY(%s)', (list(marker_set_ids),))
        sets = dict(self.temp_cur.fetchall())
        unknown = sorted(set(marker_set_ids) - set(sets))
        if unknown:
            raise AlignmentError('The database holds no marker set {}.'.format(', '.join(map(str, unknown))))
        for set_id in marker_set_ids:
            self.logger.info('Marker set {}: {}.'.format(set_id, sets[set_id]))

        self.temp_cur.execute(
            'SELECT DISTINCT m.id, m.id_in_database, m.marker_file_location, m.size, md.external_id_prefix '
            'FROM markers m JOIN marker_set_contents s ON s.marker_id = m.id '
            'JOIN marker_databases md ON md.id = m.marker_database_id WHERE s.set_id = ANY(%s)',
            (list(marker_set_ids),))
        markers = sorted((Marker(*row) for row in self.temp_cur.fetchall()), key=lambda m: m.accession)
        if not markers:
            raise AlignmentError('Marker set(s) {} hold no markers.'.format(', '.join(map(str, marker_set_ids))))

        problems = []
        for marker in markers:
            if marker.database not in MARKER_DATABASES:
                problems.append('{} is of marker database {}, which this command does not know'.format(
                    marker.accession, marker.database))
                continue
            name, version = MARKER_DATABASES[marker.database]
            if '/{}/'.format(version) not in marker.hmm:
                problems.append('{} is {}, not of {} {} as config.py says'.format(
                    marker.accession, marker.hmm, name, version))
                continue
            try:
                with open(marker.hmm, 'rb') as handle:
                    head = handle.read(4096)
            except OSError as exc:
                problems.append('{}: {}'.format(marker.accession, exc))
                continue
            if not head.startswith(HMM_MAGIC):
                problems.append('{} ({}) is not an HMM'.format(marker.accession, marker.hmm))
                continue
            leng = [line.split()[1] for line in head.decode('ascii', 'replace').splitlines()
                    if line.startswith('LENG ')]
            if not leng or int(leng[0]) != marker.size:
                problems.append('{} has {} match states in {}, where the database says {}'.format(
                    marker.accession, leng[0] if leng else 'no', marker.hmm, marker.size))
        if problems:
            raise AlignmentError('{:,} marker(s) cannot be aligned to: {}.'.format(
                len(problems), '; '.join(problems)))

        counts = {database: sum(1 for m in markers if m.database == database) for database in MARKER_DATABASES}
        self.logger.info('Aligning {:,} markers: {}.'.format(
            len(markers), ', '.join('{:,} {}'.format(n, MARKER_DATABASES[d][0]) for d, n in counts.items())))
        return markers

    def select_genomes(self, markers: Sequence[Marker], all_genomes: bool) -> Dict[str, int]:
        """The NCBI genomes to align, by accession.

        @return: genomes.id_at_source -> genomes.id.
        """

        if all_genomes:
            self.temp_cur.execute('SELECT g.id_at_source, g.id FROM genomes g WHERE g.genome_source_id <> %s',
                                  (USER_GENOME_SOURCE,))
        else:
            self.temp_cur.execute(
                'SELECT g.id_at_source, g.id FROM genomes g WHERE g.genome_source_id <> %s AND NOT EXISTS '
                '(SELECT 1 FROM aligned_markers am WHERE am.genome_id = g.id AND am.marker_id = ANY(%s))',
                (USER_GENOME_SOURCE, [m.db_id for m in markers]))
        genomes = dict(self.temp_cur.fetchall())
        self.temp_con.rollback()
        return genomes

    def run(self, marker_set_ids: Sequence[int], all_genomes: bool, gtdb_genome_path_file: str,
            out_dir: str) -> bool:
        """Align the marker genes of the genomes asked for, a batch at a time.

        Parameters
        ----------
        marker_set_ids : sequence of int
            The marker sets whose markers are aligned (marker_sets.id).
        all_genomes : bool
            Every NCBI genome of the database (--all_genomes), else only those
            with no row for any of the markers (--new_genomes).
        gtdb_genome_path_file : str
            genome_dirs file of the release, which says where each genome is.
        out_dir : str
            Directory the batches and the run's files are written to.

        @return: True where every batch this machine took finished.

        Raises
        ------
        AlignmentError
            The markers cannot be aligned to, or a genome to be aligned is not in
            the genome_dirs file; nothing is aligned.
        """

        check_dependencies([HMMALIGN])
        make_sure_path_exists(self.tmp_dir)
        make_sure_path_exists(out_dir)
        markers = self.read_markers(marker_set_ids)
        self.hmmalign_version = record_program_version(HMMALIGN)

        genomes = self.select_genomes(markers, all_genomes)
        self.logger.info('{:,} NCBI genome(s) of the database to align ({}).'.format(
            len(genomes), '--all_genomes' if all_genomes else 'with no row for any of the markers'))

        located = {accession for accession, _ in read_genome_dirs(gtdb_genome_path_file)}
        unlocated = sorted(set(genomes) - located)
        if unlocated:
            path = os.path.join(out_dir, 'not_in_genome_dirs.tsv')
            write_table([(g,) for g in unlocated], path, ('genome_id',))
            raise AlignmentError('{:,} genome(s) to align are not in {}, e.g. {}; it is not the genome_dirs '
                                 'file of the database\'s release. Every one is in {}.'.format(
                                     len(unlocated), gtdb_genome_path_file, ', '.join(unlocated[:EXAMPLES]),
                                     path))
        if not genomes:
            self.logger.info('There are no genomes to align.')
            return True

        # one set of batches for each choice of sets and genomes: the SUCCESS of a
        # batch of other markers, or of new genomes only, says nothing of this run's
        state_dir = os.path.join(out_dir, 'marker_sets_{}_{}'.format(
            '_'.join(str(i) for i in sorted(set(marker_set_ids))), 'all' if all_genomes else 'new'))
        batches = plan_batches(gtdb_genome_path_file, state_dir, self.batch_size, LAYOUT, self.logger,
                               genome_file=protein_fasta, accessions=genomes)

        done, held, failed = 0, 0, 0
        for index, batch_dir in enumerate(batches, start=1):
            label = 'Batch {:,} of {:,} ({})'.format(index, len(batches), os.path.basename(batch_dir))
            if batch_state(batch_dir) == STATE_SUCCESS:
                self.logger.info('{}: already finished, skipping.'.format(label))
                continue
            if not claim_batch(batch_dir, self.reclaim, self.lease):
                owner = read_canary(os.path.join(batch_dir, RUNNING_CANARY))
                held += 1
                self.logger.info('{}: held by {}, last heard from {}, skipping.'.format(
                    label, owner.get('host', 'another machine'),
                    age_phrase(claim_age(os.path.join(batch_dir, RUNNING_CANARY)))))
                continue

            with batch_log(batch_dir, self.logger, LAYOUT):
                self.logger.info('{}: starting.'.format(label))
                try:
                    with Heartbeat(os.path.join(batch_dir, RUNNING_CANARY), self.heartbeat):
                        aligned, missing = self.align_batch(batch_dir, markers, genomes)
                except KeyboardInterrupt:
                    self.temp_con.rollback()
                    release_claim(batch_dir)
                    self.logger.error('{}: interrupted; nothing of it was written, and the claim is given '
                                      'up.'.format(label))
                    raise
                except Exception as exc:
                    self.temp_con.rollback()
                    failed += 1
                    fail_batch(batch_dir, str(exc))
                    self.logger.error('{}: failed, nothing of it written, and will be retried by a later run: '
                                      '{}'.format(label, exc))
                    continue
                finish_batch(batch_dir, aligned=aligned, missing_protein_file=missing)
                done += 1
                self.logger.info('{}: done, {:,} genome(s) aligned, {:,} without a protein file.'.format(
                    label, aligned, missing))

        self.logger.info('{:,} batch(es) finished here, {:,} held by another machine, {:,} failed.'.format(
            done, held, failed))
        self.gather_missing(batches, out_dir)

        if failed:
            self.logger.error('{:,} batch(es) failed; they are the directories holding a FAILED file and are '
                              'retried by running the command again.'.format(failed))
            return False
        return True

    def align_batch(self, batch_dir: str, markers: Sequence[Marker], genomes: Dict[str, int]) -> Tuple[int, int]:
        """Align one batch's genomes and write their rows, in one transaction.

        @return: (genomes aligned, genomes without a protein file).
        """

        rows = read_batchfile(os.path.join(batch_dir, BATCHFILE_NAME))
        by_database = tuple((database, tuple(m.accession for m in markers if m.database == database))
                            for database in MARKER_DATABASES if any(m.database == database for m in markers))
        jobs = [(accession, proteins, by_database) for proteins, accession in rows]

        chosen: Dict[str, Dict[str, ChosenGene]] = {}
        missing: List[str] = []
        with mp.Pool(processes=self.cpus) as pool:
            for accession, genes in pool.imap_unordered(read_genome, jobs, chunksize=8):
                if genes is None:
                    missing.append(accession)
                else:
                    chosen[accession] = genes
        missing.sort()
        write_table([(g,) for g in missing], os.path.join(batch_dir, MISSING_PROTEINS_NAME),
                    MISSING_PROTEINS_HEADER)

        # one hmmalign for each marker, of every gene of the batch chosen for it,
        # each named by its genome's place in the batch
        order = sorted(chosen)
        alignments: Dict[str, Dict[str, str]] = {}
        work_dir = tempfile.mkdtemp(prefix='align_marker_genes.', dir=self.tmp_dir)
        try:
            def align(marker: Marker):
                sequences = {str(i): chosen[accession][marker.accession].sequence
                             for i, accession in enumerate(order) if marker.accession in chosen[accession]}
                return marker.accession, align_marker(marker, sequences, work_dir) if sequences else {}

            with ThreadPoolExecutor(max_workers=self.cpus) as executor:
                for accession, aligned in executor.map(align, markers):
                    alignments[accession] = aligned
        finally:
            shutil.rmtree(work_dir, ignore_errors=True)

        written = self.write_rows(markers, genomes, order, chosen, alignments)
        self.temp_con.commit()
        write_version_file(batch_dir, HMMALIGN, self.hmmalign_version)
        self.logger.info('Wrote {:,} rows of aligned_markers for {:,} genome(s).'.format(written, len(order)))
        return len(order), len(missing)

    def write_rows(self, markers: Sequence[Marker], genomes: Dict[str, int], order: Sequence[str],
                   chosen: Dict[str, Dict[str, ChosenGene]], alignments: Dict[str, Dict[str, str]]) -> int:
        """Write a batch's rows of aligned_markers, in the caller's transaction.

        Every marker of every genome is given a row, a marker the genome has no
        gene for a row of gaps.

        @return: the number of rows written.
        """

        def lines():
            for i, accession in enumerate(order):
                for marker in markers:
                    gene = chosen[accession].get(marker.accession)
                    if gene is None:
                        values = (genomes[accession], marker.db_id, GAP * marker.size, False,
                                  None, None, None, None)
                    else:
                        values = (genomes[accession], marker.db_id, alignments[marker.accession][str(i)],
                                  gene.multiple_hits, gene.evalue, gene.bitscore, gene.hit_number,
                                  gene.unique_genes)
                    yield '\t'.join(copy_text(v) for v in values) + '\n'

        columns = ', '.join(ALIGNED_COLUMNS)
        self.temp_cur.execute('DROP TABLE IF EXISTS {}'.format(NEW_ROWS_TABLE))
        self.temp_cur.execute('CREATE TEMPORARY TABLE {} (LIKE aligned_markers INCLUDING DEFAULTS) '
                              'ON COMMIT DROP'.format(NEW_ROWS_TABLE))
        self.temp_cur.copy_expert('COPY {} ({}) FROM STDIN'.format(NEW_ROWS_TABLE, columns),
                                  io.StringIO(''.join(lines())))
        self.temp_cur.execute(
            'INSERT INTO aligned_markers ({0}) SELECT {0} FROM {1} ON CONFLICT (genome_id, marker_id) DO UPDATE SET '
            '{2}'.format(columns, NEW_ROWS_TABLE,
                         ', '.join('{0} = EXCLUDED.{0}'.format(c) for c in ALIGNED_COLUMNS[2:])))
        written = self.temp_cur.rowcount
        self.temp_cur.execute('DROP TABLE {}'.format(NEW_ROWS_TABLE))
        return written

    def gather_missing(self, batches: Sequence[str], out_dir: str) -> None:
        """Gather every batch's genomes without a protein file into MISSING_PROTEINS_NAME.

        @return: None
        """

        missing: List[str] = []
        for batch_dir in batches:
            path = os.path.join(batch_dir, MISSING_PROTEINS_NAME)
            if os.path.exists(path):
                with open(path) as handle:
                    handle.readline()
                    missing.extend(line.strip() for line in handle if line.strip())
        path = os.path.join(out_dir, MISSING_PROTEINS_NAME)
        write_table([(g,) for g in sorted(missing)], path, MISSING_PROTEINS_HEADER)
        if missing:
            self.logger.warning('{:,} genome(s) have no protein file and could not be aligned; each is listed in '
                                '{}.'.format(len(missing), path))
        else:
            self.logger.info('Every genome of the batches had a protein file.')
