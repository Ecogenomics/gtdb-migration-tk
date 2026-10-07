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
import datetime
from collections import Counter, deque
from concurrent.futures import ThreadPoolExecutor
import logging
import multiprocessing as mp
import ntpath
from typing import Callable, Dict, Iterable, Iterator, List, NamedTuple, Optional, Sequence, Set, TextIO, Tuple, TypeVar

from tqdm import tqdm

from gtdb_migration_tk.batching import write_table
from gtdb_migration_tk.biolib_lite.common import make_sure_path_exists, get_num_lines
from gtdb_migration_tk.biolib_lite.seq_io import read_fasta
from gtdb_migration_tk.genometk_lite.metadata_genes import MetadataGenes
from gtdb_migration_tk.genometk_lite.metadata_nucleotide import MetadataNucleotide


# What a genome needed and did not have. The metadata of a release is generated
# while the gene calling of its last genomes is still finishing, and a genome
# whose files are not there yet is a straggler rather than a broken run: it is
# named here and passed over, where check_file_exists() would have ended the
# command on the first one and thrown away the million genomes already done. The
# file is what makes that safe -- a warning in a log of a million lines is not
# something anyone will find, and what the reader needs is the list.
MISSING_FILES_NAME = 'metadata_missing_files.tsv'
MISSING_FILES_HEADER = ('genome_id', 'missing', 'path')

# The two files a genome is expected to have, and what the log calls each. The
# genomic FASTA comes from the mirror; the GFF is the one prodigal wrote.
MISSING_GENOMIC_FASTA = 'genomic_fasta'
MISSING_PROTEIN_GFF = 'protein_gff'
MISSING_LABEL: Dict[str, str] = {MISSING_GENOMIC_FASTA: 'genomic FASTA',
                                 MISSING_PROTEIN_GFF: 'called genes (GFF)'}

# One genome as _producer() receives it, the tuple being what
# mp.Pool.imap_unordered can carry: accession, its genomic FASTA, and the GFF of
# the genes prodigal called on it.
MetadataJob = Tuple[str, str, str]

# One file a genome should have had: accession, which of the two it is, and
# where it was looked for. A row of MISSING_FILES_NAME.
MissingFile = Tuple[str, str, str]

# What gene() writes into a genome directory, and what is removed from one whose
# GFF has gone. create_metadata_tables() reads metadata.genome_gene.tsv wherever
# it is, so a table an earlier run wrote outlives the genes it was calculated
# from: r237 called 15 genomes as empty proteomes, generate_metadata wrote each
# a protein count of 0, and the patch that called them again left 3 with no GFF
# at all -- and with the 0, which create_metadata_tables() would have loaded as
# though Prodigal had found no genes, where nothing is known.
GENE_METADATA_FILES = ('metadata.genome_gene.tsv', 'metadata.genome_gene.desc.tsv')

# What _producer() hands back: the files the genome was missing, and whether gene
# metadata of an earlier run was removed for want of a GFF.
ProducerResult = Tuple[List[MissingFile], bool]

# The tables create_tables writes that give a genome a row only where it has
# something to put in one. A genome without the file a table is read from is
# passed over by its parser without a word, which across a million genomes is
# the one thing worth knowing about the run, so each is counted and the run
# ends by saying how many (MetadataTable.log_summary()). The three *_count.tsv
# tables give every genome a row, 0 where nothing was found, and are not here.
NT_TABLE = 'metadata_nt.tsv'
GENE_TABLE = 'metadata_gene.tsv'
LSU_5S_TABLE = 'metadata_lsu_5S.tsv'
TRNA_TABLE = 'metadata_trna_count.tsv'

# The tables every genome of a release should have a row in. Each genome has a
# genomic FASTA, so genomic_metadata gives each its nucleotide metadata, and one
# without it is a genome that command did not get to; the summary warns of it.
# The others are rightly missing for some genomes: no proteins called, no rRNA
# gene found.
EVERY_GENOME_TABLES = frozenset({NT_TABLE})


# How many genomes create_tables reads at once is --cpus, which defaults to 1 as
# it does for every command. Its work is opening a dozen small files in each
# genome directory, one NFS round trip after another with nothing computed
# between them, so read one genome at a time it waits on the file server for
# nearly all of a run. Threads overlap the round trips. Measured over r237's
# genome directories, on fresh samples so that nothing was cached: 75 ms a
# genome on one thread (a day and more for 1.35M genomes), 22 ms on 4, and
# 9-10 ms on 8, where it levels off -- 16, 32 and 64 threads were no faster, the
# file server being what is waited on by then. So --cpus 8, as for
# ncbi_genome_sync --nfs-jobs and list_genomes: past the knee, and no more load
# on a server others are using than buys anything.

# How many genomes are queued for each thread at once: enough that none sits
# idle between one genome and the next, and few enough that a release is not
# held in memory as futures before the first genome is written.
QUEUE_DEPTH = 4

T = TypeVar('T')
R = TypeVar('R')


def ordered_map(func: Callable[[T], R], items: Iterable[T], threads: int,
                depth: int = QUEUE_DEPTH) -> Iterator[R]:
    """func over items on threads, the results yielded in the order of items.

    At most threads * depth items are in hand at once, the next submitted as
    the oldest is yielded. An exception raised by func is raised here, for the
    item it was raised on, and the items queued behind it are cancelled.

    Parameters
    ----------
    func : callable
        Applied to each item, on a worker thread.
    items : iterable
        Read as the work is queued, so it may be a generator over a file.
    threads : int
        Worker threads; fewer than one is taken as one.
    depth : int
        Items queued per thread.

    @return: iterator over func(item), in the order of items.
    """

    threads = max(1, threads)
    with ThreadPoolExecutor(max_workers=threads) as pool:
        pending = deque()
        try:
            for item in items:
                pending.append(pool.submit(func, item))
                if len(pending) >= threads * depth:
                    yield pending.popleft().result()
            while pending:
                yield pending.popleft().result()
        finally:
            for future in pending:
                future.cancel()


class TableRow(NamedTuple):
    """One table's part of one genome, as a reader hands it back.

    The header is what the table opens with were this the first genome to
    have the file it is read from; row is the genome's line, None where the
    file named nothing to report. Both end in a newline.
    """

    header: str
    row: Optional[str]


class GenomeTables(NamedTuple):
    """Everything create_tables writes for one genome, read on a worker thread.

    tables maps each table of MetadataTable.table_sources to the genome's
    TableRow, or to None where the genome has no file that table is read
    from; the three counts are its rows of the *_count.tsv tables.
    """

    gid: str
    tables: Dict[str, Optional[TableRow]]
    ssu_count: int
    lsu_23s_count: int
    lsu_5S_count: int


class EmptyGenomeDirs(ValueError):
    """The genome_dirs file create_tables was given names no genomes."""


def taxonomy_table(prefix: str) -> str:
    """The table of one of the rRNA genes _read_taxonomy_file() reads.

    Parameters
    ----------
    prefix : str
        ssu_gg, ssu_silva or lsu_silva_23s.

    @return: the table's file name, metadata_<prefix>.tsv.
    """

    return 'metadata_{}.tsv'.format(prefix)


class MetadataTable(object):
    """Gather the metadata of every genome of a release into tables.

    Calculates nothing: genomic_metadata, rna_silva and trnascan have written
    their results into each genome directory, and this collects them into the
    ten tables update_metadata_db loads. A genome missing a file is given no
    row in that table.
    """

    def __init__(self, silva_version: str) -> None:
        """Initialization.

        Parameters
        ----------
        silva_version : str
            Version of SILVA the rRNA genes were classified against. It names
            the directory within each genome directory those results were
            written to, and must match config.SILVA_VERSION.

        @return: None
        """

        silva_folder = f'rna_silva_{silva_version}'

        # every path here is relative to one genome's directory
        self.metadata_nt_file: str = 'metadata.genome_nt.tsv'
        self.metadata_gene_file: str = 'metadata.genome_gene.tsv'
        self.ssu_gg_taxonomy_file: str = os.path.join('ssu_gg', 'ssu.taxonomy.tsv')
        self.ssu_gg_fna_file: str = os.path.join('ssu_gg', 'ssu.fna')
        self.ssu_silva_taxonomy_file: str = os.path.join(
            silva_folder, 'ssu.taxonomy.tsv')
        self.ssu_silva_fna_file: str = os.path.join(silva_folder, 'ssu.fna')
        self.ssu_silva_summary_file: str = os.path.join(
            silva_folder, 'ssu.hmm_summary.tsv')
        self.lsu_silva_23s_taxonomy_file: str = os.path.join(
            silva_folder, 'lsu_23S.taxonomy.tsv')
        self.lsu_silva_23s_fna_file: str = os.path.join(
            silva_folder, 'lsu_23S.fna')
        self.lsu_silva_23s_summary_file: str = os.path.join(
            silva_folder, 'lsu_23S.hmm_summary.tsv')

        self.lsu_5S_fna_file: str = os.path.join(silva_folder, 'lsu_5S.fna')
        self.lsu_5S_summary_file: str = os.path.join(
            silva_folder, 'lsu_5S.hmm_summary.tsv')

        # the file in a genome's directory each table's rows are read from, in
        # the order the closing summary names them
        self.table_sources: Dict[str, str] = {
            NT_TABLE: self.metadata_nt_file,
            GENE_TABLE: self.metadata_gene_file,
            taxonomy_table('ssu_gg'): self.ssu_gg_taxonomy_file,
            taxonomy_table('ssu_silva'): self.ssu_silva_taxonomy_file,
            taxonomy_table('lsu_silva_23s'): self.lsu_silva_23s_taxonomy_file,
            LSU_5S_TABLE: self.lsu_5S_fna_file,
            TRNA_TABLE: os.path.join('trna', '<gid>_trna_stats.tsv')}

        # genomes given a row in each table, and genomes without the file it
        # is read from
        self.rows: Counter = Counter()
        self.absent: Counter = Counter()

        self.logger: logging.Logger = logging.getLogger('timestamp')

    @staticmethod
    def _read_field_table(genome_id: str, path: str) -> Optional[TableRow]:
        """Read a two column field/value table generate_metadata wrote for one genome.

        metadata.genome_nt.tsv and metadata.genome_gene.tsv are both of this
        shape. The table they are gathered into is headed by the fields of the
        first genome that has one, and every genome contributes one row of
        values.

        Parameters
        ----------
        genome_id : str
            Unique identifier of genome.
        path : str
            Full path to the genome's field/value table.

        @return: the header and the genome's row, or None where the genome has
                 no such file.
        """

        try:
            with open(path) as handle:
                lines = handle.readlines()
        except FileNotFoundError:
            return None

        header = 'genome_id' + ''.join('\t' + line.split('\t')[0].strip() for line in lines) + '\n'
        row = genome_id + ''.join('\t' + line.split('\t')[1].strip() for line in lines) + '\n'
        return TableRow(header, row)

    @staticmethod
    def _read_taxonomy_file(genome_id: str,
                            metadata_taxonomy_file: str,
                            prefix: str,
                            fna_file: str,
                            summary_file: Optional[str] = None) -> Tuple[Optional[TableRow], int]:
        """Read the taxonomic information of one genome's rRNA genes.

        One method over the three rRNA tables -- ssu_gg, ssu_silva and
        lsu_silva_23s -- which differ in the prefix their fields carry and in
        whether a summary file accompanies them.

        Parameters
        ----------
        genome_id : str
            Unique identifier of genome.
        metadata_taxonomy_file : str
            Full path to file containing rRNA metadata.
        prefix : str
            Prefix to append to metadata fields.
        fna_file : str
            FASTA of the identified rRNA genes, read for the sequence of the
            hit reported.
        summary_file : str, optional
            HMM summary of the same genes, read for the length of the contig
            the hit sits on. The greengenes table has none.

        @return: the header and the genome's row, its row None where no hit
                 was reported, or None where the genome has no such table;
                 and the number of rRNA genes identified in the genome, which
                 is zero where it has no such table.
        """

        try:
            f = open(metadata_taxonomy_file)
        except FileNotFoundError:
            return None, 0

        with f:
            header_line = f.readline()  # consume header line
            headers = [prefix + '_' + x.strip().replace('ssu_', '')
                       for x in header_line.split('\t')]
            headers.append("{0}_sequence".format(prefix))
            headers.append("{0}_contig_len".format(prefix))

            if prefix == 'lsu_silva_23s':
                for n, i in enumerate(headers):
                    if i == 'lsu_silva_23s_sequence':
                        headers[n] = 'lsu_23s_sequence'
                    elif i == 'lsu_silva_23s_query_id':
                        headers[n] = 'lsu_23s_query_id'
                    elif i == 'lsu_silva_23s_length':
                        headers[n] = 'lsu_23s_length'
                    elif i == 'lsu_silva_23s_contig_len':
                        headers[n] = 'lsu_23s_contig_len'
            elif prefix == 'ssu_silva':
                for n, i in enumerate(headers):
                    if i == 'ssu_silva_sequence':
                        headers[n] = 'ssu_sequence'
                    elif i == 'ssu_silva_query_id':
                        headers[n] = 'ssu_query_id'
                    elif i == 'ssu_silva_length':
                        headers[n] = 'ssu_length'
                    elif i == 'ssu_silva_contig_len':
                        headers[n] = 'ssu_contig_len'

            header = 'genome_id' + '\t' + '\t'.join(headers) + "\n"

            # Check the CheckM headers are consistent
            split_headers = header_line.rstrip().split("\t")
            for pos in range(0, len(split_headers)):
                if split_headers[pos] == 'query_id':
                    query_id_pos = pos
                    break

            # Report hit to longest 16S rRNA gene. It is possible that
            # the HMMs identified a putative 16S rRNA gene, but that
            # there was no valid BLAST hit.
            longest_query_len = 0
            longest_ssu_hit_info = None
            identified_ssu_genes = 0
            for line in f:
                line_split = line.strip().split('\t')
                query_len = int(line_split[2])
                if query_len > longest_query_len:
                    longest_query_len = query_len
                    longest_ssu_hit_info = line_split
                    ssu_query_id = line_split[query_id_pos]

        if not longest_ssu_hit_info:
            return TableRow(header, None), identified_ssu_genes

        row = [genome_id, '\t' + '\t'.join(longest_ssu_hit_info)]
        all_genes_dict = read_fasta(fna_file, False)
        row.append('\t{0}'.format(all_genes_dict[ssu_query_id]))
        if summary_file is not None and os.path.exists(summary_file):
            with open(summary_file) as fsum:
                header_list = [x.strip() for x in fsum.readline().split('\t')]
                idx_seq = header_list.index("Sequence length")
                for line in fsum:
                    identified_ssu_genes += 1
                    sum_list = [x.strip() for x in line.split('\t')]
                    if sum_list[0] == ssu_query_id:
                        row.append("\t{0}".format(sum_list[idx_seq]))
        row.append('\n')

        return TableRow(header, ''.join(row)), identified_ssu_genes

    @staticmethod
    def _read_lsu_5S_files(accession: str,
                           fna_file: str,
                           summary_file: str) -> Tuple[Optional[TableRow], int]:
        """Read one genome's 5S LSU genes.

        The 5S genes are not classified, so there is no taxonomy table to read
        as the other rRNA genes have: the longest sequence found is reported
        from the FASTA and the HMM summary alone.

        Parameters
        ----------
        accession : str
            Unique identifier of genome.
        fna_file : str
            FASTA of the identified 5S genes.
        summary_file : str
            HMM summary of the same genes, read for the length of the contig
            each sits on.

        @return: the header and the genome's row, its row None where no gene
                 was reported, or None where no 5S sequence was identified;
                 and the number of 5S genes identified, which is zero where
                 none were.
        """

        # check if a 5S sequence was identified
        if not os.path.exists(fna_file):
            return None, 0

        header = 'genome_id\tlsu_5s_query_id\tlsu_5s_length\tlsu_5s_contig_len\tlsu_5s_sequence\n'
        seqs = read_fasta(fna_file)

        identified_genes = 0
        longest_seq = 0
        longest_seq_id = None
        longest_contig_len = None
        if os.path.exists(summary_file):
            with open(summary_file) as fsum:
                header_line = fsum.readline()  # consume header line
                header_list = [x.strip() for x in header_line.split('\t')]
                idx_seq_len = header_list.index("Sequence length")
                for line in fsum:
                    identified_genes += 1

                    line_split = list(map(str.strip, line.strip().split('\t')))
                    seq_id = line_split[0]
                    contig_len = int(line_split[idx_seq_len])
                    seq_len = len(seqs[seq_id])

                    if seq_len > longest_seq:
                        longest_seq_id = seq_id
                        longest_seq = seq_len
                        longest_contig_len = contig_len

        row = None
        if longest_seq_id:
            row = '%s\t%s\t%d\t%d\t%s\n' % (accession, longest_seq_id,
                                              longest_seq, longest_contig_len, seqs[longest_seq_id])

        return TableRow(header, row), identified_genes

    @staticmethod
    def _read_trna_file(genome_id: str, trna_file: str) -> Optional[TableRow]:
        """Read tRNA information.

        Reads the statistics tRNAscan-SE wrote for one genome, counting the
        tRNAs, the amino acids they decode and the selenocysteine tRNAs among
        them.

        Parameters
        ----------
        genome_id : str
            Unique identifier of genome.
        trna_file : str
            Full path to the genome's tRNA statistics file.

        @return: the header and the genome's row of counts, or None where the
                 genome has no such file.
        """

        try:
            with open(trna_file) as handle:
                lines = handle.readlines()
        except FileNotFoundError:
            return None

        # parse tRNA summary file
        trna_count = 0
        trna_selenocysteine_count = 0
        trna_aa_count = 0
        read_aa = False
        for line in lines:
            if line.startswith('tRNAs decoding Standard 20 AA'):
                trna_count = int(line.split(':')[1])
            elif line.startswith('Selenocysteine tRNAs (TCA)'):
                trna_selenocysteine_count = int(line.split(':')[1])
                trna_count += trna_selenocysteine_count
            elif line.startswith('Isotype / Anticodon Counts:'):
                read_aa = True

            if read_aa:
                # Parsing lines with the following format:
                # Ala   : 2	  AGC:         GGC:         CGC: 1       TGC: 1
                if ' : ' in line:
                    line_split = line.split('\t')[0]
                    line_split = list(map(str.strip, line_split.split(':')))
                    if len(line_split[0]) == 3:  # this is an amino acid
                        if int(line_split[1]) > 0:
                            trna_aa_count += 1

        return TableRow('genome_id\ttrna_count\ttrna_aa_count\ttrna_selenocysteine_count\n',
                        '%s\t%d\t%d\t%d\n' % (genome_id, trna_count, trna_aa_count,
                                              trna_selenocysteine_count))

    def _read_genome(self, genome: Tuple[str, str]) -> GenomeTables:
        """Read everything create_tables writes for one genome.

        Run on a worker thread: it reads the genome's files and nothing else,
        holding no state of the run, so that the order the tables are written
        in, and the genome each is headed by, are decided where they are
        written (create_metadata_tables()).

        Parameters
        ----------
        genome : (str, str)
            The genome's accession and directory.

        @return: the genome's part of every table.
        """

        gid, gpath = genome
        ssu_silva, ssu_count = self._read_taxonomy_file(
            gid, os.path.join(gpath, self.ssu_silva_taxonomy_file), 'ssu_silva',
            os.path.join(gpath, self.ssu_silva_fna_file),
            os.path.join(gpath, self.ssu_silva_summary_file))
        lsu_23s, lsu_23s_count = self._read_taxonomy_file(
            gid, os.path.join(gpath, self.lsu_silva_23s_taxonomy_file), 'lsu_silva_23s',
            os.path.join(gpath, self.lsu_silva_23s_fna_file),
            os.path.join(gpath, self.lsu_silva_23s_summary_file))
        lsu_5S, lsu_5S_count = self._read_lsu_5S_files(
            gid, os.path.join(gpath, self.lsu_5S_fna_file),
            os.path.join(gpath, self.lsu_5S_summary_file))
        ssu_gg, _ = self._read_taxonomy_file(
            gid, os.path.join(gpath, self.ssu_gg_taxonomy_file), 'ssu_gg',
            os.path.join(gpath, self.ssu_gg_fna_file))

        tables = {
            NT_TABLE: self._read_field_table(gid, os.path.join(gpath, self.metadata_nt_file)),
            GENE_TABLE: self._read_field_table(gid, os.path.join(gpath, self.metadata_gene_file)),
            taxonomy_table('ssu_gg'): ssu_gg,
            taxonomy_table('ssu_silva'): ssu_silva,
            taxonomy_table('lsu_silva_23s'): lsu_23s,
            LSU_5S_TABLE: lsu_5S,
            TRNA_TABLE: self._read_trna_file(gid, os.path.join(gpath, 'trna', gid + '_trna_stats.tsv'))}

        return GenomeTables(gid, tables, ssu_count, lsu_23s_count, lsu_5S_count)

    def create_metadata_tables(self, gtdb_genome_path_file: str, output_dir: str,
                               cpus: int = 1) -> None:
        """Create metadata tables.

        One pass over the release, gathering what every earlier command wrote
        into each genome directory into the ten tables the database is loaded
        from. The genomes are read cpus at a time on threads (_read_genome()),
        and written here in the order of the genome_dirs file, so the tables
        are the same whatever cpus is.

        Parameters
        ----------
        gtdb_genome_path_file : str
            genome_dirs file of the release: accession, directory, canonical
            accession, one genome per line.
        output_dir : str
            Directory the tables are written to, made if it does not exist.
        cpus : int
            Genomes read at once.

        @return: nothing; the tables are written under output_dir.

        Raises
        ------
        EmptyGenomeDirs
            The genome_dirs file names no genomes. Nothing is written: ten
            tables of no rows would replace the release's in --out_dir, and
            update_metadata_db would load them.
        """

        numlines = get_num_lines(gtdb_genome_path_file)
        if numlines == 0:
            raise EmptyGenomeDirs('{} names no genomes; no tables were written.'.format(
                gtdb_genome_path_file))

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        fouts: Dict[str, TextIO] = {table: open(os.path.join(output_dir, table), 'w')
                                    for table in self.table_sources}
        fout_ssu_silva_count = open(os.path.join(
            output_dir, 'metadata_ssu_silva_count.tsv'), 'w')
        fout_lsu_silva_23s_count = open(os.path.join(
            output_dir, 'metadata_lsu_silva_23s_count.tsv'), 'w')
        fout_lsu_5S_count = open(os.path.join(
            output_dir, 'metadata_lsu_5S_count.tsv'), 'w')

        fout_ssu_silva_count.write('%s\t%s\n' % ('genome_id', 'ssu_count'))
        fout_lsu_silva_23s_count.write(
            '%s\t%s\n' % ('genome_id', 'lsu_23s_count'))
        fout_lsu_5S_count.write('%s\t%s\n' % ('genome_id', 'lsu_5s_count'))

        def genomes(handle: TextIO) -> Iterator[Tuple[str, str]]:
            for line in handle:
                line_split = line.strip().split('\t')
                yield line_split[0], line_split[1]

        # each table is headed by the first genome, in the order of the
        # genome_dirs file, that has the file it is read from
        headed: Set[str] = set()

        self.logger.info('Gathering the metadata of {:,} genomes into {}, reading {:,} at once.'.format(
            numlines, output_dir, max(1, cpus)))
        genome_count = 0
        with open(gtdb_genome_path_file) as ggpf:
            for genome in tqdm(ordered_map(self._read_genome, genomes(ggpf), cpus),
                               ncols=100, total=numlines, smoothing=50/numlines):
                genome_count += 1

                for table, part in genome.tables.items():
                    if part is None:
                        self.absent[table] += 1
                        continue
                    if table not in headed:
                        headed.add(table)
                        fouts[table].write(part.header)
                    if part.row is not None:
                        fouts[table].write(part.row)
                        self.rows[table] += 1

                fout_ssu_silva_count.write(
                    '%s\t%d\n' % (genome.gid, genome.ssu_count))
                fout_lsu_silva_23s_count.write(
                    '%s\t%d\n' % (genome.gid, genome.lsu_23s_count))
                fout_lsu_5S_count.write(
                    '%s\t%d\n' % (genome.gid, genome.lsu_5S_count))

        for fout in fouts.values():
            fout.close()
        fout_ssu_silva_count.close()
        fout_lsu_silva_23s_count.close()
        fout_lsu_5S_count.close()

        self.log_summary(genome_count)

    def log_summary(self, genome_count: int) -> None:
        """Say how many genomes each table holds, and how many had no file.

        A parser passes over a genome without the file its table is read from,
        and nothing else says so: this is where a release finds that
        genomic_metadata, rna_silva or trnascan did not get to every genome. A
        genome with the file and still no row is one whose file named nothing
        to report, such as an rRNA table of no hits.

        Parameters
        ----------
        genome_count : int
            Genomes of the genome_dirs file, every one of them read.

        @return: nothing; a line is logged for each table of table_sources.
        """

        self.logger.info('Rows written for {:,} genomes:'.format(genome_count))
        for table, source in self.table_sources.items():
            rows = self.rows[table]
            absent = self.absent[table]
            nothing = genome_count - rows - absent

            parts = ['{:,} with a row'.format(rows)]
            if absent:
                parts.append('{:,} had no {}'.format(absent, source))
            if nothing:
                parts.append('{:,} had nothing to report'.format(nothing))
            message = '  {}: {}.'.format(table, '; '.join(parts))
            if absent and table in EVERY_GENOME_TABLES:
                self.logger.warning(message)
            else:
                self.logger.info(message)


class MetadataManager(object):
    """Create file indicating directory of each genome."""

    def __init__(self, cpus: int = 1) -> None:
        """Initialization.

        Parameters
        ----------
        cpus : int
            How many genomes have their metadata calculated at once.

        @return: None
        """

        self.cpus: int = cpus

        # a run of at least this many ambiguous bases breaks one contig from
        # the next when the scaffolds are split for the nucleotide statistics
        self.contig_break: int = 10

        self.logger: logging.Logger = logging.getLogger('timestamp')
        self.starttime: Optional[datetime.datetime] = None

########### GENERATE METADATA ######
    def generate_metadata(self, gtdb_genome_path_file: str, out_dir: str) -> None:
        """Calculate the nucleotide and gene metadata of every genome of a release.

        The results are written into each genome directory rather than gathered
        here; create_metadata_tables() is what gathers them afterwards. What
        --out_dir holds is the account of the run: which genomes could not be
        done, and why.

        Parameters
        ----------
        gtdb_genome_path_file : str
            genome_dirs file of the release: accession, directory, canonical
            accession, one genome per line.
        out_dir : str
            Directory MISSING_FILES_NAME is written to, made if it does not
            exist. No metadata is written here.

        @return: nothing; two tables and their descriptions are written into
                 each genome directory, and the genomes missing a file are
                 written to out_dir.
        """

        make_sure_path_exists(out_dir)

        self.starttime = datetime.datetime.utcnow().replace(microsecond=0)
        input_files: List[MetadataJob] = []
        with open(gtdb_genome_path_file) as ggpf:
            for line in tqdm(ggpf,total=get_num_lines(gtdb_genome_path_file)):
                gid,gpath,*_ = line.strip().split('\t')
                assembly_id = os.path.basename(os.path.normpath(gpath))

                genome_file = os.path.join(gpath, assembly_id + '_genomic.fna.gz')
                gff_file = os.path.join(gpath, 'prodigal', gid + '_protein.gff.gz')
                input_files.append((gid, genome_file, gff_file))

        # process each genome. imap_unordered rather than the vendored Parallel
        # class: a producer that raised there killed its worker silently, and the
        # run went on to report success having processed fewer genomes than it was
        # given. Here the exception reaches this loop and stops the command.
        self.logger.info('Generating metadata for {:,} genomes:'.format(
            len(input_files)))
        missing: List[MissingFile] = []
        removed = 0
        with mp.Pool(processes=self.cpus) as pool:
            for result, removed_gene in tqdm(pool.imap_unordered(self._producer, input_files),
                                             total=len(input_files), ncols=100, unit='genome'):
                # warned here rather than in the worker: several processes
                # appending to one log file interleave, and the parent is
                # reading every result anyway
                for gid, what, missing_file in result:
                    self.logger.warning('{} has no {}: {}{}'.format(
                        gid, MISSING_LABEL[what], missing_file,
                        '; removed the gene metadata an earlier run wrote'
                        if removed_gene and what == MISSING_PROTEIN_GFF else ''))
                missing.extend(result)
                removed += removed_gene

        self.report_missing(missing, len(input_files), out_dir)
        if removed:
            self.logger.warning(
                'Removed the gene metadata an earlier run wrote for {:,} genome(s) '
                'that no longer have called genes; create_tables gives them no gene '
                'row.'.format(removed))

    def report_missing(self,
                       missing: Sequence[MissingFile],
                       genome_count: int,
                       out_dir: str) -> str:
        """Name the genomes that were missing a file, once, for the release.

        Written whether or not anything was missing, so that a release with
        nothing missing says so rather than leaving the reader to wonder
        whether the run got that far.

        Parameters
        ----------
        missing : sequence of MissingFile
            What every genome was missing, gathered from the workers.
        genome_count : int
            Genomes the run was given, for the tally.
        out_dir : str
            Directory the report is written to.

        @return: the file written.
        """

        report = os.path.join(out_dir, MISSING_FILES_NAME)
        write_table(missing, report, MISSING_FILES_HEADER)

        if missing:
            genomes = len({gid for gid, _, _ in missing})
            self.logger.warning(
                '{:,} of {:,} genomes were missing a file; they are named in {}.'.format(
                    genomes, genome_count, report))
        else:
            self.logger.info(
                'Every one of {:,} genomes had both of its files; {} is empty.'.format(
                    genome_count, report))

        return report

    def _producer(self, job: MetadataJob) -> ProducerResult:
        """Process each genome.

        A genome missing a file is reported and passed over rather than ending
        the run. Which files are there decides how much of the genome can be
        done: the gene metadata needs the GFF and the genome size both, so it
        needs the two files, while the nucleotide metadata needs only the
        FASTA and is written whenever the FASTA is there. The two are separate
        files that create_metadata_tables() reads independently, so half a
        genome is worth having and is not done again when prodigal catches up.

        A genome with its FASTA and no GFF has the gene metadata of an earlier
        run removed (GENE_METADATA_FILES): it was calculated from genes the
        genome no longer has, and create_metadata_tables() would load it.

        Parameters
        ----------
        job : MetadataJob
            The genome's accession, its genomic FASTA and the GFF of the genes
            Prodigal called.

        @return: the files this genome should have had and did not, which is
                 empty for a genome that was processed in full, and whether gene
                 metadata of an earlier run was removed.
        """

        gid, genome_file, gff_file = job
        full_genome_dir, _ = ntpath.split(genome_file)

        missing: List[MissingFile] = []
        if not os.path.isfile(genome_file):
            missing.append((gid, MISSING_GENOMIC_FASTA, genome_file))
        if not os.path.isfile(gff_file):
            missing.append((gid, MISSING_PROTEIN_GFF, gff_file))

        # without the sequences there is nothing to calculate at all, and the
        # genome directory is left as it was found -- including its old log
        if any(what == MISSING_GENOMIC_FASTA for _, what, _ in missing):
            return missing, False

        # clean up old log files
        log_file = os.path.join(full_genome_dir, 'genometk.log')
        if os.path.exists(log_file):
            os.remove(log_file)

        # calculate metadata
        self.nucleotide(genome_file,full_genome_dir)
        if not missing:
            self.gene(genome_file,gff_file,full_genome_dir)
            return missing, False

        # no GFF: gene metadata already here is of genes the genome no longer has
        removed = False
        for name in GENE_METADATA_FILES:
            path = os.path.join(full_genome_dir, name)
            if os.path.exists(path):
                os.remove(path)
                removed = True

        return missing, removed

    def nucleotide(self, genome_file: str, output_dir: str) -> None:
        """Calculate metadata derived from one genome's nucleotide sequences.

        Parameters
        ----------
        genome_file : str
            The genome's genomic FASTA.
        output_dir : str
            Directory the statistics and their descriptions are written to,
            which is the genome directory.

        @return: nothing; metadata.genome_nt.tsv and its .desc.tsv are written
                 under output_dir.
        """

        # whether the file is there is _producer()'s to decide: it reports a
        # genome that has not got one rather than ending the release on it
        make_sure_path_exists(output_dir)

        meta_nuc = MetadataNucleotide()
        metadata_values, metadata_desc = meta_nuc.generate(genome_file,
                                                           self.contig_break)

        # write statistics to file
        output_file = os.path.join(output_dir, 'metadata.genome_nt.tsv')
        fout = open(output_file, 'w')
        for field in sorted(metadata_values.keys()):
            fout.write('%s\t%s\n' % (field, str(metadata_values[field])))
        fout.close()

        # write description to file
        output_file = os.path.join(output_dir, 'metadata.genome_nt.desc.tsv')
        fout = open(output_file, 'w')
        for field in sorted(metadata_desc.keys()):
            fout.write('%s\t%s\t%s\n' % (field,
                                         metadata_desc[field],
                                         type(metadata_values[field]).__name__.upper()))
        fout.close()

    def gene(self, genome_file: str, gff_file: str, output_dir: str) -> None:
        """Calculate metadata derived from one genome's called genes.

        Parameters
        ----------
        genome_file : str
            The genome's genomic FASTA.
        gff_file : str
            The GFF of the genes Prodigal called on it.
        output_dir : str
            Directory the statistics and their descriptions are written to,
            which is the genome directory.

        @return: nothing; metadata.genome_gene.tsv and its .desc.tsv are
                 written under output_dir.
        """

        # as nucleotide(): _producer() calls this only for a genome that has
        # both of its files
        make_sure_path_exists(output_dir)

        meta_genes = MetadataGenes()
        metadata_values, metadata_desc = meta_genes.generate(genome_file,
                                                                gff_file)

        # write statistics to file
        output_file = os.path.join(output_dir, 'metadata.genome_gene.tsv')
        fout = open(output_file, 'w')
        for field in sorted(metadata_values.keys()):
            fout.write('%s\t%s\n' % (field, str(metadata_values[field])))
        fout.close()

        # write description to file
        output_file = os.path.join(output_dir, 'metadata.genome_gene.desc.tsv')
        fout = open(output_file, 'w')
        for field in sorted(metadata_desc.keys()):
            fout.write('%s\t%s\t%s\n' % (field,
                                         metadata_desc[field],
                                         type(metadata_values[field]).__name__.upper()))
        fout.close()



