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
import logging
import multiprocessing as mp
import ntpath
from typing import Dict, List, Optional, Sequence, Set, TextIO, Tuple

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


class MetadataTable(object):
    """Create metadata table for all NCBI and user genomes.

    This script assumes the scripts metadata_generate.py
    and ssu.py have been run in order to create the
    required metadata. Four tables are generated which
    specific nucleotide derived, gene derived, and SSU
    derived metadata.
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

        # each output table is headed by the first genome that has anything to
        # put in it, so these say whether that genome has been seen yet
        self.write_nt_header: bool = True
        self.write_gene_header: bool = True
        self.write_trna_header: bool = True
        self.write_lsu_5S_header: bool = True

        # the rRNA tables share one method over three prefixes, so which of them
        # have been headed is kept by prefix rather than by a flag each
        self.taxonomy_headers: Set[str] = set()

    def _parse_nt(self, genome_id: str, metadata_nt_file: str, fout: TextIO) -> None:
        """Parse metadata file with information derived from nucleotide sequences.

        The file is the two column field/value table generate_metadata wrote for
        one genome. The output table is headed from the first genome that has
        one, and every genome contributes one row of values thereafter.

        Parameters
        ----------
        genome_id : str
            Unique identifier of genome.
        metadata_nt_file : str
            Full path to the genome's nucleotide metadata file.
        fout : TextIO
            Output stream to populate with metadata.

        @return: nothing; a row is written to fout, and a genome without the
                 file contributes none.
        """

        if not os.path.exists(metadata_nt_file):
            return

        if self.write_nt_header:
            self.write_nt_header = False

            fout.write('genome_id')
            for line in open(metadata_nt_file):
                line_split = line.split('\t')
                fout.write('\t' + line_split[0].strip())
            fout.write('\n')

        fout.write(genome_id)
        for line in open(metadata_nt_file):
            line_split = line.split('\t')
            fout.write('\t' + line_split[1].strip())
        fout.write('\n')

    def _parse_gene(self, genome_id: str, metadata_gene_file: str, fout: TextIO) -> None:
        """Parse metadata file with information derived from called genes.

        As _parse_nt, over the table generate_metadata wrote from the genes
        Prodigal called.

        Parameters
        ----------
        genome_id : str
            Unique identifier of genome.
        metadata_gene_file : str
            Full path to the genome's gene metadata file.
        fout : TextIO
            Output stream to populate with metadata.

        @return: nothing; a row is written to fout, and a genome without the
                 file contributes none.
        """

        if not os.path.exists(metadata_gene_file):
            return

        if self.write_gene_header:
            self.write_gene_header = False

            fout.write('genome_id')
            for line in open(metadata_gene_file):
                line_split = line.split('\t')
                fout.write('\t' + line_split[0].strip())
            fout.write('\n')
        fout.write(genome_id)
        for line in open(metadata_gene_file):
            line_split = line.split('\t')
            fout.write('\t' + line_split[1].strip())
        fout.write('\n')

    def _parse_taxonomy_file(self,
                             genome_id: str,
                             metadata_taxonomy_file: str,
                             fout: TextIO,
                             prefix: str,
                             fna_file: str,
                             summary_file: Optional[str] = None) -> int:
        """Parse metadata file with taxonomic information for rRNA genes.

        One method over the three rRNA tables -- ssu_gg, ssu_silva and
        lsu_silva_23s -- which differ in the prefix their fields carry and in
        whether a summary file accompanies them.

        Parameters
        ----------
        genome_id : str
            Unique identifier of genome.
        metadata_taxonomy_file : str
            Full path to file containing rRNA metadata.
        fout : TextIO
            Output stream to populate with metadata.
        prefix : str
            Prefix to append to metadata fields.
        fna_file : str
            FASTA of the identified rRNA genes, read for the sequence of the
            hit reported.
        summary_file : str, optional
            HMM summary of the same genes, read for the length of the contig
            the hit sits on. The greengenes table has none.

        @return: number of rRNA genes identified in the genome, which is zero
                 where the genome has no such table.
        """

        if not os.path.exists(metadata_taxonomy_file):
            return 0

        with open(metadata_taxonomy_file) as f:
            header_line = f.readline()  # consume header line
            if prefix not in self.taxonomy_headers:
                self.taxonomy_headers.add(prefix)

                fout.write('genome_id')
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

                fout.write('\t' + '\t'.join(headers) + "\n")

            # Check the CheckM headers are consistent
            split_headers = header_line.rstrip().split("\t")
            for pos in range(0, len(split_headers)):
                header = split_headers[pos]
                if header == 'query_id':
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

            if longest_ssu_hit_info:
                fout.write(genome_id)
                fout.write('\t' + '\t'.join(longest_ssu_hit_info))
                all_genes_dict = read_fasta(fna_file, False)
                sequence = all_genes_dict[ssu_query_id]
                fout.write('\t{0}'.format(sequence))
                if summary_file is not None and os.path.exists(summary_file):
                    with open(summary_file) as fsum:
                        header_line = fsum.readline()  # consume header line
                        header_list = [x.strip()
                                       for x in header_line.split('\t')]
                        idx_seq = header_list.index("Sequence length")
                        for line in fsum:
                            identified_ssu_genes += 1
                            sum_list = [x.strip() for x in line.split('\t')]
                            if sum_list[0] == ssu_query_id:
                                fout.write("\t{0}".format(sum_list[idx_seq]))

                fout.write('\n')

            return identified_ssu_genes

    def _parse_lsu_5S_files(self,
                            accession: str,
                            fout: TextIO,
                            fna_file: str,
                            summary_file: str) -> int:
        """Parse information from 5S LSU files.

        The 5S genes are not classified, so there is no taxonomy table to read
        as the other rRNA genes have: the longest sequence found is reported
        from the FASTA and the HMM summary alone.

        Parameters
        ----------
        accession : str
            Unique identifier of genome.
        fout : TextIO
            Output stream to populate with metadata.
        fna_file : str
            FASTA of the identified 5S genes.
        summary_file : str
            HMM summary of the same genes, read for the length of the contig
            each sits on.

        @return: number of 5S genes identified in the genome, which is zero
                 where none were.
        """

        # check if a 5S sequence was identified
        if not os.path.exists(fna_file):
            return 0

        # write header
        if self.write_lsu_5S_header:
            fout.write(
                'genome_id\tlsu_5s_query_id\tlsu_5s_length\tlsu_5s_contig_len\tlsu_5s_sequence\n')
            self.write_lsu_5S_header = False

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

        if longest_seq_id:
            fout.write('%s\t%s\t%d\t%d\t%s\n' % (accession, longest_seq_id,
                                                 longest_seq, longest_contig_len, seqs[longest_seq_id]))

        return identified_genes

    def _parse_trna_file(self, genome_id: str, trna_file: str, fout_trna_count: TextIO) -> None:
        """Parse tRNA information.

        Reads the statistics tRNAscan-SE wrote for one genome, counting the
        tRNAs, the amino acids they decode and the selenocysteine tRNAs among
        them.

        Parameters
        ----------
        genome_id : str
            Unique identifier of genome.
        trna_file : str
            Full path to the genome's tRNA statistics file.
        fout_trna_count : TextIO
            Output stream to populate with the counts.

        @return: nothing; a row is written to fout_trna_count, and a genome
                 without the file contributes none.
        """

        if not os.path.exists(trna_file):
            return

        # write header
        if self.write_trna_header:
            fout_trna_count.write(
                'genome_id\ttrna_count\ttrna_aa_count\ttrna_selenocysteine_count\n')
            self.write_trna_header = False

        # parse tRNA summary file
        trna_count = 0
        trna_selenocysteine_count = 0
        trna_aa_count = 0
        read_aa = False
        for line in open(trna_file):
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

        fout_trna_count.write('%s\t%d\t%d\t%d\n' % (
            genome_id, trna_count, trna_aa_count, trna_selenocysteine_count))

    def create_metadata_tables(self, gtdb_genome_path_file: str, output_dir: str) -> None:
        """Create metadata tables.

        One pass over the release, gathering what every earlier command wrote
        into each genome directory into the ten tables the database is loaded
        from.

        Parameters
        ----------
        gtdb_genome_path_file : str
            genome_dirs file of the release: accession, directory, canonical
            accession, one genome per line.
        output_dir : str
            Directory the tables are written to, made if it does not exist.

        @return: nothing; the tables are written under output_dir.
        """

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        fout_nt = open(os.path.join(output_dir, 'metadata_nt.tsv'), 'w')
        fout_gene = open(os.path.join(output_dir, 'metadata_gene.tsv'), 'w')
        fout_gg_taxonomy = open(os.path.join(
            output_dir, 'metadata_ssu_gg.tsv'), 'w')
        fout_ssu_silva_taxonomy = open(os.path.join(
            output_dir, 'metadata_ssu_silva.tsv'), 'w')
        fout_lsu_silva_23s_taxonomy = open(os.path.join(
            output_dir, 'metadata_lsu_silva_23s.tsv'), 'w')
        fout_lsu_5S = open(os.path.join(
            output_dir, 'metadata_lsu_5S.tsv'), 'w')
        fout_ssu_silva_count = open(os.path.join(
            output_dir, 'metadata_ssu_silva_count.tsv'), 'w')
        fout_lsu_silva_23s_count = open(os.path.join(
            output_dir, 'metadata_lsu_silva_23s_count.tsv'), 'w')
        fout_lsu_5S_count = open(os.path.join(
            output_dir, 'metadata_lsu_5S_count.tsv'), 'w')
        fout_trna_count = open(os.path.join(
            output_dir, 'metadata_trna_count.tsv'), 'w')

        fout_ssu_silva_count.write('%s\t%s\n' % ('genome_id', 'ssu_count'))
        fout_lsu_silva_23s_count.write(
            '%s\t%s\n' % ('genome_id', 'lsu_23s_count'))
        fout_lsu_5S_count.write('%s\t%s\n' % ('genome_id', 'lsu_5s_count'))

        # generate metadata for NCBI assemblies
        numlines = get_num_lines(gtdb_genome_path_file)
        with open(gtdb_genome_path_file) as ggpf:
            for line in tqdm(ggpf,ncols=100,total=numlines,smoothing=50/numlines):

                line_split = line.strip().split('\t')

                gid = line_split[0]
                gpath = line_split[1]
                assembly_id = os.path.basename(os.path.normpath(gpath))
                metadata_nt_file = os.path.join(
                    gpath, self.metadata_nt_file)
                self._parse_nt(
                    gid, metadata_nt_file, fout_nt)

                metadata_gene_file = os.path.join(
                    gpath, self.metadata_gene_file)
                self._parse_gene(
                    gid, metadata_gene_file, fout_gene)

                ssu_gg_taxonomy_file = os.path.join(
                    gpath, self.ssu_gg_taxonomy_file)
                ssu_gg_fna_file = os.path.join(
                    gpath, self.ssu_gg_fna_file)
                self._parse_taxonomy_file(
                    gid, ssu_gg_taxonomy_file, fout_gg_taxonomy, 'ssu_gg', ssu_gg_fna_file)

                ssu_silva_taxonomy_file = os.path.join(
                    gpath, self.ssu_silva_taxonomy_file)
                ssu_silva_fna_file = os.path.join(
                    gpath, self.ssu_silva_fna_file)
                ssu_silva_summary_file = os.path.join(
                    gpath, self.ssu_silva_summary_file)
                ssu_count = self._parse_taxonomy_file(gid,
                                                      ssu_silva_taxonomy_file,
                                                      fout_ssu_silva_taxonomy,
                                                      'ssu_silva',
                                                      ssu_silva_fna_file,
                                                      ssu_silva_summary_file)

                lsu_silva_23s_taxonomy_file = os.path.join(
                    gpath, self.lsu_silva_23s_taxonomy_file)
                lsu_silva_23s_fna_file = os.path.join(
                    gpath, self.lsu_silva_23s_fna_file)
                lsu_silva_23s_summary_file = os.path.join(
                    gpath, self.lsu_silva_23s_summary_file)
                lsu_23s_count = self._parse_taxonomy_file(
                    gid, lsu_silva_23s_taxonomy_file, fout_lsu_silva_23s_taxonomy, 'lsu_silva_23s', lsu_silva_23s_fna_file, lsu_silva_23s_summary_file)

                lsu_5S_fna_file = os.path.join(
                    gpath, self.lsu_5S_fna_file)
                lsu_5S_summary_file = os.path.join(
                    gpath, self.lsu_5S_summary_file)
                lsu_5S_count = self._parse_lsu_5S_files(
                    gid, fout_lsu_5S, lsu_5S_fna_file, lsu_5S_summary_file)

                fout_ssu_silva_count.write(
                    '%s\t%d\n' % (gid, ssu_count))
                fout_lsu_silva_23s_count.write(
                    '%s\t%d\n' % (gid, lsu_23s_count))
                fout_lsu_5S_count.write(
                    '%s\t%d\n' % (gid, lsu_5S_count))

                trna_file = os.path.join(
                    gpath, 'trna', gid + '_trna_stats.tsv')
                self._parse_trna_file(
                    gid, trna_file, fout_trna_count)

        fout_nt.close()
        fout_gene.close()
        fout_gg_taxonomy.close()
        fout_ssu_silva_taxonomy.close()
        fout_lsu_silva_23s_taxonomy.close()
        fout_lsu_5S.close()
        fout_ssu_silva_count.close()
        fout_lsu_silva_23s_count.close()
        fout_lsu_5S_count.close()
        fout_trna_count.close()


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
        with mp.Pool(processes=self.cpus) as pool:
            for result in tqdm(pool.imap_unordered(self._producer, input_files),
                               total=len(input_files), ncols=100, unit='genome'):
                # warned here rather than in the worker: several processes
                # appending to one log file interleave, and the parent is
                # reading every result anyway
                for gid, what, missing_file in result:
                    self.logger.warning('{} has no {}: {}'.format(
                        gid, MISSING_LABEL[what], missing_file))
                missing.extend(result)

        self.report_missing(missing, len(input_files), out_dir)

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

    def _producer(self, job: MetadataJob) -> List[MissingFile]:
        """Process each genome.

        A genome missing a file is reported and passed over rather than ending
        the run. Which files are there decides how much of the genome can be
        done: the gene metadata needs the GFF and the genome size both, so it
        needs the two files, while the nucleotide metadata needs only the
        FASTA and is written whenever the FASTA is there. The two are separate
        files that create_metadata_tables() reads independently, so half a
        genome is worth having and is not done again when prodigal catches up.

        Parameters
        ----------
        job : MetadataJob
            The genome's accession, its genomic FASTA and the GFF of the genes
            Prodigal called.

        @return: the files this genome should have had and did not, which is
                 empty for a genome that was processed in full.
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
            return missing

        # clean up old log files
        log_file = os.path.join(full_genome_dir, 'genometk.log')
        if os.path.exists(log_file):
            os.remove(log_file)

        # calculate metadata
        self.nucleotide(genome_file,full_genome_dir)
        if not missing:
            self.gene(genome_file,gff_file,full_genome_dir)

        return missing

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



