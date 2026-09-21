import os
import csv
import gzip
import logging
import multiprocessing as mp
import shutil
from typing import Dict, List, Optional, Tuple

from tqdm import tqdm

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

    def __init__(self, tmp_dir: str = '/tmp/', cpus: int = 1) -> None:
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

        @return: None
        """

        self.tmp_dir = tmp_dir
        self.cpus = cpus

        check_dependencies(['prodigal'])

        self.logger = logging.getLogger('timestamp')

    def run(self,
            gtdb_genome_path_file: str,
            trans_table_file: str,
            tt_override_file: Optional[str] = None,
            all_genomes: bool = False) -> bool:
        """Call genes for every genome of a release that still needs them.

        The translation table of each genome is given rather than chosen: trans_table
        predicted it with gTranslate, and --tt_override corrects the predictions that
        were wrong. Prodigal is handed that table, so it calls the genes once instead
        of calling them under tables 4 and 11 and keeping whichever coded more of the
        genome, which is a rule that cannot express table 25 at all.

        A genome with no table is not called. gTranslate returns no prediction for a
        few genomes of a release -- eight of r237's 1.35M, some of which have no
        genomic FASTA to predict from at all -- and there is nothing to call their
        genes under: Prodigal choosing a table by coding density is the very thing
        the summary is handed over to prevent. They are warned about and left, and
        --tt_override is how one is given a table and called after all.

        Which genomes those are is settled before any genes are called rather than
        discovered one at a time, because a run that meets the gap genome by genome
        meets it hours in, having already called everything ahead of it.

        A summary covering NO genome of the release is a different thing and still
        stops the run. That is the wrong file rather than a few unpredictable
        genomes, and carrying on would call nothing at all and report that the run
        had finished.

        Two passes then, because deciding is cheap and calling is not. The first
        spreads prodigal_parser() over the pool to sort the release into genomes
        whose proteins are already there and vouched for and genomes that are not;
        the second hands what is left to run_prodigal(). A release is mostly carried
        over from the one before, so the first pass is what keeps a run proportional
        to the genomes that are actually new.

        Parameters
        ----------
        gtdb_genome_path_file : str
            genome_dirs file of the release, accession and genome directory per line.
        trans_table_file : str
            gtranslate.translation_table_summary.tsv, as trans_table writes it.
        tt_override_file : str
            Corrections to those predictions, or None.
        all_genomes : bool
            Call genes for every genome, discarding the previous results, rather than
            only for those that need it (default = False).

        @return: True, unconditionally. It reports that the run finished rather than
                 that every genome succeeded; a genome that could not be called is
                 warned about as it is met.
        """

        tables, sources = read_translation_tables(trans_table_file, tt_override_file)
        corrected = sum(1 for source in sources.values() if source == SOURCE_OVERRIDE)
        self.logger.info('Read translation tables for {:,} genomes from {}{}.'.format(
            len(tables), trans_table_file,
            ', {:,} of them corrected by {}'.format(corrected, tt_override_file)
            if tt_override_file else ''))

        # get all path to genome data files for all genomes in GTDB release
        list_genome_tuples = []
        with open(gtdb_genome_path_file,'r') as fh :
            for line in tqdm(fh, ncols=100, leave=False, desc='Reading genomes'):
                tokens = line.strip().split('\t')
                gid = tokens[0]
                gpath = tokens[1]
                list_genome_tuples.append((gid, gpath, all_genomes))

        total = len(list_genome_tuples)
        uncalled = sorted(gid for gid, _, _ in list_genome_tuples if gid not in tables)

        # no genome of the release having a table is the wrong file, not a release
        # with a few genomes nothing could be predicted for; carrying on would call
        # nothing and report that the run had finished
        if uncalled and len(uncalled) == total:
            raise RuntimeError(
                'None of the {:,} genomes of this release has a translation table '
                'in {}. That file is for another release, or trans_table has not '
                'been run over this one.'.format(total, trans_table_file))

        if uncalled:
            list_genome_tuples = [row for row in list_genome_tuples if row[0] in tables]
            self.logger.warning(
                'warning: {:,} of {:,} genome(s) have no translation table in {} and '
                'their genes are NOT being called: {}{}. They are the genomes '
                'trans_table could not predict, which it names in full in '
                'gtranslate_no_prediction.tsv; give one a table with --tt_override '
                'to have it called.'.format(
                    len(uncalled), total, trans_table_file,
                    ', '.join(uncalled[:MISSING_TABLES_LOGGED]),
                    ' ...' if len(uncalled) > MISSING_TABLES_LOGGED else ''))

        # determine which genomes need to be processed by Prodigal
        self.logger.info('Running prodigal on genomes.')
        with mp.Pool(processes=self.cpus) as pool:
            genome_paths = list(tqdm(pool.imap_unordered(self.prodigal_parser, list_genome_tuples),
                                    total=len(list_genome_tuples), unit='genome'))


        # run Prodigal on genomes requiring gene calling
        genome_paths = [x for x in genome_paths if x != ('null', 'null')]

        # what the first pass was for: a release is mostly carried over from the
        # one before, so the two counts are how much of this run is work and how
        # much of it was already done. A run over a whole release that reports
        # every genome as requiring gene calling has not carried anything across
        skipped = len(list_genome_tuples) - len(genome_paths)
        self.logger.info(
            '{:,} genome(s) require gene calling; {:,} already have valid '
            'Prodigal results.'.format(len(genome_paths), skipped))

        self.run_prodigal(genome_paths, tables, sources)

        return True

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
                # verify checksum
                checksum = sha256_rb(gzip.GzipFile(fileobj=open(aa_gene_file, 'rb')))
                cur_checksum = open(checksum_file).readline().strip()
                if checksum == cur_checksum:
                    return ('null', 'null')

            return (gid, gpath)

    def run_prodigal(self,
                     genome_paths: List[Tuple[str, str]],
                     tables: Dict[str, int],
                     sources: Dict[str, str]) -> None:
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

        @return: None; the results are written into the genome directories.
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
        summary_stats = Prodigal(cpus=self.cpus).run(tasks)

        # everything else a genome's directory holds was written by the worker
        # that called its genes; this is the one record the wrapper has no
        # business knowing about
        self.logger.info('Recording the translation table of each genome.')
        for task in tasks:
            table_file = os.path.join(os.path.dirname(task.aa_gene_file),
                                      'prodigal_translation_table.tsv')
            with open(table_file, 'w') as handle:
                handle.write('{}\t{}\t{}\n'.format(
                    'best_translation_table',
                    summary_stats[task.genome_id].best_translation_table,
                    sources.get(task.genome_id, SOURCE_PREDICTED)))
