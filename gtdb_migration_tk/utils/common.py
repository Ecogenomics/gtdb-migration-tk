import csv
import os
import gzip
import logging
import re
import subprocess
from collections import namedtuple
from typing import Dict, Optional, Tuple

from tqdm import tqdm

from gtdb_migration_tk.biolib_lite.common import canonical_gid


# Columns of the translation table summary gTranslate writes, which trans_table
# produces and prodigal consumes. The file is read by COLUMN NAME: gTranslate has
# changed what it reports between releases, and a column taken by position is a
# column that can quietly become a different one.
TT_SUMMARY_GENOME = 'user_genome'
TT_SUMMARY_TABLE = 'best_tln_table'
TT_SUMMARY_DENSITY_4 = 'coding_density_4'
TT_SUMMARY_DENSITY_11 = 'coding_density_11'

# gTranslate measures these of every genome it predicts for, so the comparison
# carries them rather than reading 1.3M FASTA files to learn what has already
# been read. The output files name them as gTranslate does.
TT_SUMMARY_GC = 'gc_percent'
TT_SUMMARY_N50 = 'n50'
TT_SUMMARY_GENOME_SIZE = 'genome_size'

# The coding densities gTranslate reports are PERCENTAGES, where the rule below
# was written against fractions, so the thresholds are scaled rather than the
# values: 0.05 of 1 is 5 of 100, and 0.7 of 1 is 70 of 100.
CHECKM_DENSITY_MARGIN = 5.0
CHECKM_DENSITY_FLOOR = 70.0


# The first two bytes of a gzip member. The summary is read by what the file IS
# and not by what it is called: gTranslate writes a batch's plain and trans_table
# writes the release's compressed, and whoever hands one to prodigal should not
# have to have kept track of which -- nor be caught out by a release summary
# somebody gunzipped, or one renamed on the way to another machine.
GZIP_MAGIC = b'\x1f\x8b'

# The genome assembly sizes a command processes by default: select_genomes selects
# none outside them, and trnascan, rna_silva and rna_ltp process none larger. The
# largest bacterial and archaeal genomes are under 20 Mbp, so nothing GTDB would
# keep comes near the upper bound; what does is a metagenome deposited as one
# genome, such as GCA_964261755.1 at 9,529 Mbp, which held an r237 rna_silva batch
# for a day on a single blastn. The lower bound is in kbp and the upper in Mbp, as
# --min_genome_size and --max_genome_size take them.
DEFAULT_MIN_GENOME_SIZE = 10.0
DEFAULT_MAX_GENOME_SIZE = 100.0
KBP = 1000
MBP = 1000000

# Where prodigal leaves a genome's called proteins, and what it calls them. The
# directory and the extension are one fact about a genome directory, written by
# prodigal and read by everything downstream of it -- hmmsearch, top_hit, checkm,
# the metadata -- so protein_fasta() is what names the file rather than each of
# them joining the same three pieces.
PRODIGAL_DIR = 'prodigal'
PROTEIN_FASTA_EXT = '_protein.faa.gz'


def protein_fasta(accession: str, genome_dir: str) -> str:
    """The called proteins of a genome, as prodigal filed them.

    Named for the ACCESSION rather than for the assembly, which is what the
    genomic FASTA is named for: prodigal writes its results under the genome's
    accession, so the two files of one genome do not share a stem. That is also
    why the accession comes first, against the reading of it: this is one of the
    functions batching.plan_batches() takes to name the file a batch is planned
    around, and they are all (accession, genome directory).

    Parameters
    ----------
    accession : str
        Accession of the genome, which names the file.
    genome_dir : str
        Genome directory of the release.

    @return: path of the protein FASTA, which may not exist.
    """

    return os.path.join(genome_dir, PRODIGAL_DIR, accession + PROTEIN_FASTA_EXT)


# How each external program the toolkit runs is asked its version, and what the
# answer looks like. Every program a command runs is recorded twice: once in the
# run's log, and once in a <program>.version file beside the results it made --
# the results outlive the run, a genome's being carried across from release to
# release while its sequences are unchanged, so the version that made them is not
# the version of the run that last looked at them and only the file can say.
# Asked of the program rather than of conda or of a path, since what matters is
# what ran and a machine of a shared run can have another build first on PATH. A
# program the toolkit runs goes here, and is asked by the command running it.
HMMER_VERSION = r'HMMER \d\S*(?: \([^)]*\))?'
VERSION_QUERIES: Dict[str, Tuple[Tuple[str, ...], str]] = {
    'tRNAscan-SE': (('tRNAscan-SE', '-h'), r'tRNAscan-SE \d\S*(?: \([^)]*\))?'),
    'prodigal': (('prodigal', '-v'), r'Prodigal V\d\S*(?: .*)?'),
    'hmmsearch': (('hmmsearch', '-h'), HMMER_VERSION),
    'nhmmer': (('nhmmer', '-h'), HMMER_VERSION),
    'blastn': (('blastn', '-version'), r'blastn: \d\S*'),
    'makeblastdb': (('makeblastdb', '-version'), r'makeblastdb: \d\S*'),
    'gtranslate': (('gtranslate', '--version'), r'gtranslate: version \d\S*'),
    'checkm2': (('checkm2', '--version'), r'(?m)^\d+\.\d+\S*'),
    'checkm': (('checkm', '-h'), r'CheckM v\d\S*'),
    'busco': (('busco', '--version'), r'BUSCO \d\S*'),
}
VERSION_EXT = '.version'


def program_version(program: str) -> str:
    """The version of an external program, as the program itself states it.

    tRNAscan-SE and Prodigal say it on stderr and the rest on stdout, so both
    are read. The exit status is not: a version printed is a version, whatever
    the program thought of the rest of what it was asked.

    Parameters
    ----------
    program : str
        A key of VERSION_QUERIES, as the program is called on the command line.

    @return: the version as the program printed it, e.g.
             'Prodigal V2.6.3: February, 2016'.

    Raises
    ------
    RuntimeError
        The program could not be run, or said nothing that looks like a
        version. A version file saying something other than what ran is worse
        than none, so this is not guessed at.
    """

    command, pattern = VERSION_QUERIES[program]
    try:
        proc = subprocess.run(list(command), stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, timeout=300)
    except (OSError, subprocess.SubprocessError) as error:
        raise RuntimeError('Could not run {} to learn its version: {}'.format(
            ' '.join(command), error))

    output = proc.stdout.decode('utf-8', 'replace')
    match = re.search(pattern, output)
    if not match:
        raise RuntimeError(
            '{} did not state a version matching {!r}; it said: {}'.format(
                ' '.join(command), pattern,
                ' '.join(output.split())[:200] or 'nothing'))

    return match.group(0).strip()


def record_program_version(program: str) -> str:
    """Ask an external program its version, and say it in the run's log.

    Called once per run, where the command starts, rather than by each worker:
    it is one binary for the whole run, and a program that will not say what it
    is is met before any work is claimed rather than inside it.

    Parameters
    ----------
    program : str
        A key of VERSION_QUERIES.

    @return: the version, for write_version_file() to put beside each result.
    """

    version = program_version(program)
    logging.getLogger('timestamp').info('Using {}: {}.'.format(program, version))
    return version


def version_file(directory: str, program: str) -> str:
    """Where a program's version is recorded, beside the results it made.

    Parameters
    ----------
    directory : str
        Directory holding the results.
    program : str
        The program, as it is called; the file is named for it in lower case.

    @return: e.g. <directory>/trnascan-se.version.
    """

    return os.path.join(directory, program.lower() + VERSION_EXT)


def write_version_file(directory: str, program: str, version: str) -> None:
    """Record the version of the program that made some results, beside them.

    Parameters
    ----------
    directory : str
        Directory holding the results.
    program : str
        The program, as it is called.
    version : str
        What record_program_version() returned.

    @return: None
    """

    with open(version_file(directory, program), 'w') as handle:
        handle.write('{}\n'.format(version))


def open_text(path: str):
    """Open a text file whether or not it is gzipped.

    Parameters
    ----------
    path : str
        File to open.

    @return: a file object of text lines, to be used as a context manager.
    """

    with open(path, 'rb') as handle:
        compressed = handle.read(len(GZIP_MAGIC)) == GZIP_MAGIC

    return gzip.open(path, 'rt') if compressed else open(path)


def read_translation_table_summary(summary_file: str) -> Dict[str, Dict[str, str]]:
    """Read the translation table summary gTranslate writes.

    Shared by the command that writes it and the command that acts on it, so that
    the two never disagree about which column holds what. Gzipped or not is
    decided by the file's own first two bytes, so a batch's summary and a
    release's read the same way.

    Parameters
    ----------
    summary_file : str
        gtranslate.translation_table_summary.tsv of one batch, or the release's
        gtranslate.translation_table_summary.tsv.gz.

    @return: genome ID to its row, keyed by column name.
    """

    predictions = {}
    with open_text(summary_file) as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            genome = row.get(TT_SUMMARY_GENOME)
            if genome:
                predictions[genome] = row

    return predictions


def checkm_translation_table(density_4: str, density_11: str) -> Optional[int]:
    """The translation table the coding density rule alone would choose.

    This is what CheckM and Prodigal do when nothing tells them the table: call
    the genes under tables 4 and 11 and take 4 only where it codes appreciably
    more of the genome. It is reported beside the gTranslate prediction so that
    the two can be compared -- a classifier trained on GTDB against a threshold
    on two numbers -- without calling any genes again.

    Parameters
    ----------
    density_4 : str
        Coding density under table 4, as a percentage.
    density_11 : str
        Coding density under table 11, as a percentage.

    @return: 4 or 11, or None where either density is missing or unreadable.
    """

    try:
        coding_4, coding_11 = float(density_4), float(density_11)
    except (TypeError, ValueError):
        return None

    if (coding_4 - coding_11 > CHECKM_DENSITY_MARGIN
            and coding_4 > CHECKM_DENSITY_FLOOR):
        return 4

    return 11




def read_gtdb_metadata(metadata_file, fields):
    """Parse genome quality from GTDB metadata.
    Parameters
    ----------
    metadata_file : str
        Metadata for all genomes in CSV file.
    fields : iterable
        Fields  to read.
    Return
    ------
    dict : d[genome_id] -> namedtuple
        Value for fields indicted by genome IDs.
    """

    gtdb_metadata = namedtuple('gtdb_metadata', ' '.join(fields))
    m = {}

    with open(metadata_file) as f:
        headers = f.readline().strip().split('\t')

        genome_index = headers.index('accession')

        indices = []
        for field in fields:
            indices.append(headers.index(field))

        for line in f:
            line_split = line.strip().split('\t')
            genome_id = line_split[genome_index]

            values = []
            for i in indices:
                # save values as floats or strings
                v = line_split[i]
                try:
                    values.append(float(v))
                except ValueError:
                    if v is None or v == '' or v == 'none':
                        values.append(None)
                    elif v == 'f' or v.lower() == 'false':
                        values.append(False)
                    elif v == 't' or v.lower() == 'true':
                        values.append(True)
                    else:
                        values.append(v)
            m[genome_id] = gtdb_metadata._make(values)

    return m


def count_lines(file_path: str) -> int:
    """Count the lines in a file, in order to size a progress bar.

    Gzipped files are counted too: GTDB stores the NCBI assembly summaries
    compressed, and reading one as text would decode gzip bytes as UTF-8 and
    fail rather than merely miscount.

    Newlines are counted in binary blocks rather than by iterating lines, which
    for a file of this size is several times faster and needs no decoding at all
    -- the caller only wants a number to size a bar with.

    Parameters
    ----------
    file_path : str
        File to read, optionally gzipped.

    @return: number of lines in the file.
    """

    opener = gzip.open if file_path.endswith('.gz') else open

    with opener(file_path, 'rb') as check_file:
        return sum(block.count(b'\n')
                   for block in iter(lambda: check_file.read(1024 * 1024), b''))


def read_taxonomy(taxonomy_file: str) -> Dict[str, str]:
    """Read the standardised NCBI taxonomy of the genomes.

    The file is the two-column accession / semicolon-separated lineage TSV
    ncbi_metadata_sync writes. Each genome is recorded under the accession as
    given AND under its canonical form, so that a GenBank genome of a release
    finds the lineage the taxonomy holds against its RefSeq counterpart; an
    accession given exactly is preferred to a canonical match.

    Here rather than in either command that reads it: trans_table reads the
    lineage to decide which genomes a translation table may be doubted for, and
    trnascan reads the domain out of it to choose tRNAscan-SE's model. The
    command modules do not import one another.

    Parameters
    ----------
    taxonomy_file : str
        Standardised NCBI taxonomy file.

    @return: accession, and canonical accession, to lineage.
    """

    exact, canonical = {}, {}
    with open(taxonomy_file) as handle:
        for line in tqdm(handle, ncols=100, leave=False, desc='Reading taxonomy'):
            line = line.rstrip('\n')
            if not line:
                continue
            accession, _, lineage = line.partition('\t')
            if not lineage:
                continue
            exact[accession] = lineage
            canonical.setdefault(canonical_gid(accession), lineage)

    canonical.update(exact)

    return canonical


# The two domains a genome is searched as, by the commands whose models differ
# between them: tRNAscan-SE's bacterial and archaeal models, and the bac_ and ar_
# rRNA HMMs of rna_silva. Both files a domain is read from spell it with the GTDB
# rank prefix -- the domain file in its Predicted domain column, the taxonomy in
# the first rank of each lineage -- so one table reads both.
DOMAIN_ARCHAEA = 'Archaea'
DOMAIN_BACTERIA = 'Bacteria'
DOMAIN_OF_TAXON = {'d__Archaea': DOMAIN_ARCHAEA, 'd__Bacteria': DOMAIN_BACTERIA}

# Columns of the GTDB domain file, read by name. The prediction is made from the
# genome's marker genes, and is the string below for a genome it could not be
# made for -- those fall through to the NCBI taxonomy.
DOMAIN_FILE_GENOME = 'Genome Id'
DOMAIN_FILE_DOMAIN = 'Predicted domain'
NO_PREDICTION = 'None'

# GTDB names a genome for the database it came from, e.g. GB_GCA_000009065.1,
# where a genome_dirs file names it GCA_000009065.1.
GTDB_ID_PREFIXES = ('GB_', 'RS_')


def read_domains(gtdb_domain_file: str, taxonomy_file: str) -> Dict[str, str]:
    """Read the domain of each genome from GTDB's prediction and NCBI's taxonomy.

    The NCBI lineages are read first and GTDB's predictions written over them,
    so a genome GTDB has a prediction for is searched on that and one it has
    none for falls through to where NCBI filed it. GTDB's call is preferred
    because it is made from the genome rather than from where NCBI filed it, and
    it is the one that catches a genome under the wrong domain at NCBI. Both are
    held under the accession as given and under its canonical form, so a GenBank
    genome finds what is recorded against its RefSeq counterpart.

    Here rather than in either command that reads it: trnascan chooses
    tRNAscan-SE's model by it and rna_silva its rRNA HMM, a genome told the
    wrong domain gets a worse answer from both rather than an error, and the
    command modules do not import one another.

    Parameters
    ----------
    gtdb_domain_file : str
        GTDB domain report: Genome Id and Predicted domain, by column name.
    taxonomy_file : str
        Standardised NCBI taxonomy, accession and lineage per line.

    @return: accession, and canonical accession, to DOMAIN_ARCHAEA or
             DOMAIN_BACTERIA.
    """

    # counted by canonical accession rather than by entry: each genome is held
    # under the accession as given and under its canonical form, so counting
    # the table would report every genome twice
    domains: Dict[str, str] = {}
    from_ncbi = set()
    for key, lineage in read_taxonomy(taxonomy_file).items():
        domain = DOMAIN_OF_TAXON.get(lineage.split(';')[0])
        if domain:
            domains[key] = domain
            from_ncbi.add(canonical_gid(key))

    predicted = set()
    with open(gtdb_domain_file) as handle:
        header = handle.readline().rstrip('\n').split('\t')
        genome_idx = header.index(DOMAIN_FILE_GENOME)
        domain_idx = header.index(DOMAIN_FILE_DOMAIN)

        for line in handle:
            fields = line.rstrip('\n').split('\t')
            if len(fields) <= max(genome_idx, domain_idx):
                continue

            gid = fields[genome_idx]
            for prefix in GTDB_ID_PREFIXES:
                if gid.startswith(prefix):
                    gid = gid[len(prefix):]
                    break

            # 'None' is what the column holds for a genome the prediction could
            # not be made for; NCBI's taxonomy answers for those
            domain = DOMAIN_OF_TAXON.get(fields[domain_idx])
            if not domain:
                continue

            domains[gid] = domain
            domains[canonical_gid(gid)] = domain
            predicted.add(canonical_gid(gid))

    logging.getLogger('timestamp').info(
        'Read the domain of {:,} genome(s): {:,} predicted by GTDB and {:,} '
        'taken from the NCBI taxonomy.'.format(
            len(predicted | from_ncbi), len(predicted),
            len(from_ncbi - predicted)))

    return domains


def domain_of(domains: Dict[str, str], accession: str) -> Optional[str]:
    """The domain read_domains() found for a genome.

    Parameters
    ----------
    domains : dict
        What read_domains() returned.
    accession : str
        Genome accession, as the genome_dirs file names it.

    @return: DOMAIN_ARCHAEA or DOMAIN_BACTERIA, or None where neither file
             answers for the genome -- which each caller decides what to do
             with, and says how many there were.
    """

    domain = domains.get(accession)
    if domain is None:
        domain = domains.get(canonical_gid(accession))

    return domain
