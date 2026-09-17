import csv
import gzip
from collections import namedtuple
from typing import Dict, Optional


# Columns of the translation table summary gTranslate writes, which trans_table
# produces and prodigal consumes. The file is read by COLUMN NAME: gTranslate has
# changed what it reports between releases, and a column taken by position is a
# column that can quietly become a different one.
TT_SUMMARY_GENOME = 'user_genome'
TT_SUMMARY_TABLE = 'best_tln_table'
TT_SUMMARY_DENSITY_4 = 'coding_density_4'
TT_SUMMARY_DENSITY_11 = 'coding_density_11'

# The coding densities gTranslate reports are PERCENTAGES, where the rule below
# was written against fractions, so the thresholds are scaled rather than the
# values: 0.05 of 1 is 5 of 100, and 0.7 of 1 is 70 of 100.
CHECKM_DENSITY_MARGIN = 5.0
CHECKM_DENSITY_FLOOR = 70.0


def read_translation_table_summary(summary_file: str) -> Dict[str, Dict[str, str]]:
    """Read the translation table summary gTranslate writes.

    Shared by the command that writes it and the command that acts on it, so that
    the two never disagree about which column holds what.

    Parameters
    ----------
    summary_file : str
        gtranslate.translation_table_summary.tsv, of one batch or of a release.

    @return: genome ID to its row, keyed by column name.
    """

    predictions = {}
    with open(summary_file) as handle:
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