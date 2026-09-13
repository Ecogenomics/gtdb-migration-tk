import gzip
from collections import namedtuple


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