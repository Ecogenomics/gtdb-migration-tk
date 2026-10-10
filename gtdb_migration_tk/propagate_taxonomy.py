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

"""Carry the GTDB taxonomy of one release to the genomes of the next.

propagate_gtdb_taxonomy reads the previous release's taxonomy and
representatives from its GTDB metadata files (--gtdb_metadata_prev, the
archaeal and bacterial files of a release alike) and the genomes of the new
release from the database. It read them from a metadata file exported from the
database (--gtdb_metadata_cur) until 0.1.63, which was an export taken by hand
for the one run, and could be of another database than the one
update_propagated_tax then wrote to.

A genome inherits from the genome of the previous release with its accession or
else its canonical one (biolib_lite.common.canonical_gid()): another version,
the other database, or both. Accessions were rewritten by hand -- versions up to
four past the current, then RS_GCF and GB_GCA swapped -- which matched a genome
to every predecessor it found rather than one, and counted the genomes that
moved to RefSeq as moved to GenBank and the other way round. A canonical
accession names one genome on each side over r226's metadata and r237's
database, and a run where one names two is refused, so a genome has one
predecessor or none.

The database's taxonomy of a genome of the previous release must be the
previous release's, or not yet set: no rank below the domain holding more than
its prefix, whether NULL or the bare 'p__' ... 's__' update_propagated_tax
--truncate_taxonomy writes. The run stopped at
the first that was not, and failed with a KeyError on one the database held no
taxonomy for; now every one is listed beside the log and nothing is written.

The two files are written in --output_dir under fixed names, gzipped, for
update_propagated_tax to read from the same directory (-i), which writes them
and sets the ranks below the domain to NULL for every genome they do not name.
"""

import contextlib
import csv
import logging
import os
import re
import sys
import time

from gtdb_migration_tk.biolib_lite.common import canonical_gid, check_file_exists
from gtdb_migration_tk.biolib_lite.logger import log_directory
from gtdb_migration_tk.biolib_lite.taxonomy import Taxonomy
from gtdb_migration_tk.database_configuration import GenomeDatabaseConnectionFTPUpdate
from gtdb_migration_tk.database_configuration.GenomeDatabaseConnectionFTPUpdate import one_transaction, reset_unwritten
from gtdb_migration_tk.gtdb_lite.gtdb_importer import GTDBImporter, UNKNOWN_EXAMPLES, UnknownGenomesError, id_at_source
from gtdb_migration_tk.utils.common import GZIP_SUFFIX, open_gzip_text, open_text, read_gtdb_metadata
from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm

csv.field_size_limit(sys.maxsize)

# what propagate_gtdb_taxonomy writes in --output_dir and update_propagated_tax
# reads from -i: each genome's taxonomy (genome_id, taxonomy string), and
# whether each genome of the database is a representative (genome_id, True or
# False); neither has a header
TAXONOMY_NAME = 'gtdb_taxonomy_propagated.tsv' + GZIP_SUFFIX
REPS_NAME = 'gtdb_sp_reps_propagated.tsv' + GZIP_SUFFIX

# beside the log: each genome whose taxonomy in the database is not the previous release's
MISMATCHES_NAME = 'gtdb_taxonomy_mismatches.tsv'

# the columns read from the previous release's metadata files
PREVIOUS_COLUMNS = ('accession', 'gtdb_taxonomy', 'gtdb_representative')

# each genome of the database, as GTDB names it, and its seven GTDB ranks; a
# genome with no metadata_taxonomy row has every rank NULL
CURRENT_GENOMES = ("SELECT gs.external_id_prefix || '_' || g.id_at_source, "
                   + ', '.join('m.gtdb_' + rank for rank in Taxonomy.rank_labels)
                   + ' FROM genomes g JOIN genome_sources gs ON gs.id = g.genome_source_id '
                     'LEFT JOIN metadata_taxonomy m ON m.id = g.id')

REFSEQ_PREFIX = 'RS_'
USER_PREFIX = 'U_'

# how many genomes a line of the log or an error names
EXAMPLES = 10

# set_gtdb_domain: the marker sets a genome's domain is read from (marker_sets.id),
# and the genomes it is given one -- its id, as GTDB names it, as the database
# does, and its NCBI taxonomy, which may be NULL. Every genome without a domain
# is given one: a genome with no NCBI taxonomy was passed over until 0.1.77, and
# r237 has 1,483, their NCBI taxids deleted from the taxdump (delnodes.dmp) though
# the assembly summaries still name them. Where the markers cannot decide and NCBI
# gives no domain either, the genome is given DEFAULT_DOMAIN.
BACTERIAL_MARKER_SET = 1
ARCHAEAL_MARKER_SET = 19
DEFAULT_DOMAIN = 'd__Bacteria'
MISSING_DOMAIN = (
    "SELECT mt.id, gs.external_id_prefix || '_' || g.id_at_source, g.id_at_source, mt.ncbi_taxonomy "
    'FROM metadata_taxonomy mt JOIN genomes g ON g.id = mt.id '
    'JOIN genome_sources gs ON gs.id = g.genome_source_id '
    "WHERE mt.gtdb_domain IS NULL OR mt.gtdb_domain IN ('none', 'd__')")

# How many markers of each set each genome of a range of aligned_markers' blocks
# has a gene aligned for (a row with an e-value; a row of gaps is a marker not
# found). aligned_markers is read whole, in ALIGNED_MARKER_PARTS ranges of its
# blocks, each read once, so that a bar can say how far the read is: PostgreSQL
# says nothing of a query's progress. The reading is the cost. Over gtdb_r237_dev's
# 292,783 genomes without a domain, two queries a genome, as it was, took 21 ms a
# genome, some 1.7 hours; one query, read by a parallel scan on three processes,
# 16 minutes, with nothing to say how far it was; the ranges, on one connection,
# 19 minutes, a tenth every two, and the same counts.
# The ranges are read in one REPEATABLE READ transaction, so all of them see the
# table as it stood when the first was read, as a single query would; its size is
# measured after that, and a row the snapshot sees is in a block below it. The
# last range is bounded too: PostgreSQL reads 'ctid >= x' alone by a scan of the
# whole table.
MARKERS_FOUND = (
    'SELECT am.genome_id, count(*) FILTER (WHERE msc.set_id = %s), count(*) FILTER (WHERE msc.set_id = %s) '
    'FROM aligned_markers am '
    'JOIN marker_set_contents msc ON msc.marker_id = am.marker_id AND msc.set_id IN (%s, %s) '
    "WHERE am.ctid >= %s::tid AND am.ctid < %s::tid AND am.evalue <> '' "
    'GROUP BY am.genome_id')
ALIGNED_MARKERS_BLOCKS = ("SELECT pg_relation_size('aligned_markers') / current_setting('block_size')::int, "
                          "current_setting('block_size')::int")
ALIGNED_MARKER_PARTS = 100

# what set_gtdb_domain writes in --output_dir: the genomes whose markers give
# another domain than NCBI's, those whose markers were too few to say and were
# given NCBI's, and those whose NCBI domain has no d__, which refuse the run
DOMAIN_DISAGREEMENTS_NAME = 'gtdb_domain_disagreements.tsv'
DOMAIN_DISAGREEMENTS_HEADER = ('genome_id', 'ncbi_domain', 'gtdb_domain', 'bac120_percent', 'ar53_percent')
DOMAIN_FROM_NCBI_NAME = 'gtdb_domain_from_ncbi.tsv'
DOMAIN_FROM_NCBI_HEADER = ('genome_id', 'ncbi_domain', 'bac120_percent', 'ar53_percent')
DOMAIN_DEFAULT_NAME = 'gtdb_domain_default.tsv'
DOMAIN_DEFAULT_HEADER = ('genome_id', 'gtdb_domain', 'bac120_percent', 'ar53_percent')
NCBI_DOMAIN_ERRORS_NAME = 'ncbi_domain_errors.tsv'
NCBI_DOMAIN_ERRORS_HEADER = ('genome_id', 'ncbi_domain')


class PropagationError(ValueError):
    """The taxonomy cannot be propagated as the files and the database stand."""


def shared_canonical(accessions):
    """The canonical accessions held by more than one of the accessions given.

    @return: ['G000000001: GB_GCA_000000001.1, RS_GCF_000000001.1', ...], sorted.
    """

    held = {}
    for accession in accessions:
        held.setdefault(canonical_gid(accession), []).append(accession)
    return ['{}: {}'.format(canonical, ', '.join(sorted(names)))
            for canonical, names in sorted(held.items()) if len(names) > 1]


def normalised(taxonomy):
    """A taxonomy string without a trailing ';' or space around a rank.

    @return: str.
    """

    return ';'.join(taxon.strip() for taxon in taxonomy.strip().rstrip(';').split(';'))


def write_tsv(path, header, rows):
    """Write a table with a header, the header alone where there are no rows.

    @return: path.
    """

    with open(path, 'w') as handle:
        handle.write('\t'.join(header) + '\n')
        for row in rows:
            handle.write('\t'.join(row) + '\n')
    return path


def write_gzipped(path, rows):
    """Write a headerless table gzipped, under a temporary name until it is whole.

    Parameters
    ----------
    path : str
        The table.
    rows : iterable of tuple of str
        Its rows.

    @return: the number of rows written.
    """

    partial = path + '.partial'
    written = 0
    try:
        with open_gzip_text(partial) as handle:
            for row in rows:
                handle.write('\t'.join(row) + '\n')
                written += 1
    except BaseException:
        if os.path.exists(partial):
            os.remove(partial)
        raise
    os.replace(partial, path)
    return written


class Propagate(object):
    """Propagate GTDB taxonomy between NCBI releases."""

    def __init__(self, database=None):
        self.logger = logging.getLogger('timestamp')

        self.DEFAULT_DOMAIN_THRESHOLD = 10.0


        if database is not None:
            self.temp_con = GenomeDatabaseConnectionFTPUpdate.GenomeDatabaseConnectionFTPUpdate(database)
            self.temp_con.MakePostgresConnection()
            self.temp_cur = self.temp_con.cursor()

    def previous_release(self, metadata_files):
        """The GTDB taxonomy and representatives of the previous release.

        Parameters
        ----------
        metadata_files : sequence of str
            The previous release's GTDB metadata, gzipped or not, e.g. its
            ar53_metadata and bac120_metadata files; columns are found by name.

        @return: (accession -> taxonomy string, set of representatives), each
                 accession as GTDB writes it (RS_GCF_..., GB_GCA_...).

        Raises
        ------
        PropagationError
            A genome is named twice, by one file or two, or two genomes share a
            canonical accession, so a genome of the new release could not be
            given one predecessor.
        """

        taxonomy = {}
        reps = set()
        repeated = set()
        for metadata_file in metadata_files:
            read = 0
            with open_text(metadata_file) as handle:
                rows = csv.reader(handle, delimiter='\t')
                header = next(rows)
                missing = [c for c in PREVIOUS_COLUMNS if c not in header]
                if missing:
                    raise PropagationError('{} has no {} column.'.format(metadata_file, ', '.join(missing)))
                accession_idx, taxonomy_idx, rep_idx = (header.index(c) for c in PREVIOUS_COLUMNS)
                for row in rows:
                    accession = row[accession_idx]
                    if accession in taxonomy:
                        repeated.add(accession)
                    taxonomy[accession] = row[taxonomy_idx].strip()
                    if row[rep_idx] == 't':
                        reps.add(accession)
                    read += 1
            self.logger.info('Read the GTDB taxonomy of {:,} genomes from {}.'.format(read, metadata_file))

        if repeated:
            raise PropagationError('{:,} genome(s) are named more than once by --gtdb_metadata_prev, e.g.: {}.'.format(
                len(repeated), ', '.join(sorted(repeated)[:EXAMPLES])))
        shared = shared_canonical(taxonomy)
        if shared:
            raise PropagationError('{:,} canonical accession(s) are held by more than one genome of '
                                   '--gtdb_metadata_prev, e.g.: {}.'.format(len(shared), '; '.join(shared[:EXAMPLES])))

        return taxonomy, reps

    def current_release(self):
        """The genomes of the database and the GTDB taxonomy it holds for each.

        @return: accession (RS_GCF_..., GB_GCA_..., U_...) -> the taxonomy string,
                 or None where none is set: no rank below the domain holds more
                 than its prefix, NULL as a genome without a taxonomy has, or
                 'p__' as truncate_taxonomy() writes.

        Raises
        ------
        PropagationError
            Two genomes share a canonical accession.
        """

        self.temp_cur.execute(CURRENT_GENOMES)
        current = {}
        for row in self.temp_cur.fetchall():
            accession, ranks = row[0], row[1:]
            if all(not rank or rank == prefix for rank, prefix in zip(ranks[1:], Taxonomy.rank_prefixes[1:])):
                current[accession] = None
            else:
                current[accession] = ';'.join(rank or prefix for rank, prefix in zip(ranks, Taxonomy.rank_prefixes))
        self.temp_con.rollback()

        shared = shared_canonical(current)
        if shared:
            raise PropagationError('{:,} canonical accession(s) are held by more than one genome of the database, '
                                   'e.g.: {}.'.format(len(shared), '; '.join(shared[:EXAMPLES])))
        self.logger.info('Read {:,} genomes from the database, {:,} with a GTDB taxonomy set.'.format(
            len(current), sum(1 for t in current.values() if t is not None)))
        return current

    def propagate_taxonomy(self, gtdb_metadata_prev, output_dir):
        """Carry the GTDB taxonomy and representatives of the previous release to the database's genomes.

        A genome of the database inherits from the genome of the previous release
        with its accession or, failing that, its canonical accession (another
        version, the other database, or both). The taxonomy the database already
        holds for a genome of the previous release must be the previous release's,
        or not yet set; any other refuses the run before anything is written,
        every one listed in MISMATCHES_NAME beside the log.

        Parameters
        ----------
        gtdb_metadata_prev : sequence of str
            The previous release's GTDB metadata files.
        output_dir : str
            Directory TAXONOMY_NAME and REPS_NAME are written to.

        @return: (path of the taxonomy written, path of the representatives written).

        Raises
        ------
        PropagationError
            As previous_release() and current_release(), or the database holds a
            taxonomy that is not the previous release's.
        """

        prev_taxonomy, prev_reps = self.previous_release(gtdb_metadata_prev)
        self.logger.info('The previous release has {:,} genomes, {:,} of them representatives.'.format(
            len(prev_taxonomy), len(prev_reps)))
        current = self.current_release()

        by_canonical = {canonical_gid(accession): accession for accession in prev_taxonomy}
        predecessor = {}
        kept, new_version, to_refseq, to_genbank = [], [], [], []
        for accession in current:
            if accession in prev_taxonomy:
                predecessor[accession] = accession
                kept.append(accession)
                continue
            if accession.startswith(USER_PREFIX):
                continue
            previous = by_canonical.get(canonical_gid(accession))
            if previous is None:
                continue
            predecessor[accession] = previous
            if previous[:3] == accession[:3]:
                new_version.append(previous)
            elif accession.startswith(REFSEQ_PREFIX):
                to_refseq.append(previous)
            else:
                to_genbank.append(previous)

        mismatches = [(accession, prev_taxonomy[accession], current[accession]) for accession in kept
                      if current[accession] is not None
                      and normalised(current[accession]) != normalised(prev_taxonomy[accession])]
        if mismatches:
            report = os.path.join(log_directory(), MISMATCHES_NAME)
            with open(report, 'w') as handle:
                handle.write('genome_id\tprevious_release\tdatabase\n')
                for row in sorted(mismatches):
                    handle.write('\t'.join(row) + '\n')
            raise PropagationError(
                'The database holds a GTDB taxonomy for {:,} genome(s) of the previous release that is not the '
                'previous release\'s, e.g.: {}. Each is listed in {}; nothing was written.'.format(
                    len(mismatches), ', '.join(sorted(a for a, _, _ in mismatches)[:EXAMPLES]), report))

        self.logger.info('{:,} genomes of the database are genomes of the previous release; '
                         '{:,} representatives among them.'.format(len(kept), len(prev_reps.intersection(kept))))
        for count, what in ((len(new_version), 'have a new version'),
                            (len(to_refseq), 'moved from GenBank to RefSeq'),
                            (len(to_genbank), 'moved from RefSeq to GenBank')):
            self.logger.info('{:,} genomes of the previous release {}.'.format(count, what))

        absent = sorted(set(prev_taxonomy) - set(predecessor.values()))
        if absent:
            self.logger.info('{:,} genomes of the previous release are not in the database, {:,} of them '
                             'representatives, e.g.: {}.'.format(len(absent), len(prev_reps.intersection(absent)),
                                                                 ', '.join(absent[:EXAMPLES])))

        taxonomy_file = os.path.join(output_dir, TAXONOMY_NAME)
        reps_file = os.path.join(output_dir, REPS_NAME)
        written = write_gzipped(taxonomy_file, ((accession, prev_taxonomy[previous])
                                                for accession, previous in sorted(predecessor.items())
                                                if prev_taxonomy[previous]))
        reps = sum(1 for accession in current if predecessor.get(accession) in prev_reps)
        write_gzipped(reps_file, ((accession, str(predecessor.get(accession) in prev_reps))
                                  for accession in sorted(current)))
        self.logger.info('Wrote the GTDB taxonomy of {:,} genomes to {}.'.format(written, taxonomy_file))
        self.logger.info('Wrote whether each of {:,} genomes is a representative, {:,} of them, to {}.'.format(
            len(current), reps, reps_file))

        return taxonomy_file, reps_file

    def truncate_taxonomy(self, metadata_file):
        """Truncate taxonomy string to just domain classification."""

        # get current GTDB taxonomy for all genomes
        gtdb_taxonomy = {}
        with open(metadata_file) as f:
            header = f.readline().strip().split('\t')

            gtdb_taxonomy_index = header.index('gtdb_taxonomy')

            for line in f:
                line_split = line.strip().split('\t')

                gid = line_split[0]
                gtdb_taxa = [t.strip() for t in line_split[gtdb_taxonomy_index].split(';')]
                gtdb_taxonomy[gid] = gtdb_taxa

        gtdbimporter = GTDBImporter(self.temp_cur)
        for i, rank in enumerate(Taxonomy.rank_labels):
            data_to_commit = []
            for gid, taxa in gtdb_taxonomy.items():
                if rank == 'domain':
                    rank_str = taxa[i]
                    data_to_commit.append((gid, rank_str))
                else:
                    data_to_commit.append((gid, Taxonomy.rank_prefixes[i]))

            gtdbimporter.import_metadata_to_db('metadata_taxonomy', 'gtdb_' + rank, 'TEXT', data_to_commit)

    @one_transaction
    def add_propagated_taxonomy(self, input_dir):
        """Write the taxonomy and representatives propagate_gtdb_taxonomy wrote to the database.

        Every rank of every genome TAXONOMY_NAME names is written, and every
        genome's representative status from REPS_NAME, which names every genome
        of the database. A genome the taxonomy does not name keeps its domain,
        which set_gtdb_domain gives a genome GTDB has not classified, and loses
        any rank below it (reset_unwritten()), so no genome keeps a taxonomy the
        previous release did not give it. That was --truncate_taxonomy, which set
        every genome of a metadata file exported by hand (-m) to its domain
        first, until 0.1.63; --genome_list, which kept the writes to the genomes
        of the same export, went with it, the files naming only genomes of the
        database.

        Parameters
        ----------
        input_dir : str
            The --output_dir of propagate_gtdb_taxonomy, holding TAXONOMY_NAME
            and REPS_NAME; each was named on the command line (-t, --rep_file)
            until 0.1.63.

        @return: None
        """

        taxonomy_file = os.path.join(input_dir, TAXONOMY_NAME)
        rep_id_file = os.path.join(input_dir, REPS_NAME)
        for path in (taxonomy_file, rep_id_file):
            check_file_exists(path)

        taxonomy = {}
        with open_text(taxonomy_file) as handle:
            for line in handle:
                genome_id, taxonomy_str = line.rstrip('\n').split('\t')[:2]
                taxonomy[genome_id] = [taxon.strip() for taxon in taxonomy_str.rstrip(';').split(';')]
        self.logger.info('Read the GTDB taxonomy of {:,} genomes from {}.'.format(len(taxonomy), taxonomy_file))

        rep_to_commit = []
        with open_text(rep_id_file) as repfile:
            for line in repfile:
                genome_id, isrep = line.strip().split('\t')
                rep_to_commit.append((genome_id, isrep))

        # Each rank in turn, the ranks below the domain then set to NULL for every
        # genome the taxonomy does not name, and the representatives last. Each
        # is one upsert() or UPDATE, over r237 one to three minutes long, so the
        # log has a line as each ends, with how long it took.
        gtdbimporter = GTDBImporter(self.temp_cur)
        written = [id_at_source(genome_id) for genome_id in taxonomy]
        for i, rank in enumerate(Taxonomy.rank_labels):
            field = 'gtdb_' + rank
            started = time.time()
            gtdbimporter.import_metadata_to_db('metadata_taxonomy', field, 'TEXT',
                                               [(genome_id, taxa[i]) for genome_id, taxa in taxonomy.items()])
            self.logger.info('Wrote {} for {:,} genomes in {:.0f} s.'.format(field, len(taxonomy), time.time() - started))
            if i > 0:
                started = time.time()
                cleared = reset_unwritten(self.temp_cur, 'metadata_taxonomy', field, written)
                message = 'Set {} to NULL for {:,} genome(s) the propagated taxonomy does not name, in ' \
                          '{:.0f} s.'.format(field, cleared, time.time() - started)
                if cleared:
                    self.logger.warning(message)
                else:
                    self.logger.info(message)

        started = time.time()
        gtdbimporter.import_metadata_to_db('metadata_taxonomy', 'gtdb_representative', 'BOOLEAN', rep_to_commit)
        self.logger.info('Wrote whether each of {:,} genomes is a representative, {:,} of them, in {:.0f} s.'.format(
            len(rep_to_commit), sum(1 for _, isrep in rep_to_commit if isrep == 'True'), time.time() - started))

    @one_transaction
    def set_gtdb_domain(self, output_dir):
        """Give each genome with no GTDB domain one, from its marker genes or else NCBI's.

        The domain is the marker set, bac120 or ar53, more of whose markers the
        genome has a gene aligned for, as a share of the set; where neither share
        reaches DEFAULT_DOMAIN_THRESHOLD it is the genome's NCBI domain, and where
        the genome has no NCBI domain either, DEFAULT_DOMAIN. Nothing is printed
        but a bar of how much of aligned_markers has been read: each genome whose
        markers give another domain than NCBI's is listed in
        DOMAIN_DISAGREEMENTS_NAME, each given NCBI's domain in
        DOMAIN_FROM_NCBI_NAME and each given DEFAULT_DOMAIN in
        DOMAIN_DEFAULT_NAME, all in output_dir, with a line of the log counting
        each.

        Parameters
        ----------
        output_dir : str
            Directory the lists are written to.

        @return: None
        """

        # before anything is read, so that every range of aligned_markers is read
        # from the one snapshot; it must be the transaction's first statement
        self.temp_cur.execute('SET TRANSACTION ISOLATION LEVEL REPEATABLE READ')

        self.temp_cur.execute('SELECT set_id, count(*) FROM marker_set_contents WHERE set_id IN (%s, %s) '
                              'GROUP BY set_id', (BACTERIAL_MARKER_SET, ARCHAEAL_MARKER_SET))
        set_size = dict(self.temp_cur.fetchall())

        self.temp_cur.execute(MISSING_DOMAIN)
        genomes = self.temp_cur.fetchall()

        # None where the genome has no NCBI taxonomy, or one naming no domain (d__)
        ncbi_domains = {db_id: (ncbi_taxonomy or '').split(';')[0].strip() or None
                        for db_id, _, _, ncbi_taxonomy in genomes}
        ncbi_domains = {db_id: None if domain == 'd__' else domain for db_id, domain in ncbi_domains.items()}
        self.logger.info('{:,} genome(s) have no GTDB domain, {:,} of them no NCBI domain.'.format(
            len(genomes), sum(1 for domain in ncbi_domains.values() if domain is None)))

        # refused before aligned_markers is read, which takes minutes
        bad_prefix = sorted((genome_id, ncbi_domains[db_id]) for db_id, genome_id, _, _ in genomes
                            if ncbi_domains[db_id] is not None and not ncbi_domains[db_id].startswith('d__'))
        if bad_prefix:
            path = write_tsv(os.path.join(output_dir, NCBI_DOMAIN_ERRORS_NAME), NCBI_DOMAIN_ERRORS_HEADER, bad_prefix)
            self.logger.error('{:,} genome(s) have an NCBI domain without its d__ prefix, e.g. {}; nothing was '
                              'written. Each is listed in {}.'.format(
                                  len(bad_prefix), ', '.join('{} ({})'.format(*g) for g in bad_prefix[:EXAMPLES]),
                                  path))
            sys.exit(1)

        found = self.markers_found({db_id for db_id, _, _, _ in genomes}) if genomes else {}

        domains, from_ncbi, by_default, disagreements = [], [], [], []
        for db_id, genome_id, source_id, _ in genomes:
            ncbi_domain = ncbi_domains[db_id]
            bac_count, arc_count = found.get(db_id, (0, 0))
            bac_percent = bac_count * 100.0 / set_size[BACTERIAL_MARKER_SET]
            arc_percent = arc_count * 100.0 / set_size[ARCHAEAL_MARKER_SET]
            row = ['{:.2f}'.format(bac_percent), '{:.2f}'.format(arc_percent)]
            if bac_percent < self.DEFAULT_DOMAIN_THRESHOLD and arc_percent < self.DEFAULT_DOMAIN_THRESHOLD:
                if ncbi_domain is None:
                    gtdb_domain = DEFAULT_DOMAIN
                    by_default.append([genome_id, gtdb_domain] + row)
                else:
                    gtdb_domain = ncbi_domain
                    from_ncbi.append([genome_id, ncbi_domain] + row)
            else:
                gtdb_domain = 'd__Bacteria' if bac_percent >= arc_percent else 'd__Archaea'
                # a genome NCBI gives no domain has none for the markers to disagree with
                if ncbi_domain is not None and gtdb_domain != ncbi_domain:
                    disagreements.append([genome_id, ncbi_domain, gtdb_domain] + row)
            domains.append((source_id, gtdb_domain))

        started = time.time()
        if domains:
            GTDBImporter(self.temp_cur).import_metadata_to_db('metadata_taxonomy', 'gtdb_domain', 'TEXT', domains)
        self.logger.info('Set gtdb_domain for {:,} genome(s), {:,} d__Bacteria and {:,} d__Archaea, in {:.0f} s.'.format(
            len(domains), sum(1 for _, d in domains if d == 'd__Bacteria'),
            sum(1 for _, d in domains if d == 'd__Archaea'), time.time() - started))

        path = write_tsv(os.path.join(output_dir, DOMAIN_FROM_NCBI_NAME), DOMAIN_FROM_NCBI_HEADER, sorted(from_ncbi))
        message = '{:,} genome(s) have fewer than {:.0f}% of the bac120 and of the ar53 markers aligned and were ' \
                  'given their NCBI domain; each is listed in {}.'.format(
                      len(from_ncbi), self.DEFAULT_DOMAIN_THRESHOLD, path)
        (self.logger.warning if from_ncbi else self.logger.info)(message)

        path = write_tsv(os.path.join(output_dir, DOMAIN_DEFAULT_NAME), DOMAIN_DEFAULT_HEADER, sorted(by_default))
        message = '{:,} genome(s) have fewer than {:.0f}% of the bac120 and of the ar53 markers aligned and no ' \
                  'NCBI domain, and were given {}; each is listed in {}.'.format(
                      len(by_default), self.DEFAULT_DOMAIN_THRESHOLD, DEFAULT_DOMAIN, path)
        (self.logger.warning if by_default else self.logger.info)(message)

        path = write_tsv(os.path.join(output_dir, DOMAIN_DISAGREEMENTS_NAME), DOMAIN_DISAGREEMENTS_HEADER,
                         sorted(disagreements))
        message = '{:,} genome(s) were given a GTDB domain by their markers that is not their NCBI domain; ' \
                  'each is listed in {}.'.format(len(disagreements), path)
        (self.logger.warning if disagreements else self.logger.info)(message)

    def markers_found(self, wanted):
        """How many markers of bac120 and of ar53 each genome has a gene aligned for.

        aligned_markers is read in ALIGNED_MARKER_PARTS ranges of its blocks, a
        bar on the console saying how much has been read and a line of the log
        each tenth. The transaction's snapshot is taken by then.

        Parameters
        ----------
        wanted : set of int
            The genomes (genomes.id) to count; the rest of each range is passed over.

        @return: genomes.id -> [bac120 markers found, ar53 markers found], for each
                 wanted genome with any.
        """

        self.temp_cur.execute(ALIGNED_MARKERS_BLOCKS)
        blocks, block_size = self.temp_cur.fetchone()
        parts = max(1, min(ALIGNED_MARKER_PARTS, blocks))
        bounds = [blocks * i // parts for i in range(parts + 1)]
        self.logger.info('Reading the markers aligned for each genome from aligned_markers ({:.1f} GB) in {:,} '
                         'parts.'.format(blocks * block_size / 1e9, parts))

        found = {}
        started = time.time()
        silent = getattr(self.logger, 'is_silent', False)
        # the bar only on a terminal, and not with --silent; while it is drawn the log's
        # lines go to the console through tqdm, above it rather than through it.
        # logging_redirect_tqdm() gives a logger a console handler whether it had
        # one or not, so it is not used where there is no bar
        with tqdm(total=blocks * block_size, unit='B', unit_scale=True, ncols=100, leave=False,
                  desc='Reading aligned_markers', disable=True if silent else None) as progress, \
                (contextlib.nullcontext() if progress.disable else logging_redirect_tqdm([self.logger])):
            for part in range(parts):
                low, high = bounds[part], bounds[part + 1]
                self.temp_cur.execute(MARKERS_FOUND, (BACTERIAL_MARKER_SET, ARCHAEAL_MARKER_SET,
                                                      BACTERIAL_MARKER_SET, ARCHAEAL_MARKER_SET,
                                                      '({},0)'.format(low), '({},0)'.format(high)))
                for genome, bac_count, arc_count in self.temp_cur.fetchall():
                    if genome in wanted:
                        counts = found.setdefault(genome, [0, 0])
                        counts[0] += bac_count
                        counts[1] += arc_count
                progress.update((high - low) * block_size)
                if (part + 1) * 10 // parts > part * 10 // parts and part + 1 < parts:
                    self.logger.info('Read {}% of aligned_markers in {:.1f} min.'.format(
                        (part + 1) * 100 // parts, (time.time() - started) / 60))

        self.logger.info('Read aligned_markers in {:.1f} min: {:,} of the {:,} genome(s) have a marker '
                         'aligned.'.format((time.time() - started) / 60, len(found), len(wanted)))
        return found

    def propagate_taxonomy_from_reps_to_cluster(self,taxonomy_file,metadata_file,output_file):
        """Propagate labels to all genomes in a cluster. Based on genometreetk"""

        check_file_exists(taxonomy_file)
        check_file_exists(metadata_file)


        # get representative genome information
        rep_metadata = read_gtdb_metadata(metadata_file, ['gtdb_representative',
                                                                  'gtdb_clustered_genomes'])

        rep_metadata = {canonical_gid(gid): values
                        for gid, values in rep_metadata.items()}

        explict_tax = Taxonomy().read(taxonomy_file)

        # sanity check all representatives have a taxonomy string
        rep_count = 0
        for gid in rep_metadata:
            is_rep_genome, clustered_genomes = rep_metadata.get(gid, (None, None))
            if is_rep_genome:
                rep_count += 1
                if gid not in explict_tax:
                    self.logger.error(
                        'Expected to find {} in input taxonomy as it is a GTDB representative.'.format(gid))
                    sys.exit(-1)

        self.logger.info(
            'Identified {:,} representatives in metadata file and {:,} genomes in input taxonomy file.'.format(
                rep_count,
                len(explict_tax)))

        # propagate taxonomy to genomes clustered with each representative
        fout = open(output_file, 'w')
        for rid, taxon_list in explict_tax.items():
            taxonomy_str = ';'.join(taxon_list)
            rid = canonical_gid(rid)

            is_rep_genome, clustered_genomes = rep_metadata[rid]
            if is_rep_genome:
                # assign taxonomy to representative and all genomes in the cluster
                fout.write('{}\t{}\n'.format(rid, taxonomy_str))
                for cid in [gid.strip() for gid in clustered_genomes.split(';')]:
                    cid = canonical_gid(cid)
                    if cid != rid:
                        if cid in rep_metadata:
                            fout.write('{}\t{}\n'.format(cid, taxonomy_str))
                        else:
                            self.logger.warning('Skipping {} as it is not in GTDB metadata file.'.format(cid))
            else:
                self.logger.error(
                    'Did not expected to find {} in input taxonomy as it is not a GTDB representative.'.format(rid))
                sys.exit(-1)

        self.logger.info('Taxonomy written to: {}'.format(output_file))


    @one_transaction
    def add_taxonomy_to_database(self,taxonomy_file,metadata_file,truncate_taxonomy):
        """
        Update the taxonomy in the database, this is usually used after propagate_taxonomy_from_reps_to_cluster function

        @param taxonomy_file: taxonomy file listing all genomes to update; genome id can either be canonical (G123456789),
        normal(GCF_123456789.1) or extended (RS_GCF_123456789.1)
        @param canonical_gid_mapping: TSV file in the format (canonical_id,full_genome_id)
        @return:
        """
        # read taxonomy file
        taxonomy = Taxonomy().read(taxonomy_file)

        canonical_mapping = read_gtdb_metadata(metadata_file,['formatted_accession'])
        canonical_mapping = {v.formatted_accession: k for k, v in canonical_mapping.items()}

        # a canonical ID is changed to the genome's ID before anything is written;
        # one the metadata file does not map became None, which failed the whole
        # rank's write, and the error was passed over
        genome_ids = {}
        unmapped = []
        for genome_id in taxonomy:
            if re.match(r"G\d{9}", genome_id):
                if genome_id not in canonical_mapping:
                    unmapped.append(genome_id)
                    continue
                genome_ids[genome_id] = canonical_mapping[genome_id]
            else:
                genome_ids[genome_id] = genome_id
        if unmapped:
            raise UnknownGenomesError(
                '{:,} canonical genome ID(s) of {} are not in {}, e.g. {}.'.format(
                    len(unmapped), taxonomy_file, metadata_file,
                    ', '.join(unmapped[:UNKNOWN_EXAMPLES])))

        if truncate_taxonomy:
            self.logger.info('Truncating GTDB taxonomy to domain classification.')
            self.truncate_taxonomy(metadata_file)

        # add each taxonomic rank to database
        gtdbimporter = GTDBImporter(self.temp_cur)
        for i, rank in enumerate(Taxonomy.rank_labels):
            data_to_commit = []
            for genome_id, taxa in taxonomy.items():
                rank_str = taxa[i]
                data_to_commit.append((genome_ids[genome_id], rank_str))

            gtdbimporter.import_metadata_to_db('metadata_taxonomy', 'gtdb_' + rank, 'TEXT', data_to_commit)