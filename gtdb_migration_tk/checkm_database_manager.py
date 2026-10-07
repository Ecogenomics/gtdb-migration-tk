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

"""Writing the CheckM v1 estimates of a release to metadata_genes.

update_checkm_db reads two of the release files checkm writes, checkm.profiles.tsv.gz
and checkm.qa_sh100.tsv.gz, and writes eight fields of metadata_genes from them, in
one transaction. What is written is decided by plan_checkm_import(), which reads
no database and is what the tests drive; the manager only hands the plan to the
importer.

WHICH GENOMES

Until 0.1.48 the command also took --metadata, a table exported from the
database by 'gtdb metadata export', and wrote only the genomes named in its
first column, printing every other genome as skipped. That was what kept a
genome the database does not hold out of upsert(), which fails a whole field for
one. The importer has checked the genomes against genomes.id_at_source itself
since 0.1.47 and refuses input naming one it does not hold, so the export did
nothing the importer does not, and did one thing worse: a genome missing from
the database was dropped with a line on the console rather than refused, and an
export of the wrong release or a stale one dropped genomes the run then called a
success. checkm assesses the genomes report.log says need it, all of which
update_db has put in the database, so CheckM's results are the release's own
genomes and are refused (REFUSE) where the database does not hold one.

READING THE TABLES

The columns are found by name in each table's header, as every table of the
toolkit is read, and a table missing one is refused before anything is written.
A genome is named by its bin, <accession>_protein, and is written as its
accession (ncbi_utils.assembly_accession()), which is how genomes.id_at_source
holds it; the GB_ or RS_ the command used to add only made it match the export.
A genome named twice in a table is refused rather than one value chosen: checkm
assesses each genome in one batch, so a second row says the table was put
together from more than one run. The two tables are the same genomes, the second
qa made over the same CheckM output as the first, and tables that are not are
refused, since a genome would otherwise have one strain heterogeneity of this
release and the other of the last.

ESTIMATES OF AN EARLIER GENOME

update_db keeps a genome's row in genomes, and so its row in metadata_genes, when
its sequences change or a new version takes it over, and the row then holds the
CheckM estimates of sequences the genome no longer has. checkm assesses every
such genome, and update_checkm_db overwrites them, but checkm leaves some genomes
out (no proteins, too large), naming each in checkm_not_assessed.tsv, and a
genome left out kept the estimates of its predecessor as though they were its
own. Such a genome -- named in checkm_not_assessed.tsv, holding CheckM estimates,
and changed by this update (has_changed) -- has the fields update_checkm_db
writes set to NULL, in the same transaction, with a warning, and every one is
written with the estimates it held to checkm_estimates_cleared.tsv beside the
log.

What is cleared is decided from checkm_not_assessed.tsv rather than from every
changed genome missing from the tables: CheckM run again over a few genomes, and
its tables loaded on their own, would then have cleared the estimates of every
other genome of the release. A genome named in checkm_not_assessed.tsv that this
update did not change is not cleared but warned of, since the file is then not
of this release; and one named both there and in the tables is refused.
"""

import logging
import os
from typing import Dict, List, NamedTuple, Sequence, Tuple

from gtdb_migration_tk.database_configuration import GenomeDatabaseConnectionFTPUpdate
from gtdb_migration_tk.database_configuration.GenomeDatabaseConnectionFTPUpdate import one_transaction
from gtdb_migration_tk.biolib_lite.logger import log_directory
from gtdb_migration_tk.gtdb_lite.gtdb_importer import GTDBImporter
from gtdb_migration_tk.ncbi_utils import assembly_accession
from gtdb_migration_tk.utils.common import open_text

# the table every field is written to
CHECKM_TABLE = 'metadata_genes'

# (column of checkm.profiles.tsv.gz, field of metadata_genes, type of the field)
PROFILE_FIELDS = (('Completeness', 'checkm_completeness', 'FLOAT'),
                  ('Contamination', 'checkm_contamination', 'FLOAT'),
                  ('Strain heterogeneity', 'checkm_strain_heterogeneity', 'FLOAT'),
                  ('Marker lineage', 'checkm_marker_lineage', 'TEXT'),
                  ('# genomes', 'checkm_genome_count', 'INT'),
                  ('# markers', 'checkm_marker_count', 'INT'),
                  ('# marker sets', 'checkm_marker_set_count', 'INT'))

# the same of checkm.qa_sh100.tsv.gz, strain heterogeneity at an AAI of 100%
QA_SH100_FIELDS = (('Strain heterogeneity', 'checkm_strain_heterogeneity_100', 'FLOAT'),)

# every field update_checkm_db writes, and clears for a genome checkm left out
CHECKM_FIELDS = tuple(field for _column, field, _type in PROFILE_FIELDS + QA_SH100_FIELDS)

# the column of checkm_not_assessed.tsv naming each genome checkm left out
NOT_ASSESSED_GENOME = 'genome_id'

# beside the log: each genome whose estimates were cleared, and what it held
CLEARED_NAME = 'checkm_estimates_cleared.tsv'

# how many genomes an error names
EXAMPLES = 10

# a field to write: (field, type, [(genome, value), ...])
FieldRows = Tuple[str, str, List[Tuple[str, str]]]


class CheckMTableError(ValueError):
    """A CheckM table that cannot be written to the database as it is."""


class CheckMPlan(NamedTuple):
    """What update_checkm_db writes, decided from checkm's release files."""

    # the genomes of the tables, in the order of the profile
    genomes: List[str]
    # for each field of metadata_genes: (field, type, [(genome, value), ...])
    fields: List[FieldRows]
    # the genomes checkm left out, as checkm_not_assessed.tsv names them
    not_assessed: List[str]


def read_checkm_table(path: str, fields: Sequence[Tuple[str, str, str]]) -> Tuple[List[str], List[FieldRows]]:
    """Read the fields of a CheckM table, each genome named by its accession.

    Parameters
    ----------
    path : str
        Tab-separated CheckM table, gzipped or not, its header naming the bin
        first and then the columns, e.g. checkm.profiles.tsv.gz.
    fields : sequence of (str, str, str)
        Each column to read, the field of metadata_genes it is written to and
        the field's type, e.g. PROFILE_FIELDS.

    @return: the genomes in the order of the table, e.g. ['GCA_000003645.1', ...],
             and for each field (field, type, [(genome, value), ...]).

    Raises
    ------
    CheckMTableError
        The table has no rows, is missing a column, names a genome twice, or
        has a row shorter than its header.
    """

    with open_text(path) as handle:
        header = [column.strip() for column in handle.readline().rstrip('\n').split('\t')]
        missing = [column for column, _field, _type in fields if column not in header]
        if missing:
            raise CheckMTableError('{} has no column {}; its header is: {}.'.format(
                path, ', '.join(repr(column) for column in missing), ', '.join(header)))
        indices = [header.index(column) for column, _field, _type in fields]

        genomes: List[str] = []
        values: List[List[str]] = [[] for _ in fields]
        seen = set()
        twice = []
        for line_number, line in enumerate(handle, start=2):
            if not line.strip():
                continue
            row = [value.strip() for value in line.rstrip('\n').split('\t')]
            if len(row) < len(header):
                raise CheckMTableError('{} line {:,} has {} column(s) where its header has {}.'.format(
                    path, line_number, len(row), len(header)))
            genome = assembly_accession(row[0])
            if genome in seen:
                twice.append(genome)
                continue
            seen.add(genome)
            genomes.append(genome)
            for column_values, index in zip(values, indices):
                column_values.append(row[index])

    if twice:
        raise CheckMTableError('{} names {:,} genome(s) more than once, e.g. {}; checkm assesses '
                               'each genome in one batch, so the table holds more than one run.'.format(
                                   path, len(twice), ', '.join(twice[:EXAMPLES])))
    if not genomes:
        raise CheckMTableError('{} holds no genomes.'.format(path))

    return genomes, [(field, data_type, list(zip(genomes, column_values)))
                     for (_column, field, data_type), column_values in zip(fields, values)]


def read_not_assessed(path: str) -> List[str]:
    """The genomes checkm left out, from its checkm_not_assessed.tsv.

    Parameters
    ----------
    path : str
        checkm_not_assessed.tsv, gzipped or not: genome_id, reason, detail.
        checkm writes it whether or not it left any genome out.

    @return: the genomes, in the order of the file; empty where none was left out.

    Raises
    ------
    CheckMTableError
        The file has no genome_id column.
    """

    with open_text(path) as handle:
        header = [column.strip() for column in handle.readline().rstrip('\n').split('\t')]
        if NOT_ASSESSED_GENOME not in header:
            raise CheckMTableError('{} has no column {}; its header is: {}.'.format(
                path, repr(NOT_ASSESSED_GENOME), ', '.join(header)))
        index = header.index(NOT_ASSESSED_GENOME)
        genomes = []
        for line in handle:
            if line.strip():
                genomes.append(line.rstrip('\n').split('\t')[index].strip())

    return genomes


def plan_checkm_import(checkm_profile_file: str, checkm_qa_sh100_file: str,
                       checkm_not_assessed_file: str) -> CheckMPlan:
    """Every field update_checkm_db writes and the value of each genome; no database is read.

    Parameters
    ----------
    checkm_profile_file : str
        checkm.profiles.tsv.gz: CheckM qa joined with tree_qa -o 2.
    checkm_qa_sh100_file : str
        checkm.qa_sh100.tsv.gz: CheckM qa --aai_strain 0.9999.
    checkm_not_assessed_file : str
        checkm_not_assessed.tsv: the genomes checkm left out.

    @return: the genomes of the tables, for each field of metadata_genes
             (field, type, [(genome, value), ...]), the profile's fields first,
             and the genomes checkm left out.

    Raises
    ------
    CheckMTableError
        A table cannot be read as read_checkm_table() says, the two tables do
        not name the same genomes, or a genome is named both in them and as
        not assessed.
    """

    profile_genomes, profile_fields = read_checkm_table(checkm_profile_file, PROFILE_FIELDS)
    sh100_genomes, sh100_fields = read_checkm_table(checkm_qa_sh100_file, QA_SH100_FIELDS)

    only_profile = sorted(set(profile_genomes) - set(sh100_genomes))
    only_sh100 = sorted(set(sh100_genomes) - set(profile_genomes))
    if only_profile or only_sh100:
        differences = []
        for genomes, path in ((only_profile, checkm_profile_file), (only_sh100, checkm_qa_sh100_file)):
            if genomes:
                differences.append('{:,} only in {}, e.g. {}'.format(
                    len(genomes), path, ', '.join(genomes[:EXAMPLES])))
        raise CheckMTableError('The CheckM tables are not of the same genomes: {}. Both are '
                               'written by one checkm run.'.format('; '.join(differences)))

    not_assessed = read_not_assessed(checkm_not_assessed_file)
    both = sorted(set(not_assessed) & set(profile_genomes))
    if both:
        raise CheckMTableError('{} names {:,} genome(s) the CheckM tables hold estimates of, e.g. {}; '
                               'the files are of different checkm runs.'.format(
                                   checkm_not_assessed_file, len(both), ', '.join(both[:EXAMPLES])))

    return CheckMPlan(profile_genomes, profile_fields + sh100_fields, not_assessed)


class CheckMDatabaseManager(object):

    def __init__(self, database: Dict[str, str]):
        """Initialization.

        Parameters
        ----------
        database : dict
            libpq keywords naming the database (utils.common.database_keywords()).

        @return: None
        """

        self.logger = logging.getLogger('timestamp')

        self.temp_con = GenomeDatabaseConnectionFTPUpdate.GenomeDatabaseConnectionFTPUpdate(database)
        self.temp_con.MakePostgresConnection()
        self.temp_cur = self.temp_con.cursor()

    @one_transaction
    def add_checkm_to_db(self, checkm_profile_file: str, checkm_qa_sh100_file: str,
                         checkm_not_assessed_file: str) -> None:
        """Write the CheckM estimates of a release to metadata_genes, in one transaction.

        The estimates a genome checkm left out holds of an earlier genome are
        cleared (clear_earlier_estimates()).

        Parameters
        ----------
        checkm_profile_file : str
            checkm.profiles.tsv.gz, as checkm writes it for the release.
        checkm_qa_sh100_file : str
            checkm.qa_sh100.tsv.gz, as checkm writes it for the release.
        checkm_not_assessed_file : str
            checkm_not_assessed.tsv, as checkm writes it for the release.

        @return: None

        Raises
        ------
        CheckMTableError
            The tables cannot be written as they are; nothing is written.
        UnknownGenomesError
            A genome of the tables is not in the database; nothing is written.
        """

        plan = plan_checkm_import(checkm_profile_file, checkm_qa_sh100_file, checkm_not_assessed_file)
        self.logger.info('Read the CheckM estimates of {:,} genome(s); checkm left out {:,}.'.format(
            len(plan.genomes), len(plan.not_assessed)))

        importer = GTDBImporter(self.temp_cur)
        for field, data_type, rows in plan.fields:
            importer.import_metadata_to_db(CHECKM_TABLE, field, data_type, rows)
            self.logger.info('Wrote {}.{} for {:,} genome(s).'.format(CHECKM_TABLE, field, len(rows)))

        self.clear_earlier_estimates(plan.not_assessed)

    def clear_earlier_estimates(self, not_assessed: Sequence[str]) -> List[str]:
        """Clear the CheckM estimates a genome checkm left out holds of an earlier genome.

        A genome named in checkm_not_assessed.tsv that this update changed
        (has_changed) and that holds CheckM estimates holds those of the
        sequences it had before, its row having been kept by update_db. Every
        field update_checkm_db writes is set to NULL for it, inside the caller's
        transaction, and each such genome is written with what it held to
        checkm_estimates_cleared.tsv beside the log, which is written whether or
        not there are any.

        Parameters
        ----------
        not_assessed : sequence of str
            The genomes checkm left out, e.g. ['GCA_977065575.1', ...].

        @return: the genomes whose estimates were cleared, sorted.
        """

        held = []
        if not_assessed:
            self.temp_cur.execute(
                'SELECT g.id, g.id_at_source, g.has_changed, {} '
                'FROM genomes g JOIN {} m ON m.id = g.id '
                'WHERE g.id_at_source = ANY(%s) AND ({}) '
                'ORDER BY g.id_at_source'.format(
                    ', '.join('m.' + field for field in CHECKM_FIELDS), CHECKM_TABLE,
                    ' OR '.join('m.{} IS NOT NULL'.format(field) for field in CHECKM_FIELDS)),
                (list(not_assessed),))
            held = self.temp_cur.fetchall()

        unchanged = [row[1] for row in held if not row[2]]
        if unchanged:
            self.logger.warning(
                '{:,} genome(s) checkm_not_assessed.tsv names were not changed by this update, '
                'e.g. {}; their CheckM estimates are kept. The file is not of this release, or '
                'update_db has not been run.'.format(len(unchanged), ', '.join(unchanged[:EXAMPLES])))

        cleared = [row for row in held if row[2]]
        if cleared:
            self.temp_cur.execute(
                'UPDATE {} SET {} WHERE id = ANY(%s)'.format(
                    CHECKM_TABLE, ', '.join('{} = NULL'.format(field) for field in CHECKM_FIELDS)),
                ([row[0] for row in cleared],))

        path = os.path.join(log_directory(), CLEARED_NAME)
        with open(path, 'w') as handle:
            handle.write('\t'.join(('genome_id',) + CHECKM_FIELDS) + '\n')
            for row in cleared:
                handle.write('\t'.join('' if value is None else str(value) for value in row[1:2] + row[3:]) + '\n')

        genomes = [row[1] for row in cleared]
        if genomes:
            self.logger.warning(
                'Cleared the CheckM estimates of {:,} genome(s) checkm left out, which held '
                'those of the sequences they had before this update, e.g. {}. Each is in {}, '
                'with the estimates it held.'.format(len(genomes), ', '.join(genomes[:EXAMPLES]), path))
        else:
            self.logger.info('No genome checkm left out held CheckM estimates of an earlier genome.')

        return genomes
