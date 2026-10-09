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

__author__ = 'Pierre Chaumeil'
__copyright__ = 'Copyright 2019'
__credits__ = ['Pierre Chaumeil']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Pierre Chaumeil'
__email__ = 'p.chaumeil@uq.edu.au'
__status__ = 'Development'

"""Decide which genomes are assembled from type material, from LPSN, the SeqCode and NCBI.

WHERE THE NCBI DATA COMES FROM

Everything type_table knows of a genome it reads from the files of the release,
so it can run once the release is built rather than once the database holds it:
which genomes, from the release's genome_dirs.tsv; the organism name, taxID,
strain identifiers and NCBI's type material status, from the assembly summaries;
and the species, subspecies and every name NCBI gives a taxon, from names.dmp
and nodes.dmp. It read the first five from a metadata table exported from the
database, whose NCBI taxonomy was itself made from names.dmp and nodes.dmp, and
whose strain identifiers ncbi_strains made from the same infraspecific name and
isolate the summaries carry (ncbi_utils.strain_identifiers(), shared by both):
over r237 the summaries give every strain identifier the database held for all
but 51 of 1,047,596 genomes. A species is named as ncbi_tax_manager names it,
NCBI's '(SeqCode)' removed (ncbi_utils.strip_nomenclatural_code()). A genome
whose taxID nodes.dmp does not hold -- NCBI deleted 1,288 taxIDs of r237's
genomes, nearly all of '<genus> sp.' placeholders, between the summaries and the
taxonomy -- has no NCBI species, and so is not type material, as it is in the
NCBI taxonomy of the release.

TYPE MATERIAL UNDER THE SEQCODE

A genome's gtdb_type_designation_ncbi_taxa is decided in three steps, the last
having the final word. LPSN's designation first, from the strain matching. Then
the SeqCode: a genome typing a species valid under it (seqcode download's
seqcode_table.tsv, --seqcode_table) is a type strain of species, and the type
species of its genus where the SeqCode says so. Then NCBI: a genome whose
excluded_from_refseq says it is derived from a metagenome and not used as type
is 'not used as type', with no sources and not the type species of its genus,
gtdb_type_designation_notes saying which designation the rule overrode. The
sources are those agreeing with the designation, 'LPSN' then 'SeqCode', joined
by SOURCE_SEPARATOR.

The last two steps were update_type_designation, run once the summary and the
SeqCode table were loaded, rewriting the two fields in metadata_type_material
from metadata_seqcode and metadata_ncbi. Deciding them here loads them once. It
appended 'Seqcode' to the sources LPSN had given, so a genome LPSN made the type
strain of a heterotypic synonym was 'LPSN;Seqcode' as a type strain of species;
and a genome the NCBI rule made 'not used as type' kept its sources and its
gtdb_type_species_of_genus. Over r237 the rule makes 20 genomes not used as
type, 16 of them LPSN's type strains of species, 5 of those type species of
their genus.

A FAILURE FAILS THE COMMAND

A genome the matching cannot decide is not a genome that is not type material.
Every genome absent from the LPSN summary is written to the release's summary
table as 'not type material', so one dropped there is mislabelled in
metadata_type_material. The genomes are matched on a multiprocessing Pool,
whose worker exceptions are raised in the parent, and a contradiction the
matching meets raises StrainsError rather than calling sys.exit(): the workers
were processes of their own, and a worker that raised or exited took its genome
with it while the command went on and succeeded. Each table is written beside
its final name and moved into place once whole, so a failed run leaves none
that reads as finished.

WARNINGS

What type_table has to warn of is mostly the data rather than the run: a
strain ID NCBI marks '<not considered type>' (545 in r237's names.dmp, a strain
under each spelling NCBI gives it), a genome without an NCBI species. Logged a
line apiece, they buried what else the log had to say. Each is recorded as it is
met (Strains.notice()) -- in the workers too, which hand theirs back with the
genome -- and the run ends with one WARNING per kind, its count and up to three
examples, every one written to type_table_warnings.tsv with what its kind means
and the data that bears on it.
"""

import os
import csv
import datetime
import gzip
import io
import logging
import multiprocessing as mp
import re
from collections import Counter, OrderedDict, defaultdict, namedtuple
from typing import NamedTuple

from tqdm import tqdm

from gtdb_migration_tk.batching import read_genome_dirs
from gtdb_migration_tk.ncbi_utils import (NCBI_NA, REFSEQ_PREFIX, read_summary_rows,
                                          strain_identifiers, strip_nomenclatural_code, summary_field)
from gtdb_migration_tk.taxon_utils import canonical_strain_id, check_format_strain


class StrainsError(RuntimeError):
    """Data the type material matching cannot decide a genome from."""


# The Strains instance a Pool's workers match genomes with. It is set before the
# pool is forked, so each worker inherits the metadata of every genome once,
# rather than having the instance pickled with every genome it is sent.
_MATCHER = None


def _match_genome(gid):
    """Match one genome against the strain repository; a worker of the pool."""
    return _MATCHER.match_genome(gid)


# What a species' type is called in the third column of lpsn_strains.tsv, as
# `lpsn parse_html` writes it from the LPSN web page, and the designation it
# gives. The three are never combined: on the r232 pages no species has more
# than one, each belonging to its own code -- 'Type strain' to names validly
# published under the ICNP, 'Holotype' to the ICN (Botanical Code), and
# 'Nomenclatural type' to names not validly published. A value not here, a
# combination ('Type strain;Holotype') among them, is refused rather than
# guessed at: read as no designation, as it was, it quietly made the species
# 'type strain of species'.
LPSN_TYPE_DESIGNATIONS = {'Type strain': 'type strain of species',
                          'Holotype': 'holotype of species',
                          'Nomenclatural type': 'nomenclatural type of species'}


# The columns of an assembly summary type_table reads, in this order; each is
# required, so a summary that lacks excluded_from_refseq cannot pass over NCBI's
# exclusion of a metagenome unsaid.
SUMMARY_COLUMNS = ('assembly_accession', 'taxid', 'organism_name', 'infraspecific_name',
                   'isolate', 'relation_to_type_material', 'excluded_from_refseq')

# The columns of seqcode download's seqcode_table.tsv type_table reads, by
# name: the release's genome typing each species, the species' status, and
# whether it is the type species of its genus. A species whose status holds
# SEQCODE_VALID is valid under the SeqCode, and its genome a type strain of species.
SEQCODE_GENOME = 'seqcode_type_material_accn'
SEQCODE_STATUS = 'seqcode_species_status'
SEQCODE_TYPE_SPECIES_OF_GENUS = 'seqcode_type_species_of_genus'
SEQCODE_VALID = 'Valid'

# The sources of a type designation, as gtdb_type_designation_ncbi_taxa_sources
# names them, joined by SOURCE_SEPARATOR.
LPSN_SOURCE = 'LPSN'
SEQCODE_SOURCE = 'SeqCode'
SOURCE_SEPARATOR = '; '

# NCBI's exclusion of a genome derived from a metagenome and not used as type:
# such a genome is not used as type, whatever LPSN and the SeqCode say.
NOT_USED_AS_TYPE = 'not used as type'
NCBI_METAGENOME_EXCLUSIONS = ('derived from metagenome', 'not used as type')
NCBI_METAGENOME_RULE = "NCBI's 'derived from metagenome; not used as type' rule"


# The kinds of warning type_table gathers, in the order they are reported: what
# the WARNING line counts them as, and what the warnings TSV says one means.
NOT_CONSIDERED_TYPE = 'not_considered_type'
NO_NCBI_SPECIES = 'no_ncbi_species'
SUBSPECIES_WITHOUT_SUBSP = 'subspecies_without_subsp'
MULTIPLE_TYPE_SPECIES = 'multiple_type_species'
MULTIPLE_PRIORITY_YEARS = 'multiple_priority_years'
LPSN_STRAIN_LINE_WITHOUT_IDS = 'lpsn_strain_line_without_ids'
SEQCODE_TYPE_NOT_IN_RELEASE = 'seqcode_type_not_in_release'
WARNING_KINDS = OrderedDict((
    (NOT_CONSIDERED_TYPE,
     ("strain IDs names.dmp lists as type material but marks '<not considered type>', "
      "left out of their taxa's type material",
      "NCBI's names.dmp lists the strain ID as type material of the taxon but marks it "
      "'<not considered type>'. It is not added to the taxon's type material strain IDs, so "
      "it cannot make a genome type material. NCBI lists a strain under each spelling it "
      "gives it, so one strain can give several of these.")),
    (NO_NCBI_SPECIES,
     ('genomes with no NCBI species, and so not type material',
      "The genome's NCBI taxID is not in nodes.dmp -- NCBI deleted or merged it after the "
      "assembly summaries were written -- or lies below no species node, so the genome has "
      "no NCBI species name to match against LPSN and is not type material.")),
    (SUBSPECIES_WITHOUT_SUBSP,
     ("genomes whose NCBI subspecies-rank name holds no 'subsp.', matched under their species",
      "NCBI places the genome's taxon at subspecies rank under a name with no 'subsp.', "
      "usually a strain name. Its species name is taken from that name, and it is matched "
      "against LPSN under the binomial it begins with.")),
    (MULTIPLE_TYPE_SPECIES,
     ('genera LPSN gives more than one type species, given none',
      "lpsn_species.tsv names more than one type species for the genus. None is taken as its "
      "type species, so no genome of it is flagged gtdb_type_species_of_genus.")),
    (MULTIPLE_PRIORITY_YEARS,
     ('genomes whose names have different LPSN years of priority, the first reported',
      "The names a genome's species is matched under -- its species and its own subspecies, "
      "for one -- have different years of priority in the year table. The first is "
      "reported as its priority_year.")),
    (LPSN_STRAIN_LINE_WITHOUT_IDS,
     ('lines of lpsn_strains.tsv with no strain IDs, ignored',
      "A line of lpsn_strains.tsv names a species and no strain IDs; there is nothing of it "
      "to match a genome to.")),
    (SEQCODE_TYPE_NOT_IN_RELEASE,
     ('genomes of the SeqCode table not in the release, passed over',
      "seqcode_table.tsv names the genome as the type of a species valid under the SeqCode, "
      "and the release's genome_dirs.tsv does not hold it: the table was made for another "
      "release. Nothing is written for it.")),
))
WARNINGS_NAME = 'type_table_warnings.tsv'

# The table of each genome's type material status, gzipped: a quarter of a
# gigabyte of text for r237. update_metadata_db reads it either way.
TYPE_STRAIN_SUMMARY_NAME = 'gtdb_type_strain_summary.tsv.gz'
WARNINGS_HEADER = ('warning_type', 'description', 'warning', 'extra_data')
WARNING_EXAMPLES = 3


class Notice(NamedTuple):
    """A warning type_table gathers, rather than logs, until the run ends."""

    kind: str           # a key of WARNING_KINDS
    message: str        # the warning as a log line would put it
    example: str        # how the WARNING line names it among its examples
    extra_data: str     # 'key=value; ...' of what bears on it


def gtdb_accession(accession):
    """An NCBI accession as GTDB writes it, e.g. RS_GCF_000005845.2 or GB_GCA_000005845.2.

    @return: the accession prefixed with its database.
    """

    return ('RS_' if accession.startswith(REFSEQ_PREFIX) else 'GB_') + accession


def replace_when_written(path):
    """Where a table is written before os.replace() puts it at `path`.

    @return: the temporary path, beside `path`.
    """

    return path + '.tmp'


class Strains(object):
    def __init__(self, output_dir=None, cpus=1):
        """Initialization."""
        self.year = datetime.datetime.now().year

        self.TYPE_SPECIES = 'type strain of species'
        self.NOMENCLATURAL_TYPE='nomenclatural type of species'
        self.HOLOTYPE='holotype of species'
        self.TYPE_NEOTYPE = 'type strain of neotype'
        self.TYPE_SUBSPECIES = 'type strain of subspecies'
        self.TYPE_HETERO_SYNONYM = 'type strain of heterotypic synonym'
        self.NOT_TYPE_MATERIAL = 'not type material'
        self.type_priority = [self.TYPE_SPECIES,
                              self.NOMENCLATURAL_TYPE,
                              self.HOLOTYPE,
                              self.TYPE_NEOTYPE,
                              self.TYPE_SUBSPECIES,
                              self.TYPE_HETERO_SYNONYM,
                              self.NOT_TYPE_MATERIAL]

        self.Match = namedtuple('Match', ['category',
                                          'istype',
                                          'isneotype',
                                          'gtdb_type_status',
                                          'standard_name',
                                          'strain_id',
                                          'year_date',
                                          'is_from_standard'])
        self.logger = logging.getLogger('timestamp')
        self.cpus = cpus
        self.output_dir = output_dir
        self.notices = []

    def notice(self, kind, message, example, **extra):
        """Record a warning, reported with the others of its kind when the run ends.

        @return: None
        """

        self.notices.append(Notice(kind, message, example,
                                   '; '.join('{}={}'.format(k, v) for k, v in extra.items())))

    def report_notices(self, out_dir):
        """One WARNING per kind of warning met, and every one of them in a TSV.

        The TSV is written whether or not there are any, so a run with nothing to
        warn of says so rather than leaving the question open.

        @return: None
        """

        path = os.path.join(out_dir, WARNINGS_NAME)
        by_kind = OrderedDict((kind, []) for kind in WARNING_KINDS)
        for notice in self.notices:
            by_kind[notice.kind].append(notice)

        temporary = replace_when_written(path)
        with open(temporary, 'w', encoding='utf-8') as handle:
            handle.write('\t'.join(WARNINGS_HEADER) + '\n')
            for kind, notices in by_kind.items():
                for notice in notices:
                    handle.write('\t'.join((kind, WARNING_KINDS[kind][1], notice.message,
                                            notice.extra_data)) + '\n')
        os.replace(temporary, path)

        for kind, notices in by_kind.items():
            if not notices:
                continue
            examples = []
            for notice in notices:
                if notice.example not in examples:
                    examples.append(notice.example)
                if len(examples) == WARNING_EXAMPLES:
                    break
            self.logger.warning('{:,} {} (e.g. {}); every one is in {}.'.format(
                len(notices), WARNING_KINDS[kind][0], '; '.join(examples), path))

    def load_year_dict(self, year_table):
        """Load year of priority for species as identified at LPSN."""

        dict_date = {}
        
        with open(year_table) as yt:
            for line in yt:
                infos = line.rstrip('\n').split('\t')
                sp = infos[0]
                year = int(infos[1])
                dict_date[sp] = year
                
        return dict_date

    def standardize_strain_id(self, strain_id):
        """Convert strain ID into standard format."""

        pattern = re.compile('[^A-Za-z0-9/]+')
        strain_id = strain_id.replace('strain', '')
        standardized_id = pattern.sub('', strain_id.strip()).upper()

        return standardized_id

    def fix_common_strain_id_errors(self, strain_ids):
        """Fix common erros associated with NCBI strain IDs."""

        # NCBI strain IDs sometimes contain a 'T' at the end that
        # actually designates the ID is for type material. To
        # resolve this tailing T's are removed and the new ID
        # added as a potential strain ID
        new_ids = set()
        for sid in strain_ids:
            if len(sid) > 1 and sid[-1] == 'T':
                new_ids.add(sid[0:-1])

            new_ids.add(sid)

        return new_ids

    def parse_ncbi_names_and_nodes(self, ncbi_names_file, ncbi_nodes_file, taxids_of_interest):
        """Parse NCBI names.dmp and nodes.dmp files"""

        # determine NCBI taxIDs of species and parent<->child tree
        species_taxids = set()
        subspecies_taxids = set()
        parent = {}
        with open(ncbi_nodes_file) as nodes:
            for line in nodes:
                tokens = [token.strip() for token in line.split('|')]

                cur_taxid = int(tokens[0])
                parent_taxid = int(tokens[1])
                rank = tokens[2]

                parent[cur_taxid] = parent_taxid

                if rank == 'species':
                    species_taxids.add(cur_taxid)
                elif rank == 'subspecies':
                    subspecies_taxids.add(cur_taxid)

        self.logger.info(
            'Identified {:,} NCBI taxonomy species nodes.'.format(len(species_taxids)))

        # determine species taxID of all taxa of interest
        species_of_taxid = {}
        for cur_taxid in taxids_of_interest:
            parent_taxid = cur_taxid
            while True:
                if parent_taxid in species_taxids:
                    species_of_taxid[cur_taxid] = parent_taxid
                    break

                if parent_taxid not in parent:
                    # this happens as not all genomes are defined below
                    # the rank of species and since the NCBI taxonomy and
                    # genome data are not always in sync
                    break

                parent_taxid = parent[parent_taxid]

        self.logger.info(
            'Associated {:,} NCBI taxon nodes with their parent species node.'.format(len(species_of_taxid)))

        # parse auxillary names associated with a NCBI taxID and
        # type material strain IDs for species
        category_names = {}
        type_material = defaultdict(set)
        ncbi_authority = {}
        scientific_names = {}
        with open(ncbi_names_file) as nnf:
            for line in nnf:
                tokens = [token.strip() for token in line.split('|')]
                cur_taxid = int(tokens[0])

                if tokens[3] == 'scientific name' and (cur_taxid in species_taxids
                                                       or cur_taxid in subspecies_taxids):
                    scientific_names[cur_taxid] = strip_nomenclatural_code(tokens[1])

                if '<not considered type>' in tokens[2]:
                    self.notice(NOT_CONSIDERED_TYPE,
                                f"Ignoring {tokens[1]} as it is not considered type material.",
                                '{} (taxid {})'.format(tokens[1], cur_taxid),
                                taxid=cur_taxid, unique_name=tokens[2])
                elif tokens[3] == 'type material':
                    for sid in self.fix_common_strain_id_errors([tokens[1]]):
                        type_material[cur_taxid].add(
                            self.standardize_strain_id(sid))

                if cur_taxid in taxids_of_interest:
                    if tokens[3] == 'authority':
                        ncbi_authority[cur_taxid] = tokens[1]
                    if tokens[3] in ['misspelling', 'synonym', 'equivalent name', 'scientific name']:
                        if cur_taxid not in category_names:
                            category_names[cur_taxid] = {'misspelling': [],
                                                         'synonym': [],
                                                         'equivalent name': [],
                                                         'scientific name': []}
                        category_names[cur_taxid][tokens[3]].append(tokens[1])

        self.logger.info(
            'Read auxillary species name information for {:,} NCBI taxIDs.'.format(len(category_names)))
        self.logger.info(
            'Read type material information for {:,} NCBI taxIDs.'.format(len(type_material)))

        # sanity check results
        for k, v in category_names.items():
            if len(set(v['synonym']).intersection(v.get('scientific name'))) > 0 or len(set(v['synonym']).intersection(v['equivalent name'])) > 0:
                raise StrainsError(
                    'NCBI taxID {} gives a name as both a synonym and a scientific or '
                    'equivalent name: synonyms {}, scientific names {}, equivalent names {}.'.format(
                        k, v['synonym'], v.get('scientific name'), v['equivalent name']))

        # name the species, and any subspecies, each taxon of interest lies in:
        # the subspecies nearest it, and the species above that
        lineage_names = {}
        for cur_taxid in taxids_of_interest:
            if cur_taxid not in parent:
                lineage_names[cur_taxid] = None
                continue
            species = subspecies = None
            taxid = cur_taxid
            while taxid in parent:
                if taxid in species_taxids:
                    species = scientific_names.get(taxid)
                    break
                if taxid in subspecies_taxids and subspecies is None:
                    subspecies = scientific_names.get(taxid)
                if parent[taxid] == taxid:
                    break
                taxid = parent[taxid]
            lineage_names[cur_taxid] = (species, subspecies)

        return category_names, type_material, species_of_taxid, ncbi_authority, lineage_names

    def load_genomes(self, genome_dirs_file, assembly_summary_files):
        """The NCBI data of each genome of the release, from the assembly summaries.

        See WHERE THE NCBI DATA COMES FROM. The species is added once names.dmp
        and nodes.dmp are read (name_species()).

        Parameters
        ----------
        genome_dirs_file : str
            genome_dirs.tsv of the release, naming its genomes.
        assembly_summary_files : sequence of str
            The NCBI assembly summaries the release was selected from.

        @return: (GTDB accession -> the genome's NCBI data, the genomes' taxIDs).

        Raises
        ------
        StrainsError
            A genome of the release is in none of the summaries.
        ncbi_utils.BadInput
            A summary without a column of SUMMARY_COLUMNS.
        """

        release = {accession for accession, _ in read_genome_dirs(genome_dirs_file)}

        metadata = {}
        taxids = set()
        found = set()
        for summary in assembly_summary_files:
            for _, fields, columns in read_summary_rows(summary, required=SUMMARY_COLUMNS):
                (accession, taxid, organism_name, infraspecific_name,
                 isolate, type_material, excluded) = (summary_field(fields, columns, name)
                                                      for name in SUMMARY_COLUMNS)
                if accession not in release:
                    continue
                found.add(accession)

                strain_ids = strain_identifiers(infraspecific_name, isolate)
                standard_strain_ids = {self.standardize_strain_id(sid)
                                       for sid in self.fix_common_strain_id_errors(strain_ids)} - {''}

                metadata[gtdb_accession(accession)] = {
                    'ncbi_organism_name': organism_name,
                    'ncbi_strain_ids': ';'.join(strain_ids) or 'none',
                    'ncbi_standardised_strain_ids': standard_strain_ids,
                    'ncbi_type_material_designation': type_material,
                    'ncbi_excluded_from_refseq': excluded,
                    'ncbi_taxid': int(taxid)}
                taxids.add(int(taxid))

        missing = release - found
        if missing:
            raise StrainsError(
                '{:,} genome(s) of {} are in none of the assembly summaries given, e.g. {}: '
                'the release and the summaries are of different releases, or a summary is '
                'missing.'.format(len(missing), genome_dirs_file, ', '.join(sorted(missing)[:5])))

        self.logger.info('Read the NCBI data of {:,} genomes from {:,} assembly summaries.'.format(
            len(metadata), len(assembly_summary_files)))

        return metadata, taxids

    def name_species(self, lineage_names):
        """Give each genome the species and subspecies of its NCBI taxID.

        @return: None
        """

        for gid, genome_metadata in self.metadata.items():
            taxid = genome_metadata['ncbi_taxid']
            lineage = lineage_names.get(taxid)
            species, subspecies = lineage if lineage else (None, None)
            genome_metadata['ncbi_species'] = species
            genome_metadata['ncbi_subspecies'] = subspecies

            if species is None:
                where = ('is not in nodes.dmp' if lineage is None
                         else 'lies below no species node in nodes.dmp')
                self.notice(NO_NCBI_SPECIES,
                            '{} has no NCBI species: its taxID {} {}.'.format(gid, taxid, where),
                            '{} (taxid {})'.format(gid, taxid),
                            genome=gid, taxid=taxid, in_nodes_dmp='no' if lineage is None else 'yes',
                            ncbi_organism_name=genome_metadata['ncbi_organism_name'])

            # checked here, once a genome, rather than wherever its name is asked for
            if subspecies and 'subsp.' not in subspecies.replace(' pv. ', ' subsp. '):
                self.notice(SUBSPECIES_WITHOUT_SUBSP,
                            "NCBI subspecies name without 'subsp.' definition: {}".format(subspecies),
                            '{} ({})'.format(subspecies, gid),
                            genome=gid, taxid=taxid, ncbi_species=species or '')

    def load_dsmz_strains_dictionary(self, dsmz_dir):
        # We load the dictionary of strains from DSMZ
        dsmz_strains_dic = {}
        pattern = re.compile(r'[\W_]+')
        with open(os.path.join(dsmz_dir, 'dsmz_strains.tsv'), encoding='utf-8') as dsstr:
            dsstr.readline()
            for line in dsstr:
                infos = line.rstrip('\n').split('\t')
                if len(infos) < 2:
                    self.logger.warning('Ignoring a line of dsmz_strains.tsv with no strain IDs: {}'.format(infos))
                else:
                    list_strains = [pattern.sub('', a.strip()).upper(
                    ) for a in infos[1].split('=') if (a != '' and a != 'none')]
                    dsmz_strains_dic[infos[0]] = '='.join(sorted(set(list_strains)))

        return dsmz_strains_dic

    def load_lpsn_strains_dictionary(self, lpsn_dir, lpsn_gss_file):
    
        # get co-identical strain IDs found by scraping LPSN website
        pattern = re.compile('[^A-Za-z0-9/]+')
        lpsn_strains_dic = {}
        unknown_designations = {}
        with open(os.path.join(lpsn_dir, 'lpsn_strains.tsv'), encoding='utf-8') as lpstr:
            lpstr.readline()
            for line in lpstr:
                infos = line.rstrip('\n').split('\t')
                
                sp = infos[0]

                if len(infos) == 1:
                    self.notice(LPSN_STRAIN_LINE_WITHOUT_IDS,
                                'Ignoring a line of lpsn_strains.tsv with no strain IDs: {}'.format(infos),
                                sp, species=sp)
                elif len(infos) == 2:
                    list_strains = [pattern.sub('', a.strip()).upper(
                    ) for a in infos[1].split('=') if (a != '' and a != 'none')]
                    if len(list_strains) > 0:
                        lpsn_strains_dic[sp] = {'strains': '='.join(sorted(set(list_strains))), 'neotypes': ''}
                elif len(infos) == 3:
                    list_strains = [pattern.sub('', a.strip()).upper(
                    ) for a in infos[1].split('=') if (a != '' and a != 'none')]
                    if len(list_strains) > 0:
                        designation = infos[2].strip()
                        if designation and designation not in LPSN_TYPE_DESIGNATIONS:
                            unknown_designations[sp] = designation
                        td = LPSN_TYPE_DESIGNATIONS.get(designation, '')
                        lpsn_strains_dic[sp] = {'strains': '='.join(sorted(set(list_strains))), 'neotypes': '','type_designation': td}

        if unknown_designations:
            raise StrainsError(
                '{:,} species of {} have a type designation strains does not know: {}. It '
                'knows {}, one to a species; a new designation, or a combination of them, '
                'needs deciding what it makes a genome before it is used.'.format(
                    len(unknown_designations), os.path.join(lpsn_dir, 'lpsn_strains.tsv'),
                    '; '.join('{} ({})'.format(sp, d) for sp, d in sorted(unknown_designations.items())[:10]),
                    ', '.join("'{}'".format(d) for d in LPSN_TYPE_DESIGNATIONS)))

        self.logger.info(' - identified strain ids for {:,} species on LPSN website'.format(
                            len(lpsn_strains_dic)))

        # get co-identical strain IDs in LPSN GSS file
        _, lpsn_gss_strain_ids = self.parse_lpsn_gss_metadata(lpsn_gss_file)
        new_strain_ids = 0
        website_strains_only = 0
        for sp, strain_ids in lpsn_gss_strain_ids.items():
            if sp in lpsn_strains_dic:
                scraped_strain_ids = set(lpsn_strains_dic[sp]['strains'].split('='))
                new_strain_ids += len(set(strain_ids) - scraped_strain_ids)
                website_strains_only += len(scraped_strain_ids - set(strain_ids))
                # GSS file is more reliable so defer to these co-identical strain IDs.
                #
                # The web page's type designation is deliberately NOT carried over:
                # every species in the GSS file is validly published under the ICNP,
                # whose types are type strains, so the species is 'type strain of
                # species'. Where its web page says otherwise, it is the page of a
                # different LPSN record of the same name -- not validly published,
                # a basonym, an inaccurate spelling -- which the parse took for the
                # species: 152 species of r232 (Arthrobacter pullicola, Actinoplanes
                # ferrugineus, ...), their web page giving 'Nomenclatural type', and
                # for 19 of them a different type altogether. No species of the GSS
                # file has a 'Holotype' page.
                lpsn_strains_dic[sp] = {'strains': '='.join(strain_ids), 'neotypes': lpsn_strains_dic[sp]['neotypes']}
            else:
                lpsn_strains_dic[sp] = {'strains': '='.join(strain_ids), 'neotypes': ''}



        self.logger.info(' - identified strain ids for {:,} species in LPSN GSS file; deferring to LPSN GSS data whenever possible'.format(
                            len(lpsn_gss_strain_ids)))
        self.logger.info(' - identified {:,} species exclusive to LPSN website'.format(
                            len(set(lpsn_strains_dic) - set(lpsn_gss_strain_ids))))
        self.logger.info(' - identified {:,} species exclusive to LPSN GSS file (ideally zero!)'.format(
                            len(set(lpsn_gss_strain_ids) - set(lpsn_strains_dic))))
        self.logger.info(' - identified {:,} strain IDs exclusive to LPSN GSS file (ideally zero!)'.format(
                            new_strain_ids))
        self.logger.info(' - identified {:,} strain IDs exclusive to LPSN website (ideally zero!)'.format(
                            website_strains_only))


        return lpsn_strains_dic
    
    def _read_type_species_of_genus(self, species_file):
        """Read type species of genus information from DSMZ files."""

        type_species_of_genus = {}
        genus_type_species = {}
        type_species_named = defaultdict(list)
        with open(species_file, encoding='utf-8') as lpstr:
            lpstr.readline()

            for line in lpstr:
                line_split = line.rstrip('\n').split('\t')
                sp, genus, authority , *_ = line_split
                
                if 'Type species of the genus' in authority and not genus:
                    raise StrainsError('{} reads as the type species of a genus, but {} names '
                                       'no genus for it.'.format(sp, species_file))

                if genus:
                    sp = sp.replace('s__', '')
                    type_species_of_genus[sp] = genus

                    if sp not in type_species_named[genus]:
                        type_species_named[genus].append(sp)
                    genus_type_species.setdefault(genus, sp)

        for genus, species in type_species_named.items():
            if len(species) > 1:
                del genus_type_species[genus]
                self.notice(MULTIPLE_TYPE_SPECIES,
                            'Identified multiple type species for {} in {}. Type species for this '
                            'genus will be ignored.'.format(genus, species_file),
                            genus.replace('g__', ''),
                            genus=genus, type_species=', '.join(species))

        return type_species_of_genus, genus_type_species

    def remove_brackets(self, sp_name):
        """Remove brackets from species name.

        e.g., st__[Eubacterium] siraeum 70/3
        """

        if sp_name.startswith('['):
            sp_name = sp_name.replace('[', '', 1).replace(']', '', 1)
        return sp_name

    def get_species_name(self, gid):
        """Determine species name for genome.

        This is the NCBI subspecies name if defined,
        and the species name otherwise.
        """

        ncbi_species = self.metadata[gid]['ncbi_species']
        ncbi_subspecies = self.metadata[gid]['ncbi_subspecies']
        if ncbi_subspecies:
            # fix odd designation impacting less than a dozen genomes
            ncbi_subspecies = ncbi_subspecies.replace(' pv. ', ' subsp. ')

        if ncbi_subspecies:
            return self.remove_brackets(ncbi_subspecies)

        if ncbi_species:
            return self.remove_brackets(ncbi_species)

        return None

    def get_lpsn_priority_year(self, sp):
        """Get year of priority for species according to LPSN."""

        if sp in self.lpsn_year_table:
            return self.lpsn_year_table[sp]

        return ''

    def strains_iterate(self, gid, standard_name, repository_strain_ids, raw_names, misspelling_names, synonyms, equivalent_names, isofficial, sourcest):
        """Check for matching species name and type strain IDs."""

        # search each strain ID at a given strain repository (e.g. LPSN)
        # associated with the species name
        istype = False
        year_date = ''
        category_name = ''
        matched_strain_id = None

        for repository_strain_id in repository_strain_ids.split("="):
            strain_ids = self.metadata[gid]['ncbi_expanded_standardised_strain_ids']
            if repository_strain_id in strain_ids:
                istype = True
            else:
                if len(repository_strain_id) <= 1:
                    continue  # too short to robustly identify

                # remove all white spaces and underscores, and capitalize, before
                # looking for match with standardized strain ID
                pattern = re.compile('[^A-Za-z0-9/]+')
                collapsed_names = {pattern.sub(
                    '', a).upper(): a for a in raw_names}


                for name in collapsed_names:
                    if repository_strain_id in name:
                        first_char = repository_strain_id[0]
                        p = re.compile(" {0}|'{0}".format(first_char), re.IGNORECASE)
                        # search all the matches
                        all_matches = p.findall(collapsed_names[name])
                        # Loop to find all matches using re.search
                        pos = 0
                        index_beginning = []
                        while True:
                            match = p.search(collapsed_names[name], pos)
                            if not match:
                                break
                            index_beginning.append(match.end() - 1)
                            # Update the position to start the next search after the current match
                            pos = match.end()

                        last_char = repository_strain_id[-1]
                        q = re.compile(r'{}(\s|$)'.format(last_char), re.IGNORECASE)
                        # Loop to find all matches using re.search
                        pos = 0
                        index_end = []
                        while True:
                            match = q.search(collapsed_names[name], pos)
                            if not match:
                                break
                            index_end.append(match.start() + 1)
                            # Update the position to start the next search after the current match
                            pos = match.end()

                        if index_beginning and index_end:
                            # we get the substring that matches the matches_beginning and matches_end
                            for idx_beg in index_beginning:
                                for idx_end in index_end:
                                    if idx_end > idx_beg:
                                        potential = collapsed_names[name][idx_beg:idx_end]
                                        potential = canonical_strain_id(potential)
                                        if potential == repository_strain_id:
                                            istype = True
                        # if matches_beginning and matches_end:
                        #     istype = True

            if istype:
                if sourcest == 'lpsn':
                    year_date = self.get_lpsn_priority_year(standard_name)
                else:
                    year_date = ''

                if not isofficial:
                    category_name = self.select_category_name(standard_name,
                                                              misspelling_names,
                                                              synonyms,
                                                              equivalent_names)
                else:
                    category_name = 'official_name'

                matched_strain_id = repository_strain_id
                break

        return (matched_strain_id, category_name, istype, year_date)

    def type_species_or_subspecies(self, gid):
        """Determine if genome is the 'type strain of species' or 'type strain of subspecies'."""

        sp_name = self.get_species_name(gid)
        if 'subsp.' not in sp_name:
            return self.TYPE_SPECIES
        else:
            tokens = sp_name.split()
            subsp_index = tokens.index('subsp.')
            if tokens[subsp_index - 1] == tokens[subsp_index + 1]:
                return self.TYPE_SPECIES

        return self.TYPE_SUBSPECIES

    def match_with_latinization(self, test_sp, target_sp):
        """Check for a match between the specific name of two species considering different gender suffixes."""

        # get specific name from species name
        test = test_sp.split()[1]
        target = target_sp.split()[1]
        if test == target:
            return True

        # determine gender of test name and check for match
        # with related suffix from same group of Latin adjectives
        masc = ('us', 'is', 'er')
        fem = ('a', 'is', 'eris')
        neu = ('um', 'e', 'ere')
        for s, s1, s2 in [(masc, fem, neu), (fem, masc, neu), (neu, masc, fem)]:
            for idx, suffix in enumerate(s):
                if test.endswith(suffix):
                    if test[0:-len(suffix)] + s1[idx] == target:
                        return True
                    elif test[0:-len(suffix)] + s2[idx] == target:
                        return True

        return False

    def check_heterotypic_synonym(self, spe_name, official_spe_names):
        """Check if species is a heterotypic synonym (i.e. has difference specific name)."""

        for official_name in official_spe_names:
            if self.match_with_latinization(spe_name, official_name):
                return False

        return True

    def strain_match(self,
                     gid,
                     standard_names,
                     official_standard_names,
                     misspelling_names,
                     synonyms,
                     equivalent_names,
                     strain_dictionary,
                     sourcest,
                     isofficial):
        """Match species names with stain IDs for a type source (e.g. LPSN) in order to establish if a genome is assembled from type."""

        # Match strain IDs from type sources (e.g. LPSN) associated with each standard
        # species name to strain information at NCBI. Searching is performed on the
        # raw NCBI species designations associated with a standard name as it can be
        # challenging to parse strain information from these entries.
        # The names are taken in sorted order, and a name replaces the match only
        # where its status ranks higher, so that two names matching as well as
        # each other give the same match on every run: they came from a set, in
        # an order of Python's choosing, and the first won. Each name's status
        # is its LPSN designation (designated_status()) before it is compared,
        # so a type strain of a validly published name outranks the nomenclatural
        # type of a name that is not, as they are ranked in type_priority:
        # Clostridium ramosum over Erysipelatoclostridium ramosum, which was
        # matched instead on some runs.
        match = None
        official_status = None
        repository_strain_ids = []
        for standard_name, raw_names in sorted(standard_names.items()):
            if standard_name not in strain_dictionary:
                continue

            if sourcest == 'lpsn':
                # lpsn has information for both strains and neotype strains
                repository_strain_ids = strain_dictionary.get(
                    standard_name).get('strains')

            else:
                repository_strain_ids = strain_dictionary.get(standard_name)
            matched_strain_id, category, istype, year_date = self.strains_iterate(gid,
                                                                                  standard_name,
                                                                                  repository_strain_ids,
                                                                                  raw_names,
                                                                                  misspelling_names,
                                                                                  synonyms,
                                                                                  equivalent_names,
                                                                                  isofficial,
                                                                                  sourcest)


            isneotype = False
            if not istype and sourcest == 'lpsn':
                repository_strain_ids = strain_dictionary.get(
                    standard_name).get('neotypes')
                matched_strain_id, _, isneotype, _ = self.strains_iterate(gid,
                                                                          standard_name,
                                                                          repository_strain_ids,
                                                                          raw_names,
                                                                          misspelling_names,
                                                                          synonyms,
                                                                          equivalent_names,
                                                                          isofficial,
                                                                          sourcest)

            gtdb_type_status = 'not type material'
            if istype or isneotype:
                heterotypic_synonym = False
                if not isofficial:
                    heterotypic_synonym = self.check_heterotypic_synonym(
                        standard_name, official_standard_names)

                if heterotypic_synonym:
                    gtdb_type_status = self.TYPE_HETERO_SYNONYM
                else:
                    gtdb_type_status = self.type_species_or_subspecies(gid)
                    if isneotype and gtdb_type_status == self.TYPE_SPECIES:
                        gtdb_type_status = self.TYPE_NEOTYPE

            m = self.Match(category, istype, isneotype,
                           self.designated_status(gtdb_type_status, standard_name, strain_dictionary),
                           standard_name, matched_strain_id, year_date, False)

            if category == 'official_name':
                # an official name has one type status, its LPSN designation aside
                if official_status is not None and official_status != gtdb_type_status:
                    raise StrainsError('Official species name has ambiguous type status: {}, {}, {}'.format(
                        gid,
                        gtdb_type_status,
                        official_status))
                official_status = gtdb_type_status

            if category != '':
                # it is possible for a genome to be both a 'type strain of subspecies',
                # 'type strain of heterotypic synonym', and potentially a 'type strain of species'
                # depending on the different synonyms, equivalent names, and
                # misspelling
                if (match is None or self.type_priority.index(m.gtdb_type_status)
                        < self.type_priority.index(match.gtdb_type_status)):
                    match = m

        if match :
            # add the is_from_standard field to the match named tuple to indicate if the match is from a standard name or a synonym/equivalent/misspelling name
            # This is based on is_official field and category of the match
            if isofficial and match.category == 'official_name':
                match = match._replace(is_from_standard=True)

        return match

    def designated_status(self, gtdb_type_status, standard_name, strain_dictionary):
        """A type strain of species's status as LPSN designates its type.

        LPSN gives the type of a name that is not validly published as its
        nomenclatural type, and of a name under the ICN as its holotype, which
        rank below a type strain in type_priority.

        @return: NOMENCLATURAL_TYPE or HOLOTYPE where LPSN designates the name's
                 type so, and gtdb_type_status otherwise.
        """

        if gtdb_type_status != self.TYPE_SPECIES:
            return gtdb_type_status

        entry = strain_dictionary.get(standard_name)
        designation = entry.get('type_designation') if isinstance(entry, dict) else None
        if designation == self.NOMENCLATURAL_TYPE:
            return self.NOMENCLATURAL_TYPE
        if designation == self.HOLOTYPE:
            return self.HOLOTYPE
        return gtdb_type_status

    def select_category_name(self, spe_name, misspelling_names, synonyms, equivalent_names):
        """Determine if name is a synonym, equivalent name, or misspelling."""

        # determine source of name giving highest priority to synonyms
        # and lowest priority to misspelling
        if spe_name in synonyms:
            return 'synonyms'
        elif spe_name in equivalent_names:
            return 'equivalent name'
        elif spe_name in misspelling_names:
            return 'misspelling name'

        raise StrainsError(f'Failed to identify category of name: {spe_name}')

    def standardise_names(self, potential_names):
        """Create a standard set of species names, include subsp. designations."""

        standardized = defaultdict(set)
        for raw_name in potential_names:

            # standardize the species name
            standard_name = re.sub(r'(?i)(candidatus\s)', r'', raw_name)
            standard_name = re.sub(r'(?i)(serotype.*)', r'', standard_name)
            standard_name = re.sub(r'(?i)(ser\..*)', r'', standard_name)
            standard_name = re.sub(r'\"|\'', r'', standard_name)
            standard_name = self.remove_brackets(standard_name)
            standard_name = standard_name.strip()

            name_tokens = standard_name.split(' ')

            # check if name is binomial
            if len(name_tokens) != 2:
                if len(name_tokens) == 4 and name_tokens[2] == 'subsp.':
                    standardized[standard_name].add(raw_name)

                    # if the subspecies matches the species name
                    if name_tokens[1] == name_tokens[3]:
                        standardized[' '.join(name_tokens[0:2])].add(raw_name)

                # if the name is longer than 4 words but the 3rd word still
                # subsp, we assume than the 4 first words are a subspecies name
                elif len(name_tokens) >= 4 and name_tokens[2] == 'subsp.':
                    if name_tokens[1] == name_tokens[3]:
                        standardized[' '.join(name_tokens[0:2])].add(raw_name)
                        subsp_name = '{0} {1} subsp. {1}'.format(name_tokens[0],
                                                                 name_tokens[1])
                        standardized[subsp_name].add(raw_name)
                    if name_tokens[1] != name_tokens[3]:
                        standardized[' '.join(name_tokens[0:4])].add(raw_name)
                elif len(name_tokens) >= 2:
                    standardized[' '.join(name_tokens[0:2])].add(raw_name)
                    subsp_name = '{0} {1} subsp. {1}'.format(name_tokens[0],
                                                             name_tokens[1])
                    standardized[subsp_name].add(raw_name)
            else:
                standardized[standard_name].add(standard_name)
                subsp_name = '{0} {1} subsp. {1}'.format(name_tokens[0],
                                                         name_tokens[1])
                standardized[subsp_name].add(raw_name)

        return standardized

    def parse_strains(self, sourcest, strain_dictionary, outfile):
        """Parse information for a single strain resource (e.g., LPSN or DSMZ).

        Every genome is matched, on self.cpus processes, and the table is moved
        into place only once every genome has been: an exception matching any
        one of them is raised here and fails the command (see A FAILURE FAILS
        THE COMMAND).

        @return: None
        """

        global _MATCHER

        self.sourcest = sourcest
        self.strain_dictionary = strain_dictionary
        gids = list(self.metadata)

        temporary = replace_when_written(outfile)
        matched = 0
        try:
            with open(temporary, 'w', encoding='utf-8') as fout:
                self._write_strain_header(fout, sourcest)

                if self.cpus > 1:
                    # forked, so the workers inherit this instance from _MATCHER
                    _MATCHER = self
                    with mp.get_context('fork').Pool(self.cpus) as pool:
                        for data, notices in tqdm(pool.imap(_match_genome, gids, chunksize=256),
                                                  total=len(gids), ncols=100, leave=False,
                                                  desc='Matching genomes to {}'.format(sourcest)):
                            self.notices.extend(notices)
                            if data is not None:
                                self._write_strain_row(fout, data)
                                matched += 1
                else:
                    for gid in tqdm(gids, ncols=100, leave=False,
                                    desc='Matching genomes to {}'.format(sourcest)):
                        data, notices = self.match_genome(gid)
                        self.notices.extend(notices)
                        if data is not None:
                            self._write_strain_row(fout, data)
                            matched += 1

            os.replace(temporary, outfile)
        finally:
            _MATCHER = None
            if os.path.exists(temporary):
                os.remove(temporary)

        self.logger.info('Matched {:,} of {:,} genomes to a species name at {}.'.format(
            matched, len(gids), sourcest))

    def match_genome(self, gid):
        """Determine if a genome is assembled from type material.

        @return: (the row of the strain summary for the genome, or None where it
                 has no species name or matches no name at the strain repository;
                 the warnings met, as Notices).
        """

        genome_metadata = self.metadata[gid]
        # handed back with the genome: a notice recorded in a worker would stay there
        notices = []

        species_name = self.get_species_name(gid)
        if species_name is None:
            return None, notices

        standardized_sp_names = self.standardise_names([species_name])

        # get list of misspellings, synonyms, and equivalent names associated
        # with this genome
        misspelling_names = {}
        synonyms = {}
        equivalent_names = {}
        unofficial_potential_names = set()
        if genome_metadata['ncbi_taxid'] in self.ncbi_auxiliary_names:
            unofficial_potential_names.update(self.ncbi_auxiliary_names[
                genome_metadata['ncbi_taxid']]['misspelling'])
            unofficial_potential_names.update(self.ncbi_auxiliary_names[
                genome_metadata['ncbi_taxid']]['synonym'])
            unofficial_potential_names.update(self.ncbi_auxiliary_names[
                genome_metadata['ncbi_taxid']]['equivalent name'])

            misspelling_names = self.standardise_names(self.ncbi_auxiliary_names[
                genome_metadata['ncbi_taxid']]['misspelling'])
            synonyms = self.standardise_names(self.ncbi_auxiliary_names[
                genome_metadata['ncbi_taxid']]['synonym'])
            equivalent_names = self.standardise_names(self.ncbi_auxiliary_names[
                genome_metadata['ncbi_taxid']]['equivalent name'])

        unofficial_standard_names = self.standardise_names(
            unofficial_potential_names)

        # match species and strain information from NCBI with information
        # at type repository (e.g., LPSN)
        match = self.strain_match(gid,
                                  standardized_sp_names,
                                  standardized_sp_names,
                                  None,
                                  None,
                                  None,
                                  self.strain_dictionary,
                                  self.sourcest,
                                  True)
        list_year_tables = []
        for stdname in standardized_sp_names:
            list_year_tables.append(self.get_lpsn_priority_year(stdname))
        # remove empty entries
        list_year_tables = [y for y in list_year_tables if y != '']
        if len(set(list_year_tables)) > 1:
            names = sorted(standardized_sp_names)
            notices.append(Notice(
                MULTIPLE_PRIORITY_YEARS,
                'Identified multiple different years of priority for {}: {}'.format(names, list_year_tables),
                '{} ({})'.format(gid, ', '.join(str(y) for y in list_year_tables)),
                'genome={}; names={}; years={}'.format(gid, ', '.join(names),
                                                       ', '.join(str(y) for y in list_year_tables))))

        year_date = list_year_tables[0] if len(list_year_tables) > 0 else ''

        if not match:
            # check if any of the auxillary names have a species name
            # and strain ID match with the type repository
            match = self.strain_match(gid,
                                      unofficial_standard_names,
                                      standardized_sp_names,
                                      misspelling_names,
                                      synonyms,
                                      equivalent_names,
                                      self.strain_dictionary,
                                      self.sourcest,
                                      False)

        if match:
            if self.sourcest == 'lpsn':
                # lpsn has information for both strains and neotype strains
                repository_strain_ids = self.strain_dictionary[match.standard_name].get(
                    'strains')
            else:
                repository_strain_ids = self.strain_dictionary[match.standard_name]

            return (gid,
                    species_name,
                    year_date,
                    match.istype,
                    match.isneotype,
                    match.gtdb_type_status,
                    match.category,
                    match.standard_name,
                    match.strain_id,
                    set(repository_strain_ids.split('=')),
                    match.is_from_standard), notices

        return None, notices

    def _write_strain_header(self, fout, sourcest):
        """Write the header of the strain summary of one repository."""

        fout.write(
            'genome\tncbi_organism_name\tncbi_species_name\tncbi_type_designation\tgtdb_type_designation')
        fout.write(
            '\tncbi_base_strain_ids\tncbi_canonical_strain_ids\tmatched_strain_id')
        fout.write(
            '\t{0}_match_type\t{0}_match_name\t{0}_match_strain_id\t{0}_strain_ids'.format(sourcest))
        fout.write('\tmissspellings\tequivalent_names\tsynonyms')
        fout.write('\tneotype\tpriority_year\tis_from_standard_name\n')

    def _write_strain_row(self, fout, data):
        """Report the type material status of one genome."""

        (gid,
         species_name,
         year_date,
         type_strain,
         neotype,
         gtdb_type_status,
         category_name,
         matched_sp_name,
         matched_strain_id,
         repository_strain_ids,
         is_from_standard) = data

        info_genomes = self.metadata[gid]

        misspelling = equivalent_name = synonym = ''
        if info_genomes['ncbi_taxid'] in self.ncbi_auxiliary_names:
            misspelling = '; '.join(
                self.ncbi_auxiliary_names[info_genomes['ncbi_taxid']]['misspelling'])
            equivalent_name = '; '.join(
                self.ncbi_auxiliary_names[info_genomes['ncbi_taxid']]['equivalent name'])
            synonym = '; '.join(
                self.ncbi_auxiliary_names[info_genomes['ncbi_taxid']]['synonym'])

        expanded_ids_str = '; '.join(sorted(self.metadata[gid]['ncbi_expanded_standardised_strain_ids']))
        intersect_ids_str = '; '.join(sorted(
            repository_strain_ids.intersection(self.metadata[gid]['ncbi_expanded_standardised_strain_ids'])))
        repo_ids_str = '; '.join(sorted(repository_strain_ids))

        fout.write(
            f"{gid}\t"
            f"{info_genomes['ncbi_organism_name']}\t"
            f"{species_name}\t"
            f"{info_genomes['ncbi_type_material_designation']}\t"
            f"{gtdb_type_status}\t"
            f"{self.metadata[gid]['ncbi_strain_ids']}\t"
            f"{expanded_ids_str}\t"
            f"{matched_strain_id}\t"
            f"{category_name}\t"
            f"{matched_sp_name}\t"
            f"{intersect_ids_str}\t"
            f"{repo_ids_str}\t"
            f"{misspelling}\t"
            f"{equivalent_name}\t"
            f"{synonym}\t"
            f"{neotype}\t"
            f"{year_date}\t"
            f"{is_from_standard}\n"
        )

    def _parse_strain_summary(self, strain_summary_file):
        """Parse type information from strain repository."""

        StrainInfo = namedtuple('StrainInfo', 'type_designation priority_year')

        strain_info = {}
        with open(strain_summary_file, encoding='utf-8') as f:
            header = f.readline().rstrip().split('\t')

            gid_index = header.index('genome')
            gtdb_type_designation_index = header.index('gtdb_type_designation')
            priority_year_index = header.index('priority_year')

            for line in f:
                line_split = line.rstrip('\n').split('\t')

                gid = line_split[gid_index]
                type_designation = line_split[gtdb_type_designation_index]
                priority_year = line_split[priority_year_index]

                strain_info[gid] = StrainInfo(type_designation, priority_year)

        return strain_info

    def read_seqcode_types(self, seqcode_table):
        """The genomes of the release typing a species valid under the SeqCode.

        See TYPE MATERIAL UNDER THE SEQCODE. A genome the release does not hold
        is a warning (SEQCODE_TYPE_NOT_IN_RELEASE) and passed over.

        Parameters
        ----------
        seqcode_table : str
            seqcode_table.tsv, from seqcode download.

        @return: GTDB accession -> whether the SeqCode makes its species the type
                 species of its genus.

        Raises
        ------
        StrainsError
            The table lacks a column type_table reads.
        """

        types = {}
        with open(seqcode_table, encoding='utf-8') as handle:
            header = handle.readline().rstrip('\n').split('\t')
            missing = [c for c in (SEQCODE_GENOME, SEQCODE_STATUS, SEQCODE_TYPE_SPECIES_OF_GENUS)
                       if c not in header]
            if missing:
                raise StrainsError('{} has no {} column(s): it is not a table of seqcode download.'.format(
                    seqcode_table, ', '.join(missing)))
            genome_index = header.index(SEQCODE_GENOME)
            status_index = header.index(SEQCODE_STATUS)
            type_species_index = header.index(SEQCODE_TYPE_SPECIES_OF_GENUS)

            for line in handle:
                fields = line.rstrip('\n').split('\t')
                if len(fields) < len(header) or SEQCODE_VALID not in fields[status_index]:
                    continue
                gid = gtdb_accession(fields[genome_index])
                if gid not in self.metadata:
                    self.notice(SEQCODE_TYPE_NOT_IN_RELEASE,
                                '{} types a species valid under the SeqCode and is not in the release.'.format(
                                    fields[genome_index]),
                                fields[genome_index], genome=fields[genome_index])
                    continue
                types[gid] = fields[type_species_index] == 'True'

        self.logger.info('{:,} genomes of the release type a species valid under the SeqCode, {:,} the type '
                         'species of its genus.'.format(len(types), sum(types.values())))
        return types

    def type_summary_table(self,
                           ncbi_authority,
                           lpsn_summary_file,
                           lpsn_type_species_of_genus,
                           seqcode_types,
                           summary_table_file):
        """Generate type strain summary file across all strain repositories.

        A genome's designation is LPSN's, then the SeqCode's, then NCBI's
        exclusion of a metagenome not used as type: see TYPE MATERIAL UNDER THE
        SEQCODE.
        """

        # parse strain repository files
        lpsn = self._parse_strain_summary(lpsn_summary_file)

        # write out type strain information for each genome, moved into place
        # once every genome is written
        # no time or file name in the gzip header, so that a run over the same
        # inputs writes the same bytes
        raw = open(replace_when_written(summary_table_file), 'wb')
        fout = io.TextIOWrapper(gzip.GzipFile(filename='', fileobj=raw, mode='wb', mtime=0),
                                encoding='utf-8')
        fout.write(
            "accession\tncbi_species\tncbi_organism_name\tncbi_strain_ids\tncbi_canonical_strain_ids")
        fout.write("\tncbi_taxon_authority\tncbi_type_designation")
        fout.write("\tgtdb_type_designation_ncbi_taxa\tgtdb_type_designation_ncbi_taxa_sources")
        fout.write(
            "\tlpsn_type_designation\tlpsn_priority_year")
        fout.write("\tgtdb_type_species_of_genus\tgtdb_type_designation_notes\n")

        missing_type_at_ncbi = 0
        missing_type_at_gtdb = 0
        agreed_type_of_species = 0
        agreed_type_of_subspecies = 0
        num_type_species_of_genus = 0
        seqcode_only = 0
        not_used_as_type = 0
        overridden = Counter()
        for gid, metadata in self.metadata.items():

            fout.write(gid)

            species_name = self.get_species_name(gid) or ''
            fout.write('\t{}\t{}\t{}\t{}'.format(species_name,
                                             metadata['ncbi_organism_name'],
                                             metadata['ncbi_strain_ids'],
                                             '; '.join(sorted(metadata['ncbi_expanded_standardised_strain_ids']))))

            fout.write('\t{}\t{}'.format(ncbi_authority.get(metadata['ncbi_taxid'], '').replace('"', '~'),
                                     metadata['ncbi_type_material_designation']))

            # GTDB sets the type material designation in a specific priority order
            highest_priority_designation = self.NOT_TYPE_MATERIAL
            for sr in [lpsn]:
                if gid in sr and self.type_priority.index(sr[gid].type_designation) < self.type_priority.index(highest_priority_designation):
                    highest_priority_designation = sr[gid].type_designation
                if highest_priority_designation == self.NOMENCLATURAL_TYPE or highest_priority_designation == self.HOLOTYPE:
                    highest_priority_designation = self.TYPE_SPECIES
            # a species valid under the SeqCode makes its genome a type strain of species
            in_seqcode = gid in seqcode_types
            if in_seqcode:
                if highest_priority_designation != self.TYPE_SPECIES:
                    seqcode_only += 1
                highest_priority_designation = self.TYPE_SPECIES

            type_species_of_genus = False
            canonical_sp_name = ' '.join(species_name.split()[0:2])
            if (highest_priority_designation == 'type strain of species' and
                    (species_name in lpsn_type_species_of_genus or canonical_sp_name in lpsn_type_species_of_genus
                     or seqcode_types.get(gid, False))):
                type_species_of_genus = True

            gtdb_type_sources = []
            for sr_id, sr in [(LPSN_SOURCE, lpsn)]:
                if gid in sr and sr[gid].type_designation == highest_priority_designation:
                    gtdb_type_sources.append(sr_id)
                elif (gid in sr and sr[gid].type_designation in (self.TYPE_SPECIES,self.NOMENCLATURAL_TYPE,self.HOLOTYPE) and
                    highest_priority_designation == self.TYPE_SPECIES):
                    gtdb_type_sources.append(sr_id)
            if in_seqcode:
                gtdb_type_sources.append(SEQCODE_SOURCE)

            # NCBI's exclusion of a metagenome not used as type overrides them all
            notes = ''
            excluded = metadata['ncbi_excluded_from_refseq']
            if all(exclusion in excluded for exclusion in NCBI_METAGENOME_EXCLUSIONS):
                not_used_as_type += 1
                overridden.update(gtdb_type_sources)
                notes = '{} {}.'.format(NCBI_METAGENOME_RULE, 'overrides ' + ' and '.join(gtdb_type_sources)
                                        if gtdb_type_sources else 'applies')
                highest_priority_designation = NOT_USED_AS_TYPE
                gtdb_type_sources = []
                type_species_of_genus = False
            num_type_species_of_genus += type_species_of_genus

            fout.write('\t{}'.format(highest_priority_designation))
            fout.write('\t{}'.format(SOURCE_SEPARATOR.join(gtdb_type_sources)))

            fout.write('\t{}'.format(lpsn[gid].type_designation if gid in lpsn else self.NOT_TYPE_MATERIAL))
            fout.write('\t{}'.format(lpsn[gid].priority_year if gid in lpsn else ''))
            fout.write('\t{}\t{}\n'.format(type_species_of_genus, notes))

            # NCBI's null, NCBI_NA, where it gives no type material status; these
            # counts compared with 'none', which it never holds, so counted nothing
            if metadata['ncbi_type_material_designation'] == NCBI_NA and highest_priority_designation == self.TYPE_SPECIES:
                missing_type_at_ncbi += 1

            if metadata['ncbi_type_material_designation'] == 'assembly from type material' and highest_priority_designation == self.NOT_TYPE_MATERIAL:
                missing_type_at_gtdb += 1

            if metadata['ncbi_type_material_designation'] == 'assembly from type material' and highest_priority_designation == self.TYPE_SPECIES:
                agreed_type_of_species += 1

            if (metadata['ncbi_type_material_designation'] in ['assembly from type material', 'assembly from synonym type material']
                    and highest_priority_designation == self.TYPE_SUBSPECIES):
                sp_tokens = self.get_species_name(gid).split()
                if sp_tokens[2] == 'subsp.' and sp_tokens[1] != sp_tokens[3]:
                    agreed_type_of_subspecies += 1

        self.logger.info(
            'Identified {:,} genomes designated as the type species of genus.'.format(num_type_species_of_genus))
        self.logger.info('The SeqCode makes {:,} genomes a type strain of species that LPSN does not.'.format(
            seqcode_only))
        self.logger.info('{} makes {:,} genomes not used as type, overriding {}.'.format(
            NCBI_METAGENOME_RULE, not_used_as_type,
            ', '.join('{} for {:,}'.format(source, overridden[source])
                      for source in (LPSN_SOURCE, SEQCODE_SOURCE) if overridden[source]) or 'neither LPSN nor the SeqCode'))
        self.logger.info(
            'Genomes that appear to have missing type species information at NCBI: {:,}'.format(missing_type_at_ncbi))
        self.logger.info(
            'Genomes that are only effectively published or erroneously missing type species information at GTDB: {:,}'.format(missing_type_at_gtdb))
        self.logger.info(
            'Genomes where GTDB and NCBI both designate type strain of species: {:,}'.format(agreed_type_of_species))
        self.logger.info(
            'Genomes where GTDB and NCBI both designate type strain of subspecies: {:,}'.format(agreed_type_of_subspecies))

        fout.close()
        raw.close()
        os.replace(replace_when_written(summary_table_file), summary_table_file)

    def expand_ncbi_strain_ids(self, ncbi_coidentical_strain_ids, ncbi_species_of_taxid):
        """Expand set of NCBI co-identical strain IDs associated with each genome."""

        for gid, genome_metadata in self.metadata.items():
            # determine the list of strain IDs at NCBI that are
            # associated with the genome
            strain_ids = genome_metadata['ncbi_standardised_strain_ids']
            ncbi_taxid = genome_metadata['ncbi_taxid']
            if ncbi_taxid in ncbi_coidentical_strain_ids:
                if strain_ids.intersection(ncbi_coidentical_strain_ids[ncbi_taxid]):
                    # expand list of strain IDs to include all co-identical
                    # type material strain IDs specified by the NCBI taxonomy
                    # in names.dmp for this taxon
                    strain_ids = strain_ids.union(
                        ncbi_coidentical_strain_ids[ncbi_taxid])

            # check if genome is associated with a NCBI species node which may have
            # additional relevant co-identical strain IDs
            if ncbi_taxid in ncbi_species_of_taxid:
                ncbi_sp_taxid = ncbi_species_of_taxid[ncbi_taxid]
                if ncbi_sp_taxid in ncbi_coidentical_strain_ids:
                    if strain_ids.intersection(ncbi_coidentical_strain_ids[ncbi_sp_taxid]):
                        # expand list of strain IDs to include all co-identical
                        # type material strain IDs specified by the NCBI taxonomy
                        # in names.dmp for this taxon
                        strain_ids = strain_ids.union(
                            ncbi_coidentical_strain_ids[ncbi_sp_taxid])

            self.metadata[gid]['ncbi_expanded_standardised_strain_ids'] = strain_ids

    def generate_type_strain_table(self,
                                   genome_dirs_file,
                                   assembly_summary_files,
                                   ncbi_names_file,
                                   ncbi_nodes_file,
                                   lpsn_gss_file,
                                   lpsn_dir,
                                   year_table,
                                   seqcode_table):
        """Parse multiple sources to identify genomes assembled from type material.

        seqcode_table is seqcode download's seqcode_table.tsv of the release.
        """

        # initialize data being parsed from file
        self.logger.info('Reading the genomes of the release and their NCBI assembly data.')
        self.metadata, taxids_of_interest = self.load_genomes(genome_dirs_file,
                                                              assembly_summary_files)

        self.logger.info('Reading the genomes typing species valid under the SeqCode.')
        seqcode_types = self.read_seqcode_types(seqcode_table)

        self.logger.info('Parsing year table.')
        self.lpsn_year_table = self.load_year_dict(year_table)

        self.logger.info(
            'Parsing NCBI taxonomy information from names.dmp and nodes.dmp.')
        rtn = self.parse_ncbi_names_and_nodes(
            ncbi_names_file, ncbi_nodes_file, taxids_of_interest)
        (self.ncbi_auxiliary_names,
            ncbi_coidentical_strain_ids,
            ncbi_species_of_taxid,
            ncbi_authority,
            lineage_names) = rtn
        self.name_species(lineage_names)

        # expand set of NCBI co-identical strain IDs associated with each
        # genome
        self.logger.info(
            'Expanding co-identical strain IDs associated with each genome.')
        self.expand_ncbi_strain_ids(
            ncbi_coidentical_strain_ids, ncbi_species_of_taxid)

        # identify genomes assembled from type material
        self.logger.info('Identifying genomes assembled from type material.')
        self.logger.info('Parsing information in LPSN directory.')
        lpsn_strains_dic = self.load_lpsn_strains_dictionary(lpsn_dir,
                                                                lpsn_gss_file)

        self.logger.info('Processing LPSN data.')
        lpsn_summary_file = os.path.join(self.output_dir, 'lpsn_summary.tsv')
        self.parse_strains('lpsn',
                           lpsn_strains_dic,
                           lpsn_summary_file)

        # Deprecated
        #self.logger.info('Parsing information in DSMZ directory.')
        #dsmz_strains_dic = self.load_dsmz_strains_dictionary(dsmz_dir)

        # self.logger.info('Processing DSMZ data.')
        # dsmz_summary_file = os.path.join(
        #     self.output_dir, 'dsmz_summary.tsv')
        # self.parse_strains('dsmz',
        #                    dsmz_strains_dic,
        #                    dsmz_summary_file)

        # # generate global summary file if information was generated from all
        # # sources
        self.logger.info('Reading type species of genus as defined at LPSN.')
        lpsn_type_species_of_genus, lpsn_genus_type_species = self._read_type_species_of_genus(
             os.path.join(lpsn_dir, 'lpsn_species.tsv'))
        # self.logger.info(f' - identified type species for {len(lpsn_genus_type_species):,} genera.')
        #
        # self.logger.info('Reading type species of genus as defined at BacDive.')
        # dsmz_type_species_of_genus, dsmz_genus_type_species = self._read_type_species_of_genus(
        #     os.path.join(dsmz_dir, 'dsmz_species.tsv'))
        # self.logger.info(f' - identified type species for {len(dsmz_genus_type_species):,} genera.')
        #
        # for genus in lpsn_genus_type_species:
        #     if genus in dsmz_genus_type_species:
        #         if lpsn_genus_type_species[genus] != dsmz_genus_type_species[genus]:
        #             self.logger.warning('LPSN and DSMZ disagree on type species of genus for {}: {} {}. Deferring to LPSN.'.format(
        #                                     genus,
        #                                     lpsn_genus_type_species[genus],
        #                                     dsmz_genus_type_species[genus]))
        #             del dsmz_type_species_of_genus[dsmz_genus_type_species[genus]]
        #             del dsmz_genus_type_species[genus]

        self.logger.info(
            'Generating summary type information table across all strain repositories.')
        summary_table_file = os.path.join(
            self.output_dir, TYPE_STRAIN_SUMMARY_NAME)
        self.type_summary_table(ncbi_authority,
                                lpsn_summary_file,
                                lpsn_type_species_of_genus,
                                seqcode_types,
                                summary_table_file)

        self.report_notices(self.output_dir)
        self.logger.info('Done.')
        
    def parse_lpsn_scraped_priorities(self, lpsn_scraped_species_info):
        """Parse year of priority from references scraped from LPSN."""

        priorities = {}
        dup_sp = set()
        with open(lpsn_scraped_species_info) as lsi:
            lsi.readline()
            for line in lsi:
                infos = line.rstrip('\n').split('\t')

                species_authority = infos[2]
                reference_str = species_authority.split(', ')[0]
                references = reference_str.replace('(', '').replace(')', '')
                years = re.sub(r'emend\.[^\d]*\d{4}', '', references)
                years = re.sub(r'ex [^\d]*\d{4}', ' ', years)
                years = re.findall('[1-3][0-9]{3}', years, re.DOTALL)
                years = [int(y) for y in years if int(y) <= datetime.datetime.now().year]
                
                if len(years) == 0:
                    # assume this name is validated through ICN and just take the first 
                    # date given as the year of priority
                    years = re.findall('[1-3][0-9]{3}', references, re.DOTALL)
                    years = [int(y) for y in years if int(y) <= datetime.datetime.now().year]

                sp = infos[0].replace('s__', '')
                if sp in priorities:
                    dup_sp.add(sp)
                priorities[sp.replace('s__', '')] = years[0]

        # We make sure that species and subspecies type species have the same date
        # ie Photorhabdus luminescens and Photorhabdus luminescens subsp.
        # Luminescens
        for k, v in priorities.items():
            infos_name = k.split(' ')
            if len(infos_name) == 2 and '{0} {1} subsp. {1}'.format(infos_name[0], infos_name[1]) in priorities:
                priorities[k] = min(int(v), int(priorities.get(
                    '{0} {1} subsp. {1}'.format(infos_name[0], infos_name[1]))))
            elif len(infos_name) == 4 and infos_name[1] == infos_name[3] and '{} {}'.format(infos_name[0], infos_name[1]) in priorities:
                priorities[k] = min(int(v), int(priorities.get(
                    '{} {}'.format(infos_name[0], infos_name[1]))))
                    
        return priorities, dup_sp
        
    def parse_lpsn_gss_metadata(self, lpsn_gss_file):
        """Get priority and co-identical strain IDs for species and subspecies in LPSN GSS file."""

        priorities = {}
        strain_ids = {}
        illegitimate_names = set()
        with open(lpsn_gss_file, encoding='utf-8', errors='ignore') as f:
            csv_reader = csv.reader(f)

            for line_num, tokens in enumerate(csv_reader):
                if line_num == 0:
                    genus_idx = tokens.index('genus_name')
                    specific_idx = tokens.index('sp_epithet')
                    subsp_idx = tokens.index('subsp_epithet')
                    status_idx = tokens.index('status')
                    author_idx = tokens.index('authors')
                    nom_type_idx = tokens.index('nomenclatural_type')
                else:
                    generic = tokens[genus_idx].strip().replace('"', '')
                    specific = tokens[specific_idx].strip().replace('"', '')
                    subsp = tokens[subsp_idx].strip().replace('"', '')
                    
                    if subsp:
                        taxon = '{} {} subsp. {}'.format(generic, specific, subsp)
                    elif specific:
                        taxon = '{} {}'.format(generic, specific)
                    else:
                        # skip genus entries
                        continue

                    status = tokens[status_idx].strip().replace('"', '')
                    status_tokens = [t.strip() for t in status.split(';')]
                    status_tokens = [tt.strip() for t in status_tokens for tt in t.split(',') ]
                    
                    if 'illegitimate name' in status_tokens:
                        illegitimate_names.add(taxon)
                        if taxon in priorities:
                            continue

                    # get priority references, ignoring references if they are
                    # marked as being a revied name as indicated by a 'ex' or 'emend'
                    # (e.g. Holospora (ex Hafkine 1890) Gromov and Ossipov 1981)
                    ref_str = tokens[author_idx]
                    references = ref_str.replace('(', '').replace(')', '')
                    years = re.sub(r'emend\.[^\d]*\d{4}', '', references)
                    years = re.sub(r'ex [^\d]*\d{4}', ' ', years)
                    years = re.findall('[1-3][0-9]{3}', years, re.DOTALL)
                    years = [int(y) for y in years if int(y) <= datetime.datetime.now().year]

                    if (taxon not in illegitimate_names
                        and taxon in priorities 
                        and years[0] != priorities[taxon]):
                            # conflict that can't be attributed to one of the entries being
                            # considered an illegitimate name
                            self.logger.error('Conflicting priority references for {}: {} {}'.format(
                                                taxon, years, priorities[taxon]))

                    priorities[taxon] = years[0]
                    strain_ids[taxon] = [canonical_strain_id(strain_id) 
                                            for strain_id in tokens[nom_type_idx].split(';')
                                            if check_format_strain(canonical_strain_id(strain_id),)]
        
        return priorities, strain_ids

    def generate_date_table(self, 
                                lpsn_scraped_species_info,
                                lpsn_gss_file, 
                                output_file):
        """Parse priority year from LPSN data."""
        
        self.logger.info('Reading priority references scrapped from LPSN.')
        scraped_sp_priority, dup_scraped_sp = self.parse_lpsn_scraped_priorities(lpsn_scraped_species_info)
        self.logger.info(' - read priority for {:,} species.'.format(len(scraped_sp_priority)))
        if dup_scraped_sp:
            self.logger.info(' - identified {:,} species with duplicate entries. A small number is expected.'.format(
                                    len(dup_scraped_sp)))
        
        self.logger.info('Reading priority references from LPSN GSS file.')
        gss_sp_priority, _ = self.parse_lpsn_gss_metadata(lpsn_gss_file)
        self.logger.info(' - read priority for {:,} species.'.format(len(gss_sp_priority)))
        if dup_scraped_sp:
            self.logger.info(' - {:,} of duplicated scraped species resolved in GSS file.'.format(
                            len(dup_scraped_sp.intersection(gss_sp_priority))))
        
        self.logger.info('Scrapped priority information for {:,} species not in GSS file.'.format(
                            len(set(scraped_sp_priority) - set(gss_sp_priority))))
        self.logger.info('Parsed priority information for {:,} species not on LPSN website.'.format(
                            len(set(gss_sp_priority) - set(scraped_sp_priority))))
                            
        self.logger.info('Writing out year of priority for species giving preference to GSS file.')
        output_file = open(output_file, 'w')
        same_year = 0
        diff_year = 0
        for sp in sorted(set(scraped_sp_priority).union(gss_sp_priority)):
            if sp in gss_sp_priority:
                output_file.write('{}\t{}\n'.format(sp, gss_sp_priority[sp]))
            else:
                output_file.write('{}\t{}\n'.format(sp, scraped_sp_priority[sp]))
                
            if sp in gss_sp_priority and sp in scraped_sp_priority:
                if gss_sp_priority[sp] == scraped_sp_priority[sp]:
                    same_year += 1
                else:
                    diff_year += 1
                    
        self.logger.info(' - same priority year in GSS file and website: {:,}'.format(same_year))
        self.logger.info(' - different priority year in GSS file and website: {:,}'.format(diff_year))
            
        output_file.close()
