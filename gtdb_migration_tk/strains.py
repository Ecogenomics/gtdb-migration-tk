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

"""Decide which genomes are assembled from type material, from LPSN and the NCBI taxonomy.

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
"""

import os
import csv
import datetime
import logging
import multiprocessing as mp
import re
from collections import defaultdict, namedtuple

from tqdm import tqdm

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
        with open(ncbi_names_file) as nnf:
            for line in nnf:
                tokens = [token.strip() for token in line.split('|')]
                cur_taxid = int(tokens[0])

                if '<not considered type>' in tokens[2]:
                    self.logger.warning(f"Ignoring {tokens[1]} as it is not considered type material.")
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

        return category_names, type_material, species_of_taxid, ncbi_authority

    def load_metadata(self, metadata_file):
        """Parse data from GTDB metadata file."""

        metadata = {}
        with open(metadata_file, encoding='utf-8', newline='') as metaf:
            # a tab-separated export is split on tabs, as it always was; a
            # comma-separated one is read as CSV, so a quoted value holding a
            # comma (an organism name, a strain list) stays one field
            headers_line = metaf.readline()
            if '\t' in headers_line:
                headers = headers_line.rstrip('\r\n').split('\t')
                rows = (line.rstrip('\r\n').split('\t') for line in metaf)
            else:
                headers = next(csv.reader([headers_line]))
                rows = csv.reader(metaf)

            gtdb_ncbi_organism_name_index = headers.index('ncbi_organism_name')
            gtdb_ncbi_type_material_designation_index = headers.index(
                'ncbi_type_material_designation')
            gtdb_accession_index = headers.index('accession')
            gtdb_strain_identifiers_index = headers.index(
                'ncbi_strain_identifiers')
            gtdb_ncbi_taxonomy_unfiltered_index = headers.index(
                'ncbi_taxonomy_unfiltered')
            gtdb_ncbi_taxid_index = headers.index('ncbi_taxid')

            taxids = set()
            for infos in rows:
                if not infos:
                    continue

                gid = infos[gtdb_accession_index]

                if not gid.startswith('U_'):
                    # standardize NCBI strain IDs
                    standard_strain_ids = []
                    if infos[gtdb_strain_identifiers_index] != 'none':
                        created_list = [
                            sid.strip() for sid in infos[gtdb_strain_identifiers_index].split(';')]
                        created_list = self.fix_common_strain_id_errors(
                            created_list)
                        standard_strain_ids = [self.standardize_strain_id(sid)
                                               for sid in created_list
                                               if (sid != '' and sid != 'none')]
                    metadata[gid] = {
                        'ncbi_organism_name': infos[gtdb_ncbi_organism_name_index],
                        'ncbi_strain_ids': infos[gtdb_strain_identifiers_index],
                        'ncbi_standardised_strain_ids': set(standard_strain_ids),
                        'ncbi_type_material_designation': infos[gtdb_ncbi_type_material_designation_index],
                        'ncbi_taxonomy_unfiltered': infos[gtdb_ncbi_taxonomy_unfiltered_index],
                        'ncbi_taxid': int(infos[gtdb_ncbi_taxid_index])}

                    taxids.add(int(infos[gtdb_ncbi_taxid_index]))

        return metadata, taxids

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
                    dsmz_strains_dic[infos[0]] = '='.join(set(list_strains))

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
                    self.logger.warning('Ignoring a line of lpsn_strains.tsv with no strain IDs: {}'.format(infos))
                elif len(infos) == 2:
                    list_strains = [pattern.sub('', a.strip()).upper(
                    ) for a in infos[1].split('=') if (a != '' and a != 'none')]
                    if len(list_strains) > 0:
                        lpsn_strains_dic[sp] = {'strains': '='.join(set(list_strains)), 'neotypes': ''}
                elif len(infos) == 3:
                    list_strains = [pattern.sub('', a.strip()).upper(
                    ) for a in infos[1].split('=') if (a != '' and a != 'none')]
                    if len(list_strains) > 0:
                        designation = infos[2].strip()
                        if designation and designation not in LPSN_TYPE_DESIGNATIONS:
                            unknown_designations[sp] = designation
                        td = LPSN_TYPE_DESIGNATIONS.get(designation, '')
                        lpsn_strains_dic[sp] = {'strains': '='.join(set(list_strains)), 'neotypes': '','type_designation': td}

        if unknown_designations:
            raise StrainsError(
                '{:,} species of {} have a type designation strains does not know: {}. It '
                'knows {}, one to a species; a new designation, or a combination of them, '
                'needs deciding what it makes a genome before it is used.'.format(
                    len(unknown_designations), os.path.join(lpsn_dir, 'lpsn_strains.tsv'),
                    '; '.join('{} ({})'.format(sp, d) for sp, d in sorted(unknown_designations.items())[:10]),
                    ', '.join("'{}'".format(d) for d in LPSN_TYPE_DESIGNATIONS)))

        self.logger.info(' - identified strain ids for {:,} species on LPSN website.'.format(
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



        self.logger.info(' - identified strain ids for {:,} species in LPSN GSS file; deferring to LPSN GSS data whenever possible.'.format(
                            len(lpsn_gss_strain_ids)))
        self.logger.info(' - identified {:,} species exclusive to LPSN website.'.format(
                            len(set(lpsn_strains_dic) - set(lpsn_gss_strain_ids))))
        self.logger.info(' - identified {:,} species exclusive to LPSN GSS file (ideally zero!).'.format(
                            len(set(lpsn_gss_strain_ids) - set(lpsn_strains_dic))))
        self.logger.info(' - identified {:,} strain IDs exclusive to LPSN GSS file (ideally zero!).'.format(
                            new_strain_ids))
        self.logger.info(' - identified {:,} strain IDs exclusive to LPSN website. (ideally zero!)'.format(
                            website_strains_only))


        return lpsn_strains_dic
    
    def _read_type_species_of_genus(self, species_file):
        """Read type species of genus information from DSMZ files."""

        type_species_of_genus = {}
        genus_type_species = {}
        multiple_types = set()
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

                    if genus in genus_type_species and genus_type_species[genus] != sp:
                        self.logger.warning('Identified multiple type species for {} in {}. Type species for this genus will be ignored.'.format(genus, species_file))
                        multiple_types.add(genus)
                    else:
                        genus_type_species[genus] = sp
                        
        for genus in multiple_types:
            del genus_type_species[genus]

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

        ncbi_unfiltered_taxa = [
            t.strip() for t in self.metadata[gid]['ncbi_taxonomy_unfiltered'].split(';')]
        ncbi_species = None
        ncbi_subspecies = None
        for taxon in ncbi_unfiltered_taxa:
            if taxon.startswith('s__'):
                ncbi_species = taxon[3:]
            elif taxon.startswith('sb__'):
                ncbi_subspecies = taxon[4:]

                # fix odd designation impacting less than a dozen genomes
                ncbi_subspecies = ncbi_subspecies.replace(' pv. ', ' subsp. ')

        if ncbi_subspecies:
            if 'subsp.' not in ncbi_subspecies:
                self.logger.warning(f"NCBI subspecies name without 'subsp.' definition: {ncbi_subspecies}")

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
        match = None
        gtdb_types = set()
        repository_strain_ids = []
        for standard_name, raw_names in standard_names.items():
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

            m = self.Match(category, istype, isneotype, gtdb_type_status,
                           standard_name, matched_strain_id, year_date,False)

            if category == 'official_name':
                if match:
                    prev_gtdb_type_status = match.gtdb_type_status
                    if prev_gtdb_type_status != gtdb_type_status:
                        raise StrainsError('Official species name has ambiguous type status: {}, {}, {}'.format(
                            gid,
                            gtdb_type_status,
                            prev_gtdb_type_status))
                match = m
            elif category != '':
                # it is possible for a genome to be both a 'type strain of subspecies',
                # 'type strain of heterotypic synonym', and potentially a 'type strain of species'
                # depending on the different synonyms, equivalent names, and
                # misspelling
                if match:
                    prev_gtdb_type_status = match.gtdb_type_status
                    if self.type_priority.index(gtdb_type_status) < self.type_priority.index(prev_gtdb_type_status):
                        match = m
                else:
                    match = m

        if match and match.gtdb_type_status==self.TYPE_SPECIES:
            if strain_dictionary.get(match.standard_name) and strain_dictionary.get(match.standard_name).get('type_designation'):
                td = strain_dictionary.get(match.standard_name).get('type_designation')
                if td == 'nomenclatural type of species':
                    match = match._replace(gtdb_type_status=self.NOMENCLATURAL_TYPE)
                elif td == 'holotype of species':
                    match = match._replace(gtdb_type_status=self.HOLOTYPE)

        if match :
            # add the is_from_standard field to the match named tuple to indicate if the match is from a standard name or a synonym/equivalent/misspelling name
            # This is based on is_official field and category of the match
            if isofficial and match.category == 'official_name':
                match = match._replace(is_from_standard=True)

        return match

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
                        for data in tqdm(pool.imap(_match_genome, gids, chunksize=256),
                                         total=len(gids), ncols=100, leave=False,
                                         desc='Matching genomes to {}'.format(sourcest)):
                            if data is not None:
                                self._write_strain_row(fout, data)
                                matched += 1
                else:
                    for gid in tqdm(gids, ncols=100, leave=False,
                                    desc='Matching genomes to {}'.format(sourcest)):
                        data = self.match_genome(gid)
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

        @return: the row of the strain summary for the genome, or None where it
                 has no species name or matches no name at the strain repository.
        """

        genome_metadata = self.metadata[gid]

        species_name = self.get_species_name(gid)
        if species_name is None:
            return None

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
            self.logger.warning('Identified multiple different years of priority for {}: {}'.format(
                dict(standardized_sp_names), list_year_tables))

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
                    match.is_from_standard)

        return None

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

        expanded_ids_str = '; '.join(self.metadata[gid]['ncbi_expanded_standardised_strain_ids'])
        intersect_ids_str = '; '.join(
            repository_strain_ids.intersection(self.metadata[gid]['ncbi_expanded_standardised_strain_ids']))
        repo_ids_str = '; '.join(repository_strain_ids)

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

    def type_summary_table(self,
                           ncbi_authority,
                           lpsn_summary_file,
                           lpsn_type_species_of_genus,
                           summary_table_file):
        """Generate type strain summary file across all strain repositories."""

        # parse strain repository files
        lpsn = self._parse_strain_summary(lpsn_summary_file)

        # write out type strain information for each genome, moved into place
        # once every genome is written
        fout = open(replace_when_written(summary_table_file), 'w')
        fout.write(
            "accession\tncbi_species\tncbi_organism_name\tncbi_strain_ids\tncbi_canonical_strain_ids")
        fout.write("\tncbi_taxon_authority\tncbi_type_designation")
        fout.write("\tgtdb_type_designation_ncbi_taxa\tgtdb_type_designation_ncbi_taxa_sources")
        fout.write(
            "\tlpsn_type_designation\tlpsn_priority_year")
        fout.write("\tgtdb_type_species_of_genus\tis_from_standard\n")

        missing_type_at_ncbi = 0
        missing_type_at_gtdb = 0
        agreed_type_of_species = 0
        agreed_type_of_subspecies = 0
        num_type_species_of_genus = 0
        for gid, metadata in self.metadata.items():

            fout.write(gid)

            species_name = self.get_species_name(gid)
            fout.write('\t{}\t{}\t{}\t{}'.format(species_name,
                                             metadata['ncbi_organism_name'],
                                             metadata['ncbi_strain_ids'],
                                             '; '.join(metadata['ncbi_expanded_standardised_strain_ids'])))

            fout.write('\t{}\t{}'.format(ncbi_authority.get(metadata['ncbi_taxid'], '').replace('"', '~'),
                                     metadata['ncbi_type_material_designation']))

            # GTDB sets the type material designation in a specific priority order
            highest_priority_designation = self.NOT_TYPE_MATERIAL
            for sr in [lpsn]:
                if gid in sr and self.type_priority.index(sr[gid].type_designation) < self.type_priority.index(highest_priority_designation):
                    highest_priority_designation = sr[gid].type_designation
                if highest_priority_designation == self.NOMENCLATURAL_TYPE or highest_priority_designation == self.HOLOTYPE:
                    highest_priority_designation = self.TYPE_SPECIES
            fout.write('\t{}'.format(highest_priority_designation))

            type_species_of_genus = False
            canonical_sp_name = ' '.join(species_name.split()[0:2])
            if (highest_priority_designation == 'type strain of species' and
                    (species_name in lpsn_type_species_of_genus or canonical_sp_name in lpsn_type_species_of_genus )):
                type_species_of_genus = True
                num_type_species_of_genus += 1

            gtdb_type_sources = []
            for sr_id, sr in [('LPSN', lpsn)]:
                if gid in sr and sr[gid].type_designation == highest_priority_designation:
                    gtdb_type_sources.append(sr_id)
                elif (gid in sr and sr[gid].type_designation in (self.TYPE_SPECIES,self.NOMENCLATURAL_TYPE,self.HOLOTYPE) and
                    highest_priority_designation == self.TYPE_SPECIES):
                    gtdb_type_sources.append(sr_id)
            fout.write('\t{}'.format('; '.join(gtdb_type_sources)))

            fout.write('\t{}'.format(lpsn[gid].type_designation if gid in lpsn else self.NOT_TYPE_MATERIAL))
            fout.write('\t{}'.format(lpsn[gid].priority_year if gid in lpsn else ''))
            fout.write('\t{}\n'.format(type_species_of_genus))

            if metadata['ncbi_type_material_designation'] == 'none' and highest_priority_designation == self.TYPE_SPECIES:
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
        self.logger.info(
            'Genomes that appear to have missing type species information at NCBI: {:,}'.format(missing_type_at_ncbi))
        self.logger.info(
            'Genomes that are only effectively published or erroneously missing type species information at GTDB: {:,}'.format(missing_type_at_gtdb))
        self.logger.info(
            'Genomes where GTDB and NCBI both designate type strain of species: {:,}'.format(agreed_type_of_species))
        self.logger.info(
            'Genomes where GTDB and NCBI both designate type strain of subspecies: {:,}'.format(agreed_type_of_subspecies))

        fout.close()
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
                                   metadata_file,
                                   ncbi_names_file,
                                   ncbi_nodes_file,
                                   lpsn_gss_file,
                                   lpsn_dir,
                                   year_table):
        """Parse multiple sources to identify genomes assembled from type material."""

        # initialize data being parsed from file
        self.logger.info('Parsing GTDB metadata.')
        self.metadata, taxids_of_interest = self.load_metadata(metadata_file)

        self.logger.info('Parsing year table.')
        self.lpsn_year_table = self.load_year_dict(year_table)

        self.logger.info(
            'Parsing NCBI taxonomy information from names.dmp and nodes.dmp.')
        rtn = self.parse_ncbi_names_and_nodes(
            ncbi_names_file, ncbi_nodes_file, taxids_of_interest)
        (self.ncbi_auxiliary_names,
            ncbi_coidentical_strain_ids,
            ncbi_species_of_taxid,
            ncbi_authority) = rtn

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
            self.output_dir, 'gtdb_type_strain_summary.tsv')
        self.type_summary_table(ncbi_authority,
                                lpsn_summary_file,
                                lpsn_type_species_of_genus,
                                summary_table_file)

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
