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

import gzip
import json
import os
import pickle
import subprocess
import sys
import operator
import logging
import re

import csv
import tempfile
import time
import urllib
from collections import defaultdict, Counter
from datetime import datetime

from tqdm import tqdm
import multiprocessing as mp


from gtdb_migration_tk.strains import Strains
from gtdb_migration_tk.biolib_lite.common import canonical_gid
from gtdblib.util.bio.seq_io import read_seq
from gtdblib.util.shell.filemgmt import select_delimiter, matching_brackets
from gtdblib.util.shell.gtdbshutil import make_sure_path_exists

csv.field_size_limit(sys.maxsize)

from gtdb_migration_tk.utils.prettytable import PrettyTable


class Tools(object):
    def __init__(self):
        """Initialization."""
        self.logger = logging.getLogger()

    class tqdm_log(object):
        """A basic wrapper for the tqdm progress bar. Automatically reports the
        runtime statistics after exit.
        """

        def __init__(self, iterable=None, **kwargs):
            # Setup reporting information.
            self.logger = logging.getLogger('timestamp')
            self.start_ts = None

            # Set default parameters.
            default = {'leave': False,
                       'smoothing': 0.1,
                       'bar_format': '==> Processed {n_fmt}/{total_fmt} {unit}s '
                                     '({percentage:.0f}%) |{bar:15}| [{rate_fmt}, ETA {remaining}]'}
            merged = {**default, **kwargs}
            self.args = merged

            # Instantiate tqdm
            self.tqdm = tqdm(iterable, **merged)



    def compare_metadata(self, old_meta_file, new_meta_file, only_ncbi=False, use_formatted_id = False):

        old_delimiter = select_delimiter(old_meta_file)

        old_nested_dict = {}
        with open(old_meta_file, 'r') as omf:
            old_headers = omf.readline().split(old_delimiter)
            if old_delimiter == ',':
                for line in csv.reader(omf):
                    if (only_ncbi and not line[0].startswith('U_')) or not only_ncbi:
                        id_to_add = line[0]
                        if use_formatted_id:
                            id_to_add = 'G' + id_to_add[7:16]
                        old_nested_dict[id_to_add] = {}
                        for i, j in enumerate(line):
                            old_nested_dict[id_to_add][old_headers[i]] = str(j)
            else:
                count = 0
                for raw_line in omf:
                    line = raw_line.strip('\n').split('\t')
                    count +=1
                    print(count,end='\r')
                    if (only_ncbi and not line[0].startswith('U_')) or not only_ncbi:
                        id_to_add = line[0]
                        if use_formatted_id:
                            id_to_add = 'G' + id_to_add[7:16]
                        old_nested_dict[id_to_add] = {}
                        for i, j in enumerate(line):
                            if j in ['n/a', 'na']:
                                j = 'n/a'
                            if old_headers[i] in ('ncbi_ncrna_count', 'ncbi_cds_count','ncbi_trna_count') and j == 'none':
                                j = '0'
                            if old_headers[i] in ['coding_density','gc_percentage']:
                                j=round(float(j),5)
                            if old_headers[i] == 'ncbi_date' and j !='none':
                                j = datetime.strftime(datetime.strptime(j, '%Y-%m-%d'), '%Y-%m-%d')
                            old_nested_dict[id_to_add][old_headers[i]] = str(j)

        self.logger.info('{} parsed'.format(old_meta_file))

        new_delimiter = select_delimiter(new_meta_file)

        header_summary = {}
        # in the new metadata file
        # we check if the genome id exists, and the columns names exist
        # for each common column name we compare the value for each common
        # genomes and add 1 if they are different
        number_of_genomes = 0
        with open(new_meta_file, 'r') as nmf:
            new_headers = nmf.readline().split(new_delimiter)
            if new_delimiter == ',':
                for line in csv.reader(nmf):
                    id_to_add = line[0]
                    if use_formatted_id:
                        id_to_add = 'G' + id_to_add[7:16]
                    if id_to_add in old_nested_dict:
                        number_of_genomes += 1
                        for i, j in enumerate(line):
                            if new_headers[i] in old_headers:
                                if str(j) != old_nested_dict.get(id_to_add).get(new_headers[i]):
                                    header_summary.setdefault(
                                        new_headers[i], []).append(1)
                                else:
                                    header_summary.setdefault(
                                        new_headers[i], []).append(0)
            else:
                count = 0
                for raw_line in nmf:
                    count +=1
                    print(count,end='\r')
                    line = raw_line.strip('\n').split('\t')
                    id_to_add = line[0]
                    if use_formatted_id:
                        id_to_add = 'G' + id_to_add[7:16]
                    if id_to_add in old_nested_dict:
                        number_of_genomes += 1
                        for i, j in enumerate(line):
                            if new_headers[i] in old_headers:
                                if j in ['n/a','na']:
                                    j='n/a'
                                if new_headers[i] in ('ncbi_ncrna_count','ncbi_cds_count','ncbi_trna_count') and j == 'none':
                                    j='0'
                                if new_headers[i] in ['coding_density', 'gc_percentage']:
                                    j = round(float(j), 5)
                                if new_headers[i] == 'ncbi_date' and j !='none':
                                    j = datetime.strftime(datetime.strptime(j, '%Y-%m-%d'),'%Y-%m-%d')
                                if str(j) != old_nested_dict.get(id_to_add).get(new_headers[i]):
                                    header_summary.setdefault(
                                        new_headers[i], []).append(1)
                                else:
                                    header_summary.setdefault(
                                        new_headers[i], []).append(0)

        for k, v in header_summary.items():
            header_summary[k] = round(100 * float(sum(v)) / len(v), 2)

        sorted_d = sorted(header_summary.items(), key=operator.itemgetter(1))

        # We display all common columns from less changes to more changes
        prettyt = PrettyTable()
        prettyt.field_names = ["Column name", "Difference (%)", "Note"]
        for coln in sorted_d:
            prettyt.add_row([coln[0], coln[1], ''])

        print(prettyt)

        removed_columns = set(old_headers) - \
            set(new_headers)
        new_columns = set(new_headers) - set(old_headers)

        print("Based on {} common genomes.".format(number_of_genomes))

        print("Deprecated columns:")
        for removed_column in removed_columns:
            print("\t- {}".format(removed_column))

        print("New columns:")
        for new_column in new_columns:
            print("\t- {}".format(new_column))

    def compare_selected_data(self, old_meta_file, new_meta_file, metafield, output_file, only_ncbi=False):
        old_delimiter = select_delimiter(old_meta_file)
        old_nested_dict = {}
        with open(old_meta_file, 'r') as omf:
            old_headers = omf.readline().split(old_delimiter)
            if metafield not in old_headers:
                self.logger.error(f'{metafield} is not in {old_meta_file}')
                sys.exit()

            if old_delimiter == ',':
                for line in csv.reader(omf):
                    if (only_ncbi and not line[0].startswith('U_')) or not only_ncbi:
                        old_nested_dict[line[0]] = str(
                            line[old_headers.index(metafield)])
            else:
                for raw_line in omf:
                    line = raw_line.strip('\n').split('\t')
                    if (only_ncbi and not line[0].startswith('U_')) or not only_ncbi:
                        old_nested_dict[line[0]] = str(
                            line[old_headers.index(metafield)])

        new_delimiter = select_delimiter(new_meta_file)
        new_nested_dict = {}
        with open(new_meta_file, 'r') as nmf:
            new_headers = nmf.readline().split(new_delimiter)
            if metafield not in new_headers:
                self.logger.error(f'{metafield} is not in {old_meta_file}')
                sys.exit()
            if new_delimiter == ',':
                for line in csv.reader(nmf):
                    if line[0] in old_nested_dict:
                        new_nested_dict[line[0]] = str(
                            line[new_headers.index(metafield)])
            else:
                for raw_line in nmf:
                    line = raw_line.strip('\n').split('\t')
                    if line[0] in old_nested_dict:
                        new_nested_dict[line[0]] = str(
                            line[new_headers.index(metafield)])

        outf = open(output_file, 'w')
        outf.write('genome_id\told_value\tnew_value\tsimilarity\n')
        for k, v in new_nested_dict.items():
            similarity = 'Identical'
            if v != old_nested_dict.get(k):
                similarity = "Different"
            outf.write('{}\n'.format(
                '\t'.join([k, str(old_nested_dict.get(k)), str(v), similarity])))

        self.logger.info('{} parsed'.format(old_meta_file))


    def parse_ncbi_names_and_nodes(self,ncbi_names_file, ncbi_nodes_file, metadata_file,output_file):
        tax_id_list = []
        genomes_with_taxid = {}
        strain_manager = Strains()

        with open(metadata_file, 'r') as tf:
            header_line = tf.readline()
            headers_infos = header_line.strip().split('\t')
            ncbi_taxid_index = headers_infos.index('ncbi_taxid')

            for line in tf:
                taxid = int(line.strip().split('\t')[ncbi_taxid_index])
                tax_id_list.append(taxid)
                genomes_with_taxid.setdefault(taxid, []).append(
                    line.strip().split(';')[0])

        set_tax_id = set(tax_id_list)
        print('There is {} genomes with taxid'.format(len(tax_id_list)))

        """Parse NCBI names.dmp and nodes.dmp files"""
        # determine NCBI taxIDs of species and parent<->child tree
        species_taxids = set()
        parent = {}
        rank_names = {}
        for line in open(ncbi_nodes_file):
            tokens = [token.strip() for token in line.split('|')]
            cur_taxid = int(tokens[0])
            parent_taxid = int(tokens[1])
            rank = tokens[2]
            parent[cur_taxid] = parent_taxid
            if rank == 'species':
                species_taxids.add(cur_taxid)
            rank_names[cur_taxid] = rank
        print(
            'Identified %d NCBI taxonomy species nodes.' % len(species_taxids))

        # determine species taxID of all taxa of interest
        species_of_taxid = {}
        rank_r95 = {}
        extended_taxid = []
        for it, cur_taxid in enumerate(set_tax_id):
            print('{}/{}'.format(it, len(set_tax_id)))
            parent_taxid = cur_taxid
            rank_r95.setdefault(rank_names.get(
                int(cur_taxid)), []).append(cur_taxid)
            while True:
                if parent_taxid in species_taxids:
                    species_of_taxid[cur_taxid] = parent_taxid
                    if parent_taxid != cur_taxid:
                        if rank_names.get(cur_taxid) not in ['no rank', 'subspecies']:
                            print(rank_names.get(cur_taxid))
                        extended_taxid.append(cur_taxid)
                    break
                if parent_taxid not in parent or parent_taxid == 1:
                    # this happens as not all genomes are defined below
                    # the rank of species and since the NCBI taxonomy and
                    # genome data are not always in sync
                    break
                parent_taxid = parent[parent_taxid]
        print(
            'Associated %d NCBI taxon nodes with their parent species node.' % len(species_of_taxid))
        # print(rank_r95)

        for k, v in rank_r95.items():
            gpertaxid = [genomes_with_taxid.get(i) for i in v]
            gpertaxid = [item for sublist in gpertaxid for item in sublist]
            print("{}: {} Taxids ( {} genomes ) ".format(k, len(v), len(gpertaxid)))

        # parse auxillary names associated with a NCBI taxID and
        # type material strain IDs for species
        category_names = {}
        type_material = defaultdict(set)
        ncbi_authority = {}
        with open(ncbi_names_file) as nnf:
            for line in nnf:
                tokens = [token.strip() for token in line.split('|')]
                cur_taxid = int(tokens[0])

                if tokens[3] == 'type material':
                    for sid in strain_manager.fix_common_strain_id_errors([tokens[1]]):
                        type_material[cur_taxid].add(
                            strain_manager.standardize_strain_id(sid))

                if cur_taxid in set_tax_id:
                    if tokens[3] == 'authority':
                        ncbi_authority[cur_taxid] = tokens[1]
                    if tokens[3] in ['misspelling', 'synonym', 'equivalent name', 'scientific name']:
                        if cur_taxid not in category_names:
                            category_names[cur_taxid] = {'misspelling': [],
                                                         'synonym': [],
                                                         'equivalent name': [],
                                                         'scientific name': []}
                        category_names[cur_taxid][tokens[3]].append(tokens[1])

        print(
            'Read auxillary species name information for %d NCBI taxIDs.' % len(category_names))
        print(
            'Read type material information for %d NCBI taxIDs.' % len(type_material))

        # sanity check results
        for k, v in category_names.items():
            if len(set(v['synonym']).intersection(v.get('scientific name'))) > 0 or len(set(v['synonym']).intersection(v['equivalent name'])) > 0:
                print('ERROR')
                print(v['synonym'])
                print(v.get('scientific name'))
                print(v['equivalent name'])
                sys.exit(-1)

        self.compare_ncbi_strain_ids(metadata_file,extended_taxid, type_material, species_of_taxid,output_file)


    def load_metadata(self,metadata_file):
        """Parse data from GTDB metadata file."""

        metadata = {}
        with open(metadata_file, encoding='utf-8') as metaf:
            headers_line = metaf.readline()
            separator = ','
            if '\t' in headers_line:
                separator = '\t'
            headers = headers_line.rstrip('\n').split(separator)

            gtdb_ncbi_organism_name_index = headers.index('ncbi_organism_name')
            gtdb_ncbi_type_material_designation_index = headers.index(
                'ncbi_type_material_designation')
            gtdb_accession_index = headers.index('accession')
            gtdb_taxonomy_species_name_index = headers.index('ncbi_taxonomy')
            gtdb_strain_identifiers_index = headers.index(
                'ncbi_strain_identifiers')
            gtdb_ncbi_taxonomy_unfiltered_index = headers.index(
                'ncbi_taxonomy_unfiltered')
            gtdb_ncbi_taxid_index = headers.index('ncbi_taxid')

            taxids = set()

            strain_manager = Strains()

            for line in metaf:
                infos = line.rstrip('\n').split(separator)

                if not infos[gtdb_accession_index].startswith('U_'):
                    # standardize NCBI strain IDs
                    standard_strain_ids = []
                    if infos[gtdb_strain_identifiers_index] != 'none':
                        pattern = re.compile('[\W_]+')
                        created_list = [
                            sid.strip() for sid in infos[gtdb_strain_identifiers_index].split(';')]
                        created_list = strain_manager.fix_common_strain_id_errors(
                            created_list)
                        standard_strain_ids = [strain_manager.standardize_strain_id(sid)
                                               for sid in created_list
                                               if (sid != '' and sid != 'none')]

                    if infos[gtdb_ncbi_taxid_index] != 'n/a':
                        metadata[infos[gtdb_accession_index]] = {
                            'ncbi_organism_name': infos[gtdb_ncbi_organism_name_index],
                            'taxonomy_species_name': infos[gtdb_taxonomy_species_name_index].split(';')[6].replace('s__',
                                                                                                                   ''),
                            'ncbi_strain_ids': infos[gtdb_strain_identifiers_index],
                            'ncbi_standardised_strain_ids': set(standard_strain_ids),
                            'ncbi_type_material_designation': infos[gtdb_ncbi_type_material_designation_index],
                            'ncbi_taxonomy_unfiltered': infos[gtdb_ncbi_taxonomy_unfiltered_index],
                            'ncbi_taxid': int(infos[gtdb_ncbi_taxid_index])}

                        taxids.add(int(infos[gtdb_ncbi_taxid_index]))

        return metadata, taxids


    def compare_ncbi_strain_ids(self,metadata_file,extended_taxid, ncbi_coidentical_strain_ids, ncbi_species_of_taxid,output_file):
        """Expand set of NCBI co-identical strain IDs associated with each genome."""
        metadata, taxids_of_interest = self.load_metadata(metadata_file)
        outputfile = open(output_file, 'w')
        outputfile.write(
            'genome\tncbi_taxid\tncbi_standardised_strain_ids\tncbi_expanded_standardised_strain_ids\texpanded\tcoidentical_strains\tcoidentical_intersection\tcoidentical_union\t')
        outputfile.write(
            'associated_with_species\tspecies_intersection\tspecies_union\tsame_collection\n')

        for gid, genome_metadata in metadata.items():
                # determine the list of strain IDs at NCBI that are
                # associated with the genome
            strain_ids = genome_metadata['ncbi_standardised_strain_ids']
            ncbi_taxid = genome_metadata['ncbi_taxid']
            coidentical_interbool = False
            coindentical_interlist = []
            coindentical_unionlist = []
            ncbi_sp_taxid_interbool = False
            ncbi_sp_taxid_interlist = []
            ncbi_sp_taxid_unionlist = []

            if ncbi_taxid in ncbi_coidentical_strain_ids:
                if strain_ids.intersection(ncbi_coidentical_strain_ids[ncbi_taxid]):
                        # expand list of strain IDs to include all co-identical
                        # type material strain IDs specified by the NCBI taxonomy
                        # in names.dmp for this taxon
                    coidentical_interbool = True
                    coindentical_interlist = strain_ids.intersection(
                        ncbi_coidentical_strain_ids[ncbi_taxid])
                    strain_ids = strain_ids.union(
                        ncbi_coidentical_strain_ids[ncbi_taxid])
                    coindentical_unionlist = strain_ids

            # check if genome is associated with a NCBI species node which may have
            # additional relevant co-identical strain IDs
            if ncbi_taxid in ncbi_species_of_taxid:
                ncbi_sp_taxid = ncbi_species_of_taxid[ncbi_taxid]
                if ncbi_sp_taxid in ncbi_coidentical_strain_ids:
                    if strain_ids.intersection(ncbi_coidentical_strain_ids[ncbi_sp_taxid]):
                            # expand list of strain IDs to include all co-identical
                            # type material strain IDs specified by the NCBI taxonomy
                            # in names.dmp for this taxon
                        ncbi_sp_taxid_interbool = True
                        ncbi_sp_taxid_interlist = strain_ids.intersection(
                            ncbi_coidentical_strain_ids[ncbi_sp_taxid])
                        strain_ids = strain_ids.union(
                            ncbi_coidentical_strain_ids[ncbi_sp_taxid])
                        ncbi_sp_taxid_unionlist = strain_ids

            if int(ncbi_taxid) in extended_taxid:
                expanded = False
                if len(genome_metadata['ncbi_standardised_strain_ids']) != len(strain_ids):
                    expanded = True
                list_collections = [
                    ''.join([i for i in s if not i.isdigit()]) for s in strain_ids]
                issuecollec = False
                if len(strain_ids) > 1:
                    countcoll = Counter(list_collections)
                    if countcoll.most_common()[0][1] > 1:
                        issucecollec = True
                outputfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(gid, ncbi_taxid, '='.join(
                    genome_metadata['ncbi_standardised_strain_ids']), '='.join(strain_ids), expanded,
                    coidentical_interbool, '='.join(
                        coindentical_interlist), '='.join(coindentical_unionlist),
                    ncbi_sp_taxid_interbool, '='.join(ncbi_sp_taxid_interlist), '='.join(ncbi_sp_taxid_unionlist), issuecollec))

            metadata[gid]['ncbi_expanded_standardised_strain_ids'] = strain_ids
        outputfile.close()


#    parse_ncbi_names_and_nodes('/srv/db/gtdb/metadata/release95/ncbi/taxonomy/20190725/names.dmp',
#                               '/srv/db/gtdb/metadata/release95/ncbi/taxonomy/20190725/nodes.dmp', 'r95_taxids.lst')


    def compare_markers(self,first_domain_report,second_domain_report,output_file,only_ncbi=False,use_formatted_id=False):
        first_domain_report_dict = {}
        second_domain_report_dict = {}
        outf = open(output_file,'w')
        outf.write("formatted_id\told_id\tnew_id\told_domain\tnew_domain\tsame_domain\told_ar122_marker_percent\tnew_ar122_marker_percent\t")
        outf.write("ar122_marker_percent_diff\told_bac120_marker_percent\tnew_bac120_marker_percent\t")
        outf.write("bac120_marker_percent_diff\n")


        with open(first_domain_report,'r') as fdr:
            headers = fdr.readline().strip('\n').split('\t')
            genomeid_idx = headers.index("Genome Id")
            domain_idx = headers.index("Predicted domain")
            arc_mark_idx = headers.index("Archaeal Marker Percentage")
            bac_mark_idx = headers.index("Bacterial Marker Percentage")
            for line in fdr:
                infos = line.strip('\n').split('\t')
                if only_ncbi and infos[0].startswith('U_'):
                    continue
                else:
                    id = infos[genomeid_idx]
                    if use_formatted_id:
                        id = canonical_gid(id)
                    first_domain_report_dict[id] = {"domain":infos[domain_idx],
                                                    "arc_percent":float(infos[arc_mark_idx]),
                                                    "bac_percent":float(infos[bac_mark_idx]),
                                                    "raw_id":infos[genomeid_idx]}

        with open(second_domain_report,'r') as sdr:
            headers = sdr.readline().strip('\n').split('\t')
            genomeid_idx = headers.index("Genome Id")
            domain_idx = headers.index("Predicted domain")
            arc_mark_idx = headers.index("Archaeal Marker Percentage")
            bac_mark_idx = headers.index("Bacterial Marker Percentage")
            for line in sdr:
                infos = line.strip('\n').split('\t')
                if only_ncbi and infos[0].startswith('U_'):
                    continue
                else:
                    id = infos[genomeid_idx]
                    if use_formatted_id:
                        id = canonical_gid(id)
                    if id in first_domain_report_dict:
                        outf.write(f"{id}\t{first_domain_report_dict.get(id).get('raw_id')}\t{infos[genomeid_idx]}\t"
                              f"{first_domain_report_dict.get(id).get('domain')}\t{infos[domain_idx]}\t"
                              f"{first_domain_report_dict.get(id).get('domain')==infos[domain_idx]}\t"
                              f"{first_domain_report_dict.get(id).get('arc_percent')}\t{infos[arc_mark_idx]}\t"
                              f"{abs(first_domain_report_dict.get(id).get('arc_percent')-float(infos[arc_mark_idx]))}\t"
                              f"{first_domain_report_dict.get(id).get('bac_percent')}\t{infos[bac_mark_idx]}\t"
                              f"{abs(first_domain_report_dict.get(id).get('bac_percent') - float(infos[bac_mark_idx]))}\n")

    def compare_metadata_genome_dir(self,metadata_file,genome_dir_file):
        list_genomes_in_metadata = []
        with open(metadata_file, encoding='utf-8') as metaf:
            headers_line = metaf.readline()
            separator = ','
            if '\t' in headers_line:
                separator = '\t'
            headers = headers_line.rstrip('\n').split(separator)
            gtdb_accession_index = headers.index('accession')
            for line in metaf:
                list_genomes_in_metadata.append(line.strip().split(separator)[gtdb_accession_index].replace('RS_','').replace('GB_',''))

        list_genomes_in_genomedir = {}
        with open(genome_dir_file, encoding='utf-8') as metaf:
            for line in metaf:
                infos=line.strip().split(separator)
                list_genomes_in_genomedir[infos[0]]=infos[1]

        print('Genomes in directory but not in metadata')
        for k,v in list_genomes_in_genomedir.items():
            if k not in list_genomes_in_metadata:
                print(f'{k}\t{v}')

    def get_seqcode_classification(self,rank_order,dict_classi, spe_id):
        string_tax = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
        status_tax = [None] * 7
        type_status_temp = [0] * 7
        type_status = [False] * 7

        for rank in dict_classi:
            if rank.get('rank') in rank_order:
                string_tax[rank_order.index(rank.get('rank'))] = string_tax[
                                                                     rank_order.index(rank.get('rank'))] + rank.get(
                    'name')
                status_tax[rank_order.index(rank.get('rank'))] = rank.get('status_name')
                type_status_temp[rank_order.index(rank.get('rank'))] = int(rank.get('type_accession') or 0)
            if rank.get('rank') == 'genus':
                genus_id = int(rank.get('id') or 0)

        for idx, temp_id in enumerate(type_status_temp[:-2]):
            if temp_id == genus_id:
                type_status[idx] = True
        if spe_id == type_status_temp[-2]:
            type_status[-2] = True

        return ';'.join(string_tax), status_tax, type_status

    def generate_seqcode_table(self, gtdb_genome_path_file, output_dir,cpus=1):
        mapping_dict = self.generate_seqcode_mapping(gtdb_genome_path_file, output_dir,cpus)



        rank_order = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
        raw_fields = ['description', 'formal_styling']
        skip_fields = ['register', 'parent', 'children']
        seq_code_fields = [('type_material_accn', 'type_material'),
                           ('id', 'id'),
                           ('name', 'name'),
                           ('rank', 'rank'),
                           ('species_status', 'status_name'),
                           ('priority_date', 'priority_date'),
                           ('genus_status', 'genus_status'),
                           ('family_status', 'family_status'),
                           ('order_status', 'order_status'),
                           ('class_status', 'class_status'),
                           ('phylum_status', 'phylum_status'),
                           ('type_species_of_genus', 'type_species_of_genus'),
                           ('type_genus_of_family', 'type_genus_of_family'),
                           ('type_genus_of_order', 'type_genus_of_order'),
                           ('type_genus_of_class', 'type_genus_of_class'),
                           ('type_genus_of_phylum', 'type_genus_of_phylum'),
                           ('classification', 'classification'),
                           ('proposed_by', 'proposed_by'),
                           ('created_at', 'created_at'),
                           ('updated_at', 'updated_at'),
                           ('url', 'url')]
        seq_code_order_label = [label for label,scf in seq_code_fields]
        seq_code_order_scf = [scf for label,scf in seq_code_fields]

        output_file = os.path.join(output_dir, 'seqcode_table.tsv')
        outf = open(output_file, 'w')
        outf.write('\t'.join([f'seqcode_{x}' for x in seq_code_order_label])+'\n')


        with urllib.request.urlopen("https://disc-genomics.uibk.ac.at/seqcode/type-genomes.json") as url:
            data = json.load(url)
            generic_fields = []

            for record in data.get("values"):
                if record.get("rank") == 'species':
                    spe_record = [''] * len(seq_code_order_label)
                    with urllib.request.urlopen(record.get("url")) as url_spe:
                        data_spe = json.load(url_spe)
                        #print(record)
                    string_classi, string_status, type_list = self.get_seqcode_classification(rank_order,record.get("classification"),
                                                                                 record.get("id"))
                    # print("seqcode_classification -> " + dict_classi)

                    for label, scf in seq_code_fields:
                        if scf == 'type_material':
                            if scf in record:
                                if 'assembly' in str(record.get(scf)):
                                    spe_record[seq_code_order_label.index(label)] = mapping_dict.get(str(record.get(scf).get('assembly')))
                                elif 'nuccore' in str(record.get(scf)):
                                    spe_record[seq_code_order_label.index(label)] = mapping_dict.get(str(record.get(scf).get('nuccore')))
                                else:
                                    self.logger.error(f'Unknown type material {record.get(scf)}')
                        elif scf == 'proposed_by':
                            #print(f"seqcode_{label} [{scf}]-> {str(data_spe.get(scf).get('citation'))}")
                            spe_record[seq_code_order_label.index(label)] = str(data_spe.get(scf).get('citation'))
                        elif scf == 'classification':
                            #print(f"seqcode_{label} [{scf}]-> {string_classi}")
                            spe_record[seq_code_order_label.index(label)] = string_classi
                        elif scf in ['genus_status', 'family_status', 'order_status', 'class_status', 'phylum_status']:
                            rank = scf.split('_')[0]
                            rk_status = string_status[rank_order.index(rank)]
                            #print(f"seqcode_{label} [{scf}]-> {rk_status}")
                            spe_record[seq_code_order_label.index(label)] = rk_status
                        elif scf in ['type_species_of_genus', 'type_genus_of_family', 'type_genus_of_order',
                                     'type_genus_of_class', 'type_genus_of_phylum']:
                            rank = scf.split('_')[3]
                            rk_type = type_list[rank_order.index(rank)]
                            #print(f"seqcode_{label} [{scf}]-> {rk_type}")
                            spe_record[seq_code_order_label.index(label)] = rk_type
                        else:
                            #print(f"seqcode_{label} [{scf}]-> " + str(record.get(scf)))
                            spe_record[seq_code_order_label.index(label)] = str(record.get(scf))
                    if spe_record[0] != '':
                        print('\t'.join([str(x) for x in spe_record])+'\n')
                        outf.write('\t'.join([str(x) for x in spe_record])+'\n')
        outf.close()



    def generate_seqcode_mapping(self,gtdb_genome_path_file, output_dir,cpus=1):
        sequence_infos = {}
        canonical_sequence_infos = {}
        with open(gtdb_genome_path_file) as f:
            for line in f:
                genome_id,path,canonid = line.strip().split('\t')
                canonical_sequence_infos[canonid] = genome_id

        #if file exists, ask if we should overwrite it
        generate_pkl = False
        # if os.path.exists(os.path.join(output_dir,'seq_accessions.pkl')):
        #     overwrite = input('seq_accessions.pkl already exists. Overwrite? (y/n)')
        #     if overwrite != 'y':
        #         generate_pkl = False

        if generate_pkl:
            genome_to_process = []
            with open(gtdb_genome_path_file) as f:
                for idx,line in enumerate(f):
                    genome_to_process.append((line))

            print(f"number of cpus used:{cpus}")

            # populate worker queue with data to process
            workerQueue = mp.Queue()
            writerQueue = mp.Queue()
            manager = mp.Manager()
            return_list = manager.list()

            for f in genome_to_process:
                workerQueue.put(f)

            for _ in range(cpus):
                workerQueue.put(None)

            try:
                workerProc = [mp.Process(target=self.seqcode_parser_worker,
                                         args=(workerQueue, writerQueue, return_list))
                              for _ in range(cpus)]
                writeProc = mp.Process(target=self.__writerThread,
                                       args=(len(genome_to_process), writerQueue))

                writeProc.start()

                for p in workerProc:
                    p.start()

                for p in workerProc:
                    p.join()

                writerQueue.put(None)
                writeProc.join()

            except:
                for p in workerProc:
                    p.terminate()

                writeProc.terminate()

            list_lines_to_write = [x for x in return_list if x != 'null']

            for key, val in list_lines_to_write:
                sequence_infos.setdefault(key, val)


            start = time.process_time()
            with open(os.path.join(output_dir,'seq_accessions.pkl'), 'wb') as f:
                pickle.dump(sequence_infos, f)
            self.logger.info(f'pickle dump step took {time.process_time() - start} seconds')

        start = time.process_time()
        with open('seq_accessions.pkl', 'rb') as f:
            sequence_infos = pickle.load(f)
        self.logger.info(f'pickle load step took {time.process_time() - start} seconds')

        #if file exists, ask if we should overwrite it
        generate_reverse_pkl = False
        # if os.path.exists(os.path.join(output_dir,'seq_accessions_reverse.pkl')):
        #     overwrite = input('seq_accessions_reverse.pkl already exists. Overwrite? (y/n)')
        #     if overwrite != 'y':
        #         generate_reverse_pkl = False

        if generate_reverse_pkl:
            reverse_sequence_infos = {}
            ct = 0
            length_sequence_infos = len(sequence_infos)
            for gid, infos in sequence_infos.items():
                ct += 1
                print('{}/{}'.format(ct, length_sequence_infos), end='\r')
                if infos[0] is not None:
                    reverse_sequence_infos[infos[0]] = gid
                    if '.' in infos[0]:
                        # we remove the extension
                        reverse_sequence_infos[infos[0].split('.')[0]] = gid

                for acc in infos[1]:
                    reverse_sequence_infos[acc] = gid
                    # we remove the extension
                    if '.' in acc:
                        reverse_sequence_infos[acc.split('.')[0]] = gid
                        if m := re.match(r'^([A-Z]{4}([A-Z]{2})?)([0-9]{6,})$', acc.split('.')[0]):
                            reverse_sequence_infos[m.group(1)] = gid
            print('')

            with open(os.path.join(output_dir,'seq_accessions_reverse.pkl'), 'wb') as fp:
                pickle.dump(reverse_sequence_infos, fp)
            self.logger.info(f'pickle dump step for seq_accessions_reverse took {time.process_time() - start} seconds')

        start = time.process_time()
        with open(os.path.join(output_dir,'seq_accessions_reverse.pkl'), 'rb') as f:
            reverse_sequence_infos = pickle.load(f)
        self.logger.info(f'pickle load step for seq_accessions_reverse took {time.process_time() - start} seconds')

        self.logger.info('reverse_sequence_infos loaded')
        count_found = 0
        count_not_found = 0
        number_or_records = 0

        mapping_dict = {}

        print([ x for x  in reverse_sequence_infos.keys() if 'CP' in x])


        with urllib.request.urlopen("https://disc-genomics.uibk.ac.at/seqcode/type-genomes.json") as url:
            data = json.load(url)
            for record in data.get("values"):
                number_or_records += 1
                if 'assembly' in str(record.get("type_material")):
                    seqcode_assembly_id = str(record.get("type_material").get("assembly"))
                    canonical_seqcode_id = canonical_gid(seqcode_assembly_id)
                    if canonical_seqcode_id in canonical_sequence_infos:
                        mapping_dict[seqcode_assembly_id] = canonical_sequence_infos[canonical_seqcode_id]
                        #self.logger.info(f'{seqcode_assembly_id} -> {canonical_sequence_infos[canonical_seqcode_id]}')
                        count_found += 1
                    else:
                        self.logger.info(f'{seqcode_assembly_id} -> not found')
                        count_not_found += 1

                elif 'nuccore' in str(record.get("type_material")):
                    if m := re.match(r'^([A-Z]{4}([A-Z]{2})?)([0-9]{6,})$',
                                     str(record.get("type_material").get("nuccore"))):
                        mapping_dict[record.get("type_material").get("nuccore")] = reverse_sequence_infos[m.group(1)]
                        count_found += 1
                    elif str(record.get("type_material").get("nuccore")).startswith('CP'):
                        if record.get("type_material").get("nuccore") in reverse_sequence_infos:
                            mapping_dict[record.get("type_material").get("nuccore")] = reverse_sequence_infos[
                                str(record.get("type_material").get("nuccore"))]
                            count_found += 1
                        else:
                            self.logger.info(f'{record.get("type_material").get("nuccore")} -> not found')
                            count_not_found += 1

                    else:
                        self.logger.info(f'{record.get("type_material").get("nuccore")} -> not found')
                        count_not_found += 1

        # Summary
        self.logger.info("Summary:")
        self.logger.info(f'Found {count_found} out of {number_or_records} records ({count_found / number_or_records * 100}%)')
        self.logger.info(
            f'Not found {count_not_found} out of {number_or_records} records ({count_not_found / number_or_records * 100}%)')

        self.logger.info('Done')
        return mapping_dict

    def seqcode_parser_worker(self, queueIn, writerQueue, return_list):
        while True:
            tuple_infos = queueIn.get(block=True, timeout=None)

            if tuple_infos == None:
                break
            line=tuple_infos

            infos = line.strip().split('\t')
            assembly_file = os.path.join(infos[1], os.path.basename(infos[1]) + '_assembly_report.txt')
            with open(assembly_file) as f2:
                wgs_project_id = None
                genbank_accs = []
                for line2 in f2:
                    if line2.startswith('# WGS project:'):
                        wgs_project_id = line2.split(':')[1].strip()
                    if not line2.startswith('#'):
                        gbk_acc = line2.split('\t')[4]
                        if gbk_acc != 'na':
                            genbank_accs.append(gbk_acc)
                return_list.append((infos[0],(wgs_project_id, genbank_accs)))
                writerQueue.put(infos[0])

    def __writerThread(self, numDataItems, writerQueue):
        """Store or write results of worker threads in a single thread."""

        processedItems = 0
        while True:
            a = writerQueue.get(block=True, timeout=None)
            if a == None:
                break

            processedItems += 1
            statusStr = 'Finished processing %d of %d (%.2f%%) items.' % (processedItems,
                                                                          numDataItems,
                                                                          float(processedItems) * 100 / numDataItems)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

        sys.stdout.write('\n')

    def generate_ltp_db(self, csv_file,blastdb_file, ltp_fasta_file,output_directory, output_prefix):
        """Create Living Tree Project (LTP) FASTA and taxonomy file.
        This is a mofiied version of the original script to generate the LTP database ( in /srv/db/silva ).
        Because we want to automate the creation of the LTP metadata file, we do not use the original
        Arb file to export the metadata. Instead, we use the blastdb.fasta and CSV data file downloaded
         from `LTP Website`_.

         .. _LTP Website: https://imedea.uib-csic.es/mmg/ltp/

         """

        # parse metadata file dumped from LTP ARB database
        print('Parsing metadata dumped from LTP ARB database:')

        # parse the CSV file to get fields of interest
        delim = select_delimiter(csv_file)
        info_dict = {}
        with open(csv_file, encoding='utf-8') as csvf:
            csv_reader = csv.reader(csvf, delimiter=delim)

            id_index = 0
            orgname_index = 1
            fulltax_index = 2
            type_strain_index = 4

            for row in csv_reader:
                info_dict[row[id_index]] = {'orgname': row[orgname_index],
                                            'fulltax': row[fulltax_index],
                                            'type_strain': row[type_strain_index]}

        # parse the blastdb.fasta for extra metadata
        delim_blast = select_delimiter(blastdb_file)
        with open(blastdb_file, encoding='utf-8') as blastdbf:
            for line in blastdbf:
                if line.startswith('>'):
                    infos = line.strip().split('\|')
                    print(f' - {infos}', end='\r')
                    fields = matching_brackets(infos[2])
                    strain = None
                    accession = None
                    for field in fields:
                        if field.startswith('accession='):
                            accession = field[10:]
                        if field.startswith('strain='):
                            strain = field[7:]
                    if strain is not None:
                        info_dict[accession]['strain'] = strain


        metadata = {}
        for k,v in  info_dict.items():
            seq_ids = k
            fullname_ltp = v['orgname']
            tax_ltp = v['fulltax']
            type_ltp = v['type_strain']
            strain = v.get('strain','')

            for seq_id in seq_ids.split():
                # LTP uses a single entry for identical 16S sequences with
                # the different sequence IDs seperated by a space, e.g.:
                #   SSMG01000000 SSMG01000241
                metadata[seq_id] = (
                    fullname_ltp, tax_ltp, type_ltp, strain)

        print(f' - identified {len(metadata):,} sequences')

        # create FASTA file and taxonomy file for each LTP sequence
        print('Creating FASTA file and taxonomy file for each LTP sequence:')
        make_sure_path_exists(output_directory)
        fout_fna = open(os.path.join(output_directory,output_prefix + '.fna'), 'w')
        fout_taxonomy = open(os.path.join(output_directory,output_prefix + '_taxonomy.tsv'), 'w')
        num_seqs = 0
        for seq_id, seq, annotation in read_seq(ltp_fasta_file, keep_annotation=True):
            # FASTA header line may specify additional LTP sequences
            # as part of the first token of the annotation, e.g.:
            #  >SSMG01000000 SSMG01000241      Photobacterium lucens   Bacteria;Pseudomonadota;Gammaproteobacteria;Vibrionales;Vibrionaceae;Photobacterium
            seq_ids = [seq_id]
            for token in annotation.split('\t')[0].split():
                if token in metadata:
                    seq_ids.append(token)

            # make sure all sequences have identical LTP metadata
            if len(seq_ids) > 1:
                for seq_id in seq_ids[1:]:
                    assert metadata[seq_ids[0]] == metadata[seq_id]

            for seq_id in seq_ids:
                assert seq_id in metadata

                fullname_ltp, tax_ltp, type_ltp, strain = metadata[seq_id]

                fout_fna.write('>{} {}|{}|{}|{}\n'.format(
                    seq_id,
                    fullname_ltp,
                    tax_ltp,
                    type_ltp,
                    strain))
                fout_fna.write(f'{seq}\n')

                if ';' in type_ltp:
                    print(seq_id)

                type_ltp = type_ltp.replace(';', ' ')
                strain = strain.replace(';', ' ')
                taxonomy_str = f'{tax_ltp};{fullname_ltp};{type_ltp}|{strain}'
                fout_taxonomy.write(f'{seq_id}\t{taxonomy_str}\n')

                num_seqs += 1

        print(f' - identified {num_seqs:,} sequences')

        assert len(metadata) == num_seqs

        fout_fna.close()
        fout_taxonomy.close()

        # Time to run makeblastdb on the new fasta file
        print('Running makeblastdb...')

        cmd_to_run = ['makeblastdb','-in',os.path.join(output_directory,output_prefix + '.fna'),'-dbtype','nucl']
        proc = subprocess.Popen(
            cmd_to_run, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        # print proc.returncode
        if proc.returncode != 0:
            raise RuntimeError("%r failed, status code %s stdout %r stderr %r" % (
                cmd_to_run, proc.returncode, stdout, stderr))

        print('Done.')




def symlink(target, link_name, overwrite=False):
    '''
    Create a symbolic link named link_name pointing to target.
    If link_name exists then FileExistsError is raised, unless overwrite=True.
    When trying to overwrite a directory, IsADirectoryError is raised.
    '''

    if not overwrite:
        os.symlink(target, link_name)
        return

    # os.replace() may fail if files are on different filesystems
    link_dir = os.path.dirname(link_name)

    # Create link to target with temporary filename
    while True:
        temp_link_name = tempfile.mktemp(dir=link_dir)

        # os.* functions mimic as closely as possible system functions
        # The POSIX symlink() returns EEXIST if link_name already exists
        # https://pubs.opengroup.org/onlinepubs/9699919799/functions/symlink.html
        try:
            os.symlink(target, temp_link_name)
            break
        except FileExistsError:
            pass

    # Replace link_name with temp_link_name
    try:
        # Pre-empt os.replace on a directory with a nicer message
        if not os.path.islink(link_name) and os.path.isdir(link_name):
            raise IsADirectoryError(f"Cannot symlink over existing directory: '{link_name}'")
        os.replace(temp_link_name, link_name)
    except:
        if os.path.islink(temp_link_name):
            os.remove(temp_link_name)
        raise

def openfile(filename, mode='rt'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)