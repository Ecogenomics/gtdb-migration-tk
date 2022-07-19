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

import sys
import operator
import logging
import re

import csv
from collections import defaultdict, Counter, namedtuple
from datetime import datetime

from gtdb_migration_tk.strains import Strains
from gtdb_migration_tk.biolib_lite.common import canonical_gid, select_delimiter

csv.field_size_limit(sys.maxsize)

from gtdb_migration_tk.utils.prettytable import PrettyTable


class Tools(object):
    def __init__(self):
        """Initialization."""
        self.logger = logging.getLogger()



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
        new_nested_dict = {}
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

        results = []
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
