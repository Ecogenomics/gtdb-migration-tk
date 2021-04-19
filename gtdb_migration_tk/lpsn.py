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

import os
import glob
import sys
import urllib.request
import re
import unicodedata
import string
import html
import logging
import pandas as pd
import csv

from sqlalchemy import create_engine
from bs4 import BeautifulSoup
from gtdb_migration_tk.biolib_lite.common import make_sure_path_exists, clean_html
from gtdb_migration_tk.taxon_utils import canonical_strain_id, check_format_strain


class LPSN(object):
    def __init__(self, skip_taxa_per_letter_dl, lpsn_output_dir):
        """Initialization."""

        self.outdir = lpsn_output_dir
        if self.outdir is not None and not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        self.logger = logging.getLogger('timestamp')

        self.base_url = 'https://lpsn.dsmz.de/'
        self.skip_taxa_per_letter_dl = skip_taxa_per_letter_dl

    def full_lpsn_wf(self):
        for rk in ['phylum', 'class', 'order', 'family', 'genus', 'species']:
            self.download_rank_lpsn_html(rk)
        self.download_subspecies_lpsn_html()
        self.parse_html(os.path.join(self.outdir, 'genus_html'))

    def download_rank_lpsn_html(self, rank_name):
        """

        Download all HTML pages from LPSN for a specific rank.

        """
        # Download pages listing all rank_of_interest in LPSN
        index_dir = os.path.join(self.outdir, f'{rank_name}_per_letter')
        if self.skip_taxa_per_letter_dl:
            self.logger.info(f'Skipping download of {rank_name} at LPSN and using results in: {index_dir}')
        else:
            self.logger.info(f'Beginning download of {rank_name} from LPSN.')
            make_sure_path_exists(index_dir)
            for letter in list(string.ascii_uppercase):
                url = self.base_url + f'{rank_name}?page=' + letter
                urllib.request.urlretrieve(
                    url, os.path.join(index_dir, '{}_{}.html'.format(rank_name, letter)))

        # Parse html pages lising all classes
        rank_sites_list = open(os.path.join(
            self.outdir, f'{rank_name}_list.lst'), 'w')
        link_pattern = re.compile(r'href="(/{}/[^"]+)"'.format(rank_name))
        rank_name_pattern = re.compile(r'class=\"last-child color-{0}\">\"?([^\"]*)\"?</a>'.format(rank_name))
        rank_to_download = []
        self.logger.info(f'Generating a list of all pages that contain {rank_name} information.')
        num_ranks = 0

        for letter in list(string.ascii_uppercase):
            with open(os.path.join(index_dir, '{}_{}.html'.format(rank_name, letter))) as webf:
                for line in webf:
                    if f'class="last-child color-{rank_name}"' in line:
                        line = line.replace("'", '"')
                        result = link_pattern.search(line)
                        if result:
                            rk_name = rank_name_pattern.search(line).group(1)
                            rk_name = clean_html(rk_name)

                            num_ranks += 1
                            print(' - processed {:,} names ({} level).'.format(num_ranks, rank_name), end='\r')
                            rank_sites_list.write('{}\t{}\t{}\n'.format(
                                letter, rk_name, self.base_url + result.group(1)))
                            rank_to_download.append(rk_name)
        rank_sites_list.close()

        # we remove the duplicate names with quotes
        valid_names = []
        for temp_name in rank_to_download:

            if temp_name.startswith('"'):
                if temp_name.replace('"', '').replace('Candidatus ', '') not in rank_to_download:
                    valid_names.append(temp_name)
            elif temp_name.startswith("["):
                if ', no' in temp_name:
                    potential_name = temp_name.split(', no')[0].replace('[', '')
                    if potential_name not in rank_to_download:
                        valid_names.append(temp_name)
            else:
                valid_names.append(temp_name)

        # Download individual classes html page
        self.logger.info(f'Downloading individual {rank_name} HTML pages.')
        failed_html_file = open(os.path.join(
            self.outdir, f'{rank_name}_failed.lst'), 'w')
        make_sure_path_exists(os.path.join(self.outdir, f'all_{rank_name}'))

        num_ranks = 0
        num_already_dl = 0
        with open(os.path.join(self.outdir, f'{rank_name}_list.lst')) as gsl:
            for line in gsl:
                letter, rk_name, rk_url = line.strip().split('\t')
                make_sure_path_exists(os.path.join(
                    self.outdir, f'all_{rank_name}', letter))
                url_name = os.path.basename(rk_url)
                out_file = os.path.join(self.outdir, f'all_{rank_name}', letter, url_name)

                num_already_dl = self.download_rank_name(rk_name, rk_url,
                                                         out_file, valid_names,
                                                         failed_html_file, num_already_dl)

                num_ranks += 1
                print(
                    ' - processed {:,} names ({} level), including {:,} that were previously downloaded\r'.format(
                        num_ranks,
                        rank_name,
                        num_already_dl), end='\r')

        failed_html_file.close()

    def download_rank_name(self, rk_name, rk_url, out_file, valid_names, failed_html_file, num_already_dl):

        if rk_name in valid_names:
            if not os.path.exists(out_file):
                try:
                    urllib.request.urlretrieve(os.path.join(
                        rk_url), out_file)
                except:
                    failed_html_file.write('{}\tfailed_download\n'.format(rk_url))
            else:
                num_already_dl += 1
        else:
            failed_html_file.write('{}\tduplicate_name\n'.format(rk_url))
        return num_already_dl

    def download_subspecies_lpsn_html(self):
        """

        Download all html species pages from LPSN.

        """
        # Download pages listing all subspecies in LPSN
        index_dir = os.path.join(self.outdir, 'subspecies_all_letters')
        if self.skip_taxa_per_letter_dl:
            self.logger.info('Skipping download of subspecies at LPSN and using results in: {}'.format(index_dir))
        else:
            self.logger.info('Beginning download of subspecies from LPSN.')
            search_url = 'https://lpsn.dsmz.de/advanced_search?adv%5Btaxon-name%5D=&adv%5Bcategory%5D=subspecies&adv%5Bnomenclature%5D=&adv%5Bvalid-publ%5D=&adv%5Bcandidatus%5D=&adv%5Bcorrect-name%5D=&adv%5Bauthority%5D=&adv%5Bdeposit%5D=&adv%5Betymology%5D=&adv%5Bgender%5D=&adv%5Bdate-option%5D=&adv%5Bdate%5D=&adv%5Bdate-between%5D=&adv%5Briskgroup%5D=&adv%5Bsubmit%5D=submit-adv#results'
            make_sure_path_exists(os.path.join(index_dir))
            urllib.request.urlretrieve(
                search_url, os.path.join(index_dir, 'search_subspecies.html'))

        # Parse html pages listing all subspecies
        species_sites_list = open(os.path.join(
            self.outdir, 'subspecies_list.lst'), 'w')
        link_pattern = re.compile(r'href="(/subspecies/[^"]+)"')
        subspecies_name_pattern = re.compile(r'class=\"\s?color-subspecies\">(.*)<\/a>')
        subspecies_to_download = []
        self.logger.info('Generating a list of all the pages that contain species information.')

        num_subspecies = 0
        with open(os.path.join(index_dir, 'search_subspecies.html')) as webf:
            for line in webf:
                if 'class=" color-subspecies"' in line:
                    line = line.replace("'", '"')
                    result = link_pattern.search(line)
                    if result:
                        num_subspecies += 1
                        sys.stdout.write(' - processed {:,} subspecies\r'.format(num_subspecies))
                        sys.stdout.flush()
                        subspecies_name = subspecies_name_pattern.search(line).group(1)
                        subspecies_name = clean_html(subspecies_name)
                        subspecies_to_download.append(subspecies_name)
                        species_sites_list.write(
                            '{}\t{}\n'.format(subspecies_name, self.base_url + result.group(1)))
        species_sites_list.close()
        sys.stdout.write('\n')

        # we remove the duplicate names with quotes
        valid_subspecies = []
        for subspname in subspecies_to_download:

            if subspname.startswith('"'):
                if subspname.replace('"', '').replace('Candidatus ', '') not in subspecies_to_download:
                    valid_subspecies.append(subspname)
            else:
                valid_subspecies.append(subspname)

        # Download individual species html page
        self.logger.info('Downloading individual subspecies HTML pages.')
        failed_html_file = open(os.path.join(
            self.outdir, 'subspecies_failed.lst'), 'w')
        make_sure_path_exists(os.path.join(self.outdir, 'all_subspecies'))

        num_subspecies = 0
        num_already_dl = 0
        with open(os.path.join(self.outdir, 'subspecies_list.lst')) as gsl:
            for line in gsl:
                subspecies_name, subspe_url = line.strip().split('\t')
                make_sure_path_exists(os.path.join(
                    self.outdir, 'all_subspecies'))
                subspenname = os.path.basename(subspe_url)
                out_file = os.path.join(self.outdir, 'all_subspecies', subspenname)
                if subspecies_name in valid_subspecies:
                    if not os.path.exists(out_file):
                        try:
                            urllib.request.urlretrieve(os.path.join(
                                subspe_url), out_file)
                        except:
                            failed_html_file.write('{}\n'.format(subspe_url))
                    else:
                        num_already_dl += 1
                else:
                    failed_html_file.write('{}\tduplicate_name\n'.format(subspe_url))
            num_subspecies += 1
            sys.stdout.write(' - processed {:,} subspecies, including {:,} that were previously downloaded\r'.format(
                num_subspecies, num_already_dl))
            sys.stdout.flush()
        failed_html_file.close()
        sys.stdout.write('\n')

    def parse_strains(self, list_strains):
        """
        Parse the HTML line listing all strains and return a list of strains
        """
        results = []
        # replace (now XXXX ) pattern
        now_pattern = re.compile(r'(\w+) \(now (\w+)\) (\d+)')
        # replace (formerly XXXX) pattern
        formerly_pattern = re.compile(r'(\w+) \(formerly (\w+)\) (\d+)')
        for item in list_strains:
            if '(' in item:
                now_res = now_pattern.search(item)
                former_res = formerly_pattern.search(item)
                if now_res:
                    results.append(now_res.group(1) + ' ' + now_res.group(3))
                    results.append(now_res.group(2) + ' ' + now_res.group(3))
                if former_res:
                    results.append(former_res.group(
                        1) + ' ' + former_res.group(3))
                    results.append(former_res.group(
                        2) + ' ' + former_res.group(3))
                if not former_res and not now_res:
                    item = re.sub(r'\([^)]*\)', '', item)
                    results.append(item.strip())

            else:
                results.append(item.strip())
        return results

    def remove_accents(self, data):
        """
        Remove all special characters from reference and species name
        """
        return ''.join(x for x in unicodedata.normalize('NFKD', html.unescape(data)) if x in string.printable)

    def check_format_three_terms_strain(self, strain):
        word_counter = 0
        number_counter = 0
        if strain.count(' ') == 2:
            for item in strain.split(' '):
                if item.isdigit():
                    number_counter += 1
                elif item.isupper() and item.isalpha():
                    word_counter += 1
        if word_counter != 0 and number_counter != 0 and word_counter + number_counter == len(strain.split(' ')):
            return True
        return False

    def parse_generic_html(self, rank_name, headers_order, all_rank, file):
        new_header_pattern = re.compile(r'(<b>([a-zA-Z\d\-_\s]+):</b>)(.*)')
        dict_info = {'Rank': rank_name}
        soup = BeautifulSoup(open(file), "html.parser")
        for it in soup.find_all('p'):
            text_item = str(it).replace('\n', ' ').replace('\t', '').strip()
            new_hder_result = new_header_pattern.search(text_item)
            if new_hder_result:
                if new_hder_result.group(2) not in headers_order:
                    headers_order.append(new_hder_result.group(2))
                all_tags = BeautifulSoup(new_hder_result.group(3), 'html.parser')
                if new_hder_result.group(2).startswith('16'):
                    for tag_to_del in all_tags.find_all('a'):
                        tag_to_del.decompose()
                dict_info[new_hder_result.group(2)] = all_tags.text.strip()
        for it in soup.find_all('ul', class_='notes-list'):
            table_header = it.find('span').find('b').text.replace(':', '')
            if table_header not in headers_order:
                headers_order.append(table_header)

            if it.find('tbody'):
                for th in it.select('thead tr'):
                    dict_info.setdefault(table_header, []).append(
                        ';'.join([th.text for th in th.select('th') if th.text != ' ']))
                for tr in it.select('tbody tr'):
                    dict_info.setdefault(table_header, []).append(
                        ';'.join([td.text.strip() for td in tr.select('td') if td.text != ' ']))
                dict_info[table_header] = '|'.join(dict_info.get(table_header))
            else:
                for ul in it.find_all('ul'):
                    dict_info.setdefault(table_header, []).append(
                        re.sub(' +', ' ', ul.find('li').text.replace('\n', ' ').replace('\t', '').strip()))
                dict_info[table_header] = '|'.join(dict_info.get(table_header))

        dict_info.pop("Linking", None)

        all_rank.append(dict_info)
        return headers_order, all_rank

    def parse_rank_html(self, rank_name, input_dir, headers_order, all_rank):

        if rank_name != 'subspecies':
            for letter in list(string.ascii_uppercase):
                print(f'letter: {letter}', end='\r')
                for file in glob.glob(os.path.join(os.path.join(input_dir, f'all_{rank_name}'), letter, "*")):
                    headers_order, all_rank = self.parse_generic_html(rank_name, headers_order, all_rank, file)
        else:
            for file in glob.glob(os.path.join(input_dir, 'all_subspecies', "*")):
                headers_order, all_rank = self.parse_generic_html(rank_name, headers_order, all_rank, file)

        return headers_order, all_rank

    def parse_all_ranks_tsv(self, raw_allranks_file):
        parsed_file = open(os.path.join(os.path.dirname(raw_allranks_file), 'full_parsing_parsed.tsv'), 'w')

        with open(raw_allranks_file) as fp:
            header_line = fp.readline()
            headers = header_line.strip().split('\t')

            rank_index = headers.index('Rank')
            name_index = headers.index('Name')
            originalpub_index = headers.index('Original publication')
            correctname_index = headers.index('Correct name')
            parenttaxon_index = headers.index('Parent taxon')
            childtaxa_index = headers.index('Child taxa')
            synonyms_index = headers.index('Synonyms')
            type_class_index = headers.index('Type class')
            type_order_index = headers.index('Type order')
            type_genus_index = headers.index('Type genus')
            type_species_index = headers.index('Type species')

            headers.insert(originalpub_index, 'Priority')

            nomenclature_index = headers.index('Nomenclatural status')

            name_pattern = re.compile(r'\[([^,]*),.*\]')

            all_infos = []
            all_infos_dict = {}
            for line in fp:
                infos = line.strip().split('\t')
                # rank_of_interest
                roi = infos[rank_index]

                # Parse "Correct Name" so it is just the taxon name without the publication
                if roi in ['species', 'subspecies'] and len(infos[correctname_index].split(" ")) >= 2:
                    if (len(infos[correctname_index].split(" ")) > 3 and
                            infos[correctname_index].split(" ")[2] == 'subsp.'):
                        infos[correctname_index] = " ".join(infos[correctname_index].split(' ')[0:4]).replace('"', '')
                    else:
                        infos[correctname_index] = " ".join(infos[correctname_index].split(' ')[0:2]).replace('"', '')
                elif roi != 'species' and len(infos[correctname_index].split(" ")) >= 1:
                    infos[correctname_index] = infos[correctname_index].split(' ')[0]

                # remove citation information from "Type <rank>" fields
                for type_index in [type_class_index, type_order_index, type_genus_index]:
                    if len(infos[type_index].split(' ')) >= 1 and infos[type_index] != 'n/a':
                        temp_name = infos[type_index].replace('"', '').replace("[", "").replace("]", "").replace(
                            ',', '')
                        if temp_name.startswith("Candidatus "):
                            temp_name = temp_name.replace("Candidatus ", "")
                        temp_name = temp_name.split(" ")[0]
                        infos[type_index] = temp_name
                # sorting "type species"
                if len(infos[type_species_index].split(' ')) >= 2:
                    temp_name = infos[type_species_index].replace('"', '').replace("[", "").replace("]", "").replace(
                        ',', '')
                    if temp_name.startswith("Candidatus "):
                        temp_name = temp_name.replace("Candidatus ", "")
                    temp_name = temp_name.split(" ")[0] + " " + temp_name.split(" ")[1]
                    infos[type_species_index] = temp_name

                # Parse the "Parent taxon" so it is just the taxon name without the publication
                if roi == 'subspecies' and len(infos[parenttaxon_index].split(" ")) >= 2:
                    parent = infos[parenttaxon_index].replace('"', '')
                    if parent.startswith('Candidatus'):
                        infos[parenttaxon_index] = " ".join(parent.split(' ')[1:3])
                    else:
                        infos[parenttaxon_index] = " ".join(parent.split(' ')[0:2])
                elif len(infos[parenttaxon_index].split(" ")) >= 1:
                    parent = infos[parenttaxon_index].replace('"', '')
                    if parent.startswith('Candidatus'):
                        infos[parenttaxon_index] = parent.split(' ')[1]
                    else:
                        infos[parenttaxon_index] = parent.split(' ')[0]
                infos[parenttaxon_index] = re.sub(r'\"|\[|\]|,','',infos[parenttaxon_index])

                # Clean up the "Child taxa" to a comma separated list(e.g.: for Acidobacteria this would be just
                # Acidobacteriia, Blastocatellia, ..., Vicinamibacteria)
                child_taxa = []
                if roi in ['phylum', 'class', 'order', 'family']:
                    temp_list_step_one = [x.split(';')[0] for x in infos[childtaxa_index].split('|')[1:]]
                    for temp_item in temp_list_step_one:
                        temp_item = temp_item.replace('"', '').replace("[", "").replace("]", "").replace(',', '')
                        if temp_item.startswith("Candidatus "):
                            temp_item = temp_item.replace("Candidatus ", "")
                        temp_item = temp_item.split(' ')[0]
                        child_taxa.append(temp_item)

                if roi in ['genus', 'species']:
                    for row in infos[childtaxa_index].split('|')[1:]:
                        raw_child = row.split(';')[0]
                        temp_item = raw_child.replace('"', '').replace("[", "").replace("]", "").replace(',', '')
                        if temp_item.startswith("Candidatus "):
                            temp_item = temp_item.replace("Candidatus ", "")
                        if len(temp_item.split(' ')) > 3 and temp_item.split(' ')[2] == 'subsp.':
                            child_taxa.append(" ".join(temp_item.split(' ')[0:4]))
                        elif len(temp_item.split(' ')) > 1:
                            child_taxa.append(" ".join(temp_item.split(' ')[0:2]).replace('"', ''))
                infos[childtaxa_index] = ','.join(set(child_taxa))

                # Clean up the "Synonym the same way"
                syns = []
                if roi in ['phylum', 'class', 'order', 'family', 'genus']:
                    temp_list_step_one = [x.split(';')[0] for x in infos[synonyms_index].split('|')[1:]]
                    for temp_item in temp_list_step_one:
                        temp_item = temp_item.replace('"', '').replace("[", "").replace("]", "")
                        if temp_item.startswith("Candidatus "):
                            temp_item = temp_item.replace("Candidatus ", "")
                        temp_item = temp_item.split(' ')[0]
                        syns.append(temp_item)
                if roi in ['species', 'subspecies']:
                    for row in infos[synonyms_index].split('|')[1:]:
                        raw_syn = row.split(';')[0]
                        temp_item = raw_syn.replace('"', '').replace("[", "").replace("]", "")
                        if temp_item.startswith("Candidatus "):
                            temp_item = temp_item.replace("Candidatus ", "")
                        if len(temp_item.split(' ')) > 3 and temp_item.split(' ')[2] == 'subsp.':
                            syns.append(" ".join(temp_item.split(' ')[0:4]))
                        elif len(temp_item.split(' ')) > 1:
                            syns.append(" ".join(temp_item.split(' ')[0:2]).replace('"', ''))
                infos[synonyms_index] = ','.join(set(syns))

                if 'Candidatus' in infos[name_index]:
                    infos[name_index] = infos[name_index].replace('Candidatus ', '').strip()

                infos[name_index] = infos[name_index].replace('"', '')
                if infos[rank_index] == 'species' and len(infos[name_index].split(" ")) < 2:
                    continue

                if not infos[name_index].startswith('['):
                    if infos[rank_index] == 'species':
                        split_name = infos[name_index].split(" ", 2)
                        infos[name_index] = split_name[0] + " " + split_name[1]
                        infos.insert(originalpub_index, split_name[2])
                    elif infos[rank_index] == 'subspecies':
                        split_name = infos[name_index].split(" ", 4)
                        infos[name_index] = " ".join(split_name[0:4])
                        infos.insert(originalpub_index, split_name[4])
                    else:
                        split_name = infos[name_index].split(" ", 1)
                        infos[name_index] = split_name[0]
                        infos.insert(originalpub_index, split_name[1])
                else:
                    name_result = name_pattern.search(infos[name_index])
                    if name_result:
                        infos[name_index] = name_result.group(1)
                    infos.insert(originalpub_index, 'n/a')

                if infos[rank_index] + '_' + infos[name_index] in all_infos_dict:
                    previous_infos = all_infos_dict.get(infos[rank_index] + '_' + infos[name_index])
                    if self.check_illegimate(previous_infos[nomenclature_index]):
                        all_infos_dict[infos[rank_index] + '_' + infos[name_index]] = infos
                else:
                    all_infos_dict[infos[rank_index] + '_' + infos[name_index]] = infos

                all_infos.append(infos)

        df = pd.DataFrame.from_records(list(all_infos_dict.values()), columns=headers)

        reordering_header = ['Rank', 'Name', 'Category']
        for type_rank in ['phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']:
            if f"Type {type_rank}" in headers:
                reordering_header.append(f"Type {type_rank}")
                df[f"Type {type_rank}"] = df[f"Type {type_rank}"].str.replace('"', '')

        remaining_columns = [x for x in headers if x not in reordering_header]
        reordering_header.extend(remaining_columns)

        df = df[reordering_header]
        df.to_csv(parsed_file, sep="\t", index=False)
        parsed_file.close()
        return parsed_file

    def add_lpsn_metadata(self, hostname, user, password, db, lpsn_file):

        engine_current = create_engine(f'postgresql://{user}:{password}@{hostname}:5432/{db}',
                                       convert_unicode=True,
                                       pool_size=5,
                                       max_overflow=20,
                                       pool_recycle=3600)
        df = pd.read_csv(lpsn_file, sep='\t')
        df.columns = [x.lower() for x in df.columns]
        df.columns = df.columns.str.replace(' ', '_')
        df.to_sql('lpsn_metadata', engine_current, if_exists='replace')

    def parse_subspecies_html(self, output_file, input_dir, full_list_type_species):

        type_proposed_pattern = re.compile(
            '<b>Proposed as:</b>(.*)</p>')

        for file in glob.glob(os.path.join(input_dir, 'all_subspecies', "*")):

            short_name = ''
            raw_list_strain = []
            strain_section = False
            correct_section = False
            skip_species = False

            raw_subspe_name = os.path.basename(file)
            subspe_nameparts = raw_subspe_name.split('-')
            subspe_name = ' '.join([subspe_nameparts[0].capitalize(
            ), subspe_nameparts[1], 'subsp.', subspe_nameparts[2]])
            if subspe_nameparts[1] == subspe_nameparts[2]:
                short_name = subspe_nameparts[0] + " " + subspe_nameparts[1]
            name_section = True

            with open(file) as f:
                for line in f:
                    line = line.replace("'", '"')
                    if 'glossary#citation-of-names-and-authors' in line:
                        subspecies_reference = clean_html(
                            line).split(subspe_name, 1)[1]
                    elif "glossary#abbreviations-in-proposals" in line:
                        type_proposed_section = True
                        type_proposed_pattern_results = type_proposed_pattern.search(
                            line)
                        if type_proposed_pattern_results:
                            subspecies_proposed_type = type_proposed_pattern_results.group(
                                1)
                    elif ('<b>Type strains:</b>' in line or '<b>Type strain:</b>' in line) and name_section:
                        strain_section = True
                    elif '<p' in line and strain_section and name_section:
                        raw_list_strain = list(
                            filter(lambda a: a.strip() != "no culture isolated", raw_list_strain))
                        raw_list_strain = list(
                            filter(lambda a: a.strip() != "no culture available", raw_list_strain))
                        raw_list_strain = list(
                            filter(lambda a: a.strip() != "no pure culture", raw_list_strain))
                        raw_list_strain = self.parse_strains(raw_list_strain)

                        raw_list_strain = [
                            x for x in raw_list_strain if (x != 'n/a' and check_format_strain(canonical_strain_id(x)))]
                        break
                        # print(raw_list_strain)

                    elif strain_section and name_section:
                        if line.strip() != '':
                            raw_list_strain.extend(
                                [x.replace('strain', '').replace('Strain', '').strip() for x in
                                 clean_html(line).split(';')])

                    elif '<p class="corr-name">' in line:
                        correct_section = True
                    elif correct_section:
                        correct_subspecies_pattern = re.compile(
                            r'color-subspecies">\"?<I>([a-zA-Z]+)</I> <I>([a-zA-Z]+)</I> subsp. <I>([a-zA-Z]+)</I>')
                        correct_subspecies_results = correct_subspecies_pattern.search(line)
                        if correct_subspecies_results:
                            correct_subspecies_name = correct_subspecies_results.group(
                                1) + " " + correct_subspecies_results.group(
                                2) + " subsp. " + correct_subspecies_results.group(3)
                            if correct_subspecies_name == subspe_name:
                                skip_species = True
                                break
                    elif correct_section and 'class="helper"' in line:
                        correct_section = False

                ref_type_combined = ','.join(
                    list(filter(None, [subspecies_reference, subspecies_proposed_type]))).strip()

                if short_name == '':
                    if not skip_species:
                        output_file.write(
                            "species\t{}\t{}\t{}\t{}\n".format(
                                short_name != '' and short_name in full_list_type_species,
                                subspe_name, ref_type_combined, "=".join(raw_list_strain)))

    def check_illegimate(self, status_raw):
        status = status_raw.strip().replace('"', '')
        status_tokens = [t.strip() for t in status.split(';')]
        status_tokens = [tt.strip() for t in status_tokens for tt in t.split(',')]

        if 'illegitimate name' in status_tokens:
            return True
        return False

    def parse_gss(self, gss_file):
        gss_species = {}
        record_mapping = {}
        type_species_genus = {}
        with open(gss_file) as lpsng:
            headers = lpsng.readline().strip().split(',')
            genus_name_idx = headers.index('genus_name')
            sp_epithet_idx = headers.index('sp_epithet')
            subsp_epithet_idx = headers.index('subsp_epithet')
            nomenclatural_type_idx = headers.index('nomenclatural_type')
            record_no_idx = headers.index('record_no')
            authors_idx = headers.index('authors')
            status_idx = headers.index('status')

        with open(gss_file, "r") as csvfile:
            csvreader = csv.reader(csvfile)
            next(csvreader)
            for row in csvreader:
                if row[sp_epithet_idx] != '':
                    if row[subsp_epithet_idx] != '':
                        temp_name = row[genus_name_idx] + ' ' + row[sp_epithet_idx] + ' subsp. ' + row[
                            subsp_epithet_idx]
                        status_list = row[status_idx]
                        is_illegitimate = self.check_illegimate(status_list)
                        if (temp_name in gss_species and is_illegitimate is False
                                and gss_species[temp_name]['illegitimate'] is True):
                            gss_species[temp_name] = {'type': 'species', 'authors': row[authors_idx],
                                                      'illegitimate': is_illegitimate,
                                                      'strains': row[nomenclatural_type_idx],
                                                      'record_no': row[record_no_idx]}
                            record_mapping[row[record_no_idx]] = temp_name
                        elif temp_name not in gss_species:
                            gss_species[temp_name] = {'type': 'species', 'authors': row[authors_idx],
                                                      'illegitimate': is_illegitimate,
                                                      'strains': row[nomenclatural_type_idx],
                                                      'record_no': row[record_no_idx]}
                            record_mapping[row[record_no_idx]] = temp_name

                    else:
                        temp_name = row[genus_name_idx] + ' ' + row[sp_epithet_idx]
                        status_list = row[status_idx]
                        is_illegitimate = self.check_illegimate(status_list)
                        if (temp_name in gss_species and is_illegitimate is False
                                and gss_species[temp_name]['illegitimate'] is True):
                            gss_species[temp_name] = {'type': 'species', 'authors': row[authors_idx],
                                                      'illegitimate': is_illegitimate,
                                                      'strains': row[nomenclatural_type_idx],
                                                      'record_no': row[record_no_idx]}
                            record_mapping[row[record_no_idx]] = temp_name
                        elif temp_name not in gss_species:
                            gss_species[temp_name] = {'type': 'species', 'authors': row[authors_idx],
                                                      'illegitimate': is_illegitimate,
                                                      'strains': row[nomenclatural_type_idx],
                                                      'record_no': row[record_no_idx]}

                else:
                    gss_species[row[genus_name_idx]] = {'type': 'species', 'authors': row[authors_idx],
                                                        'type_species': row[record_no_idx]}
                    record_mapping[row[record_no_idx]] = row[genus_name_idx]
                    type_species_genus[row[nomenclatural_type_idx]] = row[genus_name_idx]

        return gss_species, record_mapping, type_species_genus

    def parse_html_ts(self, output_all_ranks):
        results_tsog, results_tgof, results_tgoo = {}, {}, {}
        with open(output_all_ranks) as lsf:
            line = lsf.readline()
            headers = line.strip().split('\t')
            rank_index = headers.index('Rank')
            name_index = headers.index('Name')
            type_species_index = headers.index('Type species')
            type_genus_index = headers.index('Type genus')

            for line in lsf:
                infos = line.strip().split('\t')
                if infos[rank_index] == 'genus' and infos[type_species_index] != 'n/a':
                    results_tsog[infos[type_species_index]] = infos[name_index]
                if infos[rank_index] == 'family' and infos[type_genus_index] != 'n/a':
                    results_tgof[infos[type_genus_index]] = infos[name_index]
                if infos[rank_index] == 'order' and infos[type_genus_index] != 'n/a':
                    results_tgoo[infos[type_genus_index]] = infos[name_index]

        return results_tsog, results_tgof, results_tgoo

    def summarise_parsing(self, output_all_ranks, gss_file):
        # parse gss_file
        gss_dict, record_mapping, type_species_genus = self.parse_gss(gss_file)

        html_type_spe_of_gen, html_type_gen_of_fam, html_type_gen_of_ord = self.parse_html_ts(output_all_ranks)

        """Create metadata by parsing assembly stats files."""

        # identify type genera, species, and strains according to LPSN
        fout_type_genera = open(os.path.join(
            os.path.dirname(output_all_ranks), 'lpsn_genera.tsv'), 'w')
        fout_type_species = open(os.path.join(
            os.path.dirname(output_all_ranks), 'lpsn_species.tsv'), 'w')
        fout_type_strains = open(os.path.join(
            os.path.dirname(output_all_ranks), 'lpsn_strains.tsv'), 'w')

        fout_type_genera.write(
            'lpsn_genus\tlpsn_type_genus_of_family\tlpsn_type_genus_of_order\tlpsn_genus_authority\n')
        fout_type_species.write(
            'lpsn_species\tlpsn_type_species\tlpsn_species_authority\tsource\n')
        fout_type_strains.write('lpsn_strain\tco-identical strain IDs\n')

        list_processed_strains = []
        processed_genus = []
        processed_species = []

        strains = set()

        # Parse the lpsn summary file
        with open(output_all_ranks) as lsf:
            line = lsf.readline()
            headers = line.strip().split('\t')
            rank_index = headers.index('Rank')
            name_index = headers.index('Name')

            type_strains_index = headers.index('Type strain')
            priority_index = headers.index('Priority')

            for line in lsf:
                line_split = line.rstrip('\n').split('\t')

                if line_split[rank_index] == 'genus':
                    genus_name = line_split[name_index]
                    genus = 'g__' + genus_name
                    if genus_name in gss_dict:
                        desc = gss_dict.get(genus_name).get('authors')
                    else:
                        desc = line_split[priority_index]

                    family = ''
                    order = ''
                    if genus_name in html_type_gen_of_fam:
                        family = 'f__' + html_type_gen_of_fam.get(genus_name)
                    if genus_name in html_type_gen_of_ord:
                        order = 'o__' + html_type_gen_of_ord.get(genus_name)

                    if (genus, family, order, desc) not in processed_genus:
                        fout_type_genera.write(f'{genus}\t{family}\t'
                                               f'{order}\t{desc}\n')
                        processed_genus.append((genus, family, order, desc))

                if line_split[rank_index] in ['species', 'subspecies']:
                    processed_strains = []
                    spe_name = line_split[name_index]
                    species = 's__' + spe_name
                    if spe_name in gss_dict:
                        source = 'GSS'
                        desc = gss_dict.get(spe_name).get('authors')
                        genus = ''

                        if gss_dict.get(spe_name).get('record_no') in type_species_genus:
                            genus = type_species_genus.get(gss_dict.get(spe_name).get('record_no'))

                        strains = gss_dict.get(spe_name).get('strains').split(';')

                    else:
                        source = 'HTML'
                        desc = line_split[priority_index]
                        genus = html_type_spe_of_gen.get(spe_name, '')
                        strains = line_split[type_strains_index].split(';')

                    if (species, genus, desc) not in processed_species:
                        fout_type_species.write(f'{species}\t{genus}\t{desc}\t{source}\n')
                        processed_species.append((species, genus, desc))

                    # Normalise the strains
                    for i, strain in enumerate(strains):
                        if strain != 'n/a':
                            processed_strains.append(canonical_strain_id(strain))
                    processed_strain_string = '{0}\t{1}'.format(
                        spe_name, "=".join(processed_strains))
                    if processed_strain_string not in list_processed_strains:
                        fout_type_strains.write(
                            '{0}\n'.format(processed_strain_string))
                        list_processed_strains.append(processed_strain_string)

        fout_type_genera.close()
        fout_type_species.close()
        fout_type_strains.close()

    def parse_html(self, input_dir, gss_file):
        """
        Parse the html file of each genus.
        Store the type, the name, the reference, the strains for each species.
        """
        make_sure_path_exists(os.path.join(self.outdir, 'all_ranks'))

        self.logger.info('Parsing all pages.')
        headers_order = ['Rank']
        all_rank = []
        for rk in ['phylum', 'class', 'order', 'family', 'genus', 'species', 'subspecies']:
            self.logger.info(f'Parsing {rk}.')
            headers_order, all_rank = self.parse_rank_html(rk, input_dir, headers_order, all_rank)

        output_all_ranks = open(os.path.join(self.outdir, 'all_ranks', 'full_parsing_raw.tsv'), 'w')
        output_all_ranks.write('\t'.join(headers_order) + '\n')
        for item in all_rank:
            output_all_ranks.write(
                '\t'.join([item.get(potential_header, 'n/a') for potential_header in headers_order]) + '\n')

        output_all_ranks.close()
        parsed_file = self.parse_all_ranks_tsv(output_all_ranks.name)

        self.summarise_parsing(parsed_file.name, gss_file)
