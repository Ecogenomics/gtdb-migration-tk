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
import argparse
import urllib.request
import re
import unicodedata
import string
import html
from biolib.common import check_file_exists, make_sure_path_exists, check_dir_exists


class LPSN(object):
    def __init__(self, lpsn_output_dir):
        """Initialization."""
        self.outdir = lpsn_output_dir
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

    def full_lpsn_wf(self):
        # self.download_lpsn_html()
        self.parse_html(os.path.join(self.outdir, 'genus_html'))

    def summarise_parsing(self, lpsn_scrape_file):
        """Create metadata by parsing assembly stats files."""

        # identify type genera, species, and strains according to LPSN
        fout_type_genera = open(os.path.join(
            self.outdir, 'lpsn_genera.tsv'), 'w')
        fout_type_species = open(os.path.join(
            self.outdir, 'lpsn_species.tsv'), 'w')
        fout_type_strains = open(os.path.join(
            self.outdir, 'lpsn_strains.tsv'), 'w')

        fout_type_genera.write(
            'lpsn_genus\tlpsn_type_genus\tlpsn_genus_authority\n')
        fout_type_species.write(
            'lpsn_species\tlpsn_type_species\tlpsn_species_authority\n')
        fout_type_strains.write('lpsn_strain\n')

        list_processed_strains = []
        processed_genus = []
        processed_species = []

        strains = set()

        # Parse the lpsn summary file
        with open(lpsn_scrape_file) as lsf:
            for line in lsf:
                line_split = line.rstrip('\n').split('\t')

                if line_split[0] == 'genus':
                    genus = 'g__' + line_split[2]
                    desc = line_split[3].strip()

                    family = ''
                    if line_split[1] == 'True':
                        family = desc[desc.find(
                            'family ?') + len('family ?'):].strip()
                        family = 'f__' + family.split()[0]

                    desc = desc.replace(' ?', '')

                    if (genus, family, desc) not in processed_genus:
                        fout_type_genera.write('%s\t%s\t%s\n' %
                                               (genus, family, desc))
                        processed_genus.append((genus, family, desc))
                elif line_split[0] == 'species':
                    species = 's__' + line_split[2]
                    desc = line_split[3].strip()

                    genus = ''
                    if line_split[1] == 'True':
                        genus = 'g__' + line_split[2].split()[0]

                    if (species, genus, desc) not in processed_species:
                        cleanr = re.compile('<.*?>')
                        species = re.sub(cleanr, '', species)
                        genus = re.sub(cleanr, '', genus)
                        fout_type_species.write(
                            '%s\t%s\t%s\n' % (species, genus, desc))
                        processed_species.append((species, genus, desc))
                    processed_strains = []
                    strains = line_split[4].split("=")
                    pattern = re.compile('[\W_]+')
                    # Normalise the strains
                    for i, strain in enumerate(strains):
                        if i == 0 and strain.startswith('strain '):
                            strain = strain.replace("strain ", "", 1)
                        strain = re.sub(r'\(.+\)', ' ', strain)
                        strain = ' '.join(strain.split())
                        b = pattern.sub('', strain).upper()
                        processed_strains.append(b)
                    processed_neotypes = []
                    processed_strain_string = '{0}\t{1}'.format(
                        line_split[2], "=".join(processed_strains))
                    if processed_strain_string not in list_processed_strains:
                        fout_type_strains.write(
                            '{0}\n'.format(processed_strain_string))
                        list_processed_strains.append(processed_strain_string)

        fout_type_genera.close()
        fout_type_species.close()
        fout_type_strains.close

    def download_family_lpsn_html(self):
        '''

        Download all html species pages from LPSN.

        '''
        # Download pages listing all genus in LPSN
        print('Beginning file download lpsn ...')
        base_url = 'https://lpsn.dsmz.de/'
        make_sure_path_exists(os.path.join(self.outdir, 'family_per_letter'))
        for letter in list(string.ascii_uppercase):
            url = base_url + 'family?page=' + letter
            urllib.request.urlretrieve(
                url, os.path.join(self.outdir, 'family_per_letter', 'family_{}.html'.format(letter)))

        # Parse html pages lising all genus
        family_sites_list = open(os.path.join(
            self.outdir, 'family_list.lst'), 'w')
        link_pattern = re.compile(r'href="(/family/[^"]+)"')
        print('Get a list of all the pages that contain genus information ...')
        for letter in list(string.ascii_uppercase):
            with open(os.path.join(self.outdir, 'family_per_letter', 'family_{}.html'.format(letter))) as webf:
                for line in webf:
                    if 'class="last-child color-family"' in line:
                        line = line.replace("'", '"')
                        result = link_pattern.search(line)
                        if result:
                            print(line)
                            family_sites_list.write('{}\t{}\n'.format(
                                letter, base_url + result.group(1)))
        family_sites_list.close()

        # Download individual species html page
        failed_html_file = open(os.path.join(
            self.outdir, 'families_failed.lst'), 'w')
        make_sure_path_exists(os.path.join(self.outdir, 'all_families'))
        with open(os.path.join(self.outdir, 'family_list.lst')) as gsl:
            for line in gsl:
                letter, fam_url = line.strip().split('\t')
                make_sure_path_exists(os.path.join(
                    self.outdir, 'all_families', letter))
                famname = os.path.basename(fam_url)
                try:
                    urllib.request.urlretrieve(os.path.join(
                        fam_url), os.path.join(self.outdir, 'all_families', letter, famname))
                except:
                    failed_html_file.write('{}\n'.format(fam_url))
        failed_html_file.close()

    def download_genus_lpsn_html(self):
        '''

        Download all html species pages from LPSN.

        '''
        # Download pages listing all genus in LPSN
        print('Beginning file download lpsn ...')
        base_url = 'https://lpsn.dsmz.de/'
        make_sure_path_exists(os.path.join(self.outdir, 'genus_per_letter'))
        for letter in list(string.ascii_uppercase):
            url = base_url + 'genus?page=' + letter
            urllib.request.urlretrieve(
                url, os.path.join(self.outdir, 'genus_per_letter', 'genus_{}.html'.format(letter)))

        # Parse html pages lising all genus
        genus_sites_list = open(os.path.join(
            self.outdir, 'genus_list.lst'), 'w')
        link_pattern = re.compile(r'href="(/genus/[^"]+)"')
        print('Get a list of all the pages that contain genus information ...')
        for letter in list(string.ascii_uppercase):
            with open(os.path.join(self.outdir, 'genus_per_letter', 'genus_{}.html'.format(letter))) as webf:
                for line in webf:
                    if 'class="last-child color-genus"' in line:
                        line = line.replace("'", '"')
                        result = link_pattern.search(line)
                        if result:
                            print(line)
                            genus_sites_list.write('{}\t{}\n'.format(
                                letter, base_url + result.group(1)))
        genus_sites_list.close()

        # Download individual species html page
        failed_html_file = open(os.path.join(
            self.outdir, 'genera_failed.lst'), 'w')
        make_sure_path_exists(os.path.join(self.outdir, 'all_genera'))
        with open(os.path.join(self.outdir, 'genus_list.lst')) as gsl:
            for line in gsl:
                letter, gen_url = line.strip().split('\t')
                make_sure_path_exists(os.path.join(
                    self.outdir, 'all_genera', letter))
                genname = os.path.basename(gen_url)
                try:
                    urllib.request.urlretrieve(os.path.join(
                        gen_url), os.path.join(self.outdir, 'all_genera', letter, genname))
                except:
                    failed_html_file.write('{}\n'.format(gen_url))
        failed_html_file.close()

    def download_species_lpsn_html(self):
        '''

        Download all html species pages from LPSN.

        '''
        # Download pages listing all genus in LPSN
        print('Beginning file download lpsn ...')
        base_url = 'https://lpsn.dsmz.de/'
        make_sure_path_exists(os.path.join(self.outdir, 'species_per_letter'))
        for letter in list(string.ascii_uppercase):
            url = base_url + 'species?page=' + letter
            urllib.request.urlretrieve(
                url, os.path.join(self.outdir, 'species_per_letter', 'species_{}.html'.format(letter)))

        # Parse html pages lising all genus
        species_sites_list = open(os.path.join(
            self.outdir, 'species_list.lst'), 'w')
        link_pattern = re.compile(r'href="(/species/[^"]+)"')
        print('Get a list of all the pages that contain species information ...')
        for letter in list(string.ascii_uppercase):
            with open(os.path.join(self.outdir, 'species_per_letter', 'species_{}.html'.format(letter))) as webf:
                for line in webf:
                    if 'class="last-child color-species"' in line:
                        line = line.replace("'", '"')
                        result = link_pattern.search(line)
                        if result:
                            print(line)
                            species_sites_list.write('{}\t{}\n'.format(
                                letter, base_url + result.group(1)))
        species_sites_list.close()

        # Download individual species html page
        failed_html_file = open(os.path.join(
            self.outdir, 'species_failed.lst'), 'w')
        make_sure_path_exists(os.path.join(self.outdir, 'all_species'))
        with open(os.path.join(self.outdir, 'species_list.lst')) as gsl:
            for line in gsl:
                letter, spe_url = line.strip().split('\t')
                make_sure_path_exists(os.path.join(
                    self.outdir, 'all_species', letter))
                spenname = os.path.basename(spe_url)
                try:
                    urllib.request.urlretrieve(os.path.join(
                        spe_url), os.path.join(self.outdir, 'all_species', letter, spenname))
                except:
                    failed_html_file.write('{}\n'.format(spe_url))
        failed_html_file.close()

    def download_subspecies_lpsn_html(self):
        '''

        Download all html species pages from LPSN.

        '''
        # Download pages listing all genus in LPSN
        print('Beginning file download lpsn ...')
        search_url = 'https://lpsn.dsmz.de/advanced_search?adv%5Btaxon-name%5D=&adv%5Bcategory%5D=subspecies&adv%5Bnomenclature%5D=&adv%5Bvalid-publ%5D=&adv%5Bcandidatus%5D=&adv%5Bcorrect-name%5D=&adv%5Bauthority%5D=&adv%5Bdeposit%5D=&adv%5Betymology%5D=&adv%5Bgender%5D=&adv%5Bdate-option%5D=&adv%5Bdate%5D=&adv%5Bdate-between%5D=&adv%5Briskgroup%5D=&adv%5Bsubmit%5D=submit-adv#results'
        base_url = 'https://lpsn.dsmz.de/'
        make_sure_path_exists(os.path.join(
            self.outdir, 'subspecies_all_letters'))
        urllib.request.urlretrieve(
            search_url, os.path.join(self.outdir, 'subspecies_all_letters', 'search_subspecies.html'))

        # Parse html pages lising all genus
        species_sites_list = open(os.path.join(
            self.outdir, 'subspecies_list.lst'), 'w')
        link_pattern = re.compile(r'href="(/subspecies/[^"]+)"')
        print('Get a list of all the pages that contain species information ...')
        with open(os.path.join(self.outdir, 'subspecies_all_letters', 'search_subspecies.html')) as webf:
            for line in webf:
                if 'class=" color-subspecies"' in line:
                    line = line.replace("'", '"')
                    result = link_pattern.search(line)
                    if result:
                        print(line)
                        species_sites_list.write(
                            '{}\n'.format(base_url + result.group(1)))
        species_sites_list.close()

        # Download individual species html page
        failed_html_file = open(os.path.join(
            self.outdir, 'subspecies_failed.lst'), 'w')
        make_sure_path_exists(os.path.join(self.outdir, 'all_subspecies'))
        with open(os.path.join(self.outdir, 'subspecies_list.lst')) as gsl:
            for line in gsl:
                subspe_url = line.strip()
                make_sure_path_exists(os.path.join(
                    self.outdir, 'all_subspecies'))
                subspenname = os.path.basename(subspe_url)
                try:
                    urllib.request.urlretrieve(os.path.join(
                        subspe_url), os.path.join(self.outdir, 'all_subspecies', subspenname))
                except:
                    failed_html_file.write('{}\n'.format(subspe_url))
        failed_html_file.close()

    def cleanhtml(self, raw_html):
        '''
        Remove all HTML tags from HTML line
        '''
        cleanr = re.compile('<.*?>')
        cleantext = re.sub(cleanr, '', raw_html)
        return cleantext.strip()

    def clean_parenthesis(self, raw_line):
        '''
        Remove parenthesis content from HTML line
        '''
        #result = re.sub("[\(\[].*?[\)\]]", "", raw_line)
        result = re.sub("\(see[^\)]+\)", "", raw_line)
        return result

    def parse_strains(self, list_strains):
        '''
        Parse the HTML line listing all strains and return a list of strains
        '''
        results = []
        # replace (now XXXX ) pattern
        now_pattern = re.compile('(\w+) \(now (\w+)\) (\d+)')
        # replace (formerly XXXX) pattern
        formerly_pattern = re.compile('(\w+) \(formerly (\w+)\) (\d+)')
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
        '''
        Remove all special characters from reference and species name
        '''
        return ''.join(x for x in unicodedata.normalize('NFKD', html.unescape(data)) if x in string.printable)

    def check_format_strain(self, strain):
        if not any(char.isdigit() for char in strain):
            return False
        # self.check_format_three_terms_strain(strain)
        if all(c.isdigit() or c.isupper() for c in strain):
            return True

        special_characters = ['-', '.', ' ']
        processed_strain = str(strain)
        for spechar in special_characters:
            processed_strain = processed_strain.replace(spechar, '')
        print(processed_strain, strain)
        if all(c.isdigit() or c.isupper() for c in processed_strain):
            return True
        if strain.count(' ') == 0 and all(c.isdigit() or c.isupper() or c.lower() for c in processed_strain):
            return True
        return False

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

    def parse_family_html(self, output_file, input_dir):
        # Pattern for family
        type_genus_pattern = re.compile('color-genus">\"?<I>([a-zA-Z]+)</I>')

        full_list_type_genus = []
        for letter in list(string.ascii_uppercase):
            for file in glob.glob(os.path.join(os.path.join(input_dir, 'all_families'), letter, "*")):
                family_name = ''
                family_reference = ''
                name_section = False
                type_genus_section = False
                type_genus_list = []
                with open(file) as f:
                    for line in f:
                        line = line.replace("'", '"')

                        if '<b>Type genus:</b>' in line:
                            type_genus_section = True
                        elif '<p' in line and type_genus_section:
                            break
                        elif type_genus_section:
                            type_genus_pattern_results = type_genus_pattern.search(
                                line)
                            if type_genus_pattern_results:
                                type_genus_name = type_genus_pattern_results.group(
                                    1)
                                type_genus_list.append(type_genus_name)

                full_list_type_genus.extend(type_genus_list)
        return full_list_type_genus

    def parse_genus_html(self, output_file, input_dir, full_list_type_genus):
        # Patterns for genus
        genus_pattern = re.compile(
            r'<b>Name:</b> \"?\s?<I>([a-zA-Z]+)</I>')
        type_species_pattern = re.compile(
            '\"?<I>([a-zA-Z]+)</I> <I>([a-zA-Z]+)</I>')
        type_proposed_pattern = re.compile(
            '<b>Proposed as:</b>(.*)</p>')

        alt_genus_pattern_a = re.compile(
            r'<b>Name:</b>\s?\[([a-zA-Z]+)')
        alt_genus_pattern_b = re.compile(
            r'<b>Name:</b> \"?<I>Candidatus</I> ([a-zA-Z]+)')
        full_list_type_species = []
        for letter in list(string.ascii_uppercase):
            break
            for file in glob.glob(os.path.join(os.path.join(input_dir, 'all_genera'), letter, "*")):
                genus_name = ''
                genus_reference = ''
                genus_proposed_type = ''

                name_section = False
                type_proposed_section = False
                type_species_section = False

                type_species_list = []

                with open(file) as f:
                    for line in f:
                        line = line.replace("'", '"')
                        if 'glossary#citation-of-names-and-authors' in line:
                            genus_pattern_results = genus_pattern.search(
                                line)
                            if genus_pattern_results and genus_pattern_results.group(1) != 'Candidatus':
                                genus_name = genus_pattern_results.group(1)
                                name_section = True
                            else:
                                alt_genus_pattern_a_results = alt_genus_pattern_a.search(
                                    line)
                                if alt_genus_pattern_a_results:
                                    genus_name = alt_genus_pattern_a_results.group(
                                        1)
                                    name_section = True
                                else:
                                    alt_genus_pattern_b_results = alt_genus_pattern_b.search(
                                        line)
                                    if alt_genus_pattern_b_results:
                                        genus_name = alt_genus_pattern_b_results.group(
                                            1)
                                        name_section = True
                            if name_section:
                                genus_reference = self.cleanhtml(
                                    line.split(genus_name, 1)[1]).replace(']', '').replace('"', '').strip()

                        elif "glossary#abbreviations-in-proposals" in line:
                            type_proposed_section = True
                            type_proposed_pattern_results = type_proposed_pattern.search(
                                line)
                            if type_proposed_pattern_results:
                                genus_proposed_type = type_proposed_pattern_results.group(
                                    1)

                        elif '<b>Type species:</b>' in line and name_section:
                            type_species_section = True
                        elif '<p' in line and type_species_section and name_section:
                            break
                        elif type_species_section and name_section:
                            type_species_pattern_results = type_species_pattern.search(
                                line)
                            if type_species_pattern_results:
                                type_specie_name = type_species_pattern_results.group(
                                    1) + " " + type_species_pattern_results.group(2)
                                type_species_list.append(type_specie_name)

                ref_type_combined = ','.join(
                    list(filter(None, [genus_reference, genus_proposed_type]))).strip()

                output_file.write("genus\t{}\t{}\t{}\t\n".format(
                    genus_name in full_list_type_genus, genus_name, ref_type_combined))
                full_list_type_species.extend(type_species_list)
        return full_list_type_species

    def parse_species_html(self, output_file, input_dir, full_list_type_species):
        # Pattern for species
        species_pattern = re.compile(
            r'<b>Name:</b> \"?<I>([a-zA-Z]+)</I> <I>([a-zA-Z]+)</I>')
        candidatus_species_pattern = re.compile(
            r'<b>Name:</b> \"?<I>Candidatus</I> ([a-zA-Z]+) ([a-zA-Z]+)"')

        type_proposed_pattern = re.compile(
            '<b>Proposed as:</b>(.*)</p>')

        for letter in list(string.ascii_uppercase):
            for file in glob.glob(os.path.join(input_dir, 'all_species', letter, "*")):
                species_name = ''
                species_reference = ''
                species_strains = ''
                neotype_strains = ''
                species_proposed_type = ''
                name_section = False
                strain_section = False
                type_proposed_section = False
                raw_list_strain = []
                infos = []
                with open(file) as f:
                    for line in f:
                        line = line.replace("'", '"')

                        # if the line contains a genus
                        if 'glossary#citation-of-names-and-authors' in line:
                            if '"<I>Candidatus</I>' not in line:
                                species_pattern_results = species_pattern.search(
                                    line)
                                if species_pattern_results:
                                    species_name = species_pattern_results.group(
                                        1) + " " + species_pattern_results.group(2)
                                    name_section = True
                            else:
                                candidatus_species_pattern_results = candidatus_species_pattern.search(
                                    line)
                                if candidatus_species_pattern_results:
                                    species_name = candidatus_species_pattern_results.group(
                                        1) + " " + candidatus_species_pattern_results.group(2)
                                    name_section = True

                            if name_section:
                                species_reference = self.cleanhtml(
                                    line).split(species_name, 1)[1].replace(']', '').replace('"', '').strip()
                                if '"<I>Candidatus</I>' in line:
                                    print(species_reference)

                        elif "glossary#abbreviations-in-proposals" in line:
                            type_proposed_section = True
                            type_proposed_pattern_results = type_proposed_pattern.search(
                                line)
                            if type_proposed_pattern_results:
                                species_proposed_type = type_proposed_pattern_results.group(
                                    1)

                        elif '<b>Type strains:</b>' in line and name_section:
                            strain_section = True
                        elif '<p' in line and strain_section and name_section:
                            raw_list_strain = list(
                                filter(lambda a: a.strip() != "no culture isolated", raw_list_strain))
                            raw_list_strain = list(
                                filter(lambda a: a.strip() != "no culture available", raw_list_strain))
                            raw_list_strain = list(
                                filter(lambda a: a.strip() != "no pure culture", raw_list_strain))
                            raw_list_strain = self.parse_strains(
                                raw_list_strain)

                            raw_list_strain = [
                                x for x in raw_list_strain if self.check_format_strain(x)]
                            break
                            # print(raw_list_strain)

                        elif strain_section and name_section:
                            if line.strip() != '':
                                raw_list_strain.extend(
                                    [x.replace('strain', '').replace('Strain', '').strip() for x in self.cleanhtml(line).split(';')])

                ref_type_combined = ','.join(
                    list(filter(None, [species_reference, species_proposed_type]))).strip()

                output_file.write("species\t{}\t{}\t{}\t{}\n".format(species_name in full_list_type_species,
                                                                     species_name, ref_type_combined, "=".join(raw_list_strain)))

    def parse_subspecies_html(self, output_file, input_dir, full_list_type_species):
        type_proposed_pattern = re.compile(
            '<b>Proposed as:</b>(.*)</p>')

        for file in glob.glob(os.path.join(input_dir, 'all_subspecies', "*")):

            short_name = ''
            raw_list_strain = []
            strain_section = False

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
                        subspecies_reference = self.cleanhtml(
                            line).split(subspe_name, 1)[1]
                    elif "glossary#abbreviations-in-proposals" in line:
                        type_proposed_section = True
                        type_proposed_pattern_results = type_proposed_pattern.search(
                            line)
                        if type_proposed_pattern_results:
                            subspecies_proposed_type = type_proposed_pattern_results.group(
                                1)
                    elif '<b>Type strains:</b>' in line and name_section:
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
                            x for x in raw_list_strain if self.check_format_strain(x)]
                        break
                        # print(raw_list_strain)

                    elif strain_section and name_section:
                        if line.strip() != '':
                            raw_list_strain.extend(
                                [x.replace('strain', '').replace('Strain', '').strip() for x in self.cleanhtml(line).split(';')])

                ref_type_combined = ','.join(
                    list(filter(None, [subspecies_reference, subspecies_proposed_type]))).strip()

                if short_name == '':
                    output_file.write("species\t{}\t{}\t{}\t{}\n".format(short_name != '' and short_name in full_list_type_species,
                                                                         subspe_name, ref_type_combined, "=".join(raw_list_strain)))

    def parse_html(self, input_dir):
        '''
        Parse the html file of each genus.
        Store the type, the name, the reference, the strains for each species.
        '''
        processed_species = []
        output_file = open(os.path.join(
            self.outdir, 'lpsn_summary.tsv'), 'w')

        full_list_type_genus = self.parse_family_html(output_file, input_dir)
        full_list_type_species = self.parse_genus_html(
            output_file, input_dir, full_list_type_genus)
        self.parse_species_html(output_file, input_dir, full_list_type_species)
        self.parse_subspecies_html(
            output_file, input_dir, full_list_type_species)

        output_file.close()
        self.summarise_parsing(output_file.name)
    #=========================================================================
    #
    #                         genus_result = genus_pattern.search(line)
    #                         if genus_result:
    #                             genus_name = genus_result.group(1)
    #                             genus_ref_result = genus_ref_pattern.search(line)
    #                             if genus_ref_result:
    #                                 genus_reference = self.cleanhtml(
    #                                     genus_ref_result.group(1))
    #                         if 'Type genus of the order' in line:
    #                             is_type_genus = True
    #                         if genus_name != '':
    #                             # write information
    #                             output_file.write('genus\t{}\t{}\t{}\n'.format(
    #                                 is_type_genus, genus_name, self.remove_accents(genus_reference)))
    #                         genus_name = ''
    #                         genus_reference = ''
    #                         is_type_genus = False
    #                     # if the line contains a specie
    #                     elif '<span class="genusspecies">' in line:
    #                         if species_name != '':
    #                             species_name = ''
    #                             species_reference = ''
    #                             species_strains = ''
    #                             neotype_strains = ''
    #                             is_type_species = False
    #                         species_result = species_pattern.search(line)
    #
    #                         alt_species_result = alt_species_pattern.search(line)
    #                         if species_result:
    #                             species_name = species_result.group(
    #                                 1) + ' ' + species_result.group(2)
    #                             species_ref_results = species_ref_pattern.search(
    #                                 line)
    #                             if species_ref_results:
    #                                 species_reference = self.cleanhtml(species_ref_results.group(
    #                                     2))
    #                             subspecies_result = subspecies_pattern.search(line)
    #                             if subspecies_result:
    #                                 species_name = species_name + \
    #                                     " subsp. " + subspecies_result.group(1)
    #                                 subspecies_ref = subspecies_ref_pattern.search(
    #                                     line)
    #                                 if subspecies_ref:
    #                                     species_reference = self.cleanhtml(subspecies_ref.group(
    #                                         2))
    #                         elif alt_species_result:
    #                             species_name = alt_species_result.group(
    #                                 1) + ' ' + alt_species_result.group(2)
    #                             species_ref_results = species_ref_pattern.search(
    #                                 line)
    #                             if species_ref_results:
    #                                 species_reference = self.cleanhtml(species_ref_results.group(
    #                                     2))
    #                             subspecies_result = subspecies_pattern.search(line)
    #                             if subspecies_result:
    #                                 species_name = species_name + \
    #                                     " subsp. " + subspecies_result.group(1)
    #                                 subspecies_ref = subspecies_ref_pattern.search(
    #                                     line)
    #                                 if subspecies_ref:
    #                                     species_reference = self.cleanhtml(subspecies_ref.group(
    #                                         2))
    #
    #                         elif 'specificepithet' not in line and 'subgen' not in line and 'ord. nov.' not in line and 'fam. nov.' not in line:
    #                             genus_result = alt_genus_pattern.search(line)
    #                             if genus_result:
    #                                 genus_name = genus_result.group(1)
    #                                 genus_ref_result = genus_ref_pattern.search(
    #                                     line)
    #                                 if genus_ref_result:
    #                                     genus_reference = self.cleanhtml(
    #                                         genus_ref_result.group(1))
    #                             if 'Type genus of the order' in line:
    #                                 is_type_genus = True
    #                             if genus_name != '':
    #                                 output_file.write('genus\t{}\t{}\t{}\n'.format(
    #                                     is_type_genus, genus_name, genus_reference))
    #                             genus_name = ''
    #                             genus_reference = ''
    #                             is_type_genus = False
    #
    #                         if 'Type species of the genus' in line:
    #                             is_type_species = True
    #                     # if the line contains a list of strains
    #                     elif '<span class="taxon-subhead-typestrain">Type strain:</span>' in line:
    #                         cleanline = self.cleanhtml(line)
    #                         cleanline = self.clean_parenthesis(cleanline)
    #
    #                         cleanline = cleanline.replace(
    #                             "Type strain:", '').replace(
    #                             ".", '').strip()
    #                         species_strains = self.parse_strains(cleanline)
    #                     # if the line contains a list of strains (neotype)
    #                     elif 'taxon-subhead-neotype' in line:
    #                         neotype_result = neotype_pattern.search(line)
    #                         if neotype_result:
    #                             neotype_strains = self.parse_strains(
    #                                 self.cleanhtml(neotype_result.group(2)))
    #                     # write information for each new specie
    #                     elif species_name not in processed_species and species_name != '' and line.startswith('Taxon Spacer'):
    #                         output_file.write('species\t{}\t{}\t{}\t{}\t{}\n'.format(
    #                             is_type_species, self.remove_accents(species_name), self.remove_accents(species_reference), '='.join(species_strains), '='.join(neotype_strains)))
    #                         processed_species.append(species_name)
    #                         species_name = ''
    #                         species_reference = ''
    #                         species_strains = ''
    #                         neotype_strains = ''
    #                         is_type_species = False
    #                 if species_name not in processed_species and species_name != '':
    #                     output_file.write('species\t{}\t{}\t{}\t{}\t{}\n'.format(
    #                         is_type_species, species_name, species_reference, '='.join(species_strains), '='.join(neotype_strains)))
    #                     processed_species.append(species_name)
    #=========================================================================
