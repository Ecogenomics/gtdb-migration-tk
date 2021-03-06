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
from biolib.common import check_file_exists, make_sure_path_exists, check_dir_exists

from gtdb_migration_tk.taxon_utils import canonical_strain_id, check_format_strain


class LPSN(object):
    def __init__(self, skip_taxa_per_letter_dl, lpsn_output_dir):
        """Initialization."""

        self.outdir = lpsn_output_dir
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        self.logger = logging.getLogger('timestamp')

        self.base_url = 'https://lpsn.dsmz.de/'
        self.skip_taxa_per_letter_dl = skip_taxa_per_letter_dl

    def full_lpsn_wf(self):
        self.download_family_lpsn_html()
        self.download_genus_lpsn_html()
        self.download_species_lpsn_html()
        self.download_subspecies_lpsn_html()
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
        fout_type_strains.write('lpsn_strain\tco-identical strain IDs\n')

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

                    # Normalise the strains
                    for i, strain in enumerate(strains):
                        processed_strains.append(canonical_strain_id(strain))
                    processed_neotypes = []
                    processed_strain_string = '{0}\t{1}'.format(
                        line_split[2], "=".join(processed_strains))
                    if processed_strain_string not in list_processed_strains:
                        fout_type_strains.write(
                            '{0}\n'.format(processed_strain_string))
                        list_processed_strains.append(processed_strain_string)

        fout_type_genera.close()
        fout_type_species.close()
        fout_type_strains.close()

    def download_family_lpsn_html(self):
        '''

        Download all HTML family pages from LPSN.

        '''

        # Download pages listing all families in LPSN
        index_dir = os.path.join(self.outdir, 'family_per_letter')
        if self.skip_taxa_per_letter_dl:
            self.logger.info('Skipping download of families at LPSN and using results in: {}'.format(index_dir))
        else:
            self.logger.info('Beginning download of families from LPSN.')
            make_sure_path_exists(index_dir)
            for letter in list(string.ascii_uppercase):
                url = self.base_url + 'family?page=' + letter
                urllib.request.urlretrieve(
                    url, os.path.join(index_dir, 'family_{}.html'.format(letter)))


        # Parse html pages lising all families
        family_sites_list = open(os.path.join(
            self.outdir, 'family_list.lst'), 'w')
        link_pattern = re.compile(r'href="(/family/[^"]+)"')
        family_name_pattern = re.compile(r'class=\"last-child color-family\">(.*)</a>')
        families_to_download = []
        self.logger.info('Generating a list of all pages that contain family information.')
        num_families = 0
        for letter in list(string.ascii_uppercase):
            with open(os.path.join(index_dir, 'family_{}.html'.format(letter))) as webf:
                for line in webf:
                    if 'class="last-child color-family"' in line:
                        line = line.replace("'", '"')
                        result = link_pattern.search(line)
                        if result:
                            family_name = family_name_pattern.search(line).group(1)
                            family_name = self.cleanhtml(family_name)

                            num_families += 1
                            sys.stdout.write(' - processed {:,} families\r'.format(num_families))
                            sys.stdout.flush()
                            family_sites_list.write('{}\t{}\n'.format(
                                letter,family_name, self.base_url + result.group(1)))
                            families_to_download.append(family_name)
        family_sites_list.close()

        sys.stdout.write('\n')

        # we remove the duplicate names with quotes
        valid_families = []
        for famname in families_to_download:

            if famname.startswith('"'):
                if famname.replace('"', '').replace('Candidatus ', '') not in families_to_download:
                    valid_families.append(famname)
            else:
                valid_families.append(famname)

        # Download individual species html page
        self.logger.info('Downloading individual family HTML pages.')
        failed_html_file = open(os.path.join(
            self.outdir, 'families_failed.lst'), 'w')
        make_sure_path_exists(os.path.join(self.outdir, 'all_families'))

        num_families = 0
        num_already_dl = 0
        with open(os.path.join(self.outdir, 'family_list.lst')) as gsl:
            for line in gsl:
                letter, family_name, fam_url = line.strip().split('\t')
                make_sure_path_exists(os.path.join(
                    self.outdir, 'all_families', letter))
                famname = os.path.basename(fam_url)
                out_file = os.path.join(self.outdir, 'all_families', letter, famname)

                if family_name in valid_families:
                    if not os.path.exists(out_file):
                        try:
                            urllib.request.urlretrieve(os.path.join(
                                fam_url), out_file)
                        except:
                            failed_html_file.write('{}\tfailed_download\n'.format(fam_url))
                    else:
                        num_already_dl += 1
                else:
                    failed_html_file.write('{}\tduplicate_name\n'.format(fam_url))

                sys.stdout.write(
                    ' - processed {:,} families, including {:,} that were previously downloaded\r'.format(
                        num_families,
                        num_already_dl))
                sys.stdout.flush()

        failed_html_file.close()
        sys.stdout.write('\n')

    def download_genus_lpsn_html(self):
        """

        Download all HTML genus pages from LPSN.

        """
        # Download pages listing all genus in LPSN
        index_dir = os.path.join(self.outdir, 'genus_per_letter')
        if self.skip_taxa_per_letter_dl:
            self.logger.info('Skipping download of genera at LPSN and using results in: {}'.format(index_dir))
        else:
            self.logger.info('Beginning download of genera from LPSN.')
            make_sure_path_exists(index_dir)
            for letter in list(string.ascii_uppercase):
                url = self.base_url + 'genus?page=' + letter
                urllib.request.urlretrieve(
                    url, os.path.join(index_dir, 'genus_{}.html'.format(letter)))

        # Parse html pages lising all genus
        genus_sites_list = open(os.path.join(
            self.outdir, 'genus_list.lst'), 'w')
        link_pattern = re.compile(r'href="(/genus/[^"]+)"')
        genera_name_pattern = re.compile(r'class=\"last-child color-genus\">(.*)</a>')
        genera_to_download = []
        self.logger.info('Generating a list of all pages that contain genus information.')

        num_genera = 0
        for letter in list(string.ascii_uppercase):
            with open(os.path.join(index_dir, 'genus_{}.html'.format(letter))) as webf:
                for line in webf:
                    if 'class="last-child color-genus"' in line:
                        line = line.replace("'", '"')
                        result = link_pattern.search(line)
                        if result:
                            num_genera += 1
                            sys.stdout.write(' - processed {:,} genera\r'.format(num_genera))
                            sys.stdout.flush()
                            genus_name = genera_name_pattern.search(line).group(1)
                            genus_name = self.cleanhtml(genus_name)
                            genus_sites_list.write('{}\t{}\t{}\n'.format(
                                letter, genus_name, self.base_url + result.group(1)))
                            genera_to_download.append(genus_name)
        genus_sites_list.close()
        sys.stdout.write('\n')

        # we remove the duplicate names with quotes
        valid_genera = []
        for genname in genera_to_download:

            if genname.startswith('"'):
                if genname.replace('"', '').replace('Candidatus ', '') not in genera_to_download:
                    valid_genera.append(genname)
            else:
                valid_genera.append(genname)

        # Download individual species html page
        self.logger.info('Downloading individual genus HTML pages.')
        failed_html_file = open(os.path.join(
            self.outdir, 'genera_failed.lst'), 'w')
        make_sure_path_exists(os.path.join(self.outdir, 'all_genera'))
        num_genera = 0
        num_already_dl = 0
        with open(os.path.join(self.outdir, 'genus_list.lst')) as gsl:
            for line in gsl:
                letter, genus_name, gen_url = line.strip().split('\t')
                make_sure_path_exists(os.path.join(
                    self.outdir, 'all_genera', letter))
                genname = os.path.basename(gen_url)
                out_file = os.path.join(self.outdir, 'all_genera', letter, genname)
                if genus_name in valid_genera:
                    if not os.path.exists(out_file):
                        try:
                            urllib.request.urlretrieve(os.path.join(
                                gen_url), out_file)
                        except:
                            failed_html_file.write('{}\tfailed_download\n'.format(gen_url))
                    else:
                        num_already_dl += 1
                else:
                    failed_html_file.write('{}\tduplicate_name\n'.format(gen_url))
                num_genera += 1
                sys.stdout.write(' - processed {:,} genera, including {:,} that were previously downloaded\r'.format(
                    num_genera, num_already_dl))
                sys.stdout.flush()
        failed_html_file.close()
        sys.stdout.write('\n')

    def download_species_lpsn_html(self):
        """
        Download all html species pages from LPSN.

        """
        # Download pages listing all species in LPSN
        index_dir = os.path.join(self.outdir, 'species_per_letter')
        if self.skip_taxa_per_letter_dl:
            self.logger.info('Skipping download of species at LPSN and using results in: {}'.format(index_dir))
        else:
            self.logger.info('Beginning file of species from LPSN.')
            make_sure_path_exists(index_dir)
            for letter in list(string.ascii_uppercase):
                url = self.base_url + 'species?page=' + letter
                urllib.request.urlretrieve(
                    url, os.path.join(index_dir, 'species_{}.html'.format(letter)))

        # Parse html pages lising all species
        species_sites_list = open(os.path.join(
            self.outdir, 'species_list.lst'), 'w')
        link_pattern = re.compile(r'href="(/species/[^"]+)"')
        species_name_pattern = re.compile(r'class=\"last-child color-species\">(.*)</a>')
        species_to_download = []
        self.logger.info('Generating a list of all pages that contain species information.')

        num_species = 0
        for letter in list(string.ascii_uppercase):
            with open(os.path.join(index_dir, 'species_{}.html'.format(letter))) as webf:
                for line in webf:
                    if 'class="last-child color-species"' in line:
                        line = line.replace("'", '"')
                        result = link_pattern.search(line)
                        if result:
                            num_species += 1
                            sys.stdout.write(' - processed {:,} species\r'.format(num_species))
                            sys.stdout.flush()
                            species_name = species_name_pattern.search(line).group(1)
                            species_name = self.cleanhtml(species_name)
                            species_sites_list.write('{}\t{}\t{}\n'.format(
                                letter, species_name, self.base_url + result.group(1)))
                            species_to_download.append(species_name)
        species_sites_list.close()
        sys.stdout.write('\n')

        # we remove the duplicate names with quotes
        valid_species = []
        for spname in species_to_download:

            if spname.startswith('"'):
                if spname.replace('"', '').replace('Candidatus ', '') not in species_to_download:
                    valid_species.append(spname)
            else:
                valid_species.append(spname)

        # Download individual species html page
        self.logger.info('Downloading individual species HTML pages.')

        failed_html_file = open(os.path.join(
            self.outdir, 'species_failed.lst'), 'w')
        make_sure_path_exists(os.path.join(self.outdir, 'all_species'))

        num_species = 0
        num_already_dl = 0
        with open(os.path.join(self.outdir, 'species_list.lst')) as gsl:
            for line in gsl:
                letter, species_name, spe_url = line.strip().split('\t')
                make_sure_path_exists(os.path.join(
                    self.outdir, 'all_species', letter))
                spenname = os.path.basename(spe_url)
                out_file = os.path.join(self.outdir, 'all_species', letter, spenname)

                if species_name in valid_species:
                    if not os.path.exists(out_file):
                        try:
                            urllib.request.urlretrieve(os.path.join(
                                spe_url), out_file)
                        except:
                            failed_html_file.write('{}\n'.format(spe_url))
                    else:
                        num_already_dl += 1
                else:
                    failed_html_file.write('{}\tduplicate_name\n'.format(spe_url))
                num_species += 1
                sys.stdout.write(' - processed {:,} species, including {:,} that were previously downloaded\r'.format(
                                    num_species, num_already_dl))
                sys.stdout.flush()
        failed_html_file.close()
        sys.stdout.write('\n')

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
                        subspecies_name = self.cleanhtml(subspecies_name)
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

    def cleanhtml(self, raw_html):
        """
        Remove all HTML tags from HTML line
        """
        cleanr = re.compile('<.*?>')
        cleantext = re.sub(cleanr, '', raw_html)
        return cleantext.strip()

    def clean_parenthesis(self, raw_line):
        """
        Remove parenthesis content from HTML line
        """
        # result = re.sub("[\(\[].*?[\)\]]", "", raw_line)
        result = re.sub("\(see[^\)]+\)", "", raw_line)
        return result

    def parse_strains(self, list_strains):
        """
        Parse the HTML line listing all strains and return a list of strains
        """
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
                                #***if '"<I>Candidatus</I>' in line:
                                #***    print(species_reference)

                        elif "glossary#abbreviations-in-proposals" in line:
                            type_proposed_section = True
                            type_proposed_pattern_results = type_proposed_pattern.search(
                                line)
                            if type_proposed_pattern_results:
                                species_proposed_type = type_proposed_pattern_results.group(
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
                            raw_list_strain = self.parse_strains(
                                raw_list_strain)

                            raw_list_strain = [
                                x for x in raw_list_strain if check_format_strain(canonical_strain_id(x))]
                            
                            break

                        elif strain_section and name_section:
                            if line.strip() != '':
                                raw_list_strain.extend(
                                    [x.replace('strain', '').replace('Strain', '').strip() for x in
                                     self.cleanhtml(line).split(';')])

                ref_type_combined = ','.join(
                    list(filter(None, [species_reference, species_proposed_type]))).strip()
                if len(species_name.strip().split(' ')) >= 2:
                    output_file.write("species\t{}\t{}\t{}\t{}\n".format(species_name in full_list_type_species,
                                                                         species_name, ref_type_combined,
                                                                         "=".join(raw_list_strain)))
                else:
                    self.logger.info(f"species name : {species_name} contains less than 2 words")

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
                            x for x in raw_list_strain if check_format_strain(canonical_strain_id(x))]
                        break
                        # print(raw_list_strain)

                    elif strain_section and name_section:
                        if line.strip() != '':
                            raw_list_strain.extend(
                                [x.replace('strain', '').replace('Strain', '').strip() for x in
                                 self.cleanhtml(line).split(';')])

                ref_type_combined = ','.join(
                    list(filter(None, [subspecies_reference, subspecies_proposed_type]))).strip()

                if short_name == '':
                    output_file.write(
                        "species\t{}\t{}\t{}\t{}\n".format(short_name != '' and short_name in full_list_type_species,
                                                           subspe_name, ref_type_combined, "=".join(raw_list_strain)))

    def parse_html(self, input_dir):
        """
        Parse the html file of each genus.
        Store the type, the name, the reference, the strains for each species.
        """
        processed_species = []
        output_file = open(os.path.join(
            self.outdir, 'lpsn_summary.tsv'), 'w')

        self.logger.info('Parsing family pages.')
        full_list_type_genus = self.parse_family_html(output_file, input_dir)

        self.logger.info('Parsing genus pages.')
        full_list_type_species = self.parse_genus_html(
            output_file, input_dir, full_list_type_genus)

        self.logger.info('Parsing species pages.')
        self.parse_species_html(output_file, input_dir, full_list_type_species)

        self.logger.info('Parsing subspecies pages.')
        self.parse_subspecies_html(
            output_file, input_dir, full_list_type_species)

        output_file.close()
        self.summarise_parsing(output_file.name)
