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

from sqlalchemy import create_engine
from bs4 import BeautifulSoup
from gtdb_migration_tk.biolib_lite.common import make_sure_path_exists
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
        for rk in ['phylum','class','order','family','genus','species']:
            self.download_rank_lpsn_html(rk)
        #self.download_rank_lpsn_html()
        #self.download_class_lpsn_html()
        #self.download_order_lpsn_html()
        #self.download_family_lpsn_html()
        #self.download_genus_lpsn_html()
        #self.download_species_lpsn_html()
        self.download_subspecies_lpsn_html()
        self.parse_html(os.path.join(self.outdir, 'genus_html'))

    def summarise_parsing(self, lpsn_scrape_file,full_list_type_genus):
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
                        family = 'f__' + full_list_type_genus.get(line_split[2])

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

    def download_rank_lpsn_html(self,rank_name):
        '''

        Download all HTML pages from LPSN for a specific rank.

        '''
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
                    url, os.path.join(index_dir, '{}_{}.html'.format(rank_name,letter)))

        # Parse html pages lising all classes
        rank_sites_list = open(os.path.join(
            self.outdir, f'{rank_name}_list.lst'), 'w')
        link_pattern = re.compile(r'href="(/{}/[^"]+)"'.format(rank_name))
        rank_name_pattern = re.compile(r'class=\"last-child color-{0}\">\"?([^\"]*)\"?</a>'.format(rank_name))
        rank_to_download = []
        self.logger.info(f'Generating a list of all pages that contain {rank_name} information.')
        num_ranks = 0

        for letter in list(string.ascii_uppercase):
            with open(os.path.join(index_dir, '{}_{}.html'.format(rank_name,letter))) as webf:
                for line in webf:
                    if f'class="last-child color-{rank_name}"' in line:
                        line = line.replace("'", '"')
                        result = link_pattern.search(line)
                        if result:
                            rk_name = rank_name_pattern.search(line).group(1)
                            rk_name = self.cleanhtml(rk_name)

                            num_ranks += 1
                            print(' - processed {:,} names ({} level).'.format(num_ranks,rank_name),end='\r')
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
                    potential_name = temp_name.split(', no')[0].replace('[','')
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

                num_ranks += 1
                print(
                    ' - processed {:,} names ({} level), including {:,} that were previously downloaded\r'.format(
                        num_ranks,
                        rank_name,
                        num_already_dl),end='\r')

        failed_html_file.close()


    def download_class_lpsn_html(self):
        '''

        Download all HTML class pages from LPSN.

        '''
        # Download pages listing all classes in LPSN
        index_dir = os.path.join(self.outdir, 'class_per_letter')
        if self.skip_taxa_per_letter_dl:
            self.logger.info('Skipping download of classes at LPSN and using results in: {}'.format(index_dir))
        else:
            self.logger.info('Beginning download of classes from LPSN.')
            make_sure_path_exists(index_dir)
            for letter in list(string.ascii_uppercase):
                url = self.base_url + 'class?page=' + letter
                urllib.request.urlretrieve(
                    url, os.path.join(index_dir, 'class_{}.html'.format(letter)))

        # Parse html pages lising all classes
        class_sites_list = open(os.path.join(
            self.outdir, 'class_list.lst'), 'w')
        link_pattern = re.compile(r'href="(/class/[^"]+)"')
        class_name_pattern = re.compile(r'class=\"last-child color-class\">(.*)</a>')
        classes_to_download = []
        self.logger.info('Generating a list of all pages that contain class information.')
        num_classes = 0
        for letter in list(string.ascii_uppercase):
            with open(os.path.join(index_dir, 'class_{}.html'.format(letter))) as webf:
                for line in webf:
                    if 'class="last-child color-class"' in line:
                        line = line.replace("'", '"')
                        result = link_pattern.search(line)
                        if result:
                            class_name = class_name_pattern.search(line).group(1)
                            class_name = self.cleanhtml(class_name)

                            num_classes += 1
                            sys.stdout.write(' - processed {:,} classes\r'.format(num_classes))
                            sys.stdout.flush()
                            class_sites_list.write('{}\t{}\t{}\n'.format(
                                letter, class_name, self.base_url + result.group(1)))
                            classes_to_download.append(class_name)
        class_sites_list.close()

        sys.stdout.write('\n')

        # we remove the duplicate names with quotes
        valid_classes = []
        for claname in classes_to_download:

            if claname.startswith('"'):
                if claname.replace('"', '').replace('Candidatus ', '') not in classes_to_download:
                    valid_classes.append(claname)
            else:
                valid_classes.append(claname)

        # Download individual classes html page
        self.logger.info('Downloading individual class HTML pages.')
        failed_html_file = open(os.path.join(
            self.outdir, 'classes_failed.lst'), 'w')
        make_sure_path_exists(os.path.join(self.outdir, 'all_classes'))

        num_classes = 0
        num_already_dl = 0
        with open(os.path.join(self.outdir, 'class_list.lst')) as gsl:
            for line in gsl:
                letter, class_name, class_url = line.strip().split('\t')
                make_sure_path_exists(os.path.join(
                    self.outdir, 'all_classes', letter))
                claname = os.path.basename(class_url)
                out_file = os.path.join(self.outdir, 'all_classes', letter, claname)

                if class_name in valid_classes:
                    if not os.path.exists(out_file):
                        try:
                            urllib.request.urlretrieve(os.path.join(
                                class_url), out_file)
                        except:
                            failed_html_file.write('{}\tfailed_download\n'.format(class_url))
                    else:
                        num_already_dl += 1
                else:
                    failed_html_file.write('{}\tduplicate_name\n'.format(class_url))

                num_classes += 1
                sys.stdout.write(
                    ' - processed {:,} classes, including {:,} that were previously downloaded\r'.format(
                        num_classes,
                        num_already_dl))
                sys.stdout.flush()

        failed_html_file.close()
        sys.stdout.write('\n')

    def download_order_lpsn_html(self):
        '''

        Download all HTML order pages from LPSN.

        '''
        # Download pages listing all orders in LPSN
        index_dir = os.path.join(self.outdir, 'order_per_letter')
        if self.skip_taxa_per_letter_dl:
            self.logger.info('Skipping download of orders at LPSN and using results in: {}'.format(index_dir))
        else:
            self.logger.info('Beginning download of orders from LPSN.')
            make_sure_path_exists(index_dir)
            for letter in list(string.ascii_uppercase):
                url = self.base_url + 'order?page=' + letter
                urllib.request.urlretrieve(
                    url, os.path.join(index_dir, 'order_{}.html'.format(letter)))

        # Parse html pages lising all orders
        order_sites_list = open(os.path.join(
            self.outdir, 'order_list.lst'), 'w')
        link_pattern = re.compile(r'href="(/order/[^"]+)"')
        order_name_pattern = re.compile(r'class=\"last-child color-order\">(.*)</a>')
        orders_to_download = []
        self.logger.info('Generating a list of all pages that contain order information.')
        num_orders = 0
        for letter in list(string.ascii_uppercase):
            with open(os.path.join(index_dir, 'order_{}.html'.format(letter))) as webf:
                for line in webf:
                    if 'class="last-child color-order"' in line:
                        line = line.replace("'", '"')
                        result = link_pattern.search(line)
                        if result:
                            order_name = order_name_pattern.search(line).group(1)
                            order_name = self.cleanhtml(order_name)

                            num_orders += 1
                            sys.stdout.write(' - processed {:,} orders\r'.format(num_orders))
                            sys.stdout.flush()
                            order_sites_list.write('{}\t{}\t{}\n'.format(
                                letter, order_name, self.base_url + result.group(1)))
                            orders_to_download.append(order_name)
        order_sites_list.close()

        sys.stdout.write('\n')

        # we remove the duplicate names with quotes
        valid_orders = []
        for ordname in orders_to_download:

            if ordname.startswith('"'):
                if ordname.replace('"', '').replace('Candidatus ', '') not in orders_to_download:
                    valid_orders.append(ordname)
            else:
                valid_orders.append(ordname)

        # Download individual classes html page
        self.logger.info('Downloading individual order HTML pages.')
        failed_html_file = open(os.path.join(
            self.outdir, 'orders_failed.lst'), 'w')
        make_sure_path_exists(os.path.join(self.outdir, 'all_orders'))

        num_orders = 0
        num_already_dl = 0
        with open(os.path.join(self.outdir, 'order_list.lst')) as gsl:
            for line in gsl:
                letter, order_name, order_url = line.strip().split('\t')
                make_sure_path_exists(os.path.join(
                    self.outdir, 'all_orders', letter))
                ordname = os.path.basename(order_url)
                out_file = os.path.join(self.outdir, 'all_orders', letter, ordname)

                if order_name in valid_orders:
                    if not os.path.exists(out_file):
                        try:
                            urllib.request.urlretrieve(os.path.join(
                                order_url), out_file)
                        except:
                            failed_html_file.write('{}\tfailed_download\n'.format(order_url))
                    else:
                        num_already_dl += 1
                else:
                    failed_html_file.write('{}\tduplicate_name\n'.format(order_url))

                num_orders +=1
                sys.stdout.write(
                    ' - processed {:,} orders, including {:,} that were previously downloaded\r'.format(
                        num_orders,
                        num_already_dl))
                sys.stdout.flush()

        failed_html_file.close()
        sys.stdout.write('\n')

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
                            family_sites_list.write('{}\t{}\t{}\n'.format(
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
                num_families +=1
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

    def parse_generic_html(self,rank_name,headers_order,all_rank,file):
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
        return headers_order,all_rank


    def parse_rank_html(self,rank_name, input_dir,headers_order,all_rank):

        if rank_name != 'subspecies':
            for letter in list(string.ascii_uppercase):
                for file in glob.glob(os.path.join(os.path.join(input_dir, f'all_{rank_name}'), letter, "*")):
                    headers_order,all_rank=self.parse_generic_html(rank_name,headers_order,all_rank,file)
        else:
            for file in glob.glob(os.path.join(input_dir, 'all_subspecies', "*")):
                headers_order, all_rank = self.parse_generic_html(rank_name, headers_order, all_rank, file)

        return headers_order,all_rank

    def parse_all_ranks_tsv(self,raw_allranks_file):
        parsed_file = open(os.path.join(os.path.dirname(raw_allranks_file),'full_parsing_parsed.tsv'),'w')

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


            name_pattern = re.compile(r'\[([^,]*),.*\]')

            all_infos = []
            for line in fp:
                infos = line.strip().split('\t')
                # rank_of_interest
                roi = infos[rank_index]

                # Parse "Correct Name" so it is just the taxon name without the publication
                if roi in ['species', 'subspecies'] and len(infos[correctname_index].split(" ")) >= 2:
                    if len(infos[correctname_index].split(" ")) > 3 and infos[correctname_index].split(" ")[
                        2] == 'subsp.':
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
                infos[parenttaxon_index] = infos[parenttaxon_index].replace('"', '').replace("[", "").replace("]",
                                                                                                              "").replace(
                    ',', '')

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

                all_infos.append(infos)

        df = pd.DataFrame.from_records(all_infos, columns=headers)

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

    def add_lpsn_metadata(self,hostname,user,password,db,lpsn_file):

        engine_current = create_engine(f'postgresql://{user}:{password}@{hostname}:5432/{db}',
                                       convert_unicode=True,
                                       pool_size=5,
                                       max_overflow=20,
                                       pool_recycle=3600)
        df = pd.read_csv(lpsn_file, sep='\t')
        df.to_sql('lpsn_metadata', engine_current, if_exists='replace')

    def parse_class_html(self, input_dir,output_file_class):
        # Pattern for class
        type_order_pattern = re.compile('color-order">\"?<I>([a-zA-Z]+)</I>')
        class_pattern = re.compile(
            r'<b>Name:</b> \"?\s?<I>([a-zA-Z]+)</I>')
        type_proposed_pattern = re.compile(
            '<b>Proposed as:</b>(.*)</p>')

        alt_class_pattern_a = re.compile(
            r'<b>Name:</b>\s?\[([a-zA-Z]+)')
        alt_class_pattern_b = re.compile(
            r'<b>Name:</b> \"?<I>Candidatus</I> ([a-zA-Z]+)')




        full_list_type_order = {}
        string_to_write = {}
        for letter in list(string.ascii_uppercase):
            for file in glob.glob(os.path.join(os.path.join(input_dir, 'all_classes'), letter, "*")):
                class_name = ''
                class_reference = ''
                name_section = False
                type_order_section = False
                type_order_list = []
                type_order_reference = ''
                with open(file) as f:
                    for line in f:
                        line = line.replace("'", '"')
                        if 'glossary#citation-of-names-and-authors' in line:
                            class_pattern_results = class_pattern.search(
                                line)
                            if class_pattern_results and class_pattern_results.group(1) != 'Candidatus':
                                class_name = class_pattern_results.group(1)
                                name_section = True
                            else:
                                alt_class_pattern_a_results = alt_class_pattern_a.search(
                                    line)
                                if alt_class_pattern_a_results:
                                    class_name = alt_class_pattern_a_results.group(
                                        1)
                                    name_section = True
                                else:
                                    alt_class_pattern_b_results = alt_class_pattern_b.search(
                                        line)
                                    if alt_class_pattern_b_results:
                                        class_name = alt_class_pattern_b_results.group(
                                            1)
                                        name_section = True
                            if name_section:
                                class_reference = self.cleanhtml(
                                    line.split(class_name, 1)[1]).replace(']', '').replace('"', '').strip()


                        if '<b>Type order:</b>' in line:
                            type_order_section = True
                        elif '<p' in line and type_order_section:
                            break
                        elif type_order_section:
                            type_order_pattern_results = type_order_pattern.search(
                                line)
                            if type_order_pattern_results:
                                type_order_name = type_order_pattern_results.group(
                                    1)
                                if type_order_name.lower() == 'candidatus':
                                    alt_type_order_pattern_b = re.compile('color-order">\"?<I>Candidatus</I> ([a-zA-Z]+)')
                                    type_order_pattern_results = alt_type_order_pattern_b.search(
                                        line)
                                    if type_order_pattern_results:
                                        type_order_name = "Candidatus "+type_order_pattern_results.group(
                                            1)
                                type_order_reference = self.cleanhtml(line).split(type_order_name)[-1].replace('"','')
                                type_order_list.append(type_order_name)

                if class_name in string_to_write:
                    if 'not assigned' in string_to_write.get(class_name).get('ref'):
                        # we replace not assigned with full entry
                        string_to_write[class_name] = {'ref':class_reference.replace(', not assigned','not assigned'),
                                               'to_write':'Class\t{}\t{}\t{}\t{}\n'.format('c__' + class_name,
                                                                                '/'.join(['o__' + x for x in
                                                                                          type_order_list]),
                                                                                class_reference.replace(
                                                                                    ', not assigned', 'not assigned'),
                                                                                type_order_reference.strip())}
                else:
                        string_to_write[class_name] = {'ref':class_reference.replace(', not assigned','not assigned'),
                                               'to_write':'Class\t{}\t{}\t{}\t{}\n'.format('c__' + class_name,
                                                                                '/'.join(['o__' + x for x in
                                                                                          type_order_list]),
                                                                                class_reference.replace(
                                                                                    ', not assigned', 'not assigned'),
                                                                                type_order_reference.strip())

                                               }
                for it in type_order_list:
                    full_list_type_order[it] = class_name

        # We remove 'not assigned to class'  if class has another entry with more information
        for key in sorted(string_to_write.keys()):
            output_file_class.write(string_to_write.get(key).get('to_write'))
        return full_list_type_order

    def parse_family_html(self, output_file, input_dir,output_file_family):
        # Pattern for family
        type_genus_pattern = re.compile('color-genus">\"?<I>([a-zA-Z]+)</I>')
        family_pattern = re.compile(
            r'<b>Name:</b> \"?\s?<I>([a-zA-Z]+)</I>')
        type_proposed_pattern = re.compile(
            '<b>Proposed as:</b>(.*)</p>')

        alt_family_pattern_a = re.compile(
            r'<b>Name:</b>\s?\[([a-zA-Z]+)')
        alt_family_pattern_b = re.compile(
            r'<b>Name:</b> \"?<I>Candidatus</I> ([a-zA-Z]+)')




        full_list_type_genus = {}
        string_to_write = {}
        for letter in list(string.ascii_uppercase):
            for file in glob.glob(os.path.join(os.path.join(input_dir, 'all_families'), letter, "*")):
                family_name = ''
                family_reference = ''
                name_section = False
                type_genus_section = False
                type_genus_list = []
                type_genus_reference = ''

                with open(file) as f:
                    for line in f:
                        line = line.replace("'", '"')
                        if 'glossary#citation-of-names-and-authors' in line:
                            family_pattern_results = family_pattern.search(
                                line)
                            if family_pattern_results and family_pattern_results.group(1) != 'Candidatus':
                                family_name = family_pattern_results.group(1)
                                name_section = True
                            else:
                                alt_family_pattern_a_results = alt_family_pattern_a.search(
                                    line)
                                if alt_family_pattern_a_results:
                                    family_name = alt_family_pattern_a_results.group(
                                        1)
                                    name_section = True
                                else:
                                    alt_family_pattern_b_results = alt_family_pattern_b.search(
                                        line)
                                    if alt_family_pattern_b_results:
                                        family_name = alt_family_pattern_b_results.group(
                                            1)
                                        name_section = True
                            if name_section:
                                family_reference = self.cleanhtml(
                                    line.split(family_name, 1)[1]).replace(']', '').replace('"', '').strip()


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
                                if type_genus_name.lower() == 'candidatus':
                                    alt_type_genus_pattern_b = re.compile('color-genus">\"?<I>Candidatus</I> ([a-zA-Z]+)')
                                    type_genus_pattern_results = alt_type_genus_pattern_b.search(
                                        line)
                                    if type_genus_pattern_results:
                                        type_genus_name = "Candidatus "+type_genus_pattern_results.group(
                                            1)
                                type_genus_reference = self.cleanhtml(line).split(type_genus_name)[-1].replace('"','')
                                type_genus_list.append(type_genus_name)

                if family_name in string_to_write:
                    if 'not assigned' in string_to_write.get(family_name).get('ref'):
                        # we replace no assigned with full entry
                        string_to_write[family_name] = {'ref':family_reference.replace(', not assigned','not assigned'),
                                                        'to_write':'Family\t{}\t{}\t{}\t{}\n'.format('f__'+family_name,
                                                              '/'.join(['g__'+x for x in type_genus_list]),
                                                              family_reference.replace(', not assigned','not assigned'),
                                                                           type_genus_reference.strip())
                                                        }
                else:
                    string_to_write[family_name] = {'ref': family_reference.replace(', not assigned', 'not assigned'),
                                                    'to_write': 'Family\t{}\t{}\t{}\t{}\n'.format('f__' + family_name,
                                                                                                  '/'.join(
                                                                                                      ['g__' + x for x
                                                                                                       in
                                                                                                       type_genus_list]),
                                                                                                  family_reference.replace(
                                                                                                      ', not assigned',
                                                                                                      'not assigned'),
                                                                                                  type_genus_reference.strip())}

                # We remove 'not assigned to family'  if family has another entry with more information
                for it in type_genus_list:
                    full_list_type_genus[it] = family_name
        for key in sorted(string_to_write.keys()):
            output_file_family.write(string_to_write.get(key).get('to_write'))
        return full_list_type_genus

    def parse_order_html(self, output_file, input_dir,output_file_order):
        # Pattern for family
        type_genus_pattern = re.compile('color-genus">\"?<I>([a-zA-Z]+)</I>')
        order_pattern = re.compile(
            r'<b>Name:</b> \"?\s?<I>([a-zA-Z]+)</I>')
        type_proposed_pattern = re.compile(
            '<b>Proposed as:</b>(.*)</p>')

        alt_order_pattern_a = re.compile(
            r'<b>Name:</b>\s?\[([a-zA-Z]+)')
        alt_order_pattern_b = re.compile(
            r'<b>Name:</b> \"?<I>Candidatus</I> ([a-zA-Z]+)')


        full_list_type_genus = {}
        string_to_write = {}

        for letter in list(string.ascii_uppercase):
            for file in glob.glob(os.path.join(os.path.join(input_dir, 'all_orders'), letter, "*")):
                order_name = ''
                order_reference = ''
                name_section = False
                type_genus_section = False
                type_genus_list = []
                type_genus_reference = ''
                with open(file) as f:
                    for line in f:
                        line = line.replace("'", '"')
                        if 'glossary#citation-of-names-and-authors' in line:
                            order_pattern_results = order_pattern.search(
                                line)
                            if order_pattern_results and order_pattern_results.group(1) != 'Candidatus':
                                order_name = order_pattern_results.group(1)
                                name_section = True
                            else:
                                alt_order_pattern_a_results = alt_order_pattern_a.search(
                                    line)
                                if alt_order_pattern_a_results:
                                    order_name = alt_order_pattern_a_results.group(
                                        1)
                                    name_section = True
                                else:
                                    alt_order_pattern_b_results = alt_order_pattern_b.search(
                                        line)
                                    if alt_order_pattern_b_results:
                                        order_name = alt_order_pattern_b_results.group(
                                            1)
                                        name_section = True
                            if name_section:
                                order_reference = self.cleanhtml(
                                    line.split(order_name, 1)[1]).replace(']', '').replace('"', '').strip()


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
                                if type_genus_name.lower() == 'candidatus':
                                    alt_type_genus_pattern_b = re.compile(
                                        'color-genus">\"?<I>Candidatus</I> ([a-zA-Z]+)')
                                    type_genus_pattern_results = alt_type_genus_pattern_b.search(
                                        line)
                                    if type_genus_pattern_results:
                                        type_genus_name = "Candidatus " + type_genus_pattern_results.group(
                                            1)
                                type_genus_reference = self.cleanhtml(line).split(type_genus_name)[-1].replace('"', '')
                                type_genus_list.append(type_genus_name)

                if order_name in string_to_write:
                    if 'not assigned' in string_to_write.get(order_name).get('ref'):
                        # we replace not assigned with full entry
                        string_to_write[order_name] = {'ref':order_reference.replace(', not assigned','not assigned'),
                                                       'to_write':'Order\t{}\t{}\t{}\t{}\n'.format('o__'+order_name,
                                                              '/'.join(['g__'+x for x in type_genus_list]),
                                                              order_reference.replace(', not assigned','not assigned'),
                                                                         type_genus_reference.strip())}
                else:
                    string_to_write[order_name] = {'ref': order_reference.replace(', not assigned', 'not assigned'),
                                                   'to_write': 'Order\t{}\t{}\t{}\t{}\n'.format('o__' + order_name,
                                                                                                '/'.join(
                                                                                                    ['g__' + x for x in
                                                                                                     type_genus_list]),
                                                                                                order_reference.replace(
                                                                                                    ', not assigned',
                                                                                                    'not assigned'),
                                                                                                type_genus_reference.strip())}

                for it in type_genus_list:
                    full_list_type_genus[it] = order_name
        # We remove 'not assigned to class'  if class has another entry with more information
        for key in sorted(string_to_write.keys()):
            output_file_order.write(string_to_write.get(key).get('to_write'))
        return full_list_type_genus

    def parse_genus_html(self, output_file, input_dir,output_file_genus, full_list_type_genus):
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
        string_to_write = {}
        letters=['H']
        #for letter in list(string.ascii_uppercase):
        for letter in letters:

            for file in glob.glob(os.path.join(os.path.join(input_dir, 'all_genera'), letter, "*")):
                genus_name = ''
                genus_reference = ''
                genus_proposed_type = ''

                name_section = False
                type_proposed_section = False
                type_species_section = False
                type_species_reference = ''
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

                                if type_specie_name.lower().startswith('candidatus '):
                                    alt_type_species_pattern_b = re.compile(
                                        'color-species">\"?<I>Candidatus</I> ([a-zA-Z]+) ([a-zA-Z]+)')
                                    type_species_pattern_results = alt_type_species_pattern_b.search(
                                        line)
                                    if type_species_pattern_results:
                                        type_specie_name = "Candidatus " + type_species_pattern_results.group(
                                    1) + " " + type_species_pattern_results.group(2)
                                if genus_name == 'Hansschlegelia':
                                    a=1
                                type_species_reference = self.cleanhtml(line).split(type_specie_name)[-1].replace('"', '')


                ref_type_combined = ','.join(
                    list(filter(None, [genus_reference, genus_proposed_type]))).strip()

                if genus_name in string_to_write:
                    if 'not assigned' in string_to_write.get(genus_name).get('ref'):
                        string_to_write[genus_name] = {'ref':ref_type_combined.replace(', not assigned','not assigned'),
                                                       'to_write':'Genus\t{}\t{}\t{}\t{}\n'.format('g__'+genus_name,
                                                              '/'.join(['s__'+x for x in type_species_list]),
                                                              ref_type_combined.replace(', not assigned','not assigned'),
                                                                         type_species_reference.strip())}
                else:
                    string_to_write[genus_name] = {'ref': ref_type_combined.replace(', not assigned', 'not assigned'),
                                                   'to_write': 'Genus\t{}\t{}\t{}\t{}\n'.format('g__' + genus_name,
                                                                                                '/'.join(
                                                                                                    ['s__' + x for x in
                                                                                                     type_species_list]),
                                                                                                ref_type_combined.replace(
                                                                                                    ', not assigned',
                                                                                                    'not assigned'),
                                                                                                type_species_reference.strip())}

                output_file.write("genus\t{}\t{}\t{}\t\n".format(
                    genus_name in full_list_type_genus, genus_name, ref_type_combined))
                full_list_type_species.extend(type_species_list)
        # We remove 'not assigned to class'  if class has another entry with more information
        for key in sorted(string_to_write.keys()):
            output_file_genus.write(string_to_write.get(key).get('to_write'))
        return full_list_type_species

    def parse_species_html(self, output_file, input_dir,output_file_species, full_list_type_species):
        # Pattern for species
        species_pattern = re.compile(
            r'<b>Name:</b> \"?<I>([a-zA-Z]+)</I> <I>([a-zA-Z]+)</I>')
        candidatus_species_pattern = re.compile(
            r'<b>Name:</b> \"?<I>Candidatus</I> ([a-zA-Z]+) ([a-zA-Z]+)"')

        type_proposed_pattern = re.compile(
            '<b>Proposed as:</b>(.*)</p>')

        #for letter in list(string.ascii_uppercase):
        for letter in ['S']:
            for file in glob.glob(os.path.join(input_dir, 'all_species', letter, "*")):
                species_name = ''
                species_reference = ''
                species_strains = ''
                neotype_strains = ''
                species_proposed_type = ''
                name_section = False
                strain_section = False
                correct_section = False
                type_proposed_section = False
                raw_list_strain = []
                infos = []
                skip_species = False
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
                                    if species_name == 'Streptomyces netropsis':
                                        a = 1
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

                            strain_section = False


                        elif strain_section and name_section:
                            if line.strip() != '':
                                raw_list_strain.extend(
                                    [x.replace('strain', '').replace('Strain', '').strip() for x in
                                     self.cleanhtml(line).split(';')])

                        if '<p class="corr-name">' in line:
                            correct_section = True
                        elif correct_section:
                            correct_species_pattern = re.compile(r'color-species">\"?<I>([a-zA-Z]+)</I> <I>([a-zA-Z]+)</I>')
                            correct_species_results = correct_species_pattern.search(line)
                            if correct_species_results:

                                correct_species_name = correct_species_results.group(
                                    1) + " " + correct_species_results.group(2)
                                if correct_species_name == species_name:
                                    print(species_name)
                                    skip_species = True
                                    break
                        if correct_section and 'class="helper"' in line:
                            correct_section = False


                    ref_type_combined = ','.join(
                    list(filter(None, [species_reference, species_proposed_type]))).strip()
                if not skip_species and len(species_name.strip().split(' ')) >= 2:
                    output_file.write("species\t{}\t{}\t{}\t{}\n".format(species_name in full_list_type_species,
                                                                         species_name, ref_type_combined,
                                                                         "=".join(raw_list_strain)))
                    output_file_species.write('Species\t{}\t{}\t{}\tNone\n'.format('s__' + species_name,
                                                                  '/'.join([x for x in raw_list_strain]),
                                                                  ref_type_combined))

                elif not skip_species and len(species_name.strip().split(' ')) < 2:
                    self.logger.info(f"species name : {species_name} contains less than 2 words")

                else:
                    self.logger.info(f"species name : {species_name} skipped because duplicates")

    def parse_subspecies_html(self, output_file, input_dir,output_file_subspecies, full_list_type_species):

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

                    elif '<p class="corr-name">' in line:
                        correct_section = True
                    elif correct_section:
                        correct_subspecies_pattern = re.compile(r'color-subspecies">\"?<I>([a-zA-Z]+)</I> <I>([a-zA-Z]+)</I> subsp. <I>([a-zA-Z]+)</I>')
                        correct_subspecies_results = correct_subspecies_pattern.search(line)
                        if correct_subspecies_results:
                            correct_subspecies_name = correct_subspecies_results.group(
                                1) + " " + correct_subspecies_results.group(2) + " subsp. "+correct_subspecies_results.group(3)
                            if correct_subspecies_name == subspe_name:
                                print(correct_subspecies_name)
                                skip_species = True
                                break
                    elif correct_section and 'class="helper"' in line:
                        correct_section = False

                ref_type_combined = ','.join(
                    list(filter(None, [subspecies_reference, subspecies_proposed_type]))).strip()

                if short_name == '':
                    if not skip_species:
                        output_file.write(
                            "species\t{}\t{}\t{}\t{}\n".format(short_name != '' and short_name in full_list_type_species,
                                                           subspe_name, ref_type_combined, "=".join(raw_list_strain)))
                        output_file_subspecies.write("Subspecies\t{}\t{}\t{}\tNone\n".format('ss__'+subspe_name,
                                                                       "/".join(raw_list_strain),ref_type_combined))




    def parse_html(self, input_dir):
        """
        Parse the html file of each genus.
        Store the type, the name, the reference, the strains for each species.
        """
        processed_species = []
        output_file = open(os.path.join(
            self.outdir, 'lpsn_summary.tsv'), 'w')

        make_sure_path_exists(os.path.join(self.outdir, 'all_ranks'))

        self.logger.info('Parsing all pages.')
        headers_order = ['Rank']
        all_rank = []
        for rk in ['phylum','class','order','family','genus','species','subspecies']:
            self.logger.info(f'Parsing {rk}.')
            headers_order,all_rank = self.parse_rank_html(rk,input_dir, headers_order,all_rank)

        output_all_ranks = open(os.path.join(self.outdir, 'all_ranks', 'full_parsing_raw.tsv'), 'w')
        output_all_ranks.write('\t'.join(headers_order)+'\n')
        for item in all_rank:
            output_all_ranks.write('\t'.join([item.get(potential_header,'n/a') for potential_header in headers_order])+'\n')


        # self.logger.info('Parsing class pages.')
        # full_list_type_order = self.parse_class_html(input_dir)
        #
        # self.logger.info('Parsing order pages.')
        # full_list_type_genus_for_ord = self.parse_order_html(output_file, input_dir)
        #
        self.logger.info('Parsing family pages.')
        full_list_type_genus = self.parse_family_html(output_file, input_dir)

        self.logger.info('Parsing genus pages.')
        full_list_type_species = self.parse_genus_html(
            output_file, input_dir,full_list_type_genus)

        self.logger.info('Parsing species pages.')
        self.parse_species_html(output_file, input_dir, full_list_type_species)

        self.logger.info('Parsing subspecies pages.')
        self.parse_subspecies_html(
            output_file, input_dir, full_list_type_species)

        output_file.close()
        output_all_ranks.close()
        self.parse_all_ranks_data(output_all_ranks)
        self.summarise_parsing(output_file.name)
