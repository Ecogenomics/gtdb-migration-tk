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
                    neotypes = line_split[5].strip().split("=")
                    # Normalise the strains (neotype)
                    if line_split[5] != '' and len(neotypes) > 0:
                        for i, neotype in enumerate(neotypes):
                            if i == 0 and neotype.startswith('strain '):
                                neotype = neotype.replace("strain ", "", 1)
                            neotype = re.sub(r'\(.+\)', ' ', neotype)
                            neotype = ' '.join(neotype.split())
                            b = pattern.sub('', strain).upper()
                            processed_neotypes.append(b)
                    processed_strain_string = '{0}\t{1}\t{2}'.format(
                        line_split[2], "=".join(processed_strains), "=".join(processed_neotypes))
                    if processed_strain_string not in list_processed_strains:
                        fout_type_strains.write(
                            '{0}\n'.format(processed_strain_string))
                        list_processed_strains.append(processed_strain_string)

        fout_type_genera.close()
        fout_type_species.close()
        fout_type_strains.close

    def download_lpsn_html(self):
        '''

        Download all html pages from LPSN.

        '''

        # Download pages listing all genus in LPSN
        print('Beginning file download lpsn ...')
        url = 'http://www.bacterio.net/-ac.html'
        urllib.request.urlretrieve(url, os.path.join(self.outdir, 'ac.html'))
        url = 'http://www.bacterio.net/-dl.html'
        urllib.request.urlretrieve(url, os.path.join(self.outdir, 'dl.html'))
        url = 'http://www.bacterio.net/-mr.html'
        urllib.request.urlretrieve(url, os.path.join(self.outdir, 'mr.html'))
        url = 'http://www.bacterio.net/-sz.html'
        urllib.request.urlretrieve(url, os.path.join(self.outdir, 'sz.html'))

        # Parse html pages lising all genus
        genera_sites_list = open(os.path.join(
            self.outdir, 'genera_list.lst'), 'w')
        link_pattern = re.compile(r'href="(\w+\.html)"')
        replace_link_pattern = re.compile(r'>(\w+)</span>')
        print('Get a list of all the pages that contain genus information ...')
        for webpage in ['ac.html', 'dl.html', 'mr.html', 'sz.html']:
            with open(os.path.join(self.outdir, webpage)) as webf:
                for line in webf:
                    if 'genusspecies' in line:
                        line = line.replace("'", '"')
                        result = link_pattern.search(line)
                        if result:
                            genera_sites_list.write(
                                '{}\n'.format(result.group(1)))
                        else:
                            result = replace_link_pattern.search(line)
                            if result:
                                genera_sites_list.write(
                                    '{}.html\n'.format(result.group(1).lower()))
        genera_sites_list.close()

        # Download individual genus html page
        failed_html_file = open(os.path.join(
            self.outdir, 'genus_failed.lst'), 'w')
        make_sure_path_exists(os.path.join(self.outdir, 'genus_html'))
        with open(os.path.join(self.outdir, 'genera_list.lst')) as gsl:
            for line in gsl:
                genus = line.strip()
                try:
                    print(os.path.join('http://www.bacterio.net', genus))
                    urllib.request.urlretrieve(os.path.join(
                        'http://www.bacterio.net', genus), os.path.join(self.outdir, 'genus_html', genus))
                except:
                    failed_html_file.write('{}\n'.format(genus))
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

    def parse_strains(self, line):
        '''
        Parse the HTML line listing all strains and return a list of strains
        '''
        results = []
        # replace (now XXXX ) pattern
        now_pattern = re.compile('(\w+) \(now (\w+)\) (\d+)')
        # replace (formerly XXXX) pattern
        formerly_pattern = re.compile('(\w+) \(formerly (\w+)\) (\d+)')
        for item in line.split('='):
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

            else:
                results.append(item.strip())
        return results

    def remove_accents(self, data):
        '''
        Remove all special characters from reference and species name
        '''
        return ''.join(x for x in unicodedata.normalize('NFKD', html.unescape(data)) if x in string.printable)

    def parse_html(self, input_dir):
        '''
        Parse the html file of each genus.
        Store the type, the name, the reference, the strains for each species.
        '''
        processed_species = []
        output_file = open(os.path.join(
            self.outdir, 'lpsn_summary.tsv'), 'w')

        # Patterns for genus
        genus_pattern = re.compile(
            r'class="genus">(?:<a name=\"?[^\"]*\"? id=\"?[^\"]*\"?>)?(?:<\/a>)?(\w+)<')
        alt_genus_pattern = re.compile(
            r'class="genusspecies">(?:<\/a>)?(\w+)<')
        genus_ref_pattern = re.compile(r'</span>(.+)')

        # Pattern for species
        species_pattern = re.compile(
            r'<span class="genusspecie[s]?"><a name="?.*"? id="?.*"?></a>\s*(\w+)</span>\s*<span class="(?:specificepithet|genusspecies)">\s*(\w+)</span>')
        alt_species_pattern = re.compile(
            r'<span class="genusspecie[s]?"><a name="?.*"? id="?.*"?></a>\s*(\w+)\s+(\w+)</span>')
        species_ref_pattern = re.compile(
            r'class="(?:specificepithet|genusspecies)">([^<]*)<\/span>(.*?(?=<br>)|(.*))')
        subspecies_pattern = re.compile(
            r'<span class="subspecificepithet">([^<]*)</span>')
        subspecies_ref_pattern = re.compile(
            r'class="subspecificepithet">([^<]*)<\/span>(.*?(?=<br>)|(.*))')
        neotype_pattern = re.compile(
            r'<span class="taxon-subhead-neotype">Neotype strain:</span>\s+(strain)?(.*).$')

        for file in glob.glob(os.path.join(input_dir, "*.html")):
            genus_name = ''
            genus_reference = ''
            is_type_genus = False
            species_name = ''
            species_reference = ''
            species_strains = ''
            neotype_strains = ''
            is_type_species = False
            with open(file) as f:
                for line in f:
                    line = line.replace("'", '"')
                    # if the line contains a genus
                    if '<span class="genus">' in line and 'subgen' not in line and 'ord. nov.' not in line and 'fam. nov.' not in line:
                        genus_result = genus_pattern.search(line)
                        if genus_result:
                            genus_name = genus_result.group(1)
                            genus_ref_result = genus_ref_pattern.search(line)
                            if genus_ref_result:
                                genus_reference = self.cleanhtml(
                                    genus_ref_result.group(1))
                        if 'Type genus of the order' in line:
                            is_type_genus = True
                        if genus_name != '':
                            # write information
                            output_file.write('genus\t{}\t{}\t{}\n'.format(
                                is_type_genus, genus_name, self.remove_accents(genus_reference)))
                        genus_name = ''
                        genus_reference = ''
                        is_type_genus = False
                    # if the line contains a specie
                    elif '<span class="genusspecies">' in line:
                        if species_name != '':
                            species_name = ''
                            species_reference = ''
                            species_strains = ''
                            neotype_strains = ''
                            is_type_species = False
                        species_result = species_pattern.search(line)

                        alt_species_result = alt_species_pattern.search(line)
                        if species_result:
                            species_name = species_result.group(
                                1) + ' ' + species_result.group(2)
                            species_ref_results = species_ref_pattern.search(
                                line)
                            if species_ref_results:
                                species_reference = self.cleanhtml(species_ref_results.group(
                                    2))
                            subspecies_result = subspecies_pattern.search(line)
                            if subspecies_result:
                                species_name = species_name + \
                                    " subsp. " + subspecies_result.group(1)
                                subspecies_ref = subspecies_ref_pattern.search(
                                    line)
                                if subspecies_ref:
                                    species_reference = self.cleanhtml(subspecies_ref.group(
                                        2))
                        elif alt_species_result:
                            species_name = alt_species_result.group(
                                1) + ' ' + alt_species_result.group(2)
                            species_ref_results = species_ref_pattern.search(
                                line)
                            if species_ref_results:
                                species_reference = self.cleanhtml(species_ref_results.group(
                                    2))
                            subspecies_result = subspecies_pattern.search(line)
                            if subspecies_result:
                                species_name = species_name + \
                                    " subsp. " + subspecies_result.group(1)
                                subspecies_ref = subspecies_ref_pattern.search(
                                    line)
                                if subspecies_ref:
                                    species_reference = self.cleanhtml(subspecies_ref.group(
                                        2))

                        elif 'specificepithet' not in line and 'subgen' not in line and 'ord. nov.' not in line and 'fam. nov.' not in line:
                            genus_result = alt_genus_pattern.search(line)
                            if genus_result:
                                genus_name = genus_result.group(1)
                                genus_ref_result = genus_ref_pattern.search(
                                    line)
                                if genus_ref_result:
                                    genus_reference = self.cleanhtml(
                                        genus_ref_result.group(1))
                            if 'Type genus of the order' in line:
                                is_type_genus = True
                            if genus_name != '':
                                output_file.write('genus\t{}\t{}\t{}\n'.format(
                                    is_type_genus, genus_name, genus_reference))
                            genus_name = ''
                            genus_reference = ''
                            is_type_genus = False

                        if 'Type species of the genus' in line:
                            is_type_species = True
                    # if the line contains a list of strains
                    elif '<span class="taxon-subhead-typestrain">Type strain:</span>' in line:
                        cleanline = self.cleanhtml(line)
                        cleanline = self.clean_parenthesis(cleanline)

                        cleanline = cleanline.replace(
                            "Type strain:", '').replace(
                            ".", '').strip()
                        species_strains = self.parse_strains(cleanline)
                    # if the line contains a list of strains (neotype)
                    elif 'taxon-subhead-neotype' in line:
                        neotype_result = neotype_pattern.search(line)
                        if neotype_result:
                            neotype_strains = self.parse_strains(
                                self.cleanhtml(neotype_result.group(2)))
                    # write information for each new specie
                    elif species_name not in processed_species and species_name != '' and line.startswith('Taxon Spacer'):
                        output_file.write('species\t{}\t{}\t{}\t{}\t{}\n'.format(
                            is_type_species, self.remove_accents(species_name), self.remove_accents(species_reference), '='.join(species_strains), '='.join(neotype_strains)))
                        processed_species.append(species_name)
                        species_name = ''
                        species_reference = ''
                        species_strains = ''
                        neotype_strains = ''
                        is_type_species = False
                if species_name not in processed_species and species_name != '':
                    output_file.write('species\t{}\t{}\t{}\t{}\t{}\n'.format(
                        is_type_species, species_name, species_reference, '='.join(species_strains), '='.join(neotype_strains)))
                    processed_species.append(species_name)
        output_file.close()

        self.summarise_parsing(output_file.name)
