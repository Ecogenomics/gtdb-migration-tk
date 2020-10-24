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

import os
import sys
import argparse
import requests
import io
import re
import logging
import time

from requests.auth import HTTPBasicAuth
from unidecode import unidecode


class BacDive(object):
    # ===============================================================================
    # REST Client for PNU Web Services.
    #
    # The attributes of the class are:
    # * ``headers`` -- sets the content-type of the HTTP request-header to json
    # * ``credentials`` -- attaches the username and the password to HTTPBasicAuth for using with the `requests` library
    # ===============================================================================

    def __init__(self, output_dir, username, pwd):
        self.headers = {'Accept': 'application/json'}
        self.credentials = HTTPBasicAuth(username, pwd)
        self.output_dir = output_dir
        self.logger = logging.getLogger('timestamp')

    def getGenera(self, outfile,full_list_of_genera, urlreq=None):
        # to get all genera , but only for the first page.
        # to consider the other pages, you have to change the url to
        # 'https://bacdive.dsmz.de/api/pnu/genus/?page=2' etc.
        genus_type_species_dict = {}

        if urlreq is None:
            response = requests.get(
                'https://bacdive.dsmz.de/api/pnu/genus/', headers=self.headers, auth=self.credentials)
        else:
            while True:
                try:
                    print(urlreq)
                    response = requests.get(
                        urlreq, headers=self.headers, auth=self.credentials)
                except requests.exceptions.ConnectionError:
                    print('Max retries for {}'.format(urlreq))
                    time.sleep(10)
                    continue
                except Exception:
                    print(e)
                    print('Max retries for {}'.format(urlreq))
                    time.sleep(10)
                    continue
                break

        if response.status_code == 200:
            results = response.json()
            listgenus = results.get("results")
            urlreq = results.get("next")
            full_list_of_genera.extend(listgenus)
            for item in listgenus:
                if item.get("label") is not None and item.get('authors') is not None and item.get('taxon') is not None:
                    outfile.write('g__{0}\t\t{1}\n'.format(
                        item.get("label"), item.get('authors') + ',' + item.get('taxon')))
                else:
                    outfile.write('g__{0}\t\t\t\n'.format(item.get("label")))

                if item.get('type_species') is not None:
                    genus_type_species_dict[item.get(
                        'type_species')] = item.get('label')
            if results.get("next") is not None:
                temp_dict,full_list_of_genera = self.getGenera(outfile,full_list_of_genera, urlreq)
                for k, v in temp_dict.items():
                    genus_type_species_dict[k] = v
#             OUTPUT:
#             object of type 'dict' with the fields 'count', 'previous', 'results', 'next'
# the different genera in field 'results' are separated by ',' e.g.
# {genus1},{genus2},{genus3},
            return (genus_type_species_dict,full_list_of_genera)

    def write_all_metadata(self,full_list_of_entities,filename):
        list_keys = []
        for entity in full_list_of_entities:
            for k,v in entity.items():
                list_keys.append(k)
        list_metadata = self.unique_item(list_keys)
        list_metadata.insert(0, list_metadata.pop(list_metadata.index('label')))
        filout = open(os.path.join(self.output_dir,filename),'w')
        filout.write("{}\n".format("\t".join(list_metadata)))
        for entity in full_list_of_entities:
            entity_infos = []
            for key in list_metadata:
                if key in entity:
                    if isinstance(entity.get(key),list):
                        if key == 'literature':
                            reflist =[item.get('reference') for item in entity.get(key)]
                            element = "||".join(reflist)
                        elif key == 'species':
                            speclist = [item[1] for item in entity.get(key)]
                            element = "||".join(speclist)
                        elif key == 'type_strain':
                            type_strain_list = entity.get(key)
                            element = "||".join(type_strain_list)
                        elif key == 'synonyms':
                            synonymslist = []
                            for syn in entity.get(key):
                                allinfos = []
                                for k,v in syn.items():
                                    allinfos.append(f"{k}:{v}")
                                synonymslist.append(";".join(allinfos))
                            element = "||".join(synonymslist)
                        elif key == 'correct_name':
                            element = entity.get(key)[1]
                        else:
                            print('other_ref')
                        entity_infos.append(element)
                    else:
                        entity_infos.append(entity.get(key))
                else:
                    entity_infos.append("None")
            filout.write("\t".join([str(x).replace("<i>","").replace("</i>","") for x in entity_infos])+"\n")
        filout.close()
        print("done")


    def unique_item(self,seq):
        seen = set()
        seen_add = seen.add
        return [x for x in seq if not (x in seen or seen_add(x))]

    def getSpecies(self, outfile_species,full_list_of_species, outfile_strains, dictgenus, urlreq=None):
        # to get a list of all species , but only for the first page.
        # to consider the other pages, you have to change the url to
        # 'https://bacdive.dsmz.de/api/pnu/species/?page=2' etc.
        if urlreq is None:
            response = requests.get(
                'https://bacdive.dsmz.de/api/pnu/species/', headers=self.headers, auth=self.credentials)
        else:
            while True:
                try:
                    print(urlreq)
                    response = requests.get(
                        urlreq, headers=self.headers, auth=self.credentials)
                except requests.exceptions.ConnectionError:
                    print('Max retries for {}'.format(urlreq))
                    time.sleep(10)
                    continue
                except Exception:
                    print(e)
                    print('Max retries for {}'.format(urlreq))
                    time.sleep(10)
                    continue
                break

        if response.status_code == 200:
            results = response.json()
            listspe = results.get("results")
            full_list_of_species.extend(listspe)
            urlreq = results.get("next")
            for item in listspe:
                # if 'subsp.' in item.get("label"):
                #    continue
                if item.get("type_strain") is not None:
                    list_strains = item.get("type_strain")
                    for st in item.get("type_strain"):
                        if " " in st:
                            list_strains.append(st.replace(" ", ''))
                        if "-" in st:
                            list_strains.append(st.replace("-", ''))
                            list_strains.append(st.replace("-", ' '))
                        p = re.compile(
                            '(\w+)\s\(now\s(\w+)\)\s(\d+)', re.IGNORECASE)
                        matches = p.search(st)
                        if matches:
                            list_strains.append("{} {}".format(
                                matches.group(1), matches.group(3)))
                            list_strains.append("{}{}".format(
                                matches.group(1), matches.group(3)))
                            list_strains.append("{} {}".format(
                                matches.group(2), matches.group(3)))
                            list_strains.append("{}{}".format(
                                matches.group(2), matches.group(3)))

                    outfile_strains.write('{0}\t{1}\n'.format(
                        item.get("species"), "=".join(list_strains)))
#                else:
#                    outfile_strains.write('{0} \n'.format(item.get("species")))

                label = ''
                species_authority = ''
                type_spe = ''
                if item.get("label") is not None:
                    label = 's__' + item.get("label")
                if item.get('authors') is not None and item.get('taxon') is not None:
                    species_authority = item.get(
                        'authors') + ',' + item.get('taxon')
                    species_authority = unidecode(species_authority)
                if item.get("label") in dictgenus:
                    type_spe = 'g__' + dictgenus.get(item.get('label'))
                outfile_species.write('{0}\t{1}\t{2}\n'.format(
                    label, type_spe, species_authority))
            # if urlreq == 'https://bacdive.dsmz.de/api/pnu/species/?page=10':
            #     return full_list_of_species
            if results.get("next") is not None:
                full_list_of_species = self.getSpecies(outfile_species,full_list_of_species,
                                outfile_strains, dictgenus, urlreq)

#             OUTPUT:
#             object of type 'dict' with the fields 'count', 'previous', 'results', 'next'
# the species in field 'results' are separated by ',' e.g.
# {species1},{species2},{species3},etc
            return full_list_of_species

    def download_strains(self):

        # test connection
        response = requests.get(
            'https://bacdive.dsmz.de/api/pnu/genus/', headers=self.headers, auth=self.credentials)
        if response.status_code != 200:
            self.logger.info(f'Impossible to log the Bacdive server. Status code = {response.status_code}.')
            sys.exit()

        outfile_genera = open(os.path.join(
            self.output_dir, 'dsmz_genera.tsv'), 'w')
        outfile_genera.write(
            "dsmz_genus\tdsmz_type_genus\tdsmz_genus_authority\n")
        self.logger.info('Parsing genera....')
        full_list_of_genera = []
        dictgenus,full_list_of_genera = self.getGenera(outfile_genera,full_list_of_genera)
        self.write_all_metadata(full_list_of_genera,"metadata_genera.tsv")
        outfile_genera.close()
        self.logger.info('Parsing genera: done.')

        outfile_species = open(os.path.join(
            self.output_dir, 'dsmz_species.tsv'), 'w')
        outfile_species.write(
            "dsmz_species\tdsmz_type_species\tdsmz_species_authority\n")
        outfile_strains = open(os.path.join(
            self.output_dir, 'dsmz_strains.tsv'), 'w')
        outfile_strains.write("dsmz_species\tdsmz_strains\n")
        self.logger.info('Parsing species...')
        full_list_of_species = []
        full_list_of_species = self.getSpecies(outfile_species,full_list_of_species, outfile_strains, dictgenus)
        self.write_all_metadata(full_list_of_species,"metadata_species.tsv")

        outfile_species.close()
        outfile_strains.close()
        self.logger.info('Parsing species: done.')
