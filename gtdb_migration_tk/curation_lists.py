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

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2020'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import sys
import argparse
import re
import datetime
import os
import logging
import time
import math
from collections import defaultdict, namedtuple

from biolib.common import canonical_gid

class CurationLists(object):
    """Lists and pseudo-trees for new representatives, polyphyletic taxa, rogue genomes, and genomes with modified NCBI names."""
    
    def __init__(self, domain, output_dir):
        """Initialization."""
        
        self.domain = domain
        self.output_dir = output_dir
        self.logger = logging.getLogger('timestamp')
        
    def pseudo_tree(self, gids, out_tree):
        """Create pseudo-tree with the specified genome IDs."""
        
        pseudo_tree = '('
        pseudo_tree += ','.join(gids)
        pseudo_tree += ');'
        
        fout = open(out_tree, 'w')
        fout.write(pseudo_tree)
        fout.close()
        
    def new_gtdb_reps(self,
                        domain_gids,
                        gtdb_sp_clusters,
                        gtdb_prev_sp_clusters):
        """New GTDB representatives."""

        self.logger.info('Identifying previous GTDB representatives.')
        prev_rids = set()
        with open(gtdb_prev_sp_clusters) as f:
            f.readline()
            for line in f:
                tokens = line.strip().split('\t')
                rid = canonical_gid(tokens[0])
                prev_rids.add(rid)
        self.logger.info(' - identified {:,} previous GTDB representatives.'.format(
                            len(prev_rids)))

        self.logger.info('Identifying current GTDB representatives.')
        cur_rids = set()
        with open(gtdb_sp_clusters) as f:
            f.readline()
            for line in f:
                tokens = line.strip().split('\t')
                rid = canonical_gid(tokens[0])
                cur_rids.add(rid)
        self.logger.info(' - identified {:,} current GTDB representatives.'.format(
                            len(cur_rids)))

        self.logger.info('Creating curation list and pseudo-tree of new GTDB representatives.')
        out_file = os.path.join(self.output_dir, f'gids_new_reps.{self.domain}.lst')
        fout = open(out_file, 'w')
        new_rids = set()
        for rid in cur_rids:
            if rid in domain_gids and rid not in prev_rids:
                fout.write('{}\n'.format(rid))
                new_rids.add(rid)
        fout.close()
        self.logger.info(' - identified {:,} new GTDB representatives.'.format(
                            len(new_rids)))
                            
        self.pseudo_tree(new_rids, out_file.replace('.lst', '.tree'))

    def poly_rogue_gtdb_reps(self,
                                domain_gids,
                                taxa_gid_map,
                                gtdb_decorate_table):
        """Polyphyletic and rogue GTDB representatives."""
        
        self.logger.info('Identifying polyphyletic and rogue GTDB representatives.')
        poly_taxa_count = 0
        poly_gids = set()
        rogue_gids = set()
        with open(gtdb_decorate_table) as f:
            f.readline()
            for line in f:
                tokens = line.split('\t')
                
                taxon = tokens[0]
                fmeasure = float(tokens[2])
                rogue_in = tokens[7].strip()
                rogue_out = tokens[8].strip()
                if fmeasure < 1.0:
                    poly_taxa_count += 1
                    poly_gids.update(taxa_gid_map[taxon])
                    
                    if rogue_in:
                        for gid in rogue_in.split(','):
                            gid = canonical_gid(gid.strip())
                            if not gid.startswith('D-'):
                                rogue_gids.add(gid)
                            
                    if rogue_out:
                        for gid in rogue_out.split(','):
                            gid = canonical_gid(gid.strip())
                            if not gid.startswith('D-'):
                                rogue_gids.add(gid)

        self.logger.info(' - identified {:,} polyphyletic taxa spanning {:,} GTDB representatives.'.format(
                            poly_taxa_count,
                            len(poly_gids)))
        self.logger.info(' - identified {:,} rogue GTDB representatives.'.format(
                            len(rogue_gids)))

        self.logger.info('Creating curation lists and pseudo-trees of polyphyletic GTDB representatives.')
        out_file = os.path.join(self.output_dir, f'gids_poly_taxa.{self.domain}.lst')
        fout = open(out_file, 'w')
        for gid in poly_gids:
            fout.write('{}\n'.format(gid))
        fout.close()
        self.pseudo_tree(poly_gids, out_file.replace('.lst', '.tree'))
            
        self.logger.info('Creating curation lists and pseudo-trees of rogue GTDB representatives.')
        out_file = os.path.join(self.output_dir, f'gids_rogues.{self.domain}.lst')
        fout = open(out_file, 'w')
        for gid in rogue_gids:
            fout.write('{}\n'.format(gid))
        fout.close()
        self.pseudo_tree(rogue_gids, out_file.replace('.lst', '.tree'))
                                    
    def run(self,
                gtdb_init_taxonomy,
                gtdb_sp_clusters,
                gtdb_prev_sp_clusters,
                gtdb_decorate_table):
        """Create curation lists and pseudo-trees."""

        # get genomes
        self.logger.info('Identifying taxonomic assignment of genomes.')
        taxa_gid_map = defaultdict(set)
        domain_gids = set()
        for line in open(gtdb_init_taxonomy):
            tokens = line.strip().split('\t')
            gid = canonical_gid(tokens[0])
            
            taxa = [t.strip() for t in tokens[1].split(';')]
            for taxon in taxa:
                taxa_gid_map[taxon].add(gid)
            
            domain_gids.add(gid)
        self.logger.info(' - identified {:,} genomes.'.format(
                            len(domain_gids)))
        
        # new GTDB representatives
        self.new_gtdb_reps(domain_gids,
                            gtdb_sp_clusters,
                            gtdb_prev_sp_clusters)
                            
        # polyphyletic and rogue GTDB representatives
        self.poly_rogue_gtdb_reps(domain_gids,
                                    taxa_gid_map,
                                    gtdb_decorate_table)