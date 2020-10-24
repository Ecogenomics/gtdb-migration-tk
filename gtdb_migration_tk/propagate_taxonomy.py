import os
import sys
import csv
import argparse
from collections import defaultdict
import logging

from gtdb_migration_tk.biolib_lite.taxonomy import Taxonomy
from gtdb_migration_tk.database_configuration import GenomeDatabaseConnectionFTPUpdate
from gtdb_migration_tk.gtdb_lite.gtdb_importer import GTDBImporter

csv.field_size_limit(sys.maxsize)


class Propagate(object):
    """Propagate GTDB taxonomy between NCBI releases."""

    def __init__(self,hostname=None,user=None,password=None,db=None):
        self.password = password
        self.hostname = hostname
        self.user = user
        self.db = db

        self.logger = logging.getLogger('timestamp')


        if db is not None:
            self.temp_con = GenomeDatabaseConnectionFTPUpdate.GenomeDatabaseConnectionFTPUpdate(
                hostname, user, password, db)
            self.temp_con.MakePostgresConnection()
            self.temp_cur = self.temp_con.cursor()

    def propagate_taxonomy(self, gtdb_metadata_prev, gtdb_metadata_cur, taxonomy_file, rep_file):
        """Propagate GTDB taxonomy between NCBI releases."""

        # get GTDB taxonomy for genome in previous release
        self.logger.info('Reading GTDB taxonomy of genome in previous release:')
        prev_gtdb_taxonomy = {}
        prev_gtdb_genomes = set()
        prev_is_rep = set()
        header = True
        for row in csv.reader(open(gtdb_metadata_prev, "rt", encoding='utf-8'),delimiter='\t'):
            if header:
                header = False
                gtdb_taxonomy_index = row.index('gtdb_taxonomy')
                gtdb_rep_index = row.index('gtdb_representative')
            else:
                genome_id = row[0]
                prev_gtdb_genomes.add(genome_id)

                gtdb_taxonomy = row[gtdb_taxonomy_index]
                if gtdb_taxonomy:
                    prev_gtdb_taxonomy[genome_id] = gtdb_taxonomy

                is_rep = (row[gtdb_rep_index] == 't')
                if is_rep:
                    prev_is_rep.add(genome_id)

        self.logger.info('  %d of %d (%.1f%%) genomes in previous NCBI release had a GTDB taxonomy string' % (len(prev_gtdb_taxonomy),
                                                                                             len(prev_gtdb_genomes),
                                                                                             len(
                                                                                                 prev_gtdb_taxonomy) * 100.0 / len(
                                                                                                 prev_gtdb_genomes)))

        self.logger.info('  %d genomes were identified as representatives' % len(prev_is_rep))

        # identify previous representatives in new NCBI release
        self.logger.info('Identifying unchanged genomes in current NCBI release:')
        header = True
        fout = open(taxonomy_file, 'w')
        retained_genomes = set()
        current_genome_ids = []
        prev_rep_count = 0
        cur_reps = set()
        cur_gtdb_taxonomy = {}
        for row in csv.reader(open(gtdb_metadata_cur, "rt", encoding='utf-8'),delimiter='\t'):
            if header:
                header = False

                gtdb_rep_index = row.index('gtdb_representative')
                gtdb_taxonomy_index = row.index('gtdb_taxonomy')
            else:
                genome_id = row[0]
                current_genome_ids.append(genome_id)

                gtdb_taxonomy = row[gtdb_taxonomy_index]
                if gtdb_taxonomy:
                    cur_gtdb_taxonomy[genome_id] = gtdb_taxonomy

                if genome_id in prev_gtdb_genomes:
                    retained_genomes.add(genome_id)
                    if genome_id in prev_gtdb_taxonomy:
                        if prev_gtdb_taxonomy[genome_id] != cur_gtdb_taxonomy[genome_id]:
                            self.logger.info("GTDB taxonomy strings don't match in the two databases:")
                            self.logger.info(cur_gtdb_taxonomy[genome_id])
                            self.logger.info(prev_gtdb_taxonomy[genome_id])
                            sys.exit()

                        fout.write('%s\t%s\n' % (genome_id, prev_gtdb_taxonomy[genome_id]))

                    if genome_id in prev_is_rep:
                        prev_rep_count += 1
                        cur_reps.add(genome_id)

        remaining_prev_genomes = prev_gtdb_genomes - retained_genomes
        self.logger.info('  %d (%.1f%%) genomes unchanged in current NCBI release' % (len(retained_genomes),
                                                                     len(retained_genomes) * 100.0 / len(
                                                                         prev_gtdb_genomes)))
        self.logger.info('  %d (%.1f%%) genomes absent or modified in current NCBI release' % (len(remaining_prev_genomes),
                                                                              len(remaining_prev_genomes) * 100.0 / len(
                                                                                  prev_gtdb_genomes)))
        self.logger.info('  %d representatives unchanged in current GTDB release' % prev_rep_count)

        # try to identify what happened to absent representatives
        self.logger.info('Identifying genomes that have changed databases or version:')

        moved_to_refseq = set()
        moved_to_genbank = set()
        new_genome_version = set()
        for genome_id in current_genome_ids:
            if genome_id.startswith('U_'):
                continue

            # check for database or version change
            cur_version = int(genome_id.split('.')[-1])
            for new_version in range(1, cur_version + 5):
                new_version_id = genome_id.replace('.%d' % cur_version, '.%d' % new_version)
                if new_version_id in remaining_prev_genomes:
                    new_genome_version.add(new_version_id)
                    if new_version_id in prev_gtdb_taxonomy:
                        fout.write('%s\t%s\n' % (genome_id, prev_gtdb_taxonomy[new_version_id]))

                    if new_version_id in prev_is_rep:
                        cur_reps.add(genome_id)
                    continue

                gb_genome_id = new_version_id.replace('RS_GCF', 'GB_GCA')
                if gb_genome_id in remaining_prev_genomes:
                    moved_to_refseq.add(gb_genome_id)
                    if gb_genome_id in prev_gtdb_taxonomy:
                        fout.write('%s\t%s\n' % (genome_id, prev_gtdb_taxonomy[gb_genome_id]))

                    if gb_genome_id in prev_is_rep:
                        cur_reps.add(genome_id)

                    continue

                rs_genome_id = new_version_id.replace('GB_GCA', 'RS_GCF')
                if rs_genome_id in remaining_prev_genomes:
                    moved_to_genbank.add(rs_genome_id)
                    if rs_genome_id in prev_gtdb_taxonomy:
                        fout.write('%s\t%s\n' % (genome_id, prev_gtdb_taxonomy[rs_genome_id]))

                    if rs_genome_id in prev_is_rep:
                        cur_reps.add(genome_id)

                    continue
        fout.close()

        # write out reps
        fout_new_reps = open(rep_file, 'w')
        for genome_id in current_genome_ids:
            if genome_id in cur_reps:
                fout_new_reps.write('%s\t%s\n' % (genome_id, str(True)))
            else:
                fout_new_reps.write('%s\t%s\n' % (genome_id, str(False)))
        fout_new_reps.close()

        self.logger.info('  %d (%.1f%%) genomes moved from GenBank to RefSeq' % (
        len(moved_to_genbank), len(moved_to_genbank) * 100.0 / len(prev_gtdb_genomes)))
        count = 0
        for elem in iter(moved_to_genbank):
            count = count + 1
            if count == 10:
                break
            print(elem)
        self.logger.info('  %d (%.1f%%) genomes moved from RefSeq to GenBank' % (
        len(moved_to_refseq), len(moved_to_refseq) * 100.0 / len(prev_gtdb_genomes)))
        count = 0
        for elem in iter(moved_to_refseq):
            count = count + 1
            if count == 10:
                break
            print(elem)
        self.logger.info('  %d (%.1f%%) genomes have a new version number' % (
        len(new_genome_version), len(new_genome_version) * 100.0 / len(prev_gtdb_genomes)))

        remaining_prev_genomes = remaining_prev_genomes - moved_to_genbank - moved_to_refseq - new_genome_version
        self.logger.info('There are %d genomes not present in the current release.' % len(remaining_prev_genomes))
        self.logger.info('%d of these were representatives.' % len(prev_is_rep.intersection(remaining_prev_genomes)))

    def truncate_taxonomy(self, metadata_file):
        """Truncate taxonomy string to just domain classification."""

        # get current GTDB taxonomy for all genomes
        gtdb_taxonomy = {}
        with open(metadata_file) as f:
            header = f.readline().strip().split('\t')

            gtdb_taxonomy_index = header.index('gtdb_taxonomy')

            for line in f:
                line_split = line.strip().split('\t')

                gid = line_split[0]
                gtdb_taxa = [t.strip() for t in line_split[gtdb_taxonomy_index].split(';')]
                gtdb_taxonomy[gid] = gtdb_taxa

        for i, rank in enumerate(Taxonomy.rank_labels):
            data_to_commit = []
            for gid, taxa in gtdb_taxonomy.iteritems():
                if rank == 'domain':
                    rank_str = taxa[i]
                    data_to_commit.append((gid, rank_str))
                else:
                    data_to_commit.append((gid, Taxonomy.rank_prefixes[i]))

            gtdbimporter = GTDBImporter(self.temp_cur)
            gtdbimporter.importMetadata('metadata_taxonomy', 'gtdb_' + rank, 'TEXT', data_to_commit)
            self.temp_con.commit()

    def add_propagated_taxonomy(self, taxonomy_file, metadata_file, genome_list_file, truncate_taxonomy,rep_id_file):
        """Add taxonomy to database."""

        if truncate_taxonomy:
            self.logger.info('Truncating GTDB taxonomy to domain classification.')
            self.truncate_taxonomy(metadata_file)

        genome_list = set()
        if genome_list_file:
            for line in open(genome_list_file):
                if '\t' in line:
                    genome_list.add(line.rstrip().split('\t')[0])
                else:
                    genome_list.add(line.rstrip().split(',')[0])

        # read taxonomy file
        taxonomy = Taxonomy().read(taxonomy_file)

        # add each taxonomic rank to database
        for i, rank in enumerate(Taxonomy.rank_labels):
            data_to_commit = []
            for genome_id, taxa in taxonomy.items():
                if genome_list_file and genome_id not in genome_list:
                    continue

                rank_str = taxa[i]
                data_to_commit.append((genome_id, rank_str))

            gtdbimporter = GTDBImporter(self.temp_cur)
            gtdbimporter.importMetadata('metadata_taxonomy', 'gtdb_' + rank, 'TEXT', data_to_commit)
            self.temp_con.commit()

        rep_to_commit = []
        with open(rep_id_file) as repfile:
            for line in repfile:
                genome_id,isrep = line.strip().split('\t')
                rep_to_commit.append((genome_id,isrep))
        gtdbimporter = GTDBImporter(self.temp_cur)
        gtdbimporter.importMetadata('metadata_taxonomy', 'gtdb_representative', 'BOOLEAN', rep_to_commit)
        self.temp_con.commit()


    def set_gtdb_domain(self):
        """Set missing GTDB domain information to reflect NCBI domain."""

        self.logger.info('Identifying NCBI genomes with missing domain information.')


        q = ("SELECT id, ncbi_taxonomy FROM metadata_taxonomy "
             + "WHERE (gtdb_domain IS NULL or gtdb_domain = 'd__') and ncbi_taxonomy IS NOT NULL")
        self.temp_cur.execute(q)

        missing_domain_info = []
        for id, ncbi_taxonomy in self.temp_cur:
            ncbi_domain = list(map(str.strip, ncbi_taxonomy.split(';')))[0]
            if ncbi_domain[0:3] != 'd__':
                self.logger.error('NCBI domain has the incorrect prefix: %s' % ncbi_domain)
                sys.exit()

            # gtdb_taxonomy = list(Taxonomy.rank_prefixes)
            # gtdb_taxonomy[0] = ncbi_domain
            # gtdb_taxonomy = ';'.join(gtdb_taxonomy)
            missing_domain_info.append([ncbi_domain, id])

        q = "UPDATE metadata_taxonomy SET gtdb_domain = %s WHERE id = %s"
        self.temp_cur.executemany(q, missing_domain_info)

        self.temp_con.commit()
        self.temp_cur.close()

        self.logger.info('NCBI genomes that were missing GTDB domain info: %d' % len(missing_domain_info))



