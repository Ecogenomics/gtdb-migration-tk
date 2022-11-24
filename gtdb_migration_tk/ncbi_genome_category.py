import gzip
import os
import sys
import re
import argparse
from collections import defaultdict


class GenomeType(object):
    """Identify genomes marked by NCBI as being a MAG or SAG."""

    def __init__(self):
        pass

    def run(self, genbank_assembly_summary, refseq_assembly_summary, genome_file, output_file):
        """Identify genomes marked by NCBI as being a MAG or SAG."""

        # parse GenBank and RefSeq assembly summary report
        print('Parsing GenBank and RefSeq assembly reports.')
        gid_genome_type = {}
        for assembly_summary in [genbank_assembly_summary, refseq_assembly_summary]:
            with open(assembly_summary, encoding='utf-8') as f:
                for line in f:
                    if line.startswith('# assembly_accession'):
                        header = line.strip().split('\t')
                        exclude_idx = header.index('excluded_from_refseq')
                        continue
                    elif line.startswith('#'):
                        continue

                    tokens = [t.strip() for t in line.split('\t')]
                    gid = tokens[0]
                    exclude = tokens[exclude_idx]

                    if 'derived from single cell' in exclude:
                        gid_genome_type[gid] = 'SAG'
                    elif 'derived from metagenome' in exclude:
                        gid_genome_type[gid] = 'MAG'
                    elif 'derived from environmental source' in exclude:
                        gid_genome_type[gid] = 'ENV'
                    elif 'derived from surveillance project' in exclude:
                        pass
                    elif 'derived' in exclude:
                        # seems like an unhandled case
                        print('Unhandled case of derived genome: {}'.format(exclude))
                        sys.exit(-1)

        print(' - identified {:,} MAGs'.format(sum([1 for v in gid_genome_type.values() if v == 'MAG'])))
        print(' - identified {:,} SAGs'.format(sum([1 for v in gid_genome_type.values() if v == 'SAG'])))
        print(
            ' - identified {:,} environmental genomes'.format(sum([1 for v in gid_genome_type.values() if v == 'ENV'])))

        # parse GBFF files
        print('Parsing GBFF file for each genome.')
        metagenome_pattern = re.compile(r'derived from(\w|\s)*metagenome', re.I)

        fout = open(output_file, 'w')
        fout.write('genome_id\tncbi_genome_category\tsource\n')
        fout_sanity_check = open(output_file + '.raw', 'w')
        genome_count = 0
        for idx, line in enumerate(open(genome_file)):
            if idx % 100 == 0:
                sys.stdout.write('==> Processed %d genomes.\r' % idx)
                sys.stdout.flush()

            gid, genome_dir, _fgid = line.strip().split('\t')

            types = set()

            if gid.startswith('U_'):
                source = 'user genome'
                types.add('MAG')
            elif gid in gid_genome_type:
                source = 'assembly report'
                types.add(gid_genome_type[gid])
            else:
                source = 'GBFF file'
                assembly_id = os.path.basename(os.path.normpath(genome_dir))
                wgs_file = os.path.join(genome_dir, assembly_id + '_genomic.gbff.gz')

                open_file = open
                if wgs_file.endswith('.gz'):
                    open_file = gzip.open

                for line in open_file(wgs_file, 'rt'):
                    if 'metagenome' in line:
                        if ('/metagenome_source' in line
                                or re.search(metagenome_pattern, line)):
                            types.add('MAG')
                            fout_sanity_check.write('%s\t\%s' % (gid, line))
                    elif 'single cell' in line:
                        types.add('SAG')
                        fout_sanity_check.write('%s\t\%s' % (gid, line))
                    elif '/environmental_sample' in line or "derived from environmental source" in line:
                        types.add('ENV')
                        fout_sanity_check.write('%s\t\%s' % (gid, line))

            if 'MAG' in types and 'SAG' in types:
                print('[WARNING] Genome %s is annotated as both a MAG and SAG.' % gid)

            if 'MAG' in types:
                fout.write('%s\t%s\t%s\n' % (gid, 'derived from metagenome', source))
            elif 'SAG' in types:
                fout.write('%s\t%s\t%s\n' % (gid, 'derived from single cell', source))
            elif 'ENV' in types:
                fout.write('%s\t%s\t%s\n' %
                           (gid, 'derived from environmental sample', source))
        sys.stdout.flush()

        fout.close()
        fout_sanity_check.close()
