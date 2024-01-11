import gzip
import os
import sys
import re
import argparse
from collections import defaultdict
import multiprocessing as mp



class GenomeType(object):
    """Identify genomes marked by NCBI as being a MAG or SAG."""

    def __init__(self,cpus=1):
        self.cpus = cpus

    def run(self, genbank_assembly_summary, refseq_assembly_summary, genome_file, output_file):
        """Identify genomes marked by NCBI as being a MAG or SAG."""

        # parse GenBank and RefSeq assembly summary report
        print('Parsing GenBank and RefSeq assembly reports.')
        gid_genome_type = {}
        for assembly_summary in [genbank_assembly_summary, refseq_assembly_summary]:
            with open(assembly_summary, encoding='utf-8') as f:
                for line in f:
                    if line.startswith('#assembly_accession'):
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


        genomes_to_process = []
        for line in open(genome_file):
            genomes_to_process.append(line.strip())


        # parse GBFF files
        print('Parsing GBFF file for each genome.')
        metagenome_pattern = re.compile(r'derived from(\w|\s)*metagenome', re.I)

        fout = open(output_file, 'w')
        fout.write('genome_id\tncbi_genome_category\tsource\n')
        fout_sanity_check = open(output_file + '.raw', 'w')
        genome_count = 0

        print(f"number of cpus used:{self.cpus}")

        # populate worker queue with data to process
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()
        manager = mp.Manager()
        return_list = manager.list()
        sanity_check_return_list = manager.list()

        for f in genomes_to_process:
            workerQueue.put((f,gid_genome_type,metagenome_pattern))

        for _ in range(self.cpus):
            workerQueue.put(None)

        try:
            workerProc = [mp.Process(target=self.ncbi_genome_cat_worker,
                                     args=(workerQueue, writerQueue, return_list, sanity_check_return_list))
                          for _ in range(self.cpus)]
            writeProc = mp.Process(target=self.__writerThread,
                                   args=(len(genomes_to_process), writerQueue))

            writeProc.start()

            for p in workerProc:
                p.start()

            for p in workerProc:
                p.join()

            writerQueue.put(None)
            writeProc.join()

        except:
            for p in workerProc:
                p.terminate()

            writeProc.terminate()

        list_lines_to_write = [x for x in return_list if x != 'null']
        list_lines_to_write_sc = [x for x in sanity_check_return_list if x != 'null']

        for line in list_lines_to_write:
            fout.write(line)
        for line in list_lines_to_write_sc:
            fout_sanity_check.write(line)

        fout.close()
        fout_sanity_check.close()


        # for idx, line in enumerate(open(genome_file)):
        #     if idx % 100 == 0:
        #         sys.stdout.write('==> Processed %d genomes.\r' % idx)
        #         sys.stdout.flush()
        #


        fout.close()
        fout_sanity_check.close()

    def ncbi_genome_cat_worker(self, queueIn, queueOut,return_list,sc_return_list):
        while True:
            tuple_infos = queueIn.get(block=True, timeout=None)

            if tuple_infos == None:
                break
            line,gid_genome_type,metagenome_pattern=tuple_infos
            sc_value = 'null'
            cat_value = 'null'

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
                            sc_value= f'{gid}\t{line}\n'
                    elif 'single cell' in line:
                        types.add('SAG')
                        sc_value= f'{gid}\t{line}\n'
                    elif '/environmental_sample' in line or "derived from environmental source" in line:
                        types.add('ENV')
                        sc_value= f'{gid}\t{line}\n'

            if 'MAG' in types and 'SAG' in types:
                print('[WARNING] Genome %s is annotated as both a MAG and SAG.' % gid)

            if 'MAG' in types:
                cat_value='%s\t%s\t%s\n' % (gid, 'derived from metagenome', source)
            elif 'SAG' in types:
                cat_value='%s\t%s\t%s\n' % (gid, 'derived from single cell', source)
            elif 'ENV' in types:
                cat_value='%s\t%s\t%s\n' % (gid, 'derived from environmental sample', source)


            return_list.append(cat_value)
            sc_return_list.append(sc_value)
            queueOut.put(cat_value)

    def __writerThread(self, numDataItems, writerQueue):
        """Store or write results of worker threads in a single thread."""

        processedItems = 0
        while True:
            a = writerQueue.get(block=True, timeout=None)
            if a == None:
                break

            processedItems += 1
            statusStr = 'Finished processing %d of %d (%.2f%%) items.' % (processedItems,
                                                                          numDataItems,
                                                                          float(processedItems) * 100 / numDataItems)
            sys.stdout.write('%s\r' % statusStr)
            sys.stdout.flush()

        sys.stdout.write('\n')