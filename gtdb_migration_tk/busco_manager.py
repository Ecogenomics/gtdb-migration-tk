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
import logging
import json
import shutil
import tempfile
import multiprocessing as mp

from gtdblib.util.shell.execute import run_bash, check_on_path


class BuscoManager(object):
    """Estimate quality of fungal genomes using BUSCO."""

    def __init__(self,cpus=1):
        """Initialization."""

        self.BUSCO_DATASET_DIR = '/srv/db/busco/obd10/'
        self.BUSCO_FUNGI_DB = 'fungi_odb10'
        self.BUSCO_CANARY_DIR = './tmp_busco_canary'
        os.makedirs(self.BUSCO_CANARY_DIR, exist_ok=True)

        self.CPUS_PER_GENOME = 8

        self.cpus = cpus
        self.logger = logging.getLogger('timestamp')
        check_on_path('busco')

    def worker(self, queue_in, queue_out):
        """Process genomes with BUSCO in parallel."""

        while True:
            gid, genome_dir = queue_in.get(block=True, timeout=None)
            if gid == None:
                break

            accn = os.path.basename(os.path.normpath(genome_dir))
            genome_file = os.path.join(genome_dir, f'{accn}_genomic.fna.gz')

            # create canary file to indicate genome is being processed
            canary_file = os.path.join(self.BUSCO_CANARY_DIR, f'{gid}.canary')
            if os.path.exists(canary_file):
                queue_out.put(gid)
                continue

            open(canary_file, 'wt')

            with tempfile.TemporaryDirectory() as tmpdir:
                busco_out_dir = os.path.join(genome_dir, 'busco')
                specific_out_file = os.path.join(busco_out_dir, 'specific_summary.json')

                if not os.path.exists(specific_out_file):
                    # remove any previous results since they can't be
                    # trusted if a specific summary file doesn't exist
                    if os.path.exists(busco_out_dir):
                        shutil.rmtree(busco_out_dir)
                    os.makedirs(busco_out_dir)

                    # copy genomic FASTA file to temporary file and decompress
                    out_genome_file = os.path.join(tmpdir, f'{gid}.fna.gz')
                    shutil.copyfile(genome_file, out_genome_file)
                    cmd = f'pigz -f -d {out_genome_file}'
                    run_bash(cmd)
                    
                    # run BUSCO
                    input_genome = os.path.splitext(out_genome_file)[0]
                
                    cmd = 'busco'
                    cmd += f' -i {input_genome}'
                    cmd += ' -m genome'
                    cmd += ' --auto-lineage-euk'
                    cmd += f' --download_path {self.BUSCO_DATASET_DIR}'
                    cmd += ' --offline'
                    cmd += f' -c {self.CPUS_PER_GENOME}'
                    cmd += f' --out_path {tmpdir}'
                    cmd += f' -o {gid}'

                    try:
                        run_bash(cmd)
                    except:
                        # BUSCO was found to fail for a small number of genomes 
                        # when using the auto lineage mode. For such genomes,
                        # the fungal marker set is used.
                        shutil.rmtree(os.path.join(tmpdir, gid)) 

                        fout = open(os.path.join(busco_out_dir, 'fungal_db.info'), 'w')
                        fout.write(f'BUSCO ran with {self.BUSCO_FUNGI_DB} DB as auto lineage mode failed.')
                        fout.close()

                        cmd = 'busco'
                        cmd += f' -i {input_genome}'
                        cmd += ' -m genome'
                        cmd += f' -l {self.BUSCO_FUNGI_DB}'
                        cmd += f' --download_path {self.BUSCO_DATASET_DIR}'
                        cmd += ' --offline'
                        cmd += f' -c {self.CPUS_PER_GENOME}'
                        cmd += f' --out_path {tmpdir}'
                        cmd += f' -o {gid}'

                        run_bash(cmd)

                    finally:
                        # copy BUSCO logs to genome directory after removing any 
                        # previous results regardless of BUSCO running successfully
                        busco_results_dir = os.path.join(os.path.join(tmpdir, gid))
                        shutil.copytree(os.path.join(busco_results_dir, 'logs'), os.path.join(busco_out_dir, 'logs'))

                    # copy files that will only be valid if BUSCO ran without issue
                    for f in os.listdir(busco_results_dir):
                        if f.startswith('short_summary.generic.') and f.endswith('.json'):
                            in_file = os.path.join(busco_results_dir, f)
                            out_file = os.path.join(busco_out_dir, 'generic_summary.json')
                            shutil.copyfile(in_file, out_file)
                        elif f.startswith('short_summary.specific.') and f.endswith('.json'):
                            in_file = os.path.join(busco_results_dir, f)
                            shutil.copyfile(in_file, specific_out_file)

            queue_out.put(gid)

    def writer(self, num_genomes, writer_queue):
        """Track progress."""

        processed_items = 0
        while True:
            gid = writer_queue.get(block=True, timeout=None)
            if gid == None:
                break

            processed_items += 1
            statusStr = ' -> finished processing {:,} of {:,} ({:.2f}%) genomes'.format(
                processed_items, 
                num_genomes, 
                float(processed_items)*100/num_genomes)
            sys.stdout.write(f'{statusStr}\r')
            sys.stdout.flush()

        sys.stdout.write('\n')

    def run_busco(self, gtdb_genome_path_file, genome_report, output_dir, all_genomes=False):
        """Applying BUSCO to genomes."""

        if genome_report is not None:
            self.logger.error("Support for 'genome_report' not yet implemented.")
            sys.exit(1)

        # process genomes through BUSCO in parallel
        GENOMES_IN_PARALLEL = self.cpus // self.CPUS_PER_GENOME
        self.logger.info(f'Processing {GENOMES_IN_PARALLEL:,} genomes in parallel with {self.CPUS_PER_GENOME:,} CPUs per genome.')
        worker_queue = mp.Queue()
        num_genomes = 0
        with open(gtdb_genome_path_file) as f:
            for line in f:
                gid, genome_dir, _canonical_gid = line.strip().split('\t')

                num_genomes += 1
                worker_queue.put((gid, genome_dir))

        self.logger.info(f'Processing {num_genomes:,} genomes through BUSCO:')

        for _ in range(GENOMES_IN_PARALLEL):
            worker_queue.put((None, None))

        try:
            writer_queue = mp.Queue()

            worker_proc = [mp.Process(target=self.worker, args=(
                worker_queue, 
                writer_queue)) for _ in range(GENOMES_IN_PARALLEL)]
            
            write_proc = mp.Process(
                target=self.writer, args=(num_genomes, writer_queue))

            write_proc.start()

            for p in worker_proc:
                p.start()

            for p in worker_proc:
                p.join()

            writer_queue.put(None)
            write_proc.join()
        except:
            for p in worker_proc:
                p.terminate()

            write_proc.terminate()
                
        # create output file with BUSCO results
        out_file = os.path.join(output_dir, 'busco.tsv')
        self.logger.info(f'Creating BUSCO summary file: {out_file}')
        
        fout = open(out_file, 'w')
        fout.write('gid')
        fout.write('\tlineage_dataset\tnum_markers\tcomplete\tsingle_copy\tmulti_copy\tfragmented\tmissing')
        fout.write('\tgenome_size\tnum_scaffolds\tnum_contigs\t50_scaffolds\tn50_contigs\n')

        missing_summary_file = 0
        progress = 0
        with open(gtdb_genome_path_file) as f:
            for line in f:
                gid, genome_dir, _canonical_gid = line.strip().split('\t')

                busco_dir = os.path.join(genome_dir, 'busco')
                busco_summary_file = os.path.join(busco_dir, 'specific_summary.json')
                if not os.path.exists(busco_summary_file):
                    self.logger.warning(f"Missing BUSCO summary results file: {busco_summary_file}")
                    missing_summary_file += 1
                    continue
                
                data = json.load(open(busco_summary_file))
                fout.write(f'{gid}')
                fout.write('\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                    data['lineage_dataset']['name'],
                    data['results']['n_markers'],
                    data['results']['Complete'],
                    data['results']['Single copy'],
                    data['results']['Multi copy'],
                    data['results']['Fragmented'],
                    data['results']['Missing'],
                ))
                fout.write('\t{}\t{}\t{}\t{}\t{}\n'.format(
                    data['results']['Total length'],
                    data['results']['Number of scaffolds'],
                    data['results']['Number of contigs'],
                    data['results']['Scaffold N50'],
                    data['results']['Contigs N50'],
                ))

                progress += 1
                statusStr = ' -> finished processing {:,} of {:,} ({:.2f}%) genomes'.format(
                    progress, 
                    num_genomes, 
                    float(progress)*100/num_genomes)
                sys.stdout.write(f'{statusStr}\r')
                sys.stdout.flush()

            sys.stdout.write('\n')

        self.logger.warning(f'Missing summary file for {missing_summary_file:,} genomes.')
    
        self.logger.info('Done')