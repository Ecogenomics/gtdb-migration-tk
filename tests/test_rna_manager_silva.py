#!/usr/bin/env python3
"""Offline unit tests for rna_manager_silva.py -- nhmmer and blastn are never run.

Run with the interpreter that has tqdm:

    /opt/centos7/sw/miniconda3/envs/gtdb_migration_tk-r237/bin/python -m unittest discover -s tests -p test_rna_manager_silva.py

genometk_lite's RNA is replaced by a stub that writes the files the real one
writes and records the domain it was given, so what is tested here is the
bookkeeping around the search: which genomes are handed over, which domain each
is searched as, which are reported as unsearched, what lands in the genome
directory and in what order, and what a batch records about itself. The stub is
a module-level class because the workers run in forked processes.
"""

import gzip
import logging
import os
import pickle
import shutil
import tempfile
import unittest
from unittest import mock

from gtdb_migration_tk import batching as B
from gtdb_migration_tk import rna_manager_silva as S
from gtdb_migration_tk.utils import common as C


# A genome whose sequences say this is one the stub refuses, standing in for the
# genomes nhmmer cannot process.
REFUSE = 'REFUSE'

# A genome whose sequences say this has no rRNA gene, and gets only a canary.
NO_GENE = 'NOGENE'

NHMMER_VERSION = 'HMMER 3.4 (Aug 2023)'
BLASTN_VERSION = 'blastn: 2.16.0+'

SILVA_VERSION = '138.2'


class StubRNA(object):
    """Stands in for genometk_lite's RNA, writing what it writes."""

    def __init__(self, rna_name, domain, cpus):
        self.rna_name = rna_name
        self.domain = domain

    def run(self, genome_file, hmm_model_file, db, taxonomy_file, output_dir):
        with gzip.open(genome_file, 'rt') as handle:
            sequences = handle.read()
        if REFUSE in sequences:
            raise RuntimeError('nhmmer error: unknown base type')

        prefix = os.path.join(output_dir, self.rna_name)
        with open('{}.{}.txt'.format(prefix, self.domain), 'w') as handle:
            handle.write('# hmm {}\n'.format(os.path.basename(hmm_model_file)))

        if NO_GENE not in sequences:
            with open(prefix + '.hmm_summary.tsv', 'w') as handle:
                handle.write('Sequence Id\tHMM\n')
            with open(prefix + '.fna', 'w') as handle:
                handle.write('>contig_1\nACGT\n')
            if db is not None:
                with open(prefix + '.blastn.tsv', 'w') as handle:
                    handle.write('contig_1\t4\n')
                with open(prefix + '.taxonomy.tsv', 'w') as handle:
                    handle.write('query_id\ttaxonomy\n')

        with open(prefix + S.CANARY_EXT, 'w') as handle:
            handle.write('done.\n')


class QuietTqdm(object):
    """tqdm with the drawing taken out."""

    def __init__(self, iterable=None, **kwargs):
        self.iterable = iterable

    def __iter__(self):
        return iter(() if self.iterable is None else self.iterable)

    def __enter__(self):
        return self

    def __exit__(self, *exc_info):
        return False

    def update(self, n=1):
        pass

    def close(self):
        pass


class TempDirCase(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='rna_silva_test.')
        self.out_dir = os.path.join(self.dir, 'out')

        for module in (S, B, C):
            patch = mock.patch.object(module, 'tqdm', QuietTqdm)
            patch.start()
            self.addCleanup(patch.stop)

        versions = {S.NHMMER: NHMMER_VERSION, S.BLASTN: BLASTN_VERSION}
        patch = mock.patch.object(S, 'record_program_version',
                                  side_effect=versions.get)
        self.asked = patch.start()
        self.addCleanup(patch.stop)

        for name, value in (('check_dependencies', None), ('RNA', StubRNA)):
            patch = (mock.patch.object(S, name) if value is None
                     else mock.patch.object(S, name, value))
            patch.start()
            self.addCleanup(patch.stop)

        logger = logging.getLogger('timestamp')
        self.addCleanup(logger.setLevel, logger.level)
        logger.setLevel(logging.CRITICAL)

        self.silva_dir = self.silva_files()

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    # ------------------------------------------------------------- the release

    def silva_files(self):
        """The SILVA release directory the constructor insists on."""
        root = os.path.join(self.dir, 'silva')
        release = os.path.join(root, SILVA_VERSION)
        os.makedirs(release)
        for name in ('silva_taxonomy.ssu.tsv', 'silva_taxonomy.lsu.tsv',
                     'SILVA_{}_SSURef_NR99_tax_silva.fasta'.format(SILVA_VERSION),
                     'SILVA_{}_LSURef_tax_silva.fasta'.format(SILVA_VERSION)):
            with open(os.path.join(release, name), 'w') as handle:
                handle.write('x\tBacteria\n')
        return root

    def domain_file(self, predicted):
        """The GTDB domain report; predicted maps an accession to 'ar', 'bac'
        or None, None being a genome the markers gave no prediction for."""
        path = os.path.join(self.dir, 'gtdb_domain_report.tsv')
        with open(path, 'w') as handle:
            handle.write('Genome Id\tNCBI taxonomy\t{}\n'.format(C.DOMAIN_FILE_DOMAIN))
            for accession, domain in predicted.items():
                prefix = 'RS_' if accession.startswith('GCF') else 'GB_'
                handle.write('{}{}\td__Bacteria;p__\t{}\n'.format(
                    prefix, accession,
                    {'ar': 'd__Archaea', 'bac': 'd__Bacteria'}.get(domain,
                                                                   C.NO_PREDICTION)))
        return path

    def taxonomy_file(self, lineages):
        path = os.path.join(self.dir, 'ncbi_taxonomy.tsv')
        with open(path, 'w') as handle:
            for accession, domain in lineages.items():
                handle.write('{}\td__{};p__Whatever;c__;o__;f__;g__;s__\n'.format(
                    accession, 'Archaea' if domain == 'ar' else 'Bacteria'))
        return path

    def genome_dir(self, accession, fasta=True, refuse=False, no_gene=False,
                   searched=False):
        assembly = '{}_ASM{}v1'.format(accession, accession[4:10])
        gpath = os.path.join(self.dir, assembly)
        os.makedirs(gpath, exist_ok=True)

        if fasta:
            bases = REFUSE if refuse else NO_GENE if no_gene else 'ACGT' * 20
            with gzip.open(os.path.join(gpath, assembly + '_genomic.fna.gz'), 'wt') as handle:
                handle.write('>contig_1\n{}\n'.format(bases))

        if searched:
            results = self.results(gpath)
            os.makedirs(results, exist_ok=True)
            for name in ('ssu' + S.CANARY_EXT, 'ssu.taxonomy.tsv'):
                with open(os.path.join(results, name), 'w') as handle:
                    handle.write('from an earlier run\n')

        return gpath

    def genome_dirs_file(self, genomes):
        path = os.path.join(self.dir, 'genome_dirs.tsv')
        with open(path, 'w') as handle:
            for accession, gpath in genomes:
                handle.write('{}\t{}\t{}\n'.format(accession, gpath, 'G' + accession[4:-2]))
        return path

    # --------------------------------------------------------------- the run

    def manager(self, predicted=None, taxonomy=None, rna_gene='ssu', **kwargs):
        return S.RnaManagerSILVA(SILVA_VERSION, self.silva_dir, rna_gene,
                                 self.domain_file(predicted or {}),
                                 self.taxonomy_file(taxonomy or {}),
                                 cpus=1, tmp_dir=os.path.join(self.dir, 'tmp'),
                                 **kwargs)

    def run_silva(self, genomes, predicted=None, taxonomy=None, rna_gene='ssu',
                  all_genomes=False, remove_prior=False,
                  batch_size=B.DEFAULT_BATCH_SIZE):
        manager = self.manager(predicted, taxonomy, rna_gene, batch_size=batch_size)
        ok = manager.run(self.genome_dirs_file(genomes), self.out_dir,
                         all_genomes, remove_prior)
        return manager, ok

    # ------------------------------------------------------------- reading back

    def results(self, gpath):
        return os.path.join(gpath, S.RESULTS_DIR_FORMAT.format(SILVA_VERSION))

    def searched_as(self, gpath, rna_gene='ssu'):
        """The domain the genome was searched as, or None if never searched."""
        results = self.results(gpath)
        for domain in (S.RNA_ARCHAEA, S.RNA_BACTERIA):
            if os.path.exists(os.path.join(results, '{}.{}.txt'.format(rna_gene, domain))):
                return domain
        return None

    def run_dir(self, rna_gene='ssu'):
        return os.path.join(self.out_dir, S.RESULTS_DIR_FORMAT.format(SILVA_VERSION),
                            rna_gene)

    def batches(self, rna_gene='ssu'):
        return B.batch_dir_names(self.run_dir(rna_gene), S.LAYOUT)

    def release_report(self, rna_gene='ssu'):
        with open(os.path.join(self.run_dir(rna_gene), S.NOT_SEARCHED_RELEASE_NAME)) as handle:
            header = handle.readline().rstrip('\n').split('\t')
            self.assertEqual(tuple(header), S.NOT_SEARCHED_HEADER)
            return sorted(tuple(line.rstrip('\n').split('\t')) for line in handle)


# ------------------------------------------------------- where the domain comes from

class TheDomainOfEachGenome(TempDirCase):
    """The rRNA HMMs are per domain, and a genome told wrong gets a worse answer
    rather than an error."""

    def test_gtdb_s_own_prediction_is_used(self):
        manager = self.manager(predicted={'GCA_000002.1': 'ar', 'GCA_000001.1': 'bac'})

        self.assertEqual(manager.domain('GCA_000002.1'), S.RNA_ARCHAEA)
        self.assertEqual(manager.domain('GCA_000001.1'), S.RNA_BACTERIA)

    def test_the_ncbi_taxonomy_file_answers_for_a_genome_gtdb_could_not_predict(self):
        """Not the NCBI taxonomy column of the domain file, which is what this
        command fell back on until 0.1.35 and which says Bacteria here."""
        manager = self.manager(predicted={'GCA_000002.1': None},
                               taxonomy={'GCA_000002.1': 'ar'})

        self.assertEqual(manager.domain('GCA_000002.1'), S.RNA_ARCHAEA)

    def test_gtdb_s_prediction_wins_over_the_ncbi_taxonomy(self):
        manager = self.manager(predicted={'GCA_000002.1': 'ar'},
                               taxonomy={'GCA_000002.1': 'bac'})

        self.assertEqual(manager.domain('GCA_000002.1'), S.RNA_ARCHAEA)

    def test_a_genbank_genome_finds_the_taxonomy_of_its_refseq_counterpart(self):
        manager = self.manager(taxonomy={'GCF_000002.1': 'ar'})

        self.assertEqual(manager.domain('GCA_000002.1'), S.RNA_ARCHAEA)

    def test_a_genome_neither_file_answers_for_is_searched_as_a_bacterium(self):
        """It used to stop the run where the domain file named no domain."""
        manager = self.manager()

        self.assertEqual(manager.domain('GCA_000009.1'), S.RNA_BACTERIA)

    def test_each_genome_is_searched_with_the_hmm_of_its_domain(self):
        arc = self.genome_dir('GCA_000002.1')
        bac = self.genome_dir('GCA_000001.1')

        self.run_silva([('GCA_000002.1', arc), ('GCA_000001.1', bac)],
                       taxonomy={'GCA_000002.1': 'ar', 'GCA_000001.1': 'bac'})

        self.assertEqual(self.searched_as(arc), S.RNA_ARCHAEA)
        self.assertEqual(self.searched_as(bac), S.RNA_BACTERIA)

    def test_a_genome_of_no_known_domain_is_warned_about(self):
        unknown = self.genome_dir('GCA_000009.1')

        with self.assertLogs('timestamp', level='WARNING') as captured:
            self.run_silva([('GCA_000009.1', unknown)])

        self.assertTrue(any('no domain in either' in record.getMessage()
                            for record in captured.records), captured.records)


# ---------------------------------------------------------- what decides the work

class WhatDecidesTheWork(TempDirCase):
    """The canary beside a genome's own results."""

    def test_a_genome_with_a_canary_for_the_gene_is_skipped(self):
        gpath = self.genome_dir('GCA_000001.1', searched=True)

        self.run_silva([('GCA_000001.1', gpath)])

        self.assertIsNone(self.searched_as(gpath))

    def test_a_canary_for_another_gene_does_not_skip_this_one(self):
        gpath = self.genome_dir('GCA_000001.1', searched=True)

        self.run_silva([('GCA_000001.1', gpath)], rna_gene='lsu_23S')

        self.assertEqual(self.searched_as(gpath, 'lsu_23S'), S.RNA_BACTERIA)

    def test_all_searches_it_again(self):
        gpath = self.genome_dir('GCA_000001.1', searched=True)

        self.run_silva([('GCA_000001.1', gpath)], all_genomes=True)

        self.assertEqual(self.searched_as(gpath), S.RNA_BACTERIA)

    def test_a_genome_searched_again_keeps_nothing_its_last_search_found(self):
        """The earlier taxonomy of a gene this search no longer finds would
        otherwise go on into the metadata as though this search had made it."""
        gpath = self.genome_dir('GCA_000001.1', no_gene=True, searched=True)

        self.run_silva([('GCA_000001.1', gpath)], all_genomes=True)

        self.assertFalse(os.path.exists(os.path.join(self.results(gpath),
                                                     'ssu.taxonomy.tsv')))

    def test_another_gene_s_results_survive_a_search_without_remove(self):
        gpath = self.genome_dir('GCA_000001.1', searched=True)

        self.run_silva([('GCA_000001.1', gpath)], rna_gene='lsu_5S')

        self.assertTrue(os.path.exists(os.path.join(self.results(gpath),
                                                    'ssu' + S.CANARY_EXT)))

    def test_remove_empties_the_directory_of_a_genome_it_searches(self):
        gpath = self.genome_dir('GCA_000001.1', searched=True)

        self.run_silva([('GCA_000001.1', gpath)], rna_gene='lsu_5S',
                       remove_prior=True)

        self.assertFalse(os.path.exists(os.path.join(self.results(gpath),
                                                     'ssu' + S.CANARY_EXT)))

    def test_remove_does_not_fall_over_on_a_genome_never_searched(self):
        """It listed a directory that a new genome does not have yet."""
        gpath = self.genome_dir('GCA_000001.1')

        _, ok = self.run_silva([('GCA_000001.1', gpath)], remove_prior=True)

        self.assertTrue(ok)
        self.assertEqual(self.release_report(), [])

    def test_the_canary_is_copied_after_everything_it_vouches_for(self):
        """A copy stopped partway through must not leave a genome that looks
        finished with half its files."""
        gpath = self.genome_dir('GCA_000001.1')
        copied = []
        real_copy = shutil.copy

        def recording_copy(source, target):
            copied.append(os.path.basename(target))
            return real_copy(source, target)

        manager = self.manager(batch_size=1)
        with mock.patch.object(S.mp, 'Pool', InProcessPool), \
                mock.patch.object(S.shutil, 'copy', recording_copy):
            manager.run(self.genome_dirs_file([('GCA_000001.1', gpath)]), self.out_dir)

        # the canary reaches the genome directory once, and nothing after it
        results = [name for name in copied if name.startswith('ssu.')]
        self.assertEqual(results.index('ssu' + S.CANARY_EXT), len(results) - 1)
        self.assertIn('ssu.taxonomy.tsv', results)


class InProcessPool(object):
    """Runs the pool's work in this process, so that a patch in the test holds
    inside it."""

    sizes = None

    def __init__(self, processes=None):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def imap_unordered(self, func, items):
        for item in items:
            if InProcessPool.sizes is not None:
                InProcessPool.sizes.append(len(pickle.dumps((func, item))))
            yield func(item)


class WhatThePoolIsHanded(TempDirCase):
    """Each genome, and not the release's domain table with it."""

    MOST_A_TASK_SHOULD_WEIGH = 10000

    def test_neither_pass_hands_a_worker_the_domain_table(self):
        genomes = [('GCA_00000{}.1'.format(i),
                    self.genome_dir('GCA_00000{}.1'.format(i))) for i in (1, 2)]
        manager = self.manager()
        manager.domains = {'GCA_{:09d}.1'.format(i): C.DOMAIN_BACTERIA
                           for i in range(100000)}

        InProcessPool.sizes = []
        self.addCleanup(setattr, InProcessPool, 'sizes', None)
        with mock.patch.object(S.mp, 'Pool', InProcessPool):
            ok = manager.run(self.genome_dirs_file(genomes), self.out_dir)

        self.assertTrue(ok)
        self.assertEqual(len(InProcessPool.sizes), 4)
        self.assertLess(max(InProcessPool.sizes), self.MOST_A_TASK_SHOULD_WEIGH)


# --------------------------------------------------------- genomes that cannot be searched

class AGenomeThatCannotBeSearched(TempDirCase):
    """Named and left, rather than taking the other ten thousand of the batch."""

    def test_a_genome_with_no_sequences_is_named(self):
        whole = self.genome_dir('GCA_000001.1')
        gone = self.genome_dir('GCA_000003.1', fasta=False)

        _, ok = self.run_silva([('GCA_000001.1', whole), ('GCA_000003.1', gone)],
                               all_genomes=True)

        self.assertTrue(ok)
        self.assertEqual(self.release_report(),
                         [('GCA_000003.1', S.REASON_NO_GENOMIC_FASTA)])
        self.assertEqual(self.searched_as(whole), S.RNA_BACTERIA)

    def test_a_genome_the_search_fails_on_is_named_and_the_batch_succeeds(self):
        whole = self.genome_dir('GCA_000001.1')
        bad = self.genome_dir('GCA_000004.1', refuse=True)

        _, ok = self.run_silva([('GCA_000001.1', whole), ('GCA_000004.1', bad)])

        self.assertTrue(ok)
        self.assertEqual(self.release_report(),
                         [('GCA_000004.1', S.REASON_SEARCH_FAILED)])
        self.assertFalse(os.path.exists(os.path.join(self.results(bad),
                                                     'ssu' + S.CANARY_EXT)))

    def test_a_batch_in_which_every_genome_failed_is_failed(self):
        bad = self.genome_dir('GCA_000004.1', refuse=True)
        worse = self.genome_dir('GCA_000005.1', refuse=True)

        _, ok = self.run_silva([('GCA_000004.1', bad), ('GCA_000005.1', worse)])

        self.assertFalse(ok)
        self.assertEqual(B.batch_state(self.batches()[0]), B.STATE_FAILED)

    def test_a_genome_with_no_rrna_gene_is_searched_not_failed(self):
        gpath = self.genome_dir('GCA_000001.1', no_gene=True)

        _, ok = self.run_silva([('GCA_000001.1', gpath)])

        self.assertTrue(ok)
        self.assertEqual(self.release_report(), [])
        self.assertTrue(os.path.exists(os.path.join(self.results(gpath),
                                                    'ssu' + S.CANARY_EXT)))


# -------------------------------------------------- what the results were made with

class TheVersionsThatSearchedAGenome(TempDirCase):

    def version(self, gpath, program):
        path = C.version_file(self.results(gpath), program)
        if not os.path.exists(path):
            return None
        with open(path) as handle:
            return handle.read().strip()

    def test_a_searched_genome_records_nhmmer_and_blastn(self):
        gpath = self.genome_dir('GCA_000001.1')

        self.run_silva([('GCA_000001.1', gpath)])

        self.assertEqual(self.version(gpath, S.NHMMER), NHMMER_VERSION)
        self.assertEqual(self.version(gpath, S.BLASTN), BLASTN_VERSION)

    def test_blastn_is_recorded_only_where_it_ran(self):
        """5S is not classified, so blastn made nothing for it."""
        gpath = self.genome_dir('GCA_000001.1')

        self.run_silva([('GCA_000001.1', gpath)], rna_gene='lsu_5S')

        self.assertEqual(self.version(gpath, S.NHMMER), NHMMER_VERSION)
        self.assertIsNone(self.version(gpath, S.BLASTN))

    def test_each_version_is_asked_once_for_the_run(self):
        genomes = [('GCA_00000{}.1'.format(i),
                    self.genome_dir('GCA_00000{}.1'.format(i))) for i in (1, 2, 3)]

        self.run_silva(genomes)

        self.assertEqual(sorted(call.args[0] for call in self.asked.call_args_list),
                         [S.BLASTN, S.NHMMER])


# --------------------------------------------------------- batches and the release

class TheBatchesAndTheRelease(TempDirCase):
    """--out_dir holds the state of the run; the rRNA genes go to the genomes."""

    def setUp(self):
        super().setUp()
        self.genomes = [(accession, self.genome_dir(accession))
                        for accession in ('GCA_000001.1', 'GCA_000002.1',
                                          'GCA_000003.1', 'GCA_000004.1')]

    def test_the_release_is_cut_into_batches_of_batch_size(self):
        self.run_silva(self.genomes, batch_size=2)

        self.assertEqual(len(self.batches()), 2)

    def test_each_gene_has_batches_of_its_own(self):
        """An ssu SUCCESS must not tell an lsu_23S run it has nothing to do."""
        self.run_silva(self.genomes, batch_size=2)
        self.run_silva(self.genomes, rna_gene='lsu_23S', batch_size=2)

        self.assertEqual(len(self.batches('lsu_23S')), 2)
        for _, gpath in self.genomes:
            self.assertEqual(self.searched_as(gpath, 'lsu_23S'), S.RNA_BACTERIA)

    def test_each_batch_records_what_it_came_to(self):
        self.run_silva(self.genomes, batch_size=2)

        canary = B.read_canary(os.path.join(self.batches()[0], B.SUCCESS_CANARY))
        self.assertEqual(canary['searched'], '2')

    def test_no_results_are_written_to_the_out_dir(self):
        self.run_silva(self.genomes, batch_size=2)

        top = sorted(name for name in os.listdir(self.run_dir())
                     if not name.startswith(B.BATCH_DIR_PREFIX))
        self.assertEqual(top, [S.NOT_SEARCHED_RELEASE_NAME])

    def test_a_finished_batch_is_skipped_by_a_later_run(self):
        self.run_silva(self.genomes, batch_size=2)
        shutil.rmtree(self.results(self.genomes[0][1]))

        self.run_silva(self.genomes, batch_size=2)

        self.assertIsNone(self.searched_as(self.genomes[0][1]))

    def test_a_batch_held_by_another_machine_is_left_alone(self):
        self.run_silva(self.genomes, batch_size=2)
        batch = self.batches()[0]
        os.remove(os.path.join(batch, B.SUCCESS_CANARY))
        with open(os.path.join(batch, B.RUNNING_CANARY), 'w') as handle:
            handle.write(B.canary_payload(host='another-machine', pid='1'))

        self.run_silva(self.genomes, batch_size=2)

        self.assertEqual(B.batch_state(batch), B.STATE_RUNNING)


class UpdateSilva(TempDirCase):

    def test_it_needs_no_run_of_rna_silva_set_up(self):
        """It was called on an instance built with no SILVA directory, whose
        constructor raised joining None to a path."""
        ssu = os.path.join(self.dir, 'ssu.fasta')
        lsu = os.path.join(self.dir, 'lsu.fasta')
        for path, seq_id in ((ssu, 'A1'), (lsu, 'B1')):
            with open(path, 'w') as handle:
                handle.write('>{} Bacteria;Firmicutes\nACGT\n'.format(seq_id))

        S.RnaManagerSILVA.update_silva(ssu, lsu, self.dir)

        with open(os.path.join(self.dir, 'silva_taxonomy.ssu.tsv')) as handle:
            self.assertEqual(handle.read(), 'A1\tBacteria;Firmicutes\n')


if __name__ == '__main__':
    unittest.main()
