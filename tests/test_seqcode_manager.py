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

"""Offline unit tests for seqcode_manager.py -- seqcode download, which was
generate_seqcode_table of utils/tools.py.

That printed every species' record as JSON, read every genome's assembly report
on worker processes and passed over one that failed, ended with a KeyError on a
WGS prefix the release did not hold, never found a type given as an AP accession,
and wrote 'None' for a value the Registry did not give. The Registry is stood in
for by a dictionary of URLs, and NCBI's E-utilities by FakeNCBI.
"""

import contextlib
import io
import logging
import os
import shutil
import tempfile
import unittest
from unittest import mock

import requests

from gtdb_migration_tk import __main__ as main_module
from gtdb_migration_tk import main as main_py
from gtdb_migration_tk import seqcode_manager as S

SUMMARY_HEADER = '#assembly_accession\tbioproject\twgs_master\tassembly_level\n'


def taxon(rank, name, taxon_id, type_id, status='Valid (SeqCode)'):
    return {'id': taxon_id, 'rank': rank, 'name': name, 'status_name': status,
            'nomenclatural_type': {'class': 'Name', 'id': type_id}}


def species(species_id, name, material, genus_id=20, type_of_genus=None):
    """A list entry of type-genomes.json: genus 20 is the type of the family."""
    return {'id': species_id, 'name': name, 'rank': 'species', 'status_name': 'Valid (SeqCode)',
            'priority_date': None, 'nomenclatural_type': material,
            'classification': [taxon('domain', 'Bacteria', 1, 999, 'Valid (ICNP)'),
                               taxon('family', 'Exampleaceae', 10, genus_id),
                               taxon('genus', 'Examplus', genus_id, type_of_genus or species_id)],
            'created_at': '2025-01-01', 'updated_at': '2026-01-01',
            'url': 'https://api.seqco.de/v1/names/{}.json'.format(species_id)}


class FakeNCBI(object):
    """esummary and elink of E-utilities, over sequences {accession.version: assemblies},
    each assembly (GenBank accession, RefSeq accession or ''); an accession without
    its version is the sequence's latest, as NCBI takes it."""

    def __init__(self, sequences):
        self.requests = []
        self.nuccore, self.assembly, links = {}, {}, {}
        for number, (accession, assemblies) in enumerate(sorted(sequences.items()), start=1000):
            uids = []
            for genbank, refseq in assemblies:
                uid = str(5000 + len(self.assembly))
                self.assembly[uid] = {'uid': uid, 'assemblyaccession': refseq or genbank,
                                      'synonym': {'genbank': genbank, 'refseq': refseq}}
                uids.append(uid)
            self.nuccore[str(number)] = {'uid': str(number), 'caption': accession.split('.')[0],
                                         'accessionversion': accession, 'links': uids}

    def __call__(self, url, params):
        params = list(params)
        self.requests.append((url, params))
        value = dict(params)
        if url.endswith('/elink.fcgi'):
            assert (value['dbfrom'], value['db']) == ('nuccore', 'assembly'), params
            linksets = []
            for uid in [v for k, v in params if k == 'id']:
                linkset = {'dbfrom': 'nuccore', 'ids': [uid]}
                if self.nuccore[uid]['links']:
                    linkset['linksetdbs'] = [{'dbto': 'assembly', 'linkname': 'nuccore_assembly',
                                              'links': self.nuccore[uid]['links']}]
                linksets.append(linkset)
            return {'header': {'type': 'elink'}, 'linksets': linksets}
        assert url.endswith('/esummary.fcgi'), url
        result, unknown = {'uids': []}, []
        for asked in value['id'].split(','):
            if value['db'] == 'assembly':
                document = self.assembly[asked]
            else:
                held = [d for d in self.nuccore.values() if asked in (d['accessionversion'], d['caption'])]
                if not held:
                    unknown.append(asked)
                    continue
                document = {k: v for k, v in max(held, key=lambda d: d['accessionversion']).items() if k != 'links'}
            result['uids'].append(document['uid'])
            result[document['uid']] = document
        answer = {'header': {'type': 'esummary'}, 'result': result}
        if unknown:
            answer['error'] = 'Invalid uid {} at position= 0'.format(unknown[0])
        return answer

    def asked(self, tool):
        return [params for url, params in self.requests if url.endswith('/{}.fcgi'.format(tool))]


class SeqCodeCase(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='seqcode_test.')
        self.addCleanup(shutil.rmtree, self.dir, True)
        logging.getLogger('timestamp').addHandler(logging.NullHandler())

    def genome(self, accession):
        """A genome's directory, which is not made: nothing in it is to be read."""
        return os.path.join(self.dir, 'genomes', accession + '_ASM1v1')

    def summaries(self, *rows):
        """An assembly summary of rows (accession, wgs_master)."""
        path = os.path.join(self.dir, 'assembly_summary_bacteria_genbank.txt')
        with open(path, 'w') as handle:
            handle.write(SUMMARY_HEADER)
            for accession, wgs_master in rows:
                handle.write('{}\tPRJNA1\t{}\tContig\n'.format(accession, wgs_master))
        return [path]

    def genome_dirs(self, *genomes):
        path = os.path.join(self.dir, 'genome_dirs.tsv')
        with open(path, 'w') as handle:
            for accession, gdir in genomes:
                handle.write('{}\t{}\tG{}\n'.format(accession, gdir, accession[4:13]))
        return path


class FetchingTheRegistry(SeqCodeCase):
    def response(self, document=None, status=200):
        answer = mock.Mock()
        answer.json.return_value = document
        answer.raise_for_status.side_effect = (requests.HTTPError('{} Server Error'.format(status))
                                               if status >= 400 else None)
        return answer

    def test_a_request_not_answered_is_tried_again_waiting_longer_each_time(self):
        sleep = mock.Mock()
        with mock.patch.object(S.requests, 'get', side_effect=[requests.ConnectTimeout('timed out'),
                                                               self.response(status=503),
                                                               self.response({'values': []})]) as get:
            self.assertEqual(S.fetch_json('https://x/a.json', sleep=sleep), {'values': []})
        self.assertEqual(get.call_count, 3)
        self.assertEqual([c.args[0] for c in sleep.call_args_list], [S.FETCH_WAIT, S.FETCH_WAIT * 2])

    def test_a_request_failing_every_attempt_is_an_error_naming_the_url(self):
        with mock.patch.object(S.requests, 'get', side_effect=requests.ConnectTimeout('timed out')), \
                self.assertRaisesRegex(S.SeqCodeError, r'https://x/a.json could not be read in 4 attempt'):
            S.fetch_json('https://x/a.json', sleep=mock.Mock())

    def test_the_parameters_of_a_request_are_passed_on(self):
        with mock.patch.object(S.requests, 'get', return_value=self.response({'result': {}})) as get:
            S.fetch_json('https://x/esummary.fcgi', params=[('id', '1'), ('id', '2')])
        self.assertEqual(get.call_args.kwargs['params'], [('id', '1'), ('id', '2')])

    def test_the_registry_is_api_seqco_de(self):
        # disc-genomics.uibk.ac.at no longer answers
        self.assertEqual(S.TYPE_GENOMES_URL, 'https://api.seqco.de/v1/type-genomes.json')


class FindingTheTypeGenome(SeqCodeCase):
    def ncbi(self):
        return FakeNCBI({'CP030759.1': [('GCA_000000002.1', 'GCF_000000002.1')],
                         'CP030759.2': [('GCA_000000012.1', '')],
                         'AP043677.1': [('GCA_000000003.1', '')],
                         'CP000004.1': []})

    def test_ncbi_gives_the_assemblies_of_a_complete_sequence_versioned_or_not_and_of_ddbj(self):
        ncbi = self.ncbi()
        with mock.patch.object(S, 'fetch_json', side_effect=ncbi):
            found = S.sequence_assemblies(['CP030759.1', 'CP030759', 'AP043677', 'CP000004.1', 'CP999999.1'])
        self.assertEqual(found, {'CP030759.1': ('GCA_000000002.1', 'GCF_000000002.1'),
                                 'CP030759': ('GCA_000000012.1',),          # unversioned is the latest
                                 'AP043677': ('GCA_000000003.1',)})         # AP types were never found
        # one request of each kind, naming the caller, where the reports of every genome were read
        self.assertEqual([url.rsplit('/', 1)[1] for url, _ in ncbi.requests],
                         ['esummary.fcgi', 'elink.fcgi', 'esummary.fcgi'])
        for _, params in ncbi.requests:
            self.assertIn(('tool', S.EUTILS_TOOL), params)
            self.assertIn(('retmode', 'json'), params)

    def test_the_sequences_are_asked_of_ncbi_a_batch_at_a_time(self):
        ncbi = self.ncbi()
        with mock.patch.object(S, 'fetch_json', side_effect=ncbi), mock.patch.object(S, 'EUTILS_BATCH', 2):
            found = S.sequence_assemblies(['CP030759.1', 'CP030759.2', 'AP043677.1'])
        self.assertEqual(sorted(found), ['AP043677.1', 'CP030759.1', 'CP030759.2'])
        self.assertEqual(len(ncbi.asked('elink')), 2)
        for params in ncbi.requests:
            ids = [v for k, v in params[1] if k == 'id']
            self.assertLessEqual(sum(len(v.split(',')) for v in ids), 2)

    def test_nothing_is_asked_of_ncbi_for_no_sequences(self):
        with mock.patch.object(S, 'fetch_json') as fetch:
            self.assertEqual(S.sequence_assemblies([]), {})
        fetch.assert_not_called()

    def test_an_answer_without_its_result_is_an_error_naming_the_url(self):
        # a type passed over would leave its species out of the table unsaid
        for answer, tool in (({'esummaryresult': ['Unable to obtain query #1']}, 'esummary'),
                             ({'ERROR': 'Backend failed'}, 'elink')):
            ncbi = self.ncbi()

            def fetch(url, params):
                return answer if url.endswith('/{}.fcgi'.format(tool)) else ncbi(url, params)
            with self.subTest(tool=tool), mock.patch.object(S, 'fetch_json', side_effect=fetch), \
                    self.assertRaisesRegex(S.SeqCodeError, r'{}\.fcgi answered without'.format(tool)):
                S.sequence_assemblies(['CP030759.1'])

    def test_a_type_is_found_by_assembly_sequence_unversioned_sequence_or_wgs_prefix(self):
        canonical = {'G000000001': 'GCF_000000001.1'}
        index = {'CP030759': 'GCA_000000002.1', 'AP043677': 'GCA_000000003.1', 'JBBJJH': 'GCA_000000002.1'}
        cases = [({'assembly': 'GCA_000000001.1'}, 'GCF_000000001.1'),       # RefSeq copy of a GenBank type
                 ({'nuccore': 'CP030759.1'}, 'GCA_000000002.1'),
                 ({'nuccore': 'CP030759'}, 'GCA_000000002.1'),
                 ({'nuccore': 'AP043677'}, 'GCA_000000003.1'),                # AP types were never found
                 ({'nuccore': 'JBBJJH000000000.1'}, 'GCA_000000002.1'),
                 ({'nuccore': 'JBZZZZ000000000.1'}, None),                    # a KeyError ended the run
                 ({'assembly': 'GCA_000000099.1'}, None)]
        for material, genome in cases:
            with self.subTest(material=material):
                self.assertEqual(S.type_genome(material, canonical, index)[0], genome)

    def test_a_wgs_project_is_found_from_the_summaries_wgs_master_by_canonical_accession(self):
        canonical = {'G000000002': 'GCA_000000002.1', 'G000000003': 'GCA_000000003.1'}
        index = S.wgs_index(self.summaries(('GCF_000000002.1', 'JBBJJH000000000.1'),
                                           ('GCA_000000003.1', 'na'),
                                           ('GCA_000000099.1', 'ABZX00000000.1')), canonical)
        # the RefSeq row names the GenBank genome of the release; a genome not in it is not indexed
        self.assertEqual(index, {'JBBJJH': 'GCA_000000002.1'})
        self.assertEqual(S.type_genome({'nuccore': 'JBBJJH000000000.1'}, canonical, index)[0], 'GCA_000000002.1')
        self.assertEqual(S.type_genome({'nuccore': 'JBBJJH000000000'}, canonical, index)[0], 'GCA_000000002.1')

    def test_a_summary_without_wgs_master_is_refused(self):
        from gtdb_migration_tk.ncbi_utils import BadInput
        path = os.path.join(self.dir, 'old.txt')
        with open(path, 'w') as handle:
            handle.write('#assembly_accession\tbioproject\nGCA_000000001.1\tPRJNA1\n')
        with self.assertRaisesRegex(BadInput, 'wgs_master'):
            S.wgs_index([path], {})

    def test_the_columns_are_those_metadata_seqcode_desc_tsv_describes(self):
        # update_metadata_db loads the table into metadata_seqcode by these names; it
        # was written with seqcode_proposed_in, a column metadata_seqcode does not have
        path = os.path.join(os.path.dirname(S.__file__), 'data_files', 'table_description',
                            'metadata_seqcode.desc.tsv')
        with open(path) as handle:
            described = [line.split('\t')[0] for line in handle if line.strip()]
        self.assertEqual(S.TABLE_HEADER[0], 'seqcode_type_material_accn')    # the genome
        self.assertEqual(sorted(S.TABLE_HEADER[1:]), sorted(described))
        self.assertIn('seqcode_proposed_by', S.TABLE_HEADER)

    def test_the_types_are_read_from_the_classification_and_a_missing_genus_is_none(self):
        taxonomy, statuses, types = S.classification(species(30, 'Examplus a', None)['classification'], 30)
        self.assertEqual(taxonomy, 'd__Bacteria;p__;c__;o__;f__Exampleaceae;g__Examplus;s__')
        self.assertEqual(statuses[0], 'Valid (ICNP)')
        self.assertEqual((types[4], types[5]), (True, True))
        # a classification without a genus failed with an UnboundLocalError
        self.assertFalse(any(S.classification([taxon('domain', 'Bacteria', 1, 999)], 30)[2]))


class DownloadingTheTable(SeqCodeCase):
    def registry(self, records, species_records=None, ncbi=None):
        """fetch_json over a Registry of one page per two records, and NCBI."""
        pages = [records[i:i + 2] for i in range(0, len(records), 2)] or [[]]
        urls = {S.TYPE_GENOMES_URL: {'response': {'total_pages': len(pages), 'count': len(records)},
                                     'values': pages[0]}}
        for number, page in enumerate(pages[1:], start=2):
            urls['{}?page={}'.format(S.TYPE_GENOMES_URL, number)] = {'values': page}
        for record in records:
            urls[record['url']] = (species_records or {}).get(record['id'], {
                'proposed_in': {'citation': 'Doe et al., 2025, Microbial Genomics'}})
        self.fetched = []
        ncbi = ncbi or FakeNCBI({})

        def fetch(url, params=None):
            self.fetched.append(url)
            if url.startswith(S.EUTILS):
                return ncbi(url, params)
            if isinstance(urls.get(url), Exception):
                raise urls[url]
            return urls[url]
        return fetch

    def run_command(self, genome_dirs, fetch, summaries=None, cache_file=None):
        console = io.StringIO()
        manager = S.SeqCodeManager(fetchers=3)
        with mock.patch.object(S, 'fetch_json', side_effect=fetch), contextlib.redirect_stdout(console), \
                contextlib.redirect_stderr(console), self.assertLogs('timestamp', level='INFO') as logged:
            manager.run(genome_dirs, self.dir, summaries or self.summaries(), cache_file)
        return console.getvalue(), [r.getMessage() for r in logged.records]

    def table(self, name=S.TABLE_NAME):
        with open(os.path.join(self.dir, name)) as handle:
            return [line.split('\t') for line in handle.read().splitlines()]

    def test_each_species_typed_by_a_genome_of_the_release_is_a_row_and_nothing_is_printed(self):
        # the genome directories do not exist: nothing in them is read
        g1, g2 = self.genome('GCA_000000001.1'), self.genome('GCA_000000002.1')
        genome_dirs = self.genome_dirs(('GCA_000000001.1', g1), ('GCA_000000002.1', g2))
        records = [species(30, 'Examplus a', {'assembly': 'GCF_000000001.1'}),
                   species(31, 'Examplus b', {'nuccore': 'AP000002'}, type_of_genus=30),
                   species(32, 'Examplus c', {'assembly': 'GCA_000000099.1'})]
        # NCBI gives the RefSeq accession first; the release holds the GenBank genome
        ncbi = FakeNCBI({'AP000002.1': [('GCA_000000002.1', 'GCF_000000002.1')]})
        console, messages = self.run_command(genome_dirs, self.registry(records, ncbi=ncbi))

        self.assertEqual(console, '')
        rows = self.table()
        self.assertEqual(rows[0], list(S.TABLE_HEADER))
        first = dict(zip(rows[0], rows[1]))
        self.assertEqual((first['seqcode_type_material_accn'], first['seqcode_name'], first['seqcode_proposed_by'],
                          first['seqcode_type_species_of_genus'], first['seqcode_priority_date']),
                         ('GCA_000000001.1', 'Examplus a', 'Doe et al., 2025, Microbial Genomics', 'True', ''))
        self.assertEqual(dict(zip(rows[0], rows[2]))['seqcode_type_material_accn'], 'GCA_000000002.1')
        self.assertEqual(self.table(S.NOT_IN_RELEASE_NAME)[1:], [['32', 'Examplus c', 'GCA_000000099.1']])
        # the list is read once, and a species' record only where its type is in the release
        self.assertNotIn(records[2]['url'], self.fetched)
        self.assertEqual(self.fetched.count(S.TYPE_GENOMES_URL), 1)
        self.assertIn('Wrote 2 species to {}.'.format(os.path.join(self.dir, S.TABLE_NAME)), messages)

    def test_ncbi_is_not_asked_where_every_type_is_an_assembly_or_a_wgs_project(self):
        g1 = self.genome('GCA_000000001.1')
        g2 = self.genome('GCA_000000002.1')
        genome_dirs = self.genome_dirs(('GCA_000000001.1', g1), ('GCA_000000002.1', g2))
        records = [species(30, 'Examplus a', {'assembly': 'GCA_000000001.1'}),
                   species(31, 'Examplus b', {'nuccore': 'JBBJJH000000000.1'})]
        self.run_command(genome_dirs, self.registry(records),
                         self.summaries(('GCA_000000002.1', 'JBBJJH000000000.1')))
        self.assertEqual([url for url in self.fetched if url.startswith(S.EUTILS)], [])
        self.assertEqual([row[0] for row in self.table()[1:]], ['GCA_000000001.1', 'GCA_000000002.1'])

    def test_a_sequence_ncbi_places_in_no_genome_of_the_release_is_listed_not_in_the_release(self):
        genome_dirs = self.genome_dirs(('GCA_000000001.1', self.genome('GCA_000000001.1')))
        records = [species(30, 'Examplus a', {'nuccore': 'CP000001.1'}),
                   species(31, 'Examplus b', {'nuccore': 'CP000002.1'}),      # in an assembly not in the release
                   species(32, 'Examplus c', {'nuccore': 'CP000003'}),        # in no assembly
                   species(33, 'Examplus d', {'nuccore': 'CP000099.1'})]      # not known to NCBI
        ncbi = FakeNCBI({'CP000001.1': [('GCA_000000001.1', '')], 'CP000002.1': [('GCA_000000002.1', '')],
                         'CP000003.1': []})
        _, messages = self.run_command(genome_dirs, self.registry(records, ncbi=ncbi))
        self.assertEqual([row[0] for row in self.table()[1:]], ['GCA_000000001.1'])
        self.assertEqual(self.table(S.NOT_IN_RELEASE_NAME)[1:], [['31', 'Examplus b', 'CP000002.1'],
                                                                 ['32', 'Examplus c', 'CP000003'],
                                                                 ['33', 'Examplus d', 'CP000099.1']])
        self.assertTrue(any(m.startswith('NCBI gives the assembly of 2 of the 4 complete sequence(s), 1 of them '
                                         'in the release') for m in messages))

    def test_an_ncbi_answer_without_its_result_ends_the_run_with_nothing_written(self):
        genome_dirs = self.genome_dirs(('GCA_000000001.1', self.genome('GCA_000000001.1')))
        fetch = self.registry([species(30, 'Examplus a', {'nuccore': 'CP000001.1'})])

        def failing(url, params=None):
            return {'esummaryresult': ['Unable to obtain query #1']} if url.startswith(S.EUTILS) else fetch(url)
        with self.assertRaisesRegex(S.SeqCodeError, 'esummary.fcgi answered without its result'):
            self.run_command(genome_dirs, failing)
        self.assertEqual([n for n in os.listdir(self.dir) if n.startswith('seqcode_table')], [])

    def test_a_species_unchanged_since_the_cache_is_not_fetched_again(self):
        g1 = self.genome('GCA_000000001.1')
        genome_dirs = self.genome_dirs(('GCA_000000001.1', g1))
        records = [species(30, 'Examplus a', {'assembly': 'GCA_000000001.1'})]
        self.run_command(genome_dirs, self.registry(records))
        self.assertIn(records[0]['url'], self.fetched)

        self.run_command(genome_dirs, self.registry(records))
        self.assertNotIn(records[0]['url'], self.fetched)
        self.assertEqual(dict(zip(*self.table()))['seqcode_proposed_by'], 'Doe et al., 2025, Microbial Genomics')

        records[0]['updated_at'] = '2026-10-01'
        self.run_command(genome_dirs, self.registry(records))
        self.assertIn(records[0]['url'], self.fetched)

    def test_a_cache_given_is_the_one_read_and_written(self):
        g1 = self.genome('GCA_000000001.1')
        genome_dirs = self.genome_dirs(('GCA_000000001.1', g1))
        cache = os.path.join(self.dir, 'shared_cache.json')
        self.run_command(genome_dirs, self.registry([species(30, 'Examplus a', {'assembly': 'GCA_000000001.1'})]),
                         cache_file=cache)
        self.assertEqual(S.read_cache(cache)['30']['proposed_in'], 'Doe et al., 2025, Microbial Genomics')
        self.assertFalse(os.path.exists(os.path.join(self.dir, S.CACHE_NAME)))

    def test_species_are_fetched_several_at_a_time(self):
        import threading
        genomes = [('GCA_00000000{}.1'.format(i), self.genome('GCA_00000000{}.1'.format(i))) for i in range(1, 7)]
        genome_dirs = self.genome_dirs(*genomes)
        records = [species(30 + i, 'Examplus {}'.format(i), {'assembly': accession})
                   for i, (accession, _) in enumerate(genomes)]
        registry = self.registry(records)
        at_once, lock, peak = [0], threading.Lock(), [0]
        three_started = threading.Barrier(3, timeout=5)

        def fetch(url, params=None):
            if url.startswith(S.TYPE_GENOMES_URL):
                return registry(url)
            with lock:
                at_once[0] += 1
                peak[0] = max(peak[0], at_once[0])
            try:
                three_started.wait()
            except threading.BrokenBarrierError:
                pass
            with lock:
                at_once[0] -= 1
            return registry(url)
        self.run_command(genome_dirs, fetch)
        self.assertEqual(peak[0], 3)
        self.assertEqual(len(self.table()), 7)

    def test_a_request_that_fails_ends_the_run_with_nothing_written(self):
        g1 = self.genome('GCA_000000001.1')
        genome_dirs = self.genome_dirs(('GCA_000000001.1', g1))
        records = [species(30, 'Examplus a', {'assembly': 'GCA_000000001.1'})]
        fetch = self.registry(records, {30: S.SeqCodeError('https://api.seqco.de/v1/names/30.json could not be read')})
        with self.assertRaisesRegex(S.SeqCodeError, 'names/30.json'):
            self.run_command(genome_dirs, fetch)
        self.assertEqual([n for n in os.listdir(self.dir) if n.startswith('seqcode_table')], [])

    def test_a_run_stopped_by_a_failed_request_keeps_what_it_fetched(self):
        g1, g2 = self.genome('GCA_000000001.1'), self.genome('GCA_000000002.1')
        genome_dirs = self.genome_dirs(('GCA_000000001.1', g1), ('GCA_000000002.1', g2))
        records = [species(30, 'Examplus a', {'assembly': 'GCA_000000001.1'}),
                   species(31, 'Examplus b', {'assembly': 'GCA_000000002.1'})]
        fetch = self.registry(records, {31: S.SeqCodeError('names/31.json could not be read')})
        with self.assertRaises(S.SeqCodeError):
            self.run_command(genome_dirs, fetch)
        self.assertEqual(sorted(S.read_cache(os.path.join(self.dir, S.CACHE_NAME))), ['30'])


class TheCommandLine(SeqCodeCase):
    def argv(self, *extra):
        return ['seqcode', 'download', '-g', 'genome_dirs.tsv', '-n', 'a.txt', 'b.txt',
                '-o', os.path.join(self.dir, 'out')] + list(extra)

    def test_the_summaries_are_required_and_the_command_is_handed_its_options(self):
        with mock.patch('sys.stderr'), self.assertRaises(SystemExit) as ended:
            main_module.get_main_parser().parse_args(['seqcode', 'download', '-g', 'genome_dirs.tsv', '-o', 'out'])
        self.assertEqual(ended.exception.code, 2)
        self.assertIsNone(main_module.get_main_parser().parse_args(self.argv()).log)

        options = main_module.get_main_parser().parse_args(self.argv('-l', 'run.log', '--species_cache', 'cache.json'))
        with mock.patch.object(main_py, 'SeqCodeManager') as manager, mock.patch.object(main_py, 'check_file_exists'):
            main_py.OptionsParser().parse_options(options)
        manager.assert_called_once_with()
        manager.return_value.run.assert_called_once_with('genome_dirs.tsv', os.path.join(self.dir, 'out'),
                                                         ['a.txt', 'b.txt'], 'cache.json')

    def run_main(self, *extra):
        """main() over the command, its manager not run."""
        for name in ('timestamp', 'no_timestamp'):
            self.addCleanup(self.drop_handlers, logging.getLogger(name))
        genome_dirs = os.path.join(self.dir, 'genome_dirs.tsv')
        summary = os.path.join(self.dir, 'a.txt')
        for path in (genome_dirs, summary):
            open(path, 'w').close()
        argv = ['gtdb_migration_tk', 'seqcode', 'download', '-g', genome_dirs, '-n', summary,
                '-o', os.path.join(self.dir, 'out'), '--silent'] + list(extra)
        with mock.patch('sys.argv', argv), mock.patch.object(main_py, 'SeqCodeManager'), \
                contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            main_module.main()

    @staticmethod
    def drop_handlers(logger):
        for handler in list(logger.handlers):
            logger.removeHandler(handler)
            handler.close()

    def test_the_log_is_gtdb_migration_tk_log_in_the_out_dir_where_l_is_not_given(self):
        self.run_main()
        self.assertTrue(os.path.isfile(os.path.join(self.dir, 'out', main_module.FALLBACK_LOG)))

    def test_the_log_is_where_l_says_where_it_is_given(self):
        log = os.path.join(self.dir, 'run.log')
        self.run_main('-l', log)
        self.assertTrue(os.path.isfile(log))
        self.assertFalse(os.path.exists(os.path.join(self.dir, 'out', main_module.FALLBACK_LOG)))

    def test_a_registry_that_does_not_answer_exits_1(self):
        options = main_module.get_main_parser().parse_args(self.argv('-l', 'run.log'))
        with mock.patch.object(main_py, 'SeqCodeManager') as manager, mock.patch.object(main_py, 'check_file_exists'), \
                self.assertLogs('timestamp', level='ERROR'), self.assertRaises(SystemExit) as ended:
            manager.return_value.run.side_effect = S.SeqCodeError('no answer')
            main_py.OptionsParser().parse_options(options)
        self.assertEqual(ended.exception.code, 1)

    def test_cpus_is_no_longer_an_option(self):
        # it only read the assembly reports, which NCBI's E-utilities replace
        with mock.patch('sys.stderr'), self.assertRaises(SystemExit) as ended:
            main_module.get_main_parser().parse_args(self.argv('-l', 'run.log', '-c', '8'))
        self.assertEqual(ended.exception.code, 2)

    def test_generate_seqcode_table_and_download_seqcode_data_are_no_longer_commands(self):
        for old_name in ('generate_seqcode_table', 'download_seqcode_data'):
            with mock.patch('sys.stderr'), self.assertRaises(SystemExit) as ended:
                main_module.get_main_parser().parse_args([old_name, '-g', 'genome_dirs.tsv', '-n', 'a.txt',
                                                          '-o', 'out', '-l', 'run.log'])
            self.assertEqual(ended.exception.code, 2)

    def test_seqcode_alone_is_refused_before_a_log_is_opened(self):
        # lpsn, given no step, starts a log and ends in a TypeError
        with mock.patch('sys.stderr'), self.assertRaises(SystemExit) as ended:
            main_module.get_main_parser().parse_args(['seqcode'])
        self.assertEqual(ended.exception.code, 2)


if __name__ == '__main__':
    unittest.main()
