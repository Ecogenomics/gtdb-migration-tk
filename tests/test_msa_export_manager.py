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

"""Offline unit tests for msa_export_manager.py: export_msa and list_marker_sets.

export_msa writes the full MSA of every genome it is given, unfiltered and
untrimmed, concatenated as 'gtdb tree create' concatenated each genome's
markers, with the ARB metadata gtdb wrote. What would break silently is a genome
dropped or a column shifted, an ARB record whose fields no longer line up with
its table, or an export made from a genome that was never aligned. The database
is stood in for by a manager whose reads return a small database held in the
test: what is tested is what is decided from them, and what is written.
"""

import gzip
import logging
import os
import shutil
import tempfile
import unittest
from unittest import mock

from gtdb_migration_tk import __main__ as main_module
from gtdb_migration_tk import msa_export_manager as M


def marker(db_id, size, prefix='PFAM', id_in_database=None):
    id_in_database = id_in_database or 'PF{:05d}.1'.format(db_id)
    return M.Marker(db_id, '{}_{}'.format(prefix, id_in_database), id_in_database, 'marker{}'.format(db_id),
                    'Marker {}'.format(db_id), size)


class ConcatenatingAGenome(unittest.TestCase):
    MARKERS = [marker(1, 4), marker(2, 3), marker(3, 2), marker(4, 3)]

    def rows(self):
        return {1: M.AlignedRow('AC-D', False, '1e-50', 1),          # single copy
                2: M.AlignedRow('EFG', True, '1e-40', 1),            # several genes, one sequence
                3: M.AlignedRow('HI', True, '1e-30', 2),             # several different genes
                4: M.AlignedRow('---', False, None, None)}           # no gene found

    def test_markers_hit_by_different_genes_or_none_are_gaps(self):
        genome = M.concatenate(self.rows(), self.MARKERS)
        self.assertEqual(genome.sequence, 'AC-DEFG-----')
        self.assertEqual(genome.statuses, [M.SINGLE, M.MULTI_IDENTICAL, M.MULTIPLE, M.MISSING])

    def test_a_sequence_the_wrong_length_is_refused_rather_than_shifting_the_columns(self):
        rows = self.rows()
        rows[1] = M.AlignedRow('ACD', False, '1e-50', 1)
        with self.assertRaises(M.MSAExportError):
            M.concatenate(rows, self.MARKERS)


class ThePartitions(unittest.TestCase):

    def test_each_marker_is_its_run_of_columns_one_based_inclusive_covering_the_msa(self):
        markers = [marker(1, 121, 'PFAM', 'PF00380.20'), marker(2, 126, 'PFAM', 'PF00410.20'),
                   marker(3, 1, 'TIGR', 'TIGR00001')]
        self.assertEqual(M.gene_partitions(markers),
                         'AA, PFAM_PF00380.20 = 1-121\n'
                         'AA, PFAM_PF00410.20 = 122-247\n'
                         'AA, TIGR_TIGR00001 = 248-248\n')

    def test_a_marker_named_twice_or_unnameable_is_refused_as_locusmasker_would_refuse_it(self):
        with self.assertRaises(M.MSAExportError):
            M.gene_partitions([marker(1, 10, 'PFAM', 'PF1'), marker(2, 10, 'PFAM', 'PF1')])
        for bad in ('PF,1', 'PF=1', 'PF#1'):
            with self.subTest(name=bad), self.assertRaises(M.MSAExportError):
                M.gene_partitions([marker(1, 10, 'PFAM', bad)])


class NamingGenomes(unittest.TestCase):
    HELD = {('GB', 'GCA_000000001.1'): 1, ('RS', 'GCF_000000002.1'): 2,
            ('GB', 'X_1'): 3, ('RS', 'X_1'): 4}

    def test_a_genome_is_found_with_or_without_its_prefix(self):
        found, unknown = M.resolve_genomes(['GB_GCA_000000001.1', 'GCF_000000002.1'], self.HELD, ['GB', 'RS'])
        self.assertEqual(found, {'GB_GCA_000000001.1': 1, 'GCF_000000002.1': 2})
        self.assertEqual(unknown, [])

    def test_a_wrong_prefix_an_unheld_genome_and_an_ambiguous_one_are_unknown(self):
        found, unknown = M.resolve_genomes(['RS_GCA_000000001.1', 'GCA_000000009.1', 'X_1'],
                                           self.HELD, ['GB', 'RS'])
        self.assertEqual(found, {})
        self.assertEqual(unknown, ['RS_GCA_000000001.1', 'GCA_000000009.1', 'X_1'])

    def test_marker_set_ids_are_split_on_commas_given_once_or_repeated(self):
        self.assertEqual(M.flatten([[1, 19], [2], [19]]), [1, 19, 2])
        self.assertEqual(M.flatten(['a,b', ' c ,', 'a']), ['a', 'b', 'c'])
        self.assertEqual(M.flatten(None), [])

    def test_the_batchfile_is_read_as_gtdb_read_one(self):
        with tempfile.NamedTemporaryFile('w', suffix='.tsv', delete=False) as handle:
            handle.write('# genomes\nGB_GCA_000000001.1\tnote\n\n  RS_GCF_000000002.1  \nGB_GCA_000000001.1\n')
        self.addCleanup(os.remove, handle.name)
        self.assertEqual(M.read_genome_file(handle.name), ['GB_GCA_000000001.1', 'RS_GCF_000000002.1'])


class TheArbRecord(unittest.TestCase):
    COLUMNS = ['id', 'accession', 'formatted_accession', 'organism_name', 'ncbi_organism_name', 'ncbi_taxonomy',
               'gtdb_clustered_genomes', 'checkm_completeness', 'description']
    TAXONOMY = 'd__Bacteria;p__P;c__C;o__O;f__F;g__Escherichia;s__Escherichia coli'

    def values(self):
        metadata = dict(zip(self.COLUMNS, [7, 'RS_GCF_000000001.1', 'G000000001', 'GCF_000000001.1',
                                           'Escherichia coli K-12', self.TAXONOMY,
                                           'RS_GCF_000000001.1;GB_GCA_000000003.2', 98.7654, 'a=b\tc']))
        return M.genome_values(metadata, 'G000000001', 's__Escherichia coli: 2')

    def test_fields_add_ncbi_genus_and_species_after_ncbi_taxonomy_and_cluster_species_last(self):
        self.assertEqual(M.metadata_fields(self.COLUMNS),
                         ['accession', 'formatted_accession', 'organism_name', 'ncbi_organism_name', 'ncbi_taxonomy',
                          'ncbi_genus', 'ncbi_species', 'gtdb_clustered_genomes', 'checkm_completeness',
                          'description', 'gtdb_cluster_ncbi_species'])

    def test_an_ncbi_genome_is_given_ncbis_name_and_canonical_clustered_genomes(self):
        values = self.values()
        self.assertEqual(values['organism_name'], 'Escherichia coli K-12')
        self.assertEqual((values['ncbi_genus'], values['ncbi_species']), ('g__Escherichia', 's__Escherichia coli'))
        self.assertEqual(values['gtdb_clustered_genomes'], 'G000000001;G000000003')

    def test_the_alignment_is_left_blank_and_its_counts_are_of_the_msa_written(self):
        fields = M.metadata_fields(self.COLUMNS)
        record = M.arb_record('G000000001', fields, self.values(), M.msa_values('AC--D', [M.SINGLE, M.MULTIPLE]))
        self.assertTrue(record.startswith('BEGIN\ndb_name=G000000001\naccession=RS_GCF_000000001.1\n'))
        self.assertIn('\nmsa_aa_count=3\nmsa_length=5\nmultiple_homologs=1\naligned_seq=\nEND\n\n', record)
        self.assertIn('\nmsa_gene_count=1\nmsa_num_marker_genes=2\n', record)
        # as gtdb wrote them: '%.4g', '=' as '/', NULL as None
        self.assertIn('\ncheckm_completeness=98.77\n', record)
        self.assertIn('\ndescription=a/b\tc\n', record)

    def test_the_tsv_holds_the_same_fields_floats_to_three_places_equals_kept_and_null_none(self):
        fields = M.metadata_fields(self.COLUMNS)
        values = self.values()
        values['description'] = None
        header = M.tsv_header(fields).rstrip('\n').split('\t')
        row = M.tsv_row('G000000001', fields, values, M.msa_values('AC--D', [M.SINGLE])).rstrip('\n').split('\t')
        self.assertEqual(header[0], 'db_name')
        self.assertEqual(header[-5:], list(M.MSA_FIELDS))
        self.assertEqual(len(header), len(row))
        cells = dict(zip(header, row))
        self.assertEqual(cells['checkm_completeness'], '98.765')
        self.assertEqual(cells['description'], 'none')
        self.assertEqual(M.tsv_value('a=b\tc'), 'a=b c')
        self.assertEqual(cells['msa_aa_count'], '3')

    def test_cluster_species_are_counted_most_common_first_ties_in_cluster_order(self):
        species = {'A': 's__x', 'B': 's__y', 'C': 's__y', 'D': 's__z', 'E': None}
        self.assertEqual(M.cluster_ncbi_species('A;B;C;D;E;F', species), 's__y: 2; s__x: 1; s__z: 1')
        self.assertEqual(M.cluster_ncbi_species(None, species), '')

    def test_the_import_filter_names_every_field_and_reads_the_sequence_after_multiple_homologs(self):
        text = M.arb_filter(M.metadata_fields(self.COLUMNS))
        self.assertTrue(text.startswith('AUTODETECT\t"BEGIN"\n\nBEGIN\t"BEGIN*"\n\nMATCH\t"db_name\\=*"\n'
                                        '\tSRT "*\\=="\n\tWRITE "name"\n\nMATCH\t"organism_name\\=*"\n'))
        self.assertEqual(text.count('MATCH\t"organism_name'), 1)
        self.assertIn('MATCH\t"ncbi_species\\=*"', text)
        self.assertTrue(text.endswith('SEQUENCEAFTER\t"multiple_homologs*"\nSEQUENCESRT\t"*\\=="\n'
                                      'SEQUENCEEND\t"END"\n\nEND\t"END"\n'))


# ---- a database of four genomes, one of poor quality and one with almost no markers

MARKERS = [marker(11, 4, 'PFAM', 'PF00001.1'), marker(12, 3, 'TIGR', 'TIGR00001')]
VIEW = ['id', 'accession', 'formatted_accession', 'organism_name', 'ncbi_organism_name', 'ncbi_taxonomy',
        'gtdb_clustered_genomes', 'checkm_completeness', 'checkm_contamination']
TAX = 'd__Bacteria;p__P;c__C;o__O;f__F;g__G;s__G sp{}'
QUALITY = {1: (99.0, 0.5), 2: (95.0, 1.0), 3: (40.0, 30.0), 4: (97.0, 0.0)}


def genome(i, source):
    at_source = '{}_{:09d}.1'.format('GCA' if source == 'GB' else 'GCF', i)
    return {'id': i, 'prefix': source, 'at_source': at_source, 'accession': '{}_{}'.format(source, at_source),
            'formatted': 'G{:09d}'.format(i)}


class FakeDatabase(object):
    def __init__(self):
        self.genomes = {i: genome(i, 'RS' if i == 1 else 'GB') for i in range(1, 5)}
        self.sets = {1}
        self.aligned = {1: {11: M.AlignedRow('ACDE', False, '1e-9', 1), 12: M.AlignedRow('FGH', False, '1e-9', 1)},
                        2: {11: M.AlignedRow('ACDF', False, '1e-9', 1), 12: M.AlignedRow('FGH', True, '1e-9', 2)},
                        3: {11: M.AlignedRow('KLMN', False, '1e-9', 1), 12: M.AlignedRow('PQR', True, '1e-9', 1)},
                        4: {11: M.AlignedRow('A---', False, '1e-9', 1), 12: M.AlignedRow('---', False, None, None)}}


class FakeManager(M.MSAExportManager):
    def __init__(self, db, cpus=1):
        self.logger = logging.getLogger('timestamp')
        self.temp_con = mock.Mock()
        self.db = db
        self.cpus = cpus
        self.chunks_read = []

    def open_readers(self):
        pass

    def close_readers(self):
        pass

    def read_marker_sets(self, ids):
        return set(ids) & self.db.sets

    def read_markers(self, ids):
        return list(MARKERS)

    def read_source_prefixes(self):
        return ['GB', 'RS']

    def read_genomes(self, at_source):
        return [(g['id'], g['prefix'], g['at_source'], g['formatted'])
                for g in self.db.genomes.values() if g['at_source'] in at_source]

    def read_aligned_rows(self, genome_ids, marker_ids):
        self.chunks_read.append(list(genome_ids))
        return {g: dict(self.db.aligned[g]) for g in genome_ids if g in self.db.aligned}

    def read_metadata(self, genome_ids):
        metadata = {}
        for i in genome_ids:
            g = self.db.genomes[i]
            metadata[i] = dict(zip(VIEW, [i, g['accession'], g['formatted'], g['at_source'], 'Organism {}'.format(i),
                                          TAX.format(i), 'RS_GCF_000000001.1;GB_GCA_000000004.1' if i == 1 else None,
                                          QUALITY[i][0], QUALITY[i][1]]))
        return VIEW, metadata

    def read_ncbi_taxonomy(self, accessions, prefixes):
        return {g['accession']: TAX.format(i) for i, g in self.db.genomes.items() if g['accession'] in accessions}


class TempDirCase(unittest.TestCase):
    GENOMES = ['RS_GCF_000000001.1', 'GB_GCA_000000002.1', 'GCA_000000003.1', 'GB_GCA_000000004.1']

    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='msa_export_manager_test.')
        self.addCleanup(shutil.rmtree, self.dir, True)
        self.out_dir = os.path.join(self.dir, 'msa')
        self.log_dir = os.path.join(self.dir, 'logs')
        os.makedirs(self.log_dir)
        patcher = mock.patch.object(M, 'log_directory', return_value=self.log_dir)
        patcher.start()
        self.addCleanup(patcher.stop)
        logger = logging.getLogger('timestamp')
        handler = logging.NullHandler()
        logger.addHandler(handler)
        self.addCleanup(logger.removeHandler, handler)
        self.db = FakeDatabase()

    def export(self, genomes=None, marker_set_ids=(1,), cpus=1, **kwargs):
        self.manager = FakeManager(self.db, cpus)
        return self.manager.export_msa(self.GENOMES if genomes is None else genomes, list(marker_set_ids),
                                       self.out_dir, **kwargs)

    def read(self, suffix, prefix='gtdb'):
        path = os.path.join(self.out_dir, prefix + suffix)
        with (gzip.open(path, 'rt') if path.endswith('.gz') else open(path)) as handle:
            return handle.read()

    def read_bytes(self, name):
        with open(os.path.join(self.out_dir, name), 'rb') as handle:
            return handle.read()


class ExportingTheMSA(TempDirCase):

    def test_every_genome_is_exported_whatever_its_quality_or_residues(self):
        # genome 3 fails every CheckM threshold gtdb had; genome 4 has 1 residue of 7
        self.assertEqual(self.export(), 4)
        self.assertEqual(sorted(os.listdir(self.out_dir)),
                         sorted(['gtdb' + s for s in (M.MSA_SUFFIX, M.MULTI_HITS_SUFFIX, M.ARB_METADATA_SUFFIX,
                                                     M.ARB_FILTER_SUFFIX, M.TSV_METADATA_SUFFIX,
                                                     M.MARKERS_INFO_SUFFIX)] + [M.PARTITIONS_NAME]))

    def test_the_msa_is_every_column_ordered_by_accession_and_named_by_canonical_accession(self):
        self.export()
        # GenBank before RefSeq, as export_metadata orders them
        self.assertEqual(self.read(M.MSA_SUFFIX), '>G000000002\nACDF---\n>G000000003\nKLMNPQR\n'
                                                  '>G000000004\nA------\n>G000000001\nACDEFGH\n')
        self.assertEqual(self.read(M.MULTI_HITS_SUFFIX).splitlines(), [
            'Genome_ID\tPFAM_PF00001.1\tTIGR_TIGR00001', 'G000000002\tSingle\tMultiple',
            'G000000003\tSingle\tMulti_Identical', 'G000000004\tSingle\tMissing', 'G000000001\tSingle\tSingle'])

    def test_the_partitions_name_the_columns_of_each_marker_in_the_msa(self):
        self.export(prefix='bac', individual=True)
        with open(os.path.join(self.out_dir, M.PARTITIONS_NAME)) as handle:
            self.assertEqual(handle.read(), 'AA, PFAM_PF00001.1 = 1-4\nAA, TIGR_TIGR00001 = 5-7\n')
        # the columns named are those of each marker's own alignment
        msa = self.read(M.MSA_SUFFIX, 'bac').splitlines()[1]
        self.assertEqual((msa[0:4], msa[4:7]), (self.read('_PF00001.1.faa.gz', 'bac').splitlines()[1],
                                                self.read('_TIGR00001.faa.gz', 'bac').splitlines()[1]))

    def test_the_alignments_and_metadata_are_gzipped_the_same_bytes_each_time(self):
        self.export(prefix='gtdb_r237_ar53', individual=True)
        gzipped = sorted(n for n in os.listdir(self.out_dir) if n.endswith('.gz'))
        self.assertEqual(gzipped, ['gtdb_r237_ar53_PF00001.1.faa.gz', 'gtdb_r237_ar53_TIGR00001.faa.gz',
                                   'gtdb_r237_ar53_arb_metadata.txt.gz', 'gtdb_r237_ar53_concatenated.faa.gz',
                                   'gtdb_r237_ar53_metadata.tsv.gz'])
        self.assertFalse([n for n in os.listdir(self.out_dir) if n.endswith(('.faa', '_arb_metadata.txt',
                                                                              '_arb_metadata.tsv', '.partial'))])
        for name in gzipped:
            self.assertEqual(self.read_bytes(name)[:2], b'\x1f\x8b', name)
        first = {name: self.read_bytes(name) for name in gzipped}
        shutil.rmtree(self.out_dir)
        self.export(prefix='gtdb_r237_ar53', individual=True)
        self.assertEqual({name: self.read_bytes(name) for name in gzipped}, first)

    def test_a_genome_named_twice_is_exported_once(self):
        self.assertEqual(self.export(['GCA_000000002.1', 'GB_GCA_000000002.1']), 1)
        self.assertEqual(self.read(M.MSA_SUFFIX), '>G000000002\nACDF---\n')

    def test_the_arb_metadata_and_its_tsv_hold_a_genome_each_with_the_counts_of_the_msa(self):
        self.export()
        arb = self.read(M.ARB_METADATA_SUFFIX)
        self.assertEqual(arb.count('BEGIN\n'), 4)
        self.assertEqual(arb.count('\naligned_seq=\n'), 4)
        self.assertIn('organism_name=Organism 1\n', arb)
        self.assertIn('gtdb_cluster_ncbi_species=s__G sp1: 1; s__G sp4: 1\n', arb)
        rows = [line.split('\t') for line in self.read(M.TSV_METADATA_SUFFIX).splitlines()]
        self.assertEqual(len(rows), 5)
        self.assertEqual({len(row) for row in rows}, {len(rows[0])})
        cells = dict(zip(rows[0], rows[1]))
        self.assertEqual((cells['db_name'], cells['msa_aa_count'], cells['msa_length'], cells['msa_gene_count'],
                          cells['multiple_homologs'], cells['checkm_completeness']),
                         ('G000000002', '4', '7', '1', '1', '95.000'))
        self.assertIn('MATCH\t"msa_aa_count\\=*"', self.read(M.ARB_FILTER_SUFFIX))

    def test_markers_info_and_individual_alignments(self):
        self.export(prefix='bac', individual=True)
        info = self.read(M.MARKERS_INFO_SUFFIX, 'bac').splitlines()
        self.assertEqual(info[1:], ['PFAM_PF00001.1\tmarker11\tMarker 11\t4\t100.00\t100.00',
                                    'TIGR_TIGR00001\tmarker12\tMarker 12\t3\t25.00\t75.00'])
        self.assertEqual(self.read('_TIGR00001.faa.gz', 'bac'),
                         '>G000000002\n---\n>G000000003\nPQR\n>G000000004\n---\n>G000000001\nFGH\n')
        self.assertEqual(self.read('_PF00001.1.faa.gz', 'bac'),
                         '>G000000002\nACDF\n>G000000003\nKLMN\n>G000000004\nA---\n>G000000001\nACDE\n')

    def test_reading_on_several_connections_writes_the_same_files(self):
        self.export(individual=True)
        single = {name: self.read('', name) for name in os.listdir(self.out_dir) if name != M.PARTITIONS_NAME}
        shutil.rmtree(self.out_dir)
        with mock.patch.object(M, 'ALIGNED_CHUNK', 1):
            self.export(individual=True, cpus=3)
        self.assertEqual(sorted(self.manager.chunks_read), [[1], [2], [3], [4]])
        self.assertEqual({name: self.read('', name) for name in single}, single)

    def test_the_transaction_is_rolled_back_whether_or_not_the_export_succeeds(self):
        self.export()
        self.manager.temp_con.rollback.assert_called_once_with()
        with self.assertRaises(M.MSAExportError):
            self.export(marker_set_ids=(77,))
        self.manager.temp_con.rollback.assert_called_once_with()


class RefusingAnExport(TempDirCase):

    def assert_nothing_written(self):
        self.assertFalse(os.path.exists(self.out_dir))

    def test_a_genome_without_its_alignments_is_refused_and_listed(self):
        del self.db.aligned[2][12]
        del self.db.aligned[3]
        with self.assertRaises(M.MSAExportError) as raised:
            self.export()
        self.assertIn('align_marker_genes', str(raised.exception))
        self.assert_nothing_written()
        with open(os.path.join(self.log_dir, M.UNALIGNED_NAME)) as handle:
            self.assertEqual(handle.read(), 'genome_id\tmissing_markers\tmarkers\n'
                                            'GB_GCA_000000002.1\t1\tTIGR_TIGR00001\n'
                                            'GB_GCA_000000003.1\t2\tPFAM_PF00001.1,TIGR_TIGR00001\n')

    def test_unaligned_genomes_found_on_several_connections_are_listed_in_order(self):
        del self.db.aligned[3]
        del self.db.aligned[2]
        with mock.patch.object(M, 'ALIGNED_CHUNK', 1), self.assertRaises(M.MSAExportError):
            self.export(cpus=4)
        with open(os.path.join(self.log_dir, M.UNALIGNED_NAME)) as handle:
            self.assertEqual([line.split('\t')[0] for line in handle][1:],
                             ['GB_GCA_000000002.1', 'GB_GCA_000000003.1'])
        self.assert_nothing_written()

    def test_a_genome_the_database_does_not_hold_is_refused_and_listed(self):
        with self.assertRaises(M.MSAExportError):
            self.export(self.GENOMES + ['GB_GCA_000000099.1'])
        self.assert_nothing_written()
        with open(os.path.join(self.log_dir, M.UNKNOWN_GENOMES_NAME)) as handle:
            self.assertEqual(handle.read(), 'genome_id\nGB_GCA_000000099.1\n')

    def test_an_unknown_marker_set_or_no_genomes_are_refused(self):
        with self.assertRaises(M.MSAExportError):
            self.export(marker_set_ids=(1, 77))
        with self.assertRaises(M.MSAExportError):
            self.export([])
        with self.assertRaises(M.MSAExportError):
            self.export(prefix='a/b')
        self.assert_nothing_written()


class TheReaders(unittest.TestCase):
    """Every connection reading aligned_markers reads the main transaction's snapshot."""

    def test_each_reader_joins_the_exported_snapshot_read_only_before_anything_else(self):
        main_cursor = mock.Mock()
        main_cursor.fetchone.return_value = ('00000003-0000001B-1',)
        connections = []

        def connection(database):
            con = mock.Mock()
            con.cursor.return_value = mock.Mock()
            connections.append(con)
            return con

        with mock.patch.object(M, 'GenomeDatabaseConnectionFTPUpdate', side_effect=connection):
            manager = M.MSAExportManager({'service': 'x'}, cpus=3)
            manager.temp_cur = main_cursor
            manager.open_readers()
        main_cursor.execute.assert_called_once_with(M.EXPORT_SNAPSHOT)
        readers = connections[1:]
        self.assertEqual(len(readers), 3)
        for con in readers:
            con.conn.set_session.assert_called_once_with(isolation_level=M.ISOLATION_LEVEL_REPEATABLE_READ,
                                                         readonly=True)
            statements = con.cursor.return_value.execute.call_args_list
            self.assertEqual(statements[0], mock.call(M.SET_SNAPSHOT, ('00000003-0000001B-1',)))
        manager.close_readers()
        for con in readers:
            con.rollback.assert_called_once_with()
            con.ClosePostgresConnection.assert_called_once_with()
        self.assertEqual(manager.readers, [])


class TheCommandLine(unittest.TestCase):

    def parse(self, *args):
        return main_module.get_main_parser().parse_args(
            ['export_msa', '--db_service', 'x', '--genome_batchfile', 'g.tsv', '-o', 'out'] + list(args))

    def test_marker_sets_are_taken_comma_separated_as_gtdb_took_them_or_space_separated(self):
        options = self.parse('--marker_set_ids', '1,19', '2')
        self.assertEqual(M.flatten(options.marker_set_ids), [1, 19, 2])
        self.assertEqual((options.prefix, options.individual, options.cpus), ('gtdb', False, 1))

    def test_gtdbs_filtering_and_trimming_options_are_not_taken(self):
        for option in ('--quality_threshold', '--min_perc_aa', '--cols_per_gene', '--rnd_seed',
                       '--excluded_genome_ids', '--guaranteed_batchfile', '--no_trim'):
            with self.subTest(option=option), mock.patch('sys.stderr'):
                with self.assertRaises(SystemExit):
                    self.parse('--marker_set_ids', '1', option, '1')

    def test_list_marker_sets_needs_only_a_database(self):
        options = main_module.get_main_parser().parse_args(['list_marker_sets', '--db_service', 'x'])
        self.assertEqual(options.subparser_name, 'list_marker_sets')


if __name__ == '__main__':
    unittest.main()
