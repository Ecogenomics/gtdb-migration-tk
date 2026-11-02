#!/usr/bin/env python3
"""Offline unit tests for lpsn.py -- no LPSN website, no HTML.

Run with the interpreter that has pandas, bs4, requests and sqlalchemy:

    /opt/centos7/sw/miniconda3/envs/python3/bin/python -m unittest -v tests.test_lpsn

What is tested here is the bookkeeping of pull_html: a page whose name
duplicates a validly published one is skipped, recorded in <rank>_skipped.lst
rather than the failed list, and that file exists only when there is something
in it. Nothing is downloaded: the skip is decided before any request is made.

Also tested is which of a name's several LPSN pages parse_all_ranks_tsv keeps,
since that page's authority becomes the name's year of priority.
"""

import csv
import io
import os
import shutil
import tempfile
import unittest

from gtdb_migration_tk.lpsn import LPSN


class TempDirCase(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='lpsn_test.')
        self.lpsn = LPSN(skip_taxa_per_letter_dl=True, lpsn_output_dir=self.dir)

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def skipped_file(self, rank):
        return os.path.join(self.dir, '{}_skipped.lst'.format(rank))


class SkippedFile(TempDirCase):
    def test_written_only_when_something_was_skipped(self):
        self.assertIsNone(self.lpsn.write_skipped('genus', []))
        self.assertFalse(os.path.exists(self.skipped_file('genus')))

    def test_lists_each_skipped_page_with_its_name_and_reason(self):
        skipped = [('"Alcaligenes faecalis subsp. faecalis"',
                    'https://lpsn.dsmz.de/subspecies/alcaligenes-faecalis-faecalis-1')]

        path = self.lpsn.write_skipped('subspecies', skipped)

        self.assertEqual(path, self.skipped_file('subspecies'))
        with open(path) as handle:
            self.assertEqual(handle.read(),
                             '"Alcaligenes faecalis subsp. faecalis"\t'
                             'https://lpsn.dsmz.de/subspecies/alcaligenes-faecalis-faecalis-1\t'
                             'duplicate_name\n')

    def test_a_file_from_an_earlier_run_is_removed_when_nothing_is_skipped_now(self):
        # the file's presence is the signal, so a stale one would misreport this run
        with open(self.skipped_file('genus'), 'w') as handle:
            handle.write('old\n')

        self.lpsn.write_skipped('genus', [])

        self.assertFalse(os.path.exists(self.skipped_file('genus')))


class DuplicateNames(TempDirCase):
    def test_a_duplicate_name_is_skipped_not_failed(self):
        failed = io.StringIO()
        skipped = []
        valid_names = ['Alcaligenes faecalis subsp. faecalis']   # the unquoted twin

        n = self.lpsn.download_rank_name(
            '"Alcaligenes faecalis subsp. faecalis"',
            'https://lpsn.dsmz.de/subspecies/alcaligenes-faecalis-faecalis-1',
            os.path.join(self.dir, 'alcaligenes-faecalis-faecalis-1'),
            valid_names, failed, skipped, num_already_dl=0)

        self.assertEqual(skipped, [('"Alcaligenes faecalis subsp. faecalis"',
                                    'https://lpsn.dsmz.de/subspecies/alcaligenes-faecalis-faecalis-1')])
        self.assertEqual(failed.getvalue(), '')
        self.assertEqual(n, 0)
        # and nothing was fetched or written for it
        self.assertEqual(os.listdir(self.dir), [])

    def test_an_already_downloaded_valid_page_is_counted_not_refetched(self):
        out_file = os.path.join(self.dir, 'alcaligenes-faecalis-faecalis')
        open(out_file, 'w').close()
        failed = io.StringIO()
        skipped = []

        n = self.lpsn.download_rank_name(
            'Alcaligenes faecalis subsp. faecalis',
            'https://lpsn.dsmz.de/subspecies/alcaligenes-faecalis-faecalis',
            out_file, ['Alcaligenes faecalis subsp. faecalis'], failed, skipped,
            num_already_dl=0)

        self.assertEqual(n, 1)
        self.assertEqual(skipped, [])
        self.assertEqual(failed.getvalue(), '')


RAW_HEADER = ['Rank', 'Name', 'Category', 'Type genus', 'Nomenclatural status',
              'Taxonomic status', 'Parent taxon', 'Synonyms', 'Child taxa',
              'Nomenclatural pro-status', 'Original publication', 'Type class',
              'Type order', 'Type species', 'Correct name']


class OnePagePerName(TempDirCase):
    """parse_all_ranks_tsv keeps one page of a name LPSN gives several, and its
    authority is the name's year of priority: the page kept must be the one whose
    name has standing, whichever order the pages are read in."""

    def kept_authority(self, pages):
        """The 'Priority' parse_all_ranks_tsv keeps for one species name given its pages.

        Each page is (name with authority, nomenclatural status, pro-status, taxonomic status).
        """
        raw = os.path.join(self.dir, 'full_parsing_raw.tsv')
        with open(raw, 'w') as handle:
            handle.write('\t'.join(RAW_HEADER) + '\n')
            for name, status, pro_status, taxonomic_status in pages:
                row = dict.fromkeys(RAW_HEADER, 'n/a')
                row.update({'Rank': 'species', 'Name': name, 'Nomenclatural status': status,
                            'Nomenclatural pro-status': pro_status,
                            'Taxonomic status': taxonomic_status})
                handle.write('\t'.join(row[h] for h in RAW_HEADER) + '\n')

        self.lpsn.parse_all_ranks_tsv(raw)

        with open(os.path.join(self.dir, 'full_parsing_parsed.tsv')) as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual(len(rows), 1)
        return rows[0]['Priority']

    def assertKeptEitherWay(self, pages, authority):
        self.assertEqual(self.kept_authority(pages), authority)
        self.assertEqual(self.kept_authority(pages[::-1]), authority)

    def test_a_validly_published_page_is_kept_over_a_quoted_one(self):
        self.assertKeptEitherWay(
            [('"Acinetobacter pecorum" Gilroy et al. 2021', 'not validly published', 'n/a',
              'synonym (and no standing)'),
             ('Acinetobacter pecorum Gilroy et al. 2024', 'validly published under the ICNP', 'n/a',
              'correct name')],
            'Gilroy et al. 2024')

    def test_a_validly_published_page_is_kept_over_a_candidatus_one(self):
        # D. audaxviator once LPSN had validated it: the Candidatus page, read first, gave 2008
        self.assertKeptEitherWay(
            [('"Candidatus Desulforudis audaxviator" Chivian et al. 2008', 'not validly published',
              'pro-validly published under the ICNP', 'pro-correct name'),
             ('Desulforudis audaxviator (ex Chivian et al. 2008) Karnachuk et al. 2026',
              'validly published under the ICNP', 'n/a', 'correct name')],
            '(ex Chivian et al. 2008) Karnachuk et al. 2026')

    def test_a_pro_validly_published_candidatus_page_is_kept_over_its_later_synonym(self):
        # the Candidatus names of Chuvochina et al. 2023 are synonyms of earlier pro-valid ones
        self.assertKeptEitherWay(
            [('"Candidatus Binatus soli" Chuvochina et al. 2023', 'not validly published',
              'not pro-validly published', 'synonym (and no pro-standing)'),
             ('"Candidatus Binatus soli" Chuvochina et al. 2019', 'not validly published',
              'pro-validly published under the ICNP', 'pro-correct name')],
            'Chuvochina et al. 2019')

    def test_among_pages_without_standing_lpsns_taxonomic_status_decides(self):
        self.assertKeptEitherWay(
            [('"Candidatus Bartonella massiliensis" Mediannikov et al. 2014', 'not validly published',
              'not pro-validly published', 'synonym (and no pro-standing)'),
             ('"Bartonella massiliensis" Medkour et al. 2019', 'not validly published', 'n/a',
              'preferred name (not correct name)')],
            'Medkour et al. 2019')

    def test_an_illegitimate_page_is_replaced_even_by_one_not_validly_published(self):
        self.assertKeptEitherWay(
            [('Cyanospira capsulata Florenzano et al. 2015',
              'validly published under the ICN (Botanical Code), illegitimate name, later homonym',
              'n/a', 'in need of a replacement'),
             ('"Cyanospira capsulata" Florenzano et al. 1985', 'not validly published', 'n/a',
              'synonym (and no standing)')],
            'Florenzano et al. 1985')

    def test_a_page_with_no_status_is_below_one_with_any(self):
        self.assertKeptEitherWay(
            [('"Xenia rubra" Smith 2001', 'n/a', 'n/a', 'n/a'),
             ('"Xenia rubra" Jones 2005', 'not validly published', 'n/a', 'synonym (and no standing)')],
            'Jones 2005')

    def test_pages_of_equal_standing_keep_the_first_read(self):
        pages = [('"Methylopila turkiensis" Agafonova et al. 2015', 'inaccurate spelling', 'n/a',
                  'misspelling (and no standing)'),
                 ('Methylopila turkiensis Shmareva et al. 2017', 'inaccurate spelling', 'n/a',
                  'misspelling')]

        self.assertEqual(self.kept_authority(pages), 'Agafonova et al. 2015')
        self.assertEqual(self.kept_authority(pages[::-1]), 'Shmareva et al. 2017')


if __name__ == '__main__':
    unittest.main()
