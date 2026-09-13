#!/usr/bin/env python3
"""Offline unit tests for the NCBI metadata download -- no network, no NCBI.

Run with the interpreter that has tqdm:

    /opt/centos7/sw/miniconda3/envs/python3/bin/python -m unittest -v tests.test_ncbi_metadata_sync

The downloads are served by a local HTTP server on the loopback interface, so what
is exercised is the real urllib path: a real socket, a real Content-Length, real
blocks. What must not break silently is the naming of the four assembly summary
files, which NCBI publishes under one name and select_genomes reads the database
from, and the refusal to accept a taxonomy dump whose MD5 does not match.
"""

import hashlib
import http.server
import io
import os
import shutil
import tarfile
import tempfile
import threading
import unittest

from gtdb_migration_tk import ncbi_ftp_manager as F
from gtdb_migration_tk import ncbi_ftp_manager_tools as T
from gtdb_migration_tk import ncbi_tax_manager
from gtdb_migration_tk import ncbi_utils as U


# Two lineages rooted at "cellular organisms", which is where the parser stops
# walking, shaped as NCBI writes nodes.dmp: fields separated by "\t|\t".
NODES = [('131567', '1', 'no rank'),
         # Bacteria
         ('2', '131567', 'superkingdom'),
         ('1224', '2', 'phylum'),
         ('1236', '1224', 'class'),
         ('91347', '1236', 'order'),
         ('543', '91347', 'family'),
         ('561', '543', 'genus'),
         ('562', '561', 'species'),
         # Archaea
         ('2157', '131567', 'superkingdom'),
         ('28890', '2157', 'phylum'),
         ('183963', '28890', 'class'),
         ('2235', '183963', 'order'),
         ('2236', '2235', 'family'),
         ('2239', '2236', 'genus'),
         ('2242', '2239', 'species')]

NAMES = [('131567', 'cellular organisms'),
         ('2', 'Bacteria'), ('1224', 'Pseudomonadota'),
         ('1236', 'Gammaproteobacteria'), ('91347', 'Enterobacterales'),
         ('543', 'Enterobacteriaceae'), ('561', 'Escherichia'),
         ('562', 'Escherichia coli'),
         ('2157', 'Archaea'), ('28890', 'Euryarchaeota'),
         ('183963', 'Halobacteria'), ('2235', 'Halobacteriales'),
         ('2236', 'Halobacteriaceae'), ('2239', 'Halobacterium'),
         ('2242', 'Halobacterium salinarum')]


def dmp_files():
    """nodes.dmp and names.dmp holding one complete bacterial lineage."""

    nodes = ''.join('{}\t|\t{}\t|\t{}\t|\t\t|\t0\t|\t\t|\t11\t|\n'.format(*row)
                    for row in NODES)
    names = ''.join('{}\t|\t{}\t|\t\t|\tscientific name\t|\n'.format(*row)
                    for row in NAMES)

    return {'nodes.dmp': nodes.encode(), 'names.dmp': names.encode()}


def make_taxdump(members=None):
    """A gzipped tarball shaped like NCBI's taxdump.tar.gz, held in memory."""

    members = dmp_files() if members is None else members
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode='w:gz') as archive:
        for name, payload in members.items():
            info = tarfile.TarInfo(name)
            info.size = len(payload)
            archive.addfile(info, io.BytesIO(payload))

    return buffer.getvalue()


class ServedFiles(http.server.BaseHTTPRequestHandler):
    """Serve a dict of path -> bytes, and nothing else."""

    routes = {}

    def do_GET(self):
        body = self.routes.get(self.path)
        if body is None:
            self.send_error(404)
            return
        self.send_response(200)
        self.send_header('Content-Length', str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def log_message(self, *args):
        pass                                     # keep the test output readable


class HttpCase(unittest.TestCase):
    """A temporary directory and a local HTTP server serving `routes`."""

    routes = {}

    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='ncbi_metadata_sync_test.')
        handler = type('Handler', (ServedFiles,), {'routes': dict(self.routes)})
        self.server = http.server.HTTPServer(('127.0.0.1', 0), handler)
        self.url = 'http://127.0.0.1:{}'.format(self.server.server_port)
        # the default 0.5 s poll interval is paid back on every shutdown(), which
        # on a class of this size is most of the runtime
        self.thread = threading.Thread(target=self.server.serve_forever,
                                       kwargs={'poll_interval': 0.01}, daemon=True)
        self.thread.start()

    def tearDown(self):
        self.server.shutdown()
        self.server.server_close()
        self.thread.join(timeout=5)
        shutil.rmtree(self.dir, ignore_errors=True)

    def path(self, *parts):
        return os.path.join(self.dir, *parts)


# ------------------------------------------------------------------- download naming

class AssemblySummaryNamingTests(unittest.TestCase):
    """NCBI publishes all four of these under one name; GTDB cannot."""

    def setUp(self):
        self.downloads = T.assembly_summary_downloads()

    def test_all_four_databases_and_domains_are_downloaded(self):
        self.assertEqual(len(self.downloads), 4)

    def test_every_file_is_saved_under_a_distinct_name(self):
        # NCBI calls each of them assembly_summary.txt, so saving them into one
        # directory under their own names would leave a single file
        names = [name for *_, name in self.downloads]
        self.assertEqual(len(set(names)), 4)
        self.assertEqual(sorted(names),
                         ['assembly_summary_archaea_genbank.txt.gz',
                          'assembly_summary_archaea_refseq.txt.gz',
                          'assembly_summary_bacteria_genbank.txt.gz',
                          'assembly_summary_bacteria_refseq.txt.gz'])

    def test_each_name_carries_the_suffix_select_genomes_reads(self):
        # select_genomes decides a file's database from this suffix, so the two
        # must agree or the release is selected from the wrong genomes
        for database, _, _, name in self.downloads:
            self.assertTrue(name.endswith('_{}.txt.gz'.format(database)), name)

    def test_the_name_matches_the_url_it_is_downloaded_from(self):
        for database, domain, url, name in self.downloads:
            self.assertTrue(url.endswith('/assembly_summary.txt'), url)
            self.assertIn('/genomes/{}/{}/'.format(database, domain), url)
            self.assertEqual(name,
                             'assembly_summary_{}_{}.txt.gz'.format(domain, database))

    def test_fungi_are_not_downloaded(self):
        # fungal genomes are a separate procedure
        self.assertFalse([n for *_, n in self.downloads if 'fungi' in n])


# ------------------------------------------------------------------------- downloading

class DownloadFileTests(HttpCase):
    routes = {'/small.txt': b'hello ncbi\n',
              '/big.bin': bytes(range(256)) * 8192}

    def test_file_is_written_with_the_bytes_served(self):
        out = self.path('small.txt')
        written = T.download_file(self.url + '/small.txt', out, quiet=True)
        self.assertEqual(written, 11)
        with open(out, 'rb') as handle:
            self.assertEqual(handle.read(), b'hello ncbi\n')

    def test_a_file_larger_than_one_block_is_written_whole(self):
        # the GenBank bacteria summary is over a gigabyte, so the loop matters
        out = self.path('big.bin')
        written = T.download_file(self.url + '/big.bin', out, quiet=True)
        self.assertEqual(written, len(self.routes['/big.bin']))
        self.assertEqual(os.path.getsize(out), written)

    def test_a_failed_download_leaves_no_file_behind(self):
        # a truncated file that every later step reads as complete is the failure
        # this guards against
        out = self.path('missing.txt')
        with self.assertRaises(Exception):
            T.download_file(self.url + '/absent.txt', out, quiet=True)
        self.assertFalse(os.path.exists(out))
        self.assertFalse(os.path.exists(out + '.partial'))

    def test_no_partial_file_survives_a_success(self):
        out = self.path('small.txt')
        T.download_file(self.url + '/small.txt', out, quiet=True)
        self.assertEqual(os.listdir(self.dir), ['small.txt'])


class ExtractTarballTests(HttpCase):
    def test_members_are_extracted_into_the_directory(self):
        tarball = self.path('taxdump.tar.gz')
        with open(tarball, 'wb') as handle:
            handle.write(make_taxdump())
        T.extract_tarball(tarball, self.path('taxdump_20240914'))
        self.assertEqual(sorted(os.listdir(self.path('taxdump_20240914'))),
                         ['names.dmp', 'nodes.dmp'])

    def test_the_output_directory_is_created(self):
        tarball = self.path('taxdump.tar.gz')
        with open(tarball, 'wb') as handle:
            handle.write(make_taxdump())
        T.extract_tarball(tarball, self.path('does', 'not', 'exist'))
        self.assertTrue(os.path.exists(self.path('does', 'not', 'exist', 'names.dmp')))


# --------------------------------------------------------------------- invalid taxids

class InvalidTaxidReportTests(unittest.TestCase):
    """NCBI builds the summaries and the taxonomy dump independently, so a newly
    registered assembly names a taxid the dump does not carry. On a full release
    that is thousands of genomes, and a line each buries every other warning."""

    def setUp(self):
        self.warnings = []
        self.parser = ncbi_tax_manager.TaxonomyNCBI()
        self.parser.logger = type('Log', (), {
            'warning': lambda _self, message: self.warnings.append(message)})()

    def test_nothing_is_reported_when_every_taxid_resolves(self):
        self.parser._report_invalid_taxids([])
        self.assertEqual(self.warnings, [])

    def test_one_warning_is_written_however_many_are_invalid(self):
        invalid = [('GCA_{:09d}.1'.format(n), str(n)) for n in range(5000)]
        self.parser._report_invalid_taxids(invalid)
        self.assertEqual(len(self.warnings), 1)

    def test_the_warning_carries_the_count_and_three_examples(self):
        invalid = [('GCA_{:09d}.1'.format(n), str(n)) for n in range(5000)]
        self.parser._report_invalid_taxids(invalid)
        warning = self.warnings[0]
        self.assertIn('5,000 assemblies', warning)
        self.assertIn('GCA_000000000.1 (taxid 0)', warning)
        self.assertIn('GCA_000000002.1 (taxid 2)', warning)
        self.assertNotIn('GCA_000000003.1', warning)

    def test_fewer_than_three_are_all_reported(self):
        self.parser._report_invalid_taxids([('GCA_000000001.1', '987210895')])
        self.assertIn('1 assemblies', self.warnings[0])
        self.assertIn('GCA_000000001.1 (taxid 987210895)', self.warnings[0])


# ----------------------------------------------------------------- the command itself

TAXDUMP = make_taxdump()
TAXDUMP_MD5 = hashlib.md5(TAXDUMP).hexdigest()

# enough of an assembly summary for the taxonomy parser, which reads taxid and
# organism_name by column name and walks the lineage of each assembly
# One genome per file. They must differ: the taxonomy parser refuses an accession
# it has already seen, so four copies of one row would abort the run.
GENOMES = {('refseq', 'archaea'): ('GCF_000000001.1', '2242', 'Halobacterium salinarum'),
           ('refseq', 'bacteria'): ('GCF_000000002.1', '562', 'Escherichia coli'),
           ('genbank', 'archaea'): ('GCA_000000003.1', '2242', 'Halobacterium salinarum'),
           ('genbank', 'bacteria'): ('GCA_000000004.1', '562', 'Escherichia coli')}


def summary_file(database, domain):
    """An assembly summary holding one genome.

    The comment line matters: the taxonomy parser skips the first line of an
    assembly summary and reads the second as the header, as NCBI writes them.
    """

    accession, taxid, organism = GENOMES[(database, domain)]

    return ('#   See ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt\n'
            '#assembly_accession\ttaxid\torganism_name\tinfraspecific_name\tftp_path\n'
            '{}\t{}\t{}\tstrain=x\tftp://a\n'.format(accession, taxid, organism)).encode()


class MetadataSyncTests(HttpCase):
    routes = dict(
        {'/pub/taxonomy/taxdump.tar.gz': TAXDUMP,
         '/pub/taxonomy/taxdump.tar.gz.md5':
             '{}  taxdump.tar.gz\n'.format(TAXDUMP_MD5).encode()},
        **{'/genomes/{}/{}/assembly_summary.txt'.format(database, domain):
           summary_file(database, domain)
           for database in T.NCBI_DATABASES for domain in T.NCBI_DOMAINS})

    def setUp(self):
        super().setUp()
        # point the download at the local server instead of NCBI
        self.ncbi_ftp, T.NCBI_FTP = T.NCBI_FTP, self.url
        self.taxdump_url, T.TAXDUMP_URL = T.TAXDUMP_URL, self.url + '/pub/taxonomy/taxdump.tar.gz'
        F.TAXDUMP_URL = T.TAXDUMP_URL

    def tearDown(self):
        T.NCBI_FTP, T.TAXDUMP_URL = self.ncbi_ftp, self.taxdump_url
        F.TAXDUMP_URL = self.taxdump_url
        super().tearDown()

    RELEASE = 237

    def run_sync(self):
        F.MetadataSyncManager(self.dir).run(self.RELEASE)
        return sorted(os.listdir(self.dir))

    def stamp(self):
        import datetime
        return datetime.date.today().strftime('%Y%m%d')

    def test_assembly_summaries_land_in_the_root_ncbi_directory(self):
        # this is the directory select_genomes is pointed at
        written = self.run_sync()
        self.assertEqual(
            sorted(n for n in written if n.startswith('assembly_summary_')),
            ['assembly_summary_archaea_genbank.txt.gz',
             'assembly_summary_archaea_refseq.txt.gz',
             'assembly_summary_bacteria_genbank.txt.gz',
             'assembly_summary_bacteria_refseq.txt.gz'])

    def test_the_root_holds_only_the_summaries_the_log_and_the_taxonomy(self):
        written = self.run_sync()
        self.assertEqual([n for n in written if not n.startswith('assembly_summary_')],
                         ['taxonomy'])

    def test_everything_derived_from_the_taxonomy_is_under_one_directory(self):
        # the taxonomy of a release is then one directory to copy or archive
        self.run_sync()
        self.assertEqual(sorted(os.listdir(self.path('taxonomy'))),
                         ['standardised_taxonomy', 'taxdump_{0}'.format(self.stamp())])

    def test_the_spent_archive_is_removed_once_it_is_unpacked(self):
        # it has been checked against NCBI and its contents are on disk, so
        # keeping it holds 80 MB for a check that can no longer fail
        self.run_sync()
        extracted = self.path('taxonomy', 'taxdump_' + self.stamp())
        self.assertTrue(os.path.exists(os.path.join(extracted, 'names.dmp')))
        self.assertFalse(os.path.exists(extracted + '.tar.gz'))
        self.assertFalse(os.path.exists(extracted + '.tar.gz.md5'))

    def test_the_taxonomy_directory_is_named_for_the_day_it_was_downloaded(self):
        self.run_sync()
        self.assertTrue(os.path.isdir(self.path('taxonomy', 'taxdump_' + self.stamp())))

    def test_the_standardised_taxonomy_is_named_for_the_release(self):
        self.run_sync()
        written = os.listdir(self.path('taxonomy', 'standardised_taxonomy'))
        taxonomy = [n for n in written if n != 'failed_filters.tsv']
        self.assertTrue(taxonomy)
        for name in taxonomy:
            self.assertTrue(name.startswith('ncbi_r237_'), name)

    def test_the_failed_filter_report_is_written_beside_the_taxonomy(self):
        # it names a hardcoded file, so without a working directory to anchor it
        # the report lands wherever the operator happened to be standing
        self.run_sync()
        self.assertIn('failed_filters.tsv', os.listdir(self.path('taxonomy', 'standardised_taxonomy')))

    def test_the_standardised_taxonomy_resolves_the_lineage_of_each_assembly(self):
        # the point of the step: a taxid becomes a lineage
        self.run_sync()
        taxonomy = self.path('taxonomy', 'standardised_taxonomy', 'ncbi_r237_unfiltered_taxonomy.tsv')
        with open(taxonomy) as handle:
            lineages = dict(line.rstrip('\n').split('\t', 1) for line in handle)

        self.assertEqual(sorted(lineages), ['GCA_000000003.1', 'GCA_000000004.1',
                                            'GCF_000000001.1', 'GCF_000000002.1'])
        self.assertIn('Escherichia coli', lineages['GCF_000000002.1'])
        self.assertIn('Bacteria', lineages['GCF_000000002.1'])
        self.assertIn('Halobacterium salinarum', lineages['GCF_000000001.1'])
        self.assertIn('Archaea', lineages['GCF_000000001.1'])

    def test_a_corrupt_taxonomy_download_stops_the_run(self):
        # a truncated taxdump reads as a valid, smaller taxonomy, and only shows
        # up as genomes that lost their NCBI lineage several commands later
        self.server.RequestHandlerClass.routes['/pub/taxonomy/taxdump.tar.gz.md5'] = \
            b'0' * 32 + b'  taxdump.tar.gz\n'
        with self.assertRaises(SystemExit):
            F.MetadataSyncManager(self.dir).run(self.RELEASE)

    def test_a_missing_file_at_ncbi_stops_the_run(self):
        del self.server.RequestHandlerClass.routes['/genomes/genbank/bacteria/assembly_summary.txt']
        with self.assertRaises(SystemExit):
            F.MetadataSyncManager(self.dir).run(self.RELEASE)

    def test_the_summaries_are_stored_gzipped(self):
        self.run_sync()
        for name in os.listdir(self.dir):
            if not name.startswith('assembly_summary_'):
                continue
            with open(self.path(name), 'rb') as raw:
                self.assertEqual(raw.read(2), b'\x1f\x8b', name)      # gzip magic number

    def test_a_stored_summary_still_reads_as_an_assembly_summary(self):
        # compressing it must not put it beyond the reader every command uses
        self.run_sync()
        rows = list(U.read_assembly_summary(
            self.path('assembly_summary_bacteria_refseq.txt.gz'),
            'assembly_accession', 'taxid'))
        self.assertEqual(rows, [('GCF_000000002.1', '562')])

    def test_no_uncompressed_summary_is_left_behind(self):
        # the uncompressed GenBank bacteria summary is 1.5 GB, so it is never
        # written in the first place rather than written and removed
        self.run_sync()
        self.assertFalse([n for n in os.listdir(self.dir)
                          if n.startswith('assembly_summary_') and not n.endswith('.gz')])

    def test_the_downloaded_summaries_are_what_select_genomes_accepts(self):
        # the two commands meet here: whatever this writes, select_genomes must
        # be able to sort into RefSeq files and GenBank files
        self.run_sync()
        summaries = [self.path(n) for n in os.listdir(self.dir)
                     if n.startswith('assembly_summary_')]
        refseq, genbank = F.SelectedGenomesManager(self.dir).group_by_database(summaries)
        self.assertEqual(len(refseq), 2)
        self.assertEqual(len(genbank), 2)
