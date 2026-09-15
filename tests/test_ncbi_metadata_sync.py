#!/usr/bin/env python3
"""Offline unit tests for the NCBI metadata download -- no network, no NCBI.

Run with the interpreter that has tqdm:

    /opt/centos7/sw/miniconda3/envs/python3/bin/python -m unittest -v tests.test_ncbi_metadata_sync

The downloads are served by a local HTTP server on the loopback interface, so what
is exercised is the real urllib path: a real socket, a real Content-Length, real
blocks. What must not break silently is the naming of the assembly summary
files, which NCBI publishes under one name and select_genomes reads the database
from, and the refusal to accept a taxonomy dump whose MD5 does not match.
"""

import hashlib
import http.server
import io
import os
import re
import shutil
import tarfile
import tempfile
import threading
import unittest

from gtdb_migration_tk import ncbi_metadata_sync as M
from gtdb_migration_tk import ncbi_tax_manager
from gtdb_migration_tk import ncbi_utils as U
from gtdb_migration_tk import select_genomes as S


# Three lineages rooted at "cellular organisms", which is where the parser stops
# walking, shaped as NCBI writes nodes.dmp: fields separated by "\t|\t". The
# fungal one carries a subkingdom and a subphylum, the ranks the fungal taxonomy
# is generated to keep and the prokaryotic one to drop.
NODES = [('131567', '1', 'no rank'),
         # Bacteria; NCBI ranks the two prokaryotic domains 'domain', and Fungi
         # a kingdom under a superkingdom, which is why one standardises to a
         # d__ and the other reaches it through a rank the parser drops
         ('2', '131567', 'domain'),
         ('1224', '2', 'phylum'),
         ('1236', '1224', 'class'),
         ('91347', '1236', 'order'),
         ('543', '91347', 'family'),
         ('561', '543', 'genus'),
         ('562', '561', 'species'),
         # Archaea
         ('2157', '131567', 'domain'),
         ('28890', '2157', 'phylum'),
         ('183963', '28890', 'class'),
         ('2235', '183963', 'order'),
         ('2236', '2235', 'family'),
         ('2239', '2236', 'genus'),
         ('2242', '2239', 'species'),
         # A lineage NCBI holds under the SeqCode, which it says by appending the
         # code to the name; ranks published under both codes carry it and the
         # rest do not
         ('2802426', '2', 'phylum'),
         ('2802427', '2802426', 'class'),
         ('2802428', '2802427', 'order'),
         ('2802429', '2802428', 'family'),
         ('2802430', '2802429', 'genus'),
         ('2802431', '2802430', 'species'),
         # Fungi
         ('2759', '131567', 'superkingdom'),
         ('4751', '2759', 'kingdom'),
         ('451864', '4751', 'subkingdom'),
         ('4890', '451864', 'phylum'),
         ('147537', '4890', 'subphylum'),
         ('4891', '147537', 'class'),
         ('4892', '4891', 'order'),
         ('766764', '4892', 'family'),
         ('1535326', '766764', 'genus'),
         ('5476', '1535326', 'species')]

NAMES = [('131567', 'cellular organisms'),
         ('2', 'Bacteria'), ('1224', 'Pseudomonadota'),
         ('1236', 'Gammaproteobacteria'), ('91347', 'Enterobacterales'),
         ('543', 'Enterobacteriaceae'), ('561', 'Escherichia'),
         ('562', 'Escherichia coli'),
         ('2802426', 'Patescibacteria (SeqCode)'),
         ('2802427', 'Patescibacteriia'),
         ('2802428', 'Patescibacteriales'),
         ('2802429', 'Patescibacteriaceae (SeqCode)'),
         ('2802430', 'Patescibacter'),
         ('2802431', 'Patescibacter aquaticus (SeqCode)'),
         ('2157', 'Archaea'), ('28890', 'Euryarchaeota'),
         ('183963', 'Halobacteria'), ('2235', 'Halobacteriales'),
         ('2236', 'Halobacteriaceae'), ('2239', 'Halobacterium'),
         ('2242', 'Halobacterium salinarum'),
         ('2759', 'Eukaryota'), ('4751', 'Fungi'), ('451864', 'Dikarya'),
         ('4890', 'Ascomycota'), ('147537', 'Saccharomycotina'),
         ('4891', 'Saccharomycetes'), ('4892', 'Saccharomycetales'),
         ('766764', 'Debaryomycetaceae'), ('1535326', 'Candida'),
         ('5476', 'Candida albicans')]


def dmp_files():
    """nodes.dmp and names.dmp holding one lineage of each group."""

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
    """NCBI publishes every one of these under one name; GTDB cannot."""

    def setUp(self):
        self.downloads = {group: M.assembly_summary_downloads(group)
                          for group in M.NCBI_GROUPS}
        self.every = [row for rows in self.downloads.values() for row in rows]

    def test_each_group_downloads_its_own_domains_from_both_databases(self):
        for group, rows in self.downloads.items():
            self.assertEqual(len(rows),
                             len(M.NCBI_DATABASES) * len(M.NCBI_GROUP_DOMAINS[group]),
                             group)

    def test_the_prokaryotic_group_is_archaea_and_bacteria(self):
        self.assertEqual(sorted(n for *_, n in self.downloads[M.GROUP_PROK]),
                         ['assembly_summary_archaea_genbank.txt.gz',
                          'assembly_summary_archaea_refseq.txt.gz',
                          'assembly_summary_bacteria_genbank.txt.gz',
                          'assembly_summary_bacteria_refseq.txt.gz'])

    def test_the_fungal_group_is_the_fungal_summary_of_each_database(self):
        self.assertEqual(sorted(n for *_, n in self.downloads[M.GROUP_FUNGI]),
                         ['assembly_summary_fungi_genbank.txt.gz',
                          'assembly_summary_fungi_refseq.txt.gz'])

    def test_no_file_is_downloaded_by_both_groups(self):
        # the groups are run separately, and neither may overwrite the other's
        # summary files in a release directory holding both
        names = [name for *_, name in self.every]
        self.assertEqual(len(set(names)), len(names))

    def test_every_file_is_saved_under_a_distinct_name(self):
        # NCBI calls each of them assembly_summary.txt, so saving them into one
        # directory under their own names would leave a single file
        for group, rows in self.downloads.items():
            names = [name for *_, name in rows]
            self.assertEqual(len(set(names)), len(names), group)

    def test_each_name_carries_the_suffix_select_genomes_reads(self):
        # select_genomes decides a file's database from this suffix, so the two
        # must agree or the release is selected from the wrong genomes
        for database, _, _, name in self.every:
            self.assertTrue(name.endswith('_{}.txt.gz'.format(database)), name)

    def test_the_name_matches_the_url_it_is_downloaded_from(self):
        for database, domain, url, name in self.every:
            self.assertTrue(url.endswith('/assembly_summary.txt'), url)
            self.assertIn('/genomes/{}/{}/'.format(database, domain), url)
            self.assertEqual(name,
                             'assembly_summary_{}_{}.txt.gz'.format(domain, database))

    def test_an_unknown_group_downloads_nothing(self):
        # rather than quietly falling back on the prokaryotic domains
        with self.assertRaises(KeyError):
            M.assembly_summary_downloads('EUK')


class GroupTests(unittest.TestCase):
    """The two groups are handled differently after they are downloaded."""

    def test_only_the_fungal_taxonomy_keeps_subranks(self):
        # GTDB curates 7 ranks for prokaryotes; fungal classification leans on
        # NCBI's intermediate ranks
        self.assertFalse(M.GROUP_KEEP_SUBRANKS[M.GROUP_PROK])
        self.assertTrue(M.GROUP_KEEP_SUBRANKS[M.GROUP_FUNGI])

    def test_every_group_says_what_it_downloads_how_it_is_named_and_its_subranks(self):
        # a group added to NCBI_GROUPS and nowhere else fails here rather than
        # part way through a download
        for group in M.NCBI_GROUPS:
            self.assertIn(group, M.NCBI_GROUP_DOMAINS)
            self.assertIn(group, M.GROUP_KEEP_SUBRANKS)
            self.assertIn(group, M.GROUP_FILE_TAG)

    def test_the_groups_name_their_output_files_differently(self):
        # both groups are run into one release directory
        tags = [M.GROUP_FILE_TAG[group] for group in M.NCBI_GROUPS]
        self.assertEqual(len(set(tags)), len(tags))

    def test_no_group_takes_a_domain_of_another(self):
        domains = [domain for group in M.NCBI_GROUPS
                   for domain in M.NCBI_GROUP_DOMAINS[group]]
        self.assertEqual(len(set(domains)), len(domains))


# ------------------------------------------------------------------------- downloading

class DownloadFileTests(HttpCase):
    routes = {'/small.txt': b'hello ncbi\n',
              '/big.bin': bytes(range(256)) * 8192}

    def test_file_is_written_with_the_bytes_served(self):
        out = self.path('small.txt')
        written = M.download_file(self.url + '/small.txt', out, quiet=True)
        self.assertEqual(written, 11)
        with open(out, 'rb') as handle:
            self.assertEqual(handle.read(), b'hello ncbi\n')

    def test_a_file_larger_than_one_block_is_written_whole(self):
        # the GenBank bacteria summary is over a gigabyte, so the loop matters
        out = self.path('big.bin')
        written = M.download_file(self.url + '/big.bin', out, quiet=True)
        self.assertEqual(written, len(self.routes['/big.bin']))
        self.assertEqual(os.path.getsize(out), written)

    def test_a_failed_download_leaves_no_file_behind(self):
        # a truncated file that every later step reads as complete is the failure
        # this guards against
        out = self.path('missing.txt')
        with self.assertRaises(Exception):
            M.download_file(self.url + '/absent.txt', out, quiet=True)
        self.assertFalse(os.path.exists(out))
        self.assertFalse(os.path.exists(out + '.partial'))

    def test_no_partial_file_survives_a_success(self):
        out = self.path('small.txt')
        M.download_file(self.url + '/small.txt', out, quiet=True)
        self.assertEqual(os.listdir(self.dir), ['small.txt'])


class ExtractTarballTests(HttpCase):
    def test_members_are_extracted_into_the_directory(self):
        tarball = self.path('taxdump.tar.gz')
        with open(tarball, 'wb') as handle:
            handle.write(make_taxdump())
        M.extract_tarball(tarball, self.path('taxdump_20240914'))
        self.assertEqual(sorted(os.listdir(self.path('taxdump_20240914'))),
                         ['names.dmp', 'nodes.dmp'])

    def test_the_output_directory_is_created(self):
        tarball = self.path('taxdump.tar.gz')
        with open(tarball, 'wb') as handle:
            handle.write(make_taxdump())
        M.extract_tarball(tarball, self.path('does', 'not', 'exist'))
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


# ------------------------------------------------------------------- names.dmp parsing

class NameParsingTests(unittest.TestCase):
    """NCBI writes the nomenclatural code into the name of a taxon it holds under
    the SeqCode, and GTDB wants the name alone."""

    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='ncbi_names_test.')
        self.parser = ncbi_tax_manager.TaxonomyNCBI()

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def names_dmp(self, *records):
        """A names.dmp holding the given (tax_id, name, name class) records."""

        path = os.path.join(self.dir, 'names.dmp')
        with open(path, 'w') as handle:
            for tax_id, name, name_class in records:
                handle.write('{}\t|\t{}\t|\t\t|\t{}\t|\n'.format(tax_id, name, name_class))

        return path

    def read(self, *records):
        return {tax_id: record.name_txt for tax_id, record
                in self.parser._read_names(self.names_dmp(*records)).items()}

    def test_the_seqcode_annotation_is_removed_from_a_name(self):
        self.assertEqual(
            self.read(('2802429', 'Patescibacteriaceae (SeqCode)', 'scientific name')),
            {'2802429': 'Patescibacteriaceae'})

    def test_a_name_without_an_annotation_is_left_as_it_is(self):
        self.assertEqual(
            self.read(('562', 'Escherichia coli', 'scientific name')),
            {'562': 'Escherichia coli'})

    def test_an_annotation_is_removed_from_a_name_of_any_rank(self):
        # it is appended to species and phyla alike
        self.assertEqual(
            self.read(('2802426', 'Patescibacteria (SeqCode)', 'scientific name'),
                      ('2802431', 'Patescibacter aquaticus (SeqCode)', 'scientific name')),
            {'2802426': 'Patescibacteria', '2802431': 'Patescibacter aquaticus'})

    def test_an_annotation_elsewhere_in_a_name_is_left_alone(self):
        # only a trailing code is an annotation on the name
        self.assertEqual(
            self.read(('1', '(SeqCode) sensu lato', 'scientific name')),
            {'1': '(SeqCode) sensu lato'})

    def test_a_name_that_is_nothing_but_the_annotation_is_kept(self):
        # rather than putting an empty taxon into a lineage
        self.assertEqual(self.read(('1', '(SeqCode)', 'scientific name')),
                         {'1': '(SeqCode)'})

    def test_only_scientific_names_are_read(self):
        # the contract the stripping sits inside
        self.assertEqual(
            self.read(('562', 'Escherichia coli', 'scientific name'),
                      ('562', 'E. coli', 'equivalent name'),
                      ('562', 'Bacillus coli (SeqCode)', 'synonym')),
            {'562': 'Escherichia coli'})


# ------------------------------------------------------------------- log formatting

class LogFormattingTests(unittest.TestCase):
    """A release is counted in millions of records, and "Read 3013402 node
    records" cannot be read at a glance; "Read 3,013,402" can."""

    FILLER = 1234                                # enough records to need a comma

    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix='ncbi_log_test.')
        self.messages = []
        self.parser = ncbi_tax_manager.TaxonomyNCBI()
        self.parser.logger = type('Log', (), {
            'info': lambda _self, message: self.messages.append(message),
            'warning': lambda _self, message: self.messages.append(message),
            'error': lambda _self, message: self.messages.append(message)})()

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def taxonomy_dir(self):
        """A taxdump of the lineages above, padded out to a readable size."""

        # the filler taxa hang off the root and are named by no assembly; they
        # are here to be counted, not walked
        filler = [(str(9000000 + n), '1', 'no rank') for n in range(self.FILLER)]
        nodes = NODES + filler
        names = NAMES + [(tax_id, 'Filler taxon {}'.format(tax_id))
                         for tax_id, _, _ in filler]

        for name, payload in {'nodes.dmp': nodes, 'names.dmp': names}.items():
            with open(os.path.join(self.dir, name), 'w') as handle:
                for record in payload:
                    if len(record) == 3:
                        handle.write('{}\t|\t{}\t|\t{}\t|\t\t|\t0\t|\t\t|\t11\t|\n'.format(*record))
                    else:
                        handle.write('{}\t|\t{}\t|\t\t|\tscientific name\t|\n'.format(*record))

        return self.dir

    def assembly_summary(self):
        """An assembly summary holding FILLER genomes, all of one species."""

        path = os.path.join(self.dir, 'assembly_summary_bacteria_refseq.txt')
        with open(path, 'w') as handle:
            handle.write('#   See ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt\n')
            handle.write('#assembly_accession\ttaxid\torganism_name\tinfraspecific_name\tftp_path\n')
            for n in range(self.FILLER):
                handle.write('GCF_{:09d}.1\t562\tEscherichia coli\tstrain=x\tftp://a\n'.format(n))

        return path

    def parse(self):
        self.parser.parse_ncbi_taxonomy(self.taxonomy_dir(),
                                        [self.assembly_summary()],
                                        False,
                                        os.path.join(self.dir, 'ncbi_r237_prok'),
                                        os.path.join(self.dir, 'prok_failed_filters.tsv'))
        return self.messages

    def test_the_node_record_count_separates_thousands(self):
        self.assertIn('Read {:,} node records.'.format(len(NODES) + self.FILLER),
                      self.parse())

    def test_the_name_record_count_separates_thousands(self):
        self.assertIn('Read {:,} name records.'.format(len(NAMES) + self.FILLER),
                      self.parse())

    def test_the_assembly_count_separates_thousands(self):
        self.assertIn('Number of assemblies: {:,}'.format(self.FILLER), self.parse())

    def test_no_count_is_written_as_a_bare_run_of_digits(self):
        # the whole of what this class is for: a count of four digits or more
        # carries its separators wherever it is logged. Paths are not counts --
        # a release number, a date stamp and a temporary directory are all
        # written as they stand -- so only the words are read
        for message in self.parse():
            for word in message.split():
                if '/' in word:
                    continue
                for number in re.findall(r'(?<![\d,])\d+', word):
                    self.assertLess(len(number), 4, message)


# ----------------------------------------------------------------- the command itself

TAXDUMP = make_taxdump()
TAXDUMP_MD5 = hashlib.md5(TAXDUMP).hexdigest()

# enough of an assembly summary for the taxonomy parser, which reads taxid and
# organism_name by column name and walks the lineage of each assembly
# The accessions must differ: the taxonomy parser refuses one it has already
# seen, so copies of a row in two files would abort the run.
SEQCODE_GENOME = ('GCF_000000007.1', '2802431', 'Patescibacter aquaticus')

GENOMES = {('refseq', 'archaea'): [('GCF_000000001.1', '2242', 'Halobacterium salinarum')],
           ('refseq', 'bacteria'): [('GCF_000000002.1', '562', 'Escherichia coli'),
                                    SEQCODE_GENOME],
           ('refseq', 'fungi'): [('GCF_000000005.1', '5476', 'Candida albicans')],
           ('genbank', 'archaea'): [('GCA_000000003.1', '2242', 'Halobacterium salinarum')],
           ('genbank', 'bacteria'): [('GCA_000000004.1', '562', 'Escherichia coli')],
           ('genbank', 'fungi'): [('GCA_000000006.1', '5476', 'Candida albicans')]}


def summary_file(database, domain):
    """An assembly summary holding the genomes of one database and domain.

    The comment line matters: the taxonomy parser skips the first line of an
    assembly summary and reads the second as the header, as NCBI writes them.
    """

    rows = ''.join('{}\t{}\t{}\tstrain=x\tftp://a\n'.format(*genome)
                   for genome in GENOMES[(database, domain)])

    return ('#   See ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt\n'
            '#assembly_accession\ttaxid\torganism_name\tinfraspecific_name\tftp_path\n'
            + rows).encode()


class MetadataSyncTests(HttpCase):
    routes = dict(
        {'/pub/taxonomy/taxdump.tar.gz': TAXDUMP,
         '/pub/taxonomy/taxdump.tar.gz.md5':
             '{}  taxdump.tar.gz\n'.format(TAXDUMP_MD5).encode()},
        **{'/genomes/{}/{}/assembly_summary.txt'.format(database, domain):
           summary_file(database, domain)
           for database in M.NCBI_DATABASES
           for group in M.NCBI_GROUPS for domain in M.NCBI_GROUP_DOMAINS[group]})

    def setUp(self):
        super().setUp()
        # point the download at the local server instead of NCBI
        self.ncbi_ftp, M.NCBI_FTP = M.NCBI_FTP, self.url
        self.taxdump_url, M.TAXDUMP_URL = M.TAXDUMP_URL, self.url + '/pub/taxonomy/taxdump.tar.gz'

    def tearDown(self):
        M.NCBI_FTP, M.TAXDUMP_URL = self.ncbi_ftp, self.taxdump_url
        super().tearDown()

    RELEASE = 237

    def run_sync(self, group=M.GROUP_PROK):
        M.MetadataSyncManager(self.dir, group).run(self.RELEASE)
        return sorted(os.listdir(self.dir))

    def stamp(self):
        import datetime
        return datetime.date.today().strftime('%Y%m%d')

    def taxonomy(self, name):
        """One file of the standardised taxonomy, as gid -> taxonomy string."""

        with open(self.path('taxonomy', 'standardised_taxonomy', name)) as handle:
            return dict(line.rstrip('\n').split('\t', 1) for line in handle)

    def test_assembly_summaries_land_in_the_root_ncbi_directory(self):
        # this is the directory select_genomes is pointed at
        written = self.run_sync()
        self.assertEqual(
            sorted(n for n in written if n.startswith('assembly_summary_')),
            ['assembly_summary_archaea_genbank.txt.gz',
             'assembly_summary_archaea_refseq.txt.gz',
             'assembly_summary_bacteria_genbank.txt.gz',
             'assembly_summary_bacteria_refseq.txt.gz'])

    def test_a_group_downloads_its_own_summaries_and_no_others(self):
        # prokaryotes and fungi are curated separately, and a run of one must
        # not spend an hour downloading the other's tables
        written = self.run_sync(M.GROUP_FUNGI)
        self.assertEqual(
            sorted(n for n in written if n.startswith('assembly_summary_')),
            ['assembly_summary_fungi_genbank.txt.gz',
             'assembly_summary_fungi_refseq.txt.gz'])

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

    def test_the_second_group_reuses_the_taxonomy_already_extracted(self):
        # both groups need the taxdump; downloading it twice on one day costs
        # 80 MB and, worse, places the two groups against two different
        # downloads of the taxonomy
        self.run_sync()
        for route in ('/pub/taxonomy/taxdump.tar.gz', '/pub/taxonomy/taxdump.tar.gz.md5'):
            del self.server.RequestHandlerClass.routes[route]

        written = self.run_sync(M.GROUP_FUNGI)                # would fail on a download
        self.assertIn('assembly_summary_fungi_refseq.txt.gz', written)
        self.assertTrue(self.taxonomy('ncbi_r237_fungi_unfiltered_taxonomy.tsv'))

    def test_a_taxonomy_downloaded_on_an_earlier_day_is_not_reused(self):
        # a release is built from the taxonomy of the day it is built
        os.makedirs(self.path('taxonomy', 'taxdump_20200101'))
        for dmp in ('names.dmp', 'nodes.dmp'):
            open(self.path('taxonomy', 'taxdump_20200101', dmp), 'w').close()

        self.run_sync()
        self.assertTrue(os.path.exists(
            self.path('taxonomy', 'taxdump_' + self.stamp(), 'names.dmp')))

    def test_the_standardised_taxonomy_is_named_for_the_release_and_the_group(self):
        # both groups are run into one release directory, so a name carrying the
        # release alone would have the second run overwrite the first
        self.run_sync()
        self.run_sync(M.GROUP_FUNGI)
        written = sorted(os.listdir(self.path('taxonomy', 'standardised_taxonomy')))
        taxonomy = [n for n in written if not n.endswith('failed_filters.tsv')]
        self.assertTrue(taxonomy)
        for name in taxonomy:
            self.assertTrue(name.startswith('ncbi_r237_prok_')
                            or name.startswith('ncbi_r237_fungi_'), name)

    def test_neither_group_overwrites_what_the_other_wrote(self):
        self.run_sync()
        prokaryotic = sorted(os.listdir(self.path('taxonomy', 'standardised_taxonomy')))
        self.run_sync(M.GROUP_FUNGI)
        both = sorted(os.listdir(self.path('taxonomy', 'standardised_taxonomy')))

        self.assertTrue(set(prokaryotic).issubset(both))
        self.assertEqual(len(both), 2 * len(prokaryotic))

    def test_the_failed_filter_report_is_written_beside_the_taxonomy(self):
        # it named a hardcoded file, so without a working directory to anchor it
        # the report landed wherever the operator happened to be standing; it is
        # named for its group for the same reason the taxonomy is
        self.run_sync()
        self.assertIn('prok_failed_filters.tsv',
                      os.listdir(self.path('taxonomy', 'standardised_taxonomy')))
        self.run_sync(M.GROUP_FUNGI)
        self.assertIn('fungi_failed_filters.tsv',
                      os.listdir(self.path('taxonomy', 'standardised_taxonomy')))

    def test_the_standardised_taxonomy_resolves_the_lineage_of_each_assembly(self):
        # the point of the step: a taxid becomes a lineage
        self.run_sync()
        lineages = self.taxonomy('ncbi_r237_prok_unfiltered_taxonomy.tsv')

        self.assertEqual(sorted(lineages), ['GCA_000000003.1', 'GCA_000000004.1',
                                            'GCF_000000001.1', 'GCF_000000002.1',
                                            'GCF_000000007.1'])
        self.assertIn('Escherichia coli', lineages['GCF_000000002.1'])
        self.assertIn('Bacteria', lineages['GCF_000000002.1'])
        self.assertIn('Halobacterium salinarum', lineages['GCF_000000001.1'])
        self.assertIn('Archaea', lineages['GCF_000000001.1'])

    def test_a_group_places_its_own_assemblies_and_no_others(self):
        # the taxdump holds every lineage; what a group's taxonomy holds is the
        # assemblies of the summaries that group downloaded
        self.run_sync(M.GROUP_FUNGI)
        lineages = self.taxonomy('ncbi_r237_fungi_unfiltered_taxonomy.tsv')

        self.assertEqual(sorted(lineages), ['GCA_000000006.1', 'GCF_000000005.1'])
        self.assertIn('Candida albicans', lineages['GCF_000000005.1'])

    def test_the_fungal_taxonomy_keeps_the_subranks_ncbi_gives(self):
        # fungal classification leans on the intermediate ranks
        self.run_sync(M.GROUP_FUNGI)
        standardised = self.taxonomy('ncbi_r237_fungi_standardized.tsv')

        taxa = standardised['GCF_000000005.1'].split(';')
        self.assertEqual(len(taxa), 13)
        self.assertIn('sd__Dikarya', taxa)
        self.assertIn('sp__Saccharomycotina', taxa)
        self.assertIn('s__Candida albicans', taxa)

    def test_the_prokaryotic_taxonomy_is_the_seven_ranks_gtdb_curates(self):
        # a subphylum or subclass in it would be a rank GTDB has no name for
        self.run_sync()
        standardised = self.taxonomy('ncbi_r237_prok_standardized.tsv')

        self.assertEqual(len(standardised), 5)
        for gid, taxonomy in standardised.items():
            taxa = taxonomy.split(';')
            self.assertEqual(len(taxa), 7, gid)
            self.assertFalse([t for t in taxa if len(t.split('__')[0]) == 2], gid)

    def test_the_nomenclatural_code_is_not_carried_into_the_taxonomy(self):
        # NCBI appends the code to the name of a taxon it holds under the
        # SeqCode; f__Patescibacteriaceae (SeqCode) matches no other spelling of
        # the same taxon
        self.run_sync()
        lineages = self.taxonomy('ncbi_r237_prok_unfiltered_taxonomy.tsv')

        self.assertNotIn('SeqCode', lineages['GCF_000000007.1'])
        self.assertIn('f__Patescibacteriaceae', lineages['GCF_000000007.1'].split(';'))
        self.assertIn('p__Patescibacteria', lineages['GCF_000000007.1'].split(';'))

    def test_a_seqcode_species_reaches_the_standardised_taxonomy(self):
        # the brackets would otherwise fail the valid character check and drop
        # the genome from the standardised taxonomy entirely
        self.run_sync()
        standardised = self.taxonomy('ncbi_r237_prok_standardized.tsv')

        self.assertEqual(
            standardised['GCF_000000007.1'],
            'd__Bacteria;p__Patescibacteria;c__Patescibacteriia;'
            'o__Patescibacteriales;f__Patescibacteriaceae;g__Patescibacter;'
            's__Patescibacter aquaticus')

    def test_a_corrupt_taxonomy_download_stops_the_run(self):
        # a truncated taxdump reads as a valid, smaller taxonomy, and only shows
        # up as genomes that lost their NCBI lineage several commands later
        self.server.RequestHandlerClass.routes['/pub/taxonomy/taxdump.tar.gz.md5'] = \
            b'0' * 32 + b'  taxdump.tar.gz\n'
        with self.assertRaises(SystemExit):
            M.MetadataSyncManager(self.dir, M.GROUP_PROK).run(self.RELEASE)

    def test_a_missing_file_at_ncbi_stops_the_run(self):
        del self.server.RequestHandlerClass.routes['/genomes/genbank/bacteria/assembly_summary.txt']
        with self.assertRaises(SystemExit):
            M.MetadataSyncManager(self.dir, M.GROUP_PROK).run(self.RELEASE)

    def test_an_unknown_group_is_refused_before_anything_is_downloaded(self):
        with self.assertRaises(ValueError):
            M.MetadataSyncManager(self.dir, 'prok')
        self.assertEqual(os.listdir(self.dir), [])

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
        self.assertEqual(rows, [('GCF_000000002.1', '562'),
                                ('GCF_000000007.1', '2802431')])

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
        refseq, genbank = S.SelectedGenomesManager(self.dir).group_by_database(summaries)
        self.assertEqual(len(refseq), len(M.NCBI_PROK_DOMAINS))
        self.assertEqual(len(genbank), len(M.NCBI_PROK_DOMAINS))
