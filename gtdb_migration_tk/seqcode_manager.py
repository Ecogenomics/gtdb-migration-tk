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

"""download_seqcode_data: the SeqCode Registry's species, by their type genomes.

It was generate_seqcode_table, Tools.generate_seqcode_table() of utils/tools.py,
until 0.1.72.

THE REGISTRY
The Registry's API lists every species with a type genome (type-genomes.json,
50 a page; per_page is not honoured), each with its classification; a species'
own record is fetched for the publication that proposed it (proposed_in), which
the list does not give. The API moved from disc-genomics.uibk.ac.at/seqcode,
which no longer answers, to api.seqco.de/v1 (SEQCODE_API). In October 2026 it
listed 1,816 species in 37 pages: 1,750 typed by an assembly and 66 by a
nucleotide accession, 26 of a WGS project and 40 of a complete sequence (CP, AP).
The list is read once, where it was read twice, and a species' record is fetched
only where its type genome is in the release, where it was fetched for every
species and the row then dropped.

A species' record takes about a second to fetch, some 35 minutes for 1,743 one
after another, so SPECIES_FETCHERS are fetched at once (8 minutes for 1,743 in
October 2026; r237 types 1,783), and each citation is kept in a cache
(--species_cache, by default CACHE_NAME in --out_dir) with the updated_at of the
species' list entry: a species whose entry has not
changed since is not fetched again, so a run again, or the next release given the
same cache, fetches only the species new or changed. The cache is written as the
run ends, whether or not it ends well, so a run stopped by a request that failed
keeps what it fetched.

Every request is retried (FETCH_ATTEMPTS, waiting longer each time) and a
request that still fails ends the run, naming the URL, with nothing written: a
table missing a page of species, or a species' citation, would load as though it
were the Registry's. The table is written under a temporary name until whole.

WHICH GENOME IS A SPECIES' TYPE
The genome_dirs file (-g) says which genomes are the release's; only its
accessions are read, and no file in a genome's directory. An assembly is matched
to the release's genome by its canonical accession
(biolib_lite.common.canonical_gid()). A WGS accession (JBBJJH000000000.1) is
matched by its project prefix (JBBJJH) to the assembly whose wgs_master the
release's NCBI assembly summaries (-n) give it, read by column name: 1,268,107 of
r237's genomes have one. Any other nucleotide accession -- a complete sequence,
CP of GenBank or AP of DDBJ, 40 of 1,816 in October 2026 -- is asked of NCBI's
E-utilities (EUTILS): esummary of nuccore for the sequence's uid, elink from it
to assembly, and esummary of assembly for the GenBank and RefSeq accessions,
each request taking EUTILS_BATCH sequences, so three requests for the 40. The
Registry's record of a type genome names its sequence and no assembly.

The sequences were matched through the GenBank-Accn of every genome's assembly
report, 1.35M small files for r237 read on --cpus processes, about two hours for
the 40, and the only use the command made of -c and of the genome directories.
Those processes were started and joined by hand, and a failed one passed over
with the run carrying on with what the others had collected; only CP accessions
were kept, so AP types were never found; and a WGS prefix the release does not
hold ended the run with a KeyError. NCBI's E-utilities ask three requests a
second of a caller with no API key, which three requests one after another do
not reach; one refused (429) is tried again as any other. An answer without its
result ends the run as a request that failed does, since a type passed over
would leave its species out of the table unsaid. A sequence NCBI does not know,
or links to no assembly, is a type not in the release. Each type not in the
release is listed in NOT_IN_RELEASE_NAME, counted in the log.

Nothing is printed: the log says how far the run is, with a bar on a terminal
as the Registry and species are read. It printed every species' record, as JSON,
and every row it wrote.
"""

import json
import logging
import os
import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import requests
from tqdm import tqdm

from gtdb_migration_tk.biolib_lite.common import canonical_gid
from gtdb_migration_tk.ncbi_utils import NCBI_NA, read_summary_rows, summary_field

SEQCODE_API = 'https://api.seqco.de/v1'
TYPE_GENOMES_URL = SEQCODE_API + '/type-genomes.json'

# NCBI's E-utilities, which give the assembly of a complete sequence; tool names
# the caller, as NCBI asks, and a request takes EUTILS_BATCH ids
EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
EUTILS_TOOL = 'gtdb_migration_tk'
EUTILS_BATCH = 100

# each request is tried this often, waiting FETCH_WAIT seconds, doubled each
# time, between tries, and given FETCH_TIMEOUT seconds to answer
FETCH_ATTEMPTS = 4
FETCH_WAIT = 5
FETCH_TIMEOUT = 60

# species' records fetched at once
SPECIES_FETCHERS = 4

# what is written in --out_dir
TABLE_NAME = 'seqcode_table.tsv'
CACHE_NAME = 'seqcode_species_cache.json'
NOT_IN_RELEASE_NAME = 'seqcode_types_not_in_release.tsv'
NOT_IN_RELEASE_HEADER = ('seqcode_id', 'seqcode_name', 'type_material')

RANK_ORDER = ('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species')
RANK_PREFIXES = ('d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__')

# the table's columns, each seqcode_<label>, and the field of the species' list
# entry it is from; the ones derived from it are named for what they say
SEQCODE_FIELDS = (('type_material_accn', 'nomenclatural_type'),
                  ('id', 'id'),
                  ('name', 'name'),
                  ('rank', 'rank'),
                  ('species_status', 'status_name'),
                  ('priority_date', 'priority_date'),
                  ('genus_status', 'genus_status'),
                  ('family_status', 'family_status'),
                  ('order_status', 'order_status'),
                  ('class_status', 'class_status'),
                  ('phylum_status', 'phylum_status'),
                  ('type_species_of_genus', 'type_species_of_genus'),
                  ('type_genus_of_family', 'type_genus_of_family'),
                  ('type_genus_of_order', 'type_genus_of_order'),
                  ('type_genus_of_class', 'type_genus_of_class'),
                  ('type_genus_of_phylum', 'type_genus_of_phylum'),
                  ('classification', 'classification'),
                  ('proposed_in', 'proposed_in'),
                  ('created_at', 'created_at'),
                  ('updated_at', 'updated_at'),
                  ('url', 'url'))
TABLE_HEADER = tuple('seqcode_' + label for label, _ in SEQCODE_FIELDS)

# a WGS accession: a four or six letter project prefix, then its digits
WGS_ACCESSION = re.compile(r'^([A-Z]{4}([A-Z]{2})?)([0-9]{6,})$')


class SeqCodeError(RuntimeError):
    """The Registry or the release cannot be read; nothing is written."""


def fetch_json(url: str, params: Optional[Sequence[Tuple[str, str]]] = None, attempts: int = FETCH_ATTEMPTS,
               wait: float = FETCH_WAIT, timeout: float = FETCH_TIMEOUT, sleep=time.sleep):
    """A JSON document of the Registry or of NCBI, tried again where it is not answered.

    @return: the parsed document.

    Raises
    ------
    SeqCodeError
        Every attempt failed: no answer, an HTTP error, or what is not JSON.
    """

    problem = None
    for attempt in range(attempts):
        if attempt:
            sleep(wait * 2 ** (attempt - 1))
        try:
            response = requests.get(url, params=params, timeout=timeout)
            response.raise_for_status()
            return response.json()
        except (requests.RequestException, ValueError) as exc:
            problem = exc
    raise SeqCodeError('{} could not be read in {} attempt(s): {}'.format(url, attempts, problem))


def eutils(tool: str, params: Sequence[Tuple[str, str]], part: str):
    """One E-utilities request, as JSON.

    Parameters
    ----------
    tool : str
        esummary or elink.
    params : sequence of (name, value)
        The request's parameters; a name may be given more than once.
    part : str
        What the answer holds when it is one: result or linksets.

    @return: that part of the answer.

    Raises
    ------
    SeqCodeError
        The request failed every attempt, or NCBI answered without that part.
    """

    url = '{}/{}.fcgi'.format(EUTILS, tool)
    document = fetch_json(url, params=list(params) + [('retmode', 'json'), ('tool', EUTILS_TOOL)])
    if not isinstance(document, dict) or part not in document:
        problem = document.get('error') or document.get('ERROR') if isinstance(document, dict) else None
        raise SeqCodeError('{} answered without its {}: {}'.format(url, part, problem or document))
    return document[part]


def batches(items: Sequence[str]) -> List[Sequence[str]]:
    """items in runs of EUTILS_BATCH.

    @return: the runs.
    """

    return [items[i:i + EUTILS_BATCH] for i in range(0, len(items), EUTILS_BATCH)]


def sequence_assemblies(accessions: Iterable[str]) -> Dict[str, Tuple[str, ...]]:
    """The assemblies NCBI gives each nucleotide sequence.

    Parameters
    ----------
    accessions : iterable of str
        Sequence accessions as the Registry gives them, versioned or not; one
        without its version is NCBI's latest.

    @return: accession as given -> the GenBank and RefSeq accessions of the
             assemblies it is in; a sequence NCBI does not know, or links to no
             assembly, is absent.

    Raises
    ------
    SeqCodeError
        A request failed every attempt, or NCBI answered without its result.
    """

    accessions = sorted(set(accessions))
    uids = {}
    for batch in batches(accessions):
        result = eutils('esummary', [('db', 'nuccore'), ('id', ','.join(batch))], 'result')
        # a versioned accession is its own record; an unversioned one NCBI's
        # latest, the highest version of its caption where a batch asks for two
        held, latest = {}, {}
        for uid in result.get('uids') or ():
            document = result.get(uid) or {}
            version = document.get('accessionversion')
            if version:
                held[version] = str(uid)
                caption = document.get('caption') or version.split('.')[0]
                number = int(version.split('.')[1]) if version.split('.')[-1].isdigit() else 0
                if caption not in latest or number > latest[caption][0]:
                    latest[caption] = (number, str(uid))
        held.update((caption, uid) for caption, (_, uid) in latest.items() if caption not in held)
        uids.update((accession, held[accession]) for accession in batch if accession in held)

    links: Dict[str, set] = {}
    for batch in batches(sorted(set(uids.values()))):
        # an id each, for a linkset each
        for linkset in eutils('elink', [('dbfrom', 'nuccore'), ('db', 'assembly')] +
                              [('id', uid) for uid in batch], 'linksets'):
            for linksetdb in linkset.get('linksetdbs') or ():
                if linksetdb.get('linkname') == 'nuccore_assembly':
                    for uid in linkset.get('ids') or ():
                        links.setdefault(str(uid), set()).update(str(a) for a in linksetdb.get('links') or ())

    assemblies: Dict[str, set] = {}
    for batch in batches(sorted(set().union(*links.values()))):
        result = eutils('esummary', [('db', 'assembly'), ('id', ','.join(batch))], 'result')
        for uid in result.get('uids') or ():
            document = result.get(uid) or {}
            synonym = document.get('synonym') or {}
            names = {document.get('assemblyaccession'), synonym.get('genbank'), synonym.get('refseq')}
            assemblies[str(uid)] = {name for name in names if name}

    found = {}
    for accession, uid in uids.items():
        names = set().union(*(assemblies.get(a, set()) for a in links.get(uid, ())))
        if names:
            found[accession] = tuple(sorted(names))
    return found


def wgs_prefix(accession: str) -> Optional[str]:
    """The project prefix of a WGS accession or master, e.g. JBBJJH of JBBJJH000000000.1.

    @return: the prefix, or None where the accession is not a WGS one.
    """

    match = WGS_ACCESSION.match(accession.split('.')[0])
    return match.group(1) if match else None


def wgs_index(assembly_summaries: Sequence[str], canonical_ids: Dict[str, str]) -> Dict[str, str]:
    """Each WGS project prefix of the summaries' wgs_master, to the release's genome.

    A genome is matched by canonical accession, so the GenBank genome of a release
    is found by its RefSeq twin's row and the other way round.

    @return: e.g. {'JBBJJH': 'GCA_...'}.

    Raises
    ------
    ncbi_utils.BadInput
        A summary without the assembly_accession or wgs_master column.
    """

    index = {}
    for summary in assembly_summaries:
        for _, fields, columns in read_summary_rows(summary, required=('assembly_accession', 'wgs_master')):
            prefix = wgs_prefix(summary_field(fields, columns, 'wgs_master') or NCBI_NA)
            genome = canonical_ids.get(canonical_gid(summary_field(fields, columns, 'assembly_accession')))
            if prefix and genome:
                index[prefix] = genome
    return index


def read_cache(path: str) -> Dict[str, dict]:
    """The species' citations of an earlier run, by species id.

    @return: id -> {'updated_at': ..., 'proposed_in': citation}; empty where there
             is no cache, or one that cannot be read, which is fetched again.
    """

    try:
        with open(path) as handle:
            cache = json.load(handle)
    except (OSError, ValueError):
        return {}
    return cache if isinstance(cache, dict) else {}


def write_cache(path: str, cache: Dict[str, dict]) -> None:
    """Write the cache, under a temporary name until whole.

    @return: None
    """

    partial = path + '.partial'
    with open(partial, 'w') as handle:
        json.dump(cache, handle, sort_keys=True)
    os.replace(partial, path)


def type_genome(nomenclatural_type: Optional[dict], canonical_ids: Dict[str, str],
                index: Dict[str, str]) -> Tuple[Optional[str], str]:
    """The release's genome a species is typed by.

    @return: (genome accession or None, the type material as the Registry
             gives it).
    """

    nomenclatural_type = nomenclatural_type or {}
    if nomenclatural_type.get('assembly'):
        assembly = str(nomenclatural_type['assembly'])
        return canonical_ids.get(canonical_gid(assembly)), assembly
    if nomenclatural_type.get('nuccore'):
        accession = str(nomenclatural_type['nuccore'])
        for key in (accession, accession.split('.')[0], wgs_prefix(accession)):
            if key and key in index:
                return index[key], accession
        return None, accession
    return None, str(nomenclatural_type.get('display') or nomenclatural_type or 'none')


def classification(ranks: Sequence[dict], species_id) -> Tuple[str, List[Optional[str]], List[bool]]:
    """A species' classification, the status of each rank, and the types among them.

    @return: (d__...;s__..., the status of each rank in RANK_ORDER, whether the
             genus is the type of each rank above it and the species the type of
             its genus).
    """

    names = list(RANK_PREFIXES)
    statuses: List[Optional[str]] = [None] * len(RANK_ORDER)
    type_ids = [0] * len(RANK_ORDER)
    genus_id = None
    for rank in ranks or ():
        if rank.get('rank') in RANK_ORDER:
            i = RANK_ORDER.index(rank['rank'])
            names[i] += rank.get('name') or ''
            statuses[i] = rank.get('status_name')
            type_ids[i] = int((rank.get('nomenclatural_type') or {}).get('id') or 0)
        if rank.get('rank') == 'genus':
            genus_id = int(rank.get('id') or 0)

    types = [False] * len(RANK_ORDER)
    for i, type_id in enumerate(type_ids[:-2]):
        types[i] = genus_id is not None and type_id == genus_id
    types[-2] = species_id == type_ids[-2]
    return ';'.join(names), statuses, types


def species_row(record: dict, genome: str, proposed_in: Optional[str]) -> List[str]:
    """A species' row of the table, in SEQCODE_FIELDS order; a value the Registry
    does not give is empty, where it was the text 'None'.

    @return: the row.
    """

    taxonomy, statuses, types = classification(record.get('classification'), record.get('id'))
    row = []
    for _, field in SEQCODE_FIELDS:
        if field == 'nomenclatural_type':
            value = genome
        elif field == 'proposed_in':
            value = proposed_in
        elif field == 'classification':
            value = taxonomy
        elif field.endswith('_status') and field != 'status_name':
            value = statuses[RANK_ORDER.index(field.split('_')[0])]
        elif field.startswith('type_'):
            value = types[RANK_ORDER.index(field.split('_')[-1])]
        else:
            value = record.get(field)
        row.append('' if value is None else str(value))
    return row


def write_tsv(path: str, header: Sequence[str], rows: Iterable[Sequence[str]]) -> int:
    """Write a table, under a temporary name until it is whole.

    @return: the number of rows written.
    """

    partial = path + '.partial'
    written = 0
    try:
        with open(partial, 'w') as handle:
            handle.write('\t'.join(header) + '\n')
            for row in rows:
                handle.write('\t'.join(row) + '\n')
                written += 1
    except BaseException:
        if os.path.exists(partial):
            os.remove(partial)
        raise
    os.replace(partial, path)
    return written


class SeqCodeManager(object):
    """Download the SeqCode Registry's species and the release's genome typing each."""

    def __init__(self, fetchers: int = SPECIES_FETCHERS):
        self.logger = logging.getLogger('timestamp')
        self.fetchers = max(1, fetchers)

    def bar(self, **kwargs):
        """A progress bar on a terminal, and none with --silent."""
        return tqdm(ncols=100, leave=False, disable=True if getattr(self.logger, 'is_silent', False) else None,
                    **kwargs)

    def type_genome_records(self) -> List[dict]:
        """Every species the Registry lists with a type genome.

        @return: their list entries, page by page.
        """

        first = fetch_json(TYPE_GENOMES_URL)
        pages = int(first['response']['total_pages'])
        records = list(first.get('values') or ())
        with self.bar(total=pages, initial=1, desc='Reading the Registry', unit=' pages') as progress:
            for page in range(2, pages + 1):
                records.extend(fetch_json('{}?page={}'.format(TYPE_GENOMES_URL, page)).get('values') or ())
                progress.update()
        self.logger.info('The SeqCode Registry lists {:,} type genome(s) in {:,} page(s) ({:,} said).'.format(
            len(records), pages, int(first['response'].get('count') or 0)))
        return records

    def sequence_index(self, accessions: Sequence[str], canonical_ids: Dict[str, str]) -> Dict[str, str]:
        """Each complete sequence typing a species, to the release's genome holding it.

        @return: accession as given -> genome; one in no genome of the release is absent.
        """

        started = time.time()
        found = sequence_assemblies(accessions)
        index = {}
        for accession, assemblies in found.items():
            for assembly in assemblies:
                genome = canonical_ids.get(canonical_gid(assembly))
                if genome:
                    index[accession] = genome
                    break
        self.logger.info('NCBI gives the assembly of {:,} of the {:,} complete sequence(s), {:,} of them in the '
                         'release, in {:.0f} s.'.format(len(found), len(accessions), len(index),
                                                        time.time() - started))
        return index

    def citations(self, records: Sequence[dict], cache_file: str) -> Dict[str, Optional[str]]:
        """The publication proposing each species, from the cache or its record.

        Records are fetched self.fetchers at once; the cache is written as this
        ends, whether or not it ends well.

        @return: species id -> citation, or None where its record gives none.

        Raises
        ------
        SeqCodeError
            A species' record could not be fetched.
        """

        cache = read_cache(cache_file)
        citations, to_fetch = {}, []
        for record in records:
            held = cache.get(str(record.get('id')))
            if held is not None and held.get('updated_at') == record.get('updated_at'):
                citations[str(record.get('id'))] = held.get('proposed_in')
            else:
                to_fetch.append(record)
        self.logger.info('{:,} species\' citation(s) are in {} unchanged; {:,} to fetch, {:,} at a time.'.format(
            len(citations), cache_file, len(to_fetch), self.fetchers))

        executor = ThreadPoolExecutor(max_workers=self.fetchers)
        futures = {}
        try:
            with self.bar(total=len(to_fetch), desc='Fetching species', unit=' species') as progress:
                futures = {executor.submit(fetch_json, record['url']): record for record in to_fetch}
                for future in as_completed(futures):
                    record = futures[future]
                    citation = (future.result().get('proposed_in') or {}).get('citation')
                    citations[str(record.get('id'))] = citation
                    cache[str(record.get('id'))] = {'updated_at': record.get('updated_at'), 'proposed_in': citation}
                    progress.update()
        finally:
            # those not yet started are not (shutdown's cancel_futures is Python 3.9)
            for future in futures:
                future.cancel()
            executor.shutdown(wait=True)
            write_cache(cache_file, cache)
        return citations

    def run(self, gtdb_genome_path_file: str, out_dir: str, assembly_summaries: Sequence[str],
            cache_file: Optional[str] = None) -> None:
        """Write TABLE_NAME: each SeqCode species whose type genome is in the release.

        Parameters
        ----------
        gtdb_genome_path_file : str
            genome_dirs file of the release; only its accessions are read.
        out_dir : str
            Directory the table and the lists are written to.
        assembly_summaries : sequence of str
            The release's NCBI assembly summaries (-n), for each assembly's WGS project.
        cache_file : str
            The species' citations of earlier runs (--species_cache); by default
            CACHE_NAME in out_dir.

        @return: None

        Raises
        ------
        SeqCodeError
            A request to the Registry or to NCBI failed every attempt, or NCBI
            answered without its result; no table is written.
        ncbi_utils.BadInput
            A summary without the assembly_accession or wgs_master column.
        """

        # the release's genomes; nothing in their directories is read
        genomes = []
        with open(gtdb_genome_path_file) as handle:
            for line in handle:
                accession = line.rstrip('\n').split('\t')[0]
                if accession:
                    genomes.append(accession)
        canonical_ids = {canonical_gid(accession): accession for accession in genomes}
        self.logger.info('Read {:,} genome(s) from {}.'.format(len(genomes), gtdb_genome_path_file))

        records = [r for r in self.type_genome_records() if r.get('rank') == 'species']
        sequences = [str(r['nomenclatural_type']['nuccore']) for r in records
                     if (r.get('nomenclatural_type') or {}).get('nuccore')]
        by_wgs = sum(1 for accession in sequences if wgs_prefix(accession))
        self.logger.info('{:,} species are typed by an assembly, {:,} by a WGS project and {:,} by a complete '
                         'sequence.'.format(len(records) - len(sequences), by_wgs, len(sequences) - by_wgs))

        index = {}
        if by_wgs:
            started = time.time()
            index.update(wgs_index(assembly_summaries, canonical_ids))
            self.logger.info('Read the WGS project of {:,} genome(s) of the release from the assembly summaries '
                             'in {:.0f} s.'.format(len(index), time.time() - started))
        complete = sorted({accession for accession in sequences if not wgs_prefix(accession)})
        if complete:
            index.update(self.sequence_index(complete, canonical_ids))

        matched, not_found = [], []
        for record in records:
            genome, material = type_genome(record.get('nomenclatural_type'), canonical_ids, index)
            if genome is None:
                not_found.append((str(record.get('id')), str(record.get('name')), material))
            else:
                matched.append((record, genome))

        path = os.path.join(out_dir, NOT_IN_RELEASE_NAME)
        write_tsv(path, NOT_IN_RELEASE_HEADER, not_found)
        self.logger.info('{:,} of the {:,} species are typed by a genome the release does not hold; each is '
                         'listed in {}.'.format(len(not_found), len(records), path))

        citations = self.citations([record for record, _ in matched], cache_file or os.path.join(out_dir, CACHE_NAME))
        rows = [species_row(record, genome, citations.get(str(record.get('id')))) for record, genome in matched]
        path = os.path.join(out_dir, TABLE_NAME)
        write_tsv(path, TABLE_HEADER, rows)
        self.logger.info('Wrote {:,} species to {}.'.format(len(rows), path))
