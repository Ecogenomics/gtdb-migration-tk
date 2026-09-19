# gtdb-migration-tk
[![ace-internal](https://img.shields.io/conda/vn/ace-internal/gtdb_migration_tk.svg?color=green)](https://anaconda.org/ace-internal/gtdb_migration_tk)

Toolkit for updating the [GTDB](https://gtdb.ecogenomic.org/) to the next release.

## Overview

`gtdb-migration-tk` automates the steps required to move the GTDB from one NCBI
release to the next. It is an internal tool: it assumes access to the GTDB
directory structure and, for many commands, to the GTDB PostgreSQL database.

The toolkit covers the full migration cycle:

* **Mirror NCBI** — sync RefSeq/GenBank genome assemblies to a local directory
  and move new or modified genomes into the GTDB folder layout.
* **Annotate genomes** — call genes with Prodigal, search Pfam/TIGRFAM markers,
  identify rRNA genes against SILVA and LTP, find tRNAs, and estimate genome
  quality with CheckM/CheckM2/BUSCO.
* **Build metadata** — derive nucleotide and gene statistics, parse NCBI
  assembly summaries and taxonomy dumps, and load the results into the database.
* **Propagate taxonomy** — carry GTDB taxonomy forward to the new release and
  push the result back to the database.
* **Nomenclature** — pull and parse information from LPSN, BacDive and SeqCode.
* **Validate** — compare metadata between releases and check the database is
  populated as expected.

Every step is a subcommand of a single `gtdb_migration_tk` executable.

## Installation

### Requirements

* Python >= 3.8
* PostgreSQL access (only for the `*_db` commands)
* External tools on `$PATH`, required by the commands that call them:

  | Tool | Used by |
  | --- | --- |
  | `prodigal` | `prodigal`, `trans_table` (via gTranslate) |
  | `gtranslate` | `trans_table` |
  | `hmmsearch` | `hmmsearch`, `top_hit` |
  | `blastn`, `blastp`, `makeblastdb` | `rna_silva`, `rna_ltp`, `generate_ltp_db` |
  | `nhmmer` | `rna_silva`, `rna_ltp` |
  | `tRNAscan-SE` | `trnascan` |
  | `busco` | `busco` |

  CheckM and CheckM2 are run separately; the toolkit consumes their output.

Python dependencies are installed automatically: `requests`, `unidecode`,
`pandas`, `numpy`, `sqlalchemy`, `beautifulsoup4`, `dendropy`, `tqdm`,
`atpbar`, `python-dateutil` and `psycopg2-binary`.

### From conda

```bash
conda install -c ace-internal gtdb_migration_tk
```

### From source

```bash
git clone https://github.com/Ecogenomics/gtdb-migration-tk.git
cd gtdb-migration-tk
pip install .
```

For development, install in editable mode with the test dependencies:

```bash
pip install -e ".[test]"
```

## Usage

Run the executable with no arguments to list the available commands:

```bash
gtdb_migration_tk
```

Each command has its own help:

```bash
gtdb_migration_tk <command> -h
```

The toolkit can also be invoked as a module, which is equivalent:

```bash
python -m gtdb_migration_tk <command>
```

Most commands take `-l/--log` and write a run log there; without it the log is
written to `./gtdb_migration_tk.log`. Add `--silent` to suppress console output.

### Example: mirroring NCBI

`ncbi_genome_sync` keeps a mirror equal to the table `select_genomes` wrote: every
genome directory under `--root` that the table does not list is removed, then every
genome it lists is fetched or brought up to date, each file verified against the
`md5checksums.txt` that NCBI publishes alongside it:

```bash
S=/srv/db/gtdb/metadata/release237/ncbi/gtdb_selected_genomes.tsv.gz

# see what a run would do: how many directories removed, genomes added, genomes
# already present -- and the removal list in <log>.rm_dry_run. Changes nothing.
gtdb_migration_tk ncbi_genome_sync --gtdb_selected_genomes $S \
    --root /srv/db/gtdb/genomes -l ./logs/sync.log --dry-run

# remove what the selection does not list, then sync what it does
gtdb_migration_tk ncbi_genome_sync --gtdb_selected_genomes $S \
    --root /srv/db/gtdb/genomes -l ./logs/sync.log

# retry only the genomes that failed last time; a retry removes nothing it was not given.
# The round gets its own --log, and so its own outputs: it reads sync.fail and writes
# retry1.fail. Reusing -l ./logs/sync.log here is refused, since the run would truncate
# the file it is reading.
gtdb_migration_tk ncbi_genome_sync --retry ./logs/sync.fail \
    --root /srv/db/gtdb/genomes -l ./logs/retry1.log

# rebuild the genomes that failed verification: --delete removes each directory
# immediately before that genome is refetched, in the one run
gtdb_migration_tk ncbi_genome_sync --retry ./logs/sync.bad \
    --root /srv/db/gtdb/genomes -l ./logs/rebuild.log --delete
```

`-l/--log` names the log and the run alike: every file a run writes is that name with
`.log` stripped, beside the log — `-l ./logs/sync.log` gives `./logs/sync.fail`,
`sync.bad`, `sync.rm`, `sync.rm_dry_run` and `sync.extra` (`<log>.*` below). Any other
extension is kept, so `run.txt` yields `run.txt.fail`. Give each round its own log and
one round's outputs can never overwrite another's; a run that would truncate its own
`--retry` input is refused (exit 2). `--fail`/`--bad` override the two that are also
inputs. A failure file carries the
same columns as the input, so it can be fed back in with `--retry`. Never give a `.fail`
or any subset of the table to `--gtdb_selected_genomes`: the table defines what the
mirror should hold, and everything it leaves out would be removed. An empty table is
refused for that reason.

A genome NCBI lists but serves no directory for (`ftp_path` is `na`) is not selected from
either database — there would be nothing to mirror. Because RefSeq covers a genome only
when the RefSeq assembly was itself selected, an unserved RefSeq genome does not suppress
its GenBank counterpart: that copy is selected instead, with the `notes` column naming the
RefSeq assembly it stands in for and why. Without this the genome would leave the release
entirely, even though NCBI holds it under its GenBank accession.

`--verify` and `--verify-only` check both halves of "the mirror equals the selection":
every listed genome present and md5-clean (failures to `<log>.bad`), and nothing else
present (directories the selection does not list to `<log>.extra`). Either fails the
verification; `--delete` removes both. Given `--retry`, only the listed genomes are
verified.

NCBI publishes two checksum tables per assembly, `md5checksums.txt` and — for some
assemblies — `uncompressed_checksums.txt`, and they can contradict each other: about 1000
genomes have a stale `md5checksums.txt` entry for `<asm>_fcs_report.txt` while the other
table is right. Believing the first alone means those genomes can never sync at all (the
bytes NCBI serves do not match the checksum NCBI publishes for them, so the download fails
every run) and verify reports rot that is not there. So a file `md5checksums.txt` rejects
gets a second opinion from `uncompressed_checksums.txt`, in both the sync and the verify,
and passes if that table vouches for exactly those bytes under exactly that name. A file
neither table vouches for fails as before. Only uncompressed files are eligible: that table
lists the uncompressed *form* of everything, so a `.gz` appears in it under a name it does
not have on disk and can never match — it is refused before a request is made, and can
never excuse a corrupt archive. The second table is fetched only on an eligible
disagreement — one request per affected genome, none for a clean one — and is written into
that genome's directory as the record of why the bytes were accepted, so some genome
directories hold one and most do not; its absence is never a fault. Each run's closing
lines report how many files and genomes it settled.

A `<log>.bad` row carries a `failed_files` column naming every file of that genome at
fault, comma-separated, each tagged with which fault it was — `missing:NAME` (the manifest
lists it, the mirror does not have it), `mismatch:NAME` (present, bytes rotted) or
`unreadable:NAME`. The whole genome is checked before the row is written, so the column is
every fault rather than the first one found. A fault with no file to name — missing
directory, no `md5checksums.txt` — puts that reason in the column instead. The four columns
before it are what `--retry` reads, so the file still feeds straight back in.

`--delete` means "remove the directory rather than repair it". With `--retry` it rebuilds:
each listed genome's directory is removed immediately before that genome is fetched, which
is what a `<log>.bad` needs, since a plain sync trusts an intact manifest and would repair
nothing. With `--verify`/`--verify-only` it removes what fails verification and what the
selection does not list. A plain `--gtdb_selected_genomes` sync refuses it (exit 2) — there
it could only mean re-downloading the whole mirror.

`ncbi_genome_sync` returns meaningful exit codes so it can be driven from a wrapper
script:

| Code | Meaning |
| --- | --- |
| `0` | everything synced / verified clean |
| `1` | some genomes failed — see `<log>.fail` / `<log>.bad` — or verification found directories the selection does not list — see `<log>.extra` |
| `2` | usage error, or a malformed assembly summary |
| `74` | filesystem refused the write (disk full, quota, read-only) |
| `75` | NCBI is throttling this host, or another sync holds `--root` — retry later |
| `130` / `143` | interrupted (SIGINT / SIGTERM); in-flight genomes recorded in `<log>.fail` |

## Commands

Run `gtdb_migration_tk <command> -h` for the arguments of any command.

### Mirror and update genome directories

| Command | Description |
| --- | --- |
| `ncbi_metadata_sync` | Download the NCBI taxonomy and the RefSeq and GenBank assembly summary files one group (`--group PROK` or `FUNGI`) is built from, and generate its standardised NCBI taxonomy |
| `ncbi_genome_sync` | Sync NCBI data to a local directory |
| `select_genomes` | Select the NCBI genomes which will comprise the new GTDB release |
| `update_genomes` | Update RefSeq and GenBank genomes from the NCBI FTP mirror, carrying derived data across where the genome's sequences are unchanged |
| `list_genomes` | Produce file indicating the directory of each genome |

`update_genomes` writes the new release under `--new_directory`, RefSeq and
GenBank in trees of their own and each genome under NCBI's own nesting:

```
<new_directory>/
  refseq/GCF/000/006/805/GCF_000006805.1_ASM680v1/
  genbank/GCA/047/639/395/GCA_047639395.1_ASM4763939v1/
  report.log          the fate of every genome of the release
  to_review.log       genomes needing manual attention
  genome_dirs.tsv     the genome_dirs file of the new release
```

`genome_dirs.tsv` is in the format `list_genomes` writes and every later command
reads, so the new release does not need indexing with `list_genomes` afterwards;
that command is for indexing the mirror. It names only the genomes whose
directory was written, so a `--dry_run`, which writes none, writes no
`genome_dirs.tsv` either.

`--fresh` starts the release from the NCBI genome data alone: every genome of the
mirror is copied into the new release and reported as `new`, nothing is compared
to the previous release and no derived data is carried across, so everything
derived from the genomes is to be regenerated. It reads no
`--old_genome_dirs_file`, which is therefore required only without it.

`--resume` continues a run that was interrupted -- the disk filled, the job hit a
wall clock -- by reading the `genome_dirs.tsv` it left in `--new_directory` and
handling again only the genomes it does not name. A genome is written to that
file once its directory is there, so the genomes the run failed on or never
reached are exactly the ones missing from it, and the directories it was in the
middle of writing are replaced rather than refused. The reports of the
interrupted run are carried forward for the genomes it finished and rewritten for
the rest, so `report.log` describes the release once through rather than the
fragment the second run happened to do. Without `--resume` an interrupted run
cannot simply be repeated: every genome would be copied again, and the copy of
the first genome already in the release would fail rather than overwrite it.
Combined with `--dry_run` it reports what is left without building any of it, and
leaves the `genome_dirs.tsv` it read untouched. It continues a run, so it expects
the inputs of the run it continues: a genome the mirror has stopped offering
since keeps the directory and the report row the interrupted run gave it.

`--cpus` is the number of genomes compared, and copied, at once. It defaults to
16. Neither half of the work is bound by the CPU: copying a genome is round trips to the file server, and the comparison only
reads a FASTA where NCBI has reissued one. Measured against a mirror on NFS, 200
genomes of the size NCBI serves:

| `--cpus` | copying (`--fresh`) | comparing, MD5s agree | comparing, every FASTA reissued |
| --- | --- | --- | --- |
| 1 | 17.1s | 15.5s | 46.0s |
| 8 | 5.2s | 4.3s | 6.5s |
| 16 | 4.5s | 3.7s | 4.3s |
| 32 | 4.2s | 3.4s | 3.9s |
| 64 | 4.0s | -- | 4.2s |

Hence the default of 16. Raise it to 32 for a release NCBI has reissued heavily;
past that the copying is bound by the link rather than by how many genomes are in
flight, and the hashing turns back down. There is no reason to scale with the
core count: 64 was slower than 32 on the one part of this that is CPU work. Lower
it on a file server shared with other work.

### Gene calling and annotation

| Command | Description |
| --- | --- |
| `trans_table` | Predict the translation table of each genome using gTranslate |
| `prodigal` | Call genes using Prodigal, under the translation table `trans_table` predicted |
| `hmmsearch` | Run HMMER on new and modified genomes |
| `top_hit` | Generate TopHit file for TIGRFAM or Pfam |
| `metadata` | Generate metadata derived from nucleotide and protein files |
| `rna_silva` | Identify and classify 16S, 23S and 5S rRNA genes against SILVA |
| `rna_ltp` | Identify and classify 16S rRNA genes against LTP |
| `update_silva` | Update taxonomy files and BLAST database from the latest SILVA release |
| `generate_ltp_db` | Generate BLAST database from the LTP website |
| `trnascan` | Identify tRNAs in genomes |

`trans_table` runs gTranslate over the genomes of a release in batches of
`--batch_size` (default 10,000), each batch a directory of its own under
`--out_dir`:

```
<out_dir>/
  batch_000001/
    gtranslate_batchfile.tsv            the genomes of this batch
    RUNNING                             a machine is working on it (host, PID, time)
    PREDICTED                           gTranslate has run; its results are final
    SUCCESS                             it finished; its results are complete
    FAILED                              it was attempted and something went wrong
    FAILED.20260919T140328              what an earlier attempt said, kept
    trans_table.log                     what this command did to this batch
    gtranslate.log                      what gTranslate did to it
    no_prediction.tsv                   genomes gTranslate returned nothing for
    gtranslate.translation_table_summary.tsv
    ncbi_tt_comparison.tsv              this batch, compared against NCBI
  batch_000002/
  ncbi_tt_comparison.tsv                the whole release, once every batch has SUCCESS
  gtranslate.translation_table_summary.tsv
```

A batch is the unit of restart and of sharing. Several machines may be given the
same `--out_dir` and will divide the release between them, each claiming batches
no other machine holds; a machine lost mid-batch costs that batch rather than the
run. Rerunning the command skips the batches that succeeded and retries those
that failed.

**A claim is a lease.** The machine holding a batch touches its `RUNNING` file
every five minutes, and a claim untouched for `--lease` hours (2 by default) is
taken by whichever machine next comes to the batch, whatever host made it. A
claim of this host's whose process has gone is taken at once. `--reclaim` takes a
claim before its lease is up, which is only ever right when the machine holding
it is known to have stopped. The lease is measured against the file server's
clock, so the machines sharing an `--out_dir` need not agree about the time.

**Hours are not redone.** `PREDICTED` is written the moment gTranslate returns,
so a batch taken over between the prediction and the comparison is compared
rather than predicted again. Within a batch gTranslate resumes by itself, keeping
each genome's called genes with a checksum beside them and skipping a genome
whose files verify, so a batch interrupted at genome 7,000 of 10,000 carries on
from there. Nothing removes a batch directory before retrying it, for that
reason.

**One bad genome does not cost a batch.** gTranslate ends a run when a worker
dies, and a genome it cannot process -- a few hundred bases with no genes to
count codons in, an assembly Prodigal refuses for its runs of N -- would take the
other 9,999 with it on every retry. It is therefore given `--force`, which drops
such a genome and carries on; the accessions dropped are written to
`no_prediction.tsv` in the batch directory and warned about, since `prodigal`
needs a table for every genome of the release. `--no_force` restores the older
behaviour of stopping the batch.

**Each machine writes its own log.** Several machines appending to one `--log`
over NFS overwrite one another and leave the file full of holes, so pass a
different `-l` per machine. What happens to a batch is written to
`trans_table.log` in the batch's own directory as well, which no other machine
writes to.

The batches are settled before any of them is processed, from the genomes sorted
by accession, so which genomes are in batch N follows from the set of genomes and
not from the order of the genome_dirs file. Once the batchfiles exist they are
authoritative: a later run reuses them and says so, since partitioning a release
again that has gained a genome would move genomes between batches that are
already finished. Remove the batch directories to partition it afresh.

`ncbi_tt_comparison.tsv` covers the genomes for which an NCBI table exists:

| Column | |
| --- | --- |
| `genome_id` | accession, as the genome_dirs file names it |
| `gtranslate_tt` | the table gTranslate predicted |
| `ncbi_tt` | the table NCBI declares in the genomic GFF |
| `checkm_tt` | the table the coding density rule alone would choose, as Prodigal and CheckM do unaided; it cannot express table 25 |
| `result` | `agree` or `conflict` |
| `coding_density_4`, `coding_density_11` | as gTranslate measured them |
| `ncbi_taxonomy` | lineage from `--taxonomy_file`, `na` where it holds none |

A genome NCBI has not annotated declares no table, so it has nothing to compare
against and is left out of the file rather than given a row saying so.

`prodigal` takes the summary `trans_table` writes as `--trans_table` and calls each
genome's genes under the table named there, so Prodigal no longer chooses one by
coding density. `--tt_override` corrects it: a TSV of `genome_id` and
`translation_table` whose rows replace the prediction. Every genome of the release
must have a table, from one file or the other, before any genes are called — a run
that found the gap genome by genome would find it hours in. Each genome's
`prodigal/prodigal_translation_table.tsv` records the table used and where it came
from, `predicted by gTranslate` or `specified by --tt_override`.

### Genome quality

| Command | Description |
| --- | --- |
| `checkm` | Run CheckM on new and modified genomes |
| `join_checkm` | Join CheckM output across GTDB versions |
| `prepare_checkm2` | Prepare files to run CheckM2 for the new release |
| `join_checkm2` | Join CheckM2 output files for different batches |
| `busco` | Estimate quality of new fungal genomes |

### Metadata

| Command | Description |
| --- | --- |
| `create_tables` | Create metadata tables for all NCBI genomes |
| `parse_assemblies` | Parse NCBI assembly summary files to generate metadata |
| `parse_ncbi_dir` | Parse the GTDB directory for extra NCBI metadata |
| `add_names_dmp` | Parse an NCBI `names.dmp` file into a table |
| `ncbi_genome_category` | Identify genomes marked by NCBI as a MAG or SAG |
| `generate_seqcode_table` | Generate a metadata table for genomes in SeqCode |

### Taxonomy

| Command | Description |
| --- | --- |
| `propagate_gtdb_taxonomy` | Propagate GTDB taxonomy to the new release |
| `propagate_curated_taxonomy` | Propagate curated taxonomy from representatives to their clusters |
| `update_propagated_tax` | Push propagated taxonomy to the new database |
| `add_taxonomy_to_database` | Update the taxonomy in the database |
| `set_gtdb_domain` | Set missing GTDB domain information from the NCBI domain |
| `curation_lists` | Lists and pseudo-trees for curation review |

### Database

| Command | Description |
| --- | --- |
| `update_db` | Update the GTDB PostgreSQL database |
| `update_checkm_db` | Import CheckM estimates |
| `update_metadata_db` | Update metadata in the database |
| `update_reps_db` | Update species cluster representatives |
| `update_ncbitax_db` | Update NCBI organism names and taxonomy |
| `update_taxid_to_db` | Add the NCBI taxid for each rank of each genome |
| `update_type_designation` | Update `type_designation` once SeqCode, NCBI and LPSN data are loaded |
| `add_surveillance_genomes` | Add surveillance genomes to a GTDB table |

### Nomenclature resources

| Command | Description |
| --- | --- |
| `lpsn` | LPSN processing (`pull_html`, `parse_html`, `lpsn_wf`, `add_metadata`) |
| `bacdive` | BacDive processing (`download_strains`) — in development |
| `strains` | Combine LPSN/DSMZ information (`date_table`, `type_table`) |
| `ncbi_strains` | Parse NCBI assembly reports and GenBank files for strain identifiers |

### Validation

| Command | Description |
| --- | --- |
| `overview` | Compare metadata files between releases |
| `compare_field` | Compare a specific metadata field between two metadata files |
| `compare_markers` | Compare marker frequencies between two releases |
| `compare_metadata_genome_dir` | Compare genomes listed in the metadata file against `genome_dirs` |
| `check_unique_strains` | Check for conflicting strain identifiers from the same collection |
| `check_db_population` | Check the database contains the expected number of genomes |

## Testing

Tests live in `tests/` and are written with `unittest`, run under `pytest`.

Install the test dependency and run the suite from the repository root:

```bash
pip install -e ".[test]"
pytest
```

`pytest` picks up `testpaths = ["tests"]` from `pyproject.toml`, so no arguments
are needed. To run a single test file, or a single test:

```bash
pytest tests/test_ncbi_genome_sync.py
pytest tests/test_ncbi_genome_sync.py -k manifest -v
```

The suite is plain `unittest`, so it can also be run without installing pytest:

```bash
python -m unittest tests.test_ncbi_genome_sync
```

The tests are offline — they use no network access and no mirror directory.

## Reference database versions

The releases of the marker and rRNA databases GTDB annotates against are set in
[gtdb_migration_tk/config.py](gtdb_migration_tk/config.py):

```python
PFAM_VERSION = '33.1'
TIGRFAM_VERSION = '15.0'
SILVA_VERSION = '138.2'
LTP_VERSION = '10_2024'
```

Each names the directory, inside a genome directory, that the results of that
database are written to. Two names derive from them:

```python
MARKER_FOLDER_SUFFIX = {'pfam': '33.1_lite', 'tigrfam': '15.0_lite'}
GTDB_DERIVED_DIRS_TO_COPY = ('prodigal', 'rna_silva_138.2', 'trna', 'rna_ltp_10_2024')
```

A genome whose sequences NCBI has not changed keeps the derived data of the
previous release. The genomic FASTA MD5 NCBI publishes decides that, but it is
the MD5 of the whole file, and NCBI reissues a FASTA with rewritten deflines --
a renamed organism, a relabelled assembly -- and identical sequences. Where the
two published MD5s disagree, `update_genomes` therefore hashes the sequences
themselves before throwing anything away: the same contigs, under the same IDs,
with the same bases, ignoring the free text after each ID, the line wrapping and
the base case. Those genomes are reported as `genomic FASTA sequences unchanged`
and keep their derived data. A genome whose contig was renamed does not, its
gene calls naming a contig the new FASTA no longer has.

`MARKER_FOLDER_SUFFIX` is the default `--folder_suffix` of `hmmsearch` and
`top_hit`, so by default they write `prodigal/pfam_33.1_lite/` and
`prodigal/tigrfam_15.0_lite/`; pass `--folder_suffix` only to annotate against
a version other than the declared one. `GTDB_DERIVED_DIRS_TO_COPY` lists the
derived data `update_genomes` carries across from the previous release when a
genome's genomic FASTA is unchanged; the Pfam and TIGRFAM results travel inside
`prodigal/`, so they are not listed separately.

`rna_silva` and `rna_ltp` still take their database version on the command line
(`--silva_version`, `--ltp_version`), and the value passed must match
`config.py`, or the directories carried across will not be the ones later steps
read. Updating a database is a matter of editing one value here and re-running
the commands that use it; nothing else in the code carries a version number.

## Repository layout

```
bin/gtdb_migration_tk      executable wrapper; defers to gtdb_migration_tk/__main__.py
gtdb_migration_tk/
    __main__.py            command-line interface: argument definitions
    main.py                OptionsParser: dispatches each command to its manager
    config.py              Pfam/TIGRFAM/SILVA/LTP versions; settings that change per release
    *_manager.py           implementation of each pipeline step
    biolib_lite/           vendored helpers (sequence I/O, taxonomy, parallelism)
    genometk_lite/         vendored genome metadata helpers
    gtdb_lite/             database import helpers
    utils/                 shared utilities
tests/                     test suite
```

## Version history

See [gtdb_migration_tk/VERSION](gtdb_migration_tk/VERSION) for the changelog.

## Copyright

Copyright © Pierre-Alain Chaumeil, Aaron Mussig and Donovan Parks. Released
under the GNU General Public License v3 (GPLv3), as stated in the header of each
source file.

Feature requests and bug reports can be sent to Donovan Parks
(donovan.parks@gmail.com) or posted on
[GitHub](https://github.com/Ecogenomics/gtdb-migration-tk/issues).
