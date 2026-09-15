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
  | `prodigal` | `prodigal`, `prodigal_check` |
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
# already present -- and the removal list in <base>.rm_dry_run. Changes nothing.
gtdb_migration_tk ncbi_genome_sync --gtdb_selected_genomes $S \
    --root /srv/db/gtdb/genomes -l ./logs/sync.log --dry-run

# remove what the selection does not list, then sync what it does
gtdb_migration_tk ncbi_genome_sync --gtdb_selected_genomes $S \
    --root /srv/db/gtdb/genomes -l ./logs/sync.log

# retry only the genomes that failed last time; a retry never removes anything
gtdb_migration_tk ncbi_genome_sync --retry ./logs/gtdb_selected_genomes.fail \
    --root /srv/db/gtdb/genomes -l ./logs/sync.log
```

`<base>.fail`, `<base>.bad` and `<base>.rm` (what was removed, each directory recorded
before it is deleted) are written to the `--log` directory. A failure file carries the
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
every listed genome present and md5-clean (failures to `<base>.bad`), and nothing else
present (directories the selection does not list to `<base>.extra`). Either fails the
verification; `--delete` removes both. Given `--retry`, only the listed genomes are
verified.

`ncbi_genome_sync` returns meaningful exit codes so it can be driven from a wrapper
script:

| Code | Meaning |
| --- | --- |
| `0` | everything synced / verified clean |
| `1` | some genomes failed — see `<base>.fail` / `<base>.bad` — or verification found directories the selection does not list — see `<base>.extra` |
| `2` | usage error, or a malformed assembly summary |
| `74` | filesystem refused the write (disk full, quota, read-only) |
| `75` | NCBI is throttling this host, or another sync holds `--root` — retry later |
| `130` / `143` | interrupted (SIGINT / SIGTERM); in-flight genomes recorded in `<base>.fail` |

## Commands

Run `gtdb_migration_tk <command> -h` for the arguments of any command.

### Mirror and update genome directories

| Command | Description |
| --- | --- |
| `ncbi_metadata_sync` | Download the NCBI taxonomy and the RefSeq and GenBank assembly summary files (archaea, bacteria and fungi) a release is built from, and generate the 7 rank NCBI taxonomy |
| `ncbi_genome_sync` | Sync NCBI data to a local directory |
| `select_genomes` | Select the NCBI genomes which will comprise the new GTDB release |
| `update_genomes` | Update RefSeq and GenBank genomes from the NCBI FTP mirror, carrying derived data across where the genomic FASTA is unchanged |
| `list_genomes` | Produce file indicating the directory of each genome |

### Gene calling and annotation

| Command | Description |
| --- | --- |
| `prodigal` | Call genes using Prodigal |
| `prodigal_check` | Check the Prodigal translation table matches NCBI |
| `hmmsearch` | Run HMMER on new and modified genomes |
| `top_hit` | Generate TopHit file for TIGRFAM or Pfam |
| `metadata` | Generate metadata derived from nucleotide and protein files |
| `rna_silva` | Identify and classify 16S, 23S and 5S rRNA genes against SILVA |
| `rna_ltp` | Identify and classify 16S rRNA genes against LTP |
| `update_silva` | Update taxonomy files and BLAST database from the latest SILVA release |
| `generate_ltp_db` | Generate BLAST database from the LTP website |
| `trnascan` | Identify tRNAs in genomes |

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
