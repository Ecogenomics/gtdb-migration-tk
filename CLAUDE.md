# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

An internal toolkit for moving the GTDB from one NCBI release to the next. Every
step is a subcommand of one executable, `gtdb_migration_tk`. It assumes access to
the GTDB genome directory layout and, for the `*_db` commands, the GTDB
PostgreSQL database. Python 3.8 or later.

## Commands

```bash
pip install -e ".[test]"                  # editable install with pytest
gtdb_migration_tk                         # list commands
gtdb_migration_tk <command> -h            # help for one command
python -m gtdb_migration_tk <command>     # equivalent
bin/gtdb_migration_tk <command>           # from a source checkout, no install needed
```

Tests are plain `unittest`, offline, and need no mirror or database:

```bash
pytest                                                # whole suite (testpaths in pyproject.toml)
pytest tests/test_ncbi_genome_sync.py -k manifest -v  # one file / one match
python -m unittest discover -s tests                  # without pytest
python -m unittest -v tests.test_select_genomes      # one module
```

There is no linter or formatter configured. Releases are cut by publishing a
GitHub release, which builds a conda package on the `ace-internal` channel from
`conda/meta.yaml`; that recipe reads the version from `setup.py`.

**Versioning.** `gtdb_migration_tk/VERSION` is the changelog. Its first line is
the version, read by both `setup.py` and `gtdb_migration_tk/__init__.py`. Bump it
by adding a new version line and its notes at the top.

## Architecture

### The CLI is three layers, and a new command touches all three

1. `__main__.py` defines every argument once as a `__name(group, required)`
   helper (e.g. `__database_setup`, `__cpus`, `__log_file`) and composes
   subparsers from them inside `subparser()` / `arg_group()` context managers.
   Arguments go in two groups, `'required named arguments'` then
   `'options arguments'`. The `print_help()` text is maintained by hand and is
   separate from argparse.
2. `main.py` `OptionsParser.parse_options()` is an if/elif chain on
   `options.subparser_name`, each branch calling a one-method-per-command wrapper
   that builds a manager and passes the `options` fields through.
3. A module holds the implementation as a class. The four NCBI commands live in
   modules named for them (`ncbi_metadata_sync.py` `NCBIMetadataSync`,
   `select_genomes.py` `SelectGenomes`, `ncbi_genome_sync.py` `NCBIGenomeSync`,
   `update_genomes.py` `UpdateGenomes`); every other command is a `*_manager.py`
   holding a `*Manager`.

So adding a command means: an argparse block and a `print_help()` line in
`__main__.py`, a method plus an elif in `main.py`, the implementation module, and
the command table in `README.md`.

### Logging and exit codes

`main()` in `__main__.py` calls `logger_setup()` (from `biolib_lite/logger.py`)
before dispatching. That creates two named loggers, `'timestamp'` and
`'no_timestamp'`; every manager does `logging.getLogger('timestamp')` rather
than creating its own. `--log` sets the log file; without it the log goes to
`./gtdb_migration_tk.log`. `parse_options()` returns an exit code, and `main()`
only calls `sys.exit()` when it is non-zero. Today only `ncbi_genome_sync` returns a
meaningful code (see README for the table); every other command returns 0.

### `ncbi_genome_sync.py` is deliberately self-contained

It is a standalone script grafted onto the toolkit. It owns its argparse via
`add_sync_arguments()`, which both `build_parser()` (standalone) and `__main__.py`
(as a subcommand) call, so the interface has one definition, `-l/--log` included.
Orchestration lives in `NCBIGenomeSync`, whose `run()` `main.py` calls like any
other manager; `main()` remains only as the script entry point, parsing argv and
opening the log when there is no toolkit to have done it. The workers it drives
(`sync_genome`, `verify_genome`, `http_get`, `prune_mirror`) stay module-level
functions, as do the rate limiter, circuit breaker and stop flag they share:
they run on `-j` threads at once and hold no per-run state. It has its own module logger and its own exit codes. The ~340 line
module docstring is the design document, with named sections (RATE LIMITING,
RESTART AND FRESHNESS, TUNING, SHARED OPERATION) that the inline comments refer
back to. Read it before changing sync behaviour. The argparse `dest` for the
summary file is `summary`, and the tests depend on that.

### NCBI assembly summary files are read by column name, never by position

`ncbi_utils.py` is the single reader, shared by `ncbi_genome_sync.py` and
`select_genomes.py`. It finds columns from the `#assembly_accession ...` header
row and refuses a table with no header (`BadInput`, a `ValueError`). NCBI has
grown `assembly_summary.txt` from 23 to 38 columns; a positional reader would
silently mirror the wrong files or build a release from the wrong genomes. Do
not slice these tables by index anywhere.

`ncbi_utils.py` also holds what more than one NCBI command knows about NCBI's
files: `REFSEQ_PREFIX`/`GENBANK_PREFIX`, `MD5_LINE_RE` (a line of
`md5checksums.txt`) and `has_ftp_path()` (an assembly NCBI lists but does not
serve). It is a leaf: it imports nothing from the package, and
`ncbi_genome_sync.py` imports from it and from nothing else in the package. The
command modules do not import one another; anything two of them need goes here.

### Release update: deciding vs. doing

`update_genomes.py` `UpdateGenomes` sorts the genomes of one database (given
as an accession prefix, `REFSEQ_PREFIX` or `GENBANK_PREFIX`) into removed, new
and shared by comparing the mirror's and the previous release's genome_dirs
files, filtered to that prefix; it reads no summary file, since the mirror is a
copy of the selection. `update_genomes` runs it once per prefix into one output
directory, with reports named for the prefix (`report_gcf.log`,
`gcf_to_review.log`, `report_gca.log`, `gca_to_review.log`). `FTPTools`, in the
same module, does the resulting copying, comparing and reporting.

Genome IDs are compared in canonical form via
`biolib_lite.common.canonical_gid()`: `GCF_005435135.1` and `GCA_005435135.1`
both become `G005435135`, which is how a GenBank genome is matched to its RefSeq
counterpart. Use it rather than slicing accessions.

The lingua franca between commands is the **genome_dirs file**: a TSV of
`accession<TAB>path`, one genome per line. `list_genomes` writes it
(`directory_manager.py`) by walking a tree and keeping the genomes named by
`--gtdb_selected_genomes`, and the update, comparison and validation commands
consume old, new and FTP variants of it. It says where each genome of a release
is held locally; the selection table says which genomes and where NCBI serves
them. Whether a tree holds what it should is `ncbi_genome_sync --verify`.

### `config.py` is the only place a reference database version lives

`PFAM_VERSION`, `TIGRFAM_VERSION`, `SILVA_VERSION` and `LTP_VERSION` each name
the directory inside a genome directory that database's results are written to.
Two names derive from them, and `tests/test_config.py` asserts the derivation
holds: `MARKER_FOLDER_SUFFIX` (`{'pfam': '33.1_lite', 'tigrfam': '15.0_lite'}`)
is the default `--folder_suffix` of `hmmsearch` and `top_hit`, resolved in
`main.py`, so those commands write `prodigal/pfam_33.1_lite/` unless told
otherwise; `GTDB_DERIVED_DIRS_TO_COPY` is the derived data `FTPTools` carries
across from the previous release when a genome's FASTA is unchanged. The
Pfam/TIGRFAM results and the version-free symlinks to them
(`prodigal/<gid>_pfam_lite.tsv.gz -> ./pfam_33.1_lite/...`) live inside
`prodigal/`, so copying `prodigal/` with `symlinks=True` carries them intact;
they are not listed separately. `marker_manager.py` itself does not read
`config.py`; it builds `pfam_<suffix>/` from whatever suffix it is handed.
`rna_silva` and `rna_ltp` take `--silva_version` and `--ltp_version` on the
command line, and those must match `config.py`.

### Database access

Credentials are always passed on the command line (`--hostname -u -d -p`, via the
`__database_setup` helper), never read from a config file.
`database_configuration/GenomeDatabaseConnectionFTPUpdate.py` is a thin psycopg2
wrapper most `*_db` managers use; `lpsn.py` and `ncbi_tax_manager.py` use
SQLAlchemy `create_engine` directly. `gtdb_lite/gtdb_importer.py` relies on an
`upsert` stored procedure that exists in the GTDB database.

### Vendored libraries

`biolib_lite/` and `gtdb_lite/` are trimmed copies of biolib and gtdb-lib,
vendored in 0.0.7 to drop those dependencies; `genometk_lite/` holds vendored
genome metadata helpers built on `biolib_lite`. Add helpers there rather than
reintroducing the upstream packages.

### Package naming gotcha

`gtdb_migration_tk/utils/` is a package (`common.py`, `tools.py`,
`prettytable.py`). A sibling `gtdb_migration_tk/utils.py` would be silently
shadowed by it and never importable. Put small shared helpers in
`utils/common.py`.

## Conventions

- Every source file starts with the GPLv3 header block. Docstrings use
  numpy-style `Parameters` sections and end with an `@return:` line.
- Module docstrings in the refactored modules (`ncbi_metadata_sync.py`,
  `select_genomes.py`, `update_genomes.py`, `ncbi_utils.py`, `config.py`,
  `ncbi_genome_sync.py`) explain *why* the code is shaped as it is, not what it
  does. Keep that up when touching them.
- Tests live in `tests/test_<module>.py`, one `TempDirCase` base for anything
  touching disk. Test names read as sentences about the contract that would
  otherwise break silently in production.
- Commit messages use `feat:`, `fix:`, `chore:` prefixes. Work happens on
  feature branches merged to `master` by pull request.
- Capitalisation in prose and messages: GTDB, RefSeq, GenBank, NCBI.
