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
python -m unittest -v tests.test_ncbi_ftp_manager     # one module
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
3. A `*_manager.py` module holds the implementation as a class.

So adding a command means: an argparse block and a `print_help()` line in
`__main__.py`, a method plus an elif in `main.py`, the manager, and the command
table in `README.md`.

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
(as a subcommand, with `--log` injected) call, so the interface has one
definition. It has its own module logger and its own exit codes. The ~340 line
module docstring is the design document, with named sections (RATE LIMITING,
RESTART AND FRESHNESS, TUNING, SHARED OPERATION) that the inline comments refer
back to. Read it before changing sync behaviour. The argparse `dest` for the
summary file is `summary`, and the tests depend on that.

### NCBI assembly summary files are read by column name, never by position

`ncbi_utils.py` is the single reader, shared by `ncbi_genome_sync.py` and
`ncbi_ftp_manager.py`. It finds columns from the `#assembly_accession ...` header
row and refuses a table with no header (`BadInput`, a `ValueError`). NCBI has
grown `assembly_summary.txt` from 23 to 38 columns; a positional reader would
silently mirror the wrong files or build a release from the wrong genomes. Do
not slice these tables by index anywhere.

### Release update: deciding vs. doing

`ncbi_ftp_manager.py` decides which genomes belong in a release (`RefSeqManager`
wants every "latest" RefSeq assembly; `GenBankManager` wants GenBank assemblies
only where RefSeq falls short, logging each decision to `gca_selection.log`).
`ncbi_ftp_manager_tools.py` `FTPTools` does the resulting copying, comparing and
reporting, and is the only consumer of `config.py`.

Genome IDs are compared in canonical form via
`biolib_lite.common.canonical_gid()`: `GCF_005435135.1` and `GCA_005435135.1`
both become `G005435135`, which is how a GenBank genome is matched to its RefSeq
counterpart. Use it rather than slicing accessions.

The lingua franca between commands is the **genome_dirs file**: a TSV of
`accession<TAB>path`, one genome per line. `list_genomes` writes it
(`directory_manager.py`), and the update, comparison and validation commands
consume old, new and FTP variants of it.

### `config.py` is the only place a marker database version lives

`PFAM_VERSION` and `TIGRFAM_VERSION` there derive every directory name, file
suffix and symlink used for marker annotations. `tests/test_config.py` asserts
the derivation holds. `marker_manager.py` does not read `config.py`: the
`hmmsearch` and `top_hit` commands take `--folder_suffix` (e.g. `33.1_lite`)
and build `pfam_<suffix>/` and `_pfam_<suffix>.tsv` from it. That suffix must
match what `config.py` derives (`pfam_33.1_lite`), otherwise `FTPTools` creates
symlinks to annotation files that were never written. The README and the
`config.py` docstring call this flag `--hmm_version`; that name is stale, the
flag is `--folder_suffix`.

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
- Module docstrings in the refactored modules (`ncbi_ftp_manager.py`,
  `ncbi_utils.py`, `config.py`, `ncbi_genome_sync.py`) explain *why* the code is
  shaped as it is, not what it does. Keep that up when touching them.
- Tests live in `tests/test_<module>.py`, one `TempDirCase` base for anything
  touching disk. Test names read as sentences about the contract that would
  otherwise break silently in production.
- Commit messages use `feat:`, `fix:`, `chore:` prefixes. Work happens on
  feature branches merged to `master` by pull request.
- Capitalisation in prose and messages: GTDB, RefSeq, GenBank, NCBI.
