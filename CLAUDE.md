# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

An internal toolkit for moving the GTDB from one NCBI release to the next. Every
step is a subcommand of one executable, `gtdb_migration_tk`. It assumes access to
the GTDB genome directory layout and, for the `*_db` commands, the GTDB
PostgreSQL database. Python 3.8 or later.

## Commands

Run the toolkit, its tests, and anything that imports it in the
`gtdb_migration_tk-r237` conda environment. That is the environment the r237
release is being built with: it holds every `install_requires` dependency, and
the package is installed in it against the working checkout, so what runs there
is the code in the tree rather than a copy of it. A change is verified by running
it there -- not under whatever `python` happens to be on `PATH`, which is an
interpreter too old to import the package.

```bash
conda activate gtdb_migration_tk-r237     # before anything below
gtdb_migration_tk                         # list commands
gtdb_migration_tk <command> -h            # help for one command
python -m gtdb_migration_tk <command>     # equivalent
bin/gtdb_migration_tk <command>           # from a source checkout, no install needed
```

Tests are plain `unittest`, offline, and need no mirror or database. `pytest` is
not installed in that environment, so run them with `unittest`, which the suite
is written in and which needs nothing extra:

```bash
python -m unittest discover -s tests                         # whole suite
python -m unittest tests.test_ncbi_utils                     # one module
python -m unittest -v tests.test_update_genomes.ResumeTests  # one class
```

`pip install -e ".[test]"` is how a NEW environment is set up, and brings pytest
with it (`pytest`, `testpaths` in `pyproject.toml`). Do not run it against
`gtdb_migration_tk-r237`: that environment is building a release, and it is not
the place to be resolving dependencies.

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

### `batching.py` is the coordination both long commands share

A release is a million-odd genomes and the work over it takes days, so `trans_table`
and `prodigal` both cut it into batches and let several machines divide them by
claiming batch directories under one `--out_dir`. That machinery -- the batchfiles,
the RUNNING/SUCCESS/FAILED canaries, claims as leases with a heartbeat, the
per-batch log, `plan_batches()`, and `write_table()`/`concatenate()` -- lives in
`batching.py` and knows nothing about gTranslate or Prodigal. What differs between
the commands in name alone -- what the batchfile is called, what the batch's log is
called -- travels in a `BatchLayout`, one per command. Add to `batching.py` rather
than to either command when both would want it.

What stays with a command is what a batch is FOR: `trans_table` owns `PREDICTED`
(the hours are over, only the comparison is left) and the files it hands gTranslate;
`prodigal` owns the decision that a genome's proteins are already vouched for.
`prodigal`'s `--out_dir` holds only the state of the run -- the called genes go into
the genome directories, which is why two machines on different batches never write
to the same place.

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
files: the database table (`REFSEQ`, `GENBANK`, `NCBI_DATABASES` and the
prefixes derived from them), `NCBI_HOST`/`NCBI_URL`, the naming of a saved
summary file (`assembly_summary_filename()`/`assembly_summary_database()`), the
columns every GTDB genome table opens with (`GENOME_COLUMNS`, `table_header()`),
`NCBI_NA`, `has_ftp_path()`, the manifest (`MD5_MANIFEST`, `GENOMIC_FASTA_EXT`,
`read_md5_manifest()`) and the block files are read and hashed in (`CHUNK`,
`file_md5()`). It is a leaf: it imports nothing from the package, and
`ncbi_genome_sync.py` imports from it and from nothing else in the package. The
command modules do not import one another; anything two of them need goes here.

### Release update: deciding vs. doing

`update_genomes.py` `UpdateGenomes` sorts the genomes of a release into removed,
new and shared by comparing the mirror's and the previous release's genome_dirs
files; it reads no summary file, since the mirror is a copy of the selection.
RefSeq and GenBank are done in ONE pass, every genome decided on its own
accession, writing `report.log`, `to_review.log` and `genome_dirs.tsv`. It ran once per accession
prefix until 0.1.7, from when `GenBankManager` needed the RefSeq run's output;
that decision now belongs to `select_genomes.py`. What the split reported for
free is kept as the per-database breakdown on every count logged
(`database_label()`, `tally_by_database()`, `count_by_database()`, and
`ComparisonTally.by_database` for the comparison outcomes). `FTPTools`, in the
same module, does the resulting copying, comparing and reporting.

Whether a shared genome keeps its derived data is decided by
`compare_genome_directories()` in two steps. The genomic FASTA MD5 published in
each `md5checksums.txt` is compared first, costing no read. Where those differ,
`sequences_md5()` hashes what the FASTA says the GENOME is -- the contig IDs and
the bases, ignoring the free text after each ID, the line wrapping and the base
case -- because NCBI reissues a FASTA with rewritten deflines and untouched
sequences, and the published MD5, being of the whole file, changes with them.
Equal sequences give `STATUS_SEQUENCES_UNCHANGED`, a fourth outcome that carries
the derived data across exactly as `STATUS_FASTA_UNCHANGED` does. The contig ID
and the division between contigs are deliberately part of the digest: the derived
data names the contigs it was called on, so a renamed or merged contig must
regenerate however unchanged the bases. The file is hashed as it decompresses,
nothing written.

The release tree splits the databases at the top where the mirror does not.
NCBI nests every genome of both under one `all/`
(`all/GCA/047/639/395/GCA_047639395.1_ASM4763939v1`); `release_genome_dir()`
keeps the nesting and replaces `all/` with the database's `name`, giving
`genbank/GCA/047/639/395/...` and `refseq/GCF/...`. The nesting is taken from the
mirror path, never rebuilt from the accession: NCBI defines it, the sync laid it
down from NCBI's URLs, and a reshaped release is no longer what
`ncbi_genome_sync --verify` checks.

Genome IDs are compared in canonical form via
`biolib_lite.common.canonical_gid()`: `GCF_005435135.1` and `GCA_005435135.1`
both become `G005435135`, which is how a GenBank genome is matched to its RefSeq
counterpart. Use it rather than slicing accessions.

The lingua franca between commands is the **genome_dirs file**: a headerless TSV
of `accession<TAB>absolute path<TAB>canonical accession`, one genome per line.
Readers split on tabs and ignore further columns, so columns may be appended but
never reordered. Two commands write one. `list_genomes`
(`directory_manager.py`) walks a tree and keeps the genomes named by
`--gtdb_selected_genomes`; that is how the MIRROR is indexed. `update_genomes`
writes `genome_dirs.tsv` for the release it builds, from the paths it placed
(`genome_dirs_row()`, `FTPTools.record_genome_dir()`), so the new release is not
walked back afterwards — only genomes whose directory was written are in it, and
a dry run, having written none, writes no file. The update, comparison and
validation commands consume old, new and FTP variants. A genome_dirs file says
where each genome of a release is held locally; the selection table says which
genomes and where NCBI serves them. Whether a tree holds what it should is
`ncbi_genome_sync --verify`.

### `config.py` is the only place a reference database version lives

`PFAM_VERSION`, `TIGRFAM_VERSION`, `SILVA_VERSION` and `LTP_VERSION` each name
the directory inside a genome directory that database's results are written to.
Two names derive from them, and `tests/test_config.py` asserts the derivation
holds: `MARKER_FOLDER_SUFFIX` (`{'pfam': '33.1_lite', 'tigrfam': '15.0_lite'}`)
is the default `--folder_suffix` of `hmmsearch` and `top_hit`, resolved in
`main.py`, so those commands write `prodigal/pfam_33.1_lite/` unless told
otherwise; `GTDB_DERIVED_DIRS_TO_COPY` is the derived data `FTPTools` carries
across from the previous release when a genome's sequences are unchanged. The
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
- Commit messages use `feat:`, `fix:`, `chore:` prefixes: a subject line, then
  a body explaining why the change is shaped as it is, as the module docstrings
  do. Work happens on feature branches merged to `master` by pull request, opened
  with the GitHub CLI:

  ```bash
  git checkout -b <branch>             # never commit to master
  git commit                           # feat:/fix:/chore: subject, then why
  git push -u origin <branch>
  gh pr create --base master --fill    # --fill takes title and body from the commit
  gh pr view --web                     # open it in a browser
  ```

  `gh` is GitHub's CLI and does what `git` does not: pull requests, issues,
  reviews, releases, Actions. It authenticates separately from `git push`, once,
  with `gh auth login`. Useful afterwards: `gh pr list`, `gh pr checks`,
  `gh pr view <n>`, `gh pr merge <n>`.
- Capitalisation in prose and messages: GTDB, RefSeq, GenBank, NCBI.
