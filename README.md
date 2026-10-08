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
  | `checkm` | `checkm` |
  | `checkm2` | `checkm2`, `trans_table` (quality of the genomes NCBI disagrees about) |
  | `hmmsearch` | `hmmsearch`, `top_hit` |
  | `blastn`, `blastp`, `makeblastdb` | `rna_silva`, `rna_ltp`, `generate_ltp_db` |
  | `nhmmer` | `rna_silva`, `rna_ltp` |
  | `tRNAscan-SE` | `trnascan` |
  | `busco` | `busco` |

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

Most commands take `-l/--log` and write a run log there. Where that file cannot
be opened — `-l logs/run.log` where `./logs` is a file, which is one keystroke
from `-l logs` — the run says so in one line and is logged to
`gtdb_migration_tk.log` under the command's `--out_dir` instead, that being where
the rest of what the run produces goes; a command with no `--out_dir` writes it
to the current directory. Add `--silent` to suppress console output.

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

`select_genomes` passes over a genome whose `genome_size_ungapped` in its
assembly summary file is below `--min_genome_size` (10 kbp by default) or above
`--max_genome_size` (100 Mbp). Genomes under the minimum are fragments or plasmids
filed as genomes. Genomes over the maximum are metagenomes deposited as one genome.
Over r237's summary files the defaults pass over 139 GenBank genomes under 10 kbp
and 8 over 100 Mbp, seven of them between 1.1 and 9.5 Gbp. Each is named in
`genomes_filtered_by_size.tsv` in `--output_dir` (`genome_id`,
`genome_size_ungapped`, `reason` of `genome_too_small` or `genome_too_large`),
which is written even when it has no rows. A RefSeq genome passed over for its size
leaves its GenBank counterpart to be judged on its own size. A genome whose summary
file states no size (archived files predate the column) is selected, and the run
says how many there were.

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

`report.log` is a headerless TSV of `accession<TAB>outcome`, one genome per line.
The outcome is `new`, `removed`, or one of the four outcomes of a comparison:
`genomic FASTA file unchanged`, `genomic FASTA sequences unchanged`,
`genomic FASTA file changed`, or `to_curate;<exception>: <message>`. The commands
that come after the update take their `--report` from it and read it through
`genomes_to_regenerate()` and `genomes_in_release()` in `update_genomes.py`
rather than parsing it themselves: the first gives the genomes the release did
not bring derived data with (`new` and `genomic FASTA file changed`), which is
what `hmmsearch` and `checkm` work on; the second gives every genome
the release holds a directory for, which is what `--all_genomes` asks for. A
`removed` genome is not in the release and a `to_curate` genome was never copied
into it, so neither is ever handed to a later command. A genome new to NCBI that
NCBI publishes without its genomic FASTA -- a directory of reports alone, as
GCA_056491145.1 and GCF_056491165.1 are -- is `to_curate;no genomic FASTA: ...`
rather than `new`, and is added by a later release once NCBI publishes it. A report written before
0.0.9 opens with a domain column and is refused rather than read as empty.

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
1, as it does for every command; give it 16. Neither half of the work is bound by
the CPU: copying a genome is round trips to the file server, and the comparison
only reads a FASTA where NCBI has reissued one. Measured against a mirror on NFS, 200
genomes of the size NCBI serves:

| `--cpus` | copying (`--fresh`) | comparing, MD5s agree | comparing, every FASTA reissued |
| --- | --- | --- | --- |
| 1 | 17.1s | 15.5s | 46.0s |
| 8 | 5.2s | 4.3s | 6.5s |
| 16 | 4.5s | 3.7s | 4.3s |
| 32 | 4.2s | 3.4s | 3.9s |
| 64 | 4.0s | -- | 4.2s |

Hence 16. Raise it to 32 for a release NCBI has reissued heavily;
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
| `align_marker_genes` | Align the marker genes of the database's genomes into `aligned_markers`, in batches under `--out_dir` |
| `genomic_metadata` | Generate metadata derived from nucleotide and protein files |
| `rna_silva` | Identify and classify 16S, 23S and 5S rRNA genes against SILVA, in batches under `--out_dir` |
| `rna_ltp` | Classify the 16S rRNA genes `rna_silva` extracted against LTP, in batches under `--out_dir` |
| `update_silva` | Update taxonomy files and BLAST database from the latest SILVA release |
| `generate_ltp_db` | Generate BLAST database from the LTP website |
| `trnascan` | Identify tRNAs in genomes, in batches under `--out_dir` |

`trans_table` runs gTranslate over the genomes of a release in batches of
`--batch_size` (default 10,000), each batch a directory of its own under
`--out_dir`:

```
<out_dir>/
  batch_000001/
    gtranslate_batchfile.tsv.gz         the genomes of this batch
    gtranslate_batchfile_present.tsv    the copy gTranslate reads, while it runs
    RUNNING                             a machine is working on it (host, PID, time)
    PREDICTED                           gTranslate has run; its results are final
    SUCCESS                             it finished; its results are complete
    FAILED                              it was attempted and something went wrong
    FAILED.20260919T140328              what an earlier attempt said, kept
    trans_table.log                     what this command did to this batch
    gtranslate.log                      what gTranslate did to it
    no_prediction.tsv                   genomes gTranslate returned nothing for
    gtranslate.translation_table_summary.tsv
    ncbi_tt_conflict.tsv                genomes of this batch NCBI disagrees about
    gtranslate_ncbi_tt_comparison.tsv.gz  every genome of this batch both called
  batch_000002/
  checkm2/                              working space, removed once the estimates
                                        are in ncbi_tt_conflict.tsv
  ncbi_tt_conflict.tsv                  the whole release, once every batch has SUCCESS
  gtranslate_ncbi_tt_comparison.tsv.gz
  gtranslate_no_prediction.tsv          genomes the release has no table for
  gtranslate.translation_table_summary.tsv.gz
```

The release summary is gzipped -- a row per genome, 116 MB of text for r237 and
34 MB compressed -- while `ncbi_tt_conflict.tsv` is a few hundred rows meant to
be read and is not. A batch's own summary is gTranslate's output and is left as
gTranslate wrote it. `prodigal` takes either: the summary is read by its first two
bytes rather than by its name, so a release summary that has been gunzipped, or
renamed on the way to another machine, still reads.

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

**The plan is kept compressed; gTranslate is handed a copy.** gTranslate opens a
batchfile with a plain `open()`, so `gtranslate_batchfile_present.tsv` is written
before it starts and removed once it has finished: a finished batch keeps
`gtranslate_batchfile.tsv.gz` alone. A batch taken over mid-run writes that copy
again from the plan, and a batch gTranslate failed on keeps it, being what a retry
looks at to see what went in. A batch planned by an earlier version holds an
uncompressed `gtranslate_batchfile.tsv`, which is read as it stands -- a release
whose batches looked unplanned would be partitioned again with its batches already
done.

**Hours are not redone.** `PREDICTED` is written the moment gTranslate returns,
so a batch taken over between the prediction and the comparison is compared
rather than predicted again. Within a batch gTranslate resumes by itself, keeping
each genome's called genes with a checksum beside them and skipping a genome
whose files verify, so a batch interrupted at genome 7,000 of 10,000 carries on
from there. Nothing removes a batch directory before retrying it, for that
reason.

**The genomes with no table are named once, for the release.**
`gtranslate_no_prediction.tsv` lists every genome the release has no translation
table for, with the reason:

| `reason` | |
| --- | --- |
| `gtranslate_returned_no_prediction` | gTranslate was handed the genome and returned nothing for it |
| `no_genomic_fasta` | the genome had no genomic FASTA, so it was never handed over |

The two are different failures and the row says which, because reporting a missing
file as a prediction failure sends the reader looking at the wrong thing. The file
is written even when there is nothing in it, so a release with nothing missing
says so. It is worked out per batch from what the batch was asked about against
what its summary answered, not concatenated from the batches' own
`no_prediction.tsv` -- those are written by the run that *predicts* a batch, and a
batch already predicted is never predicted again. For r237 it is eight genomes of
1,346,118: six of the first kind and two of the second. `prodigal` needs a table
for every genome of the release, so this is the list to correct with
`--tt_override`.

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

`ncbi_tt_conflict.tsv` holds the genomes where the table gTranslate predicted is
not the one NCBI declares, one row each:

| Column | |
| --- | --- |
| `genome_id` | accession, as the genome_dirs file names it |
| `gtranslate_tt` | the table gTranslate predicted |
| `ncbi_tt` | the table NCBI declares in the genomic GFF |
| `checkm_tt` | the table the coding density rule alone would choose, as Prodigal and CheckM do unaided; it cannot express table 25 |
| `checkm_conflict` | `True` where gTranslate and that rule disagree about the genome being recoded at all: 11 against 4, or 4 or 25 against 11. 25 against 4 is `False` -- the rule picks between 4 and 11 alone, so 4 is the closest it can come to saying 25 |
| `coding_density_4`, `coding_density_11` | as gTranslate measured them |
| `gc_percent`, `n50`, `genome_size` | as gTranslate measured them; what tells a conflict about a real genome from one about 200 kb of something barely assembled |
| `cm2_completeness_gtranslate_tt`, `cm2_contamination_gtranslate_tt` | what CheckM2 makes of the genome with its genes called under `gtranslate_tt` |
| `pass_qc_gtranslate_tt` | `True` where those pass standard GTDB QC, `na` where CheckM2 returned nothing |
| `cm2_completeness_ncbi_tt`, `cm2_contamination_ncbi_tt` | the same under `ncbi_tt` |
| `pass_qc_ncbi_tt` | and the same verdict under `ncbi_tt` |
| `ncbi_taxonomy` | lineage from `--taxonomy_file`, `na` where it holds none |

Every row is a conflict, so there is no column saying so -- `ncbi_conflict` lives
in the comparison, where it tells the rows apart. A genome NCBI has not annotated
declares no table and is in neither file.

`gtranslate_ncbi_tt_comparison.tsv.gz` holds **every** genome the two both called,
agreements and all -- 786,144 of r237's 1.35M, which is why it is gzipped in the
batch as well as in the release:

| Column | |
| --- | --- |
| `genome_id`, `gtranslate_tt`, `ncbi_tt`, `checkm_tt` | as in the conflict file |
| `ncbi_conflict` | `True` where gTranslate and NCBI named different tables; these are the rows the conflict file holds |
| `checkm_conflict` | the narrower question, as above: whether gTranslate and the density rule disagree about the genome being recoded at all |
| `coding_density_4`, `coding_density_11`, `gc_percent`, `n50`, `genome_size` | as gTranslate measured them |
| `ncbi_taxonomy` | lineage from `--taxonomy_file` |

The agreements are there because the rate the two differ at, and whether the
genomes they differ about are unlike the ones they agree about, are questions the
agreements have to be present to answer. The conflicts are that file filtered, not
a second walk over the genomes -- every GFF has already been read once -- so the
two cannot come to disagree about which genomes conflicted. How many genomes NCBI declares a table for, how many of those the two
agree and disagree about, and what percentage of them disagree is logged for each
batch and again for the release. The rate is of the genomes that could be
compared, not of the release: a genome NCBI has not annotated is not one the two
agree or disagree about.

**The two CheckM2 columns per table are the evidence about that table.** Genes
called under the wrong genetic code are truncated at every TGA and the markers
CheckM2 counts go with them, so each conflicting genome is put to CheckM2 twice,
once with Prodigal forced to `gtranslate_tt` and once forced to `ncbi_tt`. One
run would say how good the genome is; two say which table makes it look like a
genome at all. CheckM2's own choice is not asked for -- left to itself it picks
between tables 4 and 11 by coding density, which is what `checkm_tt` already
reports.

**Standard GTDB QC** is completeness > 50%, contamination < 10%, and a quality
score of `completeness - 5 * contamination` > 50. All three must hold: the score
alone would keep a genome 96% complete and 9% contaminated, and the completeness
alone would keep anything that had been sequenced. The thresholds are exclusive,
so a genome exactly 50% complete does not pass. The verdict is given per table
because a genome can pass under one and fail under the other -- which is the case
worth looking at, the conflict having changed whether the release keeps the
genome at all and not merely by how much. A genome CheckM2 returned nothing for
is `na` rather than `False`: it was not looked at and found wanting.

These runs are made once for the release, after every batch has succeeded, and
grouped by table: the conflicts are a few hundred genomes of a million-odd, and
CheckM2 searches the whole DIAMOND database once per run whatever the run holds.
A release already predicted therefore picks them up by running the command again
-- the batches are `SUCCESS` and are skipped, and the work happens where the
release files are written. Run that pass on ONE machine: the batches are claimed
one machine at a time but the release files are not, so several machines re-run
together would each start CheckM2 in the same directories, and `checkm2 predict
--force` empties its output directory as it starts.

The `checkm2/` directory is removed once the estimates are in
`ncbi_tt_conflict.tsv`. What it holds -- the staged genomes, the called proteins,
the DIAMOND output -- is about 600 KB per genome per table and says nothing the
conflict file does not now say. Running the command again therefore makes those
runs afresh rather than reading them, which is minutes for the few hundred
genomes a release conflicts about; the batches are what must never be redone, and
they are not. A table CheckM2 produced nothing for keeps the directory, so
retrying it does not also redo the tables that worked. A genome with no FASTA, or a run that fails, leaves
`na` in those four columns and the rest of the row intact; nothing here can cost
the release the comparison, which is written and counted before CheckM2 starts.

`prodigal` takes the summary `trans_table` writes as `--trans_table` and calls each
genome's genes under the table named there, so Prodigal no longer chooses one by
coding density. `--tt_override` corrects it: a TSV of `genome_id` and
`translation_table` whose rows replace the prediction.

`prodigal` takes an `--out_dir` and cuts the release into batches under it, the
same machinery `trans_table` uses (`batching.py`). Several machines may be given
the same `--out_dir` and will divide the release between them, each claiming
batches no other machine holds; `--batch_size`, `--reclaim` and `--lease` behave
as they do there, and a batch is skipped once it has `SUCCESS`.

**`--out_dir` holds the state of the run and nothing else.** The called genes go
into each genome's own `prodigal/` directory, as they always have, which is also
why two machines on different batches never write to the same place:

```
<out_dir>/
  batch_000001/
    prodigal_batchfile.tsv.gz         the genomes of this batch
    RUNNING / SUCCESS / FAILED        as in trans_table
    prodigal.log                      what this command did to this batch
    not_called.tsv                    genomes of this batch that got no genes
    meta_fallback.tsv                 genomes of this batch called in meta mode
  prodigal_not_called.tsv             the whole release, once every batch has SUCCESS
  prodigal_meta_fallback.tsv          the whole release, once every batch has SUCCESS
```

`prodigal_not_called.tsv` names every genome the release has no genes for, with
the reason: `no_translation_table`, `no_genomic_fasta`, or `prodigal_failed`. The
next command needs to know which genomes have no proteins, and it should not have
to look in 135 directories to find out.

**A genome single mode refuses is called in meta mode.** Single mode trains
Prodigal's model on the genome itself and is what a genome big enough to train on
is called under; it refuses a draft assembly it cannot train on — `saw too many
regions of N's` — however ordinary the bases between the gaps. Meta mode uses
Prodigal's precalculated parameters instead, under the SAME translation table.
Eleven genomes of one 2014 submission of N-rich actinomycetes, 9 to 13 Mb each,
had been called under neither mode since 2020; meta mode calls eight to twelve
thousand genes for each of them. gTranslate falls back the same way.

Those genomes are named in `prodigal_meta_fallback.tsv` under `--out_dir`, with
what single mode said about each, and each one's own
`prodigal/prodigal_translation_table.tsv` gains a `prodigal_mode` line saying the
same thing. Their proteins were called from precalculated parameters rather than
from a model trained on the genome, and afterwards the proteome looks like any
other, so unless the run says which genomes those were nothing downstream can
tell. A genome called as asked gets no such line.

**A genome Prodigal refuses costs that genome and not its batch.** Where neither
mode will call it, it is named `prodigal_failed` and what both modes said is in
the log; the other ten thousand genomes of the batch are called as usual. A batch
of several genomes in which EVERY one was refused is failed instead, that being
Prodigal not working on the machine rather than a batch of difficult genomes.

A rerun skips finished batches rather than re-reading the proteins of the whole
release, which is what batching buys over the per-genome checksum alone. The
checksum still decides genome by genome inside a batch that is not finished, which
is what a batch retried after a failure leans on. `--all_genomes` does the
finished batches again, since otherwise it would skip every batch it was asked to
redo.

**A genome is skipped only where its proteins are there, vouched for, and hold
something.** The digest is of the decompressed protein file, so a genome is judged
by its amino acids rather than by the gzip container. An empty protein file agrees
with its digest exactly — Prodigal leaves one where it failed, and
`da39a3ee5e6b4b0d3255bfef95601890afd80709` is the digest of nothing — so what is
asked is whether there are proteins in it. Otherwise a failure is carried from
release to release, re-vouched for at every step and named nowhere.

**A genome with no table is not called.** gTranslate returns no prediction for a
handful of genomes of a release — eight of r237's 1.35M, two of which have no
genomic FASTA to predict from at all — and there is nothing to call their genes
under; letting Prodigal pick a table by coding density is the very thing handing it
the summary prevents. They are named in the log and left, and `--tt_override` is
how one is given a table and called after all. `trans_table` lists them in full in
`gtranslate_no_prediction.tsv`.

Which genomes those are is settled before any genes are called rather than
discovered one at a time, because a run that meets the gap genome by genome meets
it hours in. A summary covering **no** genome of the release still stops the run:
that is the wrong file rather than a few unpredictable genomes, and carrying on
would call nothing and report that the run had finished.

Each genome's `prodigal/prodigal_translation_table.tsv` records the table used and
where it came from, `predicted by gTranslate` or `specified by --tt_override`.

`hmmsearch` takes an `--out_dir` and is batched the same way, so the marker search
can be spread over several machines as the gene calling is. `--batch_size`,
`--reclaim`, `--lease` and `--all` behave as they do for `prodigal`, and a batch
is skipped once it has `SUCCESS`.

**The batches of one database live under `<out_dir>/<marker directory>/.`** One
run searches one database at one version, and a run of the other database over the
same release is different work for the same genomes: sharing batch directories
would have the `SUCCESS` of a Pfam batch tell a TIGRFAM run that batch was done.
So one `--out_dir` carries both, and a new Pfam release gets its own directory
beside them for the same reason.

```
<out_dir>/
  pfam_33.1_lite/
    batch_000001/
      hmmsearch_batchfile.tsv.gz      the genomes of this batch
      RUNNING / SUCCESS / FAILED      as in trans_table
      hmmsearch.log                   what this command did to this batch
      not_searched.tsv                genomes of this batch that got no marker table
    hmmsearch_not_searched.tsv        the whole release, once every batch has SUCCESS
  tigrfam_15.0_lite/
    batch_000001/
    ...
```

The batchfile names each genome's **protein** FASTA, `prodigal/<gid>_protein.faa.gz`,
where `trans_table` and `prodigal` name the genomic FASTA: it is the file this
command reads. A genome whose proteins are not there, or are an empty file, is
named in `not_searched.tsv` under `no_protein_file` and left; it does not cost the
other ten thousand genomes of the batch. Those are the genomes `prodigal` listed
in `prodigal_not_called.tsv`, seen from the other side.

**What decides the work is the marker table, not the report.** A genome is skipped
where `prodigal/<marker dir>/<gid>_<marker>.tsv.gz` is there and its `.sha256`
agrees. `--report` is the release's own account of which genomes did not carry
their derived data across, and it is used to flag disagreement: a genome annotated
already though the report calls it new, or unannotated though the report does not,
is searched and said so in the log. Neither answer withholds the work. `--all`
discards what is there and searches every genome again, finished batches included.

Both `--db` runs write into `prodigal/` in each genome directory, so `top_hit` and
everything downstream read them exactly where they always have. `--out_dir` holds
the state of the run and nothing else, which is why two machines on different
batches never write to the same place.

**`--hmm_db_path` is a directory for `--db pfam` and a file for `--db tigrfam`.**
Pfam is searched by `PfamScan`, which is handed the directory `Pfam-A.hmm` sits
in; TIGRFAM is searched by `hmmsearch`, which is handed the HMM file itself.

```
--db pfam     --hmm_db_path /srv/db/gtdb/marker_genes/hmms_extended_pfam33.1_tigr15
--db tigrfam  --hmm_db_path /srv/db/gtdb/marker_genes/hmms_extended_pfam33.1_tigr15/tigrfam.hmm
```

The path is checked before the release is cut into batches, and a run given the
one where the other was wanted stops there with a line saying so. It used not to
be: `hmmsearch` reads a directory as a file that "appears to be empty", wrote no
marker table, and the run met that one call later as a `FileNotFoundError` on the
table in a worker.

`align_marker_genes` aligns the markers of the marker sets it is given
(`--marker_set_ids`, `marker_sets.id`: 1 = bac120, 2 = ar122, 19 = ar53) for the NCBI
genomes of the database and writes them to `aligned_markers`, which GTDB's MSAs and
trees are built from. It was `gtdb -r power realign_updated_genomes` of the gtdb
package until 0.1.66. `--new_genomes` aligns the genomes with no row for any of the
markers -- those `update_db` added, and those whose sequences changed, whose rows
it removes; `--all_genomes` aligns every NCBI genome again, rewriting its rows,
which is rarely wanted. One of the two is required. HMMER's `hmmalign` must be on
`PATH`, 3.2 or later: over 4,400 genomes added from 2015 to 2025, hmmalign 3.2.1, 3.3,
3.3.2 and 3.4 each rebuild all 520,859 rows of `aligned_markers` compared, and
3.1b1 and 3.1b2 differ in about 250 (TIGRFAM markers alone, 5% of genomes). The
`gtdb_migration_tk-r237` environment puts HMMER 3.4 first on `PATH` when activated.

```bash
gtdb_migration_tk align_marker_genes --db_service gtdb_r237 --marker_set_ids 1 2 --new_genomes \
    -g genome_dirs.tsv -o align_marker_genes -l align_marker_genes/align.log -c 16
```

For each marker the gene of the genome's top-hit table
(`prodigal/<marker dir>/<gid>_<marker dir>_tophit.tsv.gz`) naming it with the
highest bitscore, the first of equal ones, is aligned to the marker's HMM, and the
HMM's match states kept; `multiple_hits`, `hit_number` and `unique_genes` say how
many genes named it and how many distinct sequences they were. A genome is not
expected to have every marker: one no gene names is a row of gaps as long as the
HMM. The genes of a batch are aligned to a marker by one `hmmalign` rather than one
a gene, the alignment of each being the same; over 150 r237 genomes kept from r232,
every field of all 34,950 rows is the one the gtdb package wrote, but for one
e-value whose top-hit table was made again since. On 16 processes a batch of 1,000
genomes is read and aligned in about 26 seconds.

The markers are checked first: a set the database does not hold, or an HMM that is
missing, of another length than `markers.size`, or of another Pfam or TIGRFAM
version than `config.py`, refuses the run. So does a genome to align that the
genome_dirs file does not locate, every one listed in `not_in_genome_dirs.tsv`. A
genome without called proteins (`prodigal/<gid>_protein.faa.gz`) is passed over: the
log counts them and `missing_protein_file.tsv` in `--out_dir` lists them. One with
proteins and no top-hit table, or whose table names a gene its proteins do not hold,
fails its batch.

The genomes are cut into batches of `--batch_size` (1,000) under
`<out_dir>/marker_sets_<ids>_<new|all>/`, claimed as the other batched commands
claim theirs, so a stopped run, or several machines, carry on from the batches
finished. Each batch's rows are written in one transaction as it ends, and its
SUCCESS after that commits: a genome has all its rows or none, and a batch that
failed is rolled back, said in the log and in its FAILED file, and retried by the
next run, which exits 1.

`genomic_metadata` derives each genome's nucleotide statistics (GC, genome size, N50) and
gene statistics (protein count, coding bases, coding density) and writes them into
the genome's own directory as `metadata.genome_nt.tsv` and
`metadata.genome_gene.tsv`, with a `.desc.tsv` beside each naming the fields.
`create_tables` is what gathers them afterwards (see [Metadata](#metadata)).
`genomic_metadata` itself reads two files per genome:

| File | Written by |
| --- | --- |
| `<assembly>_genomic.fna.gz` | the mirror, carried across by `update_genomes` |
| `prodigal/<gid>_protein.gff.gz` | `prodigal` |

So `genomic_metadata` waits on `prodigal` and on nothing else. It reads no marker table
and no rRNA result, and can run while `hmmsearch` is still going -- the two write
to different places inside `prodigal/`.

**A genome missing a file is reported, not fatal.** The metadata of a release is
generated while the gene calling of its last genomes is still finishing, and until
0.1.28 a genome whose files were not there ended the command: the check called
`sys.exit()` from inside a pool worker, on the first such genome, throwing away
however many of the million were already done. Such a genome is now warned about
and named in `metadata_missing_files.tsv` in `--out_dir`, one row per missing file:

| `missing` | |
| --- | --- |
| `genomic_fasta` | the genome's sequences are not in the mirror |
| `protein_gff` | `prodigal` has not called this genome's genes |

with the path that was looked for in the third column. The file is written even
when there is nothing in it, so a release with nothing missing says so. `--out_dir`
is required and holds that report alone; the metadata itself goes into the genome
directories as it always has.

How much of a genome is done follows from which files it has. The nucleotide
statistics need only the FASTA and are written whenever the FASTA is there, so a
genome whose genes are not called yet still contributes its `metadata.genome_nt.tsv`
-- `create_tables` reads the two files independently, and that half is not
calculated again once `prodigal` catches up. The gene statistics are a coding
density, which needs the genome size as well as the GFF, so a genome with no
sequences yields neither file and is left untouched.

`--cpus` is a count of genomes in flight, not of threads: each is a worker process
handling whole genomes one after another, and nothing within a genome is parallel.

`trnascan` identifies the tRNAs of each genome with tRNAscan-SE, writing them into
the genome's own `trna/` directory as `<gid>_trna.tsv`, `<gid>_trna_stats.tsv` (the
file `create_tables` reads) and `<gid>_trna.log`, with a `.sha256` beside the
table. It is batched under `--out_dir` exactly as `prodigal` and `hmmsearch` are,
with the same `--batch_size`, `--reclaim`, `--lease` and `--all`:

```
<out_dir>/
  batch_000001/
    trnascan_batchfile.tsv.gz       the genomes of this batch
    RUNNING / SUCCESS / FAILED      as in trans_table
    trnascan.log                    what this command did to this batch
    not_scanned.tsv                 genomes of this batch that got no tRNAs
  trnascan_not_scanned.tsv          the whole release, once every batch has SUCCESS
```

**What decides the work is the checksum, not a report.** A genome is skipped where
`trna/<gid>_trna.tsv` is there and its `.sha256` agrees, which is why the command
takes no `--report`: `trna` is in `GTDB_DERIVED_DIRS_TO_COPY`, so a genome whose
sequences did not change carries its tRNAs and its checksum across from the
previous release and is skipped without anything having to look up what became of
it. A table whose checksum disagrees was written by a run interrupted partway
through it, and is scanned again.

**The domain decides the model.** tRNAscan-SE searches with a bacterial or an
archaeal model and the two give different answers, so each genome's domain comes
from the same two files `rna_silva` reads it from:

| Argument | |
| --- | --- |
| `-d`, `--gtdb_domain_file` | GTDB's own `Predicted domain`, from the genome's marker genes |
| `-t`, `--taxonomy_file` | the standardised NCBI taxonomy, for the genomes GTDB has no prediction for |

GTDB's call is preferred because it is made from the genome rather than from where
NCBI filed it, and it is the one that catches a genome under the wrong domain at
NCBI; `None` in that column means the markers gave no answer, and the NCBI lineage
answers for those. The taxonomy is matched on the accession and then on its
canonical form, so a GenBank genome finds the lineage recorded against its RefSeq
counterpart. A genome neither file answers for is scanned as a bacterium, and the
run says how many of those there were.

Until 0.1.29 the domain came from which of four NCBI assembly summary files an
accession appeared in. Nothing in those files says "bacteria" -- NCBI tells them
apart by directory -- so the answer was asserted by which argument each file was
passed as, and swapping two of them on the command line would have scanned every
archaeon as a bacterium in silence.

**A genome that got no tRNAs is named once, for the release**, in
`trnascan_not_scanned.tsv`:

| `reason` | |
| --- | --- |
| `no_genomic_fasta` | the genome's sequences are not where the release says they are |
| `trnascan_failed` | tRNAscan-SE was given the genome and returned an error |
| `genome_too_large` | the assembly is larger than `--max_genome_size` |

Neither costs the other ten thousand genomes of the batch, and neither is retried
for ever by a batch that fails identically every time. A batch in which several
genomes were to be scanned and every one of them failed is failed rather than
recorded as a success: that is not a batch of difficult genomes, it is tRNAscan-SE
not working on this machine. One genome failing on its own is not, for the reason
the naming exists.

**A genome larger than `--max_genome_size` is left out** (100 Mbp by default, well
above the largest bacterial or archaeal genome). What comes near it is a
metagenome deposited as one genome: GCA_964261755.1 is 9,529 Mbp in 107,915
contigs, and held an r237 `rna_silva` batch for a day. Its size is read from NCBI's
`<assembly>_assembly_stats.txt` (the `all all all all total-length` row), and
counted from the genomic FASTA only where that file is missing. A genome whose size
cannot be told is processed. Only genomes that are still to be processed are sized,
so a genome whose results are already there is not named. Each one left out is
logged with its size and named as `genome_too_large`. It gets no results and no
canary, so a later run with a larger `--max_genome_size` processes it.
`rna_silva` and `rna_ltp` take the same option.

`rna_silva` identifies, extracts and classifies one rRNA gene per run -- `-r ssu`,
`lsu_23S` or `lsu_5S` -- writing into each genome's own
`rna_silva_<version>/` directory (`ssu.fna`,
`ssu.hmm_summary.tsv`, `ssu.taxonomy.tsv`, ...). It is batched under `--out_dir`
as `trnascan` is, with the same `--batch_size`, `--reclaim`, `--lease`, `--all`
and `--tmp_dir`. Each gene at each SILVA version has batches of its own, for the
reason `hmmsearch` keeps a directory per database: an ssu batch and an lsu_23S
batch cover the same genomes and are different work.

```
<out_dir>/
  rna_silva_138.2/
    ssu/
      batch_000001/
        rna_silva_batchfile.tsv.gz  the genomes of this batch
        RUNNING / SUCCESS / FAILED  as in trans_table
        rna_silva.log               what this command did to this batch
        not_searched.tsv            genomes of this batch that were not searched
      rna_silva_not_searched.tsv    the whole release, once every batch has SUCCESS
    lsu_23S/
      ...
```

**What decides the work is the canary.** A genome is skipped where
`<gene>.canary.txt` is in its results directory, which is written once the gene has
been searched for, whether or not one was found. The results are made in
`--tmp_dir` and copied into place with the canary last, so a run stopped mid-copy
leaves a genome that is searched again. What an earlier search of the same gene
left is removed before the new results are copied in; `--remove` goes further and
empties the directory of every genome it searches, the other genes' results
included.

**The domain decides the HMM**, `bac_16S` against `ar_16S` and so on, and is read
exactly as `trnascan` reads it: `-d/--gtdb_domain_file` for GTDB's own prediction,
`-t/--taxonomy_file` for the standardised NCBI taxonomy where GTDB has none. Until
0.1.35 the fallback was the `NCBI taxonomy` column of the domain file itself. A
genome neither file answers for is searched as a bacterium, and each batch says how
many of those it had.

**A genome that could not be searched is named** in `rna_silva_not_searched.tsv`:

| `reason` | |
| --- | --- |
| `no_genomic_fasta` | the genome's sequences are not where the release says they are |
| `rna_search_failed` | nhmmer or blastn was given the genome and failed on it |
| `genome_too_large` | the assembly is larger than `--max_genome_size`, as for `trnascan` |

A genome with no copy of the gene is not among them: it was searched, and has a
canary. As with `trnascan`, a batch in which several genomes were to be searched and
every one failed is failed rather than recorded as a success.

`rna_ltp` classifies the 16S rRNA genes `rna_silva` extracted against the LTP, and
writes into each genome's own `rna_ltp_<ltp version>/` directory (`ssu.blastn.tsv`,
`ssu.taxonomy.tsv`, `ltp.canary.txt`). It is batched as `rna_silva` is, with the
same options, under `<out_dir>/rna_ltp_<ltp version>-silva_<ssu version>/`: both
versions decide what a genome's results are.

| Argument | |
| --- | --- |
| `--ltp_version` | the LTP release classified against |
| `-v`, `--ssu_version` | the SILVA version naming the `rna_silva_<version>/` directory the genes are read from |
| `-p`, `--rnapath` | the directory holding one directory per LTP release; its default, `/srv/db/silva/`, is `rna_silva`'s, so pass `-p /srv/db/silva/ltp` |

A batch is planned around each genome's `rna_silva_<version>/ssu.fna`, the file
the command reads. A genome with none is one of two things, and `rna_silva`'s own
`ssu.canary.txt` tells them apart: where it is there, `rna_silva` searched the
genome and found no 16S gene, which is counted (`no_ssu_gene` in each batch's
SUCCESS) and not named; where it is not, `rna_silva` has not searched the genome
yet, and it is named in `rna_ltp_not_classified.tsv`:

| `reason` | |
| --- | --- |
| `ssu_not_identified` | `rna_silva` has not searched the genome for its 16S gene |
| `blastn_failed` | blastn was given the genome's genes and failed on them |
| `genome_too_large` | the assembly is larger than `--max_genome_size`, as for `trnascan` |

`rna_silva` writes no canary for a genome it left out for its size, so a genome
with neither `ssu.fna` nor a canary is sized too. One that is too large is named
`genome_too_large`, not `ssu_not_identified`, since `rna_silva` will not be
searching it.

A genome is skipped where `ltp.canary.txt` is there, and its results are copied
into place with the canary last. The command takes no domain file: the domain
decides which HMM a gene is searched for with, which is `rna_silva`'s work, and
the LTP classification of a gene already extracted does not depend on it.

### Genome quality

| Command | Description |
| --- | --- |
| `checkm` | Run CheckM on new and modified genomes |
| `checkm2` | Run CheckM2 on new and modified genomes |
| `join_checkm` | Join CheckM output across GTDB versions |
| `busco` | Estimate quality of new fungal genomes |

`checkm` and `checkm2` estimate the completeness and contamination of the
genomes `--report` says are to be regenerated -- `new` and `genomic FASTA file
changed`, as `genomes_to_regenerate()` reads it -- and of no others; a genome
whose derived data the release carried across is not assessed again. `--all`
takes every genome the release holds instead, and a `--report` of `none` every
genome of the genome_dirs file. They cut those genomes into batches under `--out_dir` with the machinery
`trans_table`, `prodigal` and `trnascan` share, so several machines given one
`--out_dir` divide the work (`--batch_size`, default 1,000 for `checkm` and
5,000 for `checkm2`, `--reclaim`, `--lease`). The batches are planned in
`--out_dir` itself, so each command is given an `--out_dir` of its own; one that
holds batch directories another command planned is refused, since that
command's `SUCCESS` would tell this one it had nothing to do.

```bash
gtdb_migration_tk checkm2 -g genome_dirs.tsv --report report.log \
    -o checkm2 -l checkm2/checkm2.$(hostname).log --tmp_dir /tmp -c 32
```

Unlike the gene calling commands, their results live in the batches and not in
the genome directories, so the genomes to assess are decided once, when the
batches are planned, and only those are planned. A plan already under
`--out_dir` is used as it stands, whatever `--report` or `--all` a later run is
given; remove the batch directories to plan again.

Both are handed the proteins the `prodigal` command called (`--genes`) and take
them as correct, so the genes of a release are called once, under the table
`trans_table` chose, and CheckM and CheckM2 assess the proteome everything else
is made from. Left to call genes themselves, both would call them under tables 4
and 11 and keep whichever codes more of the genome, which cannot express table
25.

Once every batch has succeeded, whichever machine finishes last writes the
release files beside the batches, as `trans_table` does: each table is the
batches' own concatenated in batch order under a single header, the genomes left
out are gathered into one file written even when it has no rows, and the count
of genomes assessed is added up from every batch's `SUCCESS`. A batch in which
no genome could be assessed adds no rows. The tables, a row per genome, are
gzipped; the genomes left out and the versions are not. A table an earlier run
left uncompressed, and the `checkm.alignment_file.tsv` it joined from the
batches, are removed once the gzipped tables are written.

| File | |
| --- | --- |
| `checkm.profiles.tsv.gz` | CheckM `qa` joined with `tree_qa -o 2`, one row per genome, the bin named `<accession>_protein` |
| `checkm.qa_sh100.tsv.gz` | CheckM `qa --aai_strain 0.9999` |
| `checkm2.quality_report.tsv.gz` | CheckM2's `quality_report.tsv`, one row per genome, named by accession. Handed genes, CheckM2 writes none of the statistics it takes from calling genes itself (the table used, coding density, genome size, GC, N50 and the rest); its estimates are made from the proteins either way |
| `checkm_not_assessed.tsv`, `checkm2_not_assessed.tsv` | `genome_id`, `reason`, `detail` of each genome left out |
| `checkm.version`, `checkm2.version` | the version of CheckM or CheckM2 that made the release files, gathered from the `checkm.version` or `checkm2.version` each batch was given. Where the batches were made by more than one version, each is on a line of its own and the log names the batches each made |
| `batch_*/alignment_file.tsv` | `checkm` only: the alignments of multi-copy genes `qa --aai_strain 0.9999 -a` writes, kept in each batch and not gathered for the release. Nothing reads them, and joined for r237 they were 142 GB, 137 GB of it from one batch, copied again by every run that wrote the release files |
| `pplacer.version` | `checkm` only: the version of the pplacer CheckM placed the genomes with, gathered the same way. pplacer cannot be asked (bioconda's says `dev` to `--version`), so each batch learns it from the pplacer process CheckM starts, and the conda package that executable came from; the batch's log names the executable. It is written only where every batch that ran CheckM recorded one |

| `reason` | |
| --- | --- |
| `no_protein_file` | the genome has no `prodigal/<accession>_protein.faa.gz`, or it is empty |
| `genome_too_large` | the assembly is larger than `--max_genome_size`, as for `trnascan`; `detail` is its size in bases |

A genome larger than `--max_genome_size` (100 Mbp by default, as for
`select_genomes`) is sized when its batch is assessed, not when the batches are
planned, so the limit applies to a plan made before it. A batch that finishes
without such a genome is done: to assess it later under a larger limit, remove
that batch's `SUCCESS`.

`checkm` gives pplacer at most 64 threads, whatever `--cpus` says, since pplacer
runs badly on more than 64; the rest of CheckM is given all of them, so `-c 96`
runs CheckM on 96 and pplacer on 64.

CheckM or CheckM2 exiting non-zero fails the batch, which the next run repeats.
Both programs put a Unix socket under `--tmp_dir`, and a socket path is limited
to 107 characters, so give them a short one: a `--tmp_dir` longer than about 75
characters fails every batch with `OSError: AF_UNIX path too long`.

### Metadata

| Command | Description |
| --- | --- |
| `create_tables` | Gather the metadata written into each genome directory into the tables `update_metadata_db` loads |
| `parse_ncbi_assemblies` | Parse NCBI assembly summary files to generate metadata |
| `parse_ncbi_dir` | Parse the GTDB directory for extra NCBI metadata |
| `add_names_dmp` | Parse an NCBI `names.dmp` file into a table |
| `parse_ncbi_genome_category` | Identify genomes marked by NCBI as a MAG, SAG or environmental genome |
| `generate_seqcode_table` | Generate a metadata table for genomes in SeqCode |

`create_tables` calculates nothing. It walks the genomes of the genome_dirs file
it is given (`--gtdb_genome_path_file`), which may be any release's and not only
NCBI's, and gathers what earlier commands wrote into each genome directory into
nine tables in `--out_dir`, each gzipped (`metadata_nt.tsv.gz` and so on):

| Read from the genome directory | Written by | Table |
| --- | --- | --- |
| `metadata.genome_nt.tsv` | `genomic_metadata` | `metadata_nt.tsv` |
| `metadata.genome_gene.tsv` | `genomic_metadata` | `metadata_gene.tsv` |
| `rna_silva_<ver>/ssu.*` | `rna_silva -r ssu` | `metadata_ssu_silva.tsv`, `metadata_ssu_silva_count.tsv` |
| `rna_silva_<ver>/lsu_23S.*` | `rna_silva -r lsu_23S` | `metadata_lsu_silva_23s.tsv`, `metadata_lsu_silva_23s_count.tsv` |
| `rna_silva_<ver>/lsu_5S.*` | `rna_silva -r lsu_5S` | `metadata_lsu_5S.tsv`, `metadata_lsu_5S_count.tsv` |
| `trna/<gid>_trna_stats.tsv` | `trnascan` | `metadata_trna_count.tsv` |

`--silva_version` only names the `rna_silva_<ver>` directory, and must match
`config.SILVA_VERSION`. `metadata_ssu_gg.tsv`, the Greengenes classification of
the 16S rRNA genes from `ssu_gg/`, is no longer written: nothing has written
`ssu_gg/` since `rna_silva` took over, and r237 had none. A run removes one an
earlier run left in `--out_dir`, which `update_metadata_db --input_folder` would
otherwise refuse as a table it does not know. The eight `ssu_gg_*` fields were
dropped from `metadata_rna` and `metadata_view` of `gtdb_r237_dev` on
2026-10-08, the view recreated without them; the values of the 86,077 genomes
that held them (nearly all added from 2015 to 2018), and the view's definition
before and after, are in `release237/processing/drop_ssu_gg/`. A genome without a file is given no row in that table;
the three `*_count.tsv` tables are the exception, giving every genome a row, 0
where nothing was found. The run ends by logging, for every other table, how many
genomes have a row, how many had no file to read it from, and how many had one
with nothing to report (an rRNA table of no hits). A genome with no
`metadata.genome_nt.tsv` is a WARNING, since every genome has a genomic FASTA to
calculate it from; the rest are INFO, being rightly missing for some genomes:

```
Rows written for 2 genomes:
  metadata_nt.tsv: 2 with a row.
  metadata_gene.tsv: 1 with a row; 1 had no metadata.genome_gene.tsv.
  ...
  metadata_trna_count.tsv: 1 with a row; 1 had no trna/<gid>_trna_stats.tsv.
```

That is where a release finds that `genomic_metadata`, `rna_silva` or `trnascan`
did not get to every genome. The log goes to `-l/--log`, which is required, so
that the run's account of what is missing is where it was asked to be. A
genome_dirs file naming no genomes is
refused, exiting 1, before anything is written, rather than replacing the tables
in `--out_dir` with nine of no rows. The table names are the ones
`update_metadata_db --input_folder` knows, and it refuses any other `.tsv` there.

`create_tables`, `parse_ncbi_assemblies`, `parse_ncbi_dir`,
`parse_ncbi_genome_category` and `ncbi_strains` write their tables gzipped, `<table>.tsv.gz`, with no time or file name in the
gzip header, so the same genomes give the same bytes. A table an earlier run
left uncompressed in the same `--out_dir` is removed once the gzipped one is
written, since `update_metadata_db` refuses a folder holding a table both ways.
r237's `ncbi_assembly_summary.tsv` is 73 MB gzipped, where it was 598 MB, for 11
seconds more of a 38-second run.

`create_tables` reads `--cpus` genomes at once on threads (default 1; give it
8) and writes them in the order of the genome_dirs file, so the tables are the
same whatever `--cpus` is. The work is opening a dozen small files in each genome directory,
NFS round trips with nothing computed between them. Measured on r237's genome
directories, nothing cached:

| `--cpus` | per genome | 1.35M genomes |
| --- | --- | --- |
| 1 | 75 ms | ~28 h |
| 4 | 22 ms | ~8 h |
| 8 | 9-10 ms | ~3.5 h |
| 16, 32, 64 | 9-11 ms | no faster |

Past 8 the file server is what is waited on, so more threads only add load to a
server others are using.

`parse_ncbi_assemblies` (`parse_assemblies` until 0.1.55) takes the NCBI
assembly summaries the release was selected from (`-n`), as `select_genomes`
and `strains type_table` do, rather than `--rb`, `--ra`, `--gb` and `--ga`. It
writes `ncbi_assembly_summary.tsv.gz` in `-o/--out_dir`, one row for every genome of
the summaries, and requires a log (`-l`). It takes no list of genomes: which of
them the database loads is for `update_metadata_db` to decide.

```bash
gtdb_migration_tk parse_ncbi_assemblies \
    -n ncbi/assembly_summary_{archaea,bacteria}_{refseq,genbank}.txt.gz \
    -o parse_ncbi_assemblies -l parse_ncbi_assemblies/parse_ncbi_assemblies.log
```

The table is named as `update_metadata_db --input_folder` loads it, against
`metadata_ncbi_assembly_file.desc.tsv`, whose fields are the ones written.
Columns are found by name in each summary, and a summary may be gzipped. A
genome is written once, from the first summary holding it. The log says where
the table is written, how many of each summary's genomes are written, and the
total, and warns of a genome in more than one summary. Over r237's four
summaries that is 3,874,314 genomes, a 598 MB table, in under a minute: every
genome NCBI holds, not only the release's. `update_metadata_db` loads only the
genomes the database holds (every one of r237's 1,346,116 is in the
summaries), skipping and counting the rest.

`parse_ncbi_dir` reads what NCBI's own files in each genome directory of a
genome_dirs file (`-g`) say -- `_assembly_stats.txt`, `_genomic.gff.gz` and
`_genomic.gbff.gz` -- on `--cpus` processes, and writes
`ncbi_assembly_metadata.tsv.gz` in `-o/--out_dir`, the name `update_metadata_db
--input_folder` loads against `metadata_ncbi_assembly.desc.tsv`. It requires a
log (`-l`).

```bash
gtdb_migration_tk parse_ncbi_dir -g genome_dirs.tsv \
    -o parse_ncbi_dir -l parse_ncbi_dir/parse_ncbi_dir.log --cpus 32
```

Its time is the file server's: an uncached genome directory is several round
trips before a byte is read. The GenBank file is read only as far as the end of
its first record's source feature, which holds the source qualifiers
(`isolation_source`, `lat_lon` and the rest, the BioSample's, the same on every
record), and the translation table is taken from the GFF, read once for it and
its counts; the GenBank files of r237 are 1.9 MB gzipped on average, 2.5 TB in
all, for those four qualifiers. Measured on fresh r237 genomes:

| `--cpus` | before 0.1.60 | since |
| --- | --- | --- |
| 16 | 54 genomes/s | 66 genomes/s |
| 32 | 47 genomes/s | 78 genomes/s (~4.8 h for 1.35M) |
| 64 | 58 genomes/s | 76 genomes/s |

Past 32 it is no faster: opening one small file in each of 1.35M uncached genome
directories alone runs at ~160 a second from one machine.

Every genome is given a row, written as soon as it is read, in no particular
order. A genome missing one of the three files is empty in the fields read from
it, which `update_metadata_db` loads as NULL, and the log ends with a line for
each file, `Identified <n> genomes with a missing <file> file, e.g.: ...`: a
WARNING for the statistics and the GenBank file, which every assembly has, and
INFO for the GFF, which NCBI publishes only for the assemblies it annotated (121
of 300 r237 genomes drawn at random had none). Four of its columns are not
loaded, the database having no field for them: `ncbi_contig_l50`,
`ncbi_component_count`, `ncbi_geo_loc_name` and `ncbi_metagenome_source`.
`ncbi_isolation_source` and `ncbi_lat_lon` are, since 0.1.55; loading the table
sets them, as every field it loads, to NULL for every genome it gives no value,
unless `--do_not_null_field` is given.

`parse_ncbi_genome_category` (named `ncbi_genome_category` until 0.1.62) writes
the genomes of a genome_dirs file (`-g`) that NCBI
marks as a MAG, a SAG or an environmental genome, with where that was found, to
`ncbi_genome_category.tsv.gz` in `-o/--out_dir`, as `parse_ncbi_dir` writes its
table: `genome_id`, `ncbi_genome_category` (`derived from metagenome`, `derived
from single cell`, `derived from environmental sample`) and `source`, in
accession order. `update_metadata_db --input_folder` loads it by that name
against `metadata_ncbi_genome_category.desc.tsv`, `source` not being loaded;
until 0.1.62 it went to whatever file `-o` named, loaded with `--metadata_table`. It takes
the assembly summaries the release was selected from as `-n`, as
`select_genomes` and `ncbi_strains` do (`--genbank_assembly_summary` and
`--refseq_assembly_summary`, one of each, until 0.1.62).

```bash
gtdb_migration_tk parse_ncbi_genome_category -g genome_dirs.tsv \
    -n ncbi/assembly_summary_{archaea,bacteria}_{refseq,genbank}.txt.gz \
    -o parse_ncbi_genome_category \
    -l parse_ncbi_genome_category/parse_ncbi_genome_category.log -c 16
```

The summaries' `excluded_from_refseq` decides most genomes (814,000 of r237's
1,346,119); a user genome (`U_`) is a MAG. The rest are looked for in the head
of their GenBank file's first record, its header and source feature, on `--cpus`
processes, and the lines a category was found on are written beside the table,
`ncbi_genome_category_evidence.tsv`, which `update_metadata_db` passes over. In r226 that search added 104 SAGs, every one found there.
Over 5,103 r237 genomes, those 104 among them, the table is the one 0.1.61 wrote.
Measured on 16 processes: 157 genomes/s where reading the whole file managed 88,
about an hour for r237's 532,119. Before 0.1.62 a run took 8 to 13 hours a part
in r226, most of it spent pickling the summaries' category of every genome once
for each genome handed to a worker.

An `excluded_from_refseq` naming a `derived from` the command does not know
refuses the run before any GenBank file is read, naming every such value, exiting
1; add it to `SUMMARY_PHRASES` or `NOT_A_CATEGORY`. A genome that cannot be read
stops the run. The table is written as `ncbi_genome_category.tsv.gz.partial` and renamed once
every genome is decided, so a failed run leaves no table, or the one an earlier
run wrote. A genome with no GenBank file, one in none of the summaries, and one
marked both a MAG and a SAG (written as a MAG) are each counted in a closing
WARNING.

### Taxonomy

| Command | Description |
| --- | --- |
| `propagate_gtdb_taxonomy` | Propagate GTDB taxonomy to the new release |
| `propagate_curated_taxonomy` | Propagate curated taxonomy from representatives to their clusters |
| `update_propagated_tax` | Push propagated taxonomy to the new database |
| `add_taxonomy_to_database` | Update the taxonomy in the database |
| `set_gtdb_domain` | Set missing GTDB domain information from the NCBI domain |
| `curation_lists` | Lists and pseudo-trees for curation review |

`propagate_gtdb_taxonomy` carries the previous release's GTDB taxonomy and
species representatives to the genomes of the database (`--db_service`, or
`--hostname -u -d -p`). It reads the previous release from its GTDB metadata
files (`--gtdb_metadata_prev`, one or more, gzipped or not: a release's
`ar53_metadata` and `bac120_metadata`) and the genomes and the taxonomy they
hold from the database, where it read them from a metadata file exported by
hand (`--gtdb_metadata_cur`) until 0.1.63. It writes two headerless tables,
gzipped, in `-o/--out_dir`, which `update_propagated_tax` reads from the same
directory, `-i/--input_dir` (`-t` and `--rep_file`, which named each file, until
0.1.63):

| File | Rows |
| --- | --- |
| `gtdb_taxonomy_propagated.tsv.gz` | `genome_id`, taxonomy: each genome inheriting a taxonomy |
| `gtdb_sp_reps_propagated.tsv.gz` | `genome_id`, `True` or `False`: every genome of the database, whether it inherits representative status |

```bash
gtdb_migration_tk propagate_gtdb_taxonomy --db_service gtdb_r237 \
    --gtdb_metadata_prev ar53_metadata_r232.tsv.gz bac120_metadata_r232.tsv.gz \
    -o propagate_taxonomy -l propagate_taxonomy/propagate.log
gtdb_migration_tk update_propagated_tax --db_service gtdb_r237 -i propagate_taxonomy \
    -l propagate_taxonomy/update.log
```

`update_propagated_tax` writes every rank of every genome the taxonomy names,
and every genome's representative status, in one transaction. A genome the
taxonomy does not name keeps its domain, which `set_gtdb_domain` gives a genome
GTDB has not classified, and its ranks below the domain are set to NULL, each
rank's count logged as a WARNING, so no genome keeps a taxonomy the previous
release did not give it; over r237 that is none. Until 0.1.63 it took
`--truncate_taxonomy`, which set every genome of a metadata file exported by
hand (`-m`) to its domain first, and `--genome_list`, which kept the writes to
the genomes of the same export; the files name only genomes of the database.

A genome inherits from the genome of the previous release with its accession or
else its canonical accession, `G` and the nine digits: a new version, a move
between GenBank and RefSeq, or both. The log counts each, and the genomes of
the previous release not in the database with the representatives among them.
A run where one canonical accession names two genomes of the previous release,
or of the database, is refused, as is one where the previous release names a
genome twice, so no genome has more than one predecessor. The database's
taxonomy of a genome of the previous release must be the previous release's, or
not yet set (no rank below the domain holding more than its prefix: NULL, or the
bare `p__` ... `s__` `update_propagated_tax --truncate_taxonomy` writes); every one that is not is
listed in `gtdb_taxonomy_mismatches.tsv` beside the log, nothing is written, and
the run exits 1. Against `gtdb_r237_dev` with r232's two files it takes about a
minute: 900,882 genomes inherit a taxonomy, none differs, and both tables are the
ones 0.1.62 wrote from an export of the same database.

### Database

| Command | Description |
| --- | --- |
| `update_db` | Update the GTDB PostgreSQL database |
| `update_checkm_db` | Import CheckM estimates |
| `update_checkm2_db` | Import CheckM2 estimates |
| `update_metadata_db` | Update metadata in the database |
| `update_reps_db` | Update species cluster representatives |
| `update_ncbi_tax_db` | Update NCBI organism names and taxonomy (`update_ncbitax_db` until 0.1.57) |
| `update_taxid_to_db` | Add the NCBI taxid for each rank of each genome |
| `update_type_designation` | Update `type_designation` once SeqCode, NCBI and LPSN data are loaded |
| `add_surveillance_genomes` | Add surveillance genomes to a GTDB table |

Every command here, and `add_taxonomy_to_database`, `update_propagated_tax`,
`set_gtdb_domain` and `lpsn add_metadata`, is told which database to use by a
PostgreSQL connection service, one for each GTDB release, named in
`~/.pg_service.conf`, with its password in `~/.pgpass`:

```ini
# ~/.pg_service.conf
[gtdb_r237]
host=<host>
port=5432
dbname=<db>
user=<user>
```

```
# ~/.pgpass -- chmod 600, or libpq will not read it
<host>:5432:<db>:<user>:<password>
```

and given `--db_service gtdb_r237`, which `psql service=gtdb_r237` checks. A new
release is a new section, `[gtdb_r<release>]`, and its own line in `~/.pgpass`,
so a command names the release it is working on. `--hostname`, `-u`, `-d` and
`-p` still work, without a service or overriding what it says; without `-p` the
password is read from `~/.pgpass`. A password on the command line is shown by
`ps` to every user of the machine for as long as the run lasts.
`PGSERVICE=gtdb_r237` in the environment does what `--db_service gtdb_r237` does.

Each of these commands is one transaction, committed when it is done: a run that
fails leaves the database as it was, and exits non-zero. Once a field is
written, it is set to NULL for every genome holding a value the run did not give
it, in that same transaction, unless `--do_not_null_field` is given: the same
end as setting it to NULL for every genome first, without rewriting the rows of
the genomes given a value, which `upsert()` rewrites anyway (setting
`metadata_ncbi.ncbi_organism_name` to NULL for every genome took four minutes
for r237, a new version of each of the table's rows; it now touches none).
Nothing asks first, but
`update_metadata_db` given `--genome_list` (below). A genome the database
does not hold refuses the run (a metadata table, CheckM results, a taxonomy or a
cluster file of another release, or `update_db` not yet run), except in
`update_ncbi_tax_db`, whose NCBI files cover every assembly NCBI holds and whose
other genomes are skipped. Either way each one is listed in
`unknown_genomes.<table>.<field>.tsv` beside the log.

`update_metadata_db` loads the tables other commands write: every `.tsv` and
`.tsv.gz` of `-i/--input_folder`, each a table it knows by its name without
`.gz` (`metadata_*.tsv` from `create_tables`, `ncbi_assembly_summary.tsv` from
`parse_ncbi_assemblies`, `ncbi_assembly_metadata.tsv` from `parse_ncbi_dir`,
`ncbi_genome_category.tsv` from `parse_ncbi_genome_category`,
`strain_summary_file.tsv` from `ncbi_strains`, which all write them gzipped) and
loads against its descriptions in `data_files/table_description/`, or one table,
gzipped or not, given as `--metadata_table` with `--metadata_table_desc`. A
folder holding a table it does not know, a table both gzipped and not, or a
table without its description, is refused, every one named, before anything is
written; `substrains.tsv` and `ncbi_genome_category_evidence.tsv`, the reports
`ncbi_strains` and `parse_ncbi_genome_category` write beside their tables, are
passed over. A column no description names is not loaded, and the log says which.

It loads the genomes the database holds, `genomes.id_at_source`, skipping and
counting the rest, since a table may cover more (`parse_ncbi_assemblies`' covers
every genome of NCBI's summaries). `--genome_list` loads only the genomes in the
first column of a file instead, and a genome of it the database does not hold
refuses the run. Given `--genome_list` without `--do_not_null_field`, it asks
`[y/n]` before reaching the database: every field loaded is written only for
those listed and set to NULL for every other genome, removing its metadata. A run with no terminal to answer ends there, having changed
nothing.

A table is read once, in chunks of 100,000 genomes, each chunk written before
the next is read, so what is held does not grow with the table. The fields a
table gives one database table are written together, one UPDATE a chunk, so
each row is rewritten once: they were written a field at a time through
`upsert()`, which rewrote every row once for each field, so r237's
`metadata_gene.tsv` took six minutes for three fields and the 35 fields of
`ncbi_assembly_metadata.tsv` would each have rewritten `metadata_ncbi`, 1.4 GB.
A field a genome is given no value is left as it was, and a genome with no row
in the table is given one, as `upsert()` did.

r237's 3.87M-genome `ncbi_assembly_summary.tsv` reads in 14 seconds at a peak of
600 MB, most of it the database's genomes, where it was read whole into 8.5 GB.
Every chunk is in the run's one transaction, so a run stopped part way leaves
the database as it was and can be run again. A table is refused, and nothing
written, where a row has more or fewer columns than the header, a genome is
named twice, or an INT field holds anything but a whole number (`12.0` is
written as 12, `12.5` refused).

A whole number too big for PostgreSQL's `integer` is not written: the genome's
field is left without a value, as an empty one is, so it is NULL unless
`--do_not_null_field` keeps what it held. Each is counted in a WARNING and
listed in `int_out_of_range.<table>.tsv` beside the log. It failed `upsert()`'s
cast, rolling back the run with a traceback: r237 holds GCA_964261755.1, NCBI's
`UC_feces_MAGs_combined`, 9.5 Gbp of MAGs submitted as one assembly (322,603%
contamination by CheckM), whose `genome_size`, `gc_count`, `coding_bases`,
`ncbi_total_length` and `ncbi_ungapped_length` do not fit.

`update_ncbi_tax_db` writes each genome's NCBI organism name
(`metadata_ncbi.ncbi_organism_name`), standardised taxonomy
(`metadata_taxonomy.ncbi_taxonomy`) and unfiltered taxonomy
(`ncbi_taxonomy_unfiltered`) from the three files `parse_ncbi_taxonomy` writes:

```bash
gtdb_migration_tk update_ncbi_tax_db --db_service gtdb_r237 \
    -n ncbi_r237_organism_names.tsv --filtered ncbi_r237_standardized.tsv \
    --unfiltered ncbi_r237_unfiltered_taxonomy.tsv \
    -o update_ncbi_tax_db -l update_ncbi_tax_db/update_ncbi_tax_db.log
```

It writes the genomes the database holds, or with `--genome_list` those listed,
asking first as `update_metadata_db` does; a genome is matched however a file or
the list names it (`RS_GCF_...` or `GCF_...`). Each such genome without one of
the three -- no line in its file, or an empty value -- is written to
`ncbi_tax_missing.tsv` in `-o/--out_dir` (`genome_id`, then each field, `missing`
or empty), written whether or not there are any, and counted in a WARNING per
field. For r237 that is 1,486 genomes without a taxonomy, those whose taxid NCBI
deleted. A file naming a genome twice refuses the run. Each file is written in
chunks, as `update_metadata_db` writes a table.

`update_db` brings the NCBI genomes of the `genomes` table into line with a
release, from the two files `update_genomes` wrote for it: `report.log` says what
became of each genome and `genome_dirs.tsv` where the release keeps it. One run
covers RefSeq and GenBank, each genome's source taken from its accession.

```bash
gtdb_migration_tk update_db --db_service gtdb_r237 \
    -g release237/genome_dirs.tsv --report release237/report.log -f <YYYY-MM-DD> \
    -o update_db -l update_db/update_db.log -c 32 --dry_run
```

| Outcome in `report.log` | Genome in the database | Genome not in the database |
| --- | --- | --- |
| `new` | | added, or, where it is a new version of a genome being removed, takes over that genome's row (keeping its id and genome-list memberships, losing its aligned markers); not added where its directory has no genomic FASTA (`no genomic FASTA`), a release built before 0.1.41 holding such genomes as `new` |
| `genomic FASTA file changed` | paths and hashes updated, aligned markers deleted | added |
| `genomic FASTA sequences unchanged` | genomic FASTA hash updated | added |
| `genomic FASTA file unchanged` | paths updated if moved; hashes only with `--rehash_all` | added |
| `removed`, `to_curate` | deleted, with its metadata, aligned markers and genome-list memberships | nothing to do |

A genome already in the database whose genomic FASTA has gone from the release
stops the run: that release has lost a file. The files hashed are those that are
new or may have changed; `--rehash_all`
hashes every genome of the release, which the r237 update needs once: the
database holds SHA-256 for genomes added before 2022 and SHA-1 since, and about
6.6% of its genomic hashes no longer match their file. The hash, in both
`*_sha256` columns, is SHA-1 of the decompressed file. `has_changed` is set TRUE
for the genomes added, versioned or whose sequences changed, and FALSE for every
other NCBI genome. The files are hashed first and the update then made in one
transaction; `--dry_run` makes it and rolls it back.

| File in `--out_dir` | |
| --- | --- |
| `update_db_report.tsv` | `genome_id`, `report_outcome`, `action`, `detail` for every genome of the report and of the database, for the latest run |
| `genome_lists_affected.tsv` | each curated genome list a deleted genome was in, with `run_started` and `dry_run`; appended by every run, and on disk before the commit |
| `update_db_hashes.tsv.gz` | every hash made, with the file's size and modification time, read back by the next run |

Nothing reaches the database before the one `COMMIT`, so a run stopped at any
point -- an error, Ctrl-C, a machine reset -- leaves it as it was or wholly
updated, and is restarted by running the same command again:

- the server aborts the transaction of a client that vanished once it has sat
  idle for 10 minutes (`idle_in_transaction_session_timeout`), releasing its
  locks, so a restarted run does not wait behind them;
- the files hashed before the stop are not hashed again: hashes are flushed to
  `update_db_hashes.tsv.gz` every 30 seconds and taken from it while a file's
  size and modification time are unchanged, so a run after a dry run hashes
  almost nothing;
- a run after one that committed -- when a crash left it unknown whether the
  commit landed -- changes nothing more: a `new` or changed genome the database
  already holds with the same files is `already updated`, keeping its aligned
  markers, and `has_changed` is TRUE for the genomes this release brought on
  every run, given the same `-f`.

`--checkm_profile`, `-r/--repository` and `--report_dir` are gone: the genomes to
add were the rows of a CheckM profile, which left out what CheckM did not
assess, the run was made once per database, and the report directory is
`--out_dir`.

`update_checkm_db` writes the estimates `checkm` made for the release to
`metadata_genes`, from its release files, once `update_db` has added the
genomes:

```bash
gtdb_migration_tk update_checkm_db --db_service gtdb_r237 \
    -f checkm/checkm.profiles.tsv.gz -q checkm/checkm.qa_sh100.tsv.gz \
    -n checkm/checkm_not_assessed.tsv -o update_checkm_db -l update_checkm_db.log
```

| Field of `metadata_genes` | From |
| --- | --- |
| `checkm_completeness`, `checkm_contamination`, `checkm_strain_heterogeneity`, `checkm_marker_lineage`, `checkm_genome_count`, `checkm_marker_count`, `checkm_marker_set_count` | `checkm.profiles.tsv.gz`: Completeness, Contamination, Strain heterogeneity, Marker lineage, # genomes, # markers, # marker sets |
| `checkm_strain_heterogeneity_100` | `checkm.qa_sh100.tsv.gz`: Strain heterogeneity |

Columns are found by name, and each genome is written as the accession of its
bin. Tables that are not of the same genomes, or that name a genome twice or
lack a column, are refused before the database is asked. `-m/--metadata` is
gone: it wrote only the genomes of a table exported from the database, dropping
any other with a line on the console, where a genome the database does not hold
now refuses the run.

`update_db` keeps the row of a genome whose sequences changed or that a new
version took over, and with it the CheckM estimates of what the genome was. A
genome `checkm` left out (`-n`) that this update changed (`has_changed`) and that
holds estimates has the eight fields above set to NULL in the same transaction,
with a warning; each is written, with the estimates it held, to
`checkm_estimates_cleared.tsv` beside the log, which is written, its header
alone, where none was cleared. A genome the not-assessed file names that this
update did not change is warned of and kept. Only the genomes that file names
are cleared, so CheckM tables of a few genomes loaded on their own clear nothing
of the rest.

Once the estimates are written and cleared, every RefSeq and GenBank genome in
the database with no completeness that the not-assessed file does not name is
warned of and written to `checkm_estimates_missing.tsv` in `-o/--out_dir`
(`genome_id`, `status`). The file is written, header alone, when there are none.
These are genomes `checkm` never planned, such as one `update_db` added that
`report.log` calls unchanged. `status` is read from the genome's row of `genomes`:
`new` (`has_changed`, and `date_added` on the day of `last_update`: added, or a
new version that took over its predecessor's row, as `report.log` calls it new),
`updated` (`has_changed`, added earlier: its sequences changed under its
accession) or `unchanged`.

`update_checkm2_db` does the same for the estimates `checkm2` made, from its
release files:

```bash
gtdb_migration_tk update_checkm2_db --db_service gtdb_r237 \
    -f checkm2/checkm2.quality_report.tsv.gz -n checkm2/checkm2_not_assessed.tsv \
    -o update_checkm2_db -l update_checkm2_db.log
```

| Field of `metadata_genes` | From `checkm2.quality_report.tsv.gz` |
| --- | --- |
| `checkm2_completeness` | Completeness |
| `checkm2_contamination` | Contamination |
| `checkm2_model` | Completeness_Model_Used, as `General` or `Specific` |

CheckM2 names the model it estimated completeness with, e.g. `Neural Network
(Specific Model)` or `Gradient Boost (General Model)`; the database holds
`Specific` or `General`. A report naming a model that is neither is refused,
naming the model and its genomes, and nothing is written. A genome `checkm2` left
out (`-n`) that this update changed and that holds CheckM2 estimates has the
three fields set to NULL, each written with what it held to
`checkm2_estimates_cleared.tsv` beside the log, and the genomes with no
`checkm2_completeness` that the not-assessed file does not name are written to
`checkm2_estimates_missing.tsv` in `--out_dir`, as for `update_checkm_db`.
`join_checkm2` and `prepare_checkm2`, which ran CheckM2 by hand before `checkm2`
did and wrote the table this command's fields were loaded from, are gone.

### Nomenclature resources

| Command | Description |
| --- | --- |
| `lpsn` | LPSN processing (`pull_html`, `parse_html`, `lpsn_wf`, `add_metadata`) |
| `bacdive` | BacDive processing (`download_strains`) — in development |
| `strains` | Year of priority (`date_table`) and type material status (`type_table`) of each genome, from LPSN |
| `ncbi_strains` | Parse NCBI assembly reports for strain identifiers and type material status |

`ncbi_strains` writes `strain_summary_file.tsv.gz` in `--out_dir`: each genome of a
genome_dirs file (`-g`) with its organism name and strain IDs, read from its
`_assembly_report.txt` on `--cpus` processes, and its NCBI type material status,
the `relation_to_type_material` of the assembly summaries (`-n`, as
`select_genomes` takes them; `--gb`, `--ga`, `--rb` and `--ra` until 0.1.55).
`update_metadata_db` loads it into `ncbi_strain_identifiers` and
`ncbi_type_material_designation`.

```bash
gtdb_migration_tk ncbi_strains -g genome_dirs.tsv \
    -n ncbi/assembly_summary_{archaea,bacteria}_{refseq,genbank}.txt.gz \
    -o ncbi_strains -l ncbi_strains/ncbi_strains.log -c 8
```

Rows are written as the workers finish them, in no particular order. A genome in
none of the summaries has an empty type material status, which
`update_metadata_db` loads as NULL, and one with no assembly report an empty
organism name and strain IDs; the log ends with a WARNING counting each. A
genome that cannot be read stops the run. The table is written as
`strain_summary_file.tsv.gz.partial` and renamed only once every genome is in it,
so a failed run leaves no table, or the one an earlier run wrote.

Genomes whose infraspecific name names a substrain (`substr.`) are written to
`substrains.tsv` in `--out_dir`, a TSV of `genome_id` and `strain_id` (the
infraspecific name as NCBI gives it, e.g. `strain=K-12 substr. MG1655`), sorted
by genome; the log gives their count and the file in one line. It is written
on every run, its header alone where there are none. It is not a table:
`update_metadata_db --input_folder` passes over it by name rather than refusing
the folder, so `ncbi_strains` may write into the folder it loads.

`strains type_table` decides which genomes of a release are assembled from type
material, matching each genome's NCBI species and strain IDs against LPSN's type
strains, and writes `gtdb_type_strain_summary.tsv.gz` for `metadata_type_material`,
which `update_metadata_db` loads gzipped as it is.
It reads the genomes from the release's files rather than the database, so it can
run once `ncbi_metadata_sync` and `update_genomes` are done: which genomes from
`genome_dirs.tsv` (`-g`); their organism name, taxid, strain IDs and NCBI type
material status from the assembly summaries (`-n`, as `select_genomes` takes them); their
species, and every name NCBI gives a taxon, from `names.dmp` and `nodes.dmp`.

It ends with one WARNING per kind of warning, its count and up to three examples;
every warning is in `type_table_warnings.tsv`, with its kind, what that kind
means, the warning, and the data that bears on it (a taxid, a genome, a genus).

```bash
gtdb_migration_tk strains type_table -g release237/genome_dirs.tsv \
    -n ncbi/assembly_summary_*.txt.gz \
    --ncbi_names taxdump/names.dmp --ncbi_nodes taxdump/nodes.dmp \
    --lpsn_gss_file lpsn_gss_<date>.csv --lpsn_dir lpsn/parse_html/all_ranks \
    --year_table year_table.tsv -o strain_table -c 16
```

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
MARKER_DIR_SUFFIX = {'pfam': '33.1_lite', 'tigrfam': '15.0_lite'}
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

`MARKER_DIR_SUFFIX` is the default `--dir_suffix` of `hmmsearch` and
`top_hit`, so by default they write `prodigal/pfam_33.1_lite/` and
`prodigal/tigrfam_15.0_lite/`; pass `--dir_suffix` only to annotate against
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
