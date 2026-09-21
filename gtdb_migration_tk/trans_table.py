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

"""
trans_table.py -- predict the translation table of each genome with gTranslate.

gTranslate is run as a SUBPROCESS rather than imported. It requires Python >= 3.12
where this toolkit supports 3.8, and it is installed in an environment of its own
holding pinned versions of scikit-learn, xgboost and lightgbm that the predictions
depend on; importing it would make those the toolkit's own dependencies and its
floor the toolkit's floor. What crosses between the two is a file of genomes and a
file of predictions, so a subprocess is all the coupling the work needs.

BATCHES
The release is cut into batches of --batch_size genomes, each with a directory of
its own under --out_dir, because a release takes days to predict and nothing that
takes days should have to be started again from the beginning. A batch is the unit
of work, of restart and of sharing between machines: several machines may be given
the same --out_dir and will divide the release between them without being told
which part to take, and a machine that resets loses the batch it was on rather
than the run.

The batches are settled BEFORE any of them is processed, and are the plan every
machine works from. The genomes are sorted by accession first, so which genomes
are in batch N follows from the set of genomes alone and not from the order a
genome_dirs file happens to be written in. Once the batchfiles exist they are
authoritative: a later run reuses them and says so rather than partitioning the
release again, since a second partition of a release that has gained a genome
would move genomes between batches that are already finished.

The plan is cut from the genome_dirs file and nothing else. Whether a genome's
FASTA is actually on disk is asked of each batch as that batch is run, not of
the release beforehand: gTranslate checks the paths of a batchfile itself and
refuses the whole batch if one is missing, so the check has to happen, but asked
of 1.35M genomes over NFS it is hours in which nothing is written and nothing
can be resumed, where asked of ten thousand it is a minute against a batch that
then runs for hours. A genome left out is named in the batch's own directory.

CANARIES
A batch directory carries its state in files, because a file is what two machines
can both see and what survives the process that wrote it:

    RUNNING    a machine is working on this batch, and names itself, its PID
               and when it started
    PREDICTED  gTranslate has run over this batch and its predictions are final
    SUCCESS    the batch finished and its results are complete
    FAILED     the batch was attempted and something went wrong

RUNNING is created by linking a uniquely named file to it, not by opening it with
O_EXCL: the release lives on NFS, where link() is the operation that is atomic
across machines. A batch already claimed is left to the machine that claimed it,
and reported at the end rather than waited for. A batch that FAILED is retried by
the next run without a flag, the machine that failed it having cleared its own
claim; what it wrote is kept as FAILED.<timestamp> rather than deleted, so the
retry does not erase the record of why the batch failed in the first place.

THE CLAIM IS A LEASE
The machine holding a batch touches its RUNNING file every HEARTBEAT_SECONDS
while it works, and a claim that has not been touched for CLAIM_LEASE_SECONDS is
taken by whichever machine next comes to the batch. A claim of this host's whose
PID has gone is still taken at once, that being certain rather than inferred, but
it is no longer the only way a batch is freed: PIDs are reused, so after a reset
the PID a claim names is as likely to belong to something else as to be missing,
and a batch would then be skipped as busy by every machine forever. The lease
also frees a batch held by a machine that is running but wedged -- a process hung
on a dead NFS mount cannot touch its own claim any more than a dead one can,
which is the case a liveness check by PID reads exactly backwards.

The age of a claim is measured against the FILE SERVER's clock and not against
this machine's: the lease is compared with the time a file created in the same
directory is given, so nothing depends on two machines agreeing about the time.
--reclaim remains, for taking a claim before its lease is up.

WHAT IS DONE IS NOT DONE AGAIN
A batch's plan is kept gzipped, and gTranslate cannot read it: it opens a
batchfile with a plain open(). So the file gTranslate is handed is written plain
before it starts and removed once it has finished, and what a finished batch
keeps is the compressed plan alone. A batch taken over before gTranslate finished
writes that copy again from the plan, which is why removing it costs nothing, and
a batch gTranslate failed on keeps it, being what a retry looks at to see what
went in. The plan was not always compressed, so where a batch holds the older
uncompressed one that is read instead: a finished output directory is still read,
the command being run again over one to pick up the work since added to it, and a
release whose batches look unplanned would be partitioned again with its batches
already done.

gTranslate is the hours of a batch; the comparison that follows is seconds. So
the prediction step records PREDICTED the moment gTranslate returns 0, and a
batch reclaimed after that compares what is already there rather than predicting
it again. Within a batch gTranslate resumes by itself: it writes each genome's
called genes with a checksum beside them and skips a genome whose files verify,
and it clears those intermediates only when it exits cleanly, so a batch
interrupted at genome 7,000 of 10,000 carries on from genome 7,000. Nothing here
removes a batch directory before retrying it, for that reason.

A GENOME THAT CANNOT BE PREDICTED IS NOT A REASON TO LOSE THE BATCH
gTranslate calls genes in worker processes and ends the whole run when one of
them dies, so a single genome it cannot handle -- a 255 bp fragment with no genes
to count codons in, an assembly Prodigal refuses for its runs of N -- takes the
other 9,999 with it, and takes them again on every retry. --force, which gTranslate
answers by dropping such a genome and carrying on, is therefore passed unless
--no_force says not to. What it drops is what a batch cannot say anything about,
so the accessions gTranslate returned no prediction for are written to the batch
directory as no_prediction.tsv rather than being left to be discovered by prodigal
refusing the release.

THE LOG OF A BATCH IS KEPT WITH THE BATCH
Every machine writes what it does to its own --log, and several machines sharing
one log file over NFS do not append to it, they overwrite one another and leave
the file full of holes. So what happens to a batch is also written to
trans_table.log in the batch's own directory, which no other machine writes to:
whatever became of the run that started a batch, the batch says.

THE COMPARISON
gTranslate predicts a table; NCBI declares one in the GFF it serves for a genome
that it has annotated. Each batch compares the two wherever both are known and
writes ncbi_tt_conflict.tsv, which holds the genomes they DISAGREE about and
carries for each the coding densities the prediction was made from, the genome's
NCBI taxonomy, the table the coding density rule alone would have chosen, and
whether that rule and gTranslate disagree about the genome being recoded at all.

Two files come out of it. gtranslate_ncbi_tt_comparison.tsv.gz holds every genome
the two BOTH called, agreements and all, because the rate they differ at and
whether the genomes they differ about are unlike the ones they agree about are
questions the agreements have to be present to answer. ncbi_tt_conflict.tsv holds
the disagreements alone, which is the few hundred rows of a release worth looking
at: a genome whose genes GTDB would call under a table NCBI does not agree with
is something to go and read, and it should not have to be filtered out of a
million rows first. The conflicts are the comparison FILTERED rather than a
second walk over the genomes -- every GFF has already been read once -- so the
two can never come to disagree about which genomes conflicted.

The comparison is gzipped and the conflicts are not, in the batch as well as in
the release: one is a row per genome compared, 786,144 of r237's 1.35M, and the
other is meant to be opened and read. A genome NCBI has not annotated is in
neither, having nothing to compare against and no ncbi_conflict to report; how
many there were is logged, with how many agreed and how many conflicted, for
every batch. The table NCBI declares is read by
ncbi_utils.ncbi_translation_table(), which prodigal reads it with too.

Both files carry the genome's GC, N50 and size as gTranslate measured them, taken
from its summary. They are what tells a conflict that is a real disagreement
about a real genome from one about 200 kb of something barely assembled, and
measuring them here would be reading 1.3M genomes a second time to learn what has
already been learnt.

Once every batch has succeeded the run also writes the conflicts and the
prediction summary for the whole release at the top of --out_dir, so that a
release finished across several machines is one file to read. The release summary
is gzipped, being a row per genome -- 116 MB of text for r237 and 34 MB
compressed -- while the conflicts are a few hundred rows meant to be looked at
and are not. A batch's own summary is gTranslate's output, written by gTranslate
and not this command's to compress. Which of them prodigal is handed does not
matter: read_translation_table_summary() decides by the file's first two bytes
rather than by its name, so a release summary that someone gunzipped, or renamed
on the way to another machine, still reads.

WHAT THE RATE IS OF
Every batch says how many of its genomes NCBI declares a table for and what share
of THOSE gTranslate disagrees about, and the release says the same about all of
them together. The denominator is the genomes that could be compared and not the
genomes of the release: a genome NCBI has not annotated is not one the two agree
or disagree about, and counting it in would turn the number into a measure of how
much of the release NCBI has annotated. The release figure is added up from the
batches' SUCCESS canaries rather than recomputed, because a release is predicted
by several machines and the machine running the last batch has compared none of
the others; the canary is what every machine leaves behind.

THE QUALITY OF A CONFLICTING GENOME
A conflict is a genome two callers disagree about, and the question it raises is
which of them is right. Completeness is the evidence: genes called under the
wrong code are truncated at every TGA, and the markers CheckM2 counts go with
them. So each conflicting genome is put to CheckM2 TWICE, once under gTranslate's
table and once under NCBI's, and ncbi_tt_conflict.tsv carries both answers --
one run would say how good the genome is, two say which table makes it look like
a genome at all. CheckM2's own choice is not asked for: left to itself it picks
between tables 4 and 11 by coding density, which is the rule checkm_tt already
reports and which cannot express 25.

Beside each pair of numbers is whether they pass standard GTDB QC: more than 50%
complete, less than 10% contaminated, and a quality score of completeness less
five times contamination above 50. The verdict is the one the release is actually
kept or dropped by, so a conflict that changes it is a different thing from one
that moves the numbers a little, and the row says which it is without the reader
doing the arithmetic. It is given per table because a genome can pass under one
and fail under the other, which is exactly the case worth looking at. A genome
CheckM2 returned nothing for is na and not False: it was not looked at and found
wanting.

The runs are made once for the RELEASE and grouped by table, not once per batch.
The conflicts are a few hundred genomes of a million-odd -- two or three per batch
-- and CheckM2 loads its models and searches the whole DIAMOND database once per
run whatever the run holds, so the cost is the number of runs. Grouped by table
the whole release is one run per table in dispute, which is three; per batch it
would be hundreds. That also means a release already predicted picks this up by
running the command again: every batch is SUCCESS and is skipped, and the work
happens where the release files are written.

The directory those runs work in is removed once the estimates are in the
conflict file, which is where they were wanted. What it holds is the staged
links, the called proteins and the DIAMOND output -- about 600 KB per genome per
table, a few hundred megabytes for a release -- and none of it says anything the
conflict file does not now say. What that gives up is a later run reading the
reports rather than making them again, which is minutes for a few hundred
genomes, against the days the batches protect: which is why the batches are never
removed and this is. A table CheckM2 produced nothing for keeps the directory, so
that retrying it does not also redo the tables that worked.

CheckM2 is not installed beside this toolkit -- the TensorFlow it needs wants an
icu the release environment cannot hold -- so it is an external program like
gTranslate, found on PATH, and is checked for before a run starts. It names a
result for the file it read and NCBI names the file for the assembly, so each
genome is linked as <accession>.fna.gz before the run: the link is what joins the
report back to the conflict rows. Nothing about this step can cost the release
the comparison, which is done, written and counted before it starts -- a genome
with no FASTA, a run that fails, a table CheckM2 returns nothing for all leave
their genomes with na in those four columns and the rest of the row intact.
"""

import contextlib
import datetime
import gzip
import logging
import os
import shutil
import socket
import subprocess
import tempfile
import threading
import time
import uuid
from concurrent.futures import ThreadPoolExecutor
from typing import (Dict, Iterator, List, NamedTuple, Optional, Sequence,
                    Tuple)

from tqdm import tqdm

from gtdb_migration_tk.biolib_lite.external.execute import check_dependencies
from gtdb_migration_tk.biolib_lite.common import canonical_gid
from gtdb_migration_tk.ncbi_utils import (GENOMIC_FASTA_EXT, NCBI_NA,
                                          genomic_gff, ncbi_translation_table)
from gtdb_migration_tk.utils.common import (TT_SUMMARY_DENSITY_4,
                                            TT_SUMMARY_DENSITY_11,
                                            TT_SUMMARY_GC,
                                            TT_SUMMARY_GENOME_SIZE,
                                            TT_SUMMARY_N50,
                                            TT_SUMMARY_TABLE,
                                            checkm_translation_table,
                                            open_text,
                                            read_translation_table_summary)


# The executable, looked up on PATH rather than given a path of its own: it is
# installed as a module here, and the module is what puts it on PATH along with
# prodigal and GTRANSLATE_MODEL_PATH, none of which this command can supply.
GTRANSLATE_BIN = 'gtranslate'

# The subcommand run. gTranslate's other subcommands either train models or plot
# what one predicted, neither of which belongs in a migration.
DETECT_TABLE = 'detect_table'

# Genomes per batch. Large enough that the cost of starting gTranslate and loading
# its classifiers is nothing against the genomes it then processes, small enough
# that a machine lost mid-batch costs hours and not days.
DEFAULT_BATCH_SIZE = 10000

# Batch directories are numbered rather than named for the genomes they hold: the
# accessions of a batch are in its batchfile, and a name is a thing to sort by.
BATCH_DIR_PREFIX = 'batch_'
BATCH_DIR_FORMAT = BATCH_DIR_PREFIX + '{:06d}'

# Written into the batch directory rather than a temporary one: it is the record
# of which genomes the batch is, and it is what a later run and another machine
# read to agree on that without partitioning the release again.
BATCHFILE_NAME = 'gtranslate_batchfile.tsv.gz'

# What earlier versions called it, written uncompressed. A finished output
# directory is still read -- the command is run again over one to pick up work
# added since -- and batch_dir_names() finding no plan there would repartition a
# release whose batches are already done.
LEGACY_BATCHFILE_NAME = 'gtranslate_batchfile.tsv'

# What gTranslate is handed when some genome of the batch has no FASTA to process,
# and the accessions left out of it. Written only in that case, so the file being
# there at all says a batch had something wrong with it. BATCHFILE_NAME stays the
# record of which genomes the batch IS, which is what the comparison reads.
# The batchfile actually handed to gTranslate, which reads it with a plain
# open() and cannot take the gzipped one. It is written before gTranslate starts
# and removed once it has finished, so the uncompressed copy exists only while
# the batch is being worked on and the plan kept for good is the compressed one.
PRESENT_BATCHFILE_NAME = 'gtranslate_batchfile_present.tsv'
MISSING_NAME = 'missing_genomic_fasta.tsv'

# The accessions gTranslate was given and returned no prediction for, which is
# what --force leaves behind. Written only when there are some, so the file being
# there at all says a batch holds genomes no table was predicted for.
NO_PREDICTION_NAME = 'no_prediction.tsv'

# What this command does to a batch, written in the batch's own directory as well
# as to --log: every machine sharing an --out_dir writes its own log, and one log
# file appended to from several machines over NFS holds none of them.
BATCH_LOG_NAME = 'trans_table.log'

# The state of a batch, held in files so that another machine can see it and so
# that it outlives the process that wrote it.
RUNNING_CANARY = 'RUNNING'
PREDICTED_CANARY = 'PREDICTED'
SUCCESS_CANARY = 'SUCCESS'
FAILED_CANARY = 'FAILED'

# How often the machine holding a batch says it is still there, and how long a
# claim outlives the last thing said. The interval is small against the hours a
# batch takes and the lease is large against the interval, so a claim is freed
# only by a machine that has genuinely stopped touching it, not by one whose
# heartbeat was late to reach the file server.
HEARTBEAT_SECONDS = 300
CLAIM_LEASE_SECONDS = 2 * 60 * 60

STATE_PENDING = 'pending'
STATE_RUNNING = 'running'
STATE_SUCCESS = 'success'
STATE_FAILED = 'failed'

# gTranslate names its summary for its --prefix, which this command passes through.
DEFAULT_PREFIX = 'gtranslate'
SUMMARY_SUFFIX = '.translation_table_summary.tsv'

# The release's summary is gzipped and a batch's is not. A batch's is gTranslate's
# own output, written by gTranslate into the batch directory and not this
# command's to name; the release's is the concatenation of all of them, which for
# r237 is 116 MB of text and 34 MB compressed. It is read back by
# read_translation_table_summary(), which decides by the file's first two bytes
# rather than by its name, so prodigal takes either.
GZIP_EXT = '.gz'

CONFLICT_NAME = 'ncbi_tt_conflict.tsv'

# The comparison of every genome the two BOTH called, agreements and all. It is a
# row per genome NCBI declares a table for -- 786,144 of r237's 1.35M -- which is
# why it is gzipped in the batch as well as in the release, and why the conflicts
# keep a file of their own: the few hundred rows worth looking at should not have
# to be filtered out of a million first. A genome NCBI has not annotated is not
# here at all, having nothing to be compared against and no ncbi_conflict to
# report; how many there were is logged.
COMPARISON_NAME = 'gtranslate_ncbi_tt_comparison.tsv.gz'

# ncbi_conflict is the plain inequality -- gTranslate said one table and NCBI
# another -- and is what decides whether a genome is in ncbi_tt_conflict.tsv.
# checkm_conflict beside it is the narrower question of whether the coding
# density rule and gTranslate disagree about the genome being RECODED at all:
# see checkm_conflict(). The genome statistics are gTranslate's own measurements,
# carried through rather than measured again.
COMPARISON_HEADER = ('genome_id', 'gtranslate_tt', 'ncbi_tt', 'checkm_tt',
                     'ncbi_conflict', 'checkm_conflict',
                     'coding_density_4', 'coding_density_11',
                     'gc_percent', 'n50', 'genome_size', 'ncbi_taxonomy')

# checkm_tt sits beside the other two tables rather than at the end: the three
# are the answers to one question, gtranslate_tt and ncbi_tt being the two that
# disagreed. It is the table the coding density rule alone would choose, which is
# what Prodigal and CheckM do unaided, and it cannot express table 25 at all -- a
# genome gTranslate calls 25 is one the old rule was never able to get right.
# checkm_conflict follows it, saying whether that matters: see checkm_conflict().
# There is no result column: every row of the file is a conflict.
CONFLICT_HEADER = ('genome_id', 'gtranslate_tt', 'ncbi_tt', 'checkm_tt',
                   'checkm_conflict', 'coding_density_4', 'coding_density_11',
                   'gc_percent', 'n50', 'genome_size', 'ncbi_taxonomy')

# The executable that estimates the quality of a conflicting genome, looked up on
# PATH as gTranslate is. It is not installable beside this toolkit -- the
# TensorFlow it needs wants an icu the release environment cannot hold -- so it
# lives in an environment of its own and is reached as an external program.
CHECKM2_BIN = 'checkm2'

# Where the CheckM2 runs of a release are kept, one directory per table, and the
# file each writes. The directory is per table because a genome is asked about
# under two tables and CheckM2 names its results for the genome alone.
CHECKM2_DIR = 'checkm2'
CHECKM2_TABLE_DIR = 'table_{}'
CHECKM2_REPORT = 'quality_report.tsv'

# Where the staged genomes of a table's run are put. It is NOT inside the run's
# own directory: --force empties the output directory before CheckM2 starts, so a
# staging directory under it would be deleted along with the last run's results.
CHECKM2_INPUT_DIR = 'input'

# What a staged genome is named. CheckM2 labels a result with the basename minus
# the last two extensions, so <accession>.fna.gz comes back as the accession
# while NCBI's own _genomic.fna.gz comes back as the assembly name. The extension
# is therefore part of the join between the report and the conflict rows, not a
# detail of it: see stage_checkm2_input().
CHECKM2_LINK_EXT = '.fna.gz'

# The columns of that report this reads. CheckM2 writes a dozen more -- the
# coding density, the N50, the model it chose -- which are about the run and not
# about the conflict, and are left in the report for whoever wants them.
CHECKM2_NAME = 'Name'
CHECKM2_COMPLETENESS = 'Completeness'
CHECKM2_CONTAMINATION = 'Contamination'

# A conflicting genome is asked about TWICE, once under each of the tables in
# dispute, because completeness under a table is the evidence about that table:
# genes called under the wrong code are truncated at every TGA, and the markers
# CheckM2 counts go with them. One run under one table would say how good the
# genome is; two say which table makes it look like a genome at all. Each side's
# verdict sits beside the numbers it was reached from rather than at the end, so
# that a row read by eye is two answers to one question and not four numbers.
CHECKM2_COLUMNS = ('cm2_completeness_gtranslate_tt', 'cm2_contamination_gtranslate_tt',
                   'pass_qc_gtranslate_tt',
                   'cm2_completeness_ncbi_tt', 'cm2_contamination_ncbi_tt',
                   'pass_qc_ncbi_tt')

# Standard GTDB QC: a genome is kept where it is more than half there, barely
# contaminated, and still more than half there once its contamination is charged
# against it at five times its weight. All three have to hold -- the quality
# score alone would keep a 96% complete genome carrying 9% contamination, and the
# completeness alone would keep anything that had been sequenced at all.
#
# The thresholds are exclusive as GTDB states them: a genome exactly 50%
# complete, or at exactly 10% contamination, does not pass.
QC_MIN_COMPLETENESS = 50.0
QC_MAX_CONTAMINATION = 10.0
QC_CONTAMINATION_WEIGHT = 5.0
QC_MIN_QUALITY = 50.0

# The release file carries the CheckM2 columns and a batch file does not: CheckM2
# runs once for the release, over the few hundred genomes every batch together
# found, rather than once per batch over the two or three each found on its own.
# They go before ncbi_taxonomy so that the lineage stays the last and longest
# field of the row, and so that the columns of the comparison stay together.
CONFLICT_HEADER_CHECKM2 = (CONFLICT_HEADER[:-1] + CHECKM2_COLUMNS
                           + CONFLICT_HEADER[-1:])

# The standard genetic code, and the two recoded ones gTranslate chooses between
# it and: 4 for the genomes that read TGA as tryptophan, 25 for those that read
# it as glycine. The density rule cannot express 25, only 4 and 11.
STANDARD_TABLE = 11
RECODED_TABLES = (4, 25)

# How many of the release's genomic FASTA files are asked about at once while the
# batches are planned. The question is one stat per genome and nothing else, so
# what it costs is round trips to the file server and not CPU. Measured against
# r237 over the NFS the genomes are held on, one at a time takes 12-23 ms a
# genome -- 4 to 8 hours for 1.35M genomes -- against 3-7 ms with 32 outstanding
# at once, so a few hours become one or two. Beyond 32 nothing further was
# measurable: the limit is the server and the load it is already under, not the
# number of threads asking. Threads rather than processes because a stat spends
# its time in the kernel waiting and brings back one number.
STAT_THREADS = 32

# How many genomes are handed to the pool at a time. Executor.map() submits every
# item it is given before the first result can be read, one future per genome, so
# the whole release at once builds a million-odd futures before a single answer
# comes back. A chunk is large enough that no thread waits for the next one to be
# cut and small enough to be nothing in memory.
STAT_CHUNK = 50000


def genomic_fasta(genome_dir: str) -> str:
    """The genomic FASTA NCBI serves for a genome.

    NCBI names the file for the assembly and names the genome directory the same,
    so the file is named rather than searched for -- _cds_from_genomic.fna.gz and
    _rna_from_genomic.fna.gz end the same way and are different files.

    Parameters
    ----------
    genome_dir : str
        Genome directory, of a release or of the mirror.

    @return: path of the genomic FASTA in that directory, which may not exist.
    """

    assembly = os.path.basename(os.path.normpath(genome_dir))

    return os.path.join(genome_dir, assembly + GENOMIC_FASTA_EXT)


def read_genome_dirs(gtdb_genome_path_file: str) -> List[Tuple[str, str]]:
    """Read the genomes of a release from its genome_dirs file.

    The file is the headerless accession / directory / canonical accession TSV
    list_genomes and update_genomes write. Rows are split on tabs and further
    columns ignored, as every other reader of it does, so a column appended later
    does not reach here.

    Parameters
    ----------
    gtdb_genome_path_file : str
        genome_dirs file of the release.

    @return: (accession, genome directory) for each genome, in the order read.
    """

    genomes = []
    with open(gtdb_genome_path_file) as handle:
        # leave=False: the bar is worth having while a release of half a million
        # genomes is read and worth nothing afterwards, and a bar left behind
        # sits in the middle of the log saying what has already been reported
        for line in tqdm(handle, ncols=100, leave=False, desc='Reading genomes'):
            line = line.strip()
            if not line:
                continue
            tokens = line.split('\t')
            genomes.append((tokens[0], tokens[1]))

    return genomes


def fasta_size(fasta: str) -> int:
    """Size of a genomic FASTA, and 0 where there is no file to have one.

    One stat rather than os.path.exists() and os.path.getsize(), which ask the
    file server the same question twice; over NFS, and once per genome of a
    release, that is half the cost of planning a run. A file that has gone
    between the two calls also raises from the second, where here it is simply a
    genome with nothing to process.

    Parameters
    ----------
    fasta : str
        Path of the genomic FASTA.

    @return: size in bytes, or 0 if the file is absent or cannot be read.
    """

    try:
        return os.stat(fasta).st_size
    except OSError:
        return 0


def split_by_fasta(rows: Sequence[Tuple[str, str]],
                   threads: int = STAT_THREADS) -> Tuple[List[Tuple[str, str]], List[str]]:
    """Sort a batch's genomes into those that can be asked about and those that cannot.

    A genome is asked about only where its genomic FASTA is on disk and is not
    empty. gTranslate checks the paths of a batchfile before it starts and
    refuses the WHOLE batch if one of them is missing, so a single absent file
    would cost the other ten thousand genomes of the batch -- and would cost them
    again on every retry, the batch failing identically each time. Filtering here
    leaves the batch to run and the accession to be named.

    The files are asked about many at a time, the genomes being held over NFS
    where a stat is a round trip to a server rather than a lookup in a cache.
    Nothing is computed here to be divided up: what the pool is for is having many
    round trips outstanding at once rather than one.

    The answer does not depend on how many threads asked: results are read back in
    the order the genomes were given, which is the order gTranslate is handed them.

    Parameters
    ----------
    rows : sequence of tuple
        (FASTA path, accession) for the genomes of a batch.
    threads : int
        Stat calls to keep in flight at once.

    @return: (present, missing), the rows to hand gTranslate and the accessions
             of the genomes left out.
    """

    present, missing = [], []
    # leave=False, as the bar reading the genome_dirs file is: it is worth having
    # while the batch is checked over and worth nothing once it is
    with tqdm(total=len(rows), ncols=100, leave=False,
              desc='Checking genomes') as pbar, \
            ThreadPoolExecutor(max_workers=max(1, threads)) as pool:
        for start in range(0, len(rows), STAT_CHUNK):
            chunk = rows[start:start + STAT_CHUNK]
            for (fasta, accession), size in zip(
                    chunk, pool.map(fasta_size, [fasta for fasta, _ in chunk])):
                if size > 0:
                    present.append((fasta, accession))
                else:
                    missing.append(accession)
                pbar.update()

    return present, missing


def write_batchfile(rows: Sequence[Tuple[str, str]], batchfile: str,
                    compress: bool = False) -> None:
    """Write the two-column batchfile gTranslate reads.

    The genome ID given is the accession the genome_dirs file names, so every row
    of the prediction table can be matched back to the genome directory it came
    from without canonicalising anything.

    Parameters
    ----------
    rows : sequence of tuple
        (FASTA path, accession) for each genome to ask about.
    batchfile : str
        File to write.
    compress : bool
        Write it gzipped, which the batch's own plan is and the copy handed to
        gTranslate is not: gTranslate reads a batchfile with a plain open().

    @return: None
    """

    with (gzip.open(batchfile, 'wt') if compress else open(batchfile, 'w')) as handle:
        for fasta, accession in rows:
            handle.write('{}\t{}\n'.format(fasta, accession))


def batchfile_path(batch_dir: str) -> str:
    """Where a batch's plan is, whichever version of this command wrote it.

    The plan is gzipped now and was not before, and a finished output directory
    is still read: the command is run again over one to pick up the work that has
    since been added to it. Asked of a batch that has neither, the answer is
    where the plan would be written, so that a caller's error names the file it
    was looking for.

    Parameters
    ----------
    batch_dir : str
        Batch directory.

    @return: path of the batchfile.
    """

    for name in (BATCHFILE_NAME, LEGACY_BATCHFILE_NAME):
        path = os.path.join(batch_dir, name)
        if os.path.exists(path):
            return path

    return os.path.join(batch_dir, BATCHFILE_NAME)


def read_batchfile(batchfile: str) -> List[Tuple[str, str]]:
    """Read back the genomes of a batch.

    A batch is read from its own batchfile and not from the genome_dirs file, so
    that a batch is self-contained: the machine that processes it needs to agree
    with the machine that planned it about which genomes it holds, and the
    batchfile is that agreement written down.

    Parameters
    ----------
    batchfile : str
        Batchfile of one batch.

    @return: (FASTA path, accession) for each genome of the batch.
    """

    rows = []
    with open_text(batchfile) as handle:
        for line in handle:
            line = line.rstrip('\n')
            if not line:
                continue
            fasta, _, accession = line.partition('\t')
            rows.append((fasta, accession))

    return rows


def check_batch_fastas(batch_dir: str,
                       threads: int = STAT_THREADS
                       ) -> Tuple[str, List[Tuple[str, str]], List[str]]:
    """Check over the genomes of one batch, just before gTranslate is run on it.

    The check belongs to the batch rather than to the plan. Asked of the whole
    release up front it is hours of stat calls before a single batch directory
    exists -- on r237, 1.35M genomes over NFS -- and a run stopped in it has
    nothing to resume from, the plan not yet written. Asked of a batch it is ten
    thousand stat calls against a batch that then runs for hours, it happens
    while other machines are already working, and a machine lost during it costs
    one batch. It also leaves the batch boundaries following from the genome_dirs
    file alone, rather than from what stat said on the day the plan was cut.

    The file handed to gTranslate is always written, and is always the plain one:
    gTranslate reads a batchfile with a plain open() and the batch's own plan is
    gzipped. It is removed once gTranslate has finished with it, so what a
    finished batch keeps is the compressed plan alone.

    Parameters
    ----------
    batch_dir : str
        Batch directory, holding the batchfile the plan cut.
    threads : int
        Stat calls to keep in flight at once.

    @return: (batchfile, present, missing) -- the file to hand gTranslate, the
             rows it names, and the accessions left out of it.
    """

    present, missing = split_by_fasta(
        read_batchfile(batchfile_path(batch_dir)), threads)

    handed_over = os.path.join(batch_dir, PRESENT_BATCHFILE_NAME)
    write_batchfile(present, handed_over)

    if missing:
        with open(os.path.join(batch_dir, MISSING_NAME), 'w') as handle:
            for accession in missing:
                handle.write('{}\n'.format(accession))

    return handed_over, present, missing


def summary_name(prefix: Optional[str] = None) -> str:
    """The name gTranslate gives its translation table summary.

    Parameters
    ----------
    prefix : str
        The --prefix the run passes gTranslate, or None for gTranslate's default.

    @return: filename of the summary within a batch directory.
    """

    return (prefix or DEFAULT_PREFIX) + SUMMARY_SUFFIX


def release_summary_name(prefix: Optional[str] = None) -> str:
    """The name the summary of a whole release is written under.

    Parameters
    ----------
    prefix : str
        The --prefix the run passes gTranslate, or None for gTranslate's default.

    @return: filename of the release summary at the top of the output directory.
    """

    return summary_name(prefix) + GZIP_EXT


def batch_dir_names(out_dir: str) -> List[str]:
    """The batch directories already planned under an output directory.

    Parameters
    ----------
    out_dir : str
        Output directory of the run.

    @return: paths of the batch directories holding a batchfile, in batch order.
    """

    if not os.path.isdir(out_dir):
        return []

    found = []
    for name in sorted(os.listdir(out_dir)):
        path = os.path.join(out_dir, name)
        if name.startswith(BATCH_DIR_PREFIX) and os.path.isdir(path):
            if os.path.exists(batchfile_path(path)):
                found.append(path)

    return found


def create_batches(rows: Sequence[Tuple[str, str]],
                   batch_size: int,
                   out_dir: str) -> List[str]:
    """Cut the release into batches and write the batchfile of each.

    Every batchfile is written before any batch is processed, so that the plan is
    complete the moment the first genome is worked on and a second machine
    starting later finds the same batches rather than making its own.

    Parameters
    ----------
    rows : sequence of tuple
        (FASTA path, accession) for every genome of the release, in the order the
        batches are to be cut from.
    batch_size : int
        Genomes per batch.
    out_dir : str
        Output directory of the run.

    @return: paths of the batch directories created, in batch order.
    """

    created = []
    for index, start in enumerate(range(0, len(rows), batch_size), start=1):
        batch_dir = os.path.join(out_dir, BATCH_DIR_FORMAT.format(index))
        os.makedirs(batch_dir, exist_ok=True)
        write_batchfile(rows[start:start + batch_size],
                        os.path.join(batch_dir, BATCHFILE_NAME), compress=True)
        created.append(batch_dir)

    return created


def canary_payload(**extra: object) -> str:
    """What a canary file says, beyond the fact that it exists.

    A claim has to name its owner for another machine to know whose it is, and
    for this machine to recognise a claim of its own left behind by a process
    that is no longer running.

    Parameters
    ----------
    extra : dict
        Further fields to record, one per line.

    @return: the contents of the canary file.
    """

    fields = [('host', socket.gethostname()),
              ('pid', os.getpid()),
              ('time', datetime.datetime.now().isoformat(timespec='seconds'))]
    fields += sorted(extra.items())

    return ''.join('{}\t{}\n'.format(key, value) for key, value in fields)


def read_canary(path: str) -> Dict[str, str]:
    """Read a canary file back.

    Parameters
    ----------
    path : str
        Canary file.

    @return: its fields, empty where the file has gone or cannot be read.
    """

    fields = {}
    try:
        with open(path) as handle:
            for line in handle:
                key, _, value = line.rstrip('\n').partition('\t')
                fields[key] = value
    except OSError:
        pass

    return fields


def process_alive(pid: str) -> bool:
    """Whether a process of this host is still running.

    Only ever asked about a PID this host recorded. A PID from another host says
    nothing here and is not looked up: PIDs are reused, and taking a batch from a
    machine that is still working on it costs more than leaving it.

    Parameters
    ----------
    pid : str
        Process ID, as the canary recorded it.

    @return: True if the process exists or cannot be ruled out, False if it is
             certainly gone.
    """

    try:
        os.kill(int(pid), 0)
    except ProcessLookupError:
        return False
    except (ValueError, TypeError):
        # not a PID at all, so nothing can be concluded from it
        return True
    except PermissionError:
        # running as another user, which means running
        return True

    return True


def server_time(directory: str) -> float:
    """What time it is by the clock of the file server holding a directory.

    A lease is only as good as the clock it is measured against, and the machines
    sharing an --out_dir have a clock each. What they do share is the file server,
    so a file is created in the directory and the time the server gives it is
    taken as now. Skew between the machines then cannot expire a live claim or
    hold a dead one.

    Parameters
    ----------
    directory : str
        Directory to ask about, which is the batch directory holding the claim.

    @return: the server's idea of now, as a POSIX timestamp; this machine's own
             clock where the directory cannot be written to.
    """

    handle, temp = None, None
    try:
        handle, temp = tempfile.mkstemp(prefix='.now.', dir=directory)
        return os.fstat(handle).st_mtime
    except OSError:
        return time.time()
    finally:
        if handle is not None:
            try:
                os.close(handle)
            except OSError:
                pass
        if temp is not None:
            try:
                os.unlink(temp)
            except OSError:
                pass


def claim_age(running_file: str) -> Optional[float]:
    """How long it is since the machine holding a batch last said so.

    Parameters
    ----------
    running_file : str
        The RUNNING canary of a batch.

    @return: seconds since the claim was last touched, or None if it has gone.
    """

    try:
        touched = os.stat(running_file).st_mtime
    except OSError:
        return None

    return max(0.0, server_time(os.path.dirname(running_file)) - touched)


class Heartbeat(object):
    """Touch a claim while its batch is worked on, so that it does not expire.

    The beating stops when the process holding the batch stops, whether it
    returns, is killed or wedges, which is the whole point: a claim outlives the
    process that made it by one lease and no longer.
    """

    def __init__(self, running_file: str, interval: float = HEARTBEAT_SECONDS) -> None:
        """Initialization.

        Parameters
        ----------
        running_file : str
            The RUNNING canary to keep alive.
        interval : float
            Seconds between touches.

        @return: None
        """

        self.running_file = running_file
        self.interval = interval
        self.stop = threading.Event()
        self.thread = None

    def beat(self) -> None:
        """Touch the claim until asked to stop.

        @return: None
        """

        # wait() returns True only when it was set, so the loop ends the moment
        # the batch does rather than after one more interval
        while not self.stop.wait(self.interval):
            try:
                os.utime(self.running_file, None)
            except OSError:
                # the claim has gone, which another machine taking the batch or
                # the batch finishing both look like; there is nothing to keep
                return

    def __enter__(self) -> 'Heartbeat':
        self.thread = threading.Thread(target=self.beat, daemon=True)
        self.thread.start()
        return self

    def __exit__(self, *exc_info: object) -> None:
        self.stop.set()
        if self.thread is not None:
            self.thread.join(timeout=self.interval)


@contextlib.contextmanager
def batch_log(batch_dir: str, logger: logging.Logger) -> Iterator[None]:
    """Write what happens to a batch into the batch's own directory as well.

    Parameters
    ----------
    batch_dir : str
        Batch directory, which takes trans_table.log.
    logger : logging.Logger
        The logger to tee, which is the 'timestamp' logger of the run.

    @return: a context in which the logger also writes to the batch.
    """

    handler = logging.FileHandler(os.path.join(batch_dir, BATCH_LOG_NAME), 'a')
    handler.setFormatter(logging.Formatter(
        fmt='[%(asctime)s] %(levelname)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'))
    logger.addHandler(handler)
    try:
        yield
    finally:
        logger.removeHandler(handler)
        handler.close()


def age_phrase(age: Optional[float]) -> str:
    """How long ago something was, as a log line says it.

    Parameters
    ----------
    age : float
        Seconds ago, or None where there is nothing to say.

    @return: a phrase naming the time, for a log message.
    """

    if age is None:
        return 'never'
    if age < 90:
        return '{:.0f}s ago'.format(age)
    if age < 5400:
        return '{:.0f}m ago'.format(age / 60)

    return '{:.1f}h ago'.format(age / 3600)


def batch_state(batch_dir: str) -> str:
    """What has happened to a batch.

    Parameters
    ----------
    batch_dir : str
        Batch directory.

    @return: one of the STATE_* constants.
    """

    if os.path.exists(os.path.join(batch_dir, SUCCESS_CANARY)):
        return STATE_SUCCESS
    if os.path.exists(os.path.join(batch_dir, RUNNING_CANARY)):
        return STATE_RUNNING
    if os.path.exists(os.path.join(batch_dir, FAILED_CANARY)):
        return STATE_FAILED

    return STATE_PENDING


def stale_claim(running_file: str,
                lease: float = CLAIM_LEASE_SECONDS) -> bool:
    """Whether a claim has been given up by the machine that made it.

    Two things say so. A claim of this host's whose process has gone is dead and
    known to be dead, which is what a reset leaves behind on the machine that
    reset. Any claim that has not been touched for a lease is dead as well: the
    machine holding a batch says so every HEARTBEAT_SECONDS for as long as it
    works, so silence for far longer than that is a machine that stopped, and
    whether it stopped by dying, by being killed or by wedging on a mount is
    neither knowable from here nor worth knowing.

    The second rule is what makes a batch recoverable from ANOTHER machine, and
    it is also the more reliable of the two: PIDs are reused, so after a reset
    the PID a claim names is as likely to belong to something new as to be
    missing, and a liveness check then holds a dead claim forever.

    Parameters
    ----------
    running_file : str
        The RUNNING canary of a batch.
    lease : float
        Seconds a claim survives without being touched.

    @return: True if the claim can be taken over without being asked to.
    """

    fields = read_canary(running_file)
    if (fields.get('host') == socket.gethostname()
            and not process_alive(fields.get('pid', ''))):
        return True

    age = claim_age(running_file)

    return age is not None and age > lease


def keep_failure_record(batch_dir: str) -> None:
    """Move a previous attempt's FAILED aside instead of deleting it.

    A batch is retried by claiming it, and the claim has to clear FAILED or the
    batch would still read as failed while it runs. Deleting it takes with it the
    only record of why the batch failed, which on a batch that fails the same way
    every time is the thing a person needs to read. It is kept under the time it
    was cleared, beside the batch it belongs to.

    Parameters
    ----------
    batch_dir : str
        Batch directory.

    @return: None
    """

    failed = os.path.join(batch_dir, FAILED_CANARY)
    kept = '{}.{}'.format(failed, datetime.datetime.now().strftime('%Y%m%dT%H%M%S'))
    try:
        os.rename(failed, kept)
    except OSError:
        pass


def release_claim(batch_dir: str) -> None:
    """Give up a claim without saying anything about how the batch went.

    What an interrupted run leaves: the batch was neither finished nor tried and
    found wanting, and the machine that held it is about to stop. Releasing it
    has the next run take it up rather than wait out the lease.

    Parameters
    ----------
    batch_dir : str
        Batch directory.

    @return: None
    """

    try:
        os.unlink(os.path.join(batch_dir, RUNNING_CANARY))
    except OSError:
        pass


def claim_batch(batch_dir: str,
                reclaim: bool = False,
                lease: float = CLAIM_LEASE_SECONDS) -> bool:
    """Take a batch for this machine, if no other machine holds it.

    The claim is made by linking a uniquely named file onto RUNNING rather than by
    creating RUNNING directly: this runs against NFS, where O_EXCL has never been
    the operation that two machines can race on safely and link() is. link()
    failing is read back rather than believed, for the same reason -- an NFS
    client can be told a link failed that in fact succeeded, and the link count of
    the file it made is what settles it.

    Parameters
    ----------
    batch_dir : str
        Batch directory to claim.
    reclaim : bool
        Take a batch another machine holds before its claim has expired. Only
        ever right when that machine is known not to be working on it.
    lease : float
        Seconds a claim survives without being touched.

    @return: True if this machine now holds the batch.
    """

    running_file = os.path.join(batch_dir, RUNNING_CANARY)

    if os.path.exists(running_file):
        if not (reclaim or stale_claim(running_file, lease)):
            return False
        try:
            os.unlink(running_file)
        except OSError:
            return False

    tmp_file = os.path.join(batch_dir, '.{}.{}.{}'.format(
        RUNNING_CANARY, os.getpid(), uuid.uuid4().hex))
    with open(tmp_file, 'w') as handle:
        handle.write(canary_payload())

    try:
        os.link(tmp_file, running_file)
        claimed = True
    except OSError:
        # the link may have been made even so, and the link count says whether
        claimed = os.stat(tmp_file).st_nlink == 2
    finally:
        try:
            os.unlink(tmp_file)
        except OSError:
            pass

    # a previous attempt on this batch is no longer what happened to it, though
    # what it had to say about itself is kept
    if claimed:
        keep_failure_record(batch_dir)

    return claimed


def finish_batch(batch_dir: str, **extra: object) -> None:
    """Record that a batch finished, and give up the claim on it.

    SUCCESS is written before RUNNING is removed. The other order leaves a moment
    in which the batch looks unclaimed and unfinished, which is the one state that
    would have a second machine repeat it.

    Parameters
    ----------
    batch_dir : str
        Batch directory.
    extra : dict
        Further fields to record in the canary.

    @return: None
    """

    with open(os.path.join(batch_dir, SUCCESS_CANARY), 'w') as handle:
        handle.write(canary_payload(**extra))

    try:
        os.unlink(os.path.join(batch_dir, RUNNING_CANARY))
    except OSError:
        pass


def mark_predicted(batch_dir: str, **extra: object) -> None:
    """Record that gTranslate has run over a batch and its results are final.

    Written the moment gTranslate returns 0, which is hours of work, and read by
    whoever next takes the batch, which after a machine is lost between the
    prediction and the comparison is another machine a minute later. The
    comparison that follows is seconds and is simply done again.

    Parameters
    ----------
    batch_dir : str
        Batch directory.
    extra : dict
        Further fields to record in the canary.

    @return: None
    """

    with open(os.path.join(batch_dir, PREDICTED_CANARY), 'w') as handle:
        handle.write(canary_payload(**extra))


def already_predicted(batch_dir: str, summary: str) -> bool:
    """Whether gTranslate has already run over a batch.

    The canary alone is not enough: it is taken to mean the predictions are there
    only alongside the summary gTranslate wrote, since a canary without the
    results it speaks for would have the comparison read a file that is not
    there.

    Parameters
    ----------
    batch_dir : str
        Batch directory.
    summary : str
        Name of gTranslate's summary within the batch directory.

    @return: True if the predictions are there to be compared.
    """

    return (os.path.exists(os.path.join(batch_dir, PREDICTED_CANARY))
            and os.path.exists(os.path.join(batch_dir, summary)))


def report_no_prediction(batch_dir: str,
                         given: Sequence[str],
                         predicted: Sequence[str]) -> List[str]:
    """Name the genomes gTranslate was given and said nothing about.

    With --force gTranslate drops a genome it cannot process and carries on,
    which is what keeps one bad genome from costing a batch of ten thousand. A
    genome dropped that way is simply absent from the summary, so without this it
    is discovered by prodigal refusing to call a release for want of a table.

    Parameters
    ----------
    batch_dir : str
        Batch directory.
    given : sequence of str
        Accessions handed to gTranslate.
    predicted : sequence of str
        Accessions the summary holds a prediction for.

    @return: the accessions with no prediction, in the order they were given.
    """

    missing = [accession for accession in given if accession not in set(predicted)]
    if not missing:
        return missing

    with open(os.path.join(batch_dir, NO_PREDICTION_NAME), 'w') as handle:
        for accession in missing:
            handle.write('{}\n'.format(accession))

    return missing


def fail_batch(batch_dir: str, reason: str) -> None:
    """Record that a batch was attempted and did not finish.

    The claim is given up, so the batch is retried by the next run without
    --reclaim: this machine is known not to be working on it, which is exactly
    what --reclaim exists to assert about a machine that cannot be asked.

    Parameters
    ----------
    batch_dir : str
        Batch directory.
    reason : str
        What went wrong, recorded for whoever reads the directory.

    @return: None
    """

    with open(os.path.join(batch_dir, FAILED_CANARY), 'w') as handle:
        handle.write(canary_payload(reason=' '.join(str(reason).split())))

    try:
        os.unlink(os.path.join(batch_dir, RUNNING_CANARY))
    except OSError:
        pass


def detect_table_command(batchfile: str,
                         out_dir: str,
                         cpus: int = 1,
                         tmp_dir: Optional[str] = None,
                         force: bool = False,
                         keep_called_genes: bool = False,
                         prefix: Optional[str] = None,
                         custom_model_path: Optional[str] = None) -> List[str]:
    """Build the gtranslate detect_table command line.

    An option the run was not given is left off the command line rather than
    passed with a default of this module's choosing, so gTranslate's own defaults
    remain the defaults and do not have to be tracked here as it changes.

    Parameters
    ----------
    batchfile : str
        Batchfile naming the genomes to process.
    out_dir : str
        Directory gTranslate writes its results to, which is the batch directory.
    cpus : int
        Number of genomes processed at once.
    tmp_dir : str
        Directory for gTranslate's intermediate files, or None for its default.
    force : bool
        Carry on when a single genome fails rather than stopping the batch.
    keep_called_genes : bool
        Keep the genes called under the predicted table.
    prefix : str
        Prefix of gTranslate's output files, or None for its default.
    custom_model_path : str
        Classifiers to predict with, or None to use GTRANSLATE_MODEL_PATH.

    @return: the command as a list of arguments, ready for subprocess.
    """

    cmd = [GTRANSLATE_BIN, DETECT_TABLE,
           '--batchfile', batchfile,
           '--out_dir', out_dir,
           '--cpus', str(cpus)]

    if tmp_dir:
        cmd += ['--tmpdir', tmp_dir]
    if prefix:
        cmd += ['--prefix', prefix]
    if custom_model_path:
        cmd += ['--custom_model_path', custom_model_path]
    if force:
        cmd += ['--force']
    if keep_called_genes:
        cmd += ['--keep_called_genes']

    return cmd


def read_taxonomy(taxonomy_file: str) -> Dict[str, str]:
    """Read the standardised NCBI taxonomy of the genomes.

    The file is the two-column accession / semicolon-separated lineage TSV
    ncbi_metadata_sync writes. Each genome is recorded under the accession as
    given AND under its canonical form, so that a GenBank genome of a release
    finds the lineage the taxonomy holds against its RefSeq counterpart; an
    accession given exactly is preferred to a canonical match.

    Parameters
    ----------
    taxonomy_file : str
        Standardised NCBI taxonomy file.

    @return: accession, and canonical accession, to lineage.
    """

    exact, canonical = {}, {}
    with open(taxonomy_file) as handle:
        for line in tqdm(handle, ncols=100, leave=False, desc='Reading taxonomy'):
            line = line.rstrip('\n')
            if not line:
                continue
            accession, _, lineage = line.partition('\t')
            if not lineage:
                continue
            exact[accession] = lineage
            canonical.setdefault(canonical_gid(accession), lineage)

    canonical.update(exact)

    return canonical


def lineage_of(accession: str, taxonomy: Dict[str, str]) -> str:
    """The NCBI lineage of a genome.

    Parameters
    ----------
    accession : str
        Accession of the genome.
    taxonomy : dict
        Taxonomy as read_taxonomy() returned it.

    @return: the lineage, or NCBI's null where the taxonomy does not hold one.
    """

    if accession in taxonomy:
        return taxonomy[accession]

    return taxonomy.get(canonical_gid(accession), NCBI_NA)


def checkm_conflict(gtranslate_table: str, checkm_table: Optional[int]) -> bool:
    """Whether the density rule would have called this genome's genes recoded or
    not recoded the other way from gTranslate.

    What matters about the two tables is not that they are different numbers but
    that they are different KINDS of answer: 11 is the standard code and 4 and 25
    are recodings of it, and a genome called under the wrong kind has its genes
    truncated or run together at every TGA. So gTranslate saying 11 where the
    density rule says 4, and gTranslate saying 4 or 25 where the rule says 11,
    are the conflicts.

    25 against 4 is NOT one of them. The density rule calls genes under tables 4
    and 11 and picks between those two alone, so 4 is the only recoding it can
    ever return; a genome gTranslate calls 25 and the rule calls 4 is one the two
    agree about as far as the rule is able to say, which is what makes table 25
    the thing gTranslate is for.

    Parameters
    ----------
    gtranslate_table : str
        Table gTranslate predicted, as the summary gives it.
    checkm_table : int
        Table the density rule chose, or None where it could not be worked out.

    @return: True where the two disagree about whether the genome is recoded.
    """

    try:
        predicted = int(gtranslate_table)
    except (TypeError, ValueError):
        return False

    if checkm_table is None:
        return False

    if predicted == STANDARD_TABLE:
        return checkm_table in RECODED_TABLES

    if predicted in RECODED_TABLES:
        return checkm_table == STANDARD_TABLE

    return False


class BadConflictFile(ValueError):
    """A conflict file without the columns a conflict file has.

    Raised rather than skipped: the file is written by this command and read by
    it, so one that does not look like one is a bug or a truncated write, and
    carrying on would annotate a release with the wrong genomes' quality.
    """


class ComparisonCounts(NamedTuple):
    """What comparing a batch against NCBI found.

    The three are kept together because the release line needs all of them and
    gets them from the batches' canaries, which is the only place a batch run on
    another machine leaves them.
    """

    compared: int
    conflicts: int
    no_ncbi_table: int


def disagreement_rate(conflicts: int, compared: int) -> float:
    """What share of the genomes NCBI declares a table for gTranslate disagrees
    with.

    The denominator is the genomes that COULD be compared and not the genomes of
    the release: a genome NCBI has not annotated is not a genome the two agree or
    disagree about, and counting it would make the rate a measure of how much of
    the release NCBI has annotated rather than of how often the two differ.

    Parameters
    ----------
    conflicts : int
        Genomes the two disagree about.
    compared : int
        Genomes NCBI declares a table for, which are the ones compared.

    @return: the percentage, and 0.0 where nothing could be compared.
    """

    if compared <= 0:
        return 0.0

    return 100.0 * conflicts / compared


def batch_counts(batches: Sequence[str]) -> Tuple[int, Optional[int]]:
    """Add up what every batch of the release recorded about its comparison.

    The counts are read from the batches' SUCCESS canaries rather than
    recomputed, because a release is predicted by several machines and a machine
    running the last batch has compared none of the others; the canary is what
    every machine leaves behind.

    Only `compared` is certain to be there. The field saying how many genomes
    NCBI declares no table for was added after r237 had been predicted, so its
    batches record the one and not the other, and a release finished before the
    change can still report the number asked for and the rate. Where any batch is
    missing it, None is returned and the release line leaves that clause out
    rather than reporting a total that is short by whatever those batches found.

    Parameters
    ----------
    batches : sequence of str
        Every batch directory of the run.

    @return: (compared, no_ncbi_table), the second None where any batch's canary
             does not record it.
    """

    compared, no_ncbi_table = 0, 0
    complete = True
    for batch_dir in batches:
        canary = read_canary(os.path.join(batch_dir, SUCCESS_CANARY))
        try:
            compared += int(canary['compared'])
        except (KeyError, ValueError):
            complete = False
            continue

        try:
            no_ncbi_table += int(canary['no_ncbi_table'])
        except (KeyError, ValueError):
            complete = False

    return compared, no_ncbi_table if complete else None


def comparison_rows(predictions: Dict[str, Dict[str, str]],
                    genome_dirs: Dict[str, str],
                    taxonomy: Dict[str, str]) -> Tuple[List[List[str]], int, int]:
    """Compare what gTranslate predicted against what NCBI declares, for every
    genome the two both called.

    Only a genome NCBI declares a table for can be compared at all: a genome NCBI
    has not annotated has nothing to compare against and no ncbi_conflict to
    report, so it is counted rather than given a row saying nothing.

    Every genome that CAN be compared gets a row, agreements included, which is
    what separates this from the conflicts: the rate the two differ at, and
    whether the genomes they differ about are unlike the ones they agree about,
    are questions the agreements have to be present to answer. The conflicts are
    then this file filtered, by conflicts_from_comparison(), so that the two can
    never come to disagree about which genomes conflicted.

    The genome statistics are gTranslate's own measurements of the FASTA it read,
    carried through from the summary. Measuring them here would be reading 1.3M
    genomes a second time to learn what has already been learnt.

    Parameters
    ----------
    predictions : dict
        Predictions as read_translation_table_summary() returned them.
    genome_dirs : dict
        Accession to genome directory, for the genomes of the batch.
    taxonomy : dict
        Taxonomy as read_taxonomy() returned it.

    @return: (rows, compared, no_ncbi_table), the comparison in accession order
             and in COMPARISON_HEADER order, the number of genomes compared, and
             the number of genomes NCBI declared no table for.
    """

    rows, compared, no_ncbi_table = [], 0, 0
    for accession in sorted(predictions):
        genome_dir = genome_dirs.get(accession)
        if genome_dir is None:
            continue

        ncbi_table = ncbi_translation_table(genomic_gff(genome_dir))
        if ncbi_table is None:
            no_ncbi_table += 1
            continue

        compared += 1

        predicted = predictions[accession]
        table = predicted.get(TT_SUMMARY_TABLE, '')
        density_4 = predicted.get(TT_SUMMARY_DENSITY_4, '')
        density_11 = predicted.get(TT_SUMMARY_DENSITY_11, '')
        checkm_table = checkm_translation_table(density_4, density_11)

        rows.append([accession,
                     table,
                     str(ncbi_table),
                     str(checkm_table) if checkm_table else NCBI_NA,
                     str(table.strip() != str(ncbi_table)),
                     str(checkm_conflict(table.strip(), checkm_table)),
                     density_4,
                     density_11,
                     predicted.get(TT_SUMMARY_GC) or NCBI_NA,
                     predicted.get(TT_SUMMARY_N50) or NCBI_NA,
                     predicted.get(TT_SUMMARY_GENOME_SIZE) or NCBI_NA,
                     lineage_of(accession, taxonomy)])

    return rows, compared, no_ncbi_table


def conflicts_from_comparison(rows: Sequence[Sequence[str]]) -> List[List[str]]:
    """The genomes of a comparison that gTranslate and NCBI disagree about.

    The conflicts are the comparison filtered rather than a second walk over the
    genomes: the GFF of every genome has already been read once to make the
    comparison, and reading it again would be hours of a release spent finding
    out what is already known.

    ncbi_conflict does not come with them. Every row of the conflict file is a
    conflict, so a column saying so would say 'True' and nothing else; the column
    exists in the comparison because there it distinguishes the rows.

    Parameters
    ----------
    rows : sequence of sequence of str
        Comparison rows, in COMPARISON_HEADER order.

    @return: the conflicting rows, in CONFLICT_HEADER order.
    """

    verdict = COMPARISON_HEADER.index('ncbi_conflict')
    keep = [COMPARISON_HEADER.index(column) for column in CONFLICT_HEADER]

    return [[row[column] for column in keep]
            for row in rows if row[verdict] == 'True']


def write_table(rows: Sequence[Sequence[str]], path: str,
                header: Sequence[str] = CONFLICT_HEADER,
                compress: bool = False) -> None:
    """Write a comparison or a conflict table, of a batch or of the release.

    The file is written whether or not there are any rows: a batch that finished
    with nothing to report says so with a header and no rows, and the release
    file is then the concatenation of every batch's, however many they found.

    Parameters
    ----------
    rows : sequence of sequence of str
        Rows, as comparison_rows(), conflicts_from_comparison() or
        annotate_conflicts() returned them.
    path : str
        File to write.
    header : sequence of str
        Column names: COMPARISON_HEADER, CONFLICT_HEADER for a batch, or
        CONFLICT_HEADER_CHECKM2 for the release, which carries the CheckM2
        columns as well.
    compress : bool
        Write it gzipped, which the comparison is and the conflicts are not.

    @return: None
    """

    with (gzip.open(path, 'wt') if compress else open(path, 'w')) as handle:
        handle.write('\t'.join(header) + '\n')
        for row in rows:
            handle.write('\t'.join(row) + '\n')


def read_conflicts(path: str) -> List[List[str]]:
    """Read a conflict file back as the columns a batch writes.

    The columns are taken by name, so a release file that has ALREADY been
    annotated is read back as the unannotated row it was made from and can be
    annotated again. Without that the CheckM2 columns would be appended to a row
    that already had them every time the command was run over a finished output
    directory, which is what running it again is for.

    Parameters
    ----------
    path : str
        Conflict file, of a batch or of the release.

    @return: its rows in CONFLICT_HEADER order, each a list of fields, and an
             empty list where the file holds nothing but a header.
    """

    rows = []
    with open_text(path) as handle:
        header = handle.readline().rstrip('\n').split('\t')
        try:
            columns = [header.index(column) for column in CONFLICT_HEADER]
        except ValueError:
            raise BadConflictFile(
                '{} does not have the columns of a conflict file.'.format(path))

        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) <= columns[-1]:
                continue
            rows.append([fields[column] for column in columns])

    return rows


def conflict_tables(rows: Sequence[Sequence[str]]) -> Dict[int, List[str]]:
    """Which genomes have to have their quality estimated under which table.

    A genome is asked about under BOTH of the tables its row disputes, so every
    accession appears under two of them. Grouping by table rather than by genome
    is what makes this a handful of CheckM2 runs instead of hundreds: CheckM2
    loads its models and searches the whole DIAMOND database once per run,
    whatever the run holds, so the cost is in the number of runs and barely in
    the number of genomes.

    Parameters
    ----------
    rows : sequence of sequence of str
        Conflicting rows, in CONFLICT_HEADER order.

    @return: table to the accessions to be run under it, in accession order.
    """

    columns = (CONFLICT_HEADER.index('gtranslate_tt'),
               CONFLICT_HEADER.index('ncbi_tt'))

    tables = {}
    for row in rows:
        for column in columns:
            try:
                table = int(row[column])
            except (IndexError, TypeError, ValueError):
                # a row whose table is not a number is one nothing can be run
                # under; it keeps its other table and gets NCBI_NA for this one
                continue
            tables.setdefault(table, set()).add(row[0])

    return {table: sorted(accessions) for table, accessions in tables.items()}


def release_fastas(batches: Sequence[str],
                   accessions: Sequence[str]) -> Dict[str, str]:
    """Where the genomic FASTA of each named genome is, across every batch.

    The batchfiles are read rather than the genome_dirs file, for the reason
    compare_batch() reads them: they say what the release was run on, and a
    genome_dirs file says what it says now. Only the genomes asked about are
    kept, which is a few hundred of a million-odd lines.

    Parameters
    ----------
    batches : sequence of str
        Every batch directory of the run.
    accessions : sequence of str
        Genomes wanted.

    @return: accession to genomic FASTA, omitting any the batchfiles do not name.
    """

    wanted = set(accessions)

    fastas = {}
    for batch_dir in batches:
        batchfile = batchfile_path(batch_dir)
        if not os.path.exists(batchfile):
            continue
        for fasta, accession in read_batchfile(batchfile):
            if accession in wanted:
                fastas[accession] = fasta

    return fastas


def stage_checkm2_input(accessions: Sequence[str],
                        fastas: Dict[str, str],
                        staging: str) -> List[str]:
    """Name each genome's FASTA for its accession, so CheckM2 gives it back.

    CheckM2 names a result for the file it read, and NCBI names the file for the
    assembly and not for the accession: GCA_000238995.1_ASM23899v1_genomic.fna.gz
    comes back as GCA_000238995.1_ASM23899v1_genomic, which no conflict row is
    keyed by. A symlink named <accession>.fna.gz comes back as the accession,
    which is what joins the report to the rows. Symlinks rather than copies
    because the genomes are the release and are read, not written.

    Parameters
    ----------
    accessions : sequence of str
        Genomes to stage.
    fastas : dict
        Accession to the genomic FASTA of that genome.
    staging : str
        Directory the links are made in, created if it is not there.

    @return: the staged paths, which omit any genome whose FASTA is missing.
    """

    os.makedirs(staging, exist_ok=True)

    staged = []
    for accession in accessions:
        fasta = fastas.get(accession)
        if not fasta or not os.path.exists(fasta):
            continue

        link = os.path.join(staging, accession + CHECKM2_LINK_EXT)
        if os.path.islink(link) or os.path.exists(link):
            os.unlink(link)
        os.symlink(fasta, link)
        staged.append(link)

    return staged


def checkm2_command(fastas: Sequence[str],
                    table: int,
                    out_dir: str,
                    threads: int) -> List[str]:
    """The CheckM2 command run over the genomes disputing one table.

    --ttable is the whole point of the run: CheckM2 left to itself calls genes
    under whichever of tables 4 and 11 gives the better coding density, which is
    the rule the checkm_tt column already reports and is not what is being asked.
    Forcing the table asks what the genome looks like if that table is right.

    Parameters
    ----------
    fastas : sequence of str
        Staged genomic FASTA files, named for their accessions.
    table : int
        Translation table genes are called under.
    out_dir : str
        Directory CheckM2 writes its report and intermediates to.
    threads : int
        Threads CheckM2 is given.

    @return: the command as a list of arguments, ready for subprocess.
    """

    # --input takes the rest of the command line, so it goes last
    return [CHECKM2_BIN, 'predict',
            '--ttable', str(table),
            '--threads', str(threads),
            '--force',
            '--output-directory', out_dir,
            '--input'] + list(fastas)


def read_checkm2_report(path: str) -> Dict[str, Tuple[str, str]]:
    """Read the completeness and contamination CheckM2 estimated.

    Read by column name rather than by position, as the NCBI tables are: CheckM2
    has added columns between releases and puts the ones this wants in the middle
    of a dozen it does not.

    Parameters
    ----------
    path : str
        quality_report.tsv of one CheckM2 run.

    @return: accession to (completeness, contamination), empty where the file is
             absent or has no header to find the columns by.
    """

    quality = {}
    try:
        with open(path) as handle:
            header = handle.readline().rstrip('\n').split('\t')
            try:
                name = header.index(CHECKM2_NAME)
                completeness = header.index(CHECKM2_COMPLETENESS)
                contamination = header.index(CHECKM2_CONTAMINATION)
            except ValueError:
                return {}

            for line in handle:
                if not line.strip():
                    continue
                fields = line.rstrip('\n').split('\t')
                if len(fields) <= max(name, completeness, contamination):
                    continue
                quality[fields[name]] = (fields[completeness],
                                         fields[contamination])
    except OSError:
        return {}

    return quality


def quality_score(completeness: float, contamination: float) -> float:
    """The GTDB quality score of a genome.

    Contamination is charged at five times the weight of completeness because the
    two are not equally recoverable: a genome missing a marker is missing it, and
    a genome carrying another organism's markers reports things about that
    organism as though they were its own.

    Parameters
    ----------
    completeness : float
        CheckM2 completeness, as a percentage.
    contamination : float
        CheckM2 contamination, as a percentage.

    @return: the score, which may be negative.
    """

    return completeness - QC_CONTAMINATION_WEIGHT * contamination


def passes_qc(completeness: str, contamination: str) -> Optional[bool]:
    """Whether a genome passes standard GTDB QC.

    All three conditions have to hold. They are not the same condition said three
    ways: the score alone would keep a genome 96% complete and 9% contaminated,
    and the completeness alone would keep anything that had been sequenced.

    A genome CheckM2 returned nothing for does not fail -- nothing is known about
    it, and na is not False. Saying otherwise would put a genome whose FASTA is
    missing among the genomes that were looked at and found wanting.

    Parameters
    ----------
    completeness : str
        CheckM2 completeness as the report gives it, or NCBI_NA.
    contamination : str
        CheckM2 contamination as the report gives it, or NCBI_NA.

    @return: True where it passes, False where it fails, None where there is no
             estimate to judge it by.
    """

    try:
        estimated_completeness = float(completeness)
        estimated_contamination = float(contamination)
    except (TypeError, ValueError):
        return None

    return (estimated_completeness > QC_MIN_COMPLETENESS
            and estimated_contamination < QC_MAX_CONTAMINATION
            and quality_score(estimated_completeness,
                              estimated_contamination) > QC_MIN_QUALITY)


def remove_checkm2_dir(checkm2_dir: str) -> None:
    """Remove what the CheckM2 runs left behind.

    Called once the estimates are in ncbi_tt_conflict.tsv, which is where they
    were wanted: what is left is the staged links, the called proteins and the
    DIAMOND output, about 600 KB per genome per table, and for a release that is
    a few hundred megabytes of a metadata directory holding nothing the conflict
    file does not now hold.

    What is given up by removing it is that a later run reads the reports rather
    than making them again. That is a few minutes for the few hundred genomes a
    release conflicts about -- CheckM2 spends most of a small run loading its
    models and opening the database -- against the days of prediction the batches
    protect, which is why they are not removed and this is.

    Errors are ignored. The directory is this command's own working space, and a
    release whose conflicts are written and annotated is finished whether or not
    the leftovers could be swept up.

    Parameters
    ----------
    checkm2_dir : str
        The CheckM2 directory of the run.

    @return: None
    """

    shutil.rmtree(checkm2_dir, ignore_errors=True)


def annotate_conflicts(rows: Sequence[Sequence[str]],
                       quality: Dict[int, Dict[str, Tuple[str, str]]]
                       ) -> List[List[str]]:
    """Put each genome's two CheckM2 estimates into its conflict row.

    A genome CheckM2 returned nothing for -- one whose FASTA is missing, one a
    run failed on, one CheckM2 itself dropped -- keeps its row and gets NCBI_NA
    for what is not known. The row is the conflict, and the conflict is there
    whether or not its quality could be estimated.

    Parameters
    ----------
    rows : sequence of sequence of str
        Conflicting rows, in CONFLICT_HEADER order.
    quality : dict
        Table to the report of the run made under it, as read_checkm2_report()
        returned each.

    @return: the rows in CONFLICT_HEADER_CHECKM2 order.
    """

    gtranslate = CONFLICT_HEADER.index('gtranslate_tt')
    ncbi = CONFLICT_HEADER.index('ncbi_tt')

    annotated = []
    for row in rows:
        estimates = []
        for column in (gtranslate, ncbi):
            try:
                table = int(row[column])
            except (IndexError, TypeError, ValueError):
                table = None

            completeness, contamination = quality.get(table, {}).get(
                row[0], (NCBI_NA, NCBI_NA))
            verdict = passes_qc(completeness, contamination)
            estimates.extend((completeness, contamination,
                              NCBI_NA if verdict is None else str(verdict)))

        annotated.append(list(row[:-1]) + estimates + [row[-1]])

    return annotated


def concatenate(files: Sequence[str], path: str, compress: bool = False) -> int:
    """Join the tables of every batch into one, keeping a single header.

    Parameters
    ----------
    files : sequence of str
        Files to join, each with the same header, in batch order.
    path : str
        File to write.
    compress : bool
        Write it gzipped, which the release summary is and the conflicts are not:
        the summary is a row per genome of the release and the conflicts are a few
        hundred rows meant to be looked at.

    @return: number of rows written, the header not counted.
    """

    written = 0
    with (gzip.open(path, 'wt') if compress else open(path, 'w')) as out:
        for index, name in enumerate(files):
            # open_text: a batch's comparison is gzipped and its conflicts are
            # not, and gTranslate's summary is whatever gTranslate wrote
            with open_text(name) as handle:
                header = handle.readline()
                if index == 0:
                    out.write(header)
                for line in handle:
                    if line.strip():
                        out.write(line)
                        written += 1

    return written


class GTranslate(object):
    """Predict the translation table of each genome of a release, in batches."""

    def __init__(self,
                 cpus: int = 1,
                 batch_size: int = DEFAULT_BATCH_SIZE,
                 tmp_dir: Optional[str] = None,
                 force: bool = True,
                 keep_called_genes: bool = False,
                 prefix: Optional[str] = None,
                 custom_model_path: Optional[str] = None,
                 reclaim: bool = False,
                 lease: float = CLAIM_LEASE_SECONDS,
                 heartbeat: float = HEARTBEAT_SECONDS) -> None:
        """Initialization.

        Parameters
        ----------
        cpus : int
            Number of genomes gTranslate processes at once.
        batch_size : int
            Genomes per batch.
        tmp_dir : str
            Directory for gTranslate's intermediate files, or None for its default.
        force : bool
            Carry on when a single genome fails rather than stopping the batch.
        keep_called_genes : bool
            Keep the genes called under the predicted table.
        prefix : str
            Prefix of gTranslate's output files, or None for its default.
        custom_model_path : str
            Classifiers to predict with, or None to use GTRANSLATE_MODEL_PATH.
        reclaim : bool
            Take over a batch another machine holds before its claim has expired.
        lease : float
            Seconds a claim survives without the machine holding it saying so.
        heartbeat : float
            Seconds between this machine saying so about a batch of its own.

        @return: None
        """

        self.cpus = cpus
        self.batch_size = batch_size
        self.tmp_dir = tmp_dir
        self.force = force
        self.keep_called_genes = keep_called_genes
        self.prefix = prefix
        self.custom_model_path = custom_model_path
        self.reclaim = reclaim
        self.lease = lease
        self.heartbeat = heartbeat

        check_dependencies([GTRANSLATE_BIN, 'checkm2', 'prodigal'])

        self.logger = logging.getLogger('timestamp')

    def plan_batches(self, gtdb_genome_path_file: str, out_dir: str) -> List[str]:
        """Settle which genomes are in which batch, once for every machine.

        A plan already under the output directory is used as it stands. It is what
        another machine is working from and what the finished batches were cut
        from, and partitioning a release again that has since gained or lost a
        genome would move genomes between batches that are already done.

        Parameters
        ----------
        gtdb_genome_path_file : str
            genome_dirs file of the release.
        out_dir : str
            Output directory of the run.

        @return: paths of the batch directories, in batch order.
        """

        existing = batch_dir_names(out_dir)
        if existing:
            self.logger.warning(
                'warning: {:,} batch(es) are already planned under {}; using them '
                'and not regenerating the batchfiles. Remove the batch directories '
                'to partition the release again.'.format(len(existing), out_dir))
            return existing

        genomes = read_genome_dirs(gtdb_genome_path_file)
        self.logger.info('Read {:,} genomes from {}.'.format(
            len(genomes), gtdb_genome_path_file))

        # sorted so that which genomes are in batch N follows from the set of
        # genomes, and not from the order the genome_dirs file was written in
        genomes.sort(key=lambda genome: genome[0])

        # named, not stat-ed: whether the file is there is asked of each batch as
        # it is run, where ten thousand stat calls are nothing against the hours
        # gTranslate then spends, rather than of the whole release here, where on
        # r237 it is hours over NFS before a single batch directory exists and a
        # run stopped in it has no plan to resume from. It also leaves which
        # genomes share a batch following from the genome_dirs file alone.
        rows = [(genomic_fasta(genome_dir), accession)
                for accession, genome_dir in genomes]
        if not rows:
            raise RuntimeError(
                '{} names no genomes.'.format(gtdb_genome_path_file))

        batches = create_batches(rows, self.batch_size, out_dir)
        self.logger.info('Planned {:,} genomes as {:,} batch(es) of up to {:,}.'.format(
            len(rows), len(batches), self.batch_size))

        return batches

    def run_gtranslate(self, batch_dir: str) -> None:
        """Run gTranslate over the genomes of one batch.

        Parameters
        ----------
        batch_dir : str
            Batch directory, which holds the batchfile and takes the results.

        @return: None
        """

        # gTranslate checks the paths of a batchfile before it starts and refuses
        # the whole batch if one of them is missing, so a genome whose FASTA is
        # not there costs the other ten thousand -- and costs them again on every
        # retry, the batch failing identically each time. It is left out here and
        # named in the batch directory instead.
        batchfile, present, missing = check_batch_fastas(batch_dir)
        if missing:
            self.logger.warning(
                'warning: {:,} genome(s) of {} have no genomic FASTA and were left '
                'out; the first is {}. They are named in {}.'.format(
                    len(missing), os.path.basename(batch_dir), missing[0],
                    MISSING_NAME))
        if not present:
            raise RuntimeError(
                'None of the {:,} genomes of {} has a genomic FASTA to process.'.format(
                    len(missing), batch_dir))

        cmd = detect_table_command(batchfile,
                                   batch_dir,
                                   cpus=self.cpus,
                                   tmp_dir=self.tmp_dir,
                                   force=self.force,
                                   keep_called_genes=self.keep_called_genes,
                                   prefix=self.prefix,
                                   custom_model_path=self.custom_model_path)
        self.logger.info('Command: {}'.format(' '.join(cmd)))

        # gTranslate's output is left to the terminal rather than read back a
        # line at a time and logged again. Its progress bars redraw one line with
        # a carriage return, and a carriage return is a line ending to a reader,
        # so re-logging what it prints turns each bar into one line per redraw and
        # the bar into a wall. Nothing is lost by not capturing it: gTranslate
        # writes its own gtranslate.log into the batch directory, and writes it
        # without the bars, which is the better record of a batch anyway.
        # --silent is honoured by discarding the output rather than by showing it.
        silent = getattr(self.logger, 'is_silent', False)
        proc = subprocess.run(cmd,
                              stdout=subprocess.DEVNULL if silent else None,
                              stderr=subprocess.STDOUT if silent else None)

        if proc.returncode != 0:
            raise RuntimeError('{} returned exit code {}; {} says what it was '
                               'doing.'.format(GTRANSLATE_BIN, proc.returncode,
                                               os.path.join(batch_dir, 'gtranslate.log')))

        # the hours of the batch are over and what they produced is final, so a
        # machine that takes this batch after here compares rather than predicts
        mark_predicted(batch_dir, genomes=len(present))

        # gTranslate has read it and will not be run over this batch again, so
        # the uncompressed copy goes and the batch keeps the compressed plan
        # alone. A batch taken over before this point writes it again from that
        # plan, which is why removing it here costs nothing.
        try:
            os.unlink(os.path.join(batch_dir, PRESENT_BATCHFILE_NAME))
        except OSError:
            pass

        # --force has gTranslate drop a genome it cannot process rather than end
        # the batch, and a genome dropped that way is simply absent from the
        # summary; it is named here rather than found by prodigal later
        predicted = read_translation_table_summary(
            os.path.join(batch_dir, summary_name(self.prefix)))
        no_prediction = report_no_prediction(
            batch_dir, [accession for _, accession in present], list(predicted))
        if no_prediction:
            self.logger.warning(
                'warning: gTranslate returned no prediction for {:,} genome(s) of '
                '{}; the first is {}. They are named in {}.'.format(
                    len(no_prediction), os.path.basename(batch_dir),
                    no_prediction[0], NO_PREDICTION_NAME))

    def compare_batch(self,
                      batch_dir: str,
                      taxonomy: Dict[str, str]) -> ComparisonCounts:
        """Compare a batch's predictions against the tables NCBI declares.

        The genome directories are taken from the batch's own batchfile, so the
        comparison asks about the genomes the batch was run on rather than about
        whatever a genome_dirs file says now.

        Parameters
        ----------
        batch_dir : str
            Batch directory.
        taxonomy : dict
            Taxonomy as read_taxonomy() returned it.

        @return: what the comparison found, which run() records in the batch's
                 canary for the release line to add up.
        """

        genome_dirs = {accession: os.path.dirname(fasta) for fasta, accession
                       in read_batchfile(batchfile_path(batch_dir))}

        predictions = read_translation_table_summary(
            os.path.join(batch_dir, summary_name(self.prefix)))

        rows, compared, no_ncbi_table = comparison_rows(
            predictions, genome_dirs, taxonomy)
        write_table(rows, os.path.join(batch_dir, COMPARISON_NAME),
                    header=COMPARISON_HEADER, compress=True)

        # the conflicts are the comparison filtered, so the two cannot come to
        # disagree about which genomes conflicted
        conflicts = conflicts_from_comparison(rows)
        write_table(conflicts, os.path.join(batch_dir, CONFLICT_NAME))

        # the genomes NCBI declares no table for are in neither file, having
        # nothing to compare, so this line is the only place they are counted
        self.logger.info(
            'Compared {:,} genomes with a translation table from NCBI: {:,} agree, '
            '{:,} conflict ({:.2f}%); {:,} genome(s) have no table from NCBI to '
            'compare.'.format(compared, compared - len(conflicts), len(conflicts),
                              disagreement_rate(len(conflicts), compared),
                              no_ncbi_table))

        return ComparisonCounts(compared=compared,
                                conflicts=len(conflicts),
                                no_ncbi_table=no_ncbi_table)

    def report_comparison(self, batches: Sequence[str], conflicts: int) -> None:
        """Say for the whole release how much of it NCBI declares a table for and
        how often gTranslate disagrees.

        Every batch has said this about itself, in its own log, on whichever of
        the machines took it; this is the line that says it about the release,
        and it is written by whichever machine finishes last.

        Parameters
        ----------
        batches : sequence of str
            Every batch directory of the run.
        conflicts : int
            Rows of the release conflict file, which are the disagreements.

        @return: None
        """

        compared, no_ncbi_table = batch_counts(batches)

        line = ('Conflicts: {:,} genome(s) have a translation table from NCBI, and '
                'gTranslate disagrees with NCBI for {:,} ({:.2f}%) genomes.'.format(
                    compared, conflicts, disagreement_rate(conflicts, compared)))

        if no_ncbi_table is not None:
            line += ' {:,} genome(s) have no table from NCBI to compare.'.format(
                no_ncbi_table)

        self.logger.info(line)

    def aggregate(self, batches: Sequence[str], out_dir: str) -> None:
        """Write the conflicts and the summary for the whole release, report the
        comparison, and estimate the quality of the genomes it conflicted about.

        Written only once every batch has succeeded, so that the files at the top
        of the output directory are either the whole release or absent, and never
        a part of it that reads like the whole. This is also where the work that
        is about the release rather than about a batch belongs: the rate NCBI and
        gTranslate disagree at, and the CheckM2 runs, which are made once over the
        genomes every batch together found.

        Parameters
        ----------
        batches : sequence of str
            Every batch directory of the run, in batch order.
        out_dir : str
            Output directory of the run.

        @return: None
        """

        unfinished = [batch for batch in batches
                      if batch_state(batch) != STATE_SUCCESS]
        if unfinished:
            self.logger.info(
                '{:,} of {:,} batch(es) are done; the release files are '
                'written once they all are.'.format(
                    len(batches) - len(unfinished), len(batches)))
            return

        conflicts = concatenate(
            [os.path.join(batch, CONFLICT_NAME) for batch in batches],
            os.path.join(out_dir, CONFLICT_NAME))
        self.logger.info('Wrote {:,} rows to {}.'.format(
            conflicts, os.path.join(out_dir, CONFLICT_NAME)))

        # a row per genome the two both called, which for r237 is 786,144 of its
        # 1.35M genomes; gzipped for the same reason the summary is
        comparison = os.path.join(out_dir, COMPARISON_NAME)
        written = concatenate(
            [os.path.join(batch, COMPARISON_NAME) for batch in batches],
            comparison, compress=True)
        self.logger.info('Wrote {:,} rows to {}.'.format(written, comparison))

        # a row per genome of the release, which for r237 is 116 MB of text and
        # 34 MB compressed; prodigal reads it either way, by its first two bytes
        summary = os.path.join(out_dir, release_summary_name(self.prefix))
        written = concatenate(
            [os.path.join(batch, summary_name(self.prefix)) for batch in batches],
            summary, compress=True)
        self.logger.info('Wrote {:,} rows to {}.'.format(written, summary))

        # a run of an older version left the same summary here uncompressed, and
        # two files a genome apart of which one is stale is how the wrong table
        # gets called; the one just written supersedes it
        stale = os.path.join(out_dir, summary_name(self.prefix))
        if os.path.exists(stale):
            os.unlink(stale)
            self.logger.info(
                'Removed {}, which an earlier run wrote uncompressed and {} now '
                'replaces.'.format(stale, summary))

        self.report_comparison(batches, conflicts)

        # after the concatenation, which has just rewritten the release file from
        # the batches and so has just removed any CheckM2 columns a previous run
        # added: the annotation is applied to the file as it stands, and applying
        # it twice would otherwise double the columns
        self.estimate_conflict_quality(batches, out_dir)

    def checkm2_table(self,
                      accessions: Sequence[str],
                      fastas: Dict[str, str],
                      checkm2_dir: str,
                      table: int) -> Dict[str, Tuple[str, str]]:
        """Estimate the quality of every genome disputing one table, under it.

        A run whose report is already there is not made again, for the reason a
        batch already predicted is not predicted again: the release is finished
        by whichever machine happens to be last, and aggregate() is reached by
        every run over a finished output directory.

        Parameters
        ----------
        accessions : sequence of str
            Genomes to run under this table.
        fastas : dict
            Accession to the genomic FASTA of that genome.
        checkm2_dir : str
            Directory the runs of the release are kept in.
        table : int
            Translation table genes are called under.

        @return: accession to (completeness, contamination), empty where the run
                 could not be made or failed.
        """

        table_dir = os.path.join(checkm2_dir, CHECKM2_TABLE_DIR.format(table))
        report = os.path.join(table_dir, CHECKM2_REPORT)

        if os.path.exists(report):
            self.logger.info(
                'Table {}: {} has already run over these {:,} genome(s); reading '
                '{}.'.format(table, CHECKM2_BIN, len(accessions), report))
            return read_checkm2_report(report)

        staged = stage_checkm2_input(
            accessions, fastas,
            os.path.join(checkm2_dir, CHECKM2_INPUT_DIR,
                         CHECKM2_TABLE_DIR.format(table)))

        if len(staged) != len(accessions):
            self.logger.warning(
                'warning: {:,} of {:,} genome(s) disputing table {} have no '
                'genomic FASTA to estimate the quality of; they keep their row '
                'and are reported as {}.'.format(
                    len(accessions) - len(staged), len(accessions), table, NCBI_NA))

        if not staged:
            return {}

        cmd = checkm2_command(staged, table, table_dir, self.cpus)
        self.logger.info(
            'Table {}: estimating the quality of {:,} genome(s).'.format(
                table, len(staged)))
        self.logger.info('Command: {} ... ({:,} genomes)'.format(
            ' '.join(cmd[:cmd.index('--input') + 1]), len(staged)))

        silent = getattr(self.logger, 'is_silent', False)
        proc = subprocess.run(cmd,
                              stdout=subprocess.DEVNULL if silent else None,
                              stderr=subprocess.STDOUT if silent else None)

        # a failed run costs the release four columns for the genomes of one
        # table and nothing else: the conflicts are found, written and counted
        # before this runs, and they are what the command is for
        if proc.returncode != 0:
            self.logger.error(
                'error: {} returned exit code {} for table {}; those genomes are '
                'reported as {} and the run can be repeated by running the '
                'command again.'.format(CHECKM2_BIN, proc.returncode, table,
                                        NCBI_NA))
            return {}

        return read_checkm2_report(report)

    def estimate_conflict_quality(self, batches: Sequence[str], out_dir: str) -> None:
        """Add to the release conflict file what CheckM2 makes of each genome
        under each of the two tables in dispute.

        Run once for the release rather than once per batch. The conflicts are a
        few hundred genomes of a million-odd, which is two or three per batch,
        and CheckM2 loads its models and searches the whole DIAMOND database once
        per run whatever the run holds; grouped by table the whole release is
        three runs, and per batch it would be hundreds.

        Parameters
        ----------
        batches : sequence of str
            Every batch directory of the run.
        out_dir : str
            Output directory of the run.

        @return: None
        """

        conflict_file = os.path.join(out_dir, CONFLICT_NAME)

        # the conflicts themselves are found, written and counted before this
        # runs, so a file this cannot read costs the release four columns and not
        # the run -- which after days on five machines is the difference that
        # matters
        try:
            rows = read_conflicts(conflict_file)
        except (BadConflictFile, OSError) as exc:
            self.logger.error(
                'error: the quality of the conflicting genomes could not be '
                'estimated: {}'.format(exc))
            return

        if not rows:
            self.logger.info(
                'No genome of the release is a conflict, so there is no quality '
                'to estimate.')
            return

        tables = conflict_tables(rows)
        fastas = release_fastas(batches, [row[0] for row in rows])

        self.logger.info(
            'Estimating with {} the quality of {:,} conflicting genome(s) under '
            'each of the {:,} table(s) in dispute: {}.'.format(
                CHECKM2_BIN, len(rows), len(tables),
                ', '.join(str(table) for table in sorted(tables))))

        checkm2_dir = os.path.join(out_dir, CHECKM2_DIR)
        quality = {table: self.checkm2_table(accessions, fastas,
                                             checkm2_dir, table)
                   for table, accessions in sorted(tables.items())}

        annotated = annotate_conflicts(rows, quality)
        write_table(annotated, conflict_file, header=CONFLICT_HEADER_CHECKM2)

        first = CONFLICT_HEADER_CHECKM2.index(CHECKM2_COLUMNS[0])
        estimated = sum(1 for row in annotated
                        if NCBI_NA not in row[first:first + len(CHECKM2_COLUMNS)])
        passing = {column: sum(1 for row in annotated
                               if row[CONFLICT_HEADER_CHECKM2.index(column)] == 'True')
                   for column in ('pass_qc_gtranslate_tt', 'pass_qc_ncbi_tt')}

        self.logger.info(
            'Wrote {} with the completeness, the contamination and the GTDB QC '
            'verdict under both tables for {:,} of {:,} conflicting genome(s); '
            '{:,} pass QC under the table gTranslate predicted and {:,} under the '
            'table NCBI declares.'.format(
                conflict_file, estimated, len(rows),
                passing['pass_qc_gtranslate_tt'], passing['pass_qc_ncbi_tt']))

        # only now: the working directory goes once what it was for is in the
        # file, and never before, so that nothing is swept up that has not landed
        incomplete = sorted(table for table, estimates in quality.items()
                            if not estimates)
        if incomplete:
            self.logger.info(
                '{} is kept: {} produced nothing for table(s) {}, and the tables '
                'that did are read rather than run again when the command is run '
                'again.'.format(checkm2_dir, CHECKM2_BIN,
                                ', '.join(str(table) for table in incomplete)))
            return

        remove_checkm2_dir(checkm2_dir)
        self.logger.info(
            'Removed {}: the called proteins and the DIAMOND output it held say '
            'nothing {} does not.'.format(checkm2_dir, conflict_file))

    def run(self, gtdb_genome_path_file: str, taxonomy_file: str, out_dir: str) -> bool:
        """Predict the translation table of every genome of a release.

        Parameters
        ----------
        gtdb_genome_path_file : str
            genome_dirs file of the release, accession and genome directory per line.
        taxonomy_file : str
            Standardised NCBI taxonomy file, for the comparison.
        out_dir : str
            Directory the batches and their results are written to.

        @return: True where every batch this machine took has finished, False
                 where one failed and is left to a later run.
        """

        binary = shutil.which(GTRANSLATE_BIN)
        self.logger.info('Using {}.'.format(binary))

        taxonomy = read_taxonomy(taxonomy_file)
        self.logger.info('Read the NCBI taxonomy of {:,} genomes from {}.'.format(
            len(taxonomy), taxonomy_file))

        batches = self.plan_batches(gtdb_genome_path_file, out_dir)

        done, held, failed = 0, 0, 0
        for index, batch_dir in enumerate(batches, start=1):
            label = 'Batch {:,} of {:,} ({})'.format(
                index, len(batches), os.path.basename(batch_dir))

            if batch_state(batch_dir) == STATE_SUCCESS:
                self.logger.info('{}: already finished, skipping.'.format(label))
                continue

            if not claim_batch(batch_dir, self.reclaim, self.lease):
                owner = read_canary(os.path.join(batch_dir, RUNNING_CANARY))
                held += 1
                self.logger.info('{}: held by {} since {}, last heard from {}, '
                                 'skipping.'.format(
                                     label, owner.get('host', 'another machine'),
                                     owner.get('time', 'an unknown time'),
                                     age_phrase(claim_age(
                                         os.path.join(batch_dir, RUNNING_CANARY)))))
                continue

            # the batch has its own log from here, since this is where anything
            # happens to it and every machine of a run writes its own --log
            with batch_log(batch_dir, self.logger):
                self.logger.info('{}: starting.'.format(label))
                try:
                    with Heartbeat(os.path.join(batch_dir, RUNNING_CANARY),
                                   self.heartbeat):
                        if already_predicted(batch_dir, summary_name(self.prefix)):
                            self.logger.info(
                                '{}: gTranslate has already run over it; comparing '
                                'what it wrote.'.format(label))
                        else:
                            self.run_gtranslate(batch_dir)
                        counts = self.compare_batch(batch_dir, taxonomy)
                except KeyboardInterrupt:
                    # nothing was decided about the batch, and the machine that
                    # held it is stopping, so it is handed back rather than left
                    # to sit out its lease
                    release_claim(batch_dir)
                    self.logger.error('{}: interrupted; the claim is given up and '
                                      'the batch carries on where it stopped.'.format(label))
                    raise
                except Exception as exc:
                    failed += 1
                    fail_batch(batch_dir, str(exc))
                    self.logger.error('{}: failed and will be retried by a later '
                                      'run: {}'.format(label, exc))
                    continue

                # the counts go in the canary because the release line adds them
                # up across batches, and a batch run on another machine leaves
                # them nowhere else
                finish_batch(batch_dir,
                             compared=counts.compared,
                             conflicts=counts.conflicts,
                             no_ncbi_table=counts.no_ncbi_table)
                done += 1
                self.logger.info('{}: done.'.format(label))

        self.logger.info(
            '{:,} batch(es) finished here, {:,} held by another machine, '
            '{:,} failed.'.format(done, held, failed))

        self.aggregate(batches, out_dir)

        # a batch that failed has already said why, in its own log and in its
        # FAILED file, and the batches that succeeded are finished and staying
        # that way; there is nothing left for a traceback to add, and a run of
        # five machines over days should not end by printing one
        if failed:
            self.logger.error(
                '{:,} batch(es) failed; they are the directories holding a {} '
                'file and are retried by running the command again.'.format(
                    failed, FAILED_CANARY))
            return False

        return True
