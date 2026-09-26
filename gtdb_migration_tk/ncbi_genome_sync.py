r"""
ncbi_genome_sync.py — keep an NCBI genome mirror equal to a GTDB selection: a manifest-driven
sync that first removes what the selection does not list, with optional md5 verification.

WHY NOT wget
------------
`wget -r [-N]` costs 23 HTTP requests per genome (1 directory index + 22 files) whether
or not -N is used; -N only changes the responses from 200 to 304, saving bytes but not
requests. NCBI throttles on requests, so that is the constraint that matters.

md5checksums.txt answers "what files exist" and "did they change" in ONE request, so an
unchanged genome costs a single request here instead of 23. A cold genome costs 12 with
the whitelist below; the 23 is the old full-mirror figure, kept for comparison.

WHY NOT rsync
-------------
It is gone. NCBI retired rsync access to ftp.ncbi.nlm.nih.gov on 2026-06-01, "to ensure the
security and long-term stability of our data delivery infrastructure"; only FTP and HTTPS
remain, and the old shell pipeline's rsync form of these URLs no longer connects:
    https://ncbiinsights.ncbi.nlm.nih.gov/2026/03/25/retire-rsync-support-ftp-downloads/
rsync's single-connection file-list exchange would have been the zero-request answer to
the throttling described below. HTTPS driven by the manifest is the closest thing left.

FILE WHITELIST
--------------
Only the files in WANTED_EXACT / WANTED_SUFFIXES are mirrored; everything else NCBI
publishes is ignored -- not downloaded, not verified, not reported missing. Measured on
a GenBank batch: 26 MB per 10 genomes against ~133 MB for a full mirror.

Names are built exactly (assembly name + suffix) and matched against the whole manifest
path, never with endswith() or basename():
  * "_genomic.fna.gz" is also a suffix of "_cds_from_genomic.fna.gz" and
    "_rna_from_genomic.fna.gz", so a suffix test silently re-admits both;
  * every wanted file lives at the genome root, so comparing the full path excludes the
    _assembly_structure/ , all_assembly_versions/ and representative/ subtrees
    structurally -- no directory blacklist to keep in sync as NCBI adds layouts.
    (1353 of 42323 entries across 3000 manifests are nested; none is a wanted name.)

FAST PATH -- AND WHAT IS TRUSTED
--------------------------------
The fetched manifest is compared to the locally installed one ENTRY BY ENTRY, not as a
blob. A manifest is only ever installed after every whitelisted file it lists has been
verified against it, so for each file:

    entry unchanged since last sync, file present  -> trusted, not re-hashed
    entry changed (re-annotation), file present    -> hashed; downloaded if it differs
    file missing, whatever the manifest says       -> downloaded

An unchanged genome therefore costs zero hashing and zero downloads (just ~12 stat calls)
-- and, within --max-age of its last confirmed sync, zero requests (see RESTART AND
FRESHNESS). A re-annotated RefSeq genome hashes only the entries that moved -- usually the
gff/gbff, not the multi-MB genomic.fna -- and a missing file is ALWAYS noticed. That last
point matters: the previous blob comparison fast-pathed straight past a deleted file as
long as the manifest matched, so a partial --delete or a newly added WANTED_SUFFIXES entry
was invisible until --full. Now a whitelist addition backfills itself on the next sync.

The trust rule relies on NCBI regenerating the manifest whenever content changes, and it
cannot see local bit-rot in a file whose entry did not move. --full re-hashes everything
and repairs; --verify re-hashes everything and reports. That is the documented policy.

WHEN NCBI'S TWO CHECKSUM TABLES DISAGREE
----------------------------------------
NCBI publishes md5checksums.txt and, for some assemblies, uncompressed_checksums.txt, and
they can contradict each other. Measured on the production mirror: ~1000 genomes whose
md5checksums.txt entry for <asm>_fcs_report.txt is stale while the other table is right.
Believing md5checksums.txt alone condemns those genomes twice over. The sync cannot
install a file it has just downloaded -- the bytes NCBI serves do not match the checksum
NCBI publishes for them -- so the genome fails on every run and never enters the mirror at
all; and a verify of one that was mirrored by some other route reports rot that is not
there.

The two files look nothing alike, and that is a trap worth stating plainly, because
assuming otherwise fails SILENTLY rather than loudly -- a table read with the wrong reader
parses to nothing, which is indistinguishable from "NCBI vouches for nothing here":

    md5checksums.txt           <md5>  ./<name>          two spaces, no header
    uncompressed_checksums.txt ./<name>\t<md5>\t...     TAB-separated, name FIRST, under a
                                                       `#file md5sum crc32 size` header

read_uncompressed_checksums() reads the second by column name, never by position, and
returns None rather than guessing when the header is missing.

So an UNCOMPRESSED file md5checksums.txt rejects gets a second opinion, in both the sync
and the verify: if uncompressed_checksums.txt vouches for exactly the bytes in hand, under
exactly that name, the file passes. A file neither table vouches for fails exactly as
before.

Only uncompressed files, because the table lists the uncompressed FORM of everything: a
.gz appears there under a name it does not have on disk (<asm>_genomic.fna against
<asm>_genomic.fna.gz), so no entry for it can ever match. Asking is a request that cannot
succeed, so a compressed file is refused before one is made and a genome whose only
failures are archives spends nothing. What remains eligible is what NCBI does not compress
-- <asm>_fcs_report.txt, the assembly and ANI reports -- which means this can never excuse
a corrupt archive, structurally rather than by policy.

It is a SECOND opinion, never a first -- md5checksums.txt settles every file it agrees
with, and the other table is fetched lazily, one request on the first eligible
disagreement in a genome and none at all for a genome that is clean. A verify of a healthy
mirror therefore still makes no requests.

The table is not trusted from disk: a local copy vouching for local bytes is two halves of
the same claim, and NCBI may since have corrected either table, so the answer comes from
NCBI each time. The copy is written into the genome directory when it settles something,
as the record of why those bytes were accepted -- which is why some genome directories
hold an uncompressed_checksums.txt and most do not. It is not in the whitelist, nothing
requires it, and its absence is never a fault. Each run's closing lines say how many files
and genomes it settled, so a change in NCBI's consistency is visible rather than silent.

CORRECTNESS -- AND BEING KILLED
-------------------------------
Every download is written to a temp file and md5-checked BEFORE being moved into place,
so a throttled or truncated transfer can never leave a bad file in the mirror. The
manifest is installed last, only once every file it lists is in place. So a run killed at
ANY instruction (SIGKILL, OOM, a node reboot) leaves each genome in one of two states --
complete with a manifest, or partial without one -- and the next run takes the fast path
on the first and hashes-then-completes the second, downloading nothing that was already
correct. Stale temp files from the kill are swept when the genome is next visited.

A stop the script can see (Ctrl-C, SIGTERM, a full disk) is cooperative: every worker
aborts at its next request, in-flight genomes are recorded in <log>.fail as "stopped:",
the summary prints, and the exit code says why (74 / 130 / 143 below). A second signal
kills outright. <root>/.ncbi_sync.lock refuses a second sync on the same mirror, so a
restart cannot race a run that is merely paused.

REMOVAL -- THE SELECTION DEFINES THE MIRROR
-------------------------------------------
The table select_genomes writes is not a list of genomes to fetch but a statement of what
--root should contain. Before anything is fetched the mirror is walked and every genome
directory the table does not name is removed: a genome NCBI suppressed, a version it
superseded, an assembly it renamed. The match is on the exact directory the ftp_path maps
to (genome_relpath), so GCF_x.1_ASMold goes when the table names GCF_x.1_ASMnew rather
than surviving beside it, which a match on accession alone would allow. The digit-triplet
directories a removal empties go too, down to all/, so list_genomes never finds a hollow
path.

Each directory is written to <log>.rm BEFORE it is deleted: a kill mid-way leaves a
record of what went, and Ctrl-C / SIGTERM stops between directories. One that will not
delete (another member's file -- see SHARED OPERATION) is logged and counted and the sync
still runs, with exit 1.

Removing a genome is ~18 metadata operations -- a scandir, ~12 unlinks, and the rmdirs up
the triplets it empties -- and on NFS every one is a round trip the client waits out, so a
serial pass over half a million genomes is hours of idle time. Two things follow. The
climb up the triplets asks rmdir instead of listing first (rmdir already refuses a
directory that is not empty), which is four fewer round trips per genome AND removes the
window between the test and the act. That in turn makes the removals safe to run on
--nfs-jobs threads: two genomes under one triplet race for it, and the loser gets
ENOTEMPTY or ENOENT and stops climbing. NFS does not throttle the way NCBI does -- the
ceiling is the server's nfsd pool, not a penalty -- so --nfs-jobs is a throughput knob
rather than a safety one: 8 by default, lower it when others are working on the mirror.
Measured on the production server: 113 genomes/s serial, 316 at 4, 465 at 8, 539 at 16. Rows land in <log>.rm in
completion order rather than sorted.

Two guards, because this deletes data:
  * A table given to --gtdb_selected_genomes that lists no genomes is refused (exit 2):
    taken literally it would empty the mirror, and the likely cause is an empty .fail
    given to the wrong flag. --retry treats the same table as "nothing to do" (exit 0).
  * --retry never removes a genome it was not given. A .fail, a .bad, or any hand-cut
    subset of the selection lists only some of the genomes the mirror should hold, so
    what it omits is not "extra". With --delete it still removes the genomes it WAS
    given, one directory at a time, immediately before refetching that genome (below).

--delete means "remove the directory rather than repair it", and what it is given decides
which directories:

  * --retry --delete rebuilds. Each listed genome's directory is removed immediately
    before that genome is fetched, in the same run, by the same worker. A .bad is the
    case it exists for: verification found the bytes on disk wrong, and a sync would
    otherwise trust an intact manifest and repair nothing. Removing per genome rather
    than in an up-front pass means a Ctrl-C or a tripped breaker leaves at most -j
    directories gone and not yet rebuilt, and those are exactly the ones the run records
    as "stopped:" in <log>.fail. A directory that will not delete is recorded as a
    failure and NOT synced: fetching into it would rebuild the state being discarded.
  * --verify/--verify-only --delete repairs the mirror's shape. Genomes that fail
    verification lose their directories (so a later sync rebuilds them), and directories
    the selection does not name are removed and recorded in <log>.rm.
  * A bare --gtdb_selected_genomes sync refuses --delete (exit 2). There it could only
    mean "discard and re-download the whole mirror", which is millions of files and no
    guard; the pruning that sync does need happens anyway.

--verify and --verify-only, given the selection, check both halves of "the mirror equals
the selection": every listed genome present and md5-clean (failures to <log>.bad, as
always), AND nothing else present -- every directory the selection does not name is listed
in <log>.extra and fails the verification (exit 1); with --delete they are removed as
well, recorded in <log>.rm. Given --retry only the listed genomes are verified, a retry
file being a subset. --verify-only never syncs, and removes nothing without --delete.

--dry-run walks the mirror and reports the three counts -- directories to remove, genomes
to add, genomes already present (whether a present one needs updating is decided per
genome by the sync, from the manifest) -- writes the removal list to <log>.rm_dry_run,
and exits without fetching, removing or verifying anything. It takes no lock, so it can
run beside a live sync to show what the NEXT run would do. The walk is one os.scandir()
per directory, on the order of 1.4M of them on NFS for a full release, overlapped across
--nfs-jobs: minutes, not hours.

RESTART AND FRESHNESS -- .last_synced
-------------------------------------
Every request costs the same under --rate, so a restart that re-asks NCBI about 400,000
genomes it finished yesterday spends ~7 hours learning nothing. Each genome therefore
carries a marker, <genome>/.last_synced, holding the UTC time this script last confirmed
it complete AGAINST NCBI -- not when NCBI last changed it. It is written LAST, after the
manifest and the unlisted files, only when nothing failed, so its existence means "the
whole sync of this genome completed"; a crash anywhere before it leaves no marker and the
genome is simply visited again (one request, fast path). --verify withdraws it on any
failure, so the next sync must consult NCBI for that genome rather than skip it as fresh.
(What that sync then repairs is unchanged by the marker: a missing file is re-downloaded;
bit-rot in a file whose manifest entry did not move is still trusted, and still needs
--full, or --verify to find it and --retry <log>.bad --delete to rebuild it.)

With --max-age DAYS (default 14) a genome whose marker is younger is skipped with ZERO
requests -- but only after the local checks that cost nothing: the manifest must be
present and every whitelisted file it lists must exist. The marker saves the round-trip,
never the evidence. The summary row is still applied (assembly_status.txt is rewritten if
version_status or the anomalies moved), also without a request. --full ignores markers
and rewrites them; --max-age 0 always asks NCBI.

This is a FRESHNESS POLICY, not only a restart optimisation: a RefSeq re-annotation that
lands inside the window is not seen until it expires. For a mirror refreshed monthly that
is invisible; a weekly schedule wants --max-age 5 or so. The summary line prints how many
genomes were "fresh, no request" so a fast run is never mistaken for a thorough one.

Cost of a skip: one small read, one manifest read, ~12 stats -- milliseconds on cold NFS,
in parallel across -j with no rate limit involved. The first run after deploying this
visits everything once (no markers exist yet); it is the last full-price pass.

NOT EVERY FILE IS CHECKSUMMED
-----------------------------
Four wanted files are absent from md5checksums.txt for essentially every genome:
    md5checksums.txt              (a manifest cannot list its own checksum)
    assembly_status.txt
    <asm>_ani_report.txt
    <asm>_ani_contam_ranges.tsv

assembly_status.txt is not fetched at all when the summary carries version_status: NCBI
generates that file from the same assembly record the summary comes from, so the script
writes it locally from the row (render_status) -- "status=<version_status>" plus an
"assembly anomaly=" line for any ASSEMBLY_ANOMALIES term in excluded_from_refseq. It
used to be refetched every run, because a genome can be suppressed or replaced with no
data file changing, and that one file was half the warm fast-path cost. A table without
version_status (a hand-cut one) still gets the fetch. The ANI files are fetched only when
the manifest moved or the file is missing locally.

A fifth file is division-specific: GenBank manifests omit <asm>_fcs_report.txt while the
server still returns it (measured 200 for 7 of 8 GCA genomes that omitted it).

That rule has a cost: a wanted file NCBI simply does not have is missing locally every
run, so it is probed (404) every run. It used to bite on <asm>_genomic.gff.gz (absent for
every unannotated assembly -- about half of GenBank) and <asm>_wgsmaster.gbff.gz (absent
for every non-WGS assembly -- 17% of RefSeq too), adding a request to the fast path of
most genomes. Both are now in MANIFEST_IS_TRUTH_SUFFIXES: measured on the server iff
listed, so an unlisted one is not probed. The fast path is one request in both divisions;
the residual case is the rare genome with no <asm>_fcs_report.txt (2 of 110).

RATE LIMITING AND THE CIRCUIT BREAKER
-------------------------------------
NCBI throttles (HTTP 503) on SUSTAINED REQUEST RATE over a multi-minute window -- not on
connections, and not on -j. Measured 2026-09-10 on 1000 GenBank genomes:

    cold  -j8 / -j9   ~19-21 req/s for 10 min   clean
    cold  -j10        ~29 req/s                  clean 6 min, then 44 x 503, 2 genomes lost
    warm  -j9, 60 s   ~42 req/s                  clean (the window had not filled)
    warm  -j9, 38 min ~42 req/s attempted        2,500 x 503, 293 genomes lost, never recovered

The same -j gives wildly different rates because a warm request (unchanged genome: the
manifest alone, ~1 KB) completes ~10x faster than a cold one. Production re-syncs are
mostly warm and run for hours, so -j alone cannot be made safe. Two mechanisms replace it:

  * --rate: a token bucket shared by all workers caps requests/second directly (default 20,
    the measured clean rate; 16 keeps a 20% margin under it). -j then only sets bandwidth
    parallelism.
  * The circuit breaker: any 429/503 pauses EVERY worker (60 s, honouring Retry-After).
    Per-request backoff demonstrably cannot recover -- one worker sleeping leaves eight
    hammering, and in every throttled pass ~93% of first-minute requests were 503s until
    stacked backoffs happened to starve the aggregate rate. If a pause is answered by
    another 503 within 30 s of resuming, it did not work: the next doubles (max 300 s).
    After 4 such failed pauses in a row the host is in a penalty state that minutes will
    not clear (measured: 60/120/240/300 s pauses each answered within 1-17 s, 67 min after
    the load that caused it) and the run STOPS with exit 75 rather than crawling. Done
    genomes stay on disk, in-flight ones go to <log>.fail as "stopped:", and a re-run of
    the same list hours later is cheap because everything done takes the fast path.

In practice, then: --rate (default 20) is the safety mechanism and -j (default 8) is not.
Cold -j8 and -j9 both measured clean and -j10 throttled, so the default sits one step below
the highest value ever measured clean rather than on it, and nothing above 9 buys anything.
With --rate 0 that margin means nothing: the warm -j9 run above was at ~42 req/s and had
already lost 36 of its first 1000 genomes. Raise --rate only from a run whose log ended
with throttled=0, and re-check it the same way afterwards.

Short tests prove nothing: 200-genome bursts were clean at every -j from 6 to 10. Only a
sustained run (>= 1000 genomes, >= 10 min) after a >= 15 min rest means anything, and the
log's throttled= / trips= counters are the signal.

TUNING -- WHAT NOT TO CHANGE
---------------------------
Do NOT call setsockopt(SO_RCVBUF) on these sockets, and do not "fix" throughput by
raising a TCP buffer in Python. Setting SO_RCVBUF disables Linux receive-buffer
autotuning and pins the window to net.core.rmem_max (212992 on this host), which
measured 4.4x SLOWER on the ~206 ms NCBI path:

    autotuned (default)     kernel gave  235392 B  ->  2.31 MB/s
    explicit  4 MB request  kernel gave  425984 B  ->  0.52 MB/s
    explicit 32 MB request  kernel gave  425984 B  ->  0.52 MB/s

NCBI's published "set the TCP buffer to 32 MB" advice therefore backfires here: it
needs net.core.rmem_max raised as root FIRST. Left alone, autotuning grows toward
tcp_rmem's 6 MB ceiling, comfortably above the ~2 MB bandwidth-delay product at this
RTT, so there is most likely nothing to gain even after a sysctl change.

What does matter is connection reuse, which this script already gets from one
persistent HTTPS connection per worker thread: the same 4.4 MB file takes 1.86 s on a
fresh connection (2.35 MB/s) versus 0.42 s on a warm one (10.48 MB/s). Only 0.42 s of
that 1.44 s gap is the TLS handshake -- the rest is TCP slow-start, which restarts on
every new connection and never reaches full window on a short transfer. Anything that
reintroduces per-file connections (spawning wget/curl, a segmented downloader) gives
that back.

--verify-jobs is a different resource entirely: md5 over local files, never NCBI. Measured
on a page-cache-warm mirror it plateaus at 16-20 workers and degrades slightly beyond, so
raising it to the core count does not help; on cold NFS it is I/O-bound rather than
CPU-bound anyway. The default of 20 is the top of the flat region.

INPUT
-----
The input is an NCBI assembly_summary.txt: the ftp_path column says where each
assembly lives, and assembly_accession names it. Both are found BY NAME from the
`#assembly_accession ...` header row, never by column number -- NCBI has grown that file
from 23 columns to 38, and pinning ftp_path to field 20 breaks silently on the next one.
The reading of those tables is in ncbi_utils.py, shared with select_genomes.py.

    ftp_path "na" or empty   the assembly has no public directory (suppressed, or
                             not yet released). NOT an error: whole-domain summaries
                             normally carry some. Every such row is named in the log
                             and counted in a stderr warning, then skipped.
    version_status           "latest", "replaced" or "suppressed". A replaced/suppressed
                             genome is SYNCED ANYWAY -- NCBI still serves the directory,
                             and pinned version lists (GTDB releases) want exactly that
                             version -- but each is named in the log and counted on
                             stderr. Also what assembly_status.txt is written from.
    excluded_from_refseq     only its ASSEMBLY_ANOMALIES terms reach assembly_status.txt;
                             the rest ("derived from metagenome") are exclusion reasons
                             NCBI does not write there.
    ftp:// path              rewritten to https:// (older summaries publish ftp URLs)
    same ftp_path twice      deduplicated
    anything not under       hard error before a single request is made, because the
    https://ftp.ncbi.../     old shell pipeline's rewritten URLs would otherwise
    genomes/                 mirror into the wrong directory (see genome_relpath)

<log>.fail and <log>.bad open with the same column names -- plus version_status and
excluded_from_refseq, so a retry can still write assembly_status.txt without a request --
and feed straight back in with no cut/awk step. Each then adds one column of its own,
`reason` and `failed_files`, which the reader ignores: it reads by name, so a table may
carry extra columns. Nothing else is accepted: a bare list of URLs has no header, and a
file whose rows start before one is rejected by name.

OUTPUTS
-------
<log> below is -l/--log with .log stripped: -l r95/sync.log writes r95/sync.fail,
r95/sync.bad, r95/sync.rm, r95/sync.rm_dry_run and r95/sync.extra. The log names the run
and the run's outputs carry its name, so one round's whole account of itself shares one
stem and sits in one directory -- give each round its own --log and nothing of one round
can overwrite another's. Any extension that is not .log is kept (`run.txt` ->
`run.txt.fail`), so a log named otherwise cannot collide with a `run.log` beside it.
--fail/--bad override the two that are also inputs.

The consequence to know: --retry reads a file this scheme writes. Retrying a round's own
.fail or .bad under that same --log would have the run truncate the list it is reading, so
it is refused (exit 2) rather than risked -- give the retry its own --log, which is what
you want anyway, since that round deserves its own history.

The run history is NOT a file of this module's own: argv, every row that had no usable
ftp_path, a cumulative progress/throttling snapshot every 60 s, and the final summary all
go to the file named by -l/--log, which appends across runs the same way. NOTHING is logged per genome -- a ~2M-genome run
would otherwise bury the log in near-identical lines -- so the log stays bounded by how
long the run took rather than by how many genomes it covered. Genomes needing action are
listed in <log>.fail / <log>.bad below, and failures print to the console as they
happen; set LOG_LEVEL to DEBUG for per-genome detail when debugging.

    <log>.fail  genomes the sync could not complete, TSV under a `#assembly_accession
                 ftp_path version_status excluded_from_refseq reason` header. TRUNCATED
                 each run. With MAX_TRIES=6 a throttled run WILL populate this; feed it
                 back in with --retry.
    <log>.bad   genomes that failed --verify, TSV under the same header with
                 `failed_files` in place of `reason`: every file of that genome at fault,
                 comma-separated and each tagged with WHICH fault -- missing:NAME (the
                 manifest lists it, the mirror does not have it), mismatch:NAME (present,
                 bytes rotted) or unreadable:NAME (present, could not be hashed; the
                 errno is in the log). The whole genome is checked before the row is
                 written, so the column is every fault, not the first. A fault with no
                 file to name -- missing directory, no md5checksums.txt, unparseable
                 md5checksums.txt -- puts that reason in the column instead. TRUNCATED
                 each run. Feed it back in with --retry --delete, which removes each
                 directory and refetches it in the one run.
    <log>.rm    genome directories removed because the selection does not list them,
                 TSV under `#assembly_accession directory`. Each row is written BEFORE
                 its directory is deleted, so an interrupted run still says what went;
                 a row whose removal then failed is reported in the log. TRUNCATED each
                 run of the selection.
    <log>.rm_dry_run  the same list, written by --dry-run instead of being acted on.
    <log>.extra  written by --verify/--verify-only given the selection: the directories
                 the selection does not list, same columns as .rm. TRUNCATED each run.

    <root>/.ncbi_sync.lock  held (flock) while a run is alive; records pid, host, start.
    <genome>/.last_synced   UTC time this genome was last confirmed complete against
                            NCBI; what --max-age reads. Written last; withdrawn by a
                            failed verify.
    <genome>/uncompressed_checksums.txt  present ONLY in a genome where NCBI's two
                            checksum tables disagreed and this one settled it (see WHEN
                            NCBI'S TWO CHECKSUM TABLES DISAGREE). Most genomes have none,
                            and nothing treats its absence as a fault.

SHARED OPERATION -- PERMISSIONS
-------------------------------
The production mirror is millions of files on NFS, operated by a group (dataadmin), not a
user. This section is written for that tree; the small-tree shortcut is at the end.

TO RUN -- from the --log directory (where the .log/.fail/.bad/.rm/.extra live), as the owner of the
files, after anything was added by another route (mv, rsync -a, cp -p, an old version of
this script). Four lines; the third draws a progress bar in completed subtrees:

    FIX='\( -type d ! -perm 2775 -exec chmod 2775 {} + \) , \( -type f ! -perm 664 -exec chmod 664 {} + \) , \( ! -group dataadmin -exec chgrp dataadmin {} + \)'
    eval find genomes -maxdepth 3 $FIX
    find genomes/all -mindepth 3 -maxdepth 3 -type d -print0 | xargs -0 -P 8 -I@@ sh -c "find \"\$1\" $FIX; echo \"\$1\"" sh @@ | /opt/miniforge/bin/python3 -m tqdm --total $(find genomes/all -mindepth 3 -maxdepth 3 -type d | wc -l) --unit dir --null
    chmod 664 *.log *.fail *.bad *.rm *.extra genomes/.ncbi_sync.lock

Line 2 fixes the shallow levels (genomes, all, GC?, GC?/NNN) and the lock file; line 3
never enters them. Line 3 runs the same three clauses on every genomes/all/GC?/NNN/NNN
subtree, 8 at a time -- units are enumerated by DEPTH, so GCA and GCF (and any division
NCBI adds under all/) are covered alike, with no per-division step. The $(...) counts the
units for the bar; that is a depth-3 walk of a few thousand directories, not the tree.
Everything is silent when there is nothing to fix and a re-run costs one walk with no
changes, so this can sit in cron beside the sync. -P 8 is a safe start: the walk is bound
by NFS round-trip latency, not CPU, so raise it while the server keeps up. To preview
what would be touched, swap each `-exec ... {} +` in FIX for `-print`.

If the files have more than one owner, only root can fix them all. Put sudo on the
command that changes things, not on eval -- eval is a shell builtin, so `sudo eval ...`
fails with "command not found":
    eval sudo find genomes -maxdepth 3 $FIX
    find genomes/all -mindepth 3 -maxdepth 3 -type d -print0 | sudo xargs -0 -P 8 -I@@ sh -c "find \"\$1\" $FIX; echo \"\$1\"" sh @@ | /opt/miniforge/bin/python3 -m tqdm --total $(find genomes/all -mindepth 3 -maxdepth 3 -type d | wc -l) --unit dir --null
    sudo chmod 664 *.log *.fail *.bad *.rm *.extra genomes/.ncbi_sync.lock
$FIX and the globs are expanded by YOUR shell before sudo runs, which is what you want.

WHY. The root carries the setgid bit (drwxrwsr-x), so the kernel hands every new file and
subdirectory to the group and propagates the bit to new subdirectories; the script's part
is to NOT take group write away again. It sets umask 002 before creating anything, makes
directories 0o775 and installs files 0o664 (DIR_MODE / FILE_MODE). Files need no execute
bit: every operation on an EXISTING file another member created -- replace, unlink,
rmtree -- is a write to its directory, so group write on directories is what lets a
second operator run. The --log directory needs the same treatment (setgid,
group-writable), because the log, .fail, .bad, .rm, .extra and genomes/.ncbi_sync.lock follow the same
rule -- one member's 644 lock file would lock everyone else out with exit 2 "cannot use
--root".

WHY THIS SHAPE. Anything created before this rule keeps its old mode (open(..., "w")
reuses the inode), and anything mv'd or rsync -a'd into the tree keeps its old GROUP --
setgid applies only at creation. Each find runs all three tests on every entry in ONE
traversal (the comma operator) and fires a syscall only where something changes. Not
`chmod -R` / `chgrp -R`: directories and files need different modes (2775 on a file sets
execute and setgid on it), and both -R forms issue their syscall for EVERY entry, changed
or not -- measured 1021 chmod(2) and 1021 fchownat(2) on a 1021-entry tree that was
already (or mostly) correct, against 0 and 970 for the find. Multiplied by millions of
NFS round trips, that is the difference between a pass you can re-run nightly and one you
schedule. chgrp does not clear setgid on a directory (measured), so clause order does not
matter. Parallelising by subtree overlaps those round trips; the walk itself is what
costs, so the units are chosen to be numerous (thousands) and independent.

SMALL TREES. On a test mirror of a few hundred entries the parallel form is slower (0.035
s against 0.022 s: sixteen process spawns cost more than they save). There, the same
clauses as one serial find do the whole job and are easier to read:
    find genomes \( -type d ! -perm 2775 -exec chmod 2775 {} + \) , \( -type f ! -perm 664 -exec chmod 664 {} + \) , \( ! -group dataadmin -exec chgrp dataadmin {} + \)
    chmod 664 *.log *.fail *.bad *.rm *.extra genomes/.ncbi_sync.lock

EXIT CODES
----------
    0   everything synced / verified clean
    1   some genomes failed (see <log>.fail / <log>.bad), or verification found
        directories the selection does not list (see <log>.extra)
    0   also: a header-only table given to --retry (an empty .fail at the end of a loop)
    2   usage error or bad --jobs/--rate/--delete combination; a table with no header
        row, a malformed ftp_path, or rows but none with a usable ftp_path; a selection
        that lists no genomes (it would empty --root -- a retry file belongs to --retry)
    74  the filesystem refused -- disk full, quota, read-only (EX_IOERR). Nothing is
        half-installed; fix it and re-run the SAME summary
    75  NCBI is refusing this host and pauses did not clear it (EX_TEMPFAIL): rest for
        hours, then re-run the SAME summary. Also: another sync already holds --root
    130 SIGINT (Ctrl-C), 143 SIGTERM: in-flight genomes are in <log>.fail as
        "stopped:"; re-run the SAME summary

USAGE
-----
    S=gtdb_selected_genomes.tsv.gz                                        # from select_genomes
    ./ncbi_genome_sync.py --gtdb_selected_genomes $S --root genomes -l sync.log --dry-run
    ./ncbi_genome_sync.py --gtdb_selected_genomes $S --root genomes -l sync.log  # remove, then sync
    ./ncbi_genome_sync.py --gtdb_selected_genomes $S --root genomes -l sync.log --verify
    ./ncbi_genome_sync.py --gtdb_selected_genomes $S --root genomes -l sync.log --full
    ./ncbi_genome_sync.py --retry sync.fail --root genomes -l retry1.log
    ./ncbi_genome_sync.py --retry sync.bad  --root genomes -l rebuild.log --delete

  Each round gets its own --log, and so its own .fail/.bad: the first line above reads
  sync.fail and writes retry1.fail, whose leftovers the next round reads. Reusing
  -l sync.log there would have the run truncate the file it is reading, and is refused.

  A subset is just a subset of the table: keep the '#assembly_accession ...' header line
  and grep/awk out the rows you want -- and give it to --retry, because as
  --gtdb_selected_genomes a subset would remove every genome it leaves out.

  --rate (default 20 req/s) is the setting that keeps NCBI from throttling -- see RATE
  LIMITING above. -j/--jobs (default 8) only sets bandwidth parallelism once --rate is on.
  --max-age (default 14 days) is what makes a restart cheap -- see RESTART AND FRESHNESS.
  --nfs-jobs (default 8) sizes the mirror walk and the removal, which touch the local
  mirror rather than NCBI and so are bounded by NFS latency, not by --rate.
  --verify-jobs (default 20) sizes the verification pass separately, since md5 is local
  and CPU-bound rather than NCBI-limited -- see TUNING. --root is the mirror directory and
  has no default; "genomes" is the usual name. -l/--log names the log AND the run: the
  .fail/.bad/.rm/.rm_dry_run/.extra are that name with .log stripped, written beside it.
"""

import os
import re
import sys
import argparse
import calendar
import collections
import email.utils
import errno
import fcntl
import hashlib
import http.client
import logging
import shlex
import shutil
import signal
import ssl
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor

from tqdm import tqdm, __version__ as tqdm_version

from gtdb_migration_tk.ncbi_utils import (ASSEMBLY_STATS_EXT, CHUNK, GENOME_COLUMNS,
                                          GENOMIC_FASTA_EXT,
                                          MD5_MANIFEST, NCBI_HOST, NCBI_URL, BadInput,
                                          file_md5, has_ftp_path, read_md5_manifest,
                                          read_summary_rows, summary_field, table_header)


HOST = NCBI_HOST
URL_PREFIX = NCBI_URL + "/genomes/"
FTP_PREFIX = "ftp://" + NCBI_HOST + "/genomes/"     # older summaries; rewritten to https

RETRY_STATUS = frozenset((429, 500, 502, 503, 504))

# 6 attempts with 1.7x backoff (sleeping 2, 3.4, 5.8, 9.8, 16.7 s between them) give up
# on a stuck file after ~38 s. Higher values idle a whole worker on one bad URL for far
# longer without helping: genuinely dead URLs return 404, which is not retried at all.
#
# This is deliberately fail-fast: under heavy 503s it turns throttling into hard failures
# rather than riding them out, so a throttled run WILL leave genomes in <log>.fail that
# need a second pass. That is the intended trade -- no worker is held hostage to one URL
# -- but it means .fail must actually be reprocessed. Raise MAX_TRIES if you would rather
# a run take longer and finish complete.
MAX_TRIES = 6
MAX_BACKOFF = 60.0

# Per-genome marker: when this script last confirmed the genome complete AGAINST NCBI.
# Written last, after the manifest, only on full success; a genome whose marker is younger
# than --max-age is not asked about at all (see RESTART AND FRESHNESS in the header).
LAST_SYNCED = ".last_synced"
LAST_SYNCED_FORMAT = "%Y-%m-%dT%H:%M:%SZ"

# Group-operable mirror (see SHARED OPERATION in the header): directories rwxrwxr-x, files
# rw-rw-r--, umask 002. The setgid bit on the root hands new entries to the group; these
# make sure the script does not take group write away again. Every operation on an
# EXISTING file -- replace, unlink, rmtree -- is a directory write, so DIR_MODE is what
# lets a second operator run; data files never need execute.
DIR_MODE = 0o775
FILE_MODE = 0o664
UMASK = 0o002

# A temp file this old beside a genome was left by a killed run: a live download rewrites
# its temp's mtime every CHUNK, and a stalled socket dies at the 120 s connection timeout.
TEMP_STALE_S = 3600.0

# errno values that mean the FILESYSTEM refused, not the network. Retrying these
# re-downloads the same bytes into the same full disk MAX_TRIES times per file and then
# records every genome in <log>.fail; the right response is to stop the run (exit 74).
FATAL_ERRNO = frozenset((errno.ENOSPC, errno.EDQUOT, errno.EROFS))

# Only these files are mirrored; everything else NCBI publishes for a genome is
# ignored (not downloaded, not verified, not reported missing). The names a full mirror
# had and this one drops are listed -- and removed from genomes synced under the old
# regime -- by rm_files.py (REMOVE_EXACT / REMOVE_SUFFIXES), which asserts at import that
# its list and this whitelist are disjoint.
WANTED_EXACT = frozenset((
    "annotation_hashes.txt",
    "assembly_status.txt",
    MD5_MANIFEST
))

# Suffixes appended to the assembly name, e.g. GCF_036600855.1_ASM3660085v1 + _genomic.fna.gz.
# These are matched by building the EXACT filename, never with endswith(): the string
# "_genomic.fna.gz" is also a suffix of "_cds_from_genomic.fna.gz" and
# "_rna_from_genomic.fna.gz", so a suffix test would silently re-admit both.
WANTED_SUFFIXES = (
    "_ani_contam_ranges.tsv",
    "_ani_report.txt",
    "_assembly_report.txt",
    ASSEMBLY_STATS_EXT,
    "_fcs_report.txt",
    GENOMIC_FASTA_EXT,
    "_genomic.gbff.gz",
    "_genomic.gff.gz",
    "_wgsmaster.gbff.gz",
)

# NCBI's second checksum table. It is NOT in the whitelist above: it is never mirrored
# for its own sake, only written into a genome directory on the rare occasion it settled a
# disagreement there (see ChecksumFallback), so some genome directories hold one and most
# do not. Nothing treats its absence as a fault.
UNCOMPRESSED_MANIFEST = "uncompressed_checksums.txt"

# Wanted names the manifest is HONEST about: on the server iff listed. Measured on 110
# genomes (87/87 gff, 72/72 wgsmaster on disk == in manifest) and 24 live manifests. When
# unlisted they 404, and because a 404 leaves nothing on disk the generic "wanted but not
# listed" probe below re-requested them on EVERY warm run -- for every unannotated
# assembly (no gff; ~half of GenBank) and every non-WGS one (no wgsmaster; 17% of RefSeq).
# So for these two the manifest is the whole truth and they are never probed. Contrast
# _fcs_report.txt (omitted from 35/110 GenBank manifests while served) and the ANI files
# (never listed, always served), which the probe exists for. This is the one hardcoded
# exception to the "derived, not hardcoded" rule in sync_genome; STATS["probe_404"] is
# how you would notice if it ever stopped holding.
MANIFEST_IS_TRUTH_SUFFIXES = ("_genomic.gff.gz", "_wgsmaster.gbff.gz")

# --------------------------------------------------------------------------- tables

# <log>.fail and <log>.bad are themselves valid input (see read_assembly_summary), so
# both carry the assembly_summary column names that the reader looks for. Those are
# GENOME_COLUMNS, which the selection select_genomes writes also opens with, so the three
# tables this reads share one shape by construction.
# Each adds the one column that says what went wrong -- `failed_files` for a verification,
# which names the files, `reason` for a sync, which has no file to name.
BAD_HEADER = table_header(*GENOME_COLUMNS) + "\tfailed_files\n"
FAIL_HEADER = table_header(*GENOME_COLUMNS) + "\treason\n"

# How <log>.bad names a file that failed verification. The tag matters as much as the
# name: a file the manifest lists and the mirror does not have is repaired by any sync,
# while one whose bytes have rotted is trusted by the manifest fast path and needs the
# directory removed first (--retry <log>.bad --delete).
FAILED_MISSING = "missing"                       # listed in the manifest, not on disk
FAILED_MISMATCH = "mismatch"                     # present, md5 differs from the manifest
FAILED_UNREADABLE = "unreadable"                 # present, could not be hashed (OSError)
FAILED_SEP = ","                                 # no NCBI assembly filename holds one

# One row of GENOME_COLUMNS as this module holds it. The second field is `url` rather
# than ftp_path because the value is normalised before it is stored -- ftp:// rewritten
# to https, the trailing slash settled -- so it is no longer the column as NCBI wrote it.
Genome = collections.namedtuple(
    "Genome", "accession url version_status excluded_from_refseq")

# The second line of assembly_status.txt, "assembly anomaly=...", is drawn from NCBI's
# "Anomalous assemblies" vocabulary -- the seven terms below, per
# ncbi.nlm.nih.gov/assembly/help/anomnotrefseq/ (now JS-rendered; the 2019-09 Wayback
# snapshot is the readable copy). The summary's excluded_from_refseq column mixes these
# with ~20 plain exclusion reasons ("derived from metagenome", "fragmented assembly",
# "genome length too small") that NCBI does NOT write to the file, so only terms in this
# set are carried across. Verified 2026-09-11 against live files for every token present
# in the archaeal GenBank+RefSeq summaries (11 tokens, one genome each): byte-identical.
ASSEMBLY_ANOMALIES = frozenset((
    "chimeric",
    "contaminated",
    "hybrid",
    "misassembled",
    "mixed culture",
    "sequence duplications",
    "unverified source organism",
))

# Run log. This module does NOT configure logging: it logs into the GTDB Migration Tk
# log, set up by gtdb_migration_tk.biolib_lite.logger.logger_setup() before the command
# runs. "timestamp.ncbi_genome_sync" is a child of the 'timestamp' logger that toolkit
# configures, so records propagate to its handlers (the GTDB log file and the console)
# without this module owning, opening or closing any of them. The logging module is used
# for its thread safety: workers log concurrently from the pool.
#
# NOTHING PER-GENOME IS LOGGED. A production run covers ~2M genomes, so one line each --
# synced, failed, status changed, swept, bad -- would dominate the toolkit log and make it
# useless. Every one of those is logged at DEBUG, which this logger's own level discards
# before it reaches any handler, so the cost is a dropped call and nothing more.
#
# The information is not lost. What actually needs acting on is written to files that are
# meant to be machine-read and fed back in: <log>.fail lists every genome the sync could
# not complete, with its reason, and <log>.bad every genome that failed --verify. The
# console shows failures as they happen through Progress.write(). The log keeps the
# per-run story: argv, the genome count, a progress/throttling snapshot every 60 s, and
# the final summary -- bounded by runtime, not by list size.
#
# Set LOG_LEVEL to logging.DEBUG to get the per-genome detail back when debugging.
LOG = logging.getLogger("timestamp.ncbi_genome_sync")
LOG_LEVEL = logging.INFO
LOG.setLevel(LOG_LEVEL)

# extra=FILE_ONLY keeps a record out of the toolkit's console handler while still allowing
# it into the log file. Used for the periodic progress snapshot, which the bar already
# shows on screen, and carried on the per-genome DEBUG records so that raising LOG_LEVEL
# for debugging does not start scrolling the progress bar away.
FILE_ONLY = {"file_only": True}


class _FileOnlyFilter(logging.Filter):
    def filter(self, record):
        return not getattr(record, "file_only", False)


_CONSOLE_FILTER = _FileOnlyFilter()


def quiet_console_detail():
    """Keep FILE_ONLY records off the console, leaving them in the log file.

    Applied to the toolkit's console handlers only; a FileHandler keeps everything.
    Harmless for other commands -- nothing else sets the file_only flag.
    """
    for handler in logging.getLogger("timestamp").handlers:
        if isinstance(handler, logging.FileHandler):
            continue
        if _CONSOLE_FILTER not in handler.filters:
            handler.addFilter(_CONSOLE_FILTER)


def wanted_files(asm):
    """The exact set of filenames to mirror for assembly `asm`."""
    return WANTED_EXACT | frozenset(asm + suffix for suffix in WANTED_SUFFIXES)


_local = threading.local()
_print_lock = threading.Lock()

# Throttling telemetry. NCBI's limit is not documented, so the only way to know whether
# a given -j is too high is to watch what it sends back: 503/429 means we are pushing
# harder than it wants. Reported at the end of every run.
_stats_lock = threading.Lock()
STATS = {"requests": 0, "retry_status": {}, "conn_drops": 0, "bytes": 0,
         "breaker_trips": 0, "paused_s": 0.0,
         # Requests spent probing for an unlisted file that turned out not to exist. The
         # budget is one request per warm genome; this is the counter that shows when a
         # layout change at NCBI starts eating into it (see MANIFEST_IS_TRUTH_SUFFIXES).
         "probe_404": 0,
         # Files md5checksums.txt condemned and uncompressed_checksums.txt cleared, and
         # the genomes they belong to. NCBI's two tables disagree for a small minority of
         # genomes; this is how much, so the scale stays visible rather than silent.
         "fallback_files": 0, "fallback_genomes": 0}


def _bump(key, n=1):
    with _stats_lock:
        STATS[key] += n


def _bump_status(status):
    with _stats_lock:
        STATS["retry_status"][status] = STATS["retry_status"].get(status, 0) + 1


# --------------------------------------------------------------------------- rate control

class RateLimiter(object):
    """Token bucket shared by every worker thread: at most `rate` requests/second in the
    steady state, with a burst allowance of one second's worth.

    This is the knob that matters. NCBI throttles on SUSTAINED request rate over a
    multi-minute window (measured 2026-09-10: ~19-21 req/s clean for 10 min; ~29 req/s
    tripped after 6 min; ~42 req/s tripped inside a minute). A worker count cannot bound
    that, because a warm request (unchanged genome) completes ~10x faster than a cold one,
    so the same -j yields wildly different rates. Capping requests directly makes -j a pure
    bandwidth-parallelism setting and lets it stay high for cold downloads.
    """

    def __init__(self, rate):
        self.lock = threading.Lock()
        self.configure(rate)

    def configure(self, rate):
        """Set the rate; the bucket starts full so the first second is never throttled."""
        with self.lock:
            self.rate = float(rate)
            self.capacity = max(1.0, self.rate)
            self.tokens = self.capacity
            self.updated = time.monotonic()

    def acquire(self):
        if self.rate <= 0:
            return
        while True:
            with self.lock:
                now = time.monotonic()
                self.tokens = min(self.capacity, self.tokens + (now - self.updated) * self.rate)
                self.updated = now
                if self.tokens >= 1.0:
                    self.tokens -= 1.0
                    return
                wait = (1.0 - self.tokens) / self.rate
            time.sleep(wait)


class RunStopped(Exception):
    """Base for 'this genome did not fail, the RUN is ending': the worker records it as
    "stopped: ..." in <log>.fail so a re-run picks it up, and makes no further requests."""


class ThrottledOut(RunStopped):
    """NCBI is refusing every request and pauses are not clearing it. Stop the run."""


class Interrupted(RunStopped):
    """SIGINT/SIGTERM arrived, or another worker hit a fatal filesystem error."""


class StopFlag(object):
    """Cooperative shutdown. Set by the signal handler (or by a Fatal error); checked before
    every request, inside every download loop and by every long sleep, so a stop completes
    within about one request per worker instead of draining the whole queue.

    Without this a Ctrl-C reached only the main thread, and ThreadPoolExecutor.__exit__ then
    waited for every QUEUED genome to finish -- minutes if a breaker pause was in force --
    and SIGTERM (a scheduler, `timeout`) was not handled at all: an instant death that left
    temp files behind and no summary. Now both signals produce a complete <log>.fail, the
    STOPPED summary, and exit 128+signal (130 for SIGINT, 143 for SIGTERM). A second signal
    falls through to the default action and kills the process at once.
    """

    def __init__(self):
        self.event = threading.Event()
        self.signum = None
        self.reason = None

    def set(self, reason, signum=None):
        if not self.event.is_set():
            self.reason, self.signum = reason, signum
        self.event.set()

    def is_set(self):
        return self.event.is_set()

    def check(self):
        """Raise Interrupted if a stop has been requested."""
        if self.event.is_set():
            raise Interrupted(self.reason)

    def sleep(self, seconds):
        """time.sleep that ends early -- with Interrupted -- when a stop arrives."""
        if self.event.wait(seconds):
            raise Interrupted(self.reason)


STOP = StopFlag()


def install_signal_handlers():
    def on_signal(signum, frame):
        name = signal.Signals(signum).name
        STOP.set("interrupted by %s" % name, signum)
        signal.signal(signum, signal.SIG_DFL)    # second one kills outright
        sys.stderr.write("\n%s: finishing in-flight requests, recording the rest in "
                         "<log>.fail (again to kill)\n" % name)
    signal.signal(signal.SIGINT, on_signal)
    signal.signal(signal.SIGTERM, on_signal)


class CircuitBreaker(object):
    """On any throttling response, pause EVERY worker, not just the one that saw it.

    Per-request exponential backoff cannot recover from NCBI throttling: while one worker
    sleeps, the other eight keep the aggregate rate up and the penalty is renewed. A global
    pause starves the rate deliberately and immediately.

    Tripping on the FIRST 503 is correct here, not over-eager: 24,000 clean cold requests
    produced zero background 503s, so one 503 means everything is being refused -- and in
    every observed episode the other in-flight workers were refused in the same second.

    Escalation is keyed to what the data showed matters: whether the pause WORKED. A trip
    within `immediate` seconds of the previous pause ending means it did not, so the next
    pause doubles (cap `maximum`). A trip after a genuinely clean stretch resets to `base`.

    After `max_consecutive` pauses in a row that each failed this way, the host is in a
    penalty state that minutes will not clear (measured 2026-09-10: 60, 120, 240 and 300 s
    pauses were each answered by a 503 within 1-17 s of resuming, 67 min after the load that
    caused it). Crawling on -- one genome per 7 minutes -- only deepens it. The breaker then
    raises ThrottledOut and the run stops cleanly with exit 75 (EX_TEMPFAIL): completed
    genomes are on disk, in-flight ones are in <log>.fail, and a re-run later is cheap
    because everything done takes the fast path.
    """

    def __init__(self, base=60.0, maximum=300.0, immediate=30.0, max_consecutive=4):
        self.base = float(base)
        self.maximum = float(maximum)
        self.immediate = float(immediate)
        self.max_consecutive = int(max_consecutive)
        self.pause = self.base
        self.paused_until = 0.0
        self.trips = 0
        self.consecutive = 0                     # pauses in a row that did not work
        self.tripped_out = False
        self.lock = threading.Lock()
        self.notify = None                       # optional callable(str) for the terminal

    def wait(self):
        """Block while a pause is in force. Raises ThrottledOut once the breaker has given
        up, so a worker never sends another request."""
        while True:
            STOP.check()
            with self.lock:
                if self.tripped_out:
                    raise ThrottledOut("NCBI refusing requests; %s consecutive pauses "
                                       "failed to clear it"
                                       % format_count(self.consecutive))
                remaining = self.paused_until - time.monotonic()
            if remaining <= 0:
                return
            STOP.sleep(min(remaining, 5.0))

    def trip(self, status, retry_after=None):
        """Record a throttling response. Returns the pause length if this call opened a
        new pause, else None (a pause was already in force)."""
        with self.lock:
            now = time.monotonic()
            if self.tripped_out or now < self.paused_until:
                return None
            since_resume = now - self.paused_until      # <0 only before the first trip
            if self.trips and since_resume <= self.immediate:
                self.consecutive += 1                    # the last pause did not work
                self.pause = min(self.maximum, self.pause * 2)
            else:
                self.consecutive = 1                     # first trip, or after a clean run
                self.pause = self.base
            if retry_after:
                try:
                    self.pause = max(self.pause, min(float(retry_after), self.maximum))
                except ValueError:
                    pass
            self.trips += 1
            if self.consecutive >= self.max_consecutive:
                self.tripped_out = True
                pause, trips, cons = 0.0, self.trips, self.consecutive
            else:
                self.paused_until = now + self.pause
                pause, trips, cons = self.pause, self.trips, self.consecutive
        _bump("breaker_trips")
        _bump("paused_s", pause)
        if pause:
            message = ("circuit breaker: HTTP %d -> all workers paused %.0f s "
                       "(trip #%s, %s consecutive)"
                       % (status, pause, format_count(trips), format_count(cons)))
            LOG.warning(message)
        else:
            message = ("circuit breaker: HTTP %d after %s consecutive failed pauses -> "
                       "NCBI is refusing this host; stopping the run"
                       % (status, format_count(cons)))
            LOG.error(message)
        if self.notify:
            self.notify(message)
        return pause or None


LIMITER = RateLimiter(0)                         # configured from --rate in main()
BREAKER = CircuitBreaker()


# --------------------------------------------------------------------------- HTTP

class HttpError(Exception):
    pass


class Fatal(Exception):
    """The filesystem refused (disk full, quota, read-only). Nothing is half-installed --
    every write goes through a temp file -- but continuing would only re-download the same
    bytes into the same full disk, so the run stops with exit 74 (EX_IOERR)."""


def raise_if_fatal(exc):
    """Re-raise a filesystem-exhaustion OSError as Fatal, and tell every other worker to
    stop. Anything else returns so the caller can treat it as an ordinary failure."""
    if isinstance(exc, OSError) and exc.errno in FATAL_ERRNO:
        STOP.set("fatal: %s" % exc.strerror)
        raise Fatal("%s: %s" % (exc.strerror, exc.filename or exc))


def _connection():
    """One persistent HTTPS connection per worker thread."""
    conn = getattr(_local, "conn", None)
    if conn is None:
        conn = http.client.HTTPSConnection(HOST, timeout=120,
                                           context=ssl.create_default_context())
        _local.conn = conn
    return conn


def _drop_connection():
    conn = getattr(_local, "conn", None)
    if conn is not None:
        try:
            conn.close()
        except Exception:
            pass
        _local.conn = None


def http_get(path, sink=None, headers=None):
    """GET `path` with retries. Returns (status, payload, response_headers).

    Without `sink`, payload is the body as bytes. With `sink` (a writable binary file)
    the body is streamed into it and payload is the body's md5 hex digest instead --
    computed on the wire so the caller can verify a download without re-reading it.
    Retries on connection errors and on NCBI's throttling codes (RETRY_STATUS).
    """
    delay = 2.0
    last = None
    for attempt in range(MAX_TRIES):
        try:
            # A retry after a mid-stream failure must start the body over. Without this
            # the second attempt appends to the bytes the first one already wrote, so the
            # md5 can never match: the caller discards it and every remaining try is
            # wasted. Fails safe (the bad file is never installed) but stalls the worker.
            if sink is not None:
                sink.seek(0)
                sink.truncate()
            BREAKER.wait()                       # honour a global pause first
            LIMITER.acquire()                    # then take a request token
            conn = _connection()
            _bump("requests")
            conn.request("GET", path, headers=headers or {})
            resp = conn.getresponse()
            status = resp.status

            if status in RETRY_STATUS:
                resp.read()                      # drain so the connection stays usable
                last = "HTTP %d" % status
                _bump_status(status)
                retry_after = resp.getheader("Retry-After")
                if status in (429, 503):
                    # Throttled. Trip the GLOBAL breaker and go straight to the next
                    # attempt: BREAKER.wait() at the top of the loop is the sleep. No
                    # per-request backoff -- that is exactly what failed to recover.
                    BREAKER.trip(status, retry_after)
                    continue
                if retry_after:                  # 500/502/504: a server hiccup, back off alone
                    try:
                        delay = max(delay, min(float(retry_after), MAX_BACKOFF))
                    except ValueError:
                        pass
            elif sink is not None and status == 200:
                digest = hashlib.md5()
                while True:
                    STOP.check()                 # abort a big download within one chunk
                    buf = resp.read(CHUNK)
                    if not buf:
                        break
                    sink.write(buf)
                    digest.update(buf)
                    _bump("bytes", len(buf))
                return status, digest.hexdigest(), resp.getheaders()
            else:
                body = resp.read()
                _bump("bytes", len(body))        # manifests and unlisted files count too
                return status, body, resp.getheaders()

        except (OSError, http.client.HTTPException) as exc:
            # Deliberately NOT `except Exception`: that also swallowed programming errors
            # (AttributeError, TypeError), retried them MAX_TRIES times and re-raised them
            # as HttpError, hiding the real cause. Transport failures are OSError (which
            # covers ConnectionResetError, socket timeouts and ssl.SSLError) and
            # http.client.HTTPException (BadStatusLine, IncompleteRead); anything else is
            # a bug and should surface immediately. One OSError is not transport at all:
            # sink.write() on a full disk. That is Fatal, not a retry.
            raise_if_fatal(exc)
            _drop_connection()
            _bump("conn_drops")
            last = repr(exc)

        if attempt < MAX_TRIES - 1:              # pointless to sleep before giving up
            STOP.sleep(min(delay, MAX_BACKOFF))
            delay *= 1.7
    raise HttpError("%s after %d attempts: %s" % (path, MAX_TRIES, last))


def genome_relpath(url):
    """URL -> the mirror-relative path for that assembly, e.g. all/GCF/036/.../ASM...

    Both this and url_to_path() previously sliced by len(URL_PREFIX) without checking the
    prefix, so a line that was not an ftp.ncbi.nlm.nih.gov /genomes/ URL silently produced
    a WRONG directory rather than an error:
        ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/1/  -> genomes/l/GCF/1
        https://evil.example.com/x/../../../tmp/pwned/ -> escapes the mirror root
    The old shell pipeline rewrote these URLs to rsync form (ftp.ncbi...::), so a stray
    non-https line in a hand-edited list is a live possibility. Reject them loudly.
    """
    if not url.startswith(URL_PREFIX):
        raise BadInput("not an %s URL: %s" % (URL_PREFIX, url))
    rel = url[len(URL_PREFIX):].strip("/")
    if not rel:
        raise BadInput("no assembly path in URL: %s" % url)
    if os.path.isabs(rel) or ".." in rel.split("/"):
        raise BadInput("path traversal in URL: %s" % url)
    return rel


def url_to_path(url):
    """https://ftp.ncbi.nlm.nih.gov/genomes/all/... -> /genomes/all/..."""
    return "/genomes/" + genome_relpath(url) + "/"


# --------------------------------------------------------------------------- helpers

def parse_manifest(data, keep):
    """md5checksums.txt bytes -> [(md5, relative_path)] for the whitelisted files only.

    `keep` is the set from wanted_files(); entries outside it are dropped so they are
    neither downloaded nor treated as missing during verification.
    """
    out = []
    for checksum, name in read_md5_manifest(data.decode("utf-8", "replace").splitlines()):
        # Match the FULL manifest path, not its basename. Every wanted file lives at
        # the genome root, so a path with any directory component cannot be one -- which
        # excludes the _assembly_structure/, all_assembly_versions/ and representative/
        # subtrees structurally, with no directory blacklist to keep in sync. (Measured:
        # 1353 of 42323 entries across 3000 manifests are nested; none has a whitelisted
        # basename, so this is equivalent today and stays correct if that ever changes.)
        if name not in keep:
            continue

        out.append((checksum, name))
    return out


def temp_for(final_path):
    """A unique temp path beside `final_path`.

    os.getpid() alone was not enough: threads share a PID, so two workers handling the
    same assembly (a duplicated line in a list) would use the same temp name and corrupt
    each other. mkstemp is atomic and unique per call.
    """
    directory = os.path.dirname(final_path) or "."
    if not os.path.isdir(directory):
        os.makedirs(directory, DIR_MODE, exist_ok=True)
    fd, path = tempfile.mkstemp(dir=directory,
                                prefix="." + os.path.basename(final_path) + ".")
    os.close(fd)
    return path


def sweep_temps(genome_dir, keep):
    """Remove temp files a killed run left beside this genome's files. Returns the count.

    install() is atomic, so a SIGKILL or OOM mid-download leaves nothing half-installed --
    but it does leave the temp (".<name>.XXXXXX", a partial multi-MB .fna per worker) that
    nothing else ever removes. Only temps for names in `keep`, and only older than
    TEMP_STALE_S: a download in progress rewrites its temp's mtime every CHUNK, so a
    concurrent run's live file is never touched.
    """
    try:
        names = os.listdir(genome_dir)
    except OSError:
        return 0
    now = time.time()
    removed = 0
    for name in names:
        if not name.startswith("."):
            continue
        stem = name[1:].rsplit(".", 1)[0]       # ".<final>.XXXXXX" -> <final>
        if stem not in keep and stem not in (LAST_SYNCED, UNCOMPRESSED_MANIFEST):
            continue
        path = os.path.join(genome_dir, name)
        try:
            if now - os.path.getmtime(path) > TEMP_STALE_S:
                os.unlink(path)
                removed += 1
        except OSError:
            pass
    return removed


def install(tmp_path, final_path):
    os.chmod(tmp_path, FILE_MODE)                # mkstemp made it 0600
    dirname = os.path.dirname(final_path)
    if dirname and not os.path.isdir(dirname):
        os.makedirs(dirname, DIR_MODE, exist_ok=True)
    os.replace(tmp_path, final_path)


def set_mtime_from_headers(path, headers):
    for key, value in headers:
        if key.lower() == "last-modified":
            parsed = email.utils.parsedate_tz(value)
            if parsed:
                stamp = calendar.timegm(parsed[:6]) - (parsed[9] or 0)
                os.utime(path, (stamp, stamp))
            return


# Its columns, by name. NCBI has grown its tables before (crc32 and size are not in every
# copy of this one), so they are found by name and never by position -- the same rule
# ncbi_utils applies to the assembly summaries, for the same reason.
UNCOMPRESSED_NAME_COLUMN = "file"
UNCOMPRESSED_MD5_COLUMN = "md5sum"

# uncompressed_checksums.txt lists the UNCOMPRESSED form of every file, so a compressed
# one appears there under a name it does not have on disk -- <asm>_genomic.fna against
# <asm>_genomic.fna.gz -- and can never be matched. Asking about it is a request that
# cannot succeed, so a compressed file is never asked about at all: a genome whose only
# failures are archives spends nothing on a second opinion.
#
# This is an extension test, which the wanted-name matching elsewhere in this module
# deliberately is not: there, endswith() on "_genomic.fna.gz" also catches
# "_cds_from_genomic.fna.gz" and silently re-admits files the mirror does not carry. Here
# the question is only "is this name a compressed one", and the suffix is the whole of the
# answer. Everything compressed in the whitelist is a .gz, and a .gz added to
# WANTED_SUFFIXES later is excluded by this without a second list to keep in step.
COMPRESSED_EXT = ".gz"


def may_be_uncompressed_checksummed(name):
    """Can uncompressed_checksums.txt possibly carry an entry under this name?

    False for a compressed file, where the answer is structural rather than a matter of
    what NCBI happens to publish. Nothing else is hardcoded: which uncompressed files NCBI
    chooses to list is NCBI's to change, so any of them may be asked about.
    """
    return not name.endswith(COMPRESSED_EXT)
_MD5_HEX = re.compile(r"^[0-9a-fA-F]{32}$")
_warned_uncompressed = threading.Lock()
_warned_uncompressed_done = False


def read_uncompressed_checksums(lines):
    """{name: md5} from NCBI's uncompressed_checksums.txt, or None if it has no header.

    NOTHING about this file resembles md5checksums.txt, and assuming otherwise is a silent
    failure rather than a loud one: md5checksums.txt is `<md5>  ./<name>`, two spaces, no
    header; this is a TAB-separated table under `#file<TAB>md5sum<TAB>crc32<TAB>size`, with
    the name FIRST and the md5 SECOND. Reading it with the manifest reader yields an empty
    table, which looks exactly like "NCBI vouches for nothing here" -- so every second
    opinion came back negative and the genomes this exists for went on failing.

    Names keep NCBI's form less the leading ./, as read_md5_manifest does, so they compare
    directly against the manifest's. The table lists the UNCOMPRESSED form of every file,
    so for a compressed one the name here (<asm>_genomic.fna) is not the name on disk
    (<asm>_genomic.fna.gz) and no match is possible. That is the point rather than a
    limitation: this table can only ever vouch for files NCBI does not compress --
    <asm>_fcs_report.txt, the assembly and ANI reports -- and can never excuse a corrupt
    archive.

    Parameters
    ----------
    lines : iterable of str
        Lines of the table, from a decoded download.

    @return: {name: md5} for every row with both columns and a well-formed md5, or None
             when the header is missing, which is the one thing that must not be guessed.
    """

    columns = None
    out = {}
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith("#"):
            if columns is None:                  # the first comment line is the header
                names = [name.strip() for name in line[1:].split("\t")]
                columns = dict((name, i) for i, name in enumerate(names))
            continue
        if columns is None:
            return None                          # rows before a header: position is a guess
        name_at = columns.get(UNCOMPRESSED_NAME_COLUMN)
        md5_at = columns.get(UNCOMPRESSED_MD5_COLUMN)
        if name_at is None or md5_at is None:
            return None
        fields = line.split("\t")
        if len(fields) <= max(name_at, md5_at):
            continue                             # a short row is skipped, not fatal
        name, md5 = fields[name_at].strip(), fields[md5_at].strip()
        if not name or not _MD5_HEX.match(md5):
            continue
        if name.startswith("./"):
            name = name[2:]
        out[name] = md5.lower()
    return out if columns is not None else None


class ChecksumFallback(object):
    """NCBI's uncompressed_checksums.txt for ONE genome, asked for only when needed.

    NCBI publishes two checksum tables per assembly and they can disagree. Measured on the
    production mirror: ~1000 genomes whose md5checksums.txt entry for <asm>_fcs_report.txt
    is stale while uncompressed_checksums.txt is right. Believing md5checksums.txt alone
    condemns those genomes twice over -- the sync cannot install a file it just downloaded,
    because the bytes NCBI serves do not match the checksum NCBI publishes for them, so the
    genome fails every run and never enters the mirror; and a verify of one that did get
    mirrored reports rot that is not there. So a file that md5checksums.txt rejects gets a
    second opinion, and passes if the other table vouches for exactly the bytes on disk.

    It is a SECOND opinion, never a first: md5checksums.txt is asked first and settles
    every file it agrees with, and this is fetched lazily -- one request, on the first
    disagreement in the genome, never for a genome whose files all match. A clean mirror
    therefore verifies with no requests at all, as it always did.

    The table is not trusted from disk. A local copy vouching for local bytes is two halves
    of the same claim, and NCBI may have since corrected either table, so the answer comes
    from NCBI each time. The copy is written into the genome directory when it settles
    something, as the record of why those bytes were accepted.

    One instance per genome, used by the one worker handling it, so nothing here is shared
    or locked; only the STATS bump at the end is.
    """

    def __init__(self, url_path, genome_dir):
        self.url_path = url_path
        self.genome_dir = genome_dir
        self.entries = None                      # None until asked; {} when NCBI had none
        self.body = None
        self.vouched = []                        # names it cleared, for the log and counts

    def _fetch(self):
        """One request per genome, at most. A 404 is ordinary: most assemblies publish no
        uncompressed_checksums.txt, and that simply means no second opinion exists."""
        if self.entries is not None:
            return
        self.entries = {}
        try:
            status, body, _ = http_get(self.url_path + UNCOMPRESSED_MANIFEST)
        except HttpError as exc:
            LOG.debug("%s%s gave up: %s", self.url_path, UNCOMPRESSED_MANIFEST, exc,
                      extra=FILE_ONLY)
            return
        if status == 404:
            _bump("probe_404")
        if status != 200 or body is None:
            return
        entries = read_uncompressed_checksums(
            body.decode("utf-8", "replace").splitlines())
        if entries is None:
            # NCBI changed the format, and every second opinion from here on is worthless.
            # Once per run, at WARNING: this is the failure that hid the first time, and
            # per genome it would be a million identical lines.
            global _warned_uncompressed_done
            with _warned_uncompressed:
                if not _warned_uncompressed_done:
                    _warned_uncompressed_done = True
                    LOG.warning("warning: %s at NCBI has no '#%s ... %s' header; it cannot "
                                "be read, so files md5checksums.txt rejects will fail as "
                                "they did before it existed (first seen %s)",
                                UNCOMPRESSED_MANIFEST, UNCOMPRESSED_NAME_COLUMN,
                                UNCOMPRESSED_MD5_COLUMN, self.url_path)
            return
        self.body = body
        self.entries = entries

    def accepts(self, name, have_md5):
        """Does the other table vouch for exactly these bytes under exactly this name?

        A compressed file is refused before any request is made: the table lists the
        uncompressed form, under a different name, so no entry for it can exist (see
        may_be_uncompressed_checksummed). The name check that follows would refuse it
        anyway; this is what keeps a genome whose only failures are archives from spending
        a request to be told so.
        """
        if not may_be_uncompressed_checksummed(name):
            return False
        self._fetch()
        if self.entries.get(name) != have_md5:
            return False
        self._install()
        self.vouched.append(name)
        return True

    def _install(self):
        """Keep the table that settled it beside the files it settled. Failing to write it
        does not un-settle anything, so it is logged rather than raised -- except a full
        disk, which raise_if_fatal turns into the run-ending Fatal as everywhere else."""
        target = os.path.join(self.genome_dir, UNCOMPRESSED_MANIFEST)
        if self.body is None or os.path.exists(target):
            return
        tmp = temp_for(target)
        try:
            with open(tmp, "wb") as handle:
                handle.write(self.body)
            install(tmp, target)
        except OSError as exc:
            if os.path.exists(tmp):
                os.unlink(tmp)
            raise_if_fatal(exc)
            LOG.debug("could not write %s: %s", target, exc, extra=FILE_ONLY)

    def report(self, url):
        """Account for what it cleared. Nothing when it cleared nothing."""
        if not self.vouched:
            return
        _bump("fallback_files", len(self.vouched))
        _bump("fallback_genomes")
        LOG.debug("%s: %s accepted by %s against md5checksums.txt",
                  url, ", ".join(self.vouched), UNCOMPRESSED_MANIFEST, extra=FILE_ONLY)


def fetch_unlisted(url_path, genome_dir, name):
    """Timestamp-conditional fetch for a file with no published md5 (mimics wget -N).
    A 404 is normal -- not every genome has every one of these."""
    target = os.path.join(genome_dir, name)
    headers = {}
    if os.path.exists(target):
        headers["If-Modified-Since"] = email.utils.formatdate(
            os.path.getmtime(target), usegmt=True)
    try:
        status, body, resp_headers = http_get(url_path + name, headers=headers)
    except HttpError as exc:
        LOG.debug("unlisted %s%s gave up: %s", url_path, name, exc)
        return
    if status == 404:
        _bump("probe_404")
    if status != 200 or body is None:
        return                                   # 304 Not Modified, or 404
    tmp = temp_for(target)
    try:
        with open(tmp, "wb") as handle:
            handle.write(body)
        install(tmp, target)
    except BaseException:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise
    set_mtime_from_headers(target, resp_headers)


# --------------------------------------------------------------------------- sidecar files

def read_last_synced(genome_dir):
    """Epoch seconds from the genome's .last_synced marker, or None if absent/unreadable.

    The timestamp is the file's CONTENT, not its mtime: cp -p, rsync -a, touch and NFS
    clock skew all make mtimes untrustworthy across a mirror move, and this value decides
    whether NCBI is consulted at all."""
    try:
        with open(os.path.join(genome_dir, LAST_SYNCED)) as handle:
            return calendar.timegm(time.strptime(handle.read().strip(), LAST_SYNCED_FORMAT))
    except (OSError, ValueError, OverflowError):
        return None


def write_last_synced(genome_dir):
    """Stamp the genome as confirmed-complete now. Atomic, like every other install."""
    tmp = temp_for(os.path.join(genome_dir, LAST_SYNCED))
    with open(tmp, "w") as handle:
        handle.write(time.strftime(LAST_SYNCED_FORMAT, time.gmtime()) + "\n")
    install(tmp, os.path.join(genome_dir, LAST_SYNCED))


def clear_last_synced(genome_dir):
    """Withdraw the marker: verify found the genome wanting, so the next sync must look."""
    try:
        os.unlink(os.path.join(genome_dir, LAST_SYNCED))
    except OSError:
        pass


def fresh_enough(genome_dir, keep, max_age_s):
    """The number of whitelisted files a fresh, complete genome carries -- or None, meaning
    "ask NCBI". The marker saves exactly the network round-trip and nothing else: it must
    be younger than max_age_s, the manifest must be present, and EVERY whitelisted file it
    lists must exist. A missing file is always noticed (see FAST PATH) -- the old blob
    comparison once skipped straight past one, and a marker that did the same would be
    that bug again, with a two-week fuse."""
    stamp = read_last_synced(genome_dir)
    if stamp is None or time.time() - stamp > max_age_s:
        return None
    try:
        with open(os.path.join(genome_dir, MD5_MANIFEST), "rb") as handle:
            entries = parse_manifest(handle.read(), keep)
    except OSError:
        return None
    for _, name in entries:
        if not os.path.isfile(os.path.join(genome_dir, name)):
            return None
    return len(entries)


def render_status(genome):
    """assembly_status.txt as NCBI would write it, built from the summary row.

    NCBI generates the file from the same assembly record the summary comes from, so the
    summary can never be staler than the file -- and a warm re-sync used to spend a whole
    request per genome fetching it (half the fast-path cost). Format per
    genomes/all/README.txt: "status=<version_status>", then "assembly anomaly=<a, b>" only
    when there is one.

    Returns bytes, or None when the row has no version_status (a hand-cut table), in which
    case sync_genome falls back to fetching the file.
    """
    if not genome.version_status:
        return None
    text = "status=%s\n" % genome.version_status
    anomalies = [term.strip() for term in genome.excluded_from_refseq.split(";")
                 if term.strip() in ASSEMBLY_ANOMALIES]
    if anomalies:
        text += "assembly anomaly=%s\n" % ", ".join(anomalies)
    return text.encode()


def write_status(genome_dir, text):
    """Install assembly_status.txt if its content differs. Returns True only when an
    EXISTING file changed -- creating it on a cold sync is not news."""
    target = os.path.join(genome_dir, "assembly_status.txt")
    existed = os.path.isfile(target)
    if existed:
        with open(target, "rb") as handle:
            if handle.read() == text:
                return False
    tmp = temp_for(target)
    with open(tmp, "wb") as handle:
        handle.write(text)
    install(tmp, target)
    return existed


# --------------------------------------------------------------------------- sync

def sync_genome(url, root, full, status_text=None, max_age_s=0.0):
    """Sync ONE genome. `status_text` is assembly_status.txt rendered from the summary
    (render_status); None means fetch it from NCBI instead. `max_age_s` > 0 lets a genome
    whose .last_synced marker is younger than that skip NCBI entirely (never under --full).

    Returns (downloaded, verified, trusted, failures, unchanged, fresh):
      downloaded  files fetched -- missing locally, or present but not matching the manifest
      verified    files hashed locally and found to match the manifest
      trusted     files present whose manifest entry is identical to the entry in the
                  manifest the PREVIOUS sync installed, accepted without re-hashing
      failures    per-file reasons; the manifest is only installed when this is empty
      unchanged   nothing needed doing: every whitelisted manifest entry identical to the
                  installed manifest's and every file present (churn in files this mirror
                  does not carry is not a change)
      fresh       unchanged AND decided without a request: the marker was young and every
                  listed file was present, so the manifest was not even fetched

    The trust rule is what makes this cheap. A manifest is only ever installed after every
    whitelisted file it lists was verified against it, so "this entry is the same as last
    time and the file is still here" is as good as re-hashing -- except for local bit-rot,
    which is --full / --verify territory by documented policy. On an unchanged manifest
    every file is trusted (zero hashing, zero downloads); on a re-annotated RefSeq genome
    only the entries that moved -- typically the gff/gbff, not the multi-MB genomic.fna --
    are hashed. It also means a MISSING file is always noticed, whatever the manifest says.
    """
    url = url.rstrip("/") + "/"
    url_path = url_to_path(url)
    genome_dir = os.path.join(root, genome_relpath(url))
    asm = os.path.basename(genome_dir)
    manifest_path = os.path.join(genome_dir, MD5_MANIFEST)
    keep = wanted_files(asm)

    if max_age_s > 0 and not full:
        present = fresh_enough(genome_dir, keep, max_age_s)
        if present is not None:
            # Zero requests. The summary row may still have moved, and that is local:
            if status_text is not None and write_status(genome_dir, status_text):
                LOG.debug("status changed %s -> %s", url,
                          status_text.decode().strip().replace("\n", " | "),
                          extra=FILE_ONLY)
            return 0, 0, present, [], True, True

    status, body, _ = http_get(url_path + MD5_MANIFEST)
    if status != 200 or body is None:
        # Nothing is created on disk until the manifest is in hand: a genome removed
        # upstream (404) used to leave an empty directory behind for list_genomes to find.
        raise HttpError("manifest HTTP %d" % status)
    os.makedirs(genome_dir, DIR_MODE, exist_ok=True)
    swept = sweep_temps(genome_dir, keep)
    if swept:
        LOG.debug("swept %d stale temp file(s) %s", swept, url, extra=FILE_ONLY)

    old_body = None
    if os.path.exists(manifest_path):
        with open(manifest_path, "rb") as handle:
            old_body = handle.read()
    entries = parse_manifest(body, keep)
    old_entries = parse_manifest(old_body, keep) if old_body is not None else None

    # Two notions of "same". blob_same decides whether the manifest file is reinstalled.
    # wanted_same -- the WHITELISTED entries unchanged -- decides whether the genome is
    # "unchanged" and whether the unlisted files below are probed. NCBI regenerates a
    # manifest whenever any file in the directory changes, including the ~10 we do not
    # mirror (*_protein.faa.gz, *_feature_table.txt.gz ...), and probing on that would
    # spend 2-3 requests per genome on churn in files this mirror never sees.
    blob_same = (old_body == body)
    wanted_same = old_entries is not None and sorted(old_entries) == sorted(entries)

    # What the previous sync verified. Empty under --full so that everything is re-hashed.
    previous = {}
    if old_entries is not None and not full:
        previous = dict((name, md5) for md5, name in old_entries)

    downloaded = verified = trusted = 0
    failures = []
    fallback = ChecksumFallback(url_path, genome_dir)

    for want_md5, name in entries:
        target = os.path.join(genome_dir, name)
        if os.path.isfile(target):
            if previous.get(name) == want_md5:
                trusted += 1
                continue
            try:
                if file_md5(target) == want_md5:
                    verified += 1
                    continue
            except OSError as exc:
                failures.append("%s (%s)" % (name, exc))
                continue

        # missing, or present but stale -> fetch to a temp file, verify, then install
        tmp = temp_for(target)
        try:
            with open(tmp, "wb") as handle:
                st, got_md5, _ = http_get(url_path + name, sink=handle)
            if st != 200:
                raise HttpError("HTTP %d" % st)
            if got_md5 != want_md5 and not fallback.accepts(name, got_md5):
                # The bytes NCBI serves do not match the checksum NCBI publishes for them
                # in EITHER table. Not a transfer fault -- a download is retried on the
                # wire before it gets here -- so re-fetching would only fail the same way.
                raise HttpError("md5 mismatch")
            install(tmp, target)
            downloaded += 1
        except Exception as exc:
            if os.path.exists(tmp):
                os.unlink(tmp)
            if isinstance(exc, RunStopped):
                raise                            # the RUN is ending; not this file's fault
            raise_if_fatal(exc)
            failures.append("%s (%s)" % (name, exc))

    # Install the new manifest only once every file it lists is in place -- this is the
    # invariant the trust rule above depends on.
    if not failures and not blob_same:
        tmp = temp_for(manifest_path)
        with open(tmp, "wb") as handle:
            handle.write(body)
        install(tmp, manifest_path)

    # Derived, not hardcoded: a name is fetched here whenever it is wanted but is not
    # md5-enforced above. That covers two separate cases:
    #   * files NCBI serves but never checksums anywhere (assembly_status.txt, the ANI
    #     reports);
    #   * files listed for one division but not another -- GenBank manifests omit
    #     *_fcs_report.txt while still serving it (HTTP 200), so a hardcoded list
    #     silently skipped it for every GCA genome.
    # No md5 exists for any of them, so -N timestamping is the only check available.
    # A 404 is normal: not every genome has every file. The one exception to "derived,
    # not hardcoded" is MANIFEST_IS_TRUTH_SUFFIXES: names measured to 404 whenever they
    # are unlisted, so probing them only ever cost a request.
    not_fetched = set([MD5_MANIFEST])
    not_fetched.update(asm + suffix for suffix in MANIFEST_IS_TRUTH_SUFFIXES)
    if status_text is not None:
        not_fetched.add("assembly_status.txt")
        if write_status(genome_dir, status_text):
            LOG.debug("status changed %s -> %s", url,
                      status_text.decode().strip().replace("\n", " | "),
                      extra=FILE_ONLY)
    enforced = set(name for _, name in entries)
    for name in sorted(keep - enforced - not_fetched):
        # Without a summary status, assembly_status.txt is fetched every run: a genome can
        # be suppressed or replaced without any data file changing, so the manifest cannot
        # report it.
        if (name == "assembly_status.txt" or not wanted_same or full
                or not os.path.exists(os.path.join(genome_dir, name))):
            fetch_unlisted(url_path, genome_dir, name)

    fallback.report(url)
    unchanged = wanted_same and not full and downloaded == 0 and not failures
    if not failures:
        # LAST, after the manifest and the unlisted files: its existence means "everything
        # above completed", which is what lets --max-age trust it.
        write_last_synced(genome_dir)
    return downloaded, verified, trusted, failures, unchanged, False


# --------------------------------------------------------------------------- verify

def delete_genome(genome_dir):
    """Remove a genome directory so the next sync rebuilds it from scratch.

    The manifest goes FIRST, deliberately. The sync trusts an intact manifest, so if a
    deletion were interrupted after removing data files but before the manifest, the
    remaining state would look like a healthy genome with files missing -- which the
    existence check now catches, but removing the manifest first means nothing has to.
    Errors are logged and reported, never swallowed (this used to be ignore_errors=True).
    Returns True only on complete removal.
    """
    manifest = os.path.join(genome_dir, MD5_MANIFEST)
    try:
        if os.path.exists(manifest):
            os.unlink(manifest)
    except OSError as exc:
        LOG.debug("delete %s: could not remove manifest: %s", genome_dir, exc,
                  extra=FILE_ONLY)
        return False

    problems = []

    def on_error(func, path, exc_info):        # 3.11-: (func, path, exc_info)
        problems.append("%s: %s" % (path, exc_info[1]))

    def on_exc(func, path, exc):               # 3.12+: (func, path, exc)
        problems.append("%s: %s" % (path, exc))

    if sys.version_info >= (3, 12):
        shutil.rmtree(genome_dir, onexc=on_exc)
    else:
        shutil.rmtree(genome_dir, onerror=on_error)
    if problems:
        LOG.debug("delete %s: incomplete: %s", genome_dir, "; ".join(problems[:5]),
                  extra=FILE_ONLY)
        return False
    return True


def verify_genome(url, root, delete):
    """Check every whitelisted, checksummed file against the local manifest.

    (This is the job the verify_genomes.sh shell script used to do, before this tool
    replaced it -- that script is no longer part of the pipeline.)

    Note this can only attest to files md5checksums.txt actually covers: the four
    unlisted files have no published checksum, so "verified clean" means "every
    checksummed wanted file matches", not "every file is present and correct". A file
    md5checksums.txt rejects is put to NCBI's uncompressed_checksums.txt before being
    called bad, since the two tables disagree for a minority of genomes and that one is
    right -- see ChecksumFallback, and WHEN NCBI'S TWO CHECKSUM TABLES DISAGREE in the
    header. That is the only reason a verify makes a request, and only for a genome that
    is already failing on a file NCBI does not compress.

    A failure WITHDRAWS the genome's .last_synced marker (unless --delete removed the whole
    directory), so the next sync consults NCBI for this genome instead of skipping it as
    fresh for up to --max-age days. Whether that sync can repair the fault is the trust
    rule's business, not the marker's: a missing file, yes; rot in a trusted entry needs
    the directory gone first (--delete here, or --retry <log>.bad --delete) or --full.

    Returns (ok, detail), detail being <log>.bad's failed_files column: every file at
    fault, tagged missing:/mismatch:/unreadable: and comma-separated, or the reason when
    the fault is the genome's rather than a file's. Empty when ok.
    """
    url = url.rstrip("/") + "/"
    genome_dir = os.path.join(root, genome_relpath(url))
    fallback = ChecksumFallback(url_to_path(url), genome_dir)
    ok, detail = _verify_files(genome_dir, delete, fallback)
    fallback.report(url)
    if not ok and os.path.isdir(genome_dir):
        clear_last_synced(genome_dir)
    return ok, detail


def _tsv_safe(text):
    """Collapse whitespace, so a value cannot break the TSV it is written into.

    <log>.bad is itself --retry input, and a tab or a newline inside a field would make
    the next run read the wrong columns. Everything else is left as it stands: mangling a
    filename to fit a format would defeat the point of naming it.
    """
    return " ".join(str(text).split())


def _verify_files(genome_dir, delete, fallback):
    """md5-check one genome directory. Returns (ok, detail).

    `detail` is what <log>.bad's failed_files column carries and what the console line
    shows. EVERY manifest entry is checked before returning, so the column lists every
    file at fault rather than the first one found -- "which files are wrong, and how" is
    the question that file exists to answer, and stopping at the first would have the
    operator rebuild a genome to discover the next fault in it. A fault that is not a
    file's -- no directory, no manifest, a manifest that will not parse -- has no file to
    name, so the reason stands in the column in its place.

    This is why --delete acts at the END rather than at the first failure: the remaining
    entries have to still be on disk to be checked. A genome that fails is removed whole,
    once, after the check, exactly as before.
    """
    manifest_path = os.path.join(genome_dir, MD5_MANIFEST)

    if not os.path.isdir(genome_dir):
        return False, "missing directory"        # nothing to delete, nothing to name
    if not os.path.isfile(manifest_path):
        if delete:
            delete_genome(genome_dir)
        return False, "no md5checksums.txt"

    with open(manifest_path, "rb") as handle:
        raw = handle.read()
    entries = parse_manifest(raw, wanted_files(os.path.basename(genome_dir)))
    if not entries:
        # No WHITELISTED file is checksummed for this assembly. That is vacuously clean,
        # not a failure -- provided the manifest itself parsed. Reporting failure here
        # made sync and verify contradict each other: sync installs the manifest and
        # reports success, verify calls it "empty manifest", and --delete then removes a
        # correctly synced genome that the next sync recreates identically, looping
        # forever. (0 of 3000 production manifests hit this, so it is latent, not active.)
        if not any(read_md5_manifest(raw.decode("utf-8", "replace").splitlines())):
            if delete:
                delete_genome(genome_dir)
            return False, "unparseable md5checksums.txt"
        return True, ""

    failed = []                                  # manifest order, so a genome reads the same twice
    for want_md5, name in entries:
        target = os.path.join(genome_dir, name)
        if not os.path.isfile(target):
            failed.append("%s:%s" % (FAILED_MISSING, _tsv_safe(name)))
            continue
        try:
            have_md5 = file_md5(target)
            # md5checksums.txt first and last for everything it agrees with; only a file
            # it condemns is put to NCBI's other table (see ChecksumFallback)
            if have_md5 != want_md5 and not fallback.accepts(name, have_md5):
                failed.append("%s:%s" % (FAILED_MISMATCH, _tsv_safe(name)))
        except OSError as exc:
            # the errno belongs in the log, not in a column that is read back as input
            LOG.debug("verify %s: cannot hash %s: %s", genome_dir, name, exc,
                      extra=FILE_ONLY)
            failed.append("%s:%s" % (FAILED_UNREADABLE, _tsv_safe(name)))

    if not failed:
        return True, ""
    if delete:
        delete_genome(genome_dir)
    return False, FAILED_SEP.join(failed)


# --------------------------------------------------------------------------- driver

def bounded_map(fn, items, workers):
    """Run fn(item) for every item on `workers` threads, with a BOUNDED queue.

    ThreadPoolExecutor.map submits every item up front -- one Future per URL at ~1.6 KB
    each, so ~3 GB resident for the four production lists (1.88M genomes) before a single
    download starts. This keeps at most a few batches in flight, so memory scales with -j
    rather than with list length. Exceptions leaking out of fn are re-raised here.
    """
    inflight = threading.BoundedSemaphore(workers * 4)

    def run(item):
        try:
            fn(item)
        finally:
            inflight.release()

    with ThreadPoolExecutor(max_workers=workers) as pool:
        pending = []
        for item in items:
            if BREAKER.tripped_out or STOP.is_set():   # stop feeding; in-flight items drain
                break
            inflight.acquire()                   # blocks while the queue is full
            pending.append(pool.submit(run, item))
            if len(pending) >= workers * 8:      # prune, surfacing any leaked exception
                still = []
                for fut in pending:
                    if fut.done():
                        fut.result()
                    else:
                        still.append(fut)
                pending = still
        for fut in pending:
            fut.result()


def format_count(value):
    """A count as the log writes it, thousands separated: 1913482 -> "1,913,482".

    %-formatting has no comma flag, so counts are rendered here and passed to the
    logger as strings. A release is a few million genomes and a summary file a few
    million rows; a bare run of seven digits in a log line cannot be read at a
    glance, and two of them cannot be compared at all.

    Identifiers keep their digits: an HTTP status, an exit code, a PID, a line
    number in a table, and the values echoed back from the command line (-j, rate=,
    max_age=) are not quantities, and a comma in them would be wrong or unusable.
    """
    return "{:,}".format(value)


def format_amount(value, decimals=1):
    """A measured amount -- megabytes, seconds, a rate -- thousands separated.

    Same reasoning as format_count(); the decimals are kept because these are the
    numbers an operator compares between runs to see whether NCBI is slowing down.
    """
    return "{:,.{}f}".format(value, decimals)


def format_duration(seconds):
    """Seconds -> compact human duration, e.g. "2h 05m 30s", "7m 12s", "43s"."""
    seconds = max(0.0, seconds)
    if seconds < 10:
        return "%.1fs" % seconds
    seconds = int(round(seconds))
    hours, rem = divmod(seconds, 3600)
    minutes, secs = divmod(rem, 60)
    if hours:
        return "%dh %02dm %02ds" % (hours, minutes, secs)
    if minutes:
        return "%dm %02ds" % (minutes, secs)
    return "%ds" % secs


class Progress(object):
    """Progress for one pass (sync or verify): a tqdm bar on a terminal, periodic plain
    lines when piped, and a cumulative snapshot line to the log every 60 s. Counters are
    updated from worker threads via tick(); everything else runs on the main thread."""

    def __init__(self, label, total, silent):
        self.label = label
        self.total = total
        self.silent = silent
        self.done = 0
        self.failed = 0
        self.skipped = 0
        self.fresh = 0                           # skipped WITHOUT a request (.last_synced)
        self.started = time.time()
        self.lock = threading.Lock()
        # A bar is for a terminal. Piped to a log file the carriage returns become
        # unreadable soup, so fall back to periodic newline-terminated lines instead --
        # a multi-hour run is usually logged and tailed, and must still show progress.
        self.tty = sys.stderr.isatty()
        self.step = max(1, total // 20)          # ~20 progress lines over a piped run
        # Log snapshots are paced by TIME, not by count: `step` scales with list size, so
        # on a 1.9M-genome run it would emit one line per ~94k genomes (hours apart).
        self.snapshot_every = 60.0               # seconds
        self.last_snapshot = self.started
        self.bar = None
        if not silent and self.tty:
            self.bar = tqdm(total=total, desc=label, unit="genome", unit_scale=False,
                            dynamic_ncols=True, file=sys.stderr, leave=True,
                            smoothing=0.05)   # near-cumulative rate: NCBI throttling makes
                                              # an instantaneous ETA swing wildly

    def tick(self, failed=False, skipped=False, fresh=False):
        with self.lock:
            self.done += 1
            if failed:
                self.failed += 1
            if skipped:
                self.skipped += 1
            if fresh:
                self.fresh += 1
            now = time.time()
            if now - self.last_snapshot >= self.snapshot_every or self.done == self.total:
                self.last_snapshot = now
                self.snapshot()
            if self.silent:
                return
            if self.bar is not None:
                # tqdm rate-limits its own redraws (mininterval), so update() every tick
                # is cheap; postfix carries the counters the plain renderer used to show.
                self.bar.set_postfix(unchanged=self.skipped, 
                                     failed=self.failed,
                                     refresh=False)
                self.bar.update(1)
            elif self.done % (10 if self.tty else self.step) == 0 or self.done == self.total:
                self.render()

    def render(self):
        """Draw the progress line.

        On a terminal this forces a redraw of the tqdm bar, which supplies elapsed, rate
        and ETA itself. Piped to a log there is no bar (see __init__), so emit the same
        information as one line -- a logged run must still show progress and an ETA.
        """
        if self.silent:
            return
        if self.bar is not None:
            self.bar.set_postfix(unchanged=self.skipped, failed=self.failed, refresh=False)
            self.bar.refresh()
            return
        elapsed = time.time() - self.started
        rate = self.done / elapsed if elapsed > 0 else 0.0
        eta = format_duration((self.total - self.done) / rate) if rate > 0 else "?"
        sys.stderr.write(
            "%s  %s: %s/%s done, %s unchanged, %s failed [%s<%s, %s genome/s]    %s"
            % ("\r" if self.tty else "", self.label,
               format_count(self.done), format_count(self.total),
               format_count(self.skipped), format_count(self.failed),
               format_duration(elapsed), eta, format_amount(rate, 2),
               "" if self.tty else "\n"))
        sys.stderr.flush()

    def snapshot(self):
        """One cumulative progress line for the log.

        Written periodically so an interrupted run still shows how far it got and when
        throttling started -- the final HTTP summary only exists if the run completes.
        """
        elapsed = time.time() - self.started
        rate = self.done / elapsed if elapsed > 0 else 0.0
        eta = format_duration((self.total - self.done) / rate) if rate > 0 else "?"
        with _stats_lock:
            requests = STATS["requests"]
            megabytes = STATS["bytes"] / 1e6
            drops = STATS["conn_drops"]
            throttled = sum(STATS["retry_status"].get(code, 0) for code in (429, 503))
            trips, paused = STATS["breaker_trips"], STATS["paused_s"]
            probe_404 = STATS["probe_404"]
        req_rate = requests / elapsed if elapsed > 0 else 0.0
        LOG.info("Progress %s %s/%s (%.1f%%) unchanged=%s fresh=%s failed=%s "
                 "elapsed=%s eta=%s rate=%s/s requests=%s req/s=%s MB=%s "
                 "probe404=%s throttled=%s trips=%s paused=%ss drops=%s",
                 self.label, format_count(self.done), format_count(self.total),
                 100.0 * self.done / self.total if self.total else 100.0,
                 format_count(self.skipped), format_count(self.fresh),
                 format_count(self.failed), format_duration(elapsed), eta,
                 format_amount(rate, 2), format_count(requests), format_amount(req_rate),
                 format_amount(megabytes), format_count(probe_404),
                 format_count(throttled), format_count(trips), format_amount(paused, 0),
                 format_count(drops), extra=FILE_ONLY)

    def write(self, message):
        """Emit a message without corrupting an active bar."""
        if self.silent:
            return
        if self.bar is not None:
            self.bar.write(message, file=sys.stderr)
        else:
            sys.stderr.write("\r%s\n" % message)
            sys.stderr.flush()

    def finish(self):
        if self.silent:
            return
        if self.bar is not None:
            self.bar.close()
            self.bar = None
        else:
            self.render()
            if self.tty:
                sys.stderr.write("\n")


# --------------------------------------------------------------------------- input

# Columns this script cannot sync without: the accession names the genome, and
# ftp_path says where to fetch it from. Both <log>.fail and <log>.bad carry them too --
# five columns each, the fifth their own -- so those are read by this same code path.
SYNC_COLUMNS = GENOME_COLUMNS[:2]                # the genome, and where to fetch it from


def read_assembly_summary(path):
    """Read the genomes to sync from an NCBI assembly_summary.txt.

    Columns are located by name (see ncbi_utils) so the column count does not matter.
    version_status and excluded_from_refseq are read when present -- they are what
    assembly_status.txt is generated from (see render_status) -- and left empty otherwise,
    in which case that file is fetched from NCBI as it used to be.

    Rows whose ftp_path is "na" or empty are NOT an error: NCBI publishes them for
    suppressed and unreleased assemblies, and a whole-domain summary normally contains
    some. They are collected and returned for the caller to report, because logging is not
    configured yet when this runs. What counts as unusable is ncbi_utils.has_ftp_path(),
    the same test select_genomes applies, so a selection never promises a genome this
    reader then refuses.

    Every URL is validated here so a malformed table fails before any download rather than
    silently mirroring into the wrong directory (see genome_relpath).

    Returns (genomes, skipped): Genome(accession, url, version_status,
    excluded_from_refseq) records, and (lineno, accession, ftp_path) for every row with
    no usable ftp_path.
    """
    genomes, skipped, bad, seen, dupes = [], [], [], set(), 0
    for lineno, fields, columns in read_summary_rows(path, required=SYNC_COLUMNS):
        accession = summary_field(fields, columns, "assembly_accession")
        raw = summary_field(fields, columns, "ftp_path")
        # a few dozen distinct values between them, but one copy each per row would
        # be ~40% of the record's memory on a 1.9M-genome summary
        version_status = sys.intern(summary_field(fields, columns, "version_status"))
        excluded = sys.intern(summary_field(fields, columns, "excluded_from_refseq"))
        if not has_ftp_path(raw):
            skipped.append((lineno, accession or "?", raw or "(empty)"))
            continue
        url = raw.rstrip("/") + "/"
        if url.startswith(FTP_PREFIX):
            url = URL_PREFIX + url[len(FTP_PREFIX):]
        try:
            genome_relpath(url)
        except BadInput as exc:
            bad.append("  %s:%d: %s" % (path, lineno, exc))
            continue
        if url in seen:                          # same genome twice would sync twice
            dupes += 1
            continue
        seen.add(url)
        genomes.append(Genome(accession, url, version_status, excluded))
    if bad:
        raise BadInput("%s malformed ftp_path(s):\n%s"
                       % (format_count(len(bad)), "\n".join(bad[:10])))
    if dupes:
        LOG.warning("note: %s duplicate ftp_path(s) in %s ignored",
                    format_count(dupes), path)
    return genomes, skipped


# --------------------------------------------------------------------------- removal

# Default for --nfs-jobs: the threads the mirror walk and the removal each run on.
#
# Both are the same kind of work against the same server -- one metadata round trip after
# another, with the client otherwise idle -- so they take one setting rather than two.
# Walking costs one os.scandir() per directory; removing a genome costs ~18 operations (a
# scandir, ~12 unlinks, and the rmdirs up the triplets it empties). Measured on the
# production NFS server, removal ran at 113 genomes/s serially, 316 on 4 threads and 465 on
# 8; 16 threads reached 539, a further 16% for twice the load, and the walk was already
# documented as saturating by about four. Hence 8: past the knee, short of pointless.
#
# This is not --jobs. NCBI punishes too many requests, so --jobs is a safety mechanism;
# NFS just gets slower, and the ceiling here is the server's nfsd pool rather than a
# penalty. It is tunable because the mirror is shared (SHARED OPERATION): back off when
# others are working, raise it when they are not.
NFS_JOBS = 8

RM_HEADER = "#assembly_accession\tdirectory\n"


def accession_of(leaf):
    """Directory name -> versioned accession: GCA_000001405.28_GRCh38.p13 -> GCA_000001405.28.
    Cut at the first '_' after the GCA_/GCF_ prefix; the assembly name after it may itself
    contain underscores."""
    cut = leaf.find("_", 4)
    return leaf if cut < 0 else leaf[:cut]


def _subdirs(path):
    """Subdirectories of `path`, sorted; [] if it does not exist. Symlinks are not followed:
    a link into another tree must be neither deleted through nor counted as a genome."""
    try:
        with os.scandir(path) as entries:
            return sorted(e.path for e in entries if e.is_dir(follow_symlinks=False))
    except FileNotFoundError:
        return []


def mirror_genome_dirs(root, silent=False, workers=NFS_JOBS):
    """Every genome directory under <root>/all, as mirror-relative paths in the form
    genome_relpath() produces (all/GCF/000/006/805/GCF_000006805.1_ASM680v1), so the two can
    be compared as sets.

    Only the SHAPE of the layout is trusted -- all/<archive>/NNN/NNN/NNN/<leaf> -- never the
    names, so whatever NCBI adds under all/ is walked alike. Files at any level are ignored,
    and a leaf directory is a genome whatever it holds: a directory this script did not
    create, at that depth, is exactly the kind of stray the removal exists for.

    The second digit triplets are the unit of parallel work, as in list_genomes: tens of
    thousands of them, holding a comparable number of genomes each, walked on `workers`
    threads (--nfs-jobs).
    """
    second = []
    for archive in _subdirs(os.path.join(root, "all")):
        for first in _subdirs(archive):
            second.extend(_subdirs(first))

    def leaves(second_dir):
        return [os.path.relpath(leaf, root)
                for third in _subdirs(second_dir)
                for leaf in _subdirs(third)]

    found = set()
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for batch in tqdm(pool.map(leaves, second), total=len(second), unit="dir",
                          desc="Indexing mirror", file=sys.stderr, leave=False,
                          disable=silent or not sys.stderr.isatty()):
            found.update(batch)
    return found


def plan_mirror(root, genomes, silent=False, workers=NFS_JOBS):
    """What bringing `root` into line with `genomes` would do to its directories.

    Returns (to_remove, to_add, present): the directories on disk the selection does not
    name, sorted; and how many the selection names that are absent, and present. The match
    is on the exact directory the ftp_path maps to, so a renamed assembly or a version bump
    is removed and re-fetched rather than left beside its replacement. Only present/absent
    is decided here; whether a present genome needs updating is the sync's per-genome
    business, from the manifest.
    """
    expected = {genome_relpath(g.url) for g in genomes}
    on_disk = mirror_genome_dirs(root, silent, workers)
    return sorted(on_disk - expected), len(expected - on_disk), len(expected & on_disk)


def remove_genome_dir(root, rel):
    """Delete one genome directory, then the digit-triplet directories its removal emptied,
    stopping at <root>/all: the levels above hold the other genomes, and the sync recreates
    any it needs.

    The climb asks rmdir rather than listing first. rmdir already refuses a directory that
    is not empty, so the test and the act are one round trip instead of two -- four fewer
    per genome, on a pass that is nothing but round trips -- and there is no window between
    them. That is also what makes the climb safe to run in parallel: two genomes under one
    triplet race to delete it, and the loser is told ENOTEMPTY or ENOENT, which is the right
    answer either way, so it simply stops climbing.
    """
    top = os.path.normpath(os.path.join(root, "all"))
    path = os.path.join(root, rel)
    shutil.rmtree(path)
    parent = os.path.dirname(path)
    while os.path.normpath(parent) != top:
        try:
            os.rmdir(parent)
        except OSError:
            break                                # not empty, or another thread got there
        parent = os.path.dirname(parent)


def prune_mirror(root, to_remove, record_path, silent=False, workers=NFS_JOBS):
    """Remove every directory in `to_remove`, recording each in `record_path` BEFORE it is
    deleted, so a kill mid-way still leaves a list of what went. Stops between directories
    once STOP is set. Returns (removed, failed): a directory that will not delete (another
    member's 644 file, say -- see SHARED OPERATION) is logged and counted, not fatal, so
    the sync that follows still runs. The record is always written, header included, so a
    stale one from an earlier run cannot be mistaken for this run's.

    Removals run on `workers` threads (--nfs-jobs): the work is one round trip after another
    and the client is otherwise idle waiting for them. Rows therefore land in the record in
    completion order rather than sorted -- sort it if you want to diff two runs. The lock
    covers the write, so a row is still complete and still on disk before its directory is
    touched; and bounded_map keeps only a few batches in flight, which matters when the
    list is half a million long.
    """
    counts = {"removed": 0, "failed": 0}
    lock = threading.Lock()
    with open(record_path, "w") as record:
        record.write(RM_HEADER)
        record.flush()
        bar = tqdm(total=len(to_remove), unit="genome", desc="Removing", file=sys.stderr,
                   leave=False, disable=silent or not sys.stderr.isatty())

        def remove_one(rel):
            if STOP.is_set():
                return
            with lock:
                record.write("%s\t%s\n" % (accession_of(os.path.basename(rel)), rel))
                record.flush()
            try:
                remove_genome_dir(root, rel)
                outcome = "removed"
            except OSError as exc:
                outcome = "failed"
                LOG.error("could not remove %s: %s", rel, exc)
            with lock:
                counts[outcome] += 1
                bar.update()

        try:
            bounded_map(remove_one, to_remove, workers)
        finally:
            bar.close()
    return counts["removed"], counts["failed"]




# --------------------------------------------------------------------------- orchestration

def report_rows(headline, details, note):
    """Rows the run treats specially: each one to the log in full, a count and the first
    few names to stderr. Called from main() once logging is live -- the table is parsed
    before any file is created, including the log, so reporting cannot happen inside
    read_assembly_summary."""
    LOG.warning(headline)
    for line in details:
        LOG.warning(line)
    LOG.warning(note)


def first_names(names, limit=5):
    names = list(names)
    shown = ", ".join(names[:limit])
    if len(names) > limit:
        shown += ", ... (all %s listed in the log)" % format_count(len(names))
    return shown


def count_rows(path):
    """Data rows in one of our TSV outputs -- the '#' header line is not a failure."""
    with open(path) as handle:
        return sum(1 for line in handle
                   if line.strip() and not line.startswith("#"))





# The lock file keeps its original name though the command is now ncbi_genome_sync:
# it is mutual exclusion between processes, and two versions of this tool running on
# one mirror must contend for the SAME file. Renaming it would let a sync started
# before an upgrade and one started after hold different locks and run concurrently,
# which is the race lock_root exists to prevent.
def lock_root(root, name=".ncbi_sync.lock", what="sync"):
    """Hold <root>/<name> for the life of this process; None if another run has it.

    The restart hazard this closes: an operator re-launches a run they believe crashed, but
    it is alive in a 300 s breaker pause. Two syncs on one root are mostly harmless (every
    install is atomic) -- but a --verify --delete racing a sync is not, and both hammer
    NCBI from one host. flock is advisory and works on NFSv4; the file records who holds it.

    Each tool takes its OWN lock (rm_files.py: .rm_files.lock): a lock serialises runs of
    the same tool, not different tools, which are designed not to interfere.
    """
    os.makedirs(root, DIR_MODE, exist_ok=True)
    handle = open(os.path.join(root, name), "a+")
    try:
        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        handle.seek(0)
        holder = handle.read().strip()
        handle.close()
        LOG.error("error: another %s holds %s -- %s. If that run is really dead, "
                  "remove %s/%s", what, root, holder or "holder unknown", root, name)
        return None
    handle.seek(0)
    handle.truncate()
    handle.write("pid %d on %s since %s\n" % (os.getpid(), os.uname()[1],
                                              time.strftime("%Y-%m-%d %H:%M:%S")))
    handle.flush()
    return handle


def add_sync_arguments(parser):
    """Register the sync options on `parser`, which may be a standalone ArgumentParser or a
    subparser of a larger tool (gtdb_migration_tk ncbi_genome_sync). Kept separate from
    build_parser() so both entry points share one definition of the interface -- -l/--log
    included: the toolkit used to add --log itself, but the .fail/.bad/.rm outputs are placed
    by it, so the script has to know it too.

    Arguments go into the two groups the rest of the toolkit uses -- "required named
    arguments" first, then "options arguments". Both groups are returned.
    """
    required = parser.add_argument_group('required named arguments')
    # The two tables a run can start from are mutually exclusive because they say different
    # things about the mirror: the selection DEFINES it (what it does not list is removed);
    # a retry file only adds to it. dest stays "summary" for the selection -- the whole
    # module, and its tests, read args.summary -- and main() points it at the retry file
    # when that is what was given.
    table = required.add_mutually_exclusive_group(required=True)
    table.add_argument("-g", "--gtdb_selected_genomes", dest="summary", metavar="FILE",
                       help="The table written by select_genomes (gtdb_selected_genomes.tsv.gz). "
                            "Genomes are taken from its assembly_accession and ftp_path columns, "
                            "and it defines the mirror: every genome directory under --root that "
                            "it does not list is REMOVED before syncing. A subset, or a .fail/.bad "
                            "file, belongs to --retry.")
    table.add_argument("--retry", metavar="FILE",
                       help="A <log>.fail or <log>.bad from an earlier run (same columns): sync "
                            "only the genomes it lists and remove nothing it does not list. "
                            "With --delete each listed genome's directory is removed and "
                            "refetched, which is what a .bad needs.")
    required.add_argument("--root", required=True, metavar="DIR",
                          help="Mirror root directory, e.g. genomes")
    required.add_argument("-l", "--log", required=True, metavar="FILE",
                          help="Log file, appended across runs. It also names this run's "
                               "outputs, with .log stripped and written beside it: "
                               "sync.log -> sync.fail, sync.bad, sync.rm, sync.rm_dry_run, "
                               "sync.extra. Give each round its own log so one round "
                               "cannot overwrite another's.")

    optional = parser.add_argument_group('options arguments')
    optional.add_argument("--dry-run", action="store_true",
                          help="Report how many genome directories would be removed, added, and "
                               "checked for updates (already present), list the removals in "
                               "<log>.rm_dry_run, and exit. Nothing is downloaded, removed or "
                               "verified, and no lock is taken.")
    optional.add_argument("-j", "--jobs", type=int, default=8,
                        help="Parallel workers (default %(default)s). Bandwidth parallelism "
                             "only; --rate is what bounds the request rate.")
    optional.add_argument("--rate", type=float, default=20.0,
                        help="Global cap on requests/second across ALL workers (default "
                             "%(default)s; 0 disables). This, not -j, is what keeps NCBI from "
                             "throttling -- see RATE LIMITING in the file header.")
    optional.add_argument("--nfs-jobs", type=int, default=NFS_JOBS,
                        help="Threads for the mirror walk and the removal (default "
                             "%(default)s). These read and delete on the local mirror, not "
                             "at NCBI, so --rate does not apply: the work is one NFS round "
                             "trip after another and the threads overlap them. Measured "
                             "4.1x over serial at 8, and 16 bought a further 16%% -- lower "
                             "it when others are working on the mirror.")
    optional.add_argument("--verify-jobs", type=int, default=20,
                        help="Workers for the verification pass (default %(default)s). "
                             "Local md5, so not NCBI-limited -- see TUNING in the file "
                             "header for why more does not help. A file that fails its "
                             "md5checksums.txt entry does cost one request, for NCBI's "
                             "uncompressed_checksums.txt, which the two tables disagreeing "
                             "makes the deciding one; only uncompressed files are eligible, "
                             "and a clean mirror asks nothing.")
    optional.add_argument("--verify", action="store_true",
                        help="After syncing, md5-verify every genome against its manifest; "
                             "given the selection, also check the mirror holds nothing else "
                             "(directories it does not list go to <log>.extra)")
    optional.add_argument("--verify-only", action="store_true",
                        help="Skip the sync; only verify, both checks above. Removes nothing "
                             "without --delete.")
    optional.add_argument("--delete", action="store_true",
                        help="Remove genome directories rather than repair them. With "
                             "--retry: each listed genome's directory is removed "
                             "immediately before that genome is fetched, so one run "
                             "rebuilds the .fail/.bad from scratch. With "
                             "--verify/--verify-only: delete failing genome dirs, and "
                             "remove directories the selection does not list (recorded in "
                             "<log>.rm). Either way the manifest goes first, so a partial "
                             "delete can never look like a healthy genome. Illegal with a "
                             "bare --gtdb_selected_genomes sync, which would mean "
                             "re-downloading the whole mirror.")
    optional.add_argument("--max-age", type=float, default=14.0, metavar="DAYS",
                        help="A genome confirmed complete against NCBI within DAYS days "
                             "(its .last_synced marker) is skipped with NO request, after "
                             "a local check that every listed file is still present "
                             "(default %(default)s; 0 = always ask NCBI). Makes a crash "
                             "restart cheap -- see RESTART AND FRESHNESS in the header.")
    optional.add_argument("--full", action="store_true",
                        help="Bypass the manifest fast path AND --max-age: re-hash and "
                             "repair every file")
    optional.add_argument("--fail", help="Failure log (default: <log>.fail beside --log)")
    optional.add_argument("--bad", help="Verification failures (default: <log>.bad beside --log)")
    optional.add_argument("--silent", action="store_true", help="Suppress output")
    return required, optional


def build_parser():
    parser = argparse.ArgumentParser(
        description="Manifest-driven NCBI genome mirror sync with optional md5 verification.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    add_sync_arguments(parser)
    return parser


def validate_args(args):
    """Cross-checks argparse cannot express. Returns an error message, or None when the
    combination is legal. Risky-but-legal combinations get a warning on stderr here."""
    if args.jobs < 1 or args.verify_jobs < 1:
        return "--jobs and --verify-jobs must be >= 1"
    if args.nfs_jobs < 1:
        return "--nfs-jobs must be >= 1"
    if args.rate < 0:
        return "--rate must be >= 0"
    if args.max_age < 0:
        return "--max-age must be >= 0"
    if args.delete and not (args.verify or args.verify_only or args.retry is not None):
        return ("--delete only acts on --retry (rebuild the listed genomes) or during "
                "verification; add --retry, --verify or --verify-only")
    if args.rate <= 0 and args.jobs > 6:
        LOG.warning("warning: -j%d with --rate 0: nothing bounds the request rate. A warm "
                    "re-sync at -j9 lost 36 of the first 1,000 genomes to 503s.", args.jobs)
    if args.jobs > 9:
        LOG.warning("warning: -j%d: cold downloads at -j10 throttled even before rate "
                    "limiting existed; above 9 buys nothing measured.", args.jobs)
    return None


Outputs = collections.namedtuple("Outputs", "stem fail bad rm rm_dry_run extra")


def table_name(path):
    """The input table's name for the log and the progress bar, .gz then .txt/.tsv off.

    A LABEL, nothing more: it says which table a run is working from, in lines like
    "Syncing gtdb_selected_genomes: 1,204,331 genomes". What the run WRITES is named after
    --log (output_paths), so this never reaches a filename.
    """
    name = os.path.basename(path)
    if name.endswith(".gz"):
        name = name[:-len(".gz")]
    for suffix in (".txt", ".tsv"):
        if name.endswith(suffix):
            return name[:-len(suffix)]
    return name


def output_paths(args):
    """Where this run writes: every file named after --log, beside --log.

    The log names the run, so the run's outputs carry its name: -l r95/sync.log gives
    r95/sync.fail, .bad, .rm, .rm_dry_run and .extra. One stem for all five, because they
    are one run's account of itself and reading them means reading them together.

    Naming them after the INPUT, as this used to, tied a run's outputs to a table rather
    than to a run: two rounds against the same selection overwrote each other, while the
    .fail of one round and the .bad of another -- different inputs, same mirror, same log
    -- sat under unrelated names. Against that, the input's name had to survive a retry
    intact (x.bad -> x.bad.fail) so a retry could not overwrite the very list it was
    retrying. That is now the operator's to keep apart, and run() enforces it: give each
    round its own --log and it cannot arise, and if it does the run refuses (exit 2)
    rather than overwrite. .log is stripped; any other extension is kept, so a log
    deliberately named run.txt yields run.txt.fail rather than silently colliding with a
    run.log beside it.

    They go beside the log rather than into the working directory so that one directory,
    chosen by the operator, holds a mirror's whole run history (and the SHARED OPERATION
    chmod has one place to run). A bare --log name means that directory is the working
    directory, as before. --fail/--bad, when given, are used exactly as given."""
    stem = os.path.basename(args.log)
    if stem.endswith(".log") and stem != ".log":
        stem = stem[:-len(".log")]
    out_dir = os.path.dirname(args.log)

    def beside_log(name):
        return os.path.join(out_dir, name) if out_dir else name

    return Outputs(stem,
                   args.fail or beside_log(stem + ".fail"),
                   args.bad or beside_log(stem + ".bad"),
                   beside_log(stem + ".rm"),
                   beside_log(stem + ".rm_dry_run"),
                   beside_log(stem + ".extra"))


def standalone_logging(log_path):
    """Run as a script there is no logger_setup(), so give -l/--log its file, in the format
    the toolkit uses, and keep warnings on stderr as they always were. Under
    gtdb_migration_tk the toolkit has already configured this logger and nothing is added."""
    logger = logging.getLogger("timestamp")
    if logger.handlers:
        return
    logger.setLevel(logging.DEBUG)
    fmt = logging.Formatter("[%(asctime)s] %(levelname)s: %(message)s", "%Y-%m-%d %H:%M:%S")
    directory = os.path.dirname(log_path)
    if directory:
        os.makedirs(directory, DIR_MODE, exist_ok=True)
    to_file = logging.FileHandler(log_path)
    to_file.setFormatter(fmt)
    to_file.setLevel(logging.INFO)
    to_console = logging.StreamHandler(sys.stderr)
    to_console.setFormatter(fmt)
    to_console.setLevel(logging.WARNING)
    logger.addHandler(to_file)
    logger.addHandler(to_console)





class NCBIGenomeSync(object):
    """Bring a local NCBI genome mirror into line with a GTDB selection.

    One run per instance: the parsed arguments, the clock the run is reported against,
    and the output paths derived from them are the whole of its state. run() returns the
    process exit code rather than raising, so a caller can tell 75 (locked, retry later)
    and 130/143 (signalled) from a plain failure -- see EXIT CODES in the header.

    What lives here is the orchestration: what is removed, what is fetched, what is
    verified, in what order, and what the run was worth at the end. What deliberately
    does NOT live here is everything a worker calls -- sync_genome, verify_genome,
    http_get, prune_mirror and the rest stay module-level functions of a root and a URL.
    They run on -j threads at once and hold no run state, and making them methods would
    imply a shared object where there is none. The state those workers do share -- the
    rate limiter, the circuit breaker, the stop flag, the HTTP counters -- is
    module-level and process-wide for the same reason: there is one NCBI being asked,
    whatever is driving it.
    """

    def __init__(self, args):
        """Record the options this run works from.

        Parameters
        ----------
        args : argparse.Namespace
            Options as add_sync_arguments() defines them, already parsed.
        """

        self.args = args
        self.started = time.time()
        self.out = None                              # output_paths(), once args are valid

    def run(self):
        """Remove what the table does not list, sync what it does, verify if asked.

        @return: the process exit code (EXIT CODES in the header).
        """

        # Before ANY file is created, so nothing gets the inherited umask. 002, not 022:
        # the mirror is operated by a group (SHARED OPERATION in the header).
        os.umask(UMASK)
        if self.args.retry is not None:
            self.args.summary = self.args.retry      # one code path reads either table
        error = validate_args(self.args)
        if error:
            LOG.error("error: %s", error)
            return 2
        LIMITER.configure(self.args.rate)
        self.out = output_paths(self.args)
        for flag, path in (("--fail", self.out.fail), ("--bad", self.out.bad)):
            if os.path.realpath(path) == os.path.realpath(self.args.summary):
                # it is read fully before being truncated, so the run would work -- and a
                # crash would then have destroyed the only list of what needed retrying.
                # Outputs are named after --log, so this is what a retry of a round's own
                # .fail/.bad under that round's log looks like: give the retry its own log
                LOG.error("error: this run's %s output (%s) is the file it reads; it "
                          "would be overwritten. Give this run its own --log -- outputs "
                          "are named after it -- or place the file with %s",
                          flag, path, flag)
                return 2

        if self.args.dry_run:
            # Changes nothing under --root, so it takes no lock: an operator may well want
            # to see what the NEXT run would do while this one is still going. It also
            # creates nothing there, hence the check lock_root() would otherwise have done.
            if not os.path.isdir(self.args.root):
                LOG.error("error: --root %s is not a directory", self.args.root)
                return 2
            return self._run()

        try:
            lock = lock_root(self.args.root)         # held until the process exits
        except OSError as exc:
            LOG.error("error: cannot use --root %s: %s", self.args.root, exc)
            return 2
        if lock is None:
            return 75                                # EX_TEMPFAIL: try again when it is done

        try:
            return self._run()
        finally:
            lock.close()                             # release the flock even if run() is
                                                     # called again in-process (tests do);
                                                     # PyPy does not close dropped handles
                                                     # promptly, so this cannot be implicit

    def _prune(self, genomes):
        """Make --root hold exactly the selection's genomes, before anything is fetched.
        Returns 1 if any removal failed, else 0; the caller checks STOP itself."""
        to_remove, to_add, present = plan_mirror(self.args.root, genomes, self.args.silent,
                                                 self.args.nfs_jobs)
        LOG.info("Mirror: %s genome dir(s) to remove, %s to add, %s present (checked for "
                 "updates by the sync)", format_count(len(to_remove)),
                 format_count(to_add), format_count(present))
        removed, failed = prune_mirror(self.args.root, to_remove, self.out.rm,
                                       self.args.silent, self.args.nfs_jobs)
        LOG.info("Removed %s genome dir(s) not in the selection (listed in %s)%s",
                 format_count(removed), self.out.rm,
                 (", %s could not be removed" % format_count(failed)) if failed else "")
        return 1 if failed else 0

    def _dry_run(self, genomes):
        """Report what a real run would do to --root, write the removal list, change nothing."""
        to_remove, to_add, present = plan_mirror(self.args.root, genomes, self.args.silent,
                                                 self.args.nfs_jobs)
        if self.args.retry is not None:
            if self.args.delete:
                LOG.info("DRY RUN (--retry --delete): %s genome(s) to add, %s present and "
                         "removed before being fetched again. A retry removes nothing it "
                         "was not given.", format_count(to_add), format_count(present))
            else:
                LOG.info("DRY RUN (--retry): %s genome(s) to add, %s present and checked for "
                         "updates. A retry removes nothing.",
                         format_count(to_add), format_count(present))
            return 0
        with open(self.out.rm_dry_run, "w") as record:
            record.write("# DRY RUN: these directories would be removed; nothing was deleted\n")
            record.write(RM_HEADER)
            for rel in to_remove:
                record.write("%s\t%s\n" % (accession_of(os.path.basename(rel)), rel))
        LOG.info("DRY RUN: %s genome dir(s) would be removed (listed in %s), %s added, %s "
                 "present and checked for updates. Nothing was downloaded, removed or verified.",
                 format_count(len(to_remove)), self.out.rm_dry_run,
                 format_count(to_add), format_count(present))
        return 0

    def _sync(self, genomes):
        """Mirror every genome, recording per-genome failures in `fail_path`.

        Returns the exit code the sync alone would warrant: 75 (EX_TEMPFAIL) if the
        circuit breaker gave up on NCBI, 1 if any genome failed, else 0.
        """
        base, fail_path = table_name(self.args.summary), self.out.fail
        rc = 0
        # --delete on a retry: rebuild, not repair. Only reachable with --retry
        # (validate_args), so it can never wipe the whole mirror.
        rebuild = self.args.delete and self.args.retry is not None
        status = (f"Syncing {base}: {format_count(len(genomes))} genomes with "
                  f"{self.args.jobs}{'' if not self.args.full else ' (full re-hash)'} "
                  f"jobs in parallel{' (--delete: each directory removed before it is '
                  f'fetched)' if rebuild else ''}")
        LOG.info(status)
        prog = Progress("Downloading", len(genomes), self.args.silent)
        BREAKER.notify = prog.write
        fail_lock = threading.Lock()
        rebuild_lock = threading.Lock()
        rebuilt = 0
        max_age_s = self.args.max_age * 86400.0
        with open(fail_path, "w") as fail_log:
            fail_log.write(FAIL_HEADER)
            fail_log.flush()

            def record(genome, reason):
                # reasons come from repr(exc) and joined failure strings: a stray tab or
                # newline would corrupt a file that is itself input to the next run
                reason = " ".join(str(reason).split())
                with fail_lock:
                    fail_log.write("%s\t%s\t%s\t%s\t%s\n" % (
                        genome.accession, genome.url, genome.version_status,
                        genome.excluded_from_refseq, reason))
                    fail_log.flush()
                prog.tick(failed=True)

            def work(genome):
                nonlocal rebuilt
                url = genome.url
                if rebuild:
                    # Per genome, immediately before its fetch, so a stop leaves at most
                    # -j directories removed and not yet rebuilt -- and those are the ones
                    # recorded as "stopped:" in <log>.fail. A directory that will not go
                    # is a failure in its own right: syncing into it would rebuild the
                    # very state the operator asked to discard.
                    genome_dir = os.path.join(self.args.root, genome_relpath(url))
                    if os.path.isdir(genome_dir):
                        if not delete_genome(genome_dir):
                            prog.write("  FAILED %s: could not remove %s" % (url, genome_dir))
                            record(genome, "could not remove %s for rebuild" % genome_dir)
                            return
                        with rebuild_lock:
                            rebuilt += 1
                try:
                    downloaded, verified, trusted, failures, skipped, fresh = sync_genome(
                        url, self.args.root, self.args.full, render_status(genome), max_age_s)
                except RunStopped as exc:
                    # not this genome's fault; recorded so a re-run picks it up
                    record(genome, "stopped: %s" % exc)
                    return
                except Fatal:
                    raise                            # surfaces through bounded_map to main()
                except Exception as exc:
                    raise_if_fatal(exc)              # an OSError from install()/write_status()
                    LOG.debug("error %s %s", url, exc, extra=FILE_ONLY)
                    prog.write("  FAILED %s: %s" % (url, exc))
                    record(genome, exc)
                    return
                if failures:
                    LOG.debug("failed %s downloaded=%d verified=%d trusted=%d reasons=%s",
                              url, downloaded, verified, trusted, "; ".join(failures),
                              extra=FILE_ONLY)
                    prog.write("  FAILED %s: %s" % (url, "; ".join(failures)))
                    record(genome, "; ".join(failures))
                else:
                    if skipped:
                        # the common case on a re-sync: DEBUG so the log tracks work done
                        LOG.debug("fresh %s" if fresh else "unchanged %s", url)
                    else:
                        LOG.debug("synced %s downloaded=%d verified=%d trusted=%d",
                                  url, downloaded, verified, trusted, extra=FILE_ONLY)
                    prog.tick(skipped=skipped, fresh=fresh)

            try:
                bounded_map(work, genomes, self.args.jobs)
            finally:
                prog.finish()

        n_fail = count_rows(fail_path)
        # What a stop should be resumed with. --max-age is the whole point of re-running
        # the same table -- except under --delete, where every directory is removed before
        # it is fetched and there is no marker left to be fresh.
        resume = ("Re-run the SAME table with --delete; genomes already rebuilt are "
                  "fetched again, since --delete removes the marker --max-age reads."
                  if rebuild else
                  "Re-run the SAME summary -- completed genomes are skipped without a "
                  "request (--max-age).")
        if rebuild:
            LOG.info("Removed %s genome dir(s) listed in %s before fetching them (--delete)",
                     format_count(rebuilt), self.args.summary)
        if STOP.is_set():
            remaining = len(genomes) - prog.done
            msg = ("STOPPED (%s): %s/%s genomes done, %s in-flight recorded in %s, %s not "
                   "attempted. %s"
                   % (STOP.reason, format_count(prog.done - prog.failed),
                      format_count(len(genomes)), format_count(prog.failed),
                      fail_path, format_count(remaining), resume))
            LOG.warning(msg)
            rc = 128 + STOP.signum if STOP.signum else 74
        elif BREAKER.tripped_out:
            remaining = len(genomes) - prog.done
            msg = ("STOPPED: NCBI is refusing requests from this host and %s consecutive "
                   "pauses did not clear it. %s/%s genomes done, %s in-flight recorded in %s, "
                   "%s not attempted. Rest this host for hours, then: %s"
                   % (format_count(BREAKER.consecutive), format_count(prog.done - prog.failed),
                      format_count(len(genomes)), format_count(prog.failed),
                      fail_path, format_count(remaining), resume))
            LOG.error(msg)
            rc = 75                              # EX_TEMPFAIL: try again later
        LOG.info("Sync complete: %s genomes, %s unchanged (%s fresh, no request), %s failed, "
                 "elapsed=%s", format_count(len(genomes)), format_count(prog.skipped),
                 format_count(prog.fresh), format_count(n_fail),
                 format_duration(time.time() - prog.started))
        if n_fail:
            LOG.info("failures listed in %s", fail_path)
            if rc == 0:
                rc = 1

        return rc

    def _verify_nothing_else(self, genomes):
        """The other half of verifying against the selection: exactly these genomes, and
        nothing else. Every genome directory under --root that the selection does not name
        is listed in self.out.extra -- always written, header included, so a stale list cannot be
        mistaken for this run's -- and with --delete removed, each recorded in self.out.rm before
        it goes. Returns how many were found; the caller fails the verification on any."""
        extras, _, _ = plan_mirror(self.args.root, genomes, self.args.silent,
                                   self.args.nfs_jobs)
        with open(self.out.extra, "w") as record:
            record.write(RM_HEADER)
            for rel in extras:
                record.write("%s\t%s\n" % (accession_of(os.path.basename(rel)), rel))
        if not extras:
            return 0
        if self.args.delete:
            removed, failed = prune_mirror(self.args.root, extras, self.out.rm,
                                           self.args.silent, self.args.nfs_jobs)
            LOG.error("%s genome dir(s) not in the selection (listed in %s): %s removed%s",
                      format_count(len(extras)), self.out.extra, format_count(removed),
                      (", %s could not be removed" % format_count(failed)) if failed else "")
        else:
            LOG.error("%s genome dir(s) not in the selection -> %s. Re-run with --delete to "
                      "remove them, or run the selection sync.",
                      format_count(len(extras)), self.out.extra)
        return len(extras)

    def _verify(self, genomes):
        """Verify the mirror against the table. Every listed genome is md5-verified against its
        manifest, failures going to self.out.bad (and their directories deleted with --delete).
        Given the selection, the mirror must also hold NOTHING ELSE: directories it does not
        name go to self.out.extra and, with --delete, are removed (see verify_nothing_else). Given
        a retry file that second check does not apply -- a retry file is a subset.

        Returns 1 if any genome failed verification or the mirror held directories the
        selection does not list (whether or not --delete removed them), else 0.
        """
        base, bad_path = table_name(self.args.summary), self.out.bad
        vjobs = self.args.verify_jobs
        LOG.info("Verifying %s: %s genomes, -j%d%s", base, format_count(len(genomes)),
                 vjobs, " (delete on fail)" if self.args.delete else "")
        prog = Progress("verify " + base, len(genomes), self.args.silent)
        bad_lock = threading.Lock()
        with open(bad_path, "w") as bad_log:
            bad_log.write(BAD_HEADER)
            bad_log.flush()

            def check(genome):
                url = genome.url
                try:
                    ok, detail = verify_genome(url, self.args.root, self.args.delete)
                except RunStopped:
                    # Only a second opinion asks NCBI (ChecksumFallback), so only that can
                    # be stopped mid-verify. The genome is unjudged, not bad: recording it
                    # would put a healthy genome in <log>.bad, and the stop already says
                    # to re-run --verify-only on the same table.
                    prog.tick()
                    return
                if ok:
                    LOG.debug("verified %s", url, extra=FILE_ONLY)
                else:
                    LOG.debug("bad %s %s%s", url, detail,
                              " (deleted)" if self.args.delete else "", extra=FILE_ONLY)
                if not ok:
                    with bad_lock:
                        # failed_files last: every file at fault, tagged with which fault,
                        # so the row says what to do about the genome as well as that it
                        # is broken. The four columns before it are what --retry reads.
                        bad_log.write("%s\t%s\t%s\t%s\t%s\n" % (
                            genome.accession, url, genome.version_status,
                            genome.excluded_from_refseq, detail))
                        bad_log.flush()
                    with _print_lock:
                        prog.write("  BAD %s: %s" % (url, detail))
                prog.tick(failed=not ok)

            try:
                bounded_map(check, genomes, vjobs)
            finally:
                prog.finish()

        n_bad = count_rows(bad_path)
        LOG.info("Verify complete: %s genomes, %s failed, elapsed=%s",
                 format_count(len(genomes)), format_count(n_bad),
                 format_duration(time.time() - prog.started))

        # the selection says what the mirror holds, so verifying against it means checking
        # for what it does NOT list as well; not begun after a stop, which wants a re-run
        n_extra = 0
        if self.args.retry is None and not STOP.is_set():
            n_extra = self._verify_nothing_else(genomes)

        if STOP.is_set():
            msg = ("STOPPED (%s): verified %s/%s; %s lists failures among those. Re-run "
                   "--verify-only on the SAME table."
                   % (STOP.reason, format_count(prog.done), format_count(len(genomes)),
                      bad_path))
            LOG.warning(msg)
            return 128 + STOP.signum if STOP.signum else 74
        if n_bad:
            LOG.error("Failed verification: %s / %s  ->  %s",
                      format_count(n_bad), format_count(len(genomes)), bad_path)
            if self.args.delete:
                LOG.info("bad directories deleted; re-run sync on %s", bad_path)
            else:
                LOG.info("re-run with --retry %s --delete to remove and refetch them",
                         bad_path)
        if n_bad or n_extra:
            return 1

        LOG.info("All %s genomes verified clean%s.", format_count(len(genomes)),
                 "" if self.args.retry is not None else ", and the mirror holds nothing else")
        return 0

    def _run(self):
        retrying = self.args.retry is not None
        try:
            genomes, skipped = read_assembly_summary(self.args.summary)
        except (OSError, BadInput) as exc:
            LOG.error("error: %s", exc)
            return 2
        if not genomes and not skipped:
            if retrying:
                # a header and no rows: the normal end of a retry loop (an empty .fail), not an
                # error
                LOG.info("nothing to do: %s lists no genomes", self.args.summary)
                return 0
            # As the selection, an empty table says the mirror should hold nothing, and the
            # removal below would oblige. The likely cause is a .fail given to the wrong flag.
            LOG.error("error: %s lists no genomes. As --gtdb_selected_genomes it would remove "
                      "every genome under %s; an empty .fail belongs to --retry",
                      self.args.summary, self.args.root)
            return 2
        if not genomes:
            LOG.error("error: no genomes with a usable ftp_path in %s (%s rows, all na)",
                      self.args.summary, format_count(len(skipped)))
            return 2

        quiet_console_detail()
        LOG.info(shlex.join(sys.argv))
        stale = [g for g in genomes if g.version_status and g.version_status != "latest"]
        LOG.info("%s=%s genomes=%s no_ftp_path=%s not_latest=%s root=%s jobs=%d "
                 "rate=%.1f max_age=%gd verify_jobs=%d nfs_jobs=%d%s%s%s%s",
                 "retry" if retrying else "selection", self.args.summary,
                 format_count(len(genomes)), format_count(len(skipped)),
                 format_count(len(stale)), self.args.root, self.args.jobs, self.args.rate,
                 self.args.max_age,
                 self.args.verify_jobs, self.args.nfs_jobs,
                 " full" if self.args.full else "", " verify" if self.args.verify else "",
                 " verify_only" if self.args.verify_only else "", " dry_run" if self.args.dry_run else "")
        LOG.info("python=%s tqdm=%s host=%s", sys.version.split()[0], tqdm_version,
                 os.uname()[1])
        if skipped:
            report_rows(
                "%s row(s) in %s have no ftp_path and will not be synced"
                % (format_count(len(skipped)), self.args.summary),
                ["no ftp_path %s (%s:%d): ftp_path=%s" % (acc, self.args.summary, lineno, raw)
                 for lineno, acc, raw in skipped],
                "warning: %s genome(s) in %s have no ftp_path (na or empty) and were "
                "skipped: %s" % (format_count(len(skipped)), self.args.summary,
                                 first_names(acc for _, acc, _ in skipped)))
        if stale:
            # Policy: a replaced or suppressed genome is synced anyway. NCBI still serves the
            # directory, and pinned version lists (GTDB releases) want exactly that version.
            # The status lands in assembly_status.txt as before; this just makes it visible.
            report_rows(
                "%s genome(s) in %s are not the latest version; synced anyway"
                % (format_count(len(stale)), self.args.summary),
                ["not latest %s version_status=%s" % (g.accession, g.version_status)
                 for g in stale],
                "note: %s genome(s) in %s are replaced/suppressed; synced anyway with the "
                "status recorded in assembly_status.txt (see log)"
                % (format_count(len(stale)), self.args.summary))

        if self.args.dry_run:
            rtn_code = self._dry_run(genomes)
            self._finish(rtn_code)
            return rtn_code

        rtn_code = 0
        install_signal_handlers()
        try:
            if not retrying and not self.args.verify_only:
                # the selection defines the mirror: what it does not list goes first, so a
                # renamed or re-versioned assembly is gone before its replacement arrives
                rtn_code = self._prune(genomes)

            if STOP.is_set():
                # interrupted during removal: nothing is fetched, and the exit code says why
                LOG.warning("STOPPED (%s) during removal; what went is listed in %s and nothing "
                            "was fetched. Re-run the SAME selection.", STOP.reason, self.out.rm)
                rtn_code = 128 + STOP.signum if STOP.signum else 74
            elif not self.args.verify_only:
                # a failed removal (1) must neither mask a sync exit code nor be cleared by a
                # clean sync
                rtn_code = self._sync(genomes) or rtn_code

            if (self.args.verify or self.args.verify_only) and (BREAKER.tripped_out or STOP.is_set()):
                # Verifying a half-synced set is not informative -- every genome the stop
                # prevented would be "bad" -- and letting it set exit 1 would bury the 75/130
                # that tells a wrapper to wait rather than retry now.
                msg = ("Skipping --verify: the sync did not run to completion (exit %d)"
                       % rtn_code)
                LOG.warning(msg)
            elif self.args.verify or self.args.verify_only:
                # a sync exit code (1, 75, 130...) is never downgraded by a bad verify
                rtn_code = rtn_code or self._verify(genomes)
        except Fatal as exc:
            msg = ("FATAL: %s -- stopping. Nothing is half-installed; fix the filesystem and "
                   "re-run the SAME table." % exc)
            LOG.error(msg)
            rtn_code = 74                            # EX_IOERR

        self._finish(rtn_code)
        return rtn_code

    def _finish(self, rtn_code):
        """The two closing lines every run ends with, whatever it did: the HTTP account and the
        runtime with the exit code."""
        throttle = sum(STATS["retry_status"].get(c, 0) for c in (429, 503))
        # the status code is an identifier and keeps its digits; the count beside it
        # is a count
        detail = ", ".join("%d x%s" % (code, format_count(n))
                           for code, n in sorted(STATS["retry_status"].items()))
        http_summary = ("HTTP: %s requests (%s probe 404s), %s MB, %s connection drops, "
                        "%s throttled%s, %s breaker trips, %s s paused"
                        % (format_count(STATS["requests"]), format_count(STATS["probe_404"]),
                           format_amount(STATS["bytes"] / 1e6), format_count(STATS["conn_drops"]),
                           format_count(throttle), (" (%s)" % detail) if detail else "",
                           format_count(STATS["breaker_trips"]),
                           format_amount(STATS["paused_s"], 0)))
        runtime = "Total runtime: %s" % format_duration(time.time() - self.started)
        LOG.info(http_summary)
        if STATS["fallback_files"]:
            # NCBI's two checksum tables disagreeing is rare but real; a count each run is
            # how a change in that -- in either direction -- becomes visible
            LOG.info("Checksums: %s file(s) in %s genome(s) rejected by md5checksums.txt "
                     "and accepted by %s (the copy that settled it is in each genome "
                     "directory)", format_count(STATS["fallback_files"]),
                     format_count(STATS["fallback_genomes"]), UNCOMPRESSED_MANIFEST)
        LOG.info("%s (exit %d)", runtime, rtn_code)


def main(args=None):
    """Entry point for the standalone script and for gtdb_migration_tk's dispatcher alike.

    argv is parsed and the log opened only when called with no arguments, which is the
    standalone case; under gtdb_migration_tk both have already happened and the parsed
    options arrive as `args`. Kept as a function because the module is also a script.

    @return: the process exit code.
    """

    if args is None:
        args = build_parser().parse_args()
        os.umask(UMASK)                              # before the log file is created
        standalone_logging(args.log)

    return NCBIGenomeSync(args).run()


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        LOG.warning("interrupted by user")
        sys.stderr.write("\nInterrupted.\n")
        sys.exit(130)
