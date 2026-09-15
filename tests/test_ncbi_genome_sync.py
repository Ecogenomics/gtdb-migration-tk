#!/usr/bin/env python3
"""Offline unit tests for ncbi_genome_sync.py -- no network, no mirror.

Run with the interpreter that has tqdm (system python3 is 3.6 and cannot import the script):

    /opt/miniforge/bin/python3 -m unittest -v test_ncbi_genome_sync

Two contracts here would break silently in production and nowhere else, so they get the
most attention:

  * render_status() must be BYTE-IDENTICAL to the assembly_status.txt NCBI serves. The
    expected values below were captured from live files on 2026-09-11 for real genomes.
    If NCBI extends its anomaly vocabulary, this is what fails.
  * read_assembly_summary() must locate columns by name and refuse a table without the
    header -- the guard against a summary revision mirroring into the wrong directory.
"""

import logging
import os
import shutil
import tempfile
import types
import unittest

from gtdb_migration_tk import ncbi_genome_sync as N

P = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/"


def write(path, text):
    with open(path, "w") as handle:            # `with`: PyPy does not flush on refcount drop
        handle.write(text)
    return path


class TempDirCase(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix="ncbi_genome_sync_test.")

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def path(self, name):
        return os.path.join(self.dir, name)


# ------------------------------------------------------------------ assembly_status.txt

class RenderStatus(unittest.TestCase):
    def render(self, version_status, excluded):
        return N.render_status(N.Genome("GCA_x", "u/", version_status, excluded))

    def test_live_captures_are_byte_identical(self):
        # (version_status, excluded_from_refseq) -> live assembly_status.txt, 2026-09-11
        live = [
            ("latest", "na", b"status=latest\n"),                              # GCA_002287175.1
            ("latest", "derived from metagenome", b"status=latest\n"),         # GCA_015351695.1
            ("latest", "contaminated; derived from metagenome",
             b"status=latest\nassembly anomaly=contaminated\n"),               # GCA_009712615.1
            ("latest", "unverified source organism; derived from metagenome",
             b"status=latest\nassembly anomaly=unverified source organism\n"),  # GCA_048356915.1
        ]
        for version_status, excluded, expected in live:
            self.assertEqual(self.render(version_status, excluded), expected, excluded)

    def test_two_anomalies_comma_joined_in_summary_order(self):
        # the mirror holds files reading "assembly anomaly=contaminated, unverified source organism"
        self.assertEqual(
            self.render("latest", "contaminated; unverified source organism; derived from metagenome"),
            b"status=latest\nassembly anomaly=contaminated, unverified source organism\n")

    def test_exclusion_reasons_never_reach_the_file(self):
        for reason in ("fragmented assembly", "genome length too small",
                       "many frameshifted proteins", "annotation fails MAG completeness check",
                       "out of scope for the RefSeq project", "RefSeq annotation failed"):
            self.assertEqual(self.render("latest", reason), b"status=latest\n", reason)

    def test_every_ncbi_anomaly_term_is_carried(self):
        for term in N.ASSEMBLY_ANOMALIES:
            self.assertEqual(self.render("latest", term),
                             ("status=latest\nassembly anomaly=%s\n" % term).encode())

    def test_replaced_and_suppressed(self):
        self.assertEqual(self.render("replaced", "na"), b"status=replaced\n")
        self.assertEqual(self.render("suppressed", ""), b"status=suppressed\n")

    def test_no_version_status_means_fetch_instead(self):
        self.assertIsNone(self.render("", "na"))


class WriteStatus(TempDirCase):
    def test_creates_then_reports_only_real_changes(self):
        self.assertFalse(N.write_status(self.dir, b"status=latest\n"))       # new: not news
        self.assertFalse(N.write_status(self.dir, b"status=latest\n"))       # same: untouched
        self.assertTrue(N.write_status(self.dir, b"status=replaced\n"))      # changed
        with open(self.path("assembly_status.txt"), "rb") as handle:
            self.assertEqual(handle.read(), b"status=replaced\n")


# ------------------------------------------------------------------ input table

# Locating columns by name is ncbi_utils' job and is tested in tests/test_ncbi_utils.py.
# What is tested here is this script's own use of it: which columns it declares it cannot
# sync without, and what it makes of each row.

class ReadAssemblySummary(TempDirCase):
    HEADER = "#assembly_accession\tbioproject\tversion_status\tftp_path\texcluded_from_refseq\n"

    def read(self, text, name="t.txt"):
        return N.read_assembly_summary(write(self.path(name), text))

    def test_real_layout_with_comment_line(self):
        genomes, skipped = self.read(
            "##  See ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt\n"
            + self.HEADER
            + "GCA_1.1\tPRJNA1\tlatest\t%sGCA_1.1_A\tna\n" % P)
        self.assertEqual(skipped, [])
        self.assertEqual(genomes, [N.Genome("GCA_1.1", P + "GCA_1.1_A/", "latest", "na")])

    def test_na_and_empty_and_short_rows_are_skipped_not_fatal(self):
        genomes, skipped = self.read(
            self.HEADER
            + "GCA_2.1\tPRJNA2\tlatest\tna\tsuppressed\n"
            + "GCA_3.1\tPRJNA3\tlatest\t\t\n"
            + "GCA_5.1\tPRJNA5\n"
            + "GCA_1.1\tPRJNA1\tlatest\t%sGCA_1.1_A/\tna\n" % P)
        self.assertEqual([g.accession for g in genomes], ["GCA_1.1"])
        self.assertEqual([(a, r) for _, a, r in skipped],
                         [("GCA_2.1", "na"), ("GCA_3.1", "(empty)"), ("GCA_5.1", "(empty)")])

    def test_ftp_scheme_rewritten_and_duplicates_dropped(self):
        genomes, _ = self.read(
            self.HEADER
            + "GCA_4.1\tx\tlatest\tftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_4.1_D\tna\n"
            + "GCA_4.1\tx\tlatest\t%sGCA_4.1_D/\tna\n" % P)
        self.assertEqual([g.url for g in genomes], [P + "GCA_4.1_D/"])

    def test_no_header_is_rejected_by_name(self):
        with self.assertRaises(N.BadInput) as ctx:
            self.read("%sGCA_1.1_A/\n" % P, name="legacy.lst")
        self.assertIn("legacy.lst:1", str(ctx.exception))

    def test_header_without_an_accession_column_is_rejected(self):
        # every row would otherwise sync into a directory named for nothing
        with self.assertRaises(N.BadInput):
            self.read("#ftp_path\tbioproject\n%sGCA_1.1_A/\tPRJNA1\n" % P)

    def test_header_without_an_ftp_path_column_is_rejected(self):
        # every row would otherwise be skipped as having no directory to fetch
        with self.assertRaises(N.BadInput):
            self.read("#assembly_accession\tversion_status\nGCA_1.1\tlatest\n")

    def test_offsite_url_fails_before_any_download(self):
        with self.assertRaises(N.BadInput):
            self.read(self.HEADER + "GCA_1\tx\tlatest\thttps://evil.example.com/x/\tna\n")

    def test_fail_and_bad_outputs_round_trip(self):
        g, _ = self.read(N.FAIL_HEADER + "GCA_1.1\t%sGCA_1.1_A/\treplaced\tcontaminated\tHTTP 404\n" % P,
                         name="x.fail")
        self.assertEqual(g, [N.Genome("GCA_1.1", P + "GCA_1.1_A/", "replaced", "contaminated")])
        g, _ = self.read(N.BAD_HEADER + "GCA_2.1\t%sGCA_2.1_B/\tsuppressed\tna\n" % P, name="x.bad")
        self.assertEqual(g[0].version_status, "suppressed")

    def test_older_two_column_output_still_reads_and_falls_back_to_fetch(self):
        g, _ = self.read("#assembly_accession\tftp_path\nGCA_3.1\t%sGCA_3.1_C/\n" % P, name="old.bad")
        self.assertEqual(g[0].version_status, "")
        self.assertIsNone(N.render_status(g[0]))


# ------------------------------------------------------------------ manifest + URLs

class ParseManifest(unittest.TestCase):
    ASM = "GCA_000001405.29_GRCh38.p14"

    def test_whitelist_full_path_and_dot_slash(self):
        keep = N.wanted_files(self.ASM)
        data = ("\n".join([
            "0" * 32 + "  ./%s_genomic.fna.gz" % self.ASM,
            "1" * 32 + "  ./%s_cds_from_genomic.fna.gz" % self.ASM,       # suffix trap
            "2" * 32 + "  ./%s_protein.faa.gz" % self.ASM,                # not wanted
            "3" * 32 + "  ./%s_assembly_structure/Primary_Assembly/x" % self.ASM,  # nested
            "4" * 32 + "  ./annotation_hashes.txt",
            "not a manifest line",
        ]) + "\n").encode()
        self.assertEqual(N.parse_manifest(data, keep), [
            ("0" * 32, "%s_genomic.fna.gz" % self.ASM),
            ("4" * 32, "annotation_hashes.txt"),
        ])

    def test_manifest_is_truth_names_are_wanted(self):
        keep = N.wanted_files(self.ASM)
        for suffix in N.MANIFEST_IS_TRUTH_SUFFIXES:
            self.assertIn(self.ASM + suffix, keep)


class GenomeRelpath(unittest.TestCase):
    def test_good(self):
        self.assertEqual(N.genome_relpath(P + "GCA_1.1_A/"), "all/GCA/000/001/405/GCA_1.1_A")
        self.assertEqual(N.url_to_path(P + "GCA_1.1_A/"), "/genomes/all/GCA/000/001/405/GCA_1.1_A/")

    def test_rejects(self):
        for url in ("ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/1/",
                    "https://evil.example.com/x/../../../tmp/pwned/",
                    N.URL_PREFIX + "all/GCA/../../etc/",
                    N.URL_PREFIX):
            with self.assertRaises(N.BadInput, msg=url):
                N.genome_relpath(url)


# ------------------------------------------------------------------ main() helpers

class MainHelpers(TempDirCase):
    def args(self, **kw):
        base = dict(summary="assembly_summary_archaea_genbank.txt", retry=None, root="genomes",
                    log="gtdb_migration_tk.log", dry_run=False, jobs=8, verify_jobs=20,
                    nfs_jobs=8, rate=20.0, max_age=14.0, verify=False, verify_only=False,
                    delete=False, fail=None, bad=None)
        base.update(kw)
        return types.SimpleNamespace(**base)

    def test_output_paths_strip_txt_but_keep_fail_bad(self):
        out = N.output_paths(self.args())
        self.assertEqual((out.base, out.fail, out.bad),
                         ("assembly_summary_archaea_genbank",
                          "assembly_summary_archaea_genbank.fail",
                          "assembly_summary_archaea_genbank.bad"))
        out = N.output_paths(self.args(summary="/x/y/run.bad"))
        self.assertEqual((out.base, out.fail), ("run.bad", "run.bad.fail"))
        out = N.output_paths(self.args(fail="f", bad="b"))
        self.assertEqual((out.fail, out.bad), ("f", "b"))

    def test_outputs_are_placed_beside_the_log(self):
        # one directory, chosen by the operator, holds a mirror's whole run history
        out = N.output_paths(self.args(summary="/in/gtdb_selected_genomes.tsv.gz",
                                       log="/logs/sync.log"))
        self.assertEqual((out.fail, out.bad, out.rm, out.rm_dry_run, out.extra),
                         ("/logs/gtdb_selected_genomes.fail",
                          "/logs/gtdb_selected_genomes.bad",
                          "/logs/gtdb_selected_genomes.rm",
                          "/logs/gtdb_selected_genomes.rm_dry_run",
                          "/logs/gtdb_selected_genomes.extra"))

    def test_a_bare_log_name_means_the_working_directory(self):
        out = N.output_paths(self.args(summary="/in/x.tsv", log="sync.log"))
        self.assertEqual(out.fail, "x.fail")

    def test_explicit_fail_and_bad_ignore_the_log_directory(self):
        out = N.output_paths(self.args(fail="/elsewhere/f", bad="/elsewhere/b", log="/logs/s.log"))
        self.assertEqual((out.fail, out.bad), ("/elsewhere/f", "/elsewhere/b"))

    def test_output_paths_strip_gz(self):
        # GTDB stores the summaries compressed, and the table select_genomes
        # writes is gtdb_selected_genomes.tsv.gz; these outputs are plain text,
        # so they must not be named as though they were gzipped
        out = N.output_paths(self.args(summary="/x/y/gtdb_selected_genomes.tsv.gz"))
        self.assertEqual((out.base, out.fail, out.bad, out.rm),
                         ("gtdb_selected_genomes",
                          "gtdb_selected_genomes.fail",
                          "gtdb_selected_genomes.bad",
                          "gtdb_selected_genomes.rm"))

    def test_output_paths_strip_gz_before_the_retry_suffix(self):
        # a gzipped retry file still keeps its .bad, so its outputs sit beside it
        # rather than overwriting the run that produced it
        out = N.output_paths(self.args(summary="/x/y/run.bad.gz"))
        self.assertEqual((out.base, out.fail), ("run.bad", "run.bad.fail"))

    def test_output_paths_strip_a_bare_gz(self):
        self.assertEqual(N.output_paths(self.args(summary="/x/y/summary.gz")).base, "summary")

    def test_validate_args(self):
        self.assertIsNone(N.validate_args(self.args()))
        self.assertIn("--jobs", N.validate_args(self.args(jobs=0)))
        self.assertIn("--rate", N.validate_args(self.args(rate=-1)))
        self.assertIn("--nfs-jobs", N.validate_args(self.args(nfs_jobs=0)))
        self.assertIn("--max-age", N.validate_args(self.args(max_age=-1)))
        self.assertIsNone(N.validate_args(self.args(max_age=0)))
        self.assertIn("--delete", N.validate_args(self.args(delete=True)))
        self.assertIsNone(N.validate_args(self.args(delete=True, verify=True)))

    def test_count_rows_ignores_header_and_blank_lines(self):
        path = write(self.path("x.fail"), N.FAIL_HEADER + "a\tb\tc\td\te\n\n" + "f\tg\th\ti\tj\n")
        self.assertEqual(N.count_rows(path), 2)

    def test_first_names_truncates(self):
        self.assertEqual(N.first_names(["a", "b"]), "a, b")
        self.assertEqual(N.first_names(list("abcdefg")), "a, b, c, d, e, ... (all 7 listed in the log)")

    def test_format_duration(self):
        self.assertEqual(N.format_duration(3.14), "3.1s")
        self.assertEqual(N.format_duration(432), "7m 12s")
        self.assertEqual(N.format_duration(7530), "2h 05m 30s")


if __name__ == "__main__":
    unittest.main()


# ------------------------------------------------------------------ crash / restart safety

import errno
import sys
import time


class FatalClassification(unittest.TestCase):
    def setUp(self):
        N.STOP.event.clear(); N.STOP.reason = N.STOP.signum = None

    tearDown = setUp

    def test_disk_full_is_fatal_and_stops_the_run(self):
        for code in (errno.ENOSPC, errno.EDQUOT, errno.EROFS):
            N.STOP.event.clear()
            with self.assertRaises(N.Fatal):
                N.raise_if_fatal(OSError(code, os.strerror(code), "/mirror/x.tmp"))
            self.assertTrue(N.STOP.is_set(), code)

    def test_transport_errors_are_not_fatal(self):
        for exc in (OSError(errno.ECONNRESET, "reset"), TimeoutError(), ValueError("x")):
            N.raise_if_fatal(exc)                # returns quietly
        self.assertFalse(N.STOP.is_set())


class StopFlagBehaviour(unittest.TestCase):
    def setUp(self):
        N.STOP.event.clear(); N.STOP.reason = N.STOP.signum = None

    tearDown = setUp

    def test_check_and_sleep_raise_interrupted_which_is_a_run_stop(self):
        N.STOP.check()                           # not set: no-op
        N.STOP.set("interrupted by SIGTERM", 15)
        with self.assertRaises(N.RunStopped):
            N.STOP.check()
        t0 = time.time()
        with self.assertRaises(N.Interrupted):
            N.STOP.sleep(30)
        self.assertLess(time.time() - t0, 1.0)   # returned at once, did not sleep 30 s
        self.assertEqual(N.STOP.signum, 15)

    def test_first_reason_wins(self):
        N.STOP.set("first", 2); N.STOP.set("second", 15)
        self.assertEqual((N.STOP.reason, N.STOP.signum), ("first", 2))

    def test_breaker_wait_honours_stop(self):
        N.STOP.set("interrupted by SIGINT", 2)
        with self.assertRaises(N.Interrupted):
            N.BREAKER.wait()


class SweepTemps(TempDirCase):
    ASM = "GCA_000001405.29_GRCh38.p14"

    def touch(self, name, age_s):
        path = write(self.path(name), "partial")
        os.utime(path, (time.time() - age_s, time.time() - age_s))
        return path

    def test_only_old_temps_of_wanted_names_are_removed(self):
        keep = N.wanted_files(self.ASM)
        old = self.touch(".%s_genomic.fna.gz.k3j2h1" % self.ASM, 2 * N.TEMP_STALE_S)
        old_manifest = self.touch(".md5checksums.txt.a1b2c3", 2 * N.TEMP_STALE_S)
        fresh = self.touch(".%s_genomic.gbff.gz.z9y8x7" % self.ASM, 10)   # a live download
        other = self.touch(".not_ours.x1y2z3", 2 * N.TEMP_STALE_S)         # not a wanted name
        real = self.touch("%s_genomic.fna.gz" % self.ASM, 2 * N.TEMP_STALE_S)  # not a temp
        self.assertEqual(N.sweep_temps(self.dir, keep), 2)
        self.assertFalse(os.path.exists(old))
        self.assertFalse(os.path.exists(old_manifest))
        for path in (fresh, other, real):
            self.assertTrue(os.path.exists(path), path)

    def test_missing_dir_is_zero(self):
        self.assertEqual(N.sweep_temps(self.path("nope"), set()), 0)


class LockRoot(TempDirCase):
    def test_second_holder_is_refused_and_told_who(self):
        root = self.path("mirror")
        first = N.lock_root(root)
        self.assertIsNotNone(first)
        # the refusal is reported through the toolkit logger, not stderr
        records = []

        class _Capture(logging.Handler):
            def emit(self, record):
                records.append(record.getMessage())

        handler = _Capture()
        logger = logging.getLogger("timestamp")
        logger.addHandler(handler)
        try:
            self.assertIsNone(N.lock_root(root))
        finally:
            logger.removeHandler(handler)
        self.assertIn("pid %d" % os.getpid(), "\n".join(records))
        first.close()                            # released: can be taken again
        again = N.lock_root(root)
        self.assertIsNotNone(again)
        again.close()


class MainRunner(TempDirCase):
    """Drives main() the way the command line does; holds no tests of its own, so the
    classes built on it do not re-run one another's."""

    def run_main(self, table, *argv, retry=False):
        # The table is --gtdb_selected_genomes (the selection: removes what it does not
        # list) or --retry (syncs only what it lists). --log is required and places the
        # .fail/.bad/.rm outputs; under gtdb_migration_tk logger_setup() would create its
        # directory. ncbi_genome_sync reports through the toolkit logger, so capture that
        # too and return it alongside stderr -- the progress bar and signal handler still
        # use stderr.
        cwd = os.getcwd(); os.chdir(self.dir)
        os.makedirs(self.path("logs"), exist_ok=True)
        stderr, sys.stderr = sys.stderr, open(self.path("stderr"), "w")
        flag = "--retry" if retry else "--gtdb_selected_genomes"
        argv_saved, sys.argv = (sys.argv, ["ncbi_genome_sync.py", flag, table,
                                           "-l", os.path.join("logs", "sync.log")] + list(argv))
        records = []

        class _Capture(logging.Handler):
            def emit(self, record):
                records.append(record.getMessage())

        handler = _Capture()
        logger = logging.getLogger("timestamp")
        logger.addHandler(handler)
        try:
            rc = N.main()
        finally:
            logger.removeHandler(handler)
            sys.argv = argv_saved
            sys.stderr.close(); sys.stderr = stderr
            os.chdir(cwd)
        with open(self.path("stderr")) as handle:
            return rc, handle.read() + "\n".join(records)

class MainEndToEndOffline(MainRunner):
    """main() paths that need no network: argument guards and the empty-table exits."""

    def test_header_only_retry_is_nothing_to_do_exit_0(self):
        # the normal end of a retry loop: an empty .fail
        write(self.path("retry.fail"), N.FAIL_HEADER)
        rc, err = self.run_main("retry.fail", "--root", self.path("mirror"), retry=True)
        self.assertEqual(rc, 0)
        self.assertIn("nothing to do", err)

    def test_header_only_selection_is_refused_exit_2(self):
        # taken literally an empty selection says the mirror should hold nothing, and the
        # removal would oblige; the likely cause is an empty .fail given to the wrong flag
        write(self.path("retry.fail"), N.FAIL_HEADER)
        os.makedirs(self.path("mirror/all/GCA/000/001/405/GCA_000001405.1_A"))
        rc, err = self.run_main("retry.fail", "--root", self.path("mirror"))
        self.assertEqual(rc, 2)
        self.assertIn("--retry", err)
        self.assertTrue(os.path.isdir(self.path("mirror/all/GCA/000/001/405/GCA_000001405.1_A")))

    def test_rows_but_all_na_is_still_an_error(self):
        write(self.path("t.txt"), "#assembly_accession\tftp_path\nGCA_1.1\tna\n")
        rc, err = self.run_main("t.txt", "--root", self.path("mirror"))
        self.assertEqual(rc, 2)
        self.assertIn("all na", err)

    def test_output_equal_to_input_is_refused(self):
        write(self.path("x.fail"), N.FAIL_HEADER + "GCA_1.1\t%sGCA_1.1_A/\tlatest\tna\tr\n" % P)
        rc, err = self.run_main("x.fail", "--fail", "x.fail", "--root", self.path("mirror"))
        self.assertEqual(rc, 2)
        self.assertIn("would be overwritten", err)

    def test_main_releases_the_root_lock_on_return(self):
        # PyPy does not close a dropped handle promptly: without an explicit close the
        # flock outlived main() and a second run in the same process was refused with 75
        write(self.path("retry.fail"), N.FAIL_HEADER)
        self.assertEqual(self.run_main("retry.fail", "--root", self.path("mirror"), retry=True)[0], 0)
        self.assertEqual(self.run_main("retry.fail", "--root", self.path("mirror"), retry=True)[0], 0)

    def test_locked_root_exits_75_before_reading(self):
        write(self.path("t.txt"), N.FAIL_HEADER)
        held = N.lock_root(self.path("mirror"))
        try:
            rc, err = self.run_main("t.txt", "--root", self.path("mirror"))
        finally:
            held.close()
        self.assertEqual(rc, 75)
        self.assertIn("another sync holds", err)


# ------------------------------------------------------------------ .last_synced / --max-age

import calendar
import hashlib


def md5hex(data):
    return hashlib.md5(data).hexdigest()


class FakeNCBI(object):
    """Stand-in for http_get: serves a manifest and files from a dict, records every call.
    Lets sync_genome run cold, end to end, with no network."""

    def __init__(self, asm, files, manifest_extra=""):
        self.asm = asm
        self.files = files                       # name -> bytes
        self.manifest = ("".join("%s  ./%s\n" % (md5hex(b), n) for n, b in sorted(files.items()))
                         + manifest_extra).encode()
        self.calls = []

    def __call__(self, path, sink=None, headers=None):
        self.calls.append(path)
        name = path.rsplit("/", 1)[-1]
        if name == "md5checksums.txt":
            return 200, self.manifest, []
        if name in self.files:
            if sink is not None:
                sink.write(self.files[name])
                return 200, md5hex(self.files[name]), []
            return 200, self.files[name], []
        return 404, b"", []


class MarkerRoundTrip(TempDirCase):
    def test_write_read_clear(self):
        self.assertIsNone(N.read_last_synced(self.dir))
        before = time.time()
        N.write_last_synced(self.dir)
        stamp = N.read_last_synced(self.dir)
        self.assertLessEqual(abs(stamp - before), 2.0)
        with open(self.path(N.LAST_SYNCED)) as handle:
            self.assertRegex(handle.read(), r"^\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ\n$")
        N.clear_last_synced(self.dir)
        self.assertIsNone(N.read_last_synced(self.dir))
        N.clear_last_synced(self.dir)             # idempotent

    def test_garbage_or_mtime_only_is_not_a_timestamp(self):
        write(self.path(N.LAST_SYNCED), "yesterday\n")
        self.assertIsNone(N.read_last_synced(self.dir))
        write(self.path(N.LAST_SYNCED), "")      # e.g. a touch(1)
        self.assertIsNone(N.read_last_synced(self.dir))


class FreshEnough(TempDirCase):
    ASM = "GCA_000000001.1_ASM1v1"

    def setUp(self):
        super().setUp()
        self.keep = N.wanted_files(self.ASM)
        self.files = {self.ASM + "_genomic.fna.gz": b"ACGT", "annotation_hashes.txt": b"h"}
        for name, data in self.files.items():
            with open(self.path(name), "wb") as handle:
                handle.write(data)
        write(self.path("md5checksums.txt"),
              "".join("%s  ./%s\n" % (md5hex(b), n) for n, b in self.files.items()))
        N.write_last_synced(self.dir)

    def test_fresh_and_complete_reports_listed_count(self):
        self.assertEqual(N.fresh_enough(self.dir, self.keep, 14 * 86400), 2)

    def test_stale_marker_asks_ncbi(self):
        write(self.path(N.LAST_SYNCED),
              time.strftime(N.LAST_SYNCED_FORMAT, time.gmtime(time.time() - 15 * 86400)) + "\n")
        self.assertIsNone(N.fresh_enough(self.dir, self.keep, 14 * 86400))

    def test_missing_file_asks_ncbi_despite_fresh_marker(self):
        os.unlink(self.path(self.ASM + "_genomic.fna.gz"))
        self.assertIsNone(N.fresh_enough(self.dir, self.keep, 14 * 86400))

    def test_missing_manifest_asks_ncbi(self):
        os.unlink(self.path("md5checksums.txt"))
        self.assertIsNone(N.fresh_enough(self.dir, self.keep, 14 * 86400))

    def test_no_marker_asks_ncbi(self):
        N.clear_last_synced(self.dir)
        self.assertIsNone(N.fresh_enough(self.dir, self.keep, 14 * 86400))

    def test_verify_failure_withdraws_marker_success_keeps_it(self):
        url = P + self.ASM + "/"
        root = self.path("root")
        gdir = os.path.join(root, N.genome_relpath(url))
        shutil.copytree(self.dir, gdir)
        self.assertEqual(N.verify_genome(url, root, delete=False), (True, ""))
        self.assertIsNotNone(N.read_last_synced(gdir))          # clean: kept
        with open(os.path.join(gdir, self.ASM + "_genomic.fna.gz"), "wb") as handle:
            handle.write(b"rotten")
        ok, reason = N.verify_genome(url, root, delete=False)
        self.assertFalse(ok); self.assertIn("md5 mismatch", reason)
        self.assertIsNone(N.read_last_synced(gdir))             # withdrawn, dir kept
        self.assertTrue(os.path.isdir(gdir))


class SyncGenomeOffline(TempDirCase):
    """sync_genome end to end against FakeNCBI: the marker is written last and only on
    success, and a fresh genome makes no request."""

    ASM = "GCA_000000002.1_ASM2v1"
    URL = P + ASM + "/"

    def setUp(self):
        super().setUp()
        self.files = {self.ASM + "_genomic.fna.gz": b"ACGTACGT" * 100,
                      self.ASM + "_assembly_report.txt": b"report",
                      "annotation_hashes.txt": b"hashes"}
        self.fake = FakeNCBI(self.ASM, self.files)
        self.real_http_get = N.http_get
        N.http_get = self.fake
        self.gdir = os.path.join(self.dir, N.genome_relpath(self.URL))

    def tearDown(self):
        N.http_get = self.real_http_get
        super().tearDown()

    def sync(self, **kw):
        return N.sync_genome(self.URL, self.dir, kw.pop("full", False),
                             b"status=latest\n", kw.pop("max_age_s", 0.0))

    def test_cold_sync_writes_marker_last_and_only_on_success(self):
        downloaded, verified, trusted, failures, unchanged, fresh = self.sync()
        self.assertEqual((downloaded, failures, unchanged, fresh), (3, [], False, False))
        for name in self.files:
            self.assertTrue(os.path.isfile(os.path.join(self.gdir, name)), name)
        manifest = os.path.join(self.gdir, "md5checksums.txt")
        marker = os.path.join(self.gdir, N.LAST_SYNCED)
        self.assertTrue(os.path.isfile(manifest))
        self.assertGreaterEqual(os.stat(marker).st_mtime_ns, os.stat(manifest).st_mtime_ns)
        self.assertIsNotNone(N.read_last_synced(self.gdir))

    def test_failed_download_leaves_neither_manifest_nor_marker(self):
        self.fake.manifest = self.fake.manifest.replace(
            md5hex(self.files["annotation_hashes.txt"]).encode(), b"0" * 32)   # md5 will mismatch
        *_, failures, unchanged, fresh = self.sync()
        self.assertEqual(len(failures), 1)
        self.assertIn("md5 mismatch", failures[0])
        self.assertFalse(os.path.exists(os.path.join(self.gdir, "md5checksums.txt")))
        self.assertFalse(os.path.exists(os.path.join(self.gdir, N.LAST_SYNCED)))
        self.assertTrue(os.path.isfile(os.path.join(self.gdir, self.ASM + "_genomic.fna.gz")))

    def test_fresh_genome_makes_no_request_but_still_applies_status(self):
        self.sync()                                              # cold: creates everything
        calls_before = len(self.fake.calls)
        result = N.sync_genome(self.URL, self.dir, False, b"status=replaced\n", 14 * 86400)
        self.assertEqual(result, (0, 0, 3, [], True, True))
        self.assertEqual(len(self.fake.calls), calls_before)      # not one request
        with open(os.path.join(self.gdir, "assembly_status.txt"), "rb") as handle:
            self.assertEqual(handle.read(), b"status=replaced\n")   # summary row applied

    def test_max_age_zero_and_full_both_ask_ncbi(self):
        self.sync()
        n = len(self.fake.calls)
        self.assertFalse(self.sync(max_age_s=0.0)[5]);        self.assertGreater(len(self.fake.calls), n)
        n = len(self.fake.calls)
        self.assertFalse(self.sync(full=True, max_age_s=14 * 86400)[5]); self.assertGreater(len(self.fake.calls), n)

    def test_missing_file_inside_window_is_repaired_not_skipped(self):
        self.sync()
        os.unlink(os.path.join(self.gdir, self.ASM + "_assembly_report.txt"))
        downloaded, *_, fresh = self.sync(max_age_s=14 * 86400)
        self.assertEqual((downloaded, fresh), (1, False))
        self.assertTrue(os.path.isfile(os.path.join(self.gdir, self.ASM + "_assembly_report.txt")))


class GroupOperablePermissions(SyncGenomeOffline):
    """A cold sync under the script's umask leaves 2775 directories and 664 files, and a
    setgid root hands the bit down to every directory the sync creates."""

    def setUp(self):
        super().setUp()
        self.saved_umask = os.umask(N.UMASK)     # main() does this; sync_genome does not
        os.chmod(self.dir, 0o2775)               # like the real mirror root

    def tearDown(self):
        os.umask(self.saved_umask)
        super().tearDown()

    def test_modes_and_setgid_inheritance(self):
        self.sync()
        import stat as st
        path = self.dir
        for part in N.genome_relpath(self.URL).split("/"):
            path = os.path.join(path, part)
            mode = os.stat(path).st_mode
            self.assertEqual(st.S_IMODE(mode) & 0o777, N.DIR_MODE, path)
            self.assertTrue(mode & st.S_ISGID, "setgid not inherited by %s" % path)
        for name in list(self.files) + ["md5checksums.txt", "assembly_status.txt", N.LAST_SYNCED]:
            f = os.path.join(self.gdir, name)
            self.assertEqual(st.S_IMODE(os.stat(f).st_mode), N.FILE_MODE, name)
            self.assertFalse(os.stat(f).st_mode & 0o111, "execute bit on %s" % name)

    def test_constants_are_group_writable_and_no_file_exec(self):
        self.assertEqual(N.DIR_MODE & 0o070, 0o070)
        self.assertEqual(N.FILE_MODE & 0o060, 0o060)
        self.assertEqual(N.FILE_MODE & 0o111, 0)
        self.assertEqual(N.UMASK & 0o020, 0)        # umask must not strip group write


# ------------------------------------------------------------------ removal: the selection defines the mirror

GCA_PREFIX = N.URL_PREFIX + "all/GCA/000/001/405/"


def selection(*asms, header=N.BAD_HEADER):
    """A selected-genomes table naming the given assembly directory names."""
    return header + "".join("%s\t%s%s/\tlatest\tna\n" % (N.accession_of(asm), GCA_PREFIX, asm)
                            for asm in asms)


class AccessionOf(unittest.TestCase):
    def test_cuts_at_the_first_underscore_after_the_prefix(self):
        self.assertEqual(N.accession_of("GCA_000001405.28_GRCh38.p13"), "GCA_000001405.28")
        self.assertEqual(N.accession_of("GCF_000006805.1_ASM680v1"), "GCF_000006805.1")

    def test_assembly_names_with_underscores_do_not_confuse_it(self):
        self.assertEqual(N.accession_of("GCA_000000001.1_my_odd_name_v2"), "GCA_000000001.1")

    def test_a_name_without_an_assembly_part_is_returned_whole(self):
        self.assertEqual(N.accession_of("GCA_000000001.1"), "GCA_000000001.1")


class MirrorWalk(TempDirCase):
    """mirror_genome_dirs() trusts the shape of the layout and nothing else."""

    def leaf(self, rel):
        path = self.path(rel)
        os.makedirs(path)
        write(os.path.join(path, "md5checksums.txt"), "")
        return path

    def test_finds_every_leaf_under_every_archive(self):
        self.leaf("all/GCA/000/001/405/GCA_000001405.1_A")
        self.leaf("all/GCA/000/001/405/GCA_000001405.2_B")
        self.leaf("all/GCF/000/006/805/GCF_000006805.1_ASM680v1")
        self.assertEqual(N.mirror_genome_dirs(self.dir),
                         {"all/GCA/000/001/405/GCA_000001405.1_A",
                          "all/GCA/000/001/405/GCA_000001405.2_B",
                          "all/GCF/000/006/805/GCF_000006805.1_ASM680v1"})

    def test_files_at_any_level_are_not_genomes(self):
        self.leaf("all/GCA/000/001/405/GCA_000001405.1_A")
        write(self.path("all/README.txt"), "")
        write(self.path("all/GCA/000/001/stray.txt"), "")
        write(self.path("all/GCA/000/001/405/stray.txt"), "")
        self.assertEqual(N.mirror_genome_dirs(self.dir), {"all/GCA/000/001/405/GCA_000001405.1_A"})

    def test_a_symlinked_leaf_is_neither_counted_nor_reached(self):
        # a link into another tree must not be deleted through
        self.leaf("all/GCA/000/001/405/GCA_000001405.1_A")
        os.symlink(self.path("all/GCA/000/001/405/GCA_000001405.1_A"),
                   self.path("all/GCA/000/001/405/GCA_000001405.9_LINK"))
        self.assertEqual(N.mirror_genome_dirs(self.dir), {"all/GCA/000/001/405/GCA_000001405.1_A"})

    def test_a_root_with_no_all_directory_holds_no_genomes(self):
        # the first run against a fresh --root
        self.assertEqual(N.mirror_genome_dirs(self.dir), set())

    def test_whatever_ncbi_adds_under_all_is_walked_alike(self):
        self.leaf("all/GCX/000/000/001/GCX_000000001.1_NEW")
        self.assertEqual(N.mirror_genome_dirs(self.dir), {"all/GCX/000/000/001/GCX_000000001.1_NEW"})


class PlanMirror(TempDirCase):
    def genomes(self, *asms):
        write(self.path("sel.tsv"), selection(*asms))
        return N.read_assembly_summary(self.path("sel.tsv"))[0]

    def test_matches_on_the_exact_directory_the_ftp_path_maps_to(self):
        # a renamed assembly (same accession, new name) is to be removed and re-fetched,
        # not left beside its replacement
        for rel in ("all/GCA/000/001/405/GCA_000001405.1_OLDNAME",
                    "all/GCA/000/001/405/GCA_000001405.2_B"):
            os.makedirs(self.path(rel))
        to_remove, to_add, present = N.plan_mirror(
            self.dir, self.genomes("GCA_000001405.1_NEWNAME", "GCA_000001405.2_B"))
        self.assertEqual(to_remove, ["all/GCA/000/001/405/GCA_000001405.1_OLDNAME"])
        self.assertEqual((to_add, present), (1, 1))

    def test_an_empty_mirror_has_everything_to_add(self):
        self.assertEqual(N.plan_mirror(self.dir, self.genomes("GCA_000001405.1_A")), ([], 1, 0))


class RemoveGenomeDir(TempDirCase):
    def leaf(self, rel):
        os.makedirs(self.path(rel))
        write(self.path(rel + "/md5checksums.txt"), "")

    def test_prunes_the_triplets_it_empties_but_never_all(self):
        self.leaf("all/GCA/000/001/405/GCA_000001405.1_A")
        N.remove_genome_dir(self.dir, "all/GCA/000/001/405/GCA_000001405.1_A")
        self.assertFalse(os.path.exists(self.path("all/GCA")))
        self.assertTrue(os.path.isdir(self.path("all")))

    def test_stops_at_the_first_level_that_still_has_contents(self):
        self.leaf("all/GCA/000/001/405/GCA_000001405.1_A")
        self.leaf("all/GCA/000/001/999/GCA_000001999.1_C")
        N.remove_genome_dir(self.dir, "all/GCA/000/001/405/GCA_000001405.1_A")
        self.assertFalse(os.path.exists(self.path("all/GCA/000/001/405")))
        self.assertTrue(os.path.isdir(self.path("all/GCA/000/001/999/GCA_000001999.1_C")))


class PruneMirror(TempDirCase):
    def leaf(self, rel):
        os.makedirs(self.path(rel))
        write(self.path(rel + "/md5checksums.txt"), "")

    def test_records_each_directory_before_deleting_it(self):
        self.leaf("all/GCA/000/001/405/GCA_000001405.1_A")
        self.leaf("all/GCA/000/001/405/GCA_000001405.2_B")
        record = self.path("run.rm")
        removed, failed = N.prune_mirror(self.dir, ["all/GCA/000/001/405/GCA_000001405.1_A",
                                                    "all/GCA/000/001/405/GCA_000001405.2_B"],
                                         record, silent=True)
        self.assertEqual((removed, failed), (2, 0))
        with open(record) as handle:
            lines = handle.read().splitlines()
        # rows land in completion order, the removals running on threads
        self.assertEqual(lines[0], N.RM_HEADER.rstrip("\n"))
        self.assertEqual(sorted(lines[1:]),
                         ["GCA_000001405.1\tall/GCA/000/001/405/GCA_000001405.1_A",
                          "GCA_000001405.2\tall/GCA/000/001/405/GCA_000001405.2_B"])
        self.assertFalse(os.path.exists(self.path("all/GCA")))

    def test_an_empty_removal_still_truncates_the_record(self):
        # a stale record from an earlier run must not be mistaken for this run's
        write(self.path("run.rm"), "GCA_old\tall/x\n")
        self.assertEqual(N.prune_mirror(self.dir, [], self.path("run.rm"), silent=True), (0, 0))
        with open(self.path("run.rm")) as handle:
            self.assertEqual(handle.read(), N.RM_HEADER)

    @unittest.skipIf(os.geteuid() == 0, "root can delete anything")
    def test_a_directory_that_will_not_delete_is_counted_not_fatal(self):
        # another member's file (SHARED OPERATION): the sync that follows still runs
        self.leaf("all/GCA/000/001/405/GCA_000001405.1_A")
        self.leaf("all/GCA/000/002/405/GCA_000002405.1_B")
        os.chmod(self.path("all/GCA/000/002/405"), 0o555)      # its parent refuses the rmdir
        try:
            removed, failed = N.prune_mirror(
                self.dir, ["all/GCA/000/001/405/GCA_000001405.1_A",
                           "all/GCA/000/002/405/GCA_000002405.1_B"], self.path("run.rm"), silent=True)
        finally:
            os.chmod(self.path("all/GCA/000/002/405"), 0o755)
        self.assertEqual((removed, failed), (1, 1))
        self.assertFalse(os.path.exists(self.path("all/GCA/000/001/405")))
        with open(self.path("run.rm")) as handle:
            self.assertEqual(len(handle.read().splitlines()), 3)   # header + both rows

    def test_stops_between_directories_once_a_stop_is_requested(self):
        self.leaf("all/GCA/000/001/405/GCA_000001405.1_A")
        stop, N.STOP = N.STOP, N.StopFlag()
        N.STOP.set("interrupted by SIGINT", 2)
        try:
            removed, failed = N.prune_mirror(self.dir, ["all/GCA/000/001/405/GCA_000001405.1_A"],
                                             self.path("run.rm"), silent=True)
        finally:
            N.STOP = stop
        self.assertEqual((removed, failed), (0, 0))
        self.assertTrue(os.path.isdir(self.path("all/GCA/000/001/405/GCA_000001405.1_A")))


class RemovalEndToEndOffline(MainRunner):
    """main() with a mirror on disk and FakeNCBI standing in for the network: what the
    selection does not list goes, what it lists is fetched, and --retry and --dry-run
    remove nothing."""

    A = "GCA_000001405.1_A"                      # on disk, not in the selection
    B = "GCA_000001405.2_B"                      # in the selection, served by FakeNCBI

    def setUp(self):
        super().setUp()
        self.files = {self.B + "_genomic.fna.gz": b"ACGT" * 50,
                      self.B + "_assembly_report.txt": b"report",
                      "annotation_hashes.txt": b"hashes"}
        self.real_http_get, N.http_get = N.http_get, FakeNCBI(self.B, self.files)
        self.root = self.path("mirror")
        os.makedirs(os.path.join(self.root, "all/GCA/000/001/405", self.A))
        write(os.path.join(self.root, "all/GCA/000/001/405", self.A, "md5checksums.txt"), "")

    def tearDown(self):
        N.http_get = self.real_http_get
        super().tearDown()

    def a_exists(self):
        return os.path.isdir(os.path.join(self.root, "all/GCA/000/001/405", self.A))

    def b_synced(self):
        return os.path.isfile(os.path.join(self.root, "all/GCA/000/001/405", self.B, "md5checksums.txt"))

    def test_the_selection_removes_what_it_does_not_list_then_syncs_what_it_does(self):
        write(self.path("sel.tsv"), selection(self.B))
        rc, err = self.run_main("sel.tsv", "--root", self.root)
        self.assertEqual(rc, 0, err)
        self.assertFalse(self.a_exists())
        self.assertTrue(self.b_synced())
        self.assertIn("1 genome dir(s) to remove, 1 to add, 0 present", err)
        with open(self.path("logs/sel.rm")) as handle:
            self.assertEqual(handle.read(),
                             N.RM_HEADER + "GCA_000001405.1\tall/GCA/000/001/405/%s\n" % self.A)

    def test_a_retry_syncs_what_it_lists_and_removes_nothing(self):
        write(self.path("sel.fail"), selection(self.B, header=N.FAIL_HEADER).replace("\tna\n", "\tna\tr\n"))
        rc, err = self.run_main("sel.fail", "--root", self.root, retry=True)
        self.assertEqual(rc, 0, err)
        self.assertTrue(self.a_exists())
        self.assertTrue(self.b_synced())
        self.assertFalse(os.path.exists(self.path("logs/sel.fail.rm")))

    def test_dry_run_reports_counts_lists_removals_and_changes_nothing(self):
        write(self.path("sel.tsv"), selection(self.B))
        calls_before = len(N.http_get.calls)
        rc, err = self.run_main("sel.tsv", "--root", self.root, "--dry-run")
        self.assertEqual(rc, 0, err)
        self.assertIn("DRY RUN: 1 genome dir(s) would be removed", err)
        self.assertIn("1 added, 0 present", err)
        self.assertTrue(self.a_exists())
        self.assertFalse(self.b_synced())
        self.assertEqual(len(N.http_get.calls), calls_before)          # no request at all
        with open(self.path("logs/sel.rm_dry_run")) as handle:
            lines = handle.read().splitlines()
        self.assertTrue(lines[0].startswith("# DRY RUN"))
        self.assertEqual(lines[-1], "GCA_000001405.1\tall/GCA/000/001/405/%s" % self.A)
        for name in ("sel.rm", "sel.fail", "sel.bad"):
            self.assertFalse(os.path.exists(self.path("logs/" + name)), name)

    def test_dry_run_takes_no_lock_so_it_can_run_beside_a_live_sync(self):
        write(self.path("sel.tsv"), selection(self.B))
        held = N.lock_root(self.root)
        try:
            rc, err = self.run_main("sel.tsv", "--root", self.root, "--dry-run")
        finally:
            held.close()
        self.assertEqual(rc, 0, err)

    def test_dry_run_against_a_missing_root_is_a_usage_error(self):
        write(self.path("sel.tsv"), selection(self.B))
        rc, err = self.run_main("sel.tsv", "--root", self.path("nowhere"), "--dry-run")
        self.assertEqual(rc, 2)
        self.assertIn("not a directory", err)
        self.assertFalse(os.path.exists(self.path("nowhere")))

    def extra_rows(self, name):
        with open(self.path("logs/" + name)) as handle:
            return handle.read()

    def test_verify_only_reports_what_the_selection_does_not_list_and_removes_nothing(self):
        # both halves of "the mirror equals the selection": B listed but absent -> .bad;
        # A present but unlisted -> .extra; the mirror is untouched without --delete
        write(self.path("sel.tsv"), selection(self.B))
        rc, err = self.run_main("sel.tsv", "--root", self.root, "--verify-only")
        self.assertEqual(rc, 1)
        self.assertTrue(self.a_exists())
        self.assertFalse(os.path.exists(self.path("logs/sel.rm")))
        self.assertEqual(self.extra_rows("sel.extra"),
                         N.RM_HEADER + "GCA_000001405.1\tall/GCA/000/001/405/%s\n" % self.A)
        self.assertIn(self.B, self.extra_rows("sel.bad"))
        self.assertIn("1 genome dir(s) not in the selection", err)
        self.assertIn("missing directory", err)

    def test_verify_only_with_delete_removes_what_the_selection_does_not_list(self):
        write(self.path("sel.tsv"), selection(self.B))
        rc, err = self.run_main("sel.tsv", "--root", self.root, "--verify-only", "--delete")
        self.assertEqual(rc, 1)                                  # B is still missing
        self.assertFalse(self.a_exists())
        self.assertIn("GCA_000001405.1\tall/GCA/000/001/405/%s" % self.A, self.extra_rows("sel.rm"))
        self.assertIn("GCA_000001405.1\tall/GCA/000/001/405/%s" % self.A, self.extra_rows("sel.extra"))
        self.assertIn("1 removed", err)

    def test_a_mirror_equal_to_the_selection_verifies_clean_exit_0(self):
        write(self.path("sel.tsv"), selection(self.B))
        self.assertEqual(self.run_main("sel.tsv", "--root", self.root)[0], 0)   # removes A, fetches B
        rc, err = self.run_main("sel.tsv", "--root", self.root, "--verify-only")
        self.assertEqual(rc, 0, err)
        self.assertEqual(self.extra_rows("sel.extra"), N.RM_HEADER)
        self.assertIn("holds nothing else", err)

    def test_a_retry_verifies_only_what_it_lists(self):
        # a retry file is a subset, so what it does not list is not "extra"
        write(self.path("sel.fail"), selection(self.B, header=N.FAIL_HEADER).replace("\tna\n", "\tna\tr\n"))
        self.assertEqual(self.run_main("sel.fail", "--root", self.root, retry=True)[0], 0)  # fetches B
        rc, err = self.run_main("sel.fail", "--root", self.root, "--verify-only", retry=True)
        self.assertEqual(rc, 0, err)
        self.assertTrue(self.a_exists())
        self.assertFalse(os.path.exists(self.path("logs/sel.fail.extra")))
        self.assertNotIn("holds nothing else", err)


class ParserContract(unittest.TestCase):
    """The interface both entry points share."""

    def parse(self, *argv):
        return N.build_parser().parse_args(list(argv))

    def test_rate_defaults_to_20(self):
        self.assertEqual(self.parse("--gtdb_selected_genomes", "s", "--root", "g", "-l", "x.log").rate, 20.0)

    def test_root_and_log_are_required_and_root_has_no_default(self):
        with self.assertRaises(SystemExit):
            self.parse("--gtdb_selected_genomes", "s", "-l", "x.log")
        with self.assertRaises(SystemExit):
            self.parse("--gtdb_selected_genomes", "s", "--root", "g")

    def test_selection_and_retry_are_mutually_exclusive_and_one_is_required(self):
        with self.assertRaises(SystemExit):
            self.parse("--root", "g", "-l", "x.log")
        with self.assertRaises(SystemExit):
            self.parse("--gtdb_selected_genomes", "s", "--retry", "r", "--root", "g", "-l", "x.log")

    def test_nfs_jobs_defaults_to_the_measured_knee(self):
        args = self.parse("--gtdb_selected_genomes", "s", "--root", "g", "-l", "x.log")
        self.assertEqual(args.nfs_jobs, N.NFS_JOBS)
        self.assertEqual(args.nfs_jobs, 8)

    def test_nfs_jobs_is_tunable(self):
        # the mirror is shared, so an operator may need to back off
        args = self.parse("--gtdb_selected_genomes", "s", "--root", "g", "-l", "x.log",
                          "--nfs-jobs", "2")
        self.assertEqual(args.nfs_jobs, 2)

    def test_the_selection_lands_in_summary_as_the_module_expects(self):
        # the whole module, and these tests, read args.summary
        self.assertEqual(self.parse("--gtdb_selected_genomes", "s", "--root", "g", "-l", "x.log").summary, "s")
        args = self.parse("--retry", "r", "--root", "g", "-l", "x.log")
        self.assertEqual((args.summary, args.retry), (None, "r"))

    def test_the_selection_takes_g_as_list_genomes_gives_it(self):
        # the same table is -g to list_genomes, so it is -g here too
        self.assertEqual(self.parse("-g", "s", "--root", "g", "-l", "x.log").summary, "s")


class NfsJobsIsHonoured(TempDirCase):
    """--nfs-jobs must actually reach the walk and the removal, not just parse."""

    def leaf(self, rel):
        os.makedirs(self.path(rel))
        write(self.path(rel + "/md5checksums.txt"), "")

    def test_the_walk_runs_on_the_requested_threads(self):
        self.leaf("all/GCA/000/001/405/GCA_000001405.1_A")
        seen = []
        real = N.ThreadPoolExecutor

        def spy(max_workers=None, **kw):
            seen.append(max_workers)
            return real(max_workers=max_workers, **kw)

        N.ThreadPoolExecutor = spy
        try:
            N.mirror_genome_dirs(self.dir, silent=True, workers=3)
        finally:
            N.ThreadPoolExecutor = real
        self.assertEqual(seen, [3])

    def test_the_removal_runs_on_the_requested_threads(self):
        self.leaf("all/GCA/000/001/405/GCA_000001405.1_A")
        seen = []
        real = N.bounded_map

        def spy(fn, items, workers):
            seen.append(workers)
            return real(fn, items, workers)

        N.bounded_map = spy
        try:
            N.prune_mirror(self.dir, ["all/GCA/000/001/405/GCA_000001405.1_A"],
                           self.path("r.rm"), silent=True, workers=5)
        finally:
            N.bounded_map = real
        self.assertEqual(seen, [5])
        self.assertFalse(os.path.exists(self.path("all/GCA")))

    def test_one_thread_still_removes_everything(self):
        # --nfs-jobs 1 is the serial path an operator backs off to
        self.leaf("all/GCA/000/001/405/GCA_000001405.1_A")
        self.leaf("all/GCA/000/002/405/GCA_000002405.1_B")
        removed, failed = N.prune_mirror(
            self.dir, ["all/GCA/000/001/405/GCA_000001405.1_A",
                       "all/GCA/000/002/405/GCA_000002405.1_B"],
            self.path("r.rm"), silent=True, workers=1)
        self.assertEqual((removed, failed), (2, 0))
        self.assertFalse(os.path.exists(self.path("all/GCA")))


# --------------------------------------------------------------- thousands separators

import re


class CountFormatting(unittest.TestCase):
    """A release is millions of genomes; a bare 1913482 in a log line cannot be
    read at a glance, and two of them cannot be compared."""

    def test_a_count_separates_thousands(self):
        self.assertEqual(N.format_count(1234), "1,234")
        self.assertEqual(N.format_count(1913482), "1,913,482")
        self.assertEqual(N.format_count(0), "0")
        self.assertEqual(N.format_count(999), "999")

    def test_an_amount_separates_thousands_and_keeps_its_decimals(self):
        self.assertEqual(N.format_amount(1234.56), "1,234.6")
        self.assertEqual(N.format_amount(1234.56, 2), "1,234.56")
        self.assertEqual(N.format_amount(1234.56, 0), "1,235")

    def test_the_names_reported_carry_the_count_separated(self):
        self.assertIn("(all 1,234 listed in the log)", N.first_names(["x"] * 1234))


class ProgressSnapshotFormatting(unittest.TestCase):
    """The progress line carries most of the numbers a long run shows."""

    def snapshot(self, **counters):
        prog = N.Progress("Downloading", 5678901, silent=True)
        for name, value in counters.items():
            setattr(prog, name, value)

        records = []

        class _Capture(logging.Handler):
            def emit(self, record):
                records.append(record.getMessage())

        handler = _Capture()
        logger = logging.getLogger("timestamp")
        logger.addHandler(handler)
        try:
            prog.snapshot()
        finally:
            logger.removeHandler(handler)

        self.assertEqual(len(records), 1)
        return records[0]

    def test_the_counts_of_a_progress_line_are_separated(self):
        message = self.snapshot(done=1234567, skipped=234567, fresh=123456, failed=1234)
        self.assertIn("1,234,567/5,678,901", message)
        self.assertIn("unchanged=234,567", message)
        self.assertIn("fresh=123,456", message)
        self.assertIn("failed=1,234", message)

    def test_the_http_counters_of_a_progress_line_are_separated(self):
        with N._stats_lock:
            before = dict(N.STATS)
            N.STATS.update(requests=2345678, conn_drops=1234, probe_404=4321,
                           breaker_trips=1111, paused_s=98765.0, bytes=9876543210)
        try:
            message = self.snapshot(done=1234567)
        finally:
            with N._stats_lock:
                N.STATS.update(before)

        self.assertIn("requests=2,345,678", message)
        self.assertIn("drops=1,234", message)
        self.assertIn("probe404=4,321", message)
        self.assertIn("trips=1,111", message)
        self.assertIn("paused=98,765s", message)
        self.assertIn("MB=9,876.5", message)

    def test_no_count_is_written_as_a_bare_run_of_digits(self):
        # percentages and the fixed fields of a duration are the exceptions, and
        # neither reaches four digits
        message = self.snapshot(done=1234567, skipped=234567, fresh=123456, failed=1234)
        for word in message.split():
            for number in re.findall(r"(?<![\d,.])\d+", word):
                self.assertLess(len(number), 4, message)
