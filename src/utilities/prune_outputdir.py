#!/usr/bin/env python3
"""Delete jacobian_runs/*/OutputDir files that both processing stages are done with.

Gated entirely on the completion markers, so it scales to every target face
without inspecting any of them by hand. Nothing here reads file contents: the
markers are the statement that the contents were already verified.

Per face, all of the following must hold before a single file is removed:

  1. Exactly one overpass_complete.* and one inversion_complete.* marker
     exist at the run directory, and they carry the same start date and the
     same shared end date S.
  2. The cutoff comes from the marker's S, never from the checkpoints as they
     stand now. S grows as the Jacobian runs advance, so recomputing it here
     would delete beyond what the markers actually cover.
  3. Every Jacobian run holds an overpass file for every date in the window.
     Compared as an explicit set rather than a count, because a count still
     matches when one date is absent and another is duplicated.
  4. inversion/data_converted agrees with the manifest jacobian.py wrote:
     every granule recorded as written or cached has its pickle on disk, and
     the manifest covers the same window as the markers. Disable with
     --skip-converted-check for a face whose inversion predates manifests.

Every check is a path check. Nothing here knows or cares what the underlying
storage is: on pcluster with FSx autoexport, deleting the file propagates the
delete to S3 on its own. Uploading products to S3, and removing objects from a
bucket nothing syncs for you, is a separate concern handled by
s3_upload_and_prune.sh, which is a one-off for working from downloaded data.

What is safe to delete, and why: overpass output for local date L consumes
OutputDir for L and L+1, so once the marker covers local dates through
last_local, the earliest work any later round can do is local date
last_local + 1, which needs OutputDir from last_local + 1 onward. Files dated
last_local or earlier are therefore referenced by nothing. The inversion is
looser still: with SatDiagOperator false it reads UTC dates up to S-1, and
its next round starts at S.

That argument holds for an early-afternoon overpass. A morning one needs the
preceding UTC day as well, so the cutoff is held back a further day whenever
OverpassTime says so -- see overpass_lookback_days.

Dry run unless --execute is passed.

Usage:
    prune_outputdir.py configs_C36S10/config_T001.yml
    prune_outputdir.py configs_C36S10/config_T*.yml --execute
    prune_outputdir.py configs_C36S10/config_T*.yml --record-deleted pruned.txt

Exit codes:
    0  every face processed, or skipped for a stated reason
    1  usage or configuration error
    2  at least one face was blocked by a failed precondition
"""

import argparse
import glob
import json
import os
import re
import sys
from datetime import datetime, timedelta

import yaml

# GEOSChem.<collection>.YYYYMMDD_0000z.nc4
OUTPUTDIR_FILE_RE = re.compile(
    r"^GEOSChem\..+\.(\d{8})_0000z\.nc4$"
)

### gcchem_internal_checkpoint.YYYYMMDD_0000z.nc4 -- deliberately the same
### pattern get_shared_end_date reads, so the set pruned here is exactly the set
### S is derived from, and keeping the newest preserves S unchanged.
##CHECKPOINT_FILE_RE = re.compile(
##    r"^gcchem_internal_checkpoint\.(\d{8})_0000z\.nc4$"
##)

# Either gcchem_internal_checkpoint.YYYYMMDD_0000z.nc4 or
# GEOSChem.Restart.YYYYMMDD_0000z.cN.nc4 -- deliberately the same
# patterns get_shared_end_date reads, so the set pruned here is exactly the set
# S is derived from, and keeping the newest preserves S unchanged.
CHECKPOINT_FILE_RE = re.compile(
    r"^(?:GEOSChem\.Restart|gcchem_internal_checkpoint)"
    r"\.(\d{8})_0000z(?:\.c\d+)?\.nc4$"
)

# <stage>_complete.<StartDate>_S<shared_end_date>. The second field is the
# shared end date the stage ran against, not the end of its coverage.
MARKER_RE = re.compile(
    r"^(?P<stage>[a-z]+)_complete\.(?P<start>\d{8})_S(?P<end>\d{8})$"
)

MANIFEST_NAME = "outputdir_pruned.json"

REPO_ROOT = os.path.dirname(
    os.path.dirname(
        os.path.dirname(os.path.abspath(__file__))
    )
)


class _Tee:
    """Write to a stream and a log file at once.

    Flushed on every write so `tail -f` follows a backgrounded run live.
    """

    def __init__(self, stream, handle):
        self.stream = stream
        self.handle = handle

    def write(self, text):
        self.stream.write(text)
        self.handle.write(text)
        self.handle.flush()
        return len(text)

    def flush(self):
        self.stream.flush()
        self.handle.flush()


def parse_date(text):
    return datetime.strptime(text, "%Y%m%d")


def read_marker(run_dirs, stage):
    """Return (start_date, shared_end_date) for one stage, or raise."""
    matches = glob.glob(os.path.join(run_dirs, f"{stage}_complete.*"))

    if not matches:
        raise LookupError(f"no {stage}_complete marker")

    if len(matches) > 1:
        # write_stage_marker removes superseded markers, so more than one
        # means something wrote them by hand. Refuse rather than guess.
        raise LookupError(
            f"{len(matches)} {stage}_complete markers: "
            + ", ".join(sorted(os.path.basename(m) for m in matches))
        )

    name = os.path.basename(matches[0])
    match = MARKER_RE.match(name)

    if match is None:
        raise LookupError(f"unparsable marker name: {name}")

    return match.group("start"), match.group("end")


def overpass_lookback_days(config):
    """UTC days before a local overpass date that the operator may need to read.

    Derived from OverpassTime rather than assumed, because it follows from the
    arithmetic. The operator sets utc_offset = longitude / 15, so the offset
    reaches +12 at longitude 180, and UTC = local_hour - offset. That lands on
    the previous day exactly when local_hour < offset -- so when the local
    overpass is before 12:00, and never otherwise.

    TROPOMI's 13:30 gives 0: local date L falls on UTC L or L+1, which is what
    process_run_day reads. A 09:30 overpass gives 1, since 09:30 local at UTC+12
    is 21:30 the day before, and a cutoff that ignored this would already have
    deleted it.

    An unreadable value returns 1. Retaining a day too many costs one date of
    OutputDir per face; retaining one too few destroys data that only
    re-simulation can restore.

    This preserves the inputs only. An operator that actually reads L-1 needs
    process_run_day extended too; it currently reads L and L+1.
    """
    raw = str(config.get("OverpassTime", "13:30")).strip().strip('"').strip("'")

    try:
        fields = raw.split(":")
        local_hour = float(fields[0])

        if len(fields) > 1:
            local_hour += float(fields[1]) / 60.0

    except (ValueError, IndexError):
        print(f"  WARNING: unreadable OverpassTime {raw!r}; assuming morning")
        return 1

    return 1 if local_hour < 12.0 else 0


def local_date_list(start_date, end_date, shared_end_date):
    """Reproduce build_date_list's local overpass dates for a given S.

    Kept in step with calculate_satellite_overpass_diagnostics.build_date_list:
    the list starts one day before StartDate, and stops one day short of S
    unless S has reached EndDate.
    """
    local_start = parse_date(start_date) - timedelta(days=1)

    local_end_exclusive = parse_date(shared_end_date)

    if local_end_exclusive != parse_date(end_date):
        local_end_exclusive -= timedelta(days=1)

    n_days = (local_end_exclusive - local_start).days

    if n_days <= 0:
        return []

    return [
        (local_start + timedelta(days=i)).strftime("%Y%m%d")
        for i in range(n_days)
    ]


def missing_overpass_files(output_dir, date_list, overpass_tag, file_prefixes):
    """Return the overpass filenames one run directory should hold but does not.

    Compared as an explicit set rather than a count: a count matches even when
    one date is absent and another is duplicated, which is exactly the failure
    a prune must not act on.
    """
    try:
        present = set(os.listdir(output_dir))
    except (FileNotFoundError, NotADirectoryError, PermissionError):
        present = set()

    missing = []

    for date_str in date_list:
        for file_prefix in file_prefixes:
            name = f"{file_prefix}.overpass.{date_str}_{overpass_tag}.nc4"

            if name not in present:
                missing.append(name)

    return missing


def check_converted_manifest(run_dirs, start_date, shared_end_date):
    """Verify data_converted against the manifest jacobian.py wrote.

    The manifest names every granule the window selected and what became of
    it, which makes this an exact check: a granule recorded as written or
    cached must have its pickle on disk, and one recorded as no_valid_obs must
    not, because the operator found nothing to save. Inferring completeness
    from which dates happen to have output cannot distinguish those two.
    """
    inversion_dir = os.path.join(run_dirs, "inversion")

    manifest_path = os.path.join(inversion_dir, "data_converted_manifest.json")

    if not os.path.isfile(manifest_path):
        return [
            "no data_converted_manifest.json; rerun the inversion so it is "
            "written, or fall back to --skip-converted-check"
        ]

    try:
        with open(manifest_path) as handle:
            manifest = json.load(handle)
    except (OSError, ValueError) as error:
        return [f"data_converted_manifest.json is unreadable: {error}"]

    if (
        manifest.get("start_date") != start_date
        or manifest.get("shared_end_date") != shared_end_date
    ):
        return [
            "manifest covers "
            f"{manifest.get('start_date')}_{manifest.get('shared_end_date')}, "
            f"but the markers cover {start_date}_{shared_end_date}"
        ]

    converted_dir = os.path.join(inversion_dir, "data_converted")

    try:
        present = set(os.listdir(converted_dir))
    except (FileNotFoundError, NotADirectoryError, PermissionError):
        return [f"data_converted directory is absent: {converted_dir}"]

    problems = []

    for entry in manifest.get("granules", []):
        granule = entry.get("granule")
        status = entry.get("status")

        pickle_name = f"{granule}_GCtoTROPOMI.pkl"

        if status in ("written", "cached"):
            if pickle_name not in present:
                problems.append(f"recorded as {status} but absent: {pickle_name}")

        elif status == "no_valid_obs":
            continue

        else:
            problems.append(f"unrecognized status {status!r} for {granule}")

    return problems


def collect_stale_checkpoints(run_dirs, run_name, shared_end_date, keep):
    """Return checkpoints that are older than S and not a run's newest.

    Two rules, and both are needed.

    Older than S. OutputDir dates through S-1 were opened during processing --
    overpass reads local date L and L+1, running to L = S-2, and the inversion
    reads granules to S-1 -- so corruption there would already have surfaced.
    Dates at or after S exist only in runs that ran ahead of the others, were
    never opened, and are RETAINED by the prune, whose cutoff sits at least a
    day earlier. Keeping the checkpoints from S onward is what preserves the
    ability to re-simulate that unvalidated tail.

    Not the newest. S is the minimum across runs of each run's maximum, so the
    slowest run's newest checkpoint is dated exactly S. A rule of "date <= S"
    would strip it of every checkpoint, and get_shared_end_date would then fall
    back to start_date and raise. Strictly-older-than-S already implies this,
    since every run's maximum is >= S, but it is enforced explicitly rather
    than left as a consequence.
    """
    # keep <= 0 means do not prune checkpoints at all. Spelled out here rather
    # than left to the caller: with keep == 0 the slice below would select
    # every checkpoint, and a run whose checkpoints all predate S would lose
    # the lot, leaving get_shared_end_date nothing to read.
    if keep <= 0:
        return []

    jacobian_root = os.path.join(run_dirs, "jacobian_runs")

    stale = []

    for run_dir_name in sorted(os.listdir(jacobian_root)):
        run_dir = os.path.join(jacobian_root, run_dir_name)

        if not os.path.isdir(run_dir) or not run_dir_name.startswith(
            f"{run_name}_"
        ):
            continue

        restarts_dir = os.path.join(run_dir, "Restarts")

        try:
            entries = os.listdir(restarts_dir)
        except (FileNotFoundError, NotADirectoryError, PermissionError):
            continue

        #dated = []
        #for name in entries:
        #    match = CHECKPOINT_FILE_RE.match(name)
        #    if match:
        #        dated.append((match.group(1), name))
        dated = []
        for name in entries:
            if os.path.islink(os.path.join(restarts_dir, name)):
                # Never prune links. The initial restart GEOSChem.Restart.<StartDate>
                # is a link to a shared file, needed to rerun from the start, and
                # now matches the pattern.
                continue
            match = CHECKPOINT_FILE_RE.match(name)
            if match:
                dated.append((match.group(1), name))
        
        # Filenames sort chronologically.
        dated.sort()

        if not dated:
            continue

        # Two conditions must both hold before a checkpoint goes: it is not
        # among the newest `keep`, and it predates S.
        for date_str, name in dated[:-keep]:
            if date_str < shared_end_date:
                stale.append(os.path.join(restarts_dir, name))

    return sorted(stale)


def collect_deletable(run_dirs, run_name, cutoff_date):
    """Return OutputDir paths dated at or before cutoff_date."""
    jacobian_root = os.path.join(run_dirs, "jacobian_runs")

    deletable = []

    for run_dir_name in sorted(os.listdir(jacobian_root)):
        run_dir = os.path.join(jacobian_root, run_dir_name)

        if not os.path.isdir(run_dir) or not run_dir_name.startswith(
            f"{run_name}_"
        ):
            continue

        output_dir = os.path.join(run_dir, "OutputDir")

        try:
            entries = os.scandir(output_dir)
        except (FileNotFoundError, NotADirectoryError, PermissionError):
            continue

        with entries:
            for entry in entries:
                match = OUTPUTDIR_FILE_RE.match(entry.name)

                if match is None:
                    continue

                if match.group(1) <= cutoff_date:
                    deletable.append(entry.path)

    return sorted(deletable)


def read_manifest(run_dirs):
    path = os.path.join(run_dirs, MANIFEST_NAME)

    if not os.path.isfile(path):
        return None

    try:
        with open(path) as handle:
            return json.load(handle)
    except (OSError, ValueError):
        return None


def write_manifest(run_dirs, cutoff_date, deleted_count, timestamp):
    path = os.path.join(run_dirs, MANIFEST_NAME)

    previous = read_manifest(run_dirs) or {"history": []}

    previous["cutoff_date"] = cutoff_date
    previous["history"].append({
        "cutoff_date": cutoff_date,
        "deleted": deleted_count,
        "at": timestamp,
    })

    with open(path, "w") as handle:
        json.dump(previous, handle, indent=2)


def process_face(config_path, args, timestamp):
    """Prune one face. Returns "pruned", "skipped" or "blocked"."""
    with open(config_path) as config_file:
        config = yaml.load(config_file, Loader=yaml.FullLoader)

    run_name = config["RunName"]
    run_dirs = os.path.join(
        os.path.expandvars(config["OutputPath"]),
        run_name,
    )

    def report(status, message):
        print(f"[{status}] {run_name}: {message}")

    if not os.path.isdir(run_dirs):
        report("BLOCKED", f"run directory does not exist: {run_dirs}")
        return "blocked"

    # ---- precondition 1: both markers, in agreement ----
    try:
        overpass_start, overpass_end = read_marker(run_dirs, "overpass")
        inversion_start, inversion_end = read_marker(run_dirs, "inversion")
    except LookupError as error:
        report("BLOCKED", str(error))
        return "blocked"

    if (overpass_start, overpass_end) != (inversion_start, inversion_end):
        report(
            "BLOCKED",
            "markers disagree: overpass covers "
            f"{overpass_start}_{overpass_end}, inversion covers "
            f"{inversion_start}_{inversion_end}",
        )
        return "blocked"

    start_date = overpass_start
    shared_end_date = overpass_end

    if start_date != str(config["StartDate"]):
        report(
            "BLOCKED",
            f"marker start {start_date} does not match config StartDate "
            f"{config['StartDate']}",
        )
        return "blocked"

    # ---- precondition 2: cutoff from the marker, not from the checkpoints ----
    date_list = local_date_list(
        start_date,
        str(config["EndDate"]),
        shared_end_date,
    )

    if not date_list:
        report("BLOCKED", f"marker window {shared_end_date} covers no dates")
        return "blocked"

    # date_list[-1] is the last local date the marker covers, which is what the
    # unreferenced argument above is stated in terms of. The lookback is then
    # subtracted so an instrument reading earlier than local date L still finds
    # its inputs; see overpass_lookback_days.
    lookback = overpass_lookback_days(config)

    cutoff_date = (
        parse_date(date_list[-1]) - timedelta(days=lookback)
    ).strftime("%Y%m%d")

    # Optional further margin: hold OutputDir N days behind even that, so a
    # product that turns out wrong can still be recomputed.
    if args.retain_days:
        cutoff_date = (
            parse_date(cutoff_date) - timedelta(days=args.retain_days)
        ).strftime("%Y%m%d")

    if cutoff_date < start_date:
        margin = f", --retain-days {args.retain_days}" if args.retain_days else ""
        report(
            "SKIPPED",
            f"cutoff {cutoff_date} is before {start_date} "
            f"(lookback {lookback}{margin}); nothing safe to delete",
        )
        return "skipped"

    # ---- record only ----
    if args.record_manifest_only:
        write_manifest(run_dirs, cutoff_date, 0, timestamp)
        report("RECORDED", f"manifest set to cutoff {cutoff_date}")
        return "skipped"

    # ---- idempotent resume ----
    manifest = read_manifest(run_dirs)

    if manifest and manifest.get("cutoff_date", "") >= cutoff_date:
        report(
            "SKIPPED",
            f"already pruned through {manifest['cutoff_date']}",
        )
        return "skipped"

    # ---- precondition 3: per-run overpass counts ----
    jacobian_root = os.path.join(run_dirs, "jacobian_runs")

    disable_run_0000 = config.get("DisableRun0000", False)

    run_dir_names = sorted(
        name
        for name in os.listdir(jacobian_root)
        if os.path.isdir(os.path.join(jacobian_root, name))
        and name.startswith(f"{run_name}_")
    )

    base_run_index = 1 if disable_run_0000 else 0
    overpass_tag = config["OverpassTime"].replace(":", "")

    run_report = []

    for run_dir_name in run_dir_names:
        suffix = run_dir_name[len(run_name) + 1:]

        if not suffix.isdigit():
            continue

        run_i = int(suffix)

        if disable_run_0000 and run_i == 0:
            continue

        # Mirrors process_run_day: the base run also carries two sampled 3D
        # collections per date.
        file_prefixes = ["GEOSChem.CH4col"] if run_i != 0 else []

        if run_i == base_run_index:
            file_prefixes += [
                "GEOSChem.BaseSpeciesConc"
                if disable_run_0000
                else "GEOSChem.SpeciesConc",
                "GEOSChem.StateMetLevEdge",
            ]

        missing = missing_overpass_files(
            os.path.join(jacobian_root, run_dir_name, "OverpassDiagnostics"),
            date_list,
            overpass_tag,
            file_prefixes,
        )

        # Every run is surveyed before reporting. Deletion is all-or-nothing
        # per face, because data_converted draws on all of them, but stopping
        # at the first incomplete run would hide how much is actually left to
        # do: one absent date and one entirely unprocessed run look identical.
        run_report.append({
            "name": run_dir_name,
            "missing": len(missing),
            "expected": len(date_list) * len(file_prefixes),
            "first": missing[0] if missing else None,
        })

    incomplete = [entry for entry in run_report if entry["missing"]]

    if incomplete:
        total_missing = sum(entry["missing"] for entry in incomplete)

        report(
            "BLOCKED",
            f"{len(incomplete)} of {len(run_report)} Jacobian run(s) have "
            f"incomplete overpass output for window "
            f"{date_list[0]}..{date_list[-1]} "
            f"({total_missing} file(s) missing in total)",
        )

        for entry in run_report:
            if not entry["missing"]:
                print(f"    {entry['name']}: complete ({entry['expected']})")
            elif entry["missing"] == entry["expected"]:
                print(
                    f"    {entry['name']}: NONE of {entry['expected']} "
                    "present -- this run was never processed"
                )
            else:
                print(
                    f"    {entry['name']}: {entry['missing']} of "
                    f"{entry['expected']} missing, first {entry['first']}"
                )

        return "blocked"

    # ---- precondition 3b: data_converted matches its manifest ----
    if not args.skip_converted_check:
        problems = check_converted_manifest(
            run_dirs,
            start_date,
            shared_end_date,
        )

        if problems:
            report(
                "BLOCKED",
                f"data_converted disagrees with its manifest "
                f"({len(problems)} problem(s)); first: {problems[0]}",
            )

            for problem in problems[:20]:
                print(f"    {problem}")

            return "blocked"

    # ---- collect and act ----
    deletable = collect_deletable(run_dirs, run_name, cutoff_date)

    stale_checkpoints = []
    if args.keep_checkpoints > 0:
        stale_checkpoints = collect_stale_checkpoints(
            run_dirs, run_name, shared_end_date, args.keep_checkpoints
        )
        deletable += stale_checkpoints

    if not deletable:
        report("SKIPPED", f"nothing dated at or before {cutoff_date}")
        return "skipped"

    total_bytes = 0

    for path in deletable:
        try:
            total_bytes += os.path.getsize(path)
        except OSError:
            pass

    gib = total_bytes / (1024 ** 3)

    # Recorded before the files go, and in dry runs too, so the list can be
    # previewed and so a bucket that does not track filesystem deletes can be
    # cleaned up from it afterwards.
    if args.record_deleted:
        with open(args.record_deleted, "a") as handle:
            for path in deletable:
                handle.write(path + "\n")

    checkpoint_note = (
        f", incl. {len(stale_checkpoints)} stale checkpoint(s)"
        if stale_checkpoints else ""
    )

    if not args.execute:
        report(
            "DRY RUN",
            f"would delete {len(deletable)} file(s), {gib:.1f} GiB, "
            f"OutputDir dated <= {cutoff_date}{checkpoint_note}",
        )

        if args.list_files:
            for path in deletable:
                print(path)

        return "pruned"

    for path in deletable:
        os.remove(path)

    # The manifest must not claim more than has happened: deleting the local
    # files is only half a prune where the bucket does not follow on its own.
    # With --defer-manifest the caller records it via --record-manifest-only
    # once the objects are confirmed gone.
    if not args.defer_manifest:
        write_manifest(
            run_dirs,
            cutoff_date,
            len(deletable),
            timestamp,
        )

    report(
        "PRUNED",
        f"deleted {len(deletable)} file(s), {gib:.1f} GiB, "
        f"OutputDir dated <= {cutoff_date}{checkpoint_note}",
    )

    return "pruned"


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "Delete OutputDir files both processing stages are finished with, "
            "gated on the completion markers."
        )
    )

    parser.add_argument(
        "configs",
        nargs="+",
        help="Face config yaml files, e.g. configs_C36S10/config_T*.yml",
    )

    parser.add_argument(
        "--execute",
        action="store_true",
        help="Actually delete. Without this nothing is removed.",
    )

    parser.add_argument(
        "--skip-converted-check",
        action="store_true",
        help=(
            "Skip verifying data_converted against its manifest. Needed only "
            "for a face whose inversion predates the manifest being written."
        ),
    )

    parser.add_argument(
        "--keep-checkpoints",
        type=int,
        default=1,
        help=(
            "Newest internal checkpoints always kept per Jacobian run, on "
            "top of the rule that only checkpoints dated before S are "
            "removed at all. 0 disables checkpoint pruning entirely. "
            "Default 1"
        ),
    )

    parser.add_argument(
        "--retain-days",
        type=int,
        default=0,
        help=(
            "Hold OutputDir this many days further back than is strictly "
            "safe, as a recovery margin. Default 0"
        ),
    )

    parser.add_argument(
        "--max-faces",
        type=int,
        default=0,
        help="Stop after this many faces are pruned; 0 means no limit",
    )

    parser.add_argument(
        "--stop-file",
        default=None,
        help="Stop cleanly between faces once this file exists",
    )

    parser.add_argument(
        "--list-files",
        action="store_true",
        help="During a dry run, print every path that would be deleted",
    )

    parser.add_argument(
        "--log",
        default=None,
        help=(
            "Log file. Defaults to "
            "logs_C36S10/prune_outputdir_<timestamp>.log in the repo. "
            "Use 'none' to write to stdout only."
        ),
    )

    parser.add_argument(
        "--record-deleted",
        default=None,
        help=(
            "Append every path deleted, or that would be deleted in a dry "
            "run, to this file. Feed it to s3_upload_and_prune.sh when the "
            "bucket does not track filesystem deletes for you."
        ),
    )

    parser.add_argument(
        "--defer-manifest",
        action="store_true",
        help=(
            "Delete the files but do not write outputdir_pruned.json. The "
            "caller records it with --record-manifest-only once the objects "
            "are confirmed gone from the bucket, so the manifest never claims "
            "a prune that only happened on this filesystem."
        ),
    )

    parser.add_argument(
        "--record-manifest-only",
        action="store_true",
        help=(
            "Write outputdir_pruned.json for the cutoff that applies now, and "
            "delete nothing. Pairs with --defer-manifest."
        ),
    )

    return parser.parse_args()


def main():
    args = parse_arguments()

    if args.keep_checkpoints < 0:
        print("ERROR: --keep-checkpoints cannot be negative", file=sys.stderr)
        return 1

    timestamp = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")

    log_handle = None

    if args.log != "none":
        if args.log:
            log_path = args.log
        else:
            log_dir = os.path.join(REPO_ROOT, "logs_C36S10")
            os.makedirs(log_dir, exist_ok=True)
            log_path = os.path.join(
                log_dir,
                "prune_outputdir_"
                + datetime.now().strftime("%Y%m%d_%H%M%S")
                + ".log",
            )

        log_handle = open(log_path, "a")
        original_stdout, original_stderr = sys.stdout, sys.stderr
        sys.stdout = _Tee(sys.stdout, log_handle)
        sys.stderr = _Tee(sys.stderr, log_handle)

        print(f"Log: {log_path}")
        print(f"Started: {timestamp}")

    counts = {"pruned": 0, "skipped": 0, "blocked": 0}

    total = len(args.configs)

    for index, config_path in enumerate(args.configs, start=1):
        # Printed before the face is touched, so a tail of the log always shows
        # what is in flight rather than only what has finished.
        print(
            f"[{index}/{total}] {datetime.now().strftime('%H:%M:%S')} "
            f"{os.path.basename(config_path)}",
            flush=True,
        )

        if args.stop_file and os.path.exists(args.stop_file):
            print(f"Stop file present, halting: {args.stop_file}")
            break

        if args.max_faces and counts["pruned"] >= args.max_faces:
            print(f"Reached --max-faces {args.max_faces}, halting.")
            break

        try:
            outcome = process_face(config_path, args, timestamp)
        except Exception as error:
            print(f"[BLOCKED] {config_path}: {error}")
            outcome = "blocked"

        counts[outcome] += 1

    print(
        "\nSummary:"
        f"\n  Pruned:  {counts['pruned']}"
        f"\n  Skipped: {counts['skipped']}"
        f"\n  Blocked: {counts['blocked']}"
        + ("" if args.execute else "\n  (dry run; pass --execute to apply)")
        + f"\n  Finished: {datetime.now().strftime('%Y-%m-%dT%H:%M:%S')}"
    )

    if log_handle is not None:
        # Restore first: closing the file while sys.stdout still wraps it makes
        # the interpreter's shutdown flush fail, which exits 120 and masks the
        # real status.
        sys.stdout, sys.stderr = original_stdout, original_stderr
        log_handle.flush()
        log_handle.close()

    return 2 if counts["blocked"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
