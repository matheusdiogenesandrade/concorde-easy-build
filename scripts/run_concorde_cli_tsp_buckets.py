#!/usr/bin/env python3
"""
Run upstream concorde_cli on TSPLIB symmetric fixtures in the same order and
10-bucket partition as external/concorde_wrapper/tests/test_concorde_symmetric_suite.cpp.

**PASS (per case):** only when concorde_cli exits 0 *and* the parsed ``Optimal Solution``
cost matches the reference entry in ``solutions.txt`` for that instance stem (within 1e-3).
Exit code 0 alone is **not** PASS unless the cost matches.

Use ``-j N`` to run up to *N* ``concorde_cli`` processes in parallel within each
batch (bucket or ``--quick``); logging order still follows the suite order.
"""
from __future__ import annotations

import argparse
import math
import os
import re
import subprocess
import sys
import threading
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

# Serialize concorde_cli stdout/stderr when running with -j > 1 and verbose output.
_verbose_print_lock = threading.Lock()

# Concorde upstream driver prints a line like: "Optimal Solution: 10628.00"
_OPTIMAL_SOLUTION_RE = re.compile(r"Optimal\s+Solution:\s*([-+0-9.eE]+)", re.IGNORECASE)


def trim_inplace_chars(s: str) -> str:
    return s.strip()


def to_upper(s: str) -> str:
    return s.upper()


def header_key_before_colon(line: str) -> str:
    col = line.find(":")
    if col == -1:
        return ""
    key = line[:col]
    return to_upper(trim_inplace_chars(key))


def header_value_after_colon(line: str) -> str:
    col = line.find(":")
    if col == -1:
        return ""
    return trim_inplace_chars(line[col + 1 :])


def read_symmetric_tsplib_dimension(filepath: Path) -> int:
    dimension = 0
    with filepath.open(encoding="utf-8", errors="replace") as f:
        for line in f:
            line = trim_inplace_chars(line)
            if not line or line[0] == "#":
                continue
            key = header_key_before_colon(line)
            if key == "DIMENSION":
                try:
                    dimension = int(header_value_after_colon(line))
                except ValueError:
                    pass
    if dimension <= 0:
        raise RuntimeError(f"Invalid or missing DIMENSION in {filepath}")
    return dimension


def stem_upper(p: Path) -> str:
    return p.stem.upper()


@dataclass(frozen=True)
class TsplibInstance:
    path: Path
    stem: str
    n: int


@dataclass
class InstanceResult:
    stem: str
    n: int
    returncode: int
    optimal: float | None
    bucket_1based: int | None  # 1–10 when running bucketed smoke; None for --quick
    reference: float | None = None  # solutions.txt (uppercase stem key)
    reference_missing: bool = False  # stem absent from solutions.txt while file was loaded
    match_ok: bool | None = None  # True if concorde cost matches reference (when both known)


def load_tsplib_solution_costs(path: Path) -> dict[str, float]:
    """Same semantics as load_tsplib_solution_costs in tsplib_symmetric_matrix.cpp."""
    out: dict[str, float] = {}
    with path.open(encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = trim_inplace_chars(raw)
            if not line or line[0] == "#":
                continue
            col = line.find(":")
            if col == -1:
                continue
            name = to_upper(trim_inplace_chars(line[:col]))
            rest = trim_inplace_chars(line[col + 1 :])
            lpar = rest.find("(")
            if lpar != -1:
                rest = trim_inplace_chars(rest[:lpar])
            try:
                out[name] = float(rest)
            except ValueError as e:
                raise RuntimeError(f"Bad solutions.txt line: {line}") from e
    return out


def build_instance_result(
    inst: TsplibInstance,
    rc: int,
    optimal: float | None,
    bucket_1based: int | None,
    ref_map: dict[str, float] | None,
) -> InstanceResult:
    reference: float | None = None
    reference_missing = False
    if ref_map is not None:
        if inst.stem in ref_map:
            reference = ref_map[inst.stem]
        else:
            reference_missing = True
    match_ok: bool | None = None
    if reference is not None and optimal is not None:
        match_ok = math.isclose(optimal, reference, rel_tol=0.0, abs_tol=1e-3)
    elif reference is not None and optimal is None and rc == 0:
        match_ok = False
    return InstanceResult(
        stem=inst.stem,
        n=inst.n,
        returncode=rc,
        optimal=optimal,
        bucket_1based=bucket_1based,
        reference=reference,
        reference_missing=reference_missing,
        match_ok=match_ok,
    )


def instance_cost_passes(r: InstanceResult) -> bool:
    """True only if the solver exited cleanly and the optimum matches solutions.txt."""
    return r.returncode == 0 and r.match_ok is True


def validate_after_run(r: InstanceResult, *, ref_map: dict[str, float] | None) -> str | None:
    """Return human-readable error if this row should abort the suite, else None."""
    if ref_map is None:
        return None
    if r.returncode != 0:
        return None
    if r.reference_missing:
        return f"stem {r.stem} not found in solutions.txt"
    if r.optimal is None:
        return f"exit 0 but no parsable Optimal Solution line for {r.stem}"
    if r.match_ok is False:
        return (
            f"cost mismatch for {r.stem}: concorde_cli reported {r.optimal:g} "
            f"vs solutions.txt reference {r.reference:g}"
        )
    return None


def parse_optimal_from_output(text: str) -> float | None:
    last: float | None = None
    for line in text.splitlines():
        m = _OPTIMAL_SOLUTION_RE.search(line)
        if m:
            try:
                last = float(m.group(1))
            except ValueError:
                pass
    return last


def run_concorde_subprocess(
    concorde_bin: Path,
    seed: int,
    inst: TsplibInstance,
) -> tuple[int, float | None, str]:
    """Run concorde_cli once; return (rc, optimal_or_none, merged_stdout_stderr)."""
    cmd = [str(concorde_bin), "-s", str(seed), str(inst.path)]
    completed = subprocess.run(cmd, capture_output=True, text=True)
    merged = (completed.stdout or "") + (completed.stderr or "")
    optimal = parse_optimal_from_output(merged)
    return completed.returncode, optimal, merged


def _concorde_subprocess_worker(
    args: tuple[Path, int, TsplibInstance],
) -> tuple[int, float | None, str]:
    concorde_bin, seed, inst = args
    return run_concorde_subprocess(concorde_bin, seed, inst)


def _concorde_rc_fail_message(rc: int, inst: TsplibInstance, bucket_1based: int | None) -> str:
    if bucket_1based is not None:
        return f"concorde_cli failed (rc={rc}) bucket={bucket_1based} instance={inst.path}"
    return f"concorde_cli failed (rc={rc}) for {inst.path}"


def run_concorde_batch(
    concorde_bin: Path,
    seed: int,
    work: list[tuple[TsplibInstance, int, int | None]],
    *,
    jobs: int,
    verbose_concorde: bool,
    ref_map: dict[str, float] | None,
    ref_loaded: bool,
    total: int,
) -> tuple[list[InstanceResult], str | None]:
    """
    Run concorde on each instance in ``work`` preserving suite order for logging.

    ``work`` items are ``(inst, global_idx_1based, bucket_1based_or_none)``.

    Returns ``(results, error_message)`` where ``error_message`` is set on the first
    failing case in list order (same early-exit semantics as the sequential script).
    """
    results: list[InstanceResult] = []
    if jobs < 1:
        jobs = 1

    if jobs == 1 or len(work) <= 1:
        for inst, idx_1, bucket in work:
            cmd = [str(concorde_bin), "-s", str(seed), str(inst.path)]
            print(f"==> {' '.join(cmd)}", flush=True)
            rc, optimal, merged = run_concorde_subprocess(concorde_bin, seed, inst)
            if verbose_concorde:
                sys.stdout.write(merged)
                if merged and not merged.endswith("\n"):
                    sys.stdout.write("\n")
                sys.stdout.flush()
            result = build_instance_result(inst, rc, optimal, bucket, ref_map)
            print_instance_footer(result, idx=idx_1, total=total, ref_loaded=ref_loaded)
            results.append(result)
            if rc != 0:
                return results, _concorde_rc_fail_message(rc, inst, bucket)
            msg = validate_after_run(result, ref_map=ref_map)
            if msg:
                return results, msg
        return results, None

    worker_args = [(concorde_bin, seed, inst) for inst, _, _ in work]
    max_workers = min(jobs, len(work))
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        outputs = list(ex.map(_concorde_subprocess_worker, worker_args))

    for (inst, idx_1, bucket), (rc, optimal, merged) in zip(work, outputs):
        cmd = [str(concorde_bin), "-s", str(seed), str(inst.path)]
        print(f"==> {' '.join(cmd)}", flush=True)
        if verbose_concorde:
            with _verbose_print_lock:
                sys.stdout.write(merged)
                if merged and not merged.endswith("\n"):
                    sys.stdout.write("\n")
                sys.stdout.flush()
        result = build_instance_result(inst, rc, optimal, bucket, ref_map)
        print_instance_footer(result, idx=idx_1, total=total, ref_loaded=ref_loaded)
        results.append(result)
        if rc != 0:
            return results, _concorde_rc_fail_message(rc, inst, bucket)
        msg = validate_after_run(result, ref_map=ref_map)
        if msg:
            return results, msg
    return results, None


def print_instance_footer(r: InstanceResult, *, idx: int, total: int, ref_loaded: bool) -> None:
    bpart = f" bucket={r.bucket_1based}/10" if r.bucket_1based is not None else ""
    if r.optimal is not None:
        op = f"{r.optimal:.2f}"
    else:
        op = "?" if r.returncode == 0 else "—"
    if ref_loaded:
        if r.reference_missing:
            ref_s = "MISSING"
        elif r.reference is not None:
            ref_s = f"{r.reference:.2f}"
        else:
            ref_s = "—"
        if r.match_ok is True:
            m_s = "YES"
        elif r.match_ok is False:
            m_s = "NO"
        else:
            m_s = "—"
        ref_part = f" ref={ref_s} match={m_s}"
    else:
        ref_part = ""
    if ref_loaded:
        tag = "PASS" if instance_cost_passes(r) else "FAIL"
    else:
        tag = "EXIT_OK" if r.returncode == 0 else "FAIL"
    print(
        f"<<< {tag} [{idx}/{total}]{bpart} {r.stem} n={r.n} rc={r.returncode} "
        f"concorde={op}{ref_part} >>>",
        flush=True,
    )


def print_final_summary(
    results: list[InstanceResult],
    *,
    label: str,
    solutions_path: Path | None,
    ref_map_loaded: bool,
) -> None:
    n = len(results)
    n_exit_ok = sum(1 for r in results if r.returncode == 0)
    print(f"\n{'=' * 88}", flush=True)
    print(f"concorde_cli summary ({label}): {n_exit_ok}/{n} runs with exit code 0", flush=True)
    if ref_map_loaded and solutions_path is not None:
        print(f"reference costs from: {solutions_path}", flush=True)
    if ref_map_loaded:
        n_match = sum(1 for r in results if r.match_ok is True)
        n_mis = sum(1 for r in results if r.match_ok is False)
        n_miss = sum(1 for r in results if r.reference_missing)
        n_bad = sum(1 for r in results if r.returncode == 0 and r.optimal is None)
        print(
            f"vs solutions.txt: {n_match}/{n} matched; "
            f"mismatch={n_mis}, missing_stem={n_miss}, unparsed_optimal={n_bad}",
            flush=True,
        )
    print(
        f"{'stem':<16} {'n':>5} {'rc':>4} {'concorde':>14} {'reference':>14} {'PASS':>6}",
        flush=True,
    )
    print(
        "(PASS = exit 0 and concorde optimum matches solutions.txt for that stem)",
        flush=True,
    )
    print("-" * 60, flush=True)
    for r in results:
        if r.optimal is not None:
            op = f"{r.optimal:.2f}"
        else:
            op = "—" if r.returncode != 0 else "(unparsed)"
        if not ref_map_loaded:
            ref_col = "—"
            row_tag = "n/c"
        else:
            if r.reference_missing:
                ref_col = "MISSING"
            elif r.reference is not None:
                ref_col = f"{r.reference:.2f}"
            else:
                ref_col = "—"
            row_tag = "PASS" if instance_cost_passes(r) else "FAIL"
        print(f"{r.stem:<16} {r.n:5d} {r.returncode:4d} {op:>14} {ref_col:>14} {row_tag:>6}", flush=True)
    print_overall_verdict(results, solutions_path=solutions_path, ref_map_loaded=ref_map_loaded)
    print(f"{'=' * 88}\n", flush=True)


def print_overall_verdict(
    results: list[InstanceResult],
    *,
    solutions_path: Path | None,
    ref_map_loaded: bool,
) -> None:
    n = len(results)
    n_exit_ok = sum(1 for r in results if r.returncode == 0)
    exit_ok = n_exit_ok == n
    if ref_map_loaded:
        n_pass = sum(1 for r in results if instance_cost_passes(r))
        all_pass = n_pass == n
        verdict = "PASS" if all_pass else "FAIL"
    else:
        all_pass = False
        verdict = "INCOMPLETE" if exit_ok else "FAIL"
    print(f"\nOVERALL: {verdict}", flush=True)
    print(f"  processes: {n_exit_ok}/{n} exited with code 0", flush=True)
    if ref_map_loaded and solutions_path is not None:
        n_pass = sum(1 for r in results if instance_cost_passes(r))
        print(f"  solutions.txt: {solutions_path}", flush=True)
        print(
            f"  PASS (cost match): {n_pass}/{n} cases (PASS requires optimum == reference, within 1e-3)",
            flush=True,
        )
        if not all_pass:
            print("  at least one case is not PASS (see table)", flush=True)
    elif not ref_map_loaded:
        print(
            "  PASS/cost: not evaluated (--skip-solutions-check or missing file). "
            "A case PASS requires matching solutions.txt.",
            flush=True,
        )


def discover_sort_by_dimension(tsp_dir: Path) -> list[TsplibInstance]:
    out: list[TsplibInstance] = []
    for ent in tsp_dir.iterdir():
        if not ent.is_file():
            continue
        if ent.suffix not in (".tsp", ".TSP"):
            continue
        stem = stem_upper(ent)
        n = read_symmetric_tsplib_dimension(ent)
        out.append(TsplibInstance(path=ent.resolve(), stem=stem, n=n))
    out.sort(key=lambda x: (x.n, x.stem))
    return out


def partition_buckets(sorted_list: list[TsplibInstance], num_buckets: int) -> list[list[TsplibInstance]]:
    total = len(sorted_list)
    base = total // num_buckets
    rem = total % num_buckets
    buckets: list[list[TsplibInstance]] = [[] for _ in range(num_buckets)]
    idx = 0
    for b in range(num_buckets):
        sz = base + (1 if b < rem else 0)
        for _ in range(sz):
            buckets[b].append(sorted_list[idx])
            idx += 1
    return buckets


def print_buckets(buckets: list[list[TsplibInstance]]) -> None:
    for b, bucket in enumerate(buckets):
        stems = ", ".join(inst.stem for inst in bucket)
        print(f"bucket {b + 1}/10 (n={len(bucket)}): {stems}")


def main() -> int:
    script_dir = Path(__file__).resolve().parent
    default_tsp = script_dir.parent / "data" / "tsp"

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tsp-dir", type=Path, default=default_tsp, help="Directory of .tsp fixtures")
    p.add_argument(
        "--concorde-bin",
        type=Path,
        default=None,
        help="concorde_cli executable (default: CONCORDE_BIN env)",
    )
    p.add_argument("--seed", type=int, default=1, help="RNG seed passed as -s to concorde_cli")
    p.add_argument(
        "--list-only",
        action="store_true",
        help="Print bucket layout (stems) and exit without running the solver",
    )
    p.add_argument(
        "--quick",
        action="store_true",
        help="Run only att48, eil51, pr76 (fast local check; skips 110-instance assertion)",
    )
    p.add_argument(
        "--bucket",
        type=int,
        choices=range(1, 11),
        default=None,
        metavar="N",
        help="Run only bucket N (1–10), same numbering as test_concorde_symmetric_suite.cpp. "
        "Default: run every bucket in order.",
    )
    p.add_argument(
        "--quiet-concorde",
        action="store_true",
        help="Do not print concorde_cli stdout/stderr (cost is still parsed when an "
        '"Optimal Solution:" line is present).',
    )
    p.add_argument(
        "--solutions-txt",
        type=Path,
        default=None,
        help="Reference tour costs (STEM : value lines, like the C++ tests). "
        "Default: <tsp-dir>/solutions.txt",
    )
    p.add_argument(
        "--skip-solutions-check",
        action="store_true",
        help="Do not load solutions.txt; rows show EXIT_OK/n/c, never PASS vs reference. "
        "OVERALL cannot be PASS (cost verification disabled).",
    )
    p.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=1,
        metavar="N",
        help="Run up to N concorde_cli processes at once (default: 1). "
        "Log order matches the suite. With N>1, prefer --quiet-concorde.",
    )
    args = p.parse_args()

    if args.jobs < 1:
        print("run_concorde_cli_tsp_buckets.py: --jobs must be >= 1", file=sys.stderr)
        return 2

    if args.quick and args.bucket is not None:
        print(
            "run_concorde_cli_tsp_buckets.py: --quick and --bucket are mutually exclusive",
            file=sys.stderr,
        )
        return 2

    tsp_dir = args.tsp_dir.resolve()
    if not tsp_dir.is_dir():
        print(f"run_concorde_cli_tsp_buckets.py: TSP_DIR is not a directory: {tsp_dir}", file=sys.stderr)
        return 2

    def resolve_concorde_bin() -> Path | None:
        if args.concorde_bin is not None:
            return Path(args.concorde_bin)
        env = os.environ.get("CONCORDE_BIN")
        if not env:
            return None
        return Path(env)

    def require_concorde_bin() -> tuple[Path | None, int]:
        concorde_bin = resolve_concorde_bin()
        if concorde_bin is None:
            print("run_concorde_cli_tsp_buckets.py: set CONCORDE_BIN or pass --concorde-bin", file=sys.stderr)
            return None, 2
        cb = str(concorde_bin)
        if not concorde_bin.is_file() or not os.access(cb, os.X_OK):
            print(f"run_concorde_cli_tsp_buckets.py: not an executable: {concorde_bin}", file=sys.stderr)
            return None, 2
        return concorde_bin, 0

    sorted_list = discover_sort_by_dimension(tsp_dir)
    if len(sorted_list) != 110:
        print(
            f"run_concorde_cli_tsp_buckets.py: expected 110 .tsp instances, found {len(sorted_list)}",
            file=sys.stderr,
        )
        return 1

    buckets = partition_buckets(sorted_list, 10)
    if args.list_only:
        if args.bucket is not None:
            b = args.bucket - 1
            bucket = buckets[b]
            stems = ", ".join(inst.stem for inst in bucket)
            print(f"bucket {args.bucket}/10 (n={len(bucket)}): {stems}")
        else:
            print_buckets(buckets)
        return 0

    solutions_path_resolved: Path | None = None
    ref_map: dict[str, float] | None = None
    if not args.skip_solutions_check:
        solutions_path_resolved = (
            args.solutions_txt.resolve() if args.solutions_txt is not None else (tsp_dir / "solutions.txt").resolve()
        )
        if not solutions_path_resolved.is_file():
            print(
                f"run_concorde_cli_tsp_buckets.py: missing solutions.txt: {solutions_path_resolved} "
                "(use --skip-solutions-check for exit-code-only smoke)",
                file=sys.stderr,
            )
            return 1
        try:
            ref_map = load_tsplib_solution_costs(solutions_path_resolved)
        except RuntimeError as e:
            print(f"run_concorde_cli_tsp_buckets.py: {e}", file=sys.stderr)
            return 1
    ref_loaded = ref_map is not None

    verbose_concorde = not args.quiet_concorde

    if args.quick:
        concorde_bin, err = require_concorde_bin()
        if err != 0:
            return err
        quick_names = ("att48.tsp", "eil51.tsp", "pr76.tsp")
        total = len(quick_names)
        work: list[tuple[TsplibInstance, int, int | None]] = []
        for i, name in enumerate(quick_names, start=1):
            path = tsp_dir / name
            if not path.is_file():
                print(f"run_concorde_cli_tsp_buckets.py: missing fixture {path}", file=sys.stderr)
                return 1
            n = read_symmetric_tsplib_dimension(path)
            inst = TsplibInstance(path=path.resolve(), stem=stem_upper(path), n=n)
            work.append((inst, i, None))
        results, batch_err = run_concorde_batch(
            concorde_bin,
            args.seed,
            work,
            jobs=args.jobs,
            verbose_concorde=verbose_concorde,
            ref_map=ref_map,
            ref_loaded=ref_loaded,
            total=total,
        )
        if batch_err:
            print(f"run_concorde_cli_tsp_buckets.py: {batch_err}", file=sys.stderr)
            print_final_summary(
                results,
                label="quick (partial)",
                solutions_path=solutions_path_resolved,
                ref_map_loaded=ref_loaded,
            )
            return 1
        print_final_summary(
            results,
            label="quick",
            solutions_path=solutions_path_resolved,
            ref_map_loaded=ref_loaded,
        )
        print("run_concorde_cli_tsp_buckets.py: OK (quick)", flush=True)
        return 0

    concorde_bin, err = require_concorde_bin()
    if err != 0:
        return err

    results = []
    idx = 0

    if args.bucket is not None:
        b = args.bucket - 1
        bucket = buckets[b]
        total = len(bucket)
        print(f"--- bucket {args.bucket}/10 ({total} instances) ---", flush=True)
        work = [(inst, i, args.bucket) for i, inst in enumerate(bucket, start=1)]
        batch_results, batch_err = run_concorde_batch(
            concorde_bin,
            args.seed,
            work,
            jobs=args.jobs,
            verbose_concorde=verbose_concorde,
            ref_map=ref_map,
            ref_loaded=ref_loaded,
            total=total,
        )
        results.extend(batch_results)
        if batch_err:
            print(f"run_concorde_cli_tsp_buckets.py: {batch_err}", file=sys.stderr)
            print_final_summary(
                results,
                label=f"bucket {args.bucket}/10 (partial)",
                solutions_path=solutions_path_resolved,
                ref_map_loaded=ref_loaded,
            )
            return 1
        print_final_summary(
            results,
            label=f"bucket {args.bucket}/10",
            solutions_path=solutions_path_resolved,
            ref_map_loaded=ref_loaded,
        )
    else:
        total = 110
        for b, bucket in enumerate(buckets):
            print(f"--- bucket {b + 1}/10 ({len(bucket)} instances) ---", flush=True)
            work = []
            for inst in bucket:
                idx += 1
                work.append((inst, idx, b + 1))
            batch_results, batch_err = run_concorde_batch(
                concorde_bin,
                args.seed,
                work,
                jobs=args.jobs,
                verbose_concorde=verbose_concorde,
                ref_map=ref_map,
                ref_loaded=ref_loaded,
                total=total,
            )
            results.extend(batch_results)
            if batch_err:
                print(f"run_concorde_cli_tsp_buckets.py: {batch_err}", file=sys.stderr)
                print_final_summary(
                    results,
                    label="all buckets (partial)",
                    solutions_path=solutions_path_resolved,
                    ref_map_loaded=ref_loaded,
                )
                return 1
        print_final_summary(
            results,
            label="all buckets",
            solutions_path=solutions_path_resolved,
            ref_map_loaded=ref_loaded,
        )

    if ref_loaded:
        print(
            "run_concorde_cli_tsp_buckets.py: OK — every case PASS (optimum matches solutions.txt)",
            flush=True,
        )
    else:
        print(
            "run_concorde_cli_tsp_buckets.py: OK (concorde exit 0 for each run; "
            "PASS vs solutions.txt was not evaluated)",
            flush=True,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
