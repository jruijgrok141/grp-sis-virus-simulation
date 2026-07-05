#!/usr/bin/env python3
"""
Run NetLogo 7 BehaviorSpace experiments headlessly, then compute lambda_max per edge list.

Requires NetLogo 7.x with NetLogo_Console. Set NETLOGO_HOME to the installation folder
(e.g. C:\\Program Files\\NetLogo 7.0.3) or pass --netlogo-home.

Usage (from project root):
  python scripts/run_pipeline.py --all
  python scripts/run_pipeline.py --experiment 01-export-networks
  python scripts/run_pipeline.py --experiment 08-dynamics-trajectory-ER  # per-tick prevalence CSV (see notebooks/virus_dynamics_prevalence.ipynb)

GUI: Tools → BehaviorSpace → experiments are stored inside `virus_simulation.nlogox`
(same definitions as `netlogo/behaviorspace_experiments.xml`; after editing the XML,
re-copy `<experiments>…</experiments>` into the `.nlogox`, stripping `<!-- -->` comments
that contain `--`, which NetLogo rejects).
NetLogo 7 uses the XML loader (`LabXMLLoader`): wrap reporters in `<metrics>` and
parameter sets in `<constants>`. Chooser values need NetLogo string literals in XML
(`&quot;...&quot;`).
Headless runs do **not** pass `--setup-file`: NetLogo 7.0.3 then fails with
`java.util.NoSuchElementException: head of empty list`. Use `--behaviorspace-setup`
only if you use a model without embedded experiments.
`--all` runs the Erdos-Renyi-focused pipeline: \texttt{01} (export \texttt{edges_ER_*.csv} only),
\texttt{05}--\texttt{07} (threshold sweeps, three seeds), and \texttt{08} (per-tick trajectory,
seed \texttt{10001}). Lattice/Ring experiments \texttt{02}--\texttt{04} are optional
(\texttt{--with-lattice-ring}). After aggregation, \texttt{write_report_macros.py} updates
\texttt{report/generated_quantities.tex}. Each recovery law uses about $24\times 20\times 3$
BehaviorSpace replicates on \texttt{05}--\texttt{07}.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd

# Allow importing analysis_utils when run as script
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.analysis_utils import compute_lambda_max_for_edge_file  # noqa: E402


def resolve_model_and_setup(root: Path, model_override: Path | None) -> tuple[Path, Path | None]:
    """Prefer netlogo/ then NetLogo/ (cross-platform)."""
    if model_override is not None:
        m = model_override.resolve()
        if not m.is_file():
            raise SystemExit(f"Model not found: {m}")
        sub = m.parent
        x = sub / "behaviorspace_experiments.xml"
        return m, x if x.is_file() else None
    for sub in ("netlogo", "NetLogo"):
        m = root / sub / "virus_simulation.nlogox"
        if m.is_file():
            x = root / sub / "behaviorspace_experiments.xml"
            return m, x if x.is_file() else None
    raise SystemExit(
        f"No virus_simulation.nlogox under {root / 'netlogo'} or {root / 'NetLogo'}"
    )


def find_netlogo_console(home: Path | None) -> Path:
    if home is None:
        home_s = os.environ.get("NETLOGO_HOME", "").strip()
        home = Path(home_s) if home_s else None
    if home is None or not home.is_dir():
        raise SystemExit(
            "Set NETLOGO_HOME to your NetLogo installation directory "
            "(folder containing NetLogo_Console.exe), or pass --netlogo-home."
        )
    exe = home / "NetLogo_Console.exe"
    if not exe.is_file():
        exe = home / "app" / "NetLogo_Console.exe"
    if not exe.is_file():
        raise SystemExit(f"NetLogo_Console.exe not found under {home}")
    return exe


def run_experiment(
    console: Path,
    model: Path,
    experiment_name: str,
    table_out: Path,
    threads: int,
    cwd: Path,
    setup_file: Path | None,
) -> None:
    table_out.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(console),
        "--headless",
        "--model",
        str(model.resolve()),
    ]
    if setup_file is not None and setup_file.is_file():
        cmd.extend(["--setup-file", str(setup_file.resolve())])
    cmd.extend(
        [
            "--experiment",
            experiment_name,
            "--table",
            str(table_out.resolve()),
            "--threads",
            str(threads),
        ]
    )
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True, cwd=str(cwd))


def build_lambda_table(edges_dir: Path, out_csv: Path) -> pd.DataFrame:
    rows = []
    for p in sorted(edges_dir.glob("edges_*.csv")):
        try:
            rows.append(compute_lambda_max_for_edge_file(p))
        except Exception as e:  # noqa: BLE001
            print(f"Skip {p}: {e}")
    if not rows:
        print("No edge files found in", edges_dir)
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    try:
        df.to_csv(out_csv, index=False)
    except PermissionError as e:
        raise SystemExit(
            f"Cannot write {out_csv}: {e}\n"
            "Close the file if it is open in Excel, Cursor, or another program, then retry."
        ) from e
    print("Wrote", out_csv)
    return df


NOTEBOOK_OUTPUT_CSVS = (
    "empirical_thresholds.csv",
    "threshold_ratio_by_run.csv",
    "threshold_ratio_summary.csv",
)


def clear_pipeline_csv_outputs(root: Path, experiments: list[str], full: bool) -> None:
    """Remove CSV files from prior pipeline runs and stale notebook exports before new outputs."""
    out = root / "output"
    removed = 0
    if full:
        for sub in ("raw", "edges"):
            d = out / sub
            if d.is_dir():
                for p in d.glob("*.csv"):
                    try:
                        p.unlink()
                        removed += 1
                    except OSError as e:
                        print(f"Warning: could not remove {p}: {e}")
        lm = out / "lambda_max.csv"
        if lm.is_file():
            try:
                lm.unlink()
                removed += 1
            except OSError as e:
                print(f"Warning: could not remove {lm}: {e}")
    else:
        out_raw = out / "raw"
        for name in experiments:
            p = out_raw / f"{name.replace(' ', '_')}_table.csv"
            if p.is_file():
                try:
                    p.unlink()
                    removed += 1
                except OSError as e:
                    print(f"Warning: could not remove {p}: {e}")
        if "01-export-networks" in experiments:
            edges_dir = out / "edges"
            if edges_dir.is_dir():
                for p in edges_dir.glob("*.csv"):
                    try:
                        p.unlink()
                        removed += 1
                    except OSError as e:
                        print(f"Warning: could not remove {p}: {e}")
            lm = out / "lambda_max.csv"
            if lm.is_file():
                try:
                    lm.unlink()
                    removed += 1
                except OSError as e:
                    print(f"Warning: could not remove {lm}: {e}")
    for fname in NOTEBOOK_OUTPUT_CSVS:
        p = out / fname
        if p.is_file():
            try:
                p.unlink()
                removed += 1
            except OSError as e:
                print(f"Warning: could not remove {p}: {e}")
    if removed:
        print(f"Removed {removed} existing CSV file(s) under {out}")


def main() -> None:
    parser = argparse.ArgumentParser(description="NetLogo BehaviorSpace + lambda_max CSV pipeline")
    parser.add_argument(
        "--project-root",
        type=Path,
        default=ROOT,
        help="Project root (default: parent of scripts/)",
    )
    parser.add_argument(
        "--netlogo-home",
        type=Path,
        default=None,
        help="NetLogo install dir (else env NETLOGO_HOME)",
    )
    parser.add_argument(
        "--model",
        type=Path,
        default=None,
        help="Path to virus_simulation.nlogox (default: NetLogo/virus_simulation.nlogox)",
    )
    parser.add_argument("--threads", type=int, default=max(1, (os.cpu_count() or 4) * 3 // 4))
    parser.add_argument(
        "--experiment",
        type=str,
        default=None,
        help="Single BehaviorSpace experiment name to run",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Run ER-focused pipeline: 01 (ER edge export), 05–07, 08 + lambda_max + report macros",
    )
    parser.add_argument(
        "--with-lattice-ring",
        action="store_true",
        help="With --all: also run 02–04 (Lattice4/Ring coarse sweeps) before the ER experiments",
    )
    parser.add_argument(
        "--no-report-macros",
        action="store_true",
        help="Skip writing report/generated_quantities.tex after aggregation",
    )
    parser.add_argument(
        "--behaviorspace-setup",
        type=Path,
        default=None,
        help="Optional path to behaviorspace XML for NetLogo --setup-file (default: omit; required setup is embedded in the .nlogox)",
    )
    args = parser.parse_args()

    root = args.project_root.resolve()
    model, _bs_xml = resolve_model_and_setup(root, args.model)
    setup_file: Path | None = None
    if args.behaviorspace_setup is not None:
        setup_file = args.behaviorspace_setup.resolve()
        if not setup_file.is_file():
            raise SystemExit(f"BehaviorSpace setup file not found: {setup_file}")

    out_edges = root / "output" / "edges"
    out_raw = root / "output" / "raw"
    out_edges.mkdir(parents=True, exist_ok=True)
    out_raw.mkdir(parents=True, exist_ok=True)

    console = find_netlogo_console(args.netlogo_home)
    print("Using NetLogo console:", console)
    run_cwd = root

    experiments_er_full = [
        "01-export-networks",
        "05-baseline-empirical-threshold-ER",
        "06-baseline-ER-power-law-Tang",
        "07-baseline-ER-lognormal-Tang",
        "08-dynamics-trajectory-ER",
    ]
    experiments_lattice_ring = [
        "02-threshold-exponential",
        "03-threshold-power-law-Tang",
        "04-threshold-lognormal-Tang",
    ]

    if args.all:
        to_run = list(experiments_er_full)
        if args.with_lattice_ring:
            to_run = [experiments_er_full[0]] + experiments_lattice_ring + experiments_er_full[1:]
    elif args.experiment:
        to_run = [args.experiment]
    else:
        parser.print_help()
        raise SystemExit("Pass --all or --experiment NAME")

    clear_pipeline_csv_outputs(root, to_run, full=args.all)

    for name in to_run:
        table_path = out_raw / f"{name.replace(' ', '_')}_table.csv"
        run_experiment(console, model, name, table_path, args.threads, run_cwd, setup_file)

    if "01-export-networks" in to_run:
        build_lambda_table(out_edges, root / "output" / "lambda_max.csv")

    raw_dir = root / "output" / "raw"
    lambda_csv = root / "output" / "lambda_max.csv"
    has_threshold_tables = any(
        p.name.endswith("_table.csv") and "01-export" not in p.name and "export-networks" not in p.name
        for p in raw_dir.glob("*_table.csv")
    )
    if has_threshold_tables and lambda_csv.is_file():
        try:
            from scripts.aggregate_thresholds import build_merged, write_threshold_outputs

            merged, cols = build_merged(raw_dir, lambda_csv)
            write_threshold_outputs(merged, cols, root / "output")
        except Exception as e:  # noqa: BLE001
            print("Warning: aggregate_thresholds failed:", e)
        else:
            if not args.no_report_macros:
                try:
                    subprocess.run(
                        [sys.executable, str(ROOT / "scripts" / "write_report_macros.py")],
                        check=False,
                        cwd=str(root),
                    )
                except Exception as e:  # noqa: BLE001
                    print("Warning: write_report_macros failed:", e)

    print("Done.")


if __name__ == "__main__":
    main()
