#!/usr/bin/env python3
"""
Python-only benchmark for the spectral SIS threshold study (research_proposal.pdf).

Replicates the NetLogo grp-SIS dynamics on a fixed Erdős–Rényi graph:
same p = avg_degree/(n-1), same recovery draws (exponential / Tang power law / lognormal),
same metrics (extinction, late-window prevalence).

Usage (from project root):
  python benchmark/run_benchmark.py --preset quick
  python benchmark/run_benchmark.py --preset er_baseline --jobs 4

Outputs under benchmark/output/: raw runs, survival curves, threshold summary, spectral ref.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

# Allow `python benchmark/run_benchmark.py` without installing package
_BENCH = Path(__file__).resolve().parent
if str(_BENCH) not in sys.path:
    sys.path.insert(0, str(_BENCH))

from analysis import summarize_regime, survival_curve  # noqa: E402
from graph import build_er_netlogo_order  # noqa: E402
from sis_engine import mean_degree  # noqa: E402
from spectral import largest_eigenvalue_adjacency, spectral_reference  # noqa: E402
from worker import run_single_task  # noqa: E402

PRESETS: dict[str, dict[str, Any]] = {
    "quick": {
        "num_nodes": 2000,
        "avg_degree": 6.0,
        "graph_seed": 10001,
        "recovery_mean": 5.0,
        "initial_infected": 5,
        "max_ticks": 3000,
        "late_window_len": 200,
        "repetitions": 4,
        "beta_grid": np.arange(0.02, 0.57, 0.03),
        "power_law_lambda": 4.24,
        "lognormal_sigma": 1.0,
    },
    "er_baseline": {
        "num_nodes": 2000,
        "avg_degree": 6.0,
        "graph_seed": 10001,
        "recovery_mean": 5.0,
        "initial_infected": 5,
        "max_ticks": 10_000,
        "late_window_len": 200,
        "repetitions": 24,
        "beta_grid": np.arange(0.018, 0.048 + 1e-9, 0.002),
        "power_law_lambda": 4.24,
        "lognormal_sigma": 1.0,
    },
}

REGIME_NAMES: dict[str, str] = {
    "exponential": "exponential",
    "power_law_tang": "power_law_tang",
    "lognormal_tang": "lognormal_tang",
}


def run_experiment(
    preset: dict[str, Any],
    output_dir: Path,
    jobs: int,
    regimes: list[str],
    use_numba: bool,
) -> None:
    rng_graph = np.random.default_rng(int(preset["graph_seed"]))
    g = build_er_netlogo_order(int(preset["num_nodes"]), float(preset["avg_degree"]), rng_graph)
    lam = largest_eigenvalue_adjacency(g)
    spec = spectral_reference(float(preset["recovery_mean"]), lam)
    k_mean = mean_degree(g)

    meta = {
        "num_nodes": g.n,
        "num_edges": len(g.edges),
        "k_mean": k_mean,
        "graph_seed": int(preset["graph_seed"]),
        **spec,
        "recovery_mean": preset["recovery_mean"],
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "spectral_reference.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

    beta_grid = np.asarray(preset["beta_grid"], dtype=float)
    reps = int(preset["repetitions"])
    neigh_data = g.neighbors

    tasks: list[tuple] = []
    run_counter = 0
    for regime_key in regimes:
        regime = REGIME_NAMES[regime_key]
        for beta in beta_grid:
            for _ in range(reps):
                run_seed = int(preset["graph_seed"]) * 1_000_000 + run_counter
                run_counter += 1
                tasks.append(
                    (
                        neigh_data,
                        g.n,
                        float(beta),
                        float(preset["recovery_mean"]),
                        regime,
                        run_seed,
                        int(preset["initial_infected"]),
                        int(preset["max_ticks"]),
                        int(preset["late_window_len"]),
                        float(preset["power_law_lambda"]),
                        float(preset["lognormal_sigma"]),
                        bool(use_numba),
                    )
                )

    rows: list[dict] = []
    if jobs <= 1:
        for t in tasks:
            rows.append(run_single_task(t))
    else:
        with ProcessPoolExecutor(max_workers=jobs) as ex:
            futs = {ex.submit(run_single_task, t): i for i, t in enumerate(tasks)}
            for fut in as_completed(futs):
                rows.append(fut.result())

    df = pd.DataFrame(rows)
    df["expt_seed"] = int(preset["graph_seed"])
    df["network_type"] = "ER"
    df["interconnection_structure"] = "other"

    raw_path = output_dir / "raw_runs.csv"
    df.to_csv(raw_path, index=False)

    # Survival curves + threshold summary per regime
    summaries = []
    inf_col = "infection_prob"
    ext_col = "bs_out_extinct"
    for regime_key in regimes:
        regime = REGIME_NAMES[regime_key]
        sub = df[df["recovery_regime"] == regime]
        if sub.empty:
            continue
        curve = survival_curve(sub, inf_col, ext_col, prev_col="bs_out_late_mean_prevalence", tick_col="bs_out_final_tick")
        curve_path = output_dir / f"survival_curve_{regime_key}.csv"
        curve.to_csv(curve_path, index=False)

        summ = summarize_regime(sub, float(spec["beta_c"]), n_boot=400)
        summ["recovery_regime"] = regime_key
        summ["tau_c_pred"] = spec["tau_c"]
        summ["beta_c_pred"] = spec["beta_c"]
        summ["lambda_max"] = lam
        summaries.append(summ)

    pd.DataFrame(summaries).to_csv(output_dir / "threshold_summary.csv", index=False)

    # Optional plot
    try:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(8, 5))
        for regime_key in regimes:
            regime = REGIME_NAMES[regime_key]
            sub = df[df["recovery_regime"] == regime]
            if sub.empty:
                continue
            curve = survival_curve(sub, inf_col, ext_col)
            ax.plot(curve[inf_col], curve["p_survive"], marker="o", ms=3, label=regime_key)
        ax.axvline(spec["beta_c"], color="k", ls="--", lw=1, label=r"$\beta_c \approx 1/(\lambda_{\max} E[W])$")
        ax.set_xlabel("infection probability")
        ax.set_ylabel("P(survive)")
        ax.legend()
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(output_dir / "p_survive_vs_beta.png", dpi=150)
        plt.close(fig)
    except Exception:  # noqa: BLE001
        pass


def main() -> None:
    p = argparse.ArgumentParser(description="Python SIS spectral benchmark")
    p.add_argument("--preset", choices=list(PRESETS.keys()), default="quick")
    p.add_argument("--output-dir", type=Path, default=None)
    p.add_argument(
        "--jobs",
        type=int,
        default=max(1, min(8, (os.cpu_count() or 4))),
        help="Parallel workers (default: up to 8)",
    )
    p.add_argument(
        "--regimes",
        nargs="*",
        default=["exponential", "power_law_tang", "lognormal_tang"],
        help="Subsets of exponential, power_law_tang, lognormal_tang",
    )
    p.add_argument(
        "--no-numba",
        action="store_true",
        help="Run pure Python sis_engine (slow; bit-for-bit with reference implementation)",
    )
    args = p.parse_args()

    out = args.output_dir or (_BENCH / "output")
    preset = PRESETS[args.preset].copy()
    run_experiment(
        preset,
        out.resolve(),
        jobs=max(1, args.jobs),
        regimes=list(args.regimes),
        use_numba=not args.no_numba,
    )
    print("Wrote outputs to", out.resolve())


if __name__ == "__main__":
    main()
