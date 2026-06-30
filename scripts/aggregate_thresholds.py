#!/usr/bin/env python3
"""
Build empirical_thresholds.csv, threshold_ratio_by_run.csv, threshold_ratio_summary.csv,
and optionally threshold_ratio_bar.png from BehaviorSpace raw tables + lambda_max.csv.

Tables from experiments 02--04 (lattice/ring) are skipped so aggregation stays ER-focused even
if those CSVs remain under output/raw/.

Called automatically at the end of run_pipeline.py --all. Standalone: python scripts/aggregate_thresholds.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.analysis_utils import (  # noqa: E402
    net_label_from_behaviorspace,
    pick_column,
    read_behaviorspace_table,
)
from scripts.threshold_estimators import row_threshold_summary  # noqa: E402


def recovery_regime_from_path(path: Path) -> str:
    n = path.stem.lower()
    if "05-baseline" in n or "baseline-empirical-threshold-er" in n:
        return "exponential"
    if "06-baseline" in n or "power-law-tang" in n or "power_law" in n:
        return "power law (Tang)"
    if "07-baseline" in n or ("lognormal" in n and "tang" in n):
        return "lognormal (Tang)"
    if "exponential" in n:
        return "exponential"
    if "power-law" in n or "power_law" in n:
        return "power law (Tang)"
    if "lognormal" in n:
        return "lognormal (Tang)"
    if "weibull" in n:
        return "weibull (k<1)"
    if "gamma-heavy" in n or "05-threshold-gamma-heavy" in n:
        return "gamma (a<1)"
    if "gamma" in n:
        return "gamma Erlang k=3"
    return path.stem


def load_threshold_tables(raw_dir: Path) -> pd.DataFrame:
    frames = []
    for path in sorted(raw_dir.glob("*_table.csv")):
        if "01-export" in path.name:
            continue
        # Per-tick trajectory exports (not one row per run)
        if "08-dynamics-trajectory" in path.name:
            continue
        # Default ER-focused pipeline: ignore optional lattice/ring coarse sweeps if present
        if path.name.startswith(("02-", "03-", "04-")):
            continue
        df = read_behaviorspace_table(path)
        df["recovery_regime"] = recovery_regime_from_path(path)
        df["source_file"] = path.name
        frames.append(df)
    if not frames:
        raise FileNotFoundError(f"No threshold tables in {raw_dir}")
    return pd.concat(frames, ignore_index=True)


def empirical_threshold_legacy(
    sub_df: pd.DataFrame,
    inf_col: str,
    ext_col: str,
    prev_col: str,
    q: float = 0.5,
) -> float:
    if sub_df.empty or inf_col not in sub_df.columns:
        return float("nan")
    rows = []
    for beta, g in sub_df.groupby(inf_col, dropna=False):
        survived = (g[ext_col].astype(float) < 0.5) & (
            pd.to_numeric(g[prev_col], errors="coerce").fillna(0.0) > 1e-4
        )
        rows.append({inf_col: beta, "p_persist": float(survived.mean())})
    if not rows:
        return float("nan")
    g2 = pd.DataFrame(rows).sort_values(inf_col, na_position="last")
    above = g2[g2["p_persist"] >= q]
    if above.empty:
        mx = g2[inf_col].max()
        return float(mx) if pd.notna(mx) else float("nan")
    return float(above.iloc[0][inf_col])


def build_merged(raw_dir: Path, lambda_path: Path) -> tuple[pd.DataFrame, dict[str, str]]:
    lam = pd.read_csv(lambda_path)
    lam["net_label"] = lam["net_label"].astype(str)
    lam["random_seed"] = lam["random_seed"].astype(int)

    runs = load_threshold_tables(raw_dir)
    net_col = pick_column(
        runs,
        "interconnection_structure",
        "interconnection-structure",
        "network_type",
        "network-type",
    )
    type_col = pick_column(runs, "network_type", "network-type")
    seed_col = pick_column(
        runs, "expt_seed", "expt-seed", "run_seed", "run-seed", "random_seed", "random-seed"
    )
    inf_col = pick_column(
        runs, "infection_prob", "infection-prob", "infection prob", "[infection-prob]"
    )
    ext_col = pick_column(runs, "bs_out_extinct", "bs-out-extinct")
    prev_col = pick_column(runs, "bs_out_late_mean_prevalence", "bs-out-late-mean-prevalence")
    tick_col = pick_column(runs, "bs_out_final_tick", "bs-out-final-tick")
    rm_col = pick_column(runs, "recovery_mean", "recovery-mean")

    missing = [
        n
        for n, c in [
            ("network column", net_col),
            ("random_seed", seed_col),
            ("infection_prob", inf_col),
            ("extinct metric", ext_col),
            ("late prevalence", prev_col),
        ]
        if c is None
    ]
    if missing:
        raise KeyError(f"Could not resolve columns: {missing}. Available: {list(runs.columns)}")

    runs["net_label"] = runs.apply(
        lambda row: net_label_from_behaviorspace(str(row[net_col]), row.get(type_col) if type_col else None),
        axis=1,
    )
    runs["random_seed"] = runs[seed_col].astype(int)

    _lam_merge = ["net_label", "random_seed", "lambda_max"]
    for _c in ("k_mean", "k2_mean", "tau_pred_homogeneous_mf", "tau_pred_heterogeneous_mf"):
        if _c in lam.columns:
            _lam_merge.append(_c)
    merged = runs.merge(lam[_lam_merge], on=["net_label", "random_seed"], how="left")
    merged[inf_col] = pd.to_numeric(merged[inf_col], errors="coerce")
    ew = float(merged[rm_col].dropna().iloc[0]) if rm_col else 5.0
    merged["beta_pred"] = 1.0 / (merged["lambda_max"] * ew)
    merged["tau_pred_spectral"] = merged["beta_pred"] * ew
    merged["tau_pred"] = merged["tau_pred_spectral"]
    if "tau_pred_homogeneous_mf" in merged.columns:
        merged["beta_pred_homogeneous_mf"] = merged["tau_pred_homogeneous_mf"] / ew
    if "tau_pred_heterogeneous_mf" in merged.columns:
        merged["beta_pred_heterogeneous_mf"] = merged["tau_pred_heterogeneous_mf"] / ew

    _n = len(merged)
    merged = merged.dropna(subset=["lambda_max"])
    if len(merged) < _n:
        print(
            "aggregate_thresholds: dropped",
            _n - len(merged),
            "rows with no lambda_max match (use the same expt-seed as in 01-export-networks, e.g. 10001).",
        )

    cols = {
        "inf_col": inf_col,
        "ext_col": ext_col,
        "prev_col": prev_col,
        "tick_col": tick_col,
    }
    return merged, cols


def write_threshold_outputs(
    merged: pd.DataFrame,
    cols: dict[str, str],
    out_dir: Path,
    *,
    n_boot: int = 400,
    plot: bool = True,
) -> None:
    inf_col = cols["inf_col"]
    ext_col = cols["ext_col"]
    prev_col = cols["prev_col"]
    tick_col = cols["tick_col"]

    rows_out = []
    for (reg, net, seed), sub in merged.groupby(["recovery_regime", "net_label", "random_seed"]):
        beta_p = float(sub["beta_pred"].iloc[0])
        if not np.isfinite(beta_p):
            beta_p = float("nan")
        summ = row_threshold_summary(
            sub,
            inf_col,
            ext_col,
            prev_col,
            tick_col,
            beta_pred=beta_p,
            n_boot=n_boot,
        )
        summ["recovery_regime"] = reg
        summ["net_label"] = net
        summ["random_seed"] = seed
        summ["beta_hat_50_legacy_persist"] = empirical_threshold_legacy(
            sub, inf_col, ext_col, prev_col
        )
        summ["ratio_legacy_over_pred"] = (
            summ["beta_hat_50_legacy_persist"] / beta_p if beta_p and np.isfinite(beta_p) else float("nan")
        )
        rows_out.append(summ)

    thresh_tbl = pd.DataFrame(rows_out)
    order = [
        "recovery_regime",
        "net_label",
        "random_seed",
        "beta_pred",
        "beta_hat_surv_50",
        "beta_hat_surv_25",
        "beta_hat_surv_50_median_boot",
        "beta_hat_surv_50_ci_low",
        "beta_hat_surv_50_ci_high",
        "beta_hat_logit_ed50_surv",
        "beta_hat_persist_50",
        "beta_hat_persist_relaxed_50",
        "median_tick_if_extinct",
        "median_tick_if_survive",
        "ratio_surv_50_over_pred",
        "beta_hat_50_legacy_persist",
        "ratio_legacy_over_pred",
    ]
    thresh_tbl = thresh_tbl[[c for c in order if c in thresh_tbl.columns]]
    out_dir.mkdir(parents=True, exist_ok=True)
    emp_path = out_dir / "empirical_thresholds.csv"
    thresh_tbl.to_csv(emp_path, index=False)
    print("Wrote", emp_path)

    regimes_paper = ["exponential", "power law (Tang)", "lognormal (Tang)"]
    net_order = ["ER", "WS", "BA", "Lat4", "Ring"]
    _t = thresh_tbl[thresh_tbl["recovery_regime"].isin(regimes_paper)].copy()
    if not _t.empty and "beta_pred" in _t.columns:
        _t["ratio_ci_low"] = _t["beta_hat_surv_50_ci_low"] / _t["beta_pred"]
        _t["ratio_ci_high"] = _t["beta_hat_surv_50_ci_high"] / _t["beta_pred"]
        _cols = [
            "recovery_regime",
            "net_label",
            "random_seed",
            "beta_pred",
            "beta_hat_surv_50",
            "beta_hat_surv_50_ci_low",
            "beta_hat_surv_50_ci_high",
            "ratio_surv_50_over_pred",
            "ratio_ci_low",
            "ratio_ci_high",
            "beta_hat_logit_ed50_surv",
        ]
        _cols = [c for c in _cols if c in _t.columns]
        _ratio_by_run = _t[_cols].sort_values(["recovery_regime", "net_label", "random_seed"])
        _ratio_by_run.to_csv(out_dir / "threshold_ratio_by_run.csv", index=False)
        _agg = (
            _t.groupby(["recovery_regime", "net_label"], as_index=False)
            .agg(
                ratio_mean=("ratio_surv_50_over_pred", "mean"),
                ratio_std=("ratio_surv_50_over_pred", "std"),
                n_seeds=("random_seed", "count"),
            )
        )
        _cats = [x for x in net_order if x in set(_agg["net_label"].astype(str))]
        _rest = sorted(set(_agg["net_label"].astype(str)) - set(_cats))
        _agg["net_label"] = pd.Categorical(_agg["net_label"], categories=_cats + _rest, ordered=True)
        _agg = _agg.sort_values(["recovery_regime", "net_label"]).reset_index(drop=True)
        _agg.to_csv(out_dir / "threshold_ratio_summary.csv", index=False)
        print("Wrote", out_dir / "threshold_ratio_by_run.csv")
        print("Wrote", out_dir / "threshold_ratio_summary.csv")

        if plot:
            try:
                import matplotlib.pyplot as plt
                import seaborn as sns

                fig, ax = plt.subplots(figsize=(8.5, 4.8))
                palette = {
                    "exponential": "C0",
                    "power law (Tang)": "C1",
                    "lognormal (Tang)": "C2",
                }
                plot_nets = [n for n in net_order if n in _t["net_label"].unique()] or sorted(
                    _t["net_label"].unique()
                )
                sns.barplot(
                    data=_t,
                    x="net_label",
                    y="ratio_surv_50_over_pred",
                    hue="recovery_regime",
                    order=plot_nets,
                    hue_order=[r for r in regimes_paper if r in _t["recovery_regime"].unique()],
                    palette=[palette[r] for r in regimes_paper if r in _t["recovery_regime"].unique()],
                    errorbar="sd",
                    capsize=0.08,
                    ax=ax,
                )
                ax.axhline(
                    1.0,
                    color="0.25",
                    linestyle="--",
                    lw=1.2,
                    label=r"$\hat{\beta}/\beta_{\mathrm{pred}} = 1$",
                )
                ax.set_xlabel("Network")
                ax.set_ylabel(r"$\hat{\beta}_{\mathrm{surv}\,50} / \beta_{\mathrm{pred}}$")
                ax.set_title(r"Empirical threshold vs spectral $1/(\lambda_{\max}\mathbb{E}[W])$")
                ax.legend(title="Recovery", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
                plt.tight_layout()
                fig_path = out_dir / "threshold_ratio_bar.png"
                fig.savefig(fig_path, dpi=150, bbox_inches="tight")
                plt.close(fig)
                print("Wrote", fig_path)
            except ImportError:
                print("matplotlib/seaborn not available; skip threshold_ratio_bar.png")


def main() -> None:
    parser = argparse.ArgumentParser(description="Aggregate BehaviorSpace outputs to threshold CSVs")
    parser.add_argument("--project-root", type=Path, default=ROOT)
    parser.add_argument("--no-plot", action="store_true")
    parser.add_argument("--n-boot", type=int, default=400)
    args = parser.parse_args()
    root = args.project_root.resolve()
    raw_dir = root / "output" / "raw"
    lambda_path = root / "output" / "lambda_max.csv"
    if not lambda_path.is_file():
        raise SystemExit(f"Missing {lambda_path}. Run pipeline with 01-export-networks first.")
    merged, cols = build_merged(raw_dir, lambda_path)
    write_threshold_outputs(
        merged,
        cols,
        root / "output",
        n_boot=args.n_boot,
        plot=not args.no_plot,
    )


if __name__ == "__main__":
    main()
