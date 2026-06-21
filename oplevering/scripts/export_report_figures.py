#!/usr/bin/env python3
"""
Build PNG figures for report/report.tex from pipeline outputs under output/.

Usage (from project root):
  python scripts/export_report_figures.py

Requires: matplotlib, pandas, numpy; optional networkx for the ER schematic (falls back if missing).

Primary ER seed (figures that stay comparable to the original single-graph report): 10001.
Additional seeds 10002--10003 require a full pipeline run after updating BehaviorSpace experiments.
"""

from __future__ import annotations

import argparse
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Circle, FancyArrowPatch
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.analysis_utils import pick_column, read_behaviorspace_table  # noqa: E402
from scripts.threshold_estimators import survival_curve  # noqa: E402

REPORT_DIR = ROOT / "report"
BS_XML = ROOT / "netlogo" / "behaviorspace_experiments.xml"
RAW = ROOT / "output" / "raw"
EDGES_DIR = ROOT / "output" / "edges"
LAM = ROOT / "output" / "lambda_max.csv"
EMP = ROOT / "output" / "empirical_thresholds.csv"

T05 = RAW / "05-baseline-empirical-threshold-ER_table.csv"
T06 = RAW / "06-baseline-ER-power-law-Tang_table.csv"
T07 = RAW / "07-baseline-ER-lognormal-Tang_table.csv"
T08 = RAW / "08-dynamics-trajectory-ER_table.csv"

# Match BehaviorSpace `expt-seed` values in netlogo experiments (01 / 05--07).
PRIMARY_ER_SEED = 10001
DEFAULT_ER_SEEDS = (10001, 10002, 10003)

REGIME_ORDER = [
    ("exponential", "Exponential", "C0"),
    ("power law (Tang)", "Power law (Tang)", "C1"),
    ("lognormal (Tang)", "Lognormal (Tang)", "C2"),
]

THRESHOLD_EXPERIMENTS = [
    (T05, "exponential", "Exponential", "C0"),
    (T06, "power law (Tang)", "Power law (Tang)", "C1"),
    (T07, "lognormal (Tang)", "Lognormal (Tang)", "C2"),
]

def _cols(df: pd.DataFrame) -> tuple[str, str, str, str]:
    ic = "infection_prob" if "infection_prob" in df.columns else None
    ec = "bs_out_extinct" if "bs_out_extinct" in df.columns else None
    pc = "bs_out_late_mean_prevalence" if "bs_out_late_mean_prevalence" in df.columns else None
    tc = "bs_out_final_tick" if "bs_out_final_tick" in df.columns else None
    assert ic and ec
    return ic, ec, pc or ic, tc or ec


def _filter_seed(df: pd.DataFrame, seed: int) -> pd.DataFrame:
    sc = pick_column(df, "expt_seed", "expt-seed", "random_seed")
    if sc is None:
        return df
    return df.loc[df[sc].astype(int) == int(seed)].copy()


def _discovered_er_seeds() -> list[int]:
    seeds: set[int] = set()
    for p in EDGES_DIR.glob("edges_ER_*.csv"):
        stem = p.stem.replace("edges_ER_", "")
        if stem.isdigit():
            seeds.add(int(stem))
    if seeds:
        return sorted(seeds)
    return list(DEFAULT_ER_SEEDS)


def _edges_er_path(seed: int) -> Path:
    return EDGES_DIR / f"edges_ER_{int(seed)}.csv"


def _beta_grid_from_threshold_tables(primary_seed: int = PRIMARY_ER_SEED) -> np.ndarray:
    """Sorted infection-prob values shared across experiments 05--07."""
    betas: set[float] = set()
    for path, *_ in THRESHOLD_EXPERIMENTS:
        if not path.is_file():
            continue
        df = _filter_seed(read_behaviorspace_table(path), primary_seed)
        ic = pick_column(df, "infection_prob", "infection-prob")
        if ic is None or df.empty:
            continue
        for v in pd.to_numeric(df[ic], errors="coerce").dropna().unique():
            betas.add(float(v))
    if not betas:
        return np.array([0.026, 0.028, 0.036])
    return np.sort(list(betas))


def _snap_to_grid(target: float, grid: np.ndarray) -> float:
    if grid.size == 0:
        return float(target)
    return float(grid[np.argmin(np.abs(grid - float(target)))])


def trajectory_betas_near_pred(
    beta_pred: float,
    grid: np.ndarray | None = None,
) -> tuple[float, float, float]:
    """Three grid values near 0.9, 1.0, 1.3 x beta_pred (for experiment 08 / F9)."""
    g = grid if grid is not None and len(grid) else _beta_grid_from_threshold_tables()
    return (
        _snap_to_grid(0.9 * beta_pred, g),
        _snap_to_grid(beta_pred, g),
        _snap_to_grid(1.3 * beta_pred, g),
    )


def _late_prevalence_summary(
    df: pd.DataFrame,
    ic: str,
    ec: str,
    pc: str,
    *,
    survivors_only: bool,
) -> pd.DataFrame:
    rows: list[dict[str, float | int]] = []
    for beta, g in df.groupby(ic, dropna=False):
        beta_f = float(beta)
        if survivors_only:
            g = g.loc[g[ec].astype(float) < 0.5]
        prev = pd.to_numeric(g[pc], errors="coerce").dropna().to_numpy()
        rows.append(
            {
                ic: beta_f,
                "median": float(np.median(prev)) if len(prev) else float("nan"),
                "q25": float(np.percentile(prev, 25)) if len(prev) else float("nan"),
                "q75": float(np.percentile(prev, 75)) if len(prev) else float("nan"),
                "n": int(len(prev)),
            }
        )
    return pd.DataFrame(rows).sort_values(ic).reset_index(drop=True)


def plot_late_prevalence_vs_beta(
    beta_pred: float,
    out_path: Path,
    *,
    primary_seed: int = PRIMARY_ER_SEED,
    survivors_only: bool = True,
) -> None:
    """
    F6 (RQ2--RQ3): late-window mean prevalence vs beta from experiments 05--07.

    Uses ``bs-out-late-mean-prevalence``; by default only non-extinct runs at the horizon.
    Re-run 05--07 after the NetLogo ``bs-collect-metrics?`` setup fix for non-zero values.
    """
    fig, ax = plt.subplots(figsize=(8.5, 5.0))
    any_data = False
    for path, _regime_key, label, color in THRESHOLD_EXPERIMENTS:
        if not path.is_file():
            continue
        df = read_behaviorspace_table(path)
        df = _filter_seed(df, primary_seed)
        if df.empty:
            continue
        ic, ec, pc, _tc = _cols(df)
        if pc == ic:
            continue
        cur = _late_prevalence_summary(df, ic, ec, pc, survivors_only=survivors_only)
        cur = cur[np.isfinite(cur["median"])]
        if cur.empty:
            continue
        any_data = True
        ax.plot(cur[ic], cur["median"], "o-", ms=4, lw=1.5, color=color, label=label)
        yerr = np.vstack(
            [
                cur["median"].to_numpy() - cur["q25"].to_numpy(),
                cur["q75"].to_numpy() - cur["median"].to_numpy(),
            ]
        )
        ax.errorbar(
            cur[ic],
            cur["median"],
            yerr=yerr,
            fmt="none",
            ecolor=color,
            alpha=0.45,
            capsize=2,
        )
    if not any_data:
        _write_placeholder_figure(
            out_path,
            "No late-prevalence data in output/raw/05--07.\nRun BehaviorSpace experiments 05--07 first.",
        )
        return
    ax.axvline(beta_pred, color="0.35", ls="--", lw=1.2, label=r"$\beta_{\mathrm{pred}}$")
    ylab = (
        r"Late-window mean prevalence $\bar{\rho}_{\mathrm{late}}$ (surviving runs)"
        if survivors_only
        else r"Late-window mean prevalence $\bar{\rho}_{\mathrm{late}}$ (all runs)"
    )
    ax.set_xlabel(r"Infection probability $\beta$")
    ax.set_ylabel(ylab)
    ax.set_title(
        f"ER seed {primary_seed}: late prevalence vs $\\beta$ "
        f"({'non-extinct at horizon' if survivors_only else 'all replicates'})"
    )
    ax.set_ylim(bottom=0)
    ax.legend(loc="best", fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def _trajectory_envelope(
    df: pd.DataFrame,
    *,
    step_col: str,
    prev_col: str,
    run_col: str,
    ic: str,
    beta: float,
    rd_col: str | None,
    regime_key: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray] | None:
    sub = df.loc[(pd.to_numeric(df[ic], errors="coerce") - beta).abs() < 1e-9].copy()
    if rd_col is not None and rd_col in sub.columns:
        sub = sub.loc[sub[rd_col].astype(str).str.lower() == regime_key.lower()]
    if sub.empty:
        return None
    curves: list[np.ndarray] = []
    for _, g in sub.groupby(run_col):
        g2 = g.sort_values(step_col)
        y = pd.to_numeric(g2[prev_col], errors="coerce").dropna().to_numpy(dtype=float)
        if len(y) > 0:
            curves.append(y)
    if not curves:
        return None
    max_len = max(len(c) for c in curves)
    arr = np.full((len(curves), max_len), np.nan, dtype=float)
    for i, c in enumerate(curves):
        arr[i, : len(c)] = c
    steps = np.arange(max_len, dtype=float)
    mean = np.nanmean(arr, axis=0)
    lo = np.nanpercentile(arr, 2.5, axis=0)
    hi = np.nanpercentile(arr, 97.5, axis=0)
    return steps, mean, lo, hi


def _mean_field_reference_curve(
    df: pd.DataFrame,
    *,
    step_col: str,
    mf_col: str,
    run_col: str,
    ic: str,
    beta: float,
    rd_col: str | None,
    regime_key: str,
) -> tuple[np.ndarray, np.ndarray] | None:
    sub = df.loc[(pd.to_numeric(df[ic], errors="coerce") - beta).abs() < 1e-9].copy()
    if rd_col is not None and rd_col in sub.columns:
        sub = sub.loc[sub[rd_col].astype(str).str.lower() == regime_key.lower()]
    if sub.empty:
        return None
    best_run = None
    best_len = 0
    for _, g in sub.groupby(run_col):
        g2 = g.sort_values(step_col)
        if len(g2) > best_len:
            best_len = len(g2)
            best_run = g2
    if best_run is None or best_len == 0:
        return None
    steps = pd.to_numeric(best_run[step_col], errors="coerce").to_numpy(dtype=float)
    mf = pd.to_numeric(best_run[mf_col], errors="coerce").to_numpy(dtype=float)
    return steps, mf


def plot_prevalence_trajectory_bands(
    out_path: Path,
    *,
    primary_seed: int = PRIMARY_ER_SEED,
    beta_pred: float | None = None,
    traj_path: Path = T08,
) -> None:
    """
    F9 (RQ3): mean microscopic prevalence +/- 95% band over replicates (experiment 08).

    One column per recovery regime; curves for beta_low, beta_mid, beta_high near beta_pred.
    """
    if not traj_path.is_file():
        _write_placeholder_figure(
            out_path,
            "Missing output/raw/08-dynamics-trajectory-ER_table.csv\n"
            "Run: python scripts/run_pipeline.py --experiment 08-dynamics-trajectory-ER",
        )
        return
    df = read_behaviorspace_table(traj_path)
    df = _filter_seed(df, primary_seed)
    step_col = pick_column(df, "step")
    gcol = pick_column(df, "bs_out_prev_grp", "bs-out-prev-grp")
    mcol = pick_column(df, "bs_out_prev_mf", "bs-out-prev-mf")
    run_col = pick_column(df, "run_number", "[run number]")
    ic = pick_column(df, "infection_prob", "infection-prob")
    rd_col = pick_column(df, "recovery_distribution", "recovery-distribution")
    if not all([step_col, gcol, mcol, run_col, ic]) or df.empty:
        _write_placeholder_figure(out_path, "08-dynamics CSV missing required columns or is empty.")
        return

    bp = beta_pred if beta_pred is not None else beta_pred_from_lambda_csv(primary_seed)
    grid = _beta_grid_from_threshold_tables(primary_seed)
    b_lo, b_mid, b_hi = trajectory_betas_near_pred(bp, grid)
    beta_specs = [
        (b_lo, rf"$\beta={b_lo:.3f}$ (low)"),
        (b_mid, rf"$\beta={b_mid:.3f}$ (near $\beta_{{\mathrm{{pred}}}}$)"),
        (b_hi, rf"$\beta={b_hi:.3f}$ (high)"),
    ]
    line_styles = ["-", "--", "-."]

    regimes_present: list[tuple[str, str, str]] = []
    if rd_col is not None and rd_col in df.columns:
        known = {r.lower(): (r, L, c) for r, L, c in REGIME_ORDER}
        for val in df[rd_col].astype(str).unique():
            key = val.strip().lower()
            if key in known:
                regimes_present.append(known[key])
    else:
        regimes_present = [(r, L, c) for r, L, c in REGIME_ORDER[:1]]

    if not regimes_present:
        _write_placeholder_figure(out_path, "No recovery regimes found in experiment 08 export.")
        return

    n_reg = len(regimes_present)
    fig, axes = plt.subplots(1, n_reg, figsize=(5.4 * n_reg, 4.8), sharey=True)
    if n_reg == 1:
        axes = [axes]

    any_plotted = False
    for ax, (regime_key, reg_title, base_color) in zip(axes, regimes_present):
        for (beta, beta_lbl), ls in zip(beta_specs, line_styles):
            env = _trajectory_envelope(
                df,
                step_col=step_col,
                prev_col=gcol,
                run_col=run_col,
                ic=ic,
                beta=beta,
                rd_col=rd_col,
                regime_key=regime_key,
            )
            if env is None:
                continue
            steps, mean, lo, hi = env
            any_plotted = True
            ax.plot(steps, mean, ls, color=base_color, lw=1.6, label=beta_lbl)
            ax.fill_between(steps, lo, hi, color=base_color, alpha=0.18, linewidth=0)
            ref = _mean_field_reference_curve(
                df,
                step_col=step_col,
                mf_col=mcol,
                run_col=run_col,
                ic=ic,
                beta=beta,
                rd_col=rd_col,
                regime_key=regime_key,
            )
            if ref is not None and beta == b_mid:
                ax.plot(ref[0], ref[1], ":", color="0.25", lw=1.2, alpha=0.85, label="Mean-field (one run)")
        ax.set_title(reg_title)
        ax.set_xlabel("Time step")
        ax.grid(True, alpha=0.3)

    if not any_plotted:
        _write_placeholder_figure(
            out_path,
            "08-dynamics export has no rows for the expected\n"
            f"beta in {{{b_lo:.3f}, {b_mid:.3f}, {b_hi:.3f}}} and three recovery laws.\n"
            "Re-run experiment 08-dynamics-trajectory-ER (updated BehaviorSpace XML).",
        )
        return

    axes[0].set_ylabel("Prevalence (fraction infected)")
    handles, labels = axes[-1].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper center", ncol=min(4, len(handles)), fontsize=8, bbox_to_anchor=(0.5, 1.04))
    fig.suptitle(
        f"Experiment 08: microscopic grp-SIS prevalence (mean + 95% band, 24 reps), ER seed {primary_seed}",
        y=1.1,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def _recovery_mean_ticks() -> float:
    try:
        mu, _, _ = load_recovery_calibration_from_bs()
        return float(mu)
    except (KeyError, ValueError, OSError):
        return 5.0


def _plot_survival_curves_on_axis(
    ax: plt.Axes,
    *,
    primary_seed: int,
    x_scale: str,
    ew: float,
) -> bool:
    """Draw survival curves on ``ax``; return False if no data were plotted."""
    any_data = False
    for path, label, color in [
        (T05, "Exponential", "C0"),
        (T06, "Power law (Tang)", "C1"),
        (T07, "Lognormal (Tang)", "C2"),
    ]:
        if not path.is_file():
            continue
        df = read_behaviorspace_table(path)
        df = _filter_seed(df, primary_seed)
        if df.empty:
            continue
        ic, ec, pc, tc = _cols(df)
        cur = survival_curve(df, ic, ec, prev_col=pc, tick_col=tc)
        x = pd.to_numeric(cur[ic], errors="coerce").astype(float)
        if x_scale == "tau":
            x = x * ew
        any_data = True
        ax.plot(
            x,
            cur["p_survive"],
            "o-",
            ms=4,
            lw=1.5,
            color=color,
            label=label,
        )
    return any_data


def plot_survival_curves(beta_pred: float, out_path: Path, *, primary_seed: int = PRIMARY_ER_SEED) -> None:
    ew = _recovery_mean_ticks()
    fig, ax = plt.subplots(figsize=(8.5, 5.0))
    if not _plot_survival_curves_on_axis(ax, primary_seed=primary_seed, x_scale="beta", ew=ew):
        _write_placeholder_figure(out_path, "Missing output/raw/05--07 — run BehaviorSpace threshold experiments.")
        return
    ax.axvline(beta_pred, color="0.35", ls="--", lw=1.2, label=r"$\beta_{\mathrm{pred}}$ (spectral, $c{=}1$)")
    ax.axhline(0.5, color="0.65", ls=":", lw=1.0)
    ax.set_xlabel(r"Infection probability $\beta$")
    ax.set_ylabel("Fraction not extinct at time limit")
    ax.set_title(f"ER expt-seed {primary_seed}: survival vs $\\beta$ (24 replicates per $\\beta$)")
    ax.set_ylim(-0.02, 1.05)
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def plot_survival_curves_tau(
    beta_pred: float,
    out_path: Path,
    *,
    primary_seed: int = PRIMARY_ER_SEED,
) -> None:
    """F1b (RQ1): same survival curves with $tau = beta * E[W]$ on the horizontal axis."""
    ew = _recovery_mean_ticks()
    tau_pred = beta_pred * ew
    fig, ax = plt.subplots(figsize=(8.5, 5.0))
    if not _plot_survival_curves_on_axis(ax, primary_seed=primary_seed, x_scale="tau", ew=ew):
        _write_placeholder_figure(out_path, "Missing output/raw/05--07 — run BehaviorSpace threshold experiments.")
        return
    ax.axvline(tau_pred, color="0.35", ls="--", lw=1.2, label=r"$\tau_{\mathrm{pred}}=1/\lambda_{\max}$ ($c{=}1$)")
    ax.axhline(0.5, color="0.65", ls=":", lw=1.0)
    ax.set_xlabel(r"Effective transmission $\tau = \beta\,\mathbb{E}[W]$")
    ax.set_ylabel("Fraction not extinct at time limit")
    ax.set_title(
        f"ER expt-seed {primary_seed}: survival vs $\\tau$ "
        rf"($\mathbb{{E}}[W]={ew:g}$ ticks)"
    )
    ax.set_ylim(-0.02, 1.05)
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def plot_ratio_bar(emp_path: Path, beta_pred: float, out_path: Path, *, primary_seed: int = PRIMARY_ER_SEED) -> None:
    if not emp_path.is_file():
        return
    t = pd.read_csv(emp_path)
    t = t[t["net_label"].astype(str) == "ER"].copy()
    if "random_seed" in t.columns:
        t = t[t["random_seed"].astype(int) == int(primary_seed)]
    order = ["exponential", "power law (Tang)", "lognormal (Tang)"]
    labels = {
        "exponential": "Exponential",
        "power law (Tang)": "Power law (Tang)",
        "lognormal (Tang)": "Lognormal (Tang)",
    }
    rows = []
    for r in order:
        sub = t[t["recovery_regime"].astype(str).str.lower() == r]
        if sub.empty:
            continue
        rows.append(sub.iloc[0])
    if not rows:
        return
    d = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    x = np.arange(len(d))
    ratios = d["ratio_surv_50_over_pred"].astype(float).values
    err_lo = (d["beta_hat_surv_50"].astype(float) - d["beta_hat_surv_50_ci_low"].astype(float)) / beta_pred
    err_hi = (d["beta_hat_surv_50_ci_high"].astype(float) - d["beta_hat_surv_50"].astype(float)) / beta_pred
    err = np.vstack([err_lo.values, err_hi.values])
    colors = ["C0", "C1", "C2"][: len(d)]
    ax.bar(x, ratios, color=colors, yerr=err, capsize=4, ecolor="0.35", width=0.65)
    ax.axhline(1.0, color="0.2", ls="--", lw=1.2)
    ax.set_xticks(x)
    ax.set_xticklabels([labels.get(str(r).lower(), r) for r in d["recovery_regime"]], rotation=12, ha="right")
    ax.set_ylabel(r"$\hat{\beta}_{\mathrm{surv}\,50} / \beta_{\mathrm{pred}}$")
    ax.set_title(f"Empirical vs spectral threshold (bootstrap CI), ER seed {primary_seed}")
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def plot_median_extinct(emp_path: Path, out_path: Path, *, primary_seed: int = PRIMARY_ER_SEED) -> None:
    if not emp_path.is_file():
        return
    t = pd.read_csv(emp_path)
    t = t[t["net_label"].astype(str) == "ER"].copy()
    if "random_seed" in t.columns:
        t = t[t["random_seed"].astype(int) == int(primary_seed)]
    order = ["exponential", "power law (Tang)", "lognormal (Tang)"]
    labels = {
        "exponential": "Exponential",
        "power law (Tang)": "Power law (Tang)",
        "lognormal (Tang)": "Lognormal (Tang)",
    }
    ys = []
    xs = []
    for r in order:
        sub = t[t["recovery_regime"].astype(str).str.lower() == r]
        if sub.empty:
            continue
        xs.append(labels.get(r, r))
        ys.append(float(sub.iloc[0]["median_tick_if_extinct"]))
    if not ys:
        return
    fig, ax = plt.subplots(figsize=(6.5, 4.0))
    ax.bar(range(len(ys)), ys, color=["C0", "C1", "C2"][: len(ys)])
    ax.set_xticks(range(len(xs)))
    ax.set_xticklabels(xs, rotation=12, ha="right")
    ax.set_ylabel("Median final tick (extinct runs only)")
    ax.set_title(f"Pipeline summary: extinction time by recovery law (ER seed {primary_seed})")
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def _write_placeholder_figure(out_path: Path, message: str) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 3.2))
    ax.text(0.5, 0.5, message, ha="center", va="center", fontsize=10, ma="center")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close(fig)


def plot_threshold_ratio_by_seed(emp_path: Path, out_path: Path) -> None:
    """Grouped bars: recovery regime x expt-seed, ratio hat_beta/beta_pred (from empirical_thresholds.csv)."""
    if not emp_path.is_file():
        _write_placeholder_figure(out_path, "Missing output/empirical_thresholds.csv — run aggregate_thresholds after BehaviorSpace.")
        return
    t = pd.read_csv(emp_path)
    t = t[t["net_label"].astype(str) == "ER"].copy()
    if t.empty or "random_seed" not in t.columns:
        _write_placeholder_figure(out_path, "No ER rows in empirical_thresholds.csv.")
        return
    seeds = sorted(t["random_seed"].astype(int).unique().tolist())
    if len(seeds) < 2:
        _write_placeholder_figure(
            out_path,
            "Only one expt-seed in empirical_thresholds.csv.\n"
            "Run experiments 01 and 05–07 with expt-seed 10001–10003, then aggregate_thresholds.",
        )
        return
    known_lower = set(t["recovery_regime"].astype(str).str.lower())
    regimes = [r for r, _, _ in REGIME_ORDER if r in known_lower]
    if not regimes:
        regimes = sorted(known_lower)
    n_r, n_s = len(regimes), len(seeds)
    fig, ax = plt.subplots(figsize=(max(7.0, 1.2 * n_r * n_s), 4.5))
    x = np.arange(n_r, dtype=float)
    width = 0.8 / n_s
    cmap = plt.cm.tab10(np.linspace(0, 0.9, n_s))
    for i, s in enumerate(seeds):
        offs = x - 0.4 + width / 2 + i * width
        ratios = []
        for r in regimes:
            row = t[
                (t["recovery_regime"].astype(str).str.lower() == r)
                & (t["random_seed"].astype(int) == s)
            ]
            ratios.append(float(row.iloc[0]["ratio_surv_50_over_pred"]) if not row.empty else float("nan"))
        ax.bar(offs, ratios, width=width * 0.92, color=cmap[i], label=f"seed {s}")
    ax.axhline(1.0, color="0.2", ls="--", lw=1.2)
    ax.set_xticks(x)
    lbls = {r: L for r, L, _ in REGIME_ORDER}
    ax.set_xticklabels([lbls.get(r, r) for r in regimes], rotation=14, ha="right")
    ax.set_ylabel(r"$\hat{\beta}_{\mathrm{surv}\,50} / \beta_{\mathrm{pred}}$")
    ax.set_title("ER robustness: empirical survival threshold vs spectral prediction by graph draw")
    ax.legend(title="expt-seed", fontsize=9)
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def _median_extinct_by_beta(df: pd.DataFrame, ic: str, ec: str, tc: str) -> pd.DataFrame:
    rows = []
    for beta, g in df.groupby(ic, dropna=False):
        ext = g[ec].astype(float) >= 0.5
        ticks = pd.to_numeric(g.loc[ext, tc], errors="coerce").dropna()
        rows.append(
            {
                ic: beta,
                "median_tick_extinct": float(np.median(ticks)) if len(ticks) else float("nan"),
                "n_extinct": int(ext.sum()),
            }
        )
    return pd.DataFrame(rows).sort_values(ic).reset_index(drop=True)


def plot_extinct_median_vs_beta_by_seed(out_path: Path, seeds: list[int] | None = None) -> None:
    """
    Facet one panel per expt-seed; each panel shows median final tick (extinct runs only) vs beta
    for the three recovery experiments.
    """
    sets = [
        (T05, "Exponential", "C0"),
        (T06, "Power law (Tang)", "C1"),
        (T07, "Lognormal (Tang)", "C2"),
    ]
    if seeds is None:
        seeds = _discovered_er_seeds()
    # Only plot seeds present in all three files
    present: list[int] = []
    for s in seeds:
        ok = True
        for path, _, _ in sets:
            if not path.is_file():
                ok = False
                break
            df = read_behaviorspace_table(path)
            df_s = _filter_seed(df, s)
            if df_s.empty:
                ok = False
                break
        if ok:
            present.append(s)
    if not present:
        return
    fig, axes = plt.subplots(
        1,
        len(present),
        figsize=(5.2 * len(present), 4.6),
        sharey=True,
    )
    if len(present) == 1:
        axes = [axes]
    for ax, seed in zip(axes, present):
        for path, label, color in sets:
            if not path.is_file():
                continue
            df = read_behaviorspace_table(path)
            df = _filter_seed(df, seed)
            ic, ec, _pc, tc = _cols(df)
            cur = _median_extinct_by_beta(df, ic, ec, tc)
            cur = cur[np.isfinite(cur["median_tick_extinct"])]
            if cur.empty:
                continue
            ax.plot(
                cur[ic],
                cur["median_tick_extinct"],
                "o-",
                ms=3,
                lw=1.4,
                color=color,
                label=label,
            )
        ax.set_xlabel(r"Infection probability $\beta$")
        ax.set_title(f"expt-seed {seed}")
        ax.grid(True, alpha=0.3)
    axes[0].set_ylabel("Median final tick (extinct runs)")
    handles, labels = axes[0].get_legend_handles_labels()
    if not handles:
        for ax in axes:
            handles, labels = ax.get_legend_handles_labels()
            if handles:
                break
    fig.legend(handles, labels, loc="upper center", ncol=3, fontsize=8, bbox_to_anchor=(0.5, 1.02))
    fig.suptitle("Extinction-time summary vs $\\beta$ (finite-horizon runs that absorbed)", y=1.08)
    fig.tight_layout()
    fig.savefig(out_path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def plot_extinct_tick_violins_by_beta(
    out_path: Path,
    *,
    primary_seed: int = PRIMARY_ER_SEED,
    betas: tuple[float, ...] = (0.026, 0.032, 0.038, 0.044),
) -> None:
    """Per regime: violin of extinction times at selected beta (extinct runs, primary seed)."""
    sets = [
        (T05, "Exponential"),
        (T06, "Power law (Tang)"),
        (T07, "Lognormal (Tang)"),
    ]
    fig, axes = plt.subplots(len(sets), 1, figsize=(9.0, 2.8 * len(sets)), sharex=True)
    if len(sets) == 1:
        axes = [axes]
    for ax, (path, reg_title) in zip(axes, sets):
        if not path.is_file():
            ax.set_visible(False)
            continue
        df = read_behaviorspace_table(path)
        df = _filter_seed(df, primary_seed)
        ic, ec, _pc, tc = _cols(df)
        parts_data: list[np.ndarray] = []
        pos: list[int] = []
        labels: list[str] = []
        for j, b in enumerate(betas):
            sub = df.loc[(pd.to_numeric(df[ic], errors="coerce") - b).abs() < 1e-9]
            sub = sub.loc[sub[ec].astype(float) >= 0.5, tc]
            sub = pd.to_numeric(sub, errors="coerce").dropna().to_numpy()
            if len(sub) > 0:
                parts_data.append(sub)
                pos.append(j + 1)
                labels.append(f"{b:.3f}")
        if parts_data:
            vp = ax.violinplot(parts_data, positions=pos, showmeans=True, showmedians=False, widths=0.65)
            for b in vp["bodies"]:
                b.set_alpha(0.75)
        ax.set_ylabel("Final tick")
        ax.set_title(f"{reg_title} — ER seed {primary_seed}")
        ax.grid(True, axis="y", alpha=0.3)
    axes[-1].set_xticks(range(1, len(betas) + 1))
    axes[-1].set_xticklabels([f"{b:.3f}" for b in betas])
    axes[-1].set_xlabel(r"Infection probability $\beta$ (selected grid values)")
    fig.suptitle("Distribution of extinction times conditional on $\\beta$ (extinct replicates only)", y=1.01)
    fig.tight_layout()
    fig.savefig(out_path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def plot_trajectory(out_path: Path) -> None:
    if not T08.is_file():
        return
    df = read_behaviorspace_table(T08)
    step_col = "step" if "step" in df.columns else None
    gcol = "bs_out_prev_grp" if "bs_out_prev_grp" in df.columns else None
    mcol = "bs_out_prev_mf" if "bs_out_prev_mf" in df.columns else None
    run_col = "run_number" if "run_number" in df.columns else None
    if not all([step_col, gcol, mcol, run_col]):
        return
    best_run = None
    best_len = 0
    for rn, g in df.groupby(run_col):
        g2 = g.sort_values(step_col)
        if len(g2) > best_len:
            best_len = len(g2)
            best_run = rn
    sub = df[df[run_col] == best_run].sort_values(step_col)
    fig, ax = plt.subplots(figsize=(8.0, 4.5))
    ax.plot(sub[step_col], sub[gcol], label="Microscopic grp-SIS", lw=1.8)
    ax.plot(sub[step_col], sub[mcol], "--", label="Homogeneous mean-field", lw=1.5, alpha=0.9)
    ax.set_xlabel("Time step")
    ax.set_ylabel("Prevalence (fraction infected)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_title(f"Experiment 08 (run {int(best_run)}), $\\beta={float(sub['infection_prob'].iloc[0]):.4f}$")
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def plot_er_network(out_path: Path, seed: int = PRIMARY_ER_SEED, n_sample: int = 420, layout_seed: int = 42) -> None:
    """Induced subgraph on a random subset of nodes for a readable schematic."""
    edge_path = _edges_er_path(seed)
    if not edge_path.is_file():
        edge_path = ROOT / "output" / "edges" / "edges_ER_10001.csv"
    if not edge_path.is_file():
        return
    edges = pd.read_csv(edge_path)
    c0, c1 = edges.columns[0], edges.columns[1]
    a = pd.to_numeric(edges[c0], errors="coerce")
    b = pd.to_numeric(edges[c1], errors="coerce")
    mask = a.notna() & b.notna()
    a = a.loc[mask].to_numpy(dtype=np.int64)
    b = b.loc[mask].to_numpy(dtype=np.int64)
    rng = np.random.default_rng(layout_seed)
    nodes = np.unique(np.concatenate([a, b]))
    pick = set(rng.choice(nodes, size=min(n_sample, len(nodes)), replace=False))
    keep = np.array([(x in pick) and (y in pick) for x, y in zip(a, b)])
    aa, bb = a[keep], b[keep]
    if len(aa) == 0:
        return
    try:
        import networkx as nx
    except ImportError:
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.text(0.5, 0.5, "Install networkx to render network figure:\n  pip install networkx", ha="center", va="center")
        ax.axis("off")
        fig.savefig(out_path, dpi=120)
        plt.close(fig)
        return

    G = nx.Graph()
    G.add_edges_from(zip(aa.tolist(), bb.tolist()))
    pos = nx.spring_layout(G, seed=layout_seed, k=0.22, iterations=60)
    fig, ax = plt.subplots(figsize=(6.8, 6.8))
    nx.draw_networkx_edges(G, pos, ax=ax, alpha=0.35, width=0.6, edge_color="0.45")
    infected = set(rng.choice(list(G.nodes()), size=min(12, len(G)), replace=False))
    ncolor = ["#c0392b" if n in infected else "#95a5a6" for n in G.nodes()]
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=22, node_color=ncolor, linewidths=0)
    ax.axis("off")
    ax.set_title(f"ER subgraph (seed {seed}, {len(G)} nodes, $|E|={G.number_of_edges()}$), schematic layout")
    fig.tight_layout()
    fig.savefig(out_path, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _load_bs_numeric_constants(experiment_name: str) -> dict[str, float]:
    """Scalar numeric constants from a BehaviorSpace experiment (first value only)."""
    tree = ET.parse(BS_XML)
    root = tree.getroot()
    for expt in root.findall("experiment"):
        if expt.get("name") != experiment_name:
            continue
        consts: dict[str, float] = {}
        const_block = expt.find("constants")
        if const_block is None:
            return consts
        for vs in const_block.findall("enumeratedValueSet"):
            if not vs.findall("value"):
                continue
            raw = vs.findall("value")[0].get("value")
            var = vs.get("variable")
            if raw is None or var is None:
                continue
            try:
                consts[var] = float(raw)
            except ValueError:
                continue
        return consts
    raise KeyError(f"Experiment {experiment_name!r} not found in {BS_XML}")


def load_recovery_calibration_from_bs() -> tuple[float, float, float]:
    """(\texttt{recovery-mean}, \texttt{power-law-lambda}, \texttt{lognormal-sigma}) from expts 05--07."""
    e5 = _load_bs_numeric_constants("05-baseline-empirical-threshold-ER")
    e6 = _load_bs_numeric_constants("06-baseline-ER-power-law-Tang")
    e7 = _load_bs_numeric_constants("07-baseline-ER-lognormal-Tang")
    mu = float(e5["recovery-mean"])
    lam = float(e6["power-law-lambda"])
    sig = float(e7["lognormal-sigma"])
    return mu, lam, sig


def _node_color(state: str, *, focal: bool = False) -> str:
    palette = {
        "S": "#4daf4a",
        "I": "#e41a1c",
        "gray": "#bdbdbd",
    }
    base = palette.get(state, "#bdbdbd")
    if focal:
        return {"S": "#2d7a2d", "I": "#c0392b"}.get(state, base)
    return base


def _draw_local_network(
    ax,
    *,
    focal_state: str,
    neighbor_states: list[str],
    show_transmission: bool = False,
    title: str,
) -> None:
    """Small star network around a focal node (Tang et al. Fig.~1 style)."""
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(-1.55, 1.55)
    ax.set_ylim(-1.55, 1.55)
    angles = np.linspace(0, 2 * np.pi, len(neighbor_states), endpoint=False) + np.pi / 2
    r = 1.05
    focal_xy = (0.0, 0.0)
    node_r = 0.22
    for ang, nstate in zip(angles, neighbor_states):
        xy = (r * np.cos(ang), r * np.sin(ang))
        if show_transmission and nstate == "I" and focal_state == "S":
            arr = FancyArrowPatch(
                xy,
                focal_xy,
                arrowstyle="-|>",
                mutation_scale=12,
                lw=1.8,
                color="#c0392b",
                zorder=1,
            )
            ax.add_patch(arr)
        else:
            ax.plot([xy[0], focal_xy[0]], [xy[1], focal_xy[1]], color="#888888", ls="--", lw=1.2, zorder=0)
        ax.add_patch(
            Circle(
                xy,
                node_r,
                facecolor=_node_color(nstate),
                edgecolor="black",
                lw=1.0,
                zorder=2,
            )
        )
        ax.text(xy[0], xy[1], nstate if nstate != "gray" else "", ha="center", va="center", fontsize=11, fontweight="bold", color="white", zorder=3)
    ax.add_patch(
        Circle(
            focal_xy,
            node_r,
            facecolor=_node_color(focal_state, focal=True),
            edgecolor="black",
            lw=1.4,
            zorder=4,
        )
    )
    ax.text(0, 0, focal_state, ha="center", va="center", fontsize=12, fontweight="bold", color="white", zorder=5)
    ax.set_title(title, fontsize=10, pad=6)


def _draw_waiting_time_pdf(
    ax,
    *,
    x: np.ndarray,
    curves: list[tuple[np.ndarray, str, str]],
    event_t: float = 0.0,
    ylabel: str,
    title: str,
    shade_idx: int = 0,
    mean_line: float | None = None,
    xlabel: str = r"time $t$",
    legend_fontsize: float = 7.5,
) -> None:
    ax.set_title(title, fontsize=10, pad=6)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ymax = 0.0
    for y, label, color in curves:
        ax.plot(x, y, color=color, lw=2.0, label=label)
        ymax = max(ymax, float(np.max(y)))
    if mean_line is not None:
        ax.axvline(mean_line, color="0.35", ls="--", lw=1.3, label=rf"$\mathbb{{E}}[W]={mean_line:g}$")
    if event_t > 0:
        si = min(shade_idx, len(curves) - 1)
        y_shade = curves[si][0]
        ax.fill_between(x, 0, y_shade, where=(x <= event_t), color=curves[si][2], alpha=0.18)
        ax.axvline(event_t, color="0.25", ls=":", lw=1.2)
        ax.text(event_t + 0.05, ymax * 0.82, rf"$\tau={event_t:g}$", fontsize=9)
    ax.set_xlim(0, float(x[-1]))
    ax.set_ylim(0, ymax * 1.12 if ymax > 0 else 1.0)
    ax.grid(True, alpha=0.25)
    if len(curves) > 1 or mean_line is not None:
        ax.legend(loc="upper right", fontsize=legend_fontsize, framealpha=0.9)


def plot_grp_sis_process_schematics(out_path: Path) -> None:
    """Tang et al. (2025) Fig.~1-style schematics with this project's infection/recovery laws."""
    mu, lam, sig = load_recovery_calibration_from_bs()
    lam = max(2.01, lam)
    t0 = mu * (lam - 2.0) / (lam - 1.0)
    alpha = lam - 1.0
    mu_logn = float(np.log(mu) - (sig**2) / 2.0)

    beta_demo = 0.35
    tau_inf = 1.8
    tau_rec = 4.2

    x_inf = np.linspace(0, 8.0, 400)
    y_inf = beta_demo * np.exp(-beta_demo * x_inf)

    x_rec = np.linspace(1e-3, 45.0, 3000)
    y_exp = stats.expon.pdf(x_rec, scale=mu)
    y_pl = stats.pareto.pdf(x_rec, alpha, scale=t0)
    y_ln = stats.lognorm.pdf(x_rec, s=sig, scale=np.exp(mu_logn))
    var_exp = float(stats.expon.var(scale=mu))
    var_pl = float(stats.pareto.var(alpha, scale=t0))
    var_ln = float(stats.lognorm.var(s=sig, scale=np.exp(mu_logn)))

    fig = plt.figure(figsize=(11.5, 6.8))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1], wspace=0.28, hspace=0.42)

    ax_a0 = fig.add_subplot(gs[0, 0])
    ax_a1 = fig.add_subplot(gs[0, 1])
    ax_a2 = fig.add_subplot(gs[0, 2])
    ax_b0 = fig.add_subplot(gs[1, 0])
    ax_b1 = fig.add_subplot(gs[1, 1])
    ax_b2 = fig.add_subplot(gs[1, 2])

    _draw_local_network(
        ax_a0,
        focal_state="S",
        neighbor_states=["I", "I", "S", "I"],
        show_transmission=True,
        title=r"$t=t_0$",
    )
    _draw_waiting_time_pdf(
        ax_a1,
        x=x_inf,
        curves=[(y_inf, rf"$\beta e^{{-{beta_demo:g}t}}$ (memoryless)", "C0")],
        event_t=tau_inf,
        ylabel="density",
        title="Infection waiting time",
        shade_idx=0,
    )
    _draw_local_network(
        ax_a2,
        focal_state="I",
        neighbor_states=["I", "I", "S", "I"],
        show_transmission=False,
        title=rf"$t=t_0+\tau$",
    )

    _draw_local_network(
        ax_b0,
        focal_state="I",
        neighbor_states=["gray", "gray", "gray", "gray"],
        show_transmission=False,
        title=r"$t=t_0$",
    )
    _draw_waiting_time_pdf(
        ax_b1,
        x=x_rec,
        curves=[
            (y_exp, rf"Exponential ($\mathrm{{Var}}={var_exp:.2f}$)", "C0"),
            (y_pl, rf"Power law ($\mathrm{{Var}}={var_pl:.2f}$)", "C1"),
            (y_ln, rf"Lognormal ($\mathrm{{Var}}={var_ln:.2f}$)", "C2"),
        ],
        event_t=tau_rec,
        ylabel=r"$f_W(w)$",
        title=rf"Recovery laws ($\lambda={lam:g}$, $\sigma={sig:g}$)",
        shade_idx=2,
        mean_line=mu,
        xlabel=r"recovery time $W$ (ticks)",
        legend_fontsize=6.5,
    )
    _draw_local_network(
        ax_b2,
        focal_state="S",
        neighbor_states=["gray", "gray", "gray", "gray"],
        show_transmission=False,
        title=rf"$t=t_0+\tau$",
    )

    fig.text(0.02, 0.965, "(a) Infection process", fontsize=11, fontweight="bold", va="top")
    fig.text(
        0.02,
        0.475,
        "(b) Recovery process",
        fontsize=11,
        fontweight="bold",
        va="top",
    )
    fig.suptitle(
        "grp-SIS infection and recovery schematics (adapted from Tang et al., 2025, Fig.~1)",
        fontsize=12,
        y=1.02,
    )
    fig.savefig(out_path, dpi=170, bbox_inches="tight")
    plt.close(fig)


def beta_pred_from_lambda_csv(primary_seed: int = PRIMARY_ER_SEED) -> float:
    beta_pred = 0.0278527163754995
    if not LAM.is_file():
        return beta_pred
    lam = pd.read_csv(LAM)
    er = lam[(lam["net_label"].astype(str) == "ER") & (lam["random_seed"].astype(int) == int(primary_seed))]
    if er.empty:
        er = lam[lam["net_label"].astype(str) == "ER"]
    if not er.empty and "lambda_max" in er.columns:
        lamv = float(er.iloc[0]["lambda_max"])
        ew = 5.0
        beta_pred = 1.0 / (lamv * ew)
    return beta_pred


def main() -> None:
    parser = argparse.ArgumentParser(description="Export report figures from output/")
    parser.add_argument(
        "--primary-seed",
        type=int,
        default=PRIMARY_ER_SEED,
        help="Seed for survival/ratio/extinct summary figures (default: 10001)",
    )
    args = parser.parse_args()
    primary = int(args.primary_seed)

    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    beta_pred = beta_pred_from_lambda_csv(primary)

    plot_er_network(REPORT_DIR / "er_network_example.png", seed=primary)
    plot_grp_sis_process_schematics(REPORT_DIR / "fig_grp_sis_process_schematics.png")
    plot_survival_curves(beta_pred, REPORT_DIR / "fig_er_survival_vs_beta.png", primary_seed=primary)
    plot_survival_curves_tau(beta_pred, REPORT_DIR / "fig_er_survival_vs_tau.png", primary_seed=primary)
    plot_late_prevalence_vs_beta(
        beta_pred,
        REPORT_DIR / "fig_er_late_prevalence_vs_beta.png",
        primary_seed=primary,
        survivors_only=True,
    )
    plot_prevalence_trajectory_bands(
        REPORT_DIR / "fig_er_prevalence_trajectory_bands.png",
        primary_seed=primary,
        beta_pred=beta_pred,
    )
    plot_ratio_bar(EMP, beta_pred, REPORT_DIR / "fig_er_threshold_ratio_bar.png", primary_seed=primary)
    plot_median_extinct(EMP, REPORT_DIR / "fig_er_median_extinct_tick.png", primary_seed=primary)
    plot_trajectory(REPORT_DIR / "fig_er_prevalence_trajectory.png")
    plot_threshold_ratio_by_seed(EMP, REPORT_DIR / "fig_er_threshold_ratio_by_seed.png")
    plot_extinct_median_vs_beta_by_seed(REPORT_DIR / "fig_er_extinct_median_vs_beta_by_seed.png")
    plot_extinct_tick_violins_by_beta(
        REPORT_DIR / "fig_er_extinct_tick_violin_by_beta.png",
        primary_seed=primary,
    )
    print("Wrote figures to", REPORT_DIR)


if __name__ == "__main__":
    main()
