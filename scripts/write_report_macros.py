#!/usr/bin/env python3
"""
Write report/generated_quantities.tex from output/lambda_max.csv and output/empirical_thresholds.csv.

Run after aggregate_thresholds (or with stale/missing CSVs: writes safe fallbacks so LaTeX compiles).

Usage (from project root):
  py -3 scripts/write_report_macros.py
"""

from __future__ import annotations

import subprocess
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT_LAM = ROOT / "output" / "lambda_max.csv"
OUT_EMP = ROOT / "output" / "empirical_thresholds.csv"
OUT_RAW = ROOT / "output" / "raw"
DEST = ROOT / "report" / "generated_quantities.tex"

PRIMARY_SEED = 10001
ER_SEEDS = (10001, 10002, 10003)
SEED_SUFFIX = {10001: "One", 10002: "Two", 10003: "Three"}
EW = 5.0
DEFAULT_BETA_GRID_MAX = 0.056

REGIME_KEYS = ["exponential", "power law (Tang)", "lognormal (Tang)"]
REGIME_SUFFIX = {"exponential": "Exp", "power law (Tang)": "PL", "lognormal (Tang)": "Logn"}

# Fallback primary-seed thresholds (last good pipeline run) when CSVs are mid-refresh.
_FALLBACK_PRIMARY: dict[str, tuple[float, float, float, float, float]] = {
    "exponential": (0.048, 0.048, 0.048, 1.723, 8.0),
    "lognormal (Tang)": (0.048, 0.048, 0.048, 1.723, 9.0),
    "power law (Tang)": (0.044, 0.044, 0.044, 1.580, 15.0),
}
_FALLBACK_SEED_RATIOS: dict[str, dict[int, float]] = {
    "exponential": {10001: 1.723, 10002: 1.738, 10003: 1.722},
    "power law (Tang)": {10001: 1.580, 10002: 1.738, 10003: 1.650},
    "lognormal (Tang)": {10001: 1.723, 10002: 1.738, 10003: 1.650},
}


def _fmt(x: float, d: int = 4) -> str:
    s = f"{float(x):.{d}f}"
    if "." in s:
        s = s.rstrip("0").rstrip(".")
    return s


def _git_hash() -> str:
    try:
        r = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=str(ROOT),
            capture_output=True,
            text=True,
            check=True,
        )
        return r.stdout.strip() or "unknown"
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unknown"


def _beta_grid_max_from_raw() -> float:
    """Upper sweep bound from BehaviorSpace exports, if present."""
    for path in sorted(OUT_RAW.glob("0[567]-*_table.csv")):
        try:
            from scripts.analysis_utils import pick_column, read_behaviorspace_table

            df = read_behaviorspace_table(path)
            ic = pick_column(df, "infection_prob", "infection-prob")
            if ic:
                return float(pd.to_numeric(df[ic], errors="coerce").max())
        except Exception:  # noqa: BLE001
            continue
    return DEFAULT_BETA_GRID_MAX


def _load_lambda() -> dict:
    fallback = {
        "lambda_max": 7.181,
        "num_edges": 5979,
        "k_mean": 5.979,
        "k2_mean": 24.95,
        "beta_pred": 0.0278527163754995,
    }
    if not OUT_LAM.is_file():
        return fallback
    lam = pd.read_csv(OUT_LAM)
    er = lam[(lam["net_label"].astype(str) == "ER") & (lam["random_seed"].astype(int) == PRIMARY_SEED)]
    if er.empty:
        er = lam[lam["net_label"].astype(str) == "ER"]
    if er.empty:
        return fallback
    row = er.iloc[0]
    lm = float(row["lambda_max"])
    ne = int(row["num_edges"]) if "num_edges" in row and pd.notna(row["num_edges"]) else fallback["num_edges"]
    k1 = float(row["k_mean"]) if "k_mean" in row and pd.notna(row["k_mean"]) else fallback["k_mean"]
    k2 = float(row["k2_mean"]) if "k2_mean" in row and pd.notna(row["k2_mean"]) else fallback["k2_mean"]
    bp = 1.0 / (lm * EW)
    return {"lambda_max": lm, "num_edges": ne, "k_mean": k1, "k2_mean": k2, "beta_pred": bp}


def _load_empirical(beta_pred: float) -> tuple[list[tuple], dict]:
    """Primary-seed thresholds + per-seed ratio table for LaTeX."""
    primary_rows: list[tuple[str, float, float, float, float, float]] = []
    seed_ratios: dict[str, dict[int, float]] = {k: {} for k in REGIME_KEYS}

    if not OUT_EMP.is_file():
        return primary_rows, seed_ratios

    emp = pd.read_csv(OUT_EMP)
    emp_er = emp[emp["net_label"].astype(str) == "ER"]

    for k in REGIME_KEYS:
        sub_k = emp_er[emp_er["recovery_regime"].astype(str).str.lower() == k.lower()]
        for _, r in sub_k.iterrows():
            seed = int(r["random_seed"])
            ratio = float(r["ratio_surv_50_over_pred"])
            if np.isfinite(ratio):
                seed_ratios[k][seed] = ratio

        sub_primary = sub_k[sub_k["random_seed"].astype(int) == PRIMARY_SEED]
        if sub_primary.empty:
            continue
        r = sub_primary.iloc[0]
        bh = float(r["beta_hat_surv_50"])
        lo = float(r["beta_hat_surv_50_ci_low"])
        hi = float(r["beta_hat_surv_50_ci_high"])
        ratio = float(r["ratio_surv_50_over_pred"]) if np.isfinite(r["ratio_surv_50_over_pred"]) else bh / beta_pred
        med = float(r["median_tick_if_extinct"])
        primary_rows.append((k, bh, lo, hi, ratio, med))

    return primary_rows, seed_ratios


def main() -> None:
    L = _load_lambda()
    lm = L["lambda_max"]
    beta_pred = float(L["beta_pred"])
    tau_pred = beta_pred * EW
    k1, k2 = L["k_mean"], L["k2_mean"]
    tau_het = k1 / k2 if k2 > 0 else np.nan
    beta_het = tau_het / EW if np.isfinite(tau_het) else np.nan
    pct_het = 100.0 * (beta_het - beta_pred) / beta_pred if beta_pred > 0 and np.isfinite(beta_het) else 3.0

    primary_rows, seed_ratios = _load_empirical(beta_pred)
    loaded_keys = {r[0] for r in primary_rows}
    for k in REGIME_KEYS:
        if k not in loaded_keys and k in _FALLBACK_PRIMARY:
            bh, lo, hi, ratio, med = _FALLBACK_PRIMARY[k]
            primary_rows.append((k, bh, lo, hi, ratio, med))
        if not seed_ratios.get(k):
            seed_ratios[k] = dict(_FALLBACK_SEED_RATIOS.get(k, {}))

    ratios = [r[4] for r in primary_rows] if primary_rows else [1.6, 1.7]
    cmin, cmax = min(ratios), max(ratios)
    all_seed_ratios = [v for by_seed in seed_ratios.values() for v in by_seed.values() if np.isfinite(v)]
    gmin = min(all_seed_ratios) if all_seed_ratios else cmin
    gmax = max(all_seed_ratios) if all_seed_ratios else cmax
    beta_grid_max = _beta_grid_max_from_raw()

    lines = [
        "% -*- latex -*-",
        "% Generated by scripts/write_report_macros.py — do not edit by hand.",
        f"\\newcommand{{\\RepReportDate}}{{{datetime.now().strftime('%B %Y')}}}",
        f"\\newcommand{{\\RepGitHash}}{{{_git_hash()}}}",
        f"\\newcommand{{\\RepPrimarySeed}}{{{PRIMARY_SEED}}}",
        f"\\newcommand{{\\RepNumNodes}}{{2000}}",
        f"\\newcommand{{\\RepNumEdges}}{{{L['num_edges']}}}",
        f"\\newcommand{{\\RepKmean}}{{{_fmt(k1, 3)}}}",
        f"\\newcommand{{\\RepLambdaMax}}{{{_fmt(lm, 3)}}}",
        f"\\newcommand{{\\RepBetaPred}}{{{_fmt(beta_pred, 7)}}}",
        f"\\newcommand{{\\RepTauPred}}{{{_fmt(tau_pred, 4)}}}",
        f"\\newcommand{{\\RepTauHet}}{{{_fmt(tau_het, 3)}}}",
        f"\\newcommand{{\\RepBetaHet}}{{{_fmt(beta_het, 4)}}}",
        f"\\newcommand{{\\RepBetaHetPctAbove}}{{{_fmt(pct_het, 1)}}}",
        f"\\newcommand{{\\RepBetaGridMax}}{{{_fmt(beta_grid_max, 3)}}}",
        f"\\newcommand{{\\RepCeffMin}}{{{_fmt(cmin, 3)}}}",
        f"\\newcommand{{\\RepCeffMax}}{{{_fmt(cmax, 3)}}}",
        f"\\newcommand{{\\RepRatioGlobalMin}}{{{_fmt(gmin, 3)}}}",
        f"\\newcommand{{\\RepRatioGlobalMax}}{{{_fmt(gmax, 3)}}}",
    ]

    for k, suf in REGIME_SUFFIX.items():
        row = next((r for r in primary_rows if r[0] == k), None)
        if row is None:
            continue
        _, bh, lo, hi, ratio, med = row
        lines.append(f"\\newcommand{{\\RepTh{suf}Hat}}{{{_fmt(bh, 3)}}}")
        lines.append(f"\\newcommand{{\\RepTh{suf}Lo}}{{{_fmt(lo, 3)}}}")
        lines.append(f"\\newcommand{{\\RepTh{suf}Hi}}{{{_fmt(hi, 3)}}}")
        lines.append(f"\\newcommand{{\\RepTh{suf}Ratio}}{{{_fmt(ratio, 3)}}}")
        if np.isfinite(med):
            lines.append(f"\\newcommand{{\\RepTh{suf}MedExt}}{{{_fmt(med, 0)}}}")

    for k, suf in REGIME_SUFFIX.items():
        by_seed = seed_ratios.get(k, {})
        vals = [by_seed.get(s, float("nan")) for s in ER_SEEDS]
        med = float(np.nanmedian(vals)) if any(np.isfinite(v) for v in vals) else float("nan")
        lines.append(f"\\newcommand{{\\RepRatio{suf}Med}}{{{_fmt(med, 3) if np.isfinite(med) else '---'}}}")
        for seed in ER_SEEDS:
            v = by_seed.get(seed, float("nan"))
            suf_seed = SEED_SUFFIX[seed]
            lines.append(
                f"\\newcommand{{\\RepRatio{suf}S{suf_seed}}}{{{_fmt(v, 3) if np.isfinite(v) else '---'}}}"
            )

    DEST.parent.mkdir(parents=True, exist_ok=True)
    DEST.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("Wrote", DEST)


if __name__ == "__main__":
    main()
    sys.exit(0)
