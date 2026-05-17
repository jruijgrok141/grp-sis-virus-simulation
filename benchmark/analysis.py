"""Aggregate survival curves and empirical thresholds (aligned with scripts/threshold_estimators)."""

from __future__ import annotations

import numpy as np
import pandas as pd


def survival_curve(
    sub: pd.DataFrame,
    inf_col: str,
    ext_col: str,
    prev_col: str | None = None,
    tick_col: str | None = None,
) -> pd.DataFrame:
    rows = []
    for beta, g in sub.groupby(inf_col, dropna=False):
        ext = g[ext_col].astype(float)
        n = len(g)
        p_surv = float((ext < 0.5).mean())
        out: dict = {inf_col: beta, "n": n, "p_survive": p_surv}
        if prev_col and prev_col in g.columns:
            prev = pd.to_numeric(g[prev_col], errors="coerce").fillna(0.0)
            p_persist = float(((ext < 0.5) & (prev > 1e-4)).mean())
            out["p_persist_strict"] = p_persist
        if tick_col and tick_col in g.columns:
            ft = pd.to_numeric(g[tick_col], errors="coerce")
            ex = ext >= 0.5
            out["median_tick_if_extinct"] = float(np.nanmedian(ft[ex])) if ex.any() else np.nan
            out["median_tick_if_survive"] = float(np.nanmedian(ft[~ex])) if (~ex).any() else np.nan
        rows.append(out)
    if not rows:
        return pd.DataFrame(columns=[inf_col, "n", "p_survive"])
    return pd.DataFrame(rows).sort_values(inf_col).reset_index(drop=True)


def threshold_smallest_beta(
    curve: pd.DataFrame,
    beta_col: str,
    prob_col: str,
    q: float = 0.5,
    fallback_max: bool = True,
) -> float:
    if curve.empty or beta_col not in curve.columns or prob_col not in curve.columns:
        return float("nan")
    g = curve.sort_values(beta_col)
    above = g[g[prob_col] >= q]
    if above.empty:
        if not fallback_max:
            return float("nan")
        mx = g[beta_col].max()
        return float(mx) if pd.notna(mx) else float("nan")
    return float(above.iloc[0][beta_col])


def bootstrap_survival_threshold(
    sub: pd.DataFrame,
    inf_col: str,
    ext_col: str,
    q: float = 0.5,
    n_boot: int = 400,
    rng: np.random.Generator | None = None,
) -> tuple[float, float, float]:
    rng = rng or np.random.default_rng(42)
    betas = sorted(sub[inf_col].dropna().unique().tolist())
    if not betas:
        return float("nan"), float("nan"), float("nan")
    est: list[float] = []
    for _ in range(n_boot):
        rows = []
        for b in betas:
            g = sub[sub[inf_col] == b]
            if len(g) == 0:
                continue
            idx = rng.integers(0, len(g), size=len(g))
            samp = g.iloc[idx]
            p_surv = float((samp[ext_col].astype(float) < 0.5).mean())
            rows.append((b, p_surv))
        if not rows:
            continue
        cur = pd.DataFrame(rows, columns=[inf_col, "p_surv"]).sort_values(inf_col)
        est.append(threshold_smallest_beta(cur, inf_col, "p_surv", q=q, fallback_max=True))
    if not est:
        return float("nan"), float("nan"), float("nan")
    arr = np.array(est, dtype=float)
    return float(np.median(arr)), float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5))


def summarize_regime(
    df: pd.DataFrame,
    beta_pred: float,
    n_boot: int = 400,
) -> dict:
    inf_col = "infection_prob"
    ext_col = "bs_out_extinct"
    prev_col = "bs_out_late_mean_prevalence"
    tick_col = "bs_out_final_tick"
    cur = survival_curve(df, inf_col, ext_col, prev_col=prev_col, tick_col=tick_col)
    beta_hat = threshold_smallest_beta(cur, inf_col, "p_survive", q=0.5, fallback_max=True)
    med, lo, hi = bootstrap_survival_threshold(df, inf_col, ext_col, q=0.5, n_boot=n_boot)

    out = {
        "beta_pred": beta_pred,
        "beta_hat_surv_50": beta_hat,
        "beta_hat_surv_50_median_boot": med,
        "beta_hat_surv_50_ci_low": lo,
        "beta_hat_surv_50_ci_high": hi,
        "ratio_surv_50_over_pred": beta_hat / beta_pred if beta_pred and np.isfinite(beta_pred) else np.nan,
    }
    if "p_persist_strict" in cur.columns:
        out["beta_hat_persist_50"] = threshold_smallest_beta(
            cur, inf_col, "p_persist_strict", q=0.5, fallback_max=True
        )
    return out
