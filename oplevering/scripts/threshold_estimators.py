"""Empirical epidemic thresholds from BehaviorSpace sweeps (points 1–3, 5)."""

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
    """Aggregate replicates by infection probability: survival and persistence rates."""
    base_cols = [inf_col, "n", "p_survive"]
    if not len(sub) or inf_col not in sub.columns:
        return pd.DataFrame(columns=base_cols)
    rows = []
    # groupby drops NaN keys; if infection_prob is all-NaN, rows stays empty
    for beta, g in sub.groupby(inf_col, dropna=False):
        ext = g[ext_col].astype(float)
        n = len(g)
        p_surv = float((ext < 0.5).mean())
        out: dict = {inf_col: beta, "n": n, "p_survive": p_surv}
        if prev_col and prev_col in g.columns:
            prev = pd.to_numeric(g[prev_col], errors="coerce").fillna(0.0)
            p_persist = float(((ext < 0.5) & (prev > 1e-4)).mean())
            out["p_persist_strict"] = p_persist
            p_relaxed = float(((ext < 0.5) & (prev > 1e-6)).mean())
            out["p_persist_relaxed"] = p_relaxed
        if tick_col and tick_col in g.columns:
            ft = pd.to_numeric(g[tick_col], errors="coerce")
            ex = ext >= 0.5
            out["median_tick_if_extinct"] = float(np.nanmedian(ft[ex])) if ex.any() else np.nan
            out["median_tick_if_survive"] = float(np.nanmedian(ft[~ex])) if (~ex).any() else np.nan
        rows.append(out)
    if not rows:
        extra: list[str] = []
        if prev_col and prev_col in sub.columns:
            extra.extend(["p_persist_strict", "p_persist_relaxed"])
        if tick_col and tick_col in sub.columns:
            extra.extend(["median_tick_if_extinct", "median_tick_if_survive"])
        return pd.DataFrame(columns=base_cols + extra)
    return pd.DataFrame(rows).sort_values(inf_col).reset_index(drop=True)


def threshold_smallest_beta(
    curve: pd.DataFrame,
    beta_col: str,
    prob_col: str,
    q: float = 0.5,
    fallback_max: bool = True,
) -> float:
    """Smallest $\\beta$ with estimated probability $\\ge q$; else max $\\beta$ on grid (if fallback_max)."""
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


def logistic_ed50_linearized(beta: np.ndarray, p: np.ndarray) -> float:
    """$\\mathrm{logit}(p) \\approx b_0 + b_1 \\beta$; ED50 where $p=0.5$ is $-b_0/b_1$."""
    beta = np.asarray(beta, dtype=float)
    p = np.asarray(p, dtype=float)
    ok = np.isfinite(beta) & np.isfinite(p) & (p > 0) & (p < 1)
    beta, p = beta[ok], p[ok]
    if len(beta) < 3 or np.nanstd(beta) < 1e-15:
        return float("nan")
    p = np.clip(p, 1e-4, 1.0 - 1e-4)
    z = np.log(p / (1.0 - p))
    A = np.column_stack([np.ones(len(beta)), beta])
    coef, _, _, _ = np.linalg.lstsq(A, z, rcond=None)
    b0, b1 = float(coef[0]), float(coef[1])
    if abs(b1) < 1e-12:
        return float("nan")
    return -b0 / b1


def bootstrap_survival_threshold(
    sub: pd.DataFrame,
    inf_col: str,
    ext_col: str,
    q: float = 0.5,
    n_boot: int = 400,
    rng: np.random.Generator | None = None,
) -> tuple[float, float, float]:
    """Stratified bootstrap over $\\beta$ grid; return median and 95% interval of $\\hat\\beta$."""
    rng = rng or np.random.default_rng(42)
    if sub.empty or inf_col not in sub.columns:
        return float("nan"), float("nan"), float("nan")
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


def row_threshold_summary(
    sub: pd.DataFrame,
    inf_col: str,
    ext_col: str,
    prev_col: str | None,
    tick_col: str | None,
    beta_pred: float,
    n_boot: int = 400,
) -> dict:
    """One dict per (regime, net, seed) for CSV export."""
    if sub.empty or inf_col not in sub.columns:
        return {
            "beta_pred": beta_pred,
            "beta_hat_surv_50": float("nan"),
            "beta_hat_surv_25": float("nan"),
            "beta_hat_surv_50_median_boot": float("nan"),
            "beta_hat_surv_50_ci_low": float("nan"),
            "beta_hat_surv_50_ci_high": float("nan"),
            "beta_hat_logit_ed50_surv": float("nan"),
            "ratio_surv_50_over_pred": float("nan"),
        }
    cur = survival_curve(sub, inf_col, ext_col, prev_col=prev_col, tick_col=tick_col)
    beta_hat_surv_50 = threshold_smallest_beta(cur, inf_col, "p_survive", q=0.5, fallback_max=True)
    beta_hat_surv_25 = threshold_smallest_beta(cur, inf_col, "p_survive", q=0.25, fallback_max=True)
    med, lo, hi = bootstrap_survival_threshold(sub, inf_col, ext_col, q=0.5, n_boot=n_boot)

    ed50 = (
        logistic_ed50_linearized(cur[inf_col].to_numpy(), cur["p_survive"].to_numpy())
        if not cur.empty
        else float("nan")
    )

    out: dict = {
        "beta_pred": beta_pred,
        "beta_hat_surv_50": beta_hat_surv_50,
        "beta_hat_surv_25": beta_hat_surv_25,
        "beta_hat_surv_50_median_boot": med,
        "beta_hat_surv_50_ci_low": lo,
        "beta_hat_surv_50_ci_high": hi,
        "beta_hat_logit_ed50_surv": ed50,
        "ratio_surv_50_over_pred": beta_hat_surv_50 / beta_pred if beta_pred else np.nan,
    }
    if prev_col and "p_persist_strict" in cur.columns:
        out["beta_hat_persist_50"] = threshold_smallest_beta(
            cur, inf_col, "p_persist_strict", q=0.5, fallback_max=True
        )
    if prev_col and "p_persist_relaxed" in cur.columns:
        out["beta_hat_persist_relaxed_50"] = threshold_smallest_beta(
            cur, inf_col, "p_persist_relaxed", q=0.5, fallback_max=True
        )
    if tick_col and "median_tick_if_extinct" in cur.columns:
        exm = cur["median_tick_if_extinct"]
        svm = cur["median_tick_if_survive"]
        out["median_tick_if_extinct"] = (
            float(np.nanmedian(exm.to_numpy())) if exm.notna().any() else float("nan")
        )
        out["median_tick_if_survive"] = (
            float(np.nanmedian(svm.to_numpy())) if svm.notna().any() else float("nan")
        )
    return out
