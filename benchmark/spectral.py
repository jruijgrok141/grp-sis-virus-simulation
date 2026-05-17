"""Largest adjacency eigenvalue (spectral reference τ_c ≈ 1/λ_max)."""

from __future__ import annotations

import numpy as np
from scipy import sparse as sp
from scipy.sparse.linalg import eigsh

from graph import ERGraph


def largest_eigenvalue_adjacency(graph: ERGraph) -> float:
    n = graph.n
    if n == 0:
        return 0.0
    rows: list[int] = []
    cols: list[int] = []
    for i, js in enumerate(graph.neighbors):
        for j in js:
            rows.append(i)
            cols.append(j)
    if not rows:
        return 0.0
    data = np.ones(len(rows), dtype=np.float64)
    mat = sp.csr_matrix((data, (rows, cols)), shape=(n, n))
    mat.eliminate_zeros()
    vals, _ = eigsh(mat, k=1, which="LA")
    return float(vals[0])


def spectral_reference(recovery_mean: float, lam_max: float) -> dict[str, float]:
    """Mean-field line τ_c = 1/λ_max, β_c = 1/(λ_max E[W])."""
    ew = float(recovery_mean)
    if lam_max <= 0:
        return {"lambda_max": lam_max, "tau_c": float("nan"), "beta_c": float("nan")}
    inv_lam = 1.0 / lam_max
    return {
        "lambda_max": lam_max,
        "tau_c": inv_lam,
        "beta_c": inv_lam / ew if ew > 0 else float("nan"),
    }
