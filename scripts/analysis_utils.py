"""Shared helpers for loading BehaviorSpace CSVs and spectral summaries."""

from __future__ import annotations

import io
import re
from pathlib import Path

import numpy as np
import pandas as pd


def net_label_from_behaviorspace(interconnection: str, network_type: object | None) -> str:
    """
    Resolve net label when BehaviorSpace uses interconnection-structure = \"other\":
    then the actual topology comes from network-type (e.g. Erdős–Rényi → ER).
    """
    ics = str(interconnection).strip().strip('"')
    if ics.lower() == "other" and network_type is not None:
        nt = str(network_type).strip().strip('"')
        if nt and nt.lower() != "nan":
            return net_label_from_network_type(nt)
    return net_label_from_network_type(ics)


def net_label_from_network_type(network_type: str) -> str:
    """
    Map a topology description (from either `network-type` or `interconnection-structure`)
    to a short label used in filenames and tables.
    """
    s = str(network_type).strip()
    # Regular structures first (used when interconnection-structure = Lattice4 / Ring)
    if "Lattice4" in s:
        return "Lat4"
    if "Lattice8" in s:
        return "Lat8"
    if "Ring" == s or " Ring" in s:
        return "Ring"

    # Legacy complex-network types from the original design
    if "Erdos" in s or "random" in s.lower():
        return "ER"
    if "Watts" in s or "small-world" in s.lower():
        return "WS"
    if "Barabasi" in s or "scale-free" in s.lower():
        return "BA"
    return "NA"


def normalize_column_key(name: object) -> str:
    """
    Normalized column name — must match what normalize_bs_columns applies to headers.
    Use this when resolving logical names (e.g. infection_prob vs infection-prob vs 'infection prob').
    """
    return (
        re.sub(r"\s+", "_", str(name).strip().strip('"').lower())
        .replace("[", "")
        .replace("]", "")
        .replace("-", "_")
    )


def normalize_bs_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out.columns = [normalize_column_key(c) for c in out.columns]
    return out


def pick_column(df: pd.DataFrame, *candidates: str) -> str | None:
    """Map logical names to an actual column present after normalize_bs_columns."""
    index = {normalize_column_key(c): c for c in df.columns}
    for cand in candidates:
        k = normalize_column_key(cand)
        if k in index:
            return index[k]
    return None


def read_behaviorspace_table(path: Path) -> pd.DataFrame:
    """Load a BehaviorSpace *table* export (CSV). Skips preamble lines."""
    raw = path.read_text(encoding="utf-8", errors="replace")
    lines = raw.splitlines()
    hdr_idx = None
    for i, line in enumerate(lines):
        low = line.lower()
        if "run_number" in low.replace(" ", "_") or "[run number]" in low:
            hdr_idx = i
            break
    if hdr_idx is None:
        for i, line in enumerate(lines):
            if "run number" in line.lower():
                hdr_idx = i
                break
    if hdr_idx is None:
        df = pd.read_csv(io.StringIO(raw))
        return normalize_bs_columns(df)
    df = pd.read_csv(io.StringIO("\n".join(lines[hdr_idx:])))
    return normalize_bs_columns(df)


def _edge_endpoint_arrays(edges: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """Parse first two columns as integer endpoints; drop header rows and junk lines."""
    if "endpoint1" in edges.columns and "endpoint2" in edges.columns:
        c0, c1 = "endpoint1", "endpoint2"
    else:
        c0, c1 = edges.columns[0], edges.columns[1]
    s0 = pd.to_numeric(edges[c0], errors="coerce")
    s1 = pd.to_numeric(edges[c1], errors="coerce")
    mask = s0.notna() & s1.notna()
    a = s0.loc[mask].to_numpy(dtype=np.int64)
    b = s1.loc[mask].to_numpy(dtype=np.int64)
    if len(a) == 0:
        raise ValueError(f"No numeric edge rows in {c0},{c1} columns")
    # One row per undirected pair (handles merged exports / duplicate blocks)
    pairs = np.column_stack([a, b])
    lo = pairs.min(axis=1)
    hi = pairs.max(axis=1)
    uniq = np.unique(np.column_stack([lo, hi]), axis=0)
    return uniq[:, 0], uniq[:, 1]


def adjacency_from_edge_list(
    path: Path, n_nodes: int | None = None
) -> np.ndarray:
    """Build symmetric 0/1 adjacency from CSV with columns endpoint1,endpoint2."""
    edges = pd.read_csv(path)
    a, b = _edge_endpoint_arrays(edges)
    if n_nodes is None:
        n_nodes = int(max(a.max(), b.max()) + 1)
    adj = np.zeros((n_nodes, n_nodes), dtype=np.float64)
    adj[a, b] = 1.0
    adj[b, a] = 1.0
    return adj


def largest_adjacency_eigenvalue(adj: np.ndarray) -> float:
    """Largest eigenvalue of symmetric adjacency (dense)."""
    w = np.linalg.eigvalsh(adj)
    return float(w.max())


def largest_adjacency_eigenvalue_sparse(n: int, a: np.ndarray, b: np.ndarray) -> float:
    """Largest eigenvalue via Lanczos (undirected simple graph)."""
    from scipy import sparse as sp
    from scipy.sparse.linalg import eigsh

    rows = np.concatenate([a, b])
    cols = np.concatenate([b, a])
    data = np.ones(len(rows), dtype=np.float64)
    mat = sp.csr_matrix((data, (rows, cols)), shape=(n, n))
    mat.eliminate_zeros()
    vals, _ = eigsh(mat, k=1, which="LA")
    return float(vals[0])


def degree_moments_from_endpoints(
    a: np.ndarray, b: np.ndarray, n_nodes: int
) -> tuple[float, float, float, float]:
    """
    Mean degree ⟨k⟩, second moment ⟨k²⟩, and MF-style critical-τ surrogates (no spectrum).

    - Homogeneous / well-mixed analogue (exact for k-regular graphs): τ ≈ 1/⟨k⟩.
    - Uncorrelated heterogeneous MF (Pastor–Satorras–style): τ ≈ ⟨k⟩/⟨k²⟩.
    """
    deg = np.zeros(n_nodes, dtype=np.float64)
    np.add.at(deg, a, 1.0)
    np.add.at(deg, b, 1.0)
    k1 = float(deg.mean())
    k2 = float((deg * deg).mean())
    tau_hom = 1.0 / k1 if k1 > 0 else float("nan")
    tau_het = k1 / k2 if k2 > 0 else float("nan")
    return k1, k2, tau_hom, tau_het


def compute_lambda_max_for_edge_file(path: Path, sparse_threshold: int = 800) -> dict:
    edges = pd.read_csv(path)
    a, b = _edge_endpoint_arrays(edges)
    n = int(max(a.max(), b.max()) + 1)
    k_mean, k2_mean, tau_mf_hom, tau_mf_het = degree_moments_from_endpoints(a, b, n)
    if n >= sparse_threshold:
        lam = largest_adjacency_eigenvalue_sparse(n, a, b)
    else:
        adj = adjacency_from_edge_list(path, n_nodes=n)
        lam = largest_adjacency_eigenvalue(adj)
    stem = path.stem
    parts = stem.replace("edges_", "").rsplit("_", 1)
    label = parts[0] if len(parts) == 2 else stem
    seed = int(parts[1]) if len(parts) == 2 and parts[1].isdigit() else None
    return {
        "edge_file": str(path),
        "net_label": label,
        "random_seed": seed,
        "num_nodes": n,
        "num_edges": len(a),
        "k_mean": k_mean,
        "k2_mean": k2_mean,
        "tau_pred_homogeneous_mf": tau_mf_hom,
        "tau_pred_heterogeneous_mf": tau_mf_het,
        "lambda_max": lam,
    }


def persistence_proxy(row: pd.Series, extinct_col: str, prev_col: str, eps: float = 1e-4) -> float:
    """1 if survived with noticeable late prevalence, else 0 (single run)."""
    ex = float(row.get(extinct_col, 1))
    prev = float(row.get(prev_col, 0.0))
    if np.isnan(prev):
        prev = 0.0
    if ex >= 0.5:
        return 0.0
    return 1.0 if prev > eps else 0.0
