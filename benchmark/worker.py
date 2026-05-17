"""Picklable worker for multiprocessing (avoid __main__ issues on Windows spawn)."""

from __future__ import annotations

from dataclasses import asdict

import numpy as np

from graph import ERGraph
from sis_engine import run_sis_single
from sis_numba import _csr_from_neighbors, numba_available, run_sis_numba


def run_single_task(payload: tuple) -> dict:
    (
        neigh_data,
        n,
        beta,
        recovery_mean,
        regime,
        run_seed,
        initial_infected,
        max_ticks,
        late_window_len,
        power_law_lambda,
        lognormal_sigma,
        use_numba,
    ) = payload

    if use_numba and numba_available():
        indptr, cols = _csr_from_neighbors(neigh_data, n)
        r = run_sis_numba(
            indptr,
            cols,
            n,
            beta,
            recovery_mean,
            regime,
            int(run_seed),
            initial_infected=initial_infected,
            max_ticks=max_ticks,
            late_window_len=late_window_len,
            power_law_lambda=power_law_lambda,
            lognormal_sigma=lognormal_sigma,
        )
    else:
        graph = ERGraph(n=n, neighbors=neigh_data, edges=[])
        rng = np.random.default_rng(run_seed)
        r = run_sis_single(
            graph,
            beta,
            recovery_mean,
            regime,
            rng,
            initial_infected=initial_infected,
            max_ticks=max_ticks,
            late_window_len=late_window_len,
            power_law_lambda=power_law_lambda,
            lognormal_sigma=lognormal_sigma,
        )
    row = asdict(r)
    row["recovery_regime"] = regime
    row["run_seed"] = run_seed
    return row
