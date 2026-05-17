"""Optional Numba-accelerated grp-SIS (same discrete-time logic as sis_engine)."""

from __future__ import annotations

import numpy as np

try:
    from numba import njit
except ImportError:
    njit = None  # type: ignore

from sis_engine import SimulationResult


def _csr_from_neighbors(neighbors: list[list[int]], n: int) -> tuple[np.ndarray, np.ndarray]:
    indptr = np.zeros(n + 1, dtype=np.int32)
    for i in range(n):
        indptr[i + 1] = indptr[i] + len(neighbors[i])
    cols = np.empty(indptr[n], dtype=np.int32)
    k = 0
    for i in range(n):
        for j in neighbors[i]:
            cols[k] = j
            k += 1
    return indptr, cols


if njit is not None:

    @njit(cache=True)
    def _draw_recovery_inline(
        regime: int,
        recovery_mean: float,
        power_law_lambda: float,
        lognormal_sigma: float,
    ) -> int:
        rm = recovery_mean
        if regime == 0:
            mu = 1.0 / rm
            u = np.random.random()
            w = -np.log(max(1e-6, 1.0 - u)) / mu
            return int(max(1.0, round(w)))
        if regime == 2:
            sig = max(1e-6, lognormal_sigma)
            mu_ln = np.log(rm) - (sig * sig) / 2.0
            w = np.exp(mu_ln + sig * np.random.standard_normal())
            return int(max(1.0, round(w)))
        lam = max(2.01, power_law_lambda)
        t0 = rm * (lam - 2.0) / (lam - 1.0)
        v = max(1e-12, np.random.random())
        t = t0 * (v ** (1.0 / (1.0 - lam)))
        return int(max(1.0, round(t)))

    @njit(cache=True)
    def _run_sis_numba_core(
        indptr: np.ndarray,
        cols: np.ndarray,
        n: int,
        beta: float,
        recovery_mean: float,
        regime: int,
        initial_infected: int,
        max_ticks: int,
        late_window_len: int,
        power_law_lambda: float,
        lognormal_sigma: float,
        run_seed: int,
    ) -> tuple[int, int, float]:
        np.random.seed(run_seed)

        state = np.zeros(n, dtype=np.int8)
        inf_time = np.zeros(n, dtype=np.int32)
        rec_deadline = np.zeros(n, dtype=np.int32)

        perm = np.random.permutation(n)
        k0 = initial_infected if initial_infected < n else n
        for t in range(k0):
            i = perm[t]
            state[i] = 1
            inf_time[i] = 0
            rec_deadline[i] = _draw_recovery_inline(regime, recovery_mean, power_law_lambda, lognormal_sigma)

        prev_buf = np.zeros(max_ticks + 1, dtype=np.float64)
        n_prev = 0
        ticks = 0

        while ticks < max_ticks:
            has_inf = False
            for i in range(n):
                if state[i] == 1:
                    has_inf = True
                    break
            if not has_inf:
                break

            for i in range(n):
                if state[i] != 1:
                    continue
                for e in range(indptr[i], indptr[i + 1]):
                    j = cols[e]
                    if state[j] != 0:
                        continue
                    if np.random.random() < beta:
                        state[j] = 1
                        inf_time[j] = 0
                        rec_deadline[j] = _draw_recovery_inline(
                            regime, recovery_mean, power_law_lambda, lognormal_sigma
                        )

            for i in range(n):
                if state[i] != 1:
                    continue
                inf_time[i] += 1
                if inf_time[i] >= rec_deadline[i]:
                    state[i] = 0
                    inf_time[i] = 0
                    rec_deadline[i] = 0

            n_i = 0
            for i in range(n):
                if state[i] == 1:
                    n_i += 1
            prev_buf[n_prev] = n_i / float(n)
            n_prev += 1
            ticks += 1

            has_inf2 = False
            for i in range(n):
                if state[i] == 1:
                    has_inf2 = True
                    break
            if not has_inf2:
                break

        extinct = 1
        for i in range(n):
            if state[i] == 1:
                extinct = 0
                break

        late_prev = 0.0
        if n_prev > 0:
            w = late_window_len
            start = n_prev - w
            if start < 0:
                start = 0
            s = 0.0
            for k in range(start, n_prev):
                s += prev_buf[k]
            late_prev = s / float(n_prev - start)

        return extinct, ticks, late_prev


def regime_to_code(regime: str) -> int:
    if regime == "exponential":
        return 0
    if regime == "power_law_tang":
        return 1
    if regime == "lognormal_tang":
        return 2
    raise ValueError(regime)


def run_sis_numba(
    indptr: np.ndarray,
    cols: np.ndarray,
    n: int,
    infection_prob: float,
    recovery_mean: float,
    regime: str,
    run_seed: int,
    *,
    initial_infected: int = 5,
    max_ticks: int = 10_000,
    late_window_len: int = 200,
    power_law_lambda: float = 4.24,
    lognormal_sigma: float = 1.0,
) -> SimulationResult:
    if njit is None:
        raise RuntimeError("numba not installed")
    code = regime_to_code(regime)
    ex, ticks, late_prev = _run_sis_numba_core(
        indptr.astype(np.int32),
        cols.astype(np.int32),
        n,
        float(infection_prob),
        float(recovery_mean),
        code,
        int(initial_infected),
        int(max_ticks),
        int(late_window_len),
        float(power_law_lambda),
        float(lognormal_sigma),
        int(run_seed),
    )
    beta = float(infection_prob)
    rm = float(recovery_mean)
    return SimulationResult(
        bs_out_extinct=float(ex),
        bs_out_final_tick=float(ticks),
        bs_out_late_mean_prevalence=float(late_prev),
        infection_prob=beta,
        recovery_mean=rm,
        tau_sim=beta * rm,
    )


def numba_available() -> bool:
    return njit is not None
