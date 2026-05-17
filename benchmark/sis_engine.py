"""Discrete-time grp-SIS on a fixed graph (NetLogo-compatible ordering)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.random import Generator

from graph import ERGraph
from recovery import RecoveryRegime, draw_recovery_time


@dataclass
class SimulationResult:
    bs_out_extinct: float  # 1.0 if no infected at end, else 0.0
    bs_out_final_tick: float
    bs_out_late_mean_prevalence: float
    infection_prob: float
    recovery_mean: float
    tau_sim: float


def _initial_infected_indices(n: int, k: int, rng: Generator) -> np.ndarray:
    k = min(k, n)
    return rng.choice(n, size=k, replace=False)


def run_sis_single(
    graph: ERGraph,
    infection_prob: float,
    recovery_mean: float,
    regime: RecoveryRegime,
    rng: Generator,
    *,
    initial_infected: int = 5,
    max_ticks: int = 10_000,
    late_window_len: int = 200,
    power_law_lambda: float = 4.24,
    lognormal_sigma: float = 1.0,
) -> SimulationResult:
    """
    One SIS trajectory. Infection then recovery each tick; prevalence sampled after dynamics
    (same relative order as NetLogo `go`).
    """
    n = graph.n
    neigh = graph.neighbors

    state = np.zeros(n, dtype=np.int8)  # 0=S, 1=I
    inf_time = np.zeros(n, dtype=np.int32)
    rec_deadline = np.zeros(n, dtype=np.int32)

    init_idx = _initial_infected_indices(n, initial_infected, rng)
    state[init_idx] = 1
    for i in init_idx:
        d = draw_recovery_time(
            regime, recovery_mean, rng, power_law_lambda=power_law_lambda, lognormal_sigma=lognormal_sigma
        )
        rec_deadline[i] = d
        inf_time[i] = 0

    prev_samples: list[float] = []
    beta = float(infection_prob)
    rm = float(recovery_mean)
    ticks = 0

    while ticks < max_ticks:
        if state.sum() == 0:
            break

        # infection-step: infected ascending index, neighbors sorted
        infected_order = np.flatnonzero(state == 1)
        for i in infected_order:
            for j in neigh[i]:
                if state[j] != 0:
                    continue
                if rng.random() < beta:
                    state[j] = 1
                    inf_time[j] = 0
                    rec_deadline[j] = draw_recovery_time(
                        regime,
                        recovery_mean,
                        rng,
                        power_law_lambda=power_law_lambda,
                        lognormal_sigma=lognormal_sigma,
                    )

        # recovery-step
        infected_order = np.flatnonzero(state == 1)
        for i in infected_order:
            inf_time[i] += 1
            if inf_time[i] >= rec_deadline[i]:
                state[i] = 0
                inf_time[i] = 0
                rec_deadline[i] = 0

        n_i = int(state.sum())
        frac = n_i / float(n) if n > 0 else 0.0
        prev_samples.append(frac)
        ticks += 1

        if int(state.sum()) == 0:
            break

    extinct = 1.0 if int(state.sum()) == 0 else 0.0
    w = late_window_len
    if len(prev_samples) == 0:
        late_prev = 0.0
    else:
        start = max(0, len(prev_samples) - w)
        late_prev = float(np.mean(prev_samples[start:]))

    return SimulationResult(
        bs_out_extinct=extinct,
        bs_out_final_tick=float(ticks),
        bs_out_late_mean_prevalence=late_prev,
        infection_prob=beta,
        recovery_mean=rm,
        tau_sim=beta * rm,
    )


def mean_degree(graph: ERGraph) -> float:
    if graph.n == 0:
        return 0.0
    return 2.0 * len(graph.edges) / graph.n
