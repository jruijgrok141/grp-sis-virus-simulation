"""Recovery-time draws matching NetLogo `draw-recovery-time` (virus_simulation.nlogox)."""

from __future__ import annotations

from typing import Literal

import numpy as np
from numpy.random import Generator

RecoveryRegime = Literal["exponential", "power_law_tang", "lognormal_tang"]


def draw_recovery_time(
    regime: RecoveryRegime,
    recovery_mean: float,
    rng: Generator,
    *,
    power_law_lambda: float = 4.24,
    lognormal_sigma: float = 1.0,
) -> int:
    """
    Integer recovery horizon in ticks, >= 1.
    Mirrors NetLogo: exponential / power law (Tang) / lognormal (Tang).
    """
    rm = float(recovery_mean)
    if regime == "exponential":
        mu = 1.0 / rm
        u = rng.random()
        w = -np.log(max(1e-6, 1.0 - u)) / mu
        return int(max(1, round(w)))

    if regime == "lognormal_tang":
        sig = max(1e-6, float(lognormal_sigma))
        mu_ln = np.log(rm) - (sig * sig) / 2.0
        w = np.exp(mu_ln + sig * rng.standard_normal())
        return int(max(1, round(w)))

    if regime == "power_law_tang":
        lam = max(2.01, float(power_law_lambda))
        t0 = rm * (lam - 2.0) / (lam - 1.0)
        v = max(1e-12, rng.random())
        t = t0 * (v ** (1.0 / (1.0 - lam)))
        return int(max(1, round(t)))

    raise ValueError(f"Unknown regime: {regime}")
