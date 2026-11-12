"""One-compartment pharmacokinetic model with first-order oral absorption.

Units used throughout: time in hours, amounts in mg, volume in L,
so concentrations come out in mg/L (use ``mg_per_l_to_pg_per_ml`` to compare
with drug-label values).

State variables (all amounts, in mg):
    gut          drug still waiting to be absorbed            G
    central      drug in the blood                            A  (C = A / V)
    eliminated   drug cleared from the blood                  E
    first_pass   drug lost before reaching the blood (1 - F)  L

    dG/dt = -ka G
    dA/dt =  F ka G - k A
    dE/dt =  k A
    dL/dt = (1 - F) ka G

Tracking E and L as well as G and A means G + A + E + L always equals the total
dose taken, which gives an exact mass-balance check on the numerical solver.
"""

from collections import Counter
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

LN2 = np.log(2)


@dataclass(frozen=True)
class Drug:
    """Parameters for one hormone."""

    name: str
    half_life_h: float  # elimination half-life
    ka: float  # absorption rate constant, 1/h
    F: float  # oral bioavailability, 0-1
    V_L: float  # volume of distribution, L

    @property
    def k(self) -> float:
        """Elimination rate constant, 1/h."""
        return LN2 / self.half_life_h


@dataclass
class SimResult:
    t: np.ndarray  # h
    gut: np.ndarray  # mg
    central: np.ndarray  # mg
    eliminated: np.ndarray  # mg
    first_pass: np.ndarray  # mg
    conc: np.ndarray  # mg/L
    dosed: np.ndarray  # mg taken so far at each time point


def mg_per_l_to_pg_per_ml(c):
    """1 mg/L = 1e9 pg / 1e3 mL = 1e6 pg/mL."""
    return np.asarray(c) * 1e6


def mg_per_l_to_ng_per_ml(c):
    """1 mg/L = 1e3 ng/mL."""
    return np.asarray(c) * 1e3


def _rhs(t, y, ka, k, F):
    G, A, _E, _L = y
    absorbed = ka * G
    return [-absorbed, F * absorbed - k * A, k * A, (1 - F) * absorbed]


def simulate(drug: Drug, dose_times, dose_mg: float, t_eval, rtol=1e-10, atol=1e-14) -> SimResult:
    """Simulate blood levels for a list of pill times.

    Each pill adds ``dose_mg`` to the gut at its time. Late or missed pills are
    just different ``dose_times``. Two pills at the same time count as a double dose.

    Pills at or after the last time in ``t_eval`` are not simulated.
    """
    t_eval = np.asarray(t_eval, dtype=float)
    if t_eval.ndim != 1 or len(t_eval) == 0:
        raise ValueError("t_eval must be a non-empty 1-D array")
    if np.any(np.diff(t_eval) <= 0) or t_eval[0] < 0:
        raise ValueError("t_eval must be increasing and start at or after 0")
    dose_times = [float(d) for d in dose_times]
    if any(d < 0 for d in dose_times):
        raise ValueError("dose times must be >= 0")

    t_start, t_end = 0.0, t_eval[-1]
    doses_at = Counter(d for d in dose_times if d < t_end)
    breakpoints = sorted({t_start, t_end, *doses_at})

    y = np.zeros(4)
    total_dosed = 0.0
    out = np.full((4, len(t_eval)), np.nan)
    dosed = np.zeros(len(t_eval))

    for i, (a, b) in enumerate(zip(breakpoints[:-1], breakpoints[1:])):
        n = doses_at.get(a, 0)
        y[0] += n * dose_mg
        total_dosed += n * dose_mg

        last = i == len(breakpoints) - 2
        mask = (t_eval >= a) & ((t_eval <= b) if last else (t_eval < b))
        pts = np.unique(np.concatenate([t_eval[mask], [a, b]]))
        sol = solve_ivp(
            _rhs, (a, b), y, t_eval=pts, args=(drug.ka, drug.k, drug.F),
            method="LSODA", rtol=rtol, atol=atol,
        )
        if not sol.success:
            raise RuntimeError(f"ODE solver failed on [{a}, {b}]: {sol.message}")

        idx = np.searchsorted(pts, t_eval[mask])
        out[:, mask] = sol.y[:, idx]
        dosed[mask] = total_dosed
        y = sol.y[:, -1].copy()

    gut, central, eliminated, first_pass = out
    return SimResult(
        t=t_eval, gut=gut, central=central, eliminated=eliminated,
        first_pass=first_pass, conc=central / drug.V_L, dosed=dosed,
    )


def bateman(drug: Drug, dose_mg: float, t):
    """Exact blood concentration (mg/L) after one dose at t = 0 (needs ka != k)."""
    t = np.asarray(t, dtype=float)
    ka, k = drug.ka, drug.k
    if np.isclose(ka, k):
        raise ValueError("Bateman formula needs ka != k")
    c = drug.F * dose_mg * ka / (drug.V_L * (ka - k)) * (np.exp(-k * t) - np.exp(-ka * t))
    return np.where(t >= 0, c, 0.0)


def tmax(ka: float, k: float) -> float:
    """Time of peak concentration after one dose."""
    return np.log(ka / k) / (ka - k)


def ka_from_tmax(tmax_h: float, half_life_h: float) -> float:
    """Absorption rate that puts the single-dose peak at ``tmax_h``.

    Picks the solution with ka > k (absorption faster than elimination), which
    is the case for every drug in this project.
    """
    k = LN2 / half_life_h
    if tmax_h >= 1 / k:
        raise ValueError(f"tmax must be < 1/k = {1 / k:.2f} h for ka > k")
    return brentq(lambda ka: tmax(ka, k) - tmax_h, k * (1 + 1e-9), 1e4)


def accumulation_ratio(drug: Drug, tau_h: float) -> float:
    """Steady-state / first-dose ratio for doses every ``tau_h`` hours, 1 / (1 - e^(-k tau))."""
    return 1 / (1 - np.exp(-drug.k * tau_h))
