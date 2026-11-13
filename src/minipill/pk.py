"""Compartment pharmacokinetic model with first-order oral absorption.

One or two compartments: with k12 = k21 = 0 the peripheral compartment is
never reached and this is the one-compartment model.

Units used throughout: time in hours, amounts in mg, volume in L,
so concentrations come out in mg/L (use ``mg_per_l_to_pg_per_ml`` to compare
with drug-label values).

State variables (all amounts, in mg):
    gut          drug still waiting to be absorbed            G
    central      drug in the blood                            A1  (C = A1 / V)
    peripheral   drug that has moved into the tissues         A2
    eliminated   drug cleared from the blood                  E
    first_pass   drug lost before reaching the blood (1 - F)  L

    dG/dt  = -ka G
    dA1/dt =  F ka G - (k10 + k12) A1 + k21 A2
    dA2/dt =  k12 A1 - k21 A2
    dE/dt  =  k10 A1
    dL/dt  = (1 - F) ka G

Tracking E and L as well as G, A1 and A2 means their sum always equals the
total dose taken, which gives an exact mass-balance check on the numerical solver.

Optional SHBG effect (combined pill only), a sixth state S in nmol/L:

    dS/dt  = kS (S_target - S),   S_target = S_max while estrogen is on, else S0
    k10(t) = k10 (S0 / S)^gamma

Estrogen counts as "on" for 24 h after each active pill. With the effect on,
the model is no longer linear, so only the step-by-step ``simulate`` handles it.
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
    ka: float  # absorption rate constant, 1/h
    F: float  # oral bioavailability, 0-1
    V_L: float  # volume of the central (blood) compartment, L
    k10: float  # elimination rate from the central compartment, 1/h
    k12: float = 0.0  # central -> peripheral, 1/h
    k21: float = 0.0  # peripheral -> central, 1/h

    def __post_init__(self):
        if min(self.ka, self.F, self.V_L, self.k10) <= 0 or self.F > 1:
            raise ValueError("ka, V_L, k10 must be > 0 and 0 < F <= 1")
        if self.k12 < 0 or self.k21 < 0 or (self.k12 > 0) != (self.k21 > 0):
            raise ValueError("k12 and k21 must both be 0 (one compartment) or both > 0")

    @classmethod
    def one_compartment(cls, name, half_life_h, ka, F, V_L):
        return cls(name=name, ka=ka, F=F, V_L=V_L, k10=LN2 / half_life_h)

    @property
    def n_compartments(self) -> int:
        return 2 if self.k12 > 0 else 1

    @property
    def disposition_rates(self) -> tuple[float, ...]:
        """Rates of the exponential phases after absorption: (k10,) or (alpha, beta), alpha > beta."""
        if self.n_compartments == 1:
            return (self.k10,)
        s = self.k10 + self.k12 + self.k21
        root = np.sqrt(s * s - 4 * self.k10 * self.k21)
        return ((s + root) / 2, (s - root) / 2)

    @property
    def terminal_rate(self) -> float:
        """Slowest elimination rate, 1/h (sets the terminal half-life)."""
        return self.disposition_rates[-1]

    @property
    def terminal_half_life_h(self) -> float:
        return LN2 / self.terminal_rate

    @property
    def clearance_L_per_h(self) -> float:
        return self.k10 * self.V_L

    @property
    def vss_L(self) -> float:
        """Volume of distribution at steady state, V1 (1 + k12 / k21)."""
        return self.V_L * (1 + self.k12 / self.k21) if self.k12 > 0 else self.V_L


@dataclass(frozen=True)
class SHBGEffect:
    """SHBG rises while the pill's estrogen is on and slows the drug's clearance."""

    s0: float  # baseline SHBG, nmol/L
    s_max: float  # level SHBG approaches with daily estrogen, nmol/L
    k_per_h: float  # rate SHBG moves toward its target, 1/h
    gamma: float  # strength: clearance scales as (s0 / S)^gamma

    def __post_init__(self):
        if min(self.s0, self.s_max, self.k_per_h) <= 0 or self.gamma < 0:
            raise ValueError("s0, s_max, k_per_h must be > 0 and gamma >= 0")

    def clearance_factor(self, S):
        return (self.s0 / np.asarray(S)) ** self.gamma


@dataclass
class SimResult:
    t: np.ndarray  # h
    gut: np.ndarray  # mg
    central: np.ndarray  # mg
    peripheral: np.ndarray  # mg
    eliminated: np.ndarray  # mg
    first_pass: np.ndarray  # mg
    conc: np.ndarray  # mg/L
    dosed: np.ndarray  # mg taken so far at each time point
    shbg: np.ndarray | None = None  # nmol/L, only with an SHBG effect


def mg_per_l_to_pg_per_ml(c):
    """1 mg/L = 1e9 pg / 1e3 mL = 1e6 pg/mL."""
    return np.asarray(c) * 1e6


def mg_per_l_to_ng_per_ml(c):
    """1 mg/L = 1e3 ng/mL."""
    return np.asarray(c) * 1e3


def _rhs(t, y, ka, F, k10, k12, k21):
    G, A1, A2, _E, _L = y
    absorbed = ka * G
    return [
        -absorbed,
        F * absorbed - (k10 + k12) * A1 + k21 * A2,
        k12 * A1 - k21 * A2,
        k10 * A1,
        (1 - F) * absorbed,
    ]


def _rhs_shbg(t, y, ka, F, k10, k12, k21, shbg: SHBGEffect, estrogen_on: bool):
    G, A1, A2, _E, _L, S = y
    absorbed = ka * G
    k10_t = k10 * shbg.clearance_factor(S)
    target = shbg.s_max if estrogen_on else shbg.s0
    return [
        -absorbed,
        F * absorbed - (k10_t + k12) * A1 + k21 * A2,
        k12 * A1 - k21 * A2,
        k10_t * A1,
        (1 - F) * absorbed,
        shbg.k_per_h * (target - S),
    ]


ESTROGEN_HOURS = 24.0  # estrogen counts as "on" for this long after each active pill


def simulate(drug: Drug, dose_times, dose_mg: float, t_eval, rtol=1e-10, atol=1e-14,
             shbg: SHBGEffect | None = None) -> SimResult:
    """Simulate blood levels for a list of pill times.

    Each pill adds ``dose_mg`` to the gut at its time. Late or missed pills are
    just different ``dose_times``. Two pills at the same time count as a double dose.

    With ``shbg``, every pill also carries estrogen, which pushes SHBG up for
    24 h after it is taken; SHBG starts at baseline (``shbg.s0``).

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
    off_times = {d + ESTROGEN_HOURS for d in doses_at if d + ESTROGEN_HOURS < t_end} if shbg else set()
    breakpoints = sorted({t_start, t_end, *doses_at, *off_times})
    args = (drug.ka, drug.F, drug.k10, drug.k12, drug.k21)

    n_state = 6 if shbg else 5
    y = np.zeros(n_state)
    if shbg:
        y[5] = shbg.s0
    total_dosed = 0.0
    out = np.full((n_state, len(t_eval)), np.nan)
    dosed = np.zeros(len(t_eval))

    for i, (a, b) in enumerate(zip(breakpoints[:-1], breakpoints[1:])):
        n = doses_at.get(a, 0)
        y[0] += n * dose_mg
        total_dosed += n * dose_mg

        last = i == len(breakpoints) - 2
        mask = (t_eval >= a) & ((t_eval <= b) if last else (t_eval < b))
        pts = np.unique(np.concatenate([t_eval[mask], [a, b]]))
        if shbg:
            on = any(d <= a < d + ESTROGEN_HOURS for d in doses_at)
            sol = solve_ivp(_rhs_shbg, (a, b), y, t_eval=pts, args=(*args, shbg, on),
                            method="LSODA", rtol=rtol, atol=atol)
        else:
            sol = solve_ivp(_rhs, (a, b), y, t_eval=pts, args=args, method="LSODA", rtol=rtol, atol=atol)
        if not sol.success:
            raise RuntimeError(f"ODE solver failed on [{a}, {b}]: {sol.message}")

        idx = np.searchsorted(pts, t_eval[mask])
        out[:, mask] = sol.y[:, idx]
        dosed[mask] = total_dosed
        y = sol.y[:, -1].copy()

    gut, central, peripheral, eliminated, first_pass = out[:5]
    return SimResult(
        t=t_eval, gut=gut, central=central, peripheral=peripheral, eliminated=eliminated,
        first_pass=first_pass, conc=central / drug.V_L, dosed=dosed,
        shbg=out[5] if shbg else None,
    )


def bateman(drug: Drug, dose_mg: float, t):
    """Exact one-compartment concentration (mg/L) after one dose at t = 0 (needs ka != k10)."""
    if drug.n_compartments != 1:
        raise ValueError("bateman() is for one-compartment drugs; use exact_single_dose()")
    t = np.asarray(t, dtype=float)
    ka, k = drug.ka, drug.k10
    if np.isclose(ka, k):
        raise ValueError("Bateman formula needs ka != k10")
    c = drug.F * dose_mg * ka / (drug.V_L * (ka - k)) * (np.exp(-k * t) - np.exp(-ka * t))
    return np.where(t >= 0, c, 0.0)


def exact_single_dose(drug: Drug, dose_mg: float, t):
    """Exact concentration (mg/L) after one dose at t = 0, for one or two compartments.

    The model is linear, y' = M y, so y(t) = P exp(D t) P^-1 y0 using the
    eigendecomposition of M. Independent of the ODE solver used in ``simulate``.
    """
    t = np.asarray(t, dtype=float)
    ka, F, k10, k12, k21 = drug.ka, drug.F, drug.k10, drug.k12, drug.k21
    M = np.array([
        [-ka, 0.0, 0.0],
        [F * ka, -(k10 + k12), k21],
        [0.0, k12, -k21],
    ])
    lam, P = np.linalg.eig(M)
    if np.min(np.abs(np.subtract.outer(lam, lam)) + np.eye(3)) < 1e-9:
        raise ValueError("repeated rates: exact solution needs distinct ka, alpha, beta")
    coef = np.linalg.solve(P, [dose_mg, 0.0, 0.0])
    tt = np.clip(t, 0, None).ravel()
    central = ((P[1] * coef) @ np.exp(np.outer(lam, tt))).reshape(t.shape)
    return np.where(t >= 0, np.real(central) / drug.V_L, 0.0)


def tmax(ka: float, k: float) -> float:
    """One-compartment time of peak concentration after one dose."""
    return np.log(ka / k) / (ka - k)


def ka_from_tmax(tmax_h: float, half_life_h: float) -> float:
    """One-compartment absorption rate that puts the single-dose peak at ``tmax_h``.

    Picks the solution with ka > k (absorption faster than elimination), which
    is the case for every drug in this project.
    """
    k = LN2 / half_life_h
    if tmax_h >= 1 / k:
        raise ValueError(f"tmax must be < 1/k = {1 / k:.2f} h for ka > k")
    return brentq(lambda ka: tmax(ka, k) - tmax_h, k * (1 + 1e-9), 1e4)


def accumulation_ratio(drug: Drug, tau_h: float) -> float:
    """Steady-state / first-dose ratio for doses every ``tau_h`` hours, 1 / (1 - e^(-k tau)).

    Exact for one compartment (k = k10). For two compartments, uses the terminal
    rate, which is approximate.
    """
    return 1 / (1 - np.exp(-drug.terminal_rate * tau_h))
