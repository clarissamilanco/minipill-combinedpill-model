"""Calibrate the SHBG effect on levonorgestrel to Alesse Table I.

Step 1, SHBG over time (3 values, 3 measurements, exact):
    S(t) = S_max - (S_max - S0) exp(-kS t) while pills are taken,
    through 57, 81, 93 nmol/L on days 1, 6, 21 (t = 0, 5, 20 days).

Step 2, how strongly SHBG slows clearance (1 value, 2 measurements):
    k10(t) = k10 (S0 / S)^gamma, gamma fitted to the measured AUC over one
    dosing interval on day 6 (46.0) and day 21 (68.3 ng*h/mL).

The levonorgestrel parameters (fitted to day 1, before SHBG rises) stay fixed.
Not used in the fit, so they are checks: day-6 and day-21 peaks (4.52, 6.0 ng/mL)
and the steady-state low before the next pill (1.9 ng/mL).

Usage:  python scripts/fit_shbg.py
Copy the printed values into src/minipill/params.py.
"""

import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, least_squares

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from minipill.params import (  # noqa: E402
    ALESSE_SHBG_NMOL_L, ALESSE_STEADY_STATE, ALESSE_TABLE_I, LEVONORGESTREL, LEVONORGESTREL_DOSE_MG,
)
from minipill.pk import SHBGEffect, mg_per_l_to_ng_per_ml, simulate  # noqa: E402

DAY = 24.0


def fit_shbg_curve():
    s0, s6, s21 = (ALESSE_SHBG_NMOL_L[d] for d in (1, 6, 21))
    t6, t21 = 5.0, 20.0  # days after the first pill

    # With r = exp(-kS t6), the day-21 equation becomes one equation in S_max.
    def mismatch(s_max):
        r = (s_max - s6) / (s_max - s0)
        return (s_max - s21) / (s_max - s0) - r ** (t21 / t6)

    s_max = brentq(mismatch, s21 + 1e-9, s21 + 50)
    k_per_day = -np.log((s_max - s6) / (s_max - s0)) / t6
    return s0, s_max, k_per_day / DAY


def day_metrics(shbg, days=(1, 6, 21)):
    """AUC over the day's dosing interval, peak, and level just before the next pill."""
    t = np.linspace(0, 21 * DAY, 21 * 240 + 1)
    sim = simulate(LEVONORGESTREL, np.arange(21) * DAY, LEVONORGESTREL_DOSE_MG, t, shbg=shbg)
    c = mg_per_l_to_ng_per_ml(sim.conc)
    out = {}
    for d in days:
        m = (t >= (d - 1) * DAY) & (t <= d * DAY)
        out[d] = {"auc": np.trapezoid(c[m], t[m]), "cmax": c[m].max(), "cmin": c[m][-1]}
    return out, t, sim.shbg


def main():
    s0, s_max, k = fit_shbg_curve()
    print(f"SHBG: baseline {s0:.0f}, plateau {s_max:.2f} nmol/L, rate {k * DAY:.4f} /day "
          f"(half-way in {np.log(2) / (k * DAY):.1f} days)")

    def residuals(p):
        m, _, _ = day_metrics(SHBGEffect(s0, s_max, k, p[0]), days=(6, 21))
        return [np.log(m[d]["auc"] / ALESSE_TABLE_I[d][0]) for d in (6, 21)]

    gamma = least_squares(residuals, [1.2], bounds=(0, 5)).x[0]
    shbg = SHBGEffect(s0, s_max, k, gamma)
    print(f"Strength gamma = {gamma:.4f}\n")

    m, t, S = day_metrics(shbg)
    peaks = {1: 2.75, 6: 4.52, 21: ALESSE_STEADY_STATE["cmax_ng_ml"]}
    print("Day   SHBG model/label   AUC model/label   peak model/label   low before next pill")
    for d in (1, 6, 21):
        s_model = S[np.searchsorted(t, (d - 1) * DAY)]
        low = f"{m[d]['cmin']:.2f} / {ALESSE_STEADY_STATE['cmin_ng_ml']}" if d == 21 else f"{m[d]['cmin']:.2f}"
        auc_label = f"{ALESSE_TABLE_I[d][0]:<5}" if d != 1 else "(n/a)"  # day-1 label AUC is to infinity
        print(f"{d:3d}   {s_model:5.1f} / {ALESSE_SHBG_NMOL_L[d]:<5}      {m[d]['auc']:5.1f} / {auc_label}"
              f"     {m[d]['cmax']:4.2f} / {peaks[d]:<5}       {low}")


if __name__ == "__main__":
    main()
