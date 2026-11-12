"""Fit levonorgestrel's two-compartment parameters to Alesse day 1 (no SHBG effect).

Day 1 is a single dose, before ethinyl estradiol has raised SHBG, so it is the
only Alesse data a model without SHBG should match.

Fixed:  F = 1 (Alesse label: "bioavailability about 100%"), dose = 0.10 mg,
        ka = 0.8 /h (ASSUMED, see below).
Fitted (4 unknowns): V1, k10, k12, k21.
Targets (4, Alesse Table I, day 1):
    Cmax 2.75 ng/mL, Tmax 1.6 h, AUC 35.2 ng*h/mL,
    terminal half-life 34.3 h = ln 2 * (Vz/F) / (CL/F) = ln 2 * 2.66 / 0.0537
    (the per-kg body weight cancels, so no weight assumption is needed).

Four targets can't pin down five rates. Every ka from about 0.4 to 1.2 /h fits
all four exactly, and they all give nearly the same 24 h level and late-pill
drop; ``--profile`` prints that table. ka = 0.8 /h is a value in the middle.

Usage:  python scripts/fit_levonorgestrel.py [--profile]
Copy the printed values into src/minipill/params.py.
"""

import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares, minimize_scalar

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from minipill.params import ALESSE_DAY1, LEVONORGESTREL_DOSE_MG  # noqa: E402
from minipill.pk import LN2, Drug, exact_single_dose, mg_per_l_to_ng_per_ml, simulate  # noqa: E402

F = 1.0
KA = 0.8
DOSE = LEVONORGESTREL_DOSE_MG
TARGETS = {
    "cmax_ng_ml": ALESSE_DAY1["cmax_ng_ml"],
    "tmax_h": ALESSE_DAY1["tmax_h"],
    "auc_ng_h_ml": ALESSE_DAY1["auc_ng_h_ml"],
    "half_life_h": ALESSE_DAY1["half_life_h"],
}


def summarize(drug: Drug) -> dict:
    def c(t):
        return float(mg_per_l_to_ng_per_ml(exact_single_dose(drug, DOSE, t)))

    peak = minimize_scalar(lambda t: -c(t), bounds=(0.05, 8), method="bounded", options={"xatol": 1e-7})
    return {
        "cmax_ng_ml": -peak.fun,
        "tmax_h": peak.x,
        "auc_ng_h_ml": float(mg_per_l_to_ng_per_ml(drug.F * DOSE / drug.clearance_L_per_h)),  # exact AUC(0-inf)
        "half_life_h": drug.terminal_half_life_h,
    }


def make_drug(log_params, ka=KA) -> Drug:
    V, k10, k12, k21 = np.exp(log_params)
    return Drug(name="levonorgestrel", ka=ka, F=F, V_L=V, k10=k10, k12=k12, k21=k21)


def fit(ka=KA, n_starts=40, seed=0):
    def residuals(p):
        s = summarize(make_drug(p, ka))
        return [np.log(s[k] / v) for k, v in TARGETS.items()]

    best = None
    rng = np.random.default_rng(seed)
    for _ in range(n_starts):  # several starting points; keep the best fit
        x0 = np.log([25.0, 0.12, 0.3, 0.06]) + rng.normal(0, 0.8, 4)
        try:
            f = least_squares(residuals, x0, xtol=1e-14, ftol=1e-14, gtol=1e-14)
        except ValueError:  # repeated rates: exact solution undefined
            continue
        if best is None or f.cost < best.cost:
            best = f
    return make_drug(best.x, ka), best.cost


def late_pill_drop(drug: Drug) -> tuple[float, float]:
    """Trough before pill 21 (on time) and the level 48 h after pill 20 (pill 21 missed)."""
    doses = np.arange(20) * 24.0
    t = np.array([20 * 24.0, 21 * 24.0])
    c = mg_per_l_to_ng_per_ml(simulate(drug, doses, DOSE, t).conc)
    return float(c[0]), float(c[1])


def main():
    if "--profile" in sys.argv:
        print("  ka    V1     k10    k12    k21   | C24    day-20 trough  48 h after last pill")
        for ka in (0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.5):
            drug, cost = fit(ka)
            if cost > 1e-10:
                print(f"  {ka:3.1f}  no exact fit")
                continue
            c24 = float(mg_per_l_to_ng_per_ml(exact_single_dose(drug, DOSE, 24.0)))
            trough, late = late_pill_drop(drug)
            print(f"  {ka:3.1f}  {drug.V_L:5.1f}  {drug.k10:.3f}  {drug.k12:.3f}  {drug.k21:.3f} | {c24:.3f}  "
                  f"{trough:.3f}          {late:.3f} ({100 * (1 - late / trough):.0f}% below)")
        return

    drug, cost = fit()
    s = summarize(drug)
    print(f"Fitted parameters (ka fixed at {KA} /h, F = {F}):")
    for name in ("V_L", "k10", "k12", "k21"):
        print(f"  {name:4s} = {getattr(drug, name):.6g}")
    print(f"\nDistribution half-life (alpha): {LN2 / drug.disposition_rates[0]:.3g} h")
    print(f"Clearance: {drug.clearance_L_per_h:.2f} L/h   Vss: {drug.vss_L:.0f} L")
    print("\nTarget            label       model")
    for k, v in TARGETS.items():
        print(f"  {k:15s} {v:9.4g}   {s[k]:9.4g}")


if __name__ == "__main__":
    main()
