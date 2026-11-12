"""Fit norethindrone's two-compartment parameters to the label's single-dose curve.

Fixed from the label: F = 0.65, dose = 0.35 mg.
Fitted (5 unknowns): ka, V1, k10, k12, k21.
Targets (5, norethindrone label, single 0.35 mg dose):
    Tmax 1.2 h, Cmax 4,817 pg/mL, C(24 h) 130 pg/mL,
    Cavg over 24 h 885 pg/mL, terminal half-life 7.7 h.

Usage:  python scripts/fit_norethindrone.py
Copy the printed values into src/minipill/params.py.
"""

import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares, minimize_scalar

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from minipill.params import BODY_WEIGHT_KG, NORETHINDRONE_DOSE_MG, NORETHINDRONE_LABEL  # noqa: E402
from minipill.pk import Drug, exact_single_dose, mg_per_l_to_pg_per_ml  # noqa: E402

F = 0.65
TARGET_HALF_LIFE_H = 7.7
T24 = np.linspace(0, 24, 4801)


def summarize(drug: Drug) -> dict:
    def c(t):
        return float(mg_per_l_to_pg_per_ml(exact_single_dose(drug, NORETHINDRONE_DOSE_MG, t)))

    peak = minimize_scalar(lambda t: -c(t), bounds=(0.05, 6), method="bounded", options={"xatol": 1e-6})
    curve = mg_per_l_to_pg_per_ml(exact_single_dose(drug, NORETHINDRONE_DOSE_MG, T24))
    return {
        "tmax_h": peak.x,
        "cmax_pg_ml": -peak.fun,
        "c24_pg_ml": c(24.0),
        "cavg_24h_pg_ml": np.trapezoid(curve, T24) / 24,
        "half_life_h": drug.terminal_half_life_h,
    }


def make_drug(log_params) -> Drug:
    ka, V, k10, k12, k21 = np.exp(log_params)
    return Drug(name="norethindrone", ka=ka, F=F, V_L=V, k10=k10, k12=k12, k21=k21)


TARGETS = {**{k: NORETHINDRONE_LABEL[k] for k in ("tmax_h", "cmax_pg_ml", "c24_pg_ml", "cavg_24h_pg_ml")},
           "half_life_h": TARGET_HALF_LIFE_H}


def residuals(log_params):
    s = summarize(make_drug(log_params))
    return [np.log(s[k] / TARGETS[k]) for k in TARGETS]


def main():
    best = None
    rng = np.random.default_rng(0)
    for _ in range(20):  # several starting points; keep the best fit
        x0 = np.log([3.0, 40.0, 0.3, 0.5, 0.15]) + rng.normal(0, 0.7, 5)
        try:
            fit = least_squares(residuals, x0, xtol=1e-14, ftol=1e-14, gtol=1e-14)
        except ValueError:
            continue
        if best is None or fit.cost < best.cost:
            best = fit
    drug = make_drug(best.x)
    s = summarize(drug)

    print("Fitted parameters:")
    for name in ("ka", "V_L", "k10", "k12", "k21"):
        print(f"  {name:4s} = {getattr(drug, name):.6g}")
    print(f"\nDistribution half-life (alpha): {np.log(2) / drug.disposition_rates[0]:.3g} h")
    print(f"Clearance:  {drug.clearance_L_per_h * 24:.0f} L/day   (label: ~600 L/day)")
    print(f"Vss:        {drug.vss_L:.0f} L = {drug.vss_L / BODY_WEIGHT_KG:.2f} L/kg at {BODY_WEIGHT_KG:.0f} kg"
          "   (label: ~4 L/kg)")
    print("\nTarget            label       model")
    for k, v in TARGETS.items():
        print(f"  {k:15s} {v:9.4g}   {s[k]:9.4g}")


if __name__ == "__main__":
    main()
