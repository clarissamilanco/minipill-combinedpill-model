"""Write the model's curves into app/pill_comparison.html, for a grid of BMIs.

With SHBG the model has no simple formula, so the page can't recompute it in the
browser. This script runs the Python model and replaces the data block between
the /*DATA*/ and /*END DATA*/ markers in the page.

For each BMI on the grid (clearance scaled by ``clearance_factor_for_bmi``):
 - single-dose terms (rate, amplitude in ng/mL) for both two-compartment drugs,
   which the page sums for any pill schedule;
 - the combined pill with SHBG over a month, and for every late-pill option.
The page blends between the two nearest grid BMIs. SHBG doesn't depend on
clearance, so its curve is stored once.

Usage:  python scripts/export_comparison_data.py
"""

import json
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from minipill.params import (  # noqa: E402
    BMI_OBESE, BMI_REFERENCE, LEVONORGESTREL, LEVONORGESTREL_DOSE_MG, NORETHINDRONE, NORETHINDRONE_DOSE_MG, SHBG,
    clearance_factor_for_bmi,
)
from minipill.pk import Drug, mg_per_l_to_ng_per_ml, simulate, single_dose_terms  # noqa: E402

DAY = 24.0
PAGE = ROOT / "app" / "pill_comparison.html"
PACK = [d * DAY for d in range(21)]
DUE = 14 * DAY  # the day-15 pill
MONTH_STEP, LATE_STEP = 0.2, 0.25  # hours
LATE_WINDOW = (-24.0, 96.0)  # hours around DUE
LATE_OPTIONS = [x / 2 for x in range(48)] + ["skip", "skip2"]  # 0, 0.5, ... 23.5 h late
BMI_GRID = np.linspace(BMI_REFERENCE, BMI_OBESE, 9)  # clearance is flat outside this range


def run(drug: Drug, doses, t):
    sim = simulate(drug, sorted(doses), LEVONORGESTREL_DOSE_MG, t, shbg=SHBG)
    return mg_per_l_to_ng_per_ml(sim.conc), sim.shbg


def terms(drug: Drug, dose_mg: float):
    rates, amps = single_dose_terms(drug, dose_mg)
    order = np.argsort(rates)
    return [[float(rates[i]), float(mg_per_l_to_ng_per_ml(amps[i]))] for i in order]


def scenario(late):
    if late == "skip":
        return [d for d in PACK if d != DUE]
    if late == "skip2":  # CDC: pills 15 and 16 missed; missed pill + pill 17 taken together at +48 h
        return [d for d in PACK if d not in (DUE, DUE + DAY)] + [DUE + 2 * DAY]
    return [d for d in PACK if d != DUE] + [DUE + late]


def combined_pill(lng: Drug):
    """Month curve, usual day-15 low and every late-pill option, with SHBG."""
    t_month = np.arange(0, 28 * DAY + 1e-9, MONTH_STEP)
    c_month, s_month = run(lng, PACK, t_month)

    t_late = np.arange(0, DUE + LATE_WINDOW[1] + 1e-9, LATE_STEP)
    c_on, _ = run(lng, PACK, t_late)
    low = c_on[int(round(DUE / LATE_STEP))]
    first = int(round((DUE + LATE_WINDOW[0]) / LATE_STEP))
    late = {}
    for opt in LATE_OPTIONS:
        c, _ = run(lng, scenario(opt), t_late)
        late[opt if isinstance(opt, str) else f"{opt:g}"] = [round(100 * v / low, 1) for v in c[first:]]
    return c_month, s_month, low, late


def main():
    data = {"bmiGrid": [round(float(b), 4) for b in BMI_GRID], "factors": [], "net2": [], "lng2": [],
            "monthStep": MONTH_STEP, "month": [], "lateStep": LATE_STEP, "lateStart": LATE_WINDOW[0],
            "usualLow": [], "late": []}
    for b in BMI_GRID:
        f = float(clearance_factor_for_bmi(b))
        net, lng = NORETHINDRONE.with_clearance_scaled(f), LEVONORGESTREL.with_clearance_scaled(f)
        c_month, s_month, low, late = combined_pill(lng)
        data["factors"].append(round(f, 4))
        data["net2"].append(terms(net, NORETHINDRONE_DOSE_MG))
        data["lng2"].append(terms(lng, LEVONORGESTREL_DOSE_MG))
        data["month"].append([round(float(v), 4) for v in c_month])
        data["usualLow"].append(round(float(low), 4))
        data["late"].append(late)
        data.setdefault("shbg", [round(float(v), 2) for v in s_month])  # same at every BMI
        print(f"BMI {b:5.2f}: clearance x{f:.3f}, combined-pill usual low {low:.3f} ng/mL")

    block = "/*DATA*/ const MODEL = " + json.dumps(data, separators=(",", ":")) + "; /*END DATA*/"
    html = PAGE.read_text(encoding="utf-8")
    new, n = re.subn(r"/\*DATA\*/.*?/\*END DATA\*/", lambda _: block, html, flags=re.S)
    if n != 1:
        raise SystemExit("data markers not found in page")
    PAGE.write_text(new, encoding="utf-8")
    print(f"Wrote {len(block) / 1024:.0f} KB of data for {len(BMI_GRID)} BMIs")


if __name__ == "__main__":
    main()
