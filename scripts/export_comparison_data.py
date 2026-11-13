"""Write the combined pill's SHBG curves into app/pill_comparison.html.

With SHBG the model has no simple formula, so the page can't recompute it in the
browser. This script runs the Python model and replaces the data block between
the /*DATA*/ and /*END DATA*/ markers in the page.

Usage:  python scripts/export_comparison_data.py
"""

import json
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from minipill.params import LEVONORGESTREL, LEVONORGESTREL_DOSE_MG, SHBG  # noqa: E402
from minipill.pk import mg_per_l_to_ng_per_ml, simulate  # noqa: E402

DAY = 24.0
PAGE = ROOT / "app" / "pill_comparison.html"
PACK = [d * DAY for d in range(21)]
DUE = 14 * DAY  # the day-15 pill
MONTH_STEP, LATE_STEP = 0.2, 0.25  # hours
LATE_WINDOW = (-24.0, 96.0)  # hours around DUE
LATE_OPTIONS = [x / 2 for x in range(48)] + ["skip", "skip2"]  # 0, 0.5, ... 23.5 h late


def run(doses, t):
    sim = simulate(LEVONORGESTREL, sorted(doses), LEVONORGESTREL_DOSE_MG, t, shbg=SHBG)
    return mg_per_l_to_ng_per_ml(sim.conc), sim.shbg


def scenario(late):
    if late == "skip":
        return [d for d in PACK if d != DUE]
    if late == "skip2":  # CDC: pills 15 and 16 missed; missed pill + pill 17 taken together at +48 h
        return [d for d in PACK if d not in (DUE, DUE + DAY)] + [DUE + 2 * DAY]
    return [d for d in PACK if d != DUE] + [DUE + late]


def main():
    t_month = np.arange(0, 28 * DAY + 1e-9, MONTH_STEP)
    c_month, s_month = run(PACK, t_month)

    t_late = np.arange(0, DUE + LATE_WINDOW[1] + 1e-9, LATE_STEP)
    c_on, _ = run(PACK, t_late)
    low = c_on[int(round(DUE / LATE_STEP))]
    first = int(round((DUE + LATE_WINDOW[0]) / LATE_STEP))
    late = {}
    for opt in LATE_OPTIONS:
        c, _ = run(scenario(opt), t_late)
        late[opt if isinstance(opt, str) else f"{opt:g}"] = [round(100 * v / low, 1) for v in c[first:]]

    data = {
        "monthStep": MONTH_STEP,
        "month": [round(float(v), 4) for v in c_month],
        "shbg": [round(float(v), 2) for v in s_month],
        "lateStep": LATE_STEP,
        "lateStart": LATE_WINDOW[0],
        "usualLow": round(float(low), 4),
        "late": late,
    }
    block = "/*DATA*/ const LNG_SHBG = " + json.dumps(data, separators=(",", ":")) + "; /*END DATA*/"
    html = PAGE.read_text(encoding="utf-8")
    new, n = re.subn(r"/\*DATA\*/.*?/\*END DATA\*/", lambda _: block, html, flags=re.S)
    if n != 1:
        raise SystemExit("data markers not found in page")
    PAGE.write_text(new, encoding="utf-8")
    print(f"Wrote {len(block) / 1024:.0f} KB of data; usual low {low:.3f} ng/mL")


if __name__ == "__main__":
    main()
