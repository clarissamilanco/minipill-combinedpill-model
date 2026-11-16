"""Print what simulate() tracks, pill by pill, for three norethindrone pills.

A learning aid: shows that the drug left in the blood and tissues at the end of
one dosing interval carries over into the next (the new pill lands on top of it),
that a pill only jumps the gut amount, and that every row is mass-balanced.

Usage:  python scripts/show_simulate_state.py
Try changing PILLS, e.g. [0, 30, 48] for a 6 h late second pill.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from minipill.params import NORETHINDRONE, NORETHINDRONE_DOSE_MG  # noqa: E402
from minipill.pk import mg_per_l_to_pg_per_ml, simulate  # noqa: E402

PILLS = [0, 24, 48]  # hours
TIMES = [0, 1.2, 12, 23.999, 24, 25.2, 47.999, 48, 72]  # 23.999 = just before the next pill


def main():
    r = simulate(NORETHINDRONE, PILLS, NORETHINDRONE_DOSE_MG, TIMES)
    print(f"Norethindrone {NORETHINDRONE_DOSE_MG} mg pills at hours {PILLS}; amounts in mg\n")
    print(f"{'hour':>7} | {'gut':>7} {'blood':>7} {'tissue':>7} {'elim':>7} {'1st pass':>8} |"
          f" {'sum':>6} {'swallowed':>9} | {'conc pg/mL':>10}")
    print("-" * 92)
    for i, t in enumerate(r.t):
        total = r.gut[i] + r.central[i] + r.peripheral[i] + r.eliminated[i] + r.first_pass[i]
        print(f"{t:7.3f} | {r.gut[i]:7.4f} {r.central[i]:7.4f} {r.peripheral[i]:7.4f} {r.eliminated[i]:7.4f}"
              f" {r.first_pass[i]:8.4f} | {total:6.3f} {r.dosed[i]:9.2f} | {mg_per_l_to_pg_per_ml(r.conc[i]):10.0f}")


if __name__ == "__main__":
    main()
