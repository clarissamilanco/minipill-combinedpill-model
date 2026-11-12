import numpy as np
import pytest

from minipill.params import NORETHINDRONE, NORETHINDRONE_DOSE_MG, NORETHINDRONE_LABEL
from minipill.pk import (
    Drug, accumulation_ratio, bateman, ka_from_tmax, mg_per_l_to_pg_per_ml, simulate, tmax,
)

DAY = 24.0
# Generic levonorgestrel-like drug (30 h half-life) for checks that don't need real values.
SLOW = Drug(name="slow", half_life_h=30.0, ka=ka_from_tmax(1.5, 30.0), F=1.0, V_L=100.0)


@pytest.mark.parametrize("drug", [NORETHINDRONE, SLOW], ids=lambda d: d.name)
def test_single_dose_matches_bateman(drug):
    t = np.linspace(0, 3 * DAY, 2001)
    sim = simulate(drug, [0.0], 1.0, t)
    exact = bateman(drug, 1.0, t)
    np.testing.assert_allclose(sim.conc, exact, rtol=1e-6, atol=1e-9 * exact.max())


def test_multiple_doses_match_superposition():
    """Late and missed pills: the linear model's answer is a sum of single-dose curves."""
    doses = [0, 24, 51, 96, 120]  # pill 3 is 3 h late, pill 4 (72 h) skipped
    t = np.linspace(0, 7 * DAY, 3001)
    sim = simulate(NORETHINDRONE, doses, NORETHINDRONE_DOSE_MG, t)
    exact = sum(bateman(NORETHINDRONE, NORETHINDRONE_DOSE_MG, t - d) for d in doses)
    np.testing.assert_allclose(sim.conc, exact, rtol=1e-6, atol=1e-9 * exact.max())


@pytest.mark.parametrize("drug", [NORETHINDRONE, SLOW], ids=lambda d: d.name)
def test_mass_balance(drug):
    doses = np.arange(0, 28 * DAY, DAY)
    t = np.linspace(0, 30 * DAY, 4001)
    sim = simulate(drug, doses, 0.35, t)
    total = sim.gut + sim.central + sim.eliminated + sim.first_pass
    np.testing.assert_allclose(total, sim.dosed, rtol=1e-8, atol=1e-12)
    assert sim.dosed[-1] == pytest.approx(28 * 0.35)


def test_two_pills_at_once_is_a_double_dose():
    t = np.linspace(0, DAY, 201)
    double = simulate(NORETHINDRONE, [0, 0], 0.35, t).conc
    single = simulate(NORETHINDRONE, [0], 0.70, t).conc
    np.testing.assert_allclose(double, single, rtol=1e-10)


@pytest.mark.parametrize(
    "drug, expected", [(NORETHINDRONE, 1.14), (SLOW, 2.35)], ids=["norethindrone", "slow"]
)
def test_steady_state_accumulation(drug, expected):
    """Level just before a dose at steady state vs. 24 h after the first dose."""
    n_days = 40
    t = np.arange(0, n_days + 1) * DAY
    sim = simulate(drug, np.arange(n_days) * DAY, 1.0, t)
    ratio = sim.conc[-1] / sim.conc[1]
    assert ratio == pytest.approx(accumulation_ratio(drug, DAY), rel=1e-5)
    assert ratio == pytest.approx(expected, abs=0.01)  # README values


def test_ka_from_tmax_reproduces_tmax():
    assert tmax(NORETHINDRONE.ka, NORETHINDRONE.k) == pytest.approx(1.2, rel=1e-9)
    assert NORETHINDRONE.ka == pytest.approx(3.0, abs=0.1)


def test_peak_time_matches_label():
    t = np.linspace(0, 6, 6001)
    sim = simulate(NORETHINDRONE, [0], NORETHINDRONE_DOSE_MG, t)
    assert t[np.argmax(sim.conc)] == pytest.approx(NORETHINDRONE_LABEL["tmax_h"], abs=0.01)


def test_invalid_inputs():
    with pytest.raises(ValueError):
        simulate(NORETHINDRONE, [-1], 0.35, [0, 1])
    with pytest.raises(ValueError):
        simulate(NORETHINDRONE, [0], 0.35, [1, 0])
    with pytest.raises(ValueError):
        ka_from_tmax(100, 7.7)


@pytest.mark.xfail(
    strict=True,
    reason="Expected: one compartment gives ~8-fold peak-to-24h drop; label shows ~37-fold. "
    "Should pass once the second compartment is added.",
)
def test_norethindrone_peak_to_24h_drop_matches_label():
    t = np.linspace(0, DAY, 2401)
    c = mg_per_l_to_pg_per_ml(simulate(NORETHINDRONE, [0], NORETHINDRONE_DOSE_MG, t).conc)
    label_drop = NORETHINDRONE_LABEL["cmax_pg_ml"] / NORETHINDRONE_LABEL["c24_pg_ml"]
    assert c.max() / c[-1] == pytest.approx(label_drop, rel=0.25)
