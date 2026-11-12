import numpy as np
import pytest

from minipill.params import (
    NORETHINDRONE, NORETHINDRONE_1C, NORETHINDRONE_DOSE_MG, NORETHINDRONE_LABEL,
)
from minipill.pk import (
    Drug, accumulation_ratio, bateman, exact_single_dose, ka_from_tmax, mg_per_l_to_pg_per_ml,
    simulate, tmax,
)

DAY = 24.0
# Generic levonorgestrel-like drug (30 h half-life) for checks that don't need real values.
SLOW = Drug.one_compartment(name="slow", half_life_h=30.0, ka=ka_from_tmax(1.5, 30.0), F=1.0, V_L=100.0)
# Generic two-compartment drug with well-separated rates.
TWO_COMP = Drug(name="two-comp", ka=2.0, F=0.8, V_L=20.0, k10=0.3, k12=0.6, k21=0.1)

ONE_COMP_DRUGS = [NORETHINDRONE_1C, SLOW]
ALL_DRUGS = [NORETHINDRONE, NORETHINDRONE_1C, SLOW, TWO_COMP]


def ids(d):
    return d.name


# --- The numerical solver ---------------------------------------------------

@pytest.mark.parametrize("drug", ONE_COMP_DRUGS, ids=ids)
def test_single_dose_matches_bateman(drug):
    t = np.linspace(0, 3 * DAY, 2001)
    sim = simulate(drug, [0.0], 1.0, t)
    exact = bateman(drug, 1.0, t)
    np.testing.assert_allclose(sim.conc, exact, rtol=1e-6, atol=1e-9 * exact.max())


@pytest.mark.parametrize("drug", ONE_COMP_DRUGS, ids=ids)
def test_exact_solution_matches_bateman(drug):
    t = np.linspace(0, 3 * DAY, 501)
    np.testing.assert_allclose(exact_single_dose(drug, 1.0, t), bateman(drug, 1.0, t), rtol=1e-9, atol=1e-15)


@pytest.mark.parametrize("drug", [NORETHINDRONE, TWO_COMP], ids=ids)
def test_two_compartment_single_dose_matches_exact(drug):
    t = np.linspace(0, 3 * DAY, 2001)
    sim = simulate(drug, [0.0], 1.0, t)
    exact = exact_single_dose(drug, 1.0, t)
    np.testing.assert_allclose(sim.conc, exact, rtol=1e-6, atol=1e-9 * exact.max())


@pytest.mark.parametrize("drug", [NORETHINDRONE, NORETHINDRONE_1C], ids=ids)
def test_multiple_doses_match_superposition(drug):
    """Late and missed pills: the linear model's answer is a sum of single-dose curves."""
    doses = [0, 24, 51, 96, 120]  # pill 3 is 3 h late, pill 4 (72 h) skipped
    t = np.linspace(0, 7 * DAY, 3001)
    sim = simulate(drug, doses, NORETHINDRONE_DOSE_MG, t)
    exact = sum(exact_single_dose(drug, NORETHINDRONE_DOSE_MG, t - d) for d in doses)
    np.testing.assert_allclose(sim.conc, exact, rtol=1e-6, atol=1e-9 * exact.max())


@pytest.mark.parametrize("drug", ALL_DRUGS, ids=ids)
def test_mass_balance(drug):
    doses = np.arange(0, 28 * DAY, DAY)
    t = np.linspace(0, 30 * DAY, 4001)
    sim = simulate(drug, doses, 0.35, t)
    total = sim.gut + sim.central + sim.peripheral + sim.eliminated + sim.first_pass
    np.testing.assert_allclose(total, sim.dosed, rtol=1e-8, atol=1e-12)
    assert sim.dosed[-1] == pytest.approx(28 * 0.35)


def test_two_pills_at_once_is_a_double_dose():
    t = np.linspace(0, DAY, 201)
    double = simulate(NORETHINDRONE, [0, 0], 0.35, t).conc
    single = simulate(NORETHINDRONE, [0], 0.70, t).conc
    np.testing.assert_allclose(double, single, rtol=1e-10)


@pytest.mark.parametrize(
    "drug, expected", [(NORETHINDRONE_1C, 1.14), (SLOW, 2.35)], ids=["norethindrone-1c", "slow"]
)
def test_steady_state_accumulation(drug, expected):
    """Level just before a dose at steady state vs. 24 h after the first dose."""
    n_days = 40
    t = np.arange(0, n_days + 1) * DAY
    sim = simulate(drug, np.arange(n_days) * DAY, 1.0, t)
    ratio = sim.conc[-1] / sim.conc[1]
    assert ratio == pytest.approx(accumulation_ratio(drug, DAY), rel=1e-5)
    assert ratio == pytest.approx(expected, abs=0.01)  # README values


# --- Two-compartment rates ------------------------------------------------------

@pytest.mark.parametrize("drug", [NORETHINDRONE, TWO_COMP], ids=ids)
def test_disposition_rates(drug):
    alpha, beta = drug.disposition_rates
    assert alpha > beta > 0
    assert alpha + beta == pytest.approx(drug.k10 + drug.k12 + drug.k21)
    assert alpha * beta == pytest.approx(drug.k10 * drug.k21)


@pytest.mark.parametrize("drug", [NORETHINDRONE, TWO_COMP], ids=ids)
def test_terminal_slope_matches_beta(drug):
    t = np.array([60.0, 72.0])
    c = simulate(drug, [0], 1.0, np.linspace(0, 72, 7201)).conc[[6000, 7200]]
    slope = -np.diff(np.log(c))[0] / np.diff(t)[0]
    assert slope == pytest.approx(drug.terminal_rate, rel=1e-4)


def test_one_compartment_rates():
    assert NORETHINDRONE_1C.disposition_rates == (NORETHINDRONE_1C.k10,)
    assert NORETHINDRONE_1C.terminal_half_life_h == pytest.approx(7.7)


def test_ka_from_tmax_reproduces_tmax():
    assert tmax(NORETHINDRONE_1C.ka, NORETHINDRONE_1C.k10) == pytest.approx(1.2, rel=1e-9)
    assert NORETHINDRONE_1C.ka == pytest.approx(3.0, abs=0.1)


def test_invalid_inputs():
    with pytest.raises(ValueError):
        simulate(NORETHINDRONE, [-1], 0.35, [0, 1])
    with pytest.raises(ValueError):
        simulate(NORETHINDRONE, [0], 0.35, [1, 0])
    with pytest.raises(ValueError):
        ka_from_tmax(100, 7.7)
    with pytest.raises(ValueError):
        Drug(name="bad", ka=1, F=1, V_L=1, k10=1, k12=0.5, k21=0)
    with pytest.raises(ValueError):
        Drug(name="bad", ka=1, F=1.5, V_L=1, k10=1)
    with pytest.raises(ValueError):
        bateman(NORETHINDRONE, 1.0, [0, 1])


# --- Norethindrone vs. the label --------------------------------------------------

def _norethindrone_single_dose(drug):
    t = np.linspace(0, DAY, 24001)
    c = mg_per_l_to_pg_per_ml(simulate(drug, [0], NORETHINDRONE_DOSE_MG, t).conc)
    return t, c


def test_norethindrone_matches_label():
    t, c = _norethindrone_single_dose(NORETHINDRONE)
    label = NORETHINDRONE_LABEL
    assert t[np.argmax(c)] == pytest.approx(label["tmax_h"], abs=0.05)
    assert c.max() == pytest.approx(label["cmax_pg_ml"], rel=0.05)
    assert c[-1] == pytest.approx(label["c24_pg_ml"], rel=0.05)
    assert np.trapezoid(c, t) / DAY == pytest.approx(label["cavg_24h_pg_ml"], rel=0.05)
    assert NORETHINDRONE.terminal_half_life_h == pytest.approx(label["half_life_h"], rel=0.05)


@pytest.mark.parametrize("drug, drop", [(NORETHINDRONE, 37), (NORETHINDRONE_1C, 7.6)], ids=["2c", "1c"])
def test_norethindrone_peak_to_24h_drop(drug, drop):
    """The check that motivated the second compartment: label shows a ~37-fold drop."""
    _, c = _norethindrone_single_dose(drug)
    assert c.max() / c[-1] == pytest.approx(drop, rel=0.05)
