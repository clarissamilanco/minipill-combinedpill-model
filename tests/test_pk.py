import numpy as np
import pytest

from minipill.params import (
    ALESSE_DAY1, ALESSE_SHBG_NMOL_L, ALESSE_STEADY_STATE, ALESSE_TABLE_I, LEVONORGESTREL, LEVONORGESTREL_1C, SHBG,
    LEVONORGESTREL_DOSE_MG, NORETHINDRONE, NORETHINDRONE_1C, NORETHINDRONE_DOSE_MG, NORETHINDRONE_LABEL,
)
from minipill.pk import (
    Drug, SHBGEffect, accumulation_ratio, bateman, exact_single_dose, ka_from_tmax, mg_per_l_to_ng_per_ml,
    mg_per_l_to_pg_per_ml, simulate, tmax,
)

DAY = 24.0
# Generic levonorgestrel-like drug (30 h half-life) for checks that don't need real values.
SLOW = Drug.one_compartment(name="slow", half_life_h=30.0, ka=ka_from_tmax(1.5, 30.0), F=1.0, V_L=100.0)
# Generic two-compartment drug with well-separated rates.
TWO_COMP = Drug(name="two-comp", ka=2.0, F=0.8, V_L=20.0, k10=0.3, k12=0.6, k21=0.1)

ONE_COMP_DRUGS = [NORETHINDRONE_1C, LEVONORGESTREL_1C, SLOW]
TWO_COMP_DRUGS = [NORETHINDRONE, LEVONORGESTREL, TWO_COMP]
ALL_DRUGS = ONE_COMP_DRUGS + TWO_COMP_DRUGS


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


@pytest.mark.parametrize("drug", TWO_COMP_DRUGS, ids=ids)
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

@pytest.mark.parametrize("drug", TWO_COMP_DRUGS, ids=ids)
def test_disposition_rates(drug):
    alpha, beta = drug.disposition_rates
    assert alpha > beta > 0
    assert alpha + beta == pytest.approx(drug.k10 + drug.k12 + drug.k21)
    assert alpha * beta == pytest.approx(drug.k10 * drug.k21)


@pytest.mark.parametrize("drug", TWO_COMP_DRUGS, ids=ids)
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


# --- Levonorgestrel vs. Alesse (no SHBG effect yet) ------------------------------

def _levonorgestrel_single_dose(drug, hours=DAY):
    t = np.linspace(0, hours, int(hours * 1000) + 1)
    c = mg_per_l_to_ng_per_ml(simulate(drug, [0], LEVONORGESTREL_DOSE_MG, t).conc)
    return t, c


def test_levonorgestrel_matches_alesse_day1():
    t, c = _levonorgestrel_single_dose(LEVONORGESTREL)
    label = ALESSE_DAY1
    assert t[np.argmax(c)] == pytest.approx(label["tmax_h"], abs=0.05)
    assert c.max() == pytest.approx(label["cmax_ng_ml"], rel=0.02)
    auc = mg_per_l_to_ng_per_ml(LEVONORGESTREL.F * LEVONORGESTREL_DOSE_MG / LEVONORGESTREL.clearance_L_per_h)
    assert auc == pytest.approx(label["auc_ng_h_ml"], rel=0.02)
    assert LEVONORGESTREL.terminal_half_life_h == pytest.approx(label["half_life_h"], rel=0.02)


def test_levonorgestrel_one_compartment_underpredicts_peak():
    """Same failure as norethindrone: the label's CL/F and Vz/F give a peak ~4x too low."""
    t, c = _levonorgestrel_single_dose(LEVONORGESTREL_1C)
    assert t[np.argmax(c)] == pytest.approx(ALESSE_DAY1["tmax_h"], abs=0.01)
    assert c.max() == pytest.approx(0.69, abs=0.01)
    assert ALESSE_DAY1["cmax_ng_ml"] / c.max() == pytest.approx(4.0, abs=0.1)


def test_without_shbg_day21_auc_is_about_half_of_alesse():
    """Behavior the science requires: without SHBG, a linear model's day-21 AUC over one
    dosing interval equals the day-1 AUC (35.2), but Alesse measured 68.3 ng*h/mL.
    The missing factor (~1.9) is the SHBG effect, matching the CL/F fall 53.7 -> 28.4."""
    t = np.linspace(20 * DAY, 21 * DAY, 2401)
    sim = simulate(LEVONORGESTREL, np.arange(21) * DAY, LEVONORGESTREL_DOSE_MG, np.concatenate([[0.0], t]))
    c = mg_per_l_to_ng_per_ml(sim.conc[1:])
    auc_day21 = np.trapezoid(c, t)
    assert auc_day21 == pytest.approx(ALESSE_DAY1["auc_ng_h_ml"], rel=0.02)
    label_auc21, label_cl1, label_cl21 = ALESSE_TABLE_I[21][0], ALESSE_TABLE_I[1][1], ALESSE_TABLE_I[21][1]
    assert label_auc21 / auc_day21 == pytest.approx(label_cl1 / label_cl21, rel=0.05)
    # Trough before the next pill: ~0.84 ng/mL without SHBG vs. 1.9 ng/mL measured.
    assert c[-1] == pytest.approx(0.84, abs=0.02)
    assert c[-1] < 0.5 * ALESSE_STEADY_STATE["cmin_ng_ml"]


# --- SHBG effect (combined pill) -------------------------------------------------

PACK_T = np.linspace(0, 28 * DAY, 28 * 120 + 1)  # every 12 min over one pack


@pytest.fixture(scope="module")
def pack_with_shbg():
    return simulate(LEVONORGESTREL, np.arange(21) * DAY, LEVONORGESTREL_DOSE_MG, PACK_T, shbg=SHBG)


def _at(t):
    return int(np.searchsorted(PACK_T, t))


def _day(sim, d):
    i, j = _at((d - 1) * DAY), _at(d * DAY)
    c = mg_per_l_to_ng_per_ml(sim.conc[i:j + 1])
    return np.trapezoid(c, PACK_T[i:j + 1]), c.max(), c[-1]


def test_shbg_with_zero_strength_changes_nothing():
    t = np.linspace(0, 5 * DAY, 1201)
    doses = np.arange(5) * DAY
    plain = simulate(LEVONORGESTREL, doses, LEVONORGESTREL_DOSE_MG, t).conc
    off = SHBGEffect(SHBG.s0, SHBG.s_max, SHBG.k_per_h, gamma=0.0)
    np.testing.assert_allclose(simulate(LEVONORGESTREL, doses, LEVONORGESTREL_DOSE_MG, t, shbg=off).conc,
                               plain, rtol=1e-7, atol=1e-12)


def test_shbg_follows_alesse(pack_with_shbg):
    S = pack_with_shbg.shbg
    for day, label in ALESSE_SHBG_NMOL_L.items():
        assert S[_at((day - 1) * DAY)] == pytest.approx(label, abs=0.5)


def test_shbg_falls_back_in_placebo_week(pack_with_shbg):
    S = pack_with_shbg.shbg
    assert SHBG.s0 < S[-1] < S[_at(21 * DAY)]


def test_mass_balance_with_shbg(pack_with_shbg):
    sim = pack_with_shbg
    total = sim.gut + sim.central + sim.peripheral + sim.eliminated + sim.first_pass
    np.testing.assert_allclose(total, sim.dosed, rtol=1e-8, atol=1e-12)


def test_shbg_matches_alesse_auc(pack_with_shbg):
    """Fitted targets: AUC over one dosing interval on days 6 and 21."""
    for day in (6, 21):
        auc, _, _ = _day(pack_with_shbg, day)
        assert auc == pytest.approx(ALESSE_TABLE_I[day][0], rel=0.05)


def test_shbg_predicts_steady_state_low_and_peak(pack_with_shbg):
    """Not fitted: the day-21 low and peak land within the label's mean ± SD."""
    _, peak, low = _day(pack_with_shbg, 21)
    assert low == pytest.approx(ALESSE_STEADY_STATE["cmin_ng_ml"], rel=0.10)  # 2.04 vs. 1.9
    assert abs(peak - ALESSE_STEADY_STATE["cmax_ng_ml"]) < 2.65  # 4.94 vs. 6.0 ± 2.65
    assert peak < ALESSE_STEADY_STATE["cmax_ng_ml"]  # known shortfall: volume kept fixed


def _late_drop(drug, doses, dose_mg, late_h, shbg=None):
    """% below the usual low when the day-15 pill is ``late_h`` hours late (None = skipped)."""
    due = 14 * DAY
    t = np.linspace(0, due + 5 * DAY, int((due + 5 * DAY) * 20) + 1)
    on = simulate(drug, doses, dose_mg, t, shbg=shbg).conc
    moved = [d for d in doses if d != due] + ([] if late_h is None else [due + late_h])
    lat = simulate(drug, sorted(moved), dose_mg, t, shbg=shbg).conc
    low = on[int(np.searchsorted(t, due))]
    after = t >= due
    return 100 * (1 - lat[after].min() / low)


def test_late_pill_hurts_minipill_more_than_combined_pill():
    """The project's question, in the model: same delay, much smaller dip for the combined pill."""
    net = _late_drop(NORETHINDRONE, np.arange(28) * DAY, NORETHINDRONE_DOSE_MG, 3)
    lng = _late_drop(LEVONORGESTREL, np.arange(21) * DAY, LEVONORGESTREL_DOSE_MG, 3, shbg=SHBG)
    assert net == pytest.approx(23, abs=1)
    assert lng < net / 3
