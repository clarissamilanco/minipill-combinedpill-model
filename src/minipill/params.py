"""Drug parameters, with sources. See sources/Model Parameters .pdf for the full table."""

from minipill.pk import LN2, Drug, SHBGEffect, ka_from_tmax

# Assumed body weight to turn the label's L/kg volume into litres. Not from a source.
BODY_WEIGHT_KG = 60.0

# --- Norethindrone (minipill) ------------------------------------------------
NORETHINDRONE_DOSE_MG = 0.35  # taken every 24 h, no break

# Label values to check the model against (norethindrone label, single 0.35 mg dose).
NORETHINDRONE_LABEL = {
    "tmax_h": 1.2,  # ± 0.5
    "cmax_pg_ml": 4817,  # ± 1533
    "cavg_24h_pg_ml": 885,  # ± 250
    "c24_pg_ml": 130,  # ± 47
    "half_life_h": 7.7,  # ± 0.5
}

# Two-compartment model (the one to use).
# F is from the label (absolute bioavailability ≈ 65%). The other five values are
# FITTED by scripts/fit_norethindrone.py to the five label values above; they are
# not measured. The fit reproduces every label value within 3%.
#
# Notes on the fit:
#  - Any starting point gives the same fit, always with ka ≈ alpha (the fast
#    distribution rate, ≈ 0.91 /h). The label's summary values can't tell absorption
#    and distribution apart, so ka and alpha individually are not well determined;
#    the curve shape, 24 h level and terminal phase are.
#  - Implied clearance ≈ 234 L/day and Vss ≈ 57 L (≈ 1 L/kg at 60 kg), lower than the
#    label's ~600 L/day and ~4 L/kg. Those label figures come from different studies
#    and don't match the label's own concentration data.
NORETHINDRONE = Drug(
    name="norethindrone",
    ka=0.911,  # 1/h, fitted
    F=0.65,  # norethindrone label
    V_L=18.7,  # L, central volume, fitted
    k10=0.520,  # 1/h, fitted
    k12=0.320,  # 1/h, fitted
    k21=0.155,  # 1/h, fitted
)

# One-compartment model, kept for comparison: it gives a ~8-fold drop from peak to
# 24 h instead of the label's ~37-fold, and with the label's 4 L/kg volume it
# underpredicts the peak about 5-fold.
NORETHINDRONE_1C = Drug.one_compartment(
    name="norethindrone (1 compartment)",
    half_life_h=7.7,  # norethindrone label
    ka=ka_from_tmax(1.2, 7.7),  # ESTIMATED: chosen so Tmax = 1.2 h (label); ≈ 3.0 /h
    F=0.65,  # norethindrone label
    V_L=4.0 * BODY_WEIGHT_KG,  # norethindrone label, ≈ 4 L/kg
)

# --- Levonorgestrel (combined pill, Alesse 0.10 mg) -------------------------
LEVONORGESTREL_DOSE_MG = 0.10  # 21 active days, then 7 placebo days

# Alesse label, Table I (mean of 22 women). Day 1 is a single dose, before
# ethinyl estradiol has raised SHBG; days 6 and 21 include the SHBG effect.
ALESSE_DAY1 = {
    "cmax_ng_ml": 2.75,  # ± 0.88
    "tmax_h": 1.6,  # ± 0.9
    "auc_ng_h_ml": 35.2,  # ± 12.8 (the AUC behind CL/F, so AUC to infinity)
    # ln 2 · (Vz/F) / (CL/F) = ln 2 · 2.66 L/kg / 0.0537 L/h/kg; body weight cancels.
    "half_life_h": 34.3,
}
ALESSE_TABLE_I = {  # day: (AUC ng·h/mL, CL/F mL/h/kg, SHBG nmol/L)
    1: (35.2, 53.7, 57),
    6: (46.0, 40.8, 81),
    21: (68.3, 28.4, 93),  # day 21 AUC is over one 24 h dosing interval
}
ALESSE_STEADY_STATE = {"cmax_ng_ml": 6.0, "cmin_ng_ml": 1.9}  # from day 19 on
ALESSE_SHBG_NMOL_L = {1: 57, 6: 81, 21: 93}  # ± 18, 25, 40

# Two-compartment model. Used alone, it has NO SHBG effect; pass SHBG below to
# simulate() for the combined pill over a pack.
# F is from the label (≈ 100%). ka is ASSUMED: the day-1 data can't fix it; any
# ka from ~0.4 to 1.2 /h fits exactly and gives a 24 h level within ±5% and the
# same late-pill drop (run scripts/fit_levonorgestrel.py --profile).
# The other four are FITTED by scripts/fit_levonorgestrel.py to Alesse day 1
# (peak, time to peak, AUC, terminal half-life); not measured.
LEVONORGESTREL = Drug(
    name="levonorgestrel",
    ka=0.8,  # 1/h, assumed
    F=1.0,  # Alesse label
    V_L=16.7,  # L, central volume, fitted
    k10=0.170,  # 1/h, fitted
    k12=0.339,  # 1/h, fitted
    k21=0.066,  # 1/h, fitted
)

# One-compartment model from the label's CL/F and Vz/F, kept for comparison:
# like norethindrone, it underpredicts the day-1 peak about 4-fold.
_LNG_HALF_LIFE_H = ALESSE_DAY1["half_life_h"]
_LNG_CL_L_PER_H = LEVONORGESTREL_DOSE_MG * 1e3 / ALESSE_DAY1["auc_ng_h_ml"]  # 2.84 L/h
LEVONORGESTREL_1C = Drug.one_compartment(
    name="levonorgestrel (1 compartment)",
    half_life_h=_LNG_HALF_LIFE_H,
    ka=ka_from_tmax(ALESSE_DAY1["tmax_h"], _LNG_HALF_LIFE_H),  # ESTIMATED so Tmax = 1.6 h
    F=1.0,
    V_L=_LNG_CL_L_PER_H * _LNG_HALF_LIFE_H / LN2,  # = Vz ≈ 141 L
)

# --- SHBG effect (combined pill) ----------------------------------------------
# The pill's estrogen makes the liver produce more SHBG, a blood protein that
# binds levonorgestrel; with more SHBG, levonorgestrel is cleared more slowly.
# FITTED by scripts/fit_shbg.py to Alesse Table I:
#  - SHBG curve (s0, s_max, rate) passes exactly through 57, 81, 93 nmol/L on days 1, 6, 21.
#  - gamma fitted to the day-6 and day-21 AUCs (46.0, 68.3 ng*h/mL): model 48.0 and 66.8.
# Checks not used in the fit: low before the next pill on day 21 is 2.04 vs. 1.9
# ng/mL measured; the day-21 peak is 4.94 vs. 6.0 ± 2.65 (low, since the model
# keeps the volume fixed while the label's Vz/F falls from 2.66 to 1.43 L/kg).
# Caveat: the label's unbound data show unbound clearance also falls (2.79 -> 1.57
# L/h/kg), so SHBG here stands for the whole estrogen effect on clearance.
SHBG = SHBGEffect(
    s0=57.0,  # nmol/L, Alesse day 1
    s_max=93.5,  # nmol/L, fitted plateau
    k_per_h=0.2143 / 24,  # 1/h, fitted; SHBG gets half-way to its plateau in 3.2 days
    gamma=1.331,  # fitted
)
