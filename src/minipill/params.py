"""Drug parameters, with sources. See sources/Model Parameters .pdf for the full table."""

from minipill.pk import Drug, ka_from_tmax

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
