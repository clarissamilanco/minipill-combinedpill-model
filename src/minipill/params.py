"""Drug parameters, with sources. See sources/Model Parameters .pdf for the full table."""

from minipill.pk import Drug, ka_from_tmax

# Assumed body weight to turn the label's L/kg volume into litres. Not from a source.
BODY_WEIGHT_KG = 60.0

# --- Norethindrone (minipill) ------------------------------------------------
NORETHINDRONE_DOSE_MG = 0.35  # taken every 24 h, no break

NORETHINDRONE = Drug(
    name="norethindrone",
    half_life_h=7.7,  # norethindrone label, single dose (7.7 ± 0.5 h)
    ka=ka_from_tmax(1.2, 7.7),  # ESTIMATED: chosen so Tmax = 1.2 h (label); ≈ 3.0 /h
    F=0.65,  # norethindrone label, absolute bioavailability ≈ 65%
    # Norethindrone label: ≈ 4 L/kg. This is the whole-body volume, so a one-compartment
    # model using it underpredicts the label's peak (4,817 pg/mL) about 5-fold.
    # Provisional until the second compartment is added.
    V_L=4.0 * BODY_WEIGHT_KG,
)

# Label values to check the model against (norethindrone label, single 0.35 mg dose).
NORETHINDRONE_LABEL = {
    "tmax_h": 1.2,
    "cmax_pg_ml": 4817,
    "cavg_24h_pg_ml": 885,
    "c24_pg_ml": 130,
}
