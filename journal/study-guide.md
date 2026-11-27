# Minipill Model: Key Functions Study Guide

Nov 26, 2026 · Exported from the [Claude Docs version](https://claude.ai/artifact/Uufkm74pbhni1AeNibZkT2), which has the drawn diagram.

## The big picture

The model is `pk.py`; everything else feeds it numbers, checks it, or shows its results. Read the files in this order: `params.py`, then `Drug` and `simulate` in `pk.py`, then one test, then `fit_norethindrone.py`.

```
Drug labels, studies  -->  scripts/fit_*.py  -->  params.py
(peak, AUC, half-life,     (least_squares finds    (the drugs, SHBG, BMI curve,
 SHBG; Alesse, Edelman)     values that match)      each value with its source)
                                                         |  drug values
                                                         v
                           tests/test_pk.py  -->  pk.py: THE MODEL
                           (57 checks)             (Drug, simulate, SHBGEffect,
                                                    exact one-pill formulas)
                                                         |  runs simulate
                                                         v
You, in a browser  <--  pill_comparison.html  <--  export_comparison_data.py
(weight, height,         (blend and updateBody      (runs 9 BMIs, writes MODEL
 late pill)               draw the charts)           into the web page)
```

The fitting scripts also call `pk.py` (`exact_single_dose`, `simulate`) to predict the curves they compare with the labels.

## 1. Describing a drug: `Drug` (pk.py)

`Drug` holds one hormone's numbers; everything else in the model reads them. It is a frozen dataclass: a short class for named values that can't be changed after it's made, so parameters can't be overwritten by accident.

| Field | Meaning | Unit |
| --- | --- | --- |
| `ka` | absorption, gut to blood | 1/h |
| `F` | fraction that survives first pass (not a rate) | 0–1 |
| `V_L` | blood (central) volume | L |
| `k10` | elimination from the blood | 1/h |
| `k12`, `k21` | blood to tissues, tissues to blood (both 0 = one compartment) | 1/h |

What comes with it:

- `__post_init__` checks the inputs as soon as a drug is made (e.g. 0 < F ≤ 1).
- `@property` values are calculated when you read them: `terminal_half_life_h`, `clearance_L_per_h` (= `k10 × V_L`), `vss_L`, `disposition_rates` (α and β).
- `Drug.one_compartment(...)` builds a one-compartment drug from its half-life.
- `with_clearance_scaled(factor)` (new) returns a copy with only `k10` multiplied; see section 5.

The real drugs live in `params.py`: `NORETHINDRONE` and `LEVONORGESTREL` (two compartments, fitted to their labels) and `*_1C` versions kept for comparison.

## 2. Simulating blood levels: `simulate` (pk.py)

`simulate(drug, dose_times, dose_mg, t_eval, shbg=None)` is the main function: give it pill times, get back the amount of drug in each compartment at each time, plus the concentration. Late, missed or double pills are just a different `dose_times` list.

**The state it tracks** (all in mg): gut, central (blood), peripheral (tissues), eliminated, first pass, and SHBG (nmol/L) when the SHBG effect is on. Gut + blood + tissues + eliminated + first pass always equals the amount swallowed, which is the mass-balance check.

**The equations:** `_rhs` (or `_rhs_shbg`) returns how fast each amount is changing right now. It does not give amounts directly; SciPy's `solve_ivp` (LSODA) steps forward using those rates.

**How the main loop works:**

1. **Breakpoints:** sort every pill time (and, with SHBG, every "estrogen off" time 24 h after a pill) into a list. Pills are instant jumps, which a solver can't handle mid-step.
2. **For each segment between breakpoints:** add any pill taken at its start to the gut (`y[0] += n * dose_mg`; `Counter` makes two pills at once a double dose).
3. **Solve** the smooth stretch with `solve_ivp`, and copy out the requested times that fall in it (`mask`, `searchsorted`).
4. **Carry the state forward:** `y = sol.y[:, -1].copy()`. The end of one segment is the start of the next, so drug left from earlier pills is still there when the next pill lands. This is how levels build up over a pack.

**After the loop:** concentration is calculated once, `conc = central / V_L` (mg/L), and everything is returned in a `SimResult` (`result.conc`, `result.central`, …).

`args = (ka, F, k10, k12, k21)`: the drug's values passed to `_rhs`, in the same order as its parameters. `V_L` isn't in it because the equations work in amounts.

## 3. Exact formulas for one pill (pk.py)

These give the exact curve after a single pill, with no step-by-step solving. They check `simulate`, and they are fast enough for fitting.

**`bateman(drug, dose_mg, t)`**: one compartment only.

$$C(t) = \frac{F \, D \, k_a}{V (k_a - k)} \left( e^{-k t} - e^{-k_a t} \right)$$

- e^(−kt) falls slowly: elimination, the long tail that sets the half-life.
- e^(−ka·t) falls fast: absorption finishing.
- At t = 0 both are 1, so C = 0 (nothing absorbed yet). Needs ka ≠ k, or it divides by zero.
- Setting the slope to 0 gives `tmax = ln(ka/k) / (ka − k)`; `ka_from_tmax` runs it backwards.

**`exact_single_dose(drug, dose_mg, t)`**: one or two compartments. Writes the model as y′ = M·y (gut, blood, tissues), finds M's eigenvalues (−ka, −α, −β) with `np.linalg.eig`, and gets a sum of three exponentials:

$$C(t) = a_1 e^{-k_a t} + a_2 e^{-\alpha t} + a_3 e^{-\beta t}$$

**`single_dose_terms(drug, dose_mg)`** (new): the first half of `exact_single_dose`, split out. Returns the rates (ka, α, β) and amplitudes (a) on their own; `exact_single_dose` now just adds them up. The export script uses it to send each drug's formula to the web page.

**Why the minipill can be summed but the combined pill can't:** without SHBG the model is linear, so a month of pills is the single-pill curve added once per pill, shifted to each pill's time (superposition). With SHBG, clearance changes over time, so pills can't be added up and `simulate` must run every scenario.

## 4. The estrogen effect: `SHBGEffect` (pk.py)

The combined pill's estrogen makes the liver produce more SHBG, and more SHBG slows levonorgestrel's clearance. `SHBGEffect` holds four numbers, fitted to the Alesse label by `scripts/fit_shbg.py`:

| Field | Value | Meaning |
| --- | --- | --- |
| `s0` | 57 nmol/L | baseline SHBG (Alesse day 1) |
| `s_max` | 93.5 nmol/L | level SHBG approaches with daily estrogen |
| `k_per_h` | 0.0089 /h | how fast SHBG moves (half-way in 3.2 days) |
| `gamma` | 1.331 | how strongly SHBG slows clearance |

**`clearance_factor(S)`** returns (s0 / S)^gamma, and `_rhs_shbg` multiplies `k10` by it at every moment:

$$k_{10}(t) = k_{10} \left( \frac{s_0}{S} \right)^{\gamma}$$

- **gamma = 0:** anything to the power 0 is 1, so SHBG has no effect and levonorgestrel doesn't build up.
- **gamma = 1.33:** a 1% rise in SHBG cuts clearance about 1.33%. At the plateau clearance is about 52% of normal, matching the label's CL/F falling 53.7 → 28.4.
- Above 1 because SHBG stands for the whole estrogen effect, not only binding (the label's unbound clearance falls too).

SHBG moves toward `s_max` while estrogen is on (24 h after each active pill) and back toward `s0` in the placebo week.

## 5. Body size: `clearance_factor_for_bmi` and `with_clearance_scaled`

Body size changes only clearance. Edelman 2009 measured levonorgestrel's half-life at 25.6 h (normal BMI, 21.9) vs. 52.1 h (obese, 37.3) with the same volume; same volume and twice the half-life gives clearance × 25.6 / 52.1 = **0.491**.

**`clearance_factor_for_bmi(bmi)`** (params.py) turns a BMI into a clearance multiplier:

```python
frac = np.clip((np.asarray(bmi, dtype=float) - BMI_REFERENCE) / (BMI_OBESE - BMI_REFERENCE), 0.0, 1.0)
return OBESE_CLEARANCE_RATIO ** frac
```

- `np.asarray(..., dtype=float)` lets it take one number, a list or an array.
- `frac` = **where** the BMI sits between the groups: 0 at 21.9, 0.5 at 29.6, 1 at 37.3.
- `np.clip` keeps frac in 0–1, so the curve is flat outside the study (×1 below 21.9, ×0.491 above 37.3).
- `OBESE_CLEARANCE_RATIO` (0.491) = **how big** the full effect is. `**` combines them; each BMI point multiplies clearance by 0.955.

| BMI | 18 | 21.9 | 25 | 29.6 | 33 | 37.3 | 45 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| multiplier | 1 | 1 | 0.87 | 0.70 | 0.60 | 0.49 | 0.49 |

**`Drug.with_clearance_scaled(factor)`** (pk.py) applies it. `Drug` is frozen, so `dataclasses.replace(self, k10=self.k10 * factor)` makes a copy with only `k10` changed. Clearance (= `k10 × V_L`) scales exactly; absorption and volumes stay put.

```python
f = clearance_factor_for_bmi(37.3)                 # 0.491
obese = LEVONORGESTREL.with_clearance_scaled(f)    # clearance 2.84 -> 1.40 L/h
```

Two things to know: the model's half-life rises ×1.75 (34 → 60 h), short of the study's ×2.03; and norethindrone gets the same change with no data behind it.

## 6. Fitting, exporting and the website

**Fitting scripts** (`scripts/fit_*.py`) find the parameter values in `params.py`. Most weren't measured; the labels only give summary numbers (peak, time to peak, 24 h level, AUC).

1. Guess parameters; build a `Drug`.
2. Predict the single-pill curve with `exact_single_dose` (fast, no solver) and read off the same summary numbers.
3. `residuals` = how far each prediction misses the label, as log(model / label).
4. SciPy's `least_squares` adjusts the guess until the misses are as small as possible.

Tricks: fit in log space (`np.exp(log_params)`) so values stay positive; start from 40 random guesses and keep the best. `fit_shbg.py` uses `simulate` instead, because its data (days 6 and 21) involve many pills and SHBG.

**`export_comparison_data.py`** runs the model and writes the results into the web page:

- `main()` loops over 9 BMIs from 21.9 to 37.3; at each it scales both drugs with `with_clearance_scaled`.
- `terms()` saves each drug's single-pill formula (from `single_dose_terms`, in ng/mL).
- `combined_pill()` runs `simulate` with SHBG for the month and all 50 late-pill options.
- Everything goes into one JSON object, `MODEL`, written between the `/*DATA*/` markers in `app/pill_comparison.html`. About 70 s.

**In the page (JavaScript):**

- `single(m, t)` sums a·e^(−r·t); `multi(m, doses, t)` adds one single-pill curve per pill (superposition).
- `blend(f)` mixes the two grid BMIs nearest the one entered: (1 − gw)·f(gi) + gw·f(gi + 1). Every BMI-dependent curve goes through it.
- `updateBody()` runs on every input change: works out BMI, finds `gi` and `gw`, fills the BMI/clearance/half-life tiles, and redraws the charts. `setUnits()` switches kg/cm and lb/ft-in without changing the BMI.

## 7. Which function is used where

| Function | File | Tests | Fitting | Website |
| --- | --- | --- | --- | --- |
| `simulate` | pk.py | yes | SHBG fit | combined-pill curves (precomputed) |
| `exact_single_dose` | pk.py | yes | norethindrone, levonorgestrel fits | no |
| `single_dose_terms` | pk.py | yes | no | minipill and no-SHBG curves (via `terms`) |
| `bateman` | pk.py | **only** here | no | no |
| `SHBGEffect` | pk.py | yes | SHBG fit | combined-pill curves |
| `clearance_factor_for_bmi` | params.py | yes | no | sets the 9-BMI grid |
| `with_clearance_scaled` | pk.py | yes | no | applies each grid BMI |
| `clearance_L_per_h` | pk.py | yes | levonorgestrel AUC; printed | no |
| `blend`, `updateBody` | pill_comparison.html | headless Chrome check | no | the BMI input |

## 8. Test yourself

Cover the answers, then check. The starred ones are from past quizzes where I slipped.

| Question | Answer |
| --- | --- |
| What does `y = sol.y[:, -1].copy()` do in `simulate`, and why does it matter? ★ | Takes the state at the end of a segment as the start of the next. Drug left from earlier pills is still there when the next pill lands, so levels build up and a late pill falls further. |
| Why does `simulate` split time at pills? | A pill is an instant jump in the gut amount; the solver assumes smooth change. Solve each smooth stretch, apply the jump between. |
| Why is `V_L` not in `args`? | The equations work in amounts (mg). Volume is used once at the end: `conc = central / V_L`. |
| In Bateman's formula, which term makes the slow fall? ★ | e^(−kt), elimination. e^(−ka·t) is the fast one (absorption finishing). |
| What happens with gamma = 0? ★ | (s0/S)^0 = 1, so SHBG has no effect on clearance; day-21 AUC stays about 35 instead of 68. |
| Which function makes the website's curves, and which fits norethindrone? ★ | `simulate` (with SHBG) and `single_dose_terms` make the curves; `exact_single_dose` fits norethindrone; `bateman` is only in tests. |
| What does `clearance_factor_for_bmi(60)` return? ★ | 0.491. frac = 2.47 is clipped down to 1. ×1 only happens at BMI 21.9 or below. |
| A new study says 0.6× instead of 0.49×. What do you change? ★ | `OBESE_CLEARANCE_RATIO` (or the half-lives in its source dictionary). `frac` only says where a BMI sits. |
| If `BMI_OBESE` moved to 40, what happens at BMI 37.3? ★ | Bigger: frac = 0.851 (not clipped), 0.491^0.851 = 0.546. A smaller power of a number below 1 is bigger. |
| Why does a higher BMI make late pills dip less in the model? | Slower clearance means a longer half-life, so less drug leaves in the hours a pill is late. |

**To revisit:** carrying the state forward in `simulate`; which way `np.clip` pushes values outside the range; that a smaller power of a number below 1 is bigger; which function is used for which fit.
