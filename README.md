# Why Is the Minipill So Strict?

**Modeling hormone levels in progestin-only and combined birth control pills**

02-120 Programming for Scientists · Fall 2026 · Final Project

> **Status:** In progress. The first version of the app goes live by November 20.
>
> **Live app:** _link coming_ · **Video demo:** _link coming_

> This is an educational model built for a course project. It is not medical advice. If you have questions about your own birth control, talk to a doctor or pharmacist.

---

## The question

The progestin-only "minipill" has to be taken within a 3-hour window: CDC guidance counts a minipill taken more than 3 hours late as missed. A single combined pill, on the other hand, can be up to 48 hours late without needing backup contraception.

**Why the difference?** This project tests how much of it comes down to how fast each hormone leaves the body (its half-life), and how much comes from the two hormones in the combined pill affecting each other.

My starting hypothesis:

- **Norethindrone**, the hormone in the minipill, has a half-life of about 8 hours, so very little is left by the time the next pill is due.
- **Levonorgestrel**, the progestin in the combined pill, lasts about 30 to 36 hours.
- The **estrogen** in the combined pill (ethinyl estradiol) makes the liver produce more SHBG, a blood protein that binds levonorgestrel and slows its clearance. So levonorgestrel builds up over the month in a way norethindrone never does.

## What the program does

The app simulates hormone levels in the blood over a full pill pack, with both pills side by side. You can:

- choose a pill: the minipill (0.35 mg norethindrone) or a combined pill (levonorgestrel + ethinyl estradiol, 21 active / 7 placebo);
- take each dose on time, late, or not at all;
- watch how far hormone levels fall and how long they take to recover.

There is no single "safe" hormone level, because the minipill only stops ovulation in about half of users and works largely by thickening cervical mucus. So instead of using a made-up threshold, the app compares each pill against itself: **how far below that pill's usual daily low point do levels fall when a dose is 3 hours late versus 48 hours late?** If the drops are similar, half-life explains the timing rules. If not, something else must matter too.

## How the model works

Each hormone uses a compartment model with oral absorption. A pill puts a dose into the gut, the hormone moves from the gut into the blood, and the body clears it from the blood:

$$\frac{dG}{dt} = -k_a G$$

$$\frac{dC}{dt} = \frac{F \, k_a \, G}{V} - k \, C, \qquad k = \frac{\ln 2}{t_{1/2}}$$

| Symbol | Meaning |
| --- | --- |
| G | amount of hormone still in the gut |
| C | concentration of hormone in the blood |
| kₐ | absorption rate (chosen to match each drug's time to peak) |
| F | fraction of the dose that reaches the blood (bioavailability) |
| V | volume of distribution |
| k | elimination rate, set by the half-life |

Each pill taken adds one dose to G at the time it's taken.

**Planned additions:**

- **A second compartment** for norethindrone, if the one-compartment model can't match its measured 24-hour level (see checks below).
- **The estrogen–SHBG effect:** ethinyl estradiol slowly raises SHBG over the cycle, and higher SHBG slows levonorgestrel clearance. The strength and speed of this effect aren't reported directly, so I calibrate them to match published data and say so clearly.

## How I know it's right

| Check | What it tests | Expected result | Status |
| --- | --- | --- | --- |
| Exact single-dose solution | The numerical simulation | Matches the Bateman function | ☐ Planned |
| Mass balance | No drug appears or disappears | Absorbed + eliminated + remaining = total dose | ☐ Planned |
| Steady-state accumulation | Repeated daily dosing | 1 / (1 − e^(−kτ)): ≈ 1.14 for norethindrone, ≈ 2.4 for levonorgestrel | ☐ Planned |
| Norethindrone single dose | Published label values | Peak ≈ 4,800 pg/mL near 1.2 h; ≈ 130 pg/mL at 24 h | ☐ Planned |
| Combined pill steady state | Published label values (0.10 mg levonorgestrel) | Steady state from day 19; peak ≈ 6.0 ng/mL, low ≈ 1.9 ng/mL | ☐ Planned |
| SHBG effect | A behavior the science requires | Without SHBG, the model underpredicts levonorgestrel's 3- to 4-fold rise; with it, the gap closes | ☐ Planned |

I expect the one-compartment model to **fail** the norethindrone 24-hour check: the measured drop is about 37-fold in 24 hours, much steeper than the 8-fold an 8-hour half-life predicts. Finding out why, and whether a second compartment fixes it, is part of the project.

## Running it

**Online:** Use the live app (link above). Nothing to install.

**Locally:** _Instructions coming once the first version is built._

**Running the tests:** _Instructions coming._

## Repository layout

```
sources/    Drug labels and papers the parameters come from
src/        The model and simulation code
tests/      Checks against exact solutions and published values
app/        The web app
journal/    Claude Code session logs (my build journal)
```

_This layout will be updated as the project grows._

## Key parameters

Full tables, with sources for every value, are in [`sources/`](sources/).

| Hormone | Half-life | Notes |
| --- | --- | --- |
| Norethindrone (minipill) | 7.7 ± 0.5 h | Peak ≈ 1.2 h after a dose; taken daily with no break |
| Levonorgestrel (combined) | ≈ 30 h single dose; 36 ± 13 h at steady state | Binds SHBG; rises 3- to 4-fold over the cycle |
| Ethinyl estradiol (combined) | ≈ 15 h single dose; 18 ± 4.7 h at steady state | Raises SHBG ≈ 1.7-fold; broken down partly by CYP3A4 |

## Assumptions and limitations

- Hormone level is used as a stand-in for protection. In reality, the minipill also works through cervical mucus, and the model doesn't simulate ovulation itself.
- Parameters are averages from label studies. Real people vary, and published data show that body weight can roughly double levonorgestrel's half-life.
- The SHBG effect is calibrated to match published data, not measured directly.
- Drug interactions, food, and vomiting or diarrhea are not modeled (yet).

## Possible extensions

- **Body weight:** reproduce published differences in levonorgestrel kinetics between normal-BMI and obese users (Edelman et al., 2009), then see what the model predicts for the minipill, where no such data exist.
- **Drug interactions:** a drug like rifampin speeds up the enzyme that breaks down ethinyl estradiol, which should lower SHBG and then levonorgestrel.

## Sources

- Teorell, T. (1937). Kinetics of distribution of substances administered to the body. *Archives Internationales de Pharmacodynamie et de Thérapie*, 57.
- [Norethindrone 0.35 mg prescribing information](https://medlibrary.org/lib/rx/meds/norethindrone/)
- [Levonorgestrel/ethinyl estradiol prescribing information (Alesse)](https://www.drugs.com/pro/alesse.html)
- [Levonorgestrel/ethinyl estradiol prescribing information (Seasonale)](https://www.drugs.com/pro/seasonale.html)
- [Levonorgestrel/ethinyl estradiol package insert (Levette)](https://www.sahpra.org.za/wp-content/uploads/2020/08/Levette_PI_Actor-Pharma_MCC-FormatJuly-2015.pdf)
- [CDC recommended actions after late or missed pills (UCSF summary)](https://beyondthepill.ucsf.edu/wp-content/uploads/2025/01/recommended-action-late-missed-progestin-508.pdf)
- Edelman, A. B., et al. (2009). [Impact of obesity on oral contraceptive pharmacokinetics and hypothalamic-pituitary-ovarian activity](https://pmc.ncbi.nlm.nih.gov/articles/PMC2736633). *Contraception*, 80(2), 119–127.

## About

Built by _[your name]_ for 02-120 Programming for Scientists at Carnegie Mellon University, Fall 2026, with Claude Code as a collaborator.
