# The supported-duration endpoint: theory, arithmetic and limits

2026-10-19. Companion to [T8.0 findings](m8-t80-findings.md). Implementation in
[`src/cardiac_redundancy/duration.py`](../src/cardiac_redundancy/duration.py) (25 tests).

This note derives the duration endpoint from the model's own accounting, states
every assumption it rests on, and says which of them are conservative and which
are optimistic. No claim about survival, recovery or injury is made anywhere in
it: a duration here is the interval over which the model's declared delivery
screens can be met.

---

## 1. Why demand is an oxygen number but the limit is a glucose number

M1 expresses whole-body demand as oxygen: 255.848 mL O₂/min at rest, allocated
across organs. The suppressed, supported state used in M8 requests about
108 mL O₂/min at the T8.0 optimum.

Meeting an oxidative demand requires both oxygen *and* a substrate to oxidise.
The model charges all oxidative demand to glucose at the stoichiometry of
complete oxidation,

> C₆H₁₂O₆ + 6 O₂ → 6 CO₂ + 6 H₂O,

so one mmol of glucose consumes 6 mmol of oxygen, which at STPD is

> 6 × 22.414 = **134.48 mL O₂ per mmol glucose**  (`O2_ML_STPD_PER_MMOL_GLUCOSE`).

A demand of D mL O₂/min therefore requires **D / 134.48 mmol glucose per minute**.
At D = 108.2 that is 0.8046 mmol/min.

Oxygen itself is not the binding resource while ventilation continues: M7 found
the oxygen-limited term is zero at every saved sample, and the deficit is
glucose-limited. That is why duration is a substrate calculation rather than an
oxygen-store calculation.

## 2. Three reservoirs with very different roles

| Reservoir | Size | Role |
|---|---|---|
| Circulating glucose | plasma 5.0 mmol/L ÷ 1.11 × 5.25 L = **23.65 mmol** | immediately available; drains |
| Hepatic glycogen | 229 mmol/L × 1.44 L = **329.8 mmol** | 14× larger, but access is rate-capped |
| Oxygen stores | lung gas plus ~5 s of each organ's demand | negligible on this timescale |

The whole-blood concentration is the plasma reading divided by the
plasma-to-whole-blood ratio of 1.11, which is how the model reports it
(`minimum_plasma_equivalent_mmol_L` multiplies back by that ratio).

## 3. The rate cap is the whole story

Hepatic export is capped at **0.21 mmol/min**, which in oxygen-equivalent terms is

> 0.105 × 134.48 = **14.1 mL O₂-eq/min**, where 0.105 = 0.21 × `liver_function_fraction` 0.5.

Against a demand of 108.2 mL O₂/min, export covers **26%**. The other 74% is
drawn from the circulating pool at a net rate

> **drain = D / 134.48 − R** = 0.8046 − 0.105 = **0.6996 mmol/min**,

and the supported duration is

> **T = usable pool / drain**.

At reference physiology with no floor that is 23.65 / 0.6996 = **33.8 min**.

### The store never empties, and that is a structural fact

> store ÷ export = 329.8 / 0.105 = **3141 min**, over two days,

against 33.8 min for the pool. The glycogen store outlasts the pool by a factor
of 93. **The binding constraint is the rate at which glucose can be mobilised,
not the quantity stored**, and that conclusion is insensitive to the store size:
it would hold for any hepatic store above roughly 2.5 mmol, under 1% of the
reference value. Arguing about glycogen content changes nothing; arguing about
export rate changes everything.

## 4. The floor enters linearly

With a declared plasma floor C_floor below the initial C₀ = 5.0 mmol/L, only the
glucose above the floor is counted:

> usable pool = (C₀ − C_floor) / 1.11 × 5.25 L,
> **T(C_floor) = (C₀ − C_floor) / 1.11 × V / (D/134.48 − R)**.

T is linear in the floor and reaches zero at C_floor = C₀, which is why the
curves in the duration figure are straight lines. A floor of 3.0 mmol/L leaves
two fifths of the pool, hence exactly two fifths of the duration: 13.5 min
against 33.8.

The model declares **no floor**, so its native endpoint is complete exhaustion.
That is not a physiological endpoint, and the evidence gives a range rather than
a cutoff: cerebral glucose transport begins limiting utilisation near 3 mmol/L,
the ADA/EASD joint statement sets 3.0 mmol/L (54 mg/dL) as the level to report
in trials, measured symptom thresholds cluster at 3 mmol/L, cognitive impairment
appears at 3.1–3.4 mmol/L, and transport half-saturation estimates span about
1.5 to 6.2 mmol/L with wide uncertainty. The floor is therefore a declared axis
with 3.0 mmol/L as an anchored reference point, never a fitted parameter.

## 5. Elasticities, and an exact identity

Writing T = P / (D/K − R) with P the usable pool and K = 134.48:

> ∂lnT/∂lnD = −(D/K)/(D/K − R) = −0.8046/0.6996 = **−1.150**
> ∂lnT/∂lnR = +R/(D/K − R) = 0.105/0.6996 = **+0.150**

These sum to exactly **−1** for any D and R, since the numerators differ by the
denominator. Measured by finite difference at a 10% step: −1.03 and +0.15.

The ratio of the levers is (D/K)/R = 0.8046/0.105 = **7.7**: a given fractional
reduction in demand buys 7.7 times what the same fractional gain in export does.
The pool empties because of demand, not for want of refilling. An earlier
expectation that export would dominate was wrong, and the measurement corrected it.

## 6. Why the design barely moves the clock

Total demand splits as D = D_organs + D_pump, and the auxiliary heart is only
**7–16%** of it (1.233 to 3.110 mL O₂ per 10 s against 16.4 for the organs).
Applying the elasticity:

> scale 0.4 → 1.0 raises D by 10.6% (105.6 → 116.8), predicting −12.6% in T.
> Observed: 34.8 → 31.0 min, **−10.9%**.

So the design enters duration only through its share of total demand, and that
share is small. Across the whole grid of scenario, floor and design, duration
spans 9.0 to 95.2 min; scenario contributes about 4×, the floor about 3.3×, and
the design about 1.15×.

**Duration is a physiology-and-reporting-convention result, not a design result.**
The design question is whether perfusion can be established at all; how long it
lasts is set by substrate logistics and by where the floor is drawn.

## 7. What the calculation assumes, in both directions

**Optimistic — these make the estimate too long:**

1. **Export already at its effective cap.** Measured export over the first
   minute was 0.105 mmol/min, which *is* the effective cap (0.21 × liver function
   0.5), so this assumption held exactly and the simulated drain of 0.704
   mmol/min gives 33.7 min against the analytic 33.8. An earlier version of this
   note applied 0.21 without the liver fraction, overstating duration by 18% and
   presenting the resulting mismatch as a 34–40 min uncertainty bracket. There is
   no such bracket: the analytic and simulated values agree.
2. **Demand constant.** Prescribed organ demand does not fall as substrate
   depletes, and contraction is prescribed, so the pump cannot starve.
3. **No floor by default**, as above.
4. **A fully accessible pool**, with no distribution or transport delay between
   compartments.

**Conservative — these make the estimate too short:**

5. **All oxidative demand is charged to glucose.** The model has no fatty-acid
   or ketone oxidation. Resting metabolism runs near a respiratory quotient of
   0.8, so roughly two thirds of energy normally comes from fat. Charging a
   fraction f of demand to glucose gives 33.8 min at f=1.0, 80 min at 0.5, 173 at
   0.3, and **no substrate limit at all below f=0.131**, where export alone covers
   demand. For scale, the glycogen store is about 244 kcal at this demand, or
   8 hours, against roughly 135,000 kcal in 15 kg of adipose, about 184 days. This
   is the single largest bias in the calculation by a wide margin, and it is a
   model-form limitation, not a physiological claim. M8 reports the glucose-only
   arm as its primary result and repeats its tasks under a declared mixed-fuel
   alternative.
6. **Hepatic export only, and only its glycogenolytic part.** The 0.21 figure is
   measured fasting **net glycogenolysis**, not maximal stimulated output
   (Magnusson et al. 1992, PMID 1401068), and excludes gluconeogenesis. Measured
   whole-body glucose production in fasted humans is around 2 mg/kg/min, about
   7.4 times the model's effective export. No muscle glycogenolysis and no lactate
   shuttle either.
7. **No anaerobic ATP** contribution at all.

The two directions do not cancel and should not be presented as if they do. The
honest statement is that the figure is a **first-order bound within a
glucose-only oxidative model**, verified against a 60 s coupled replay for its
drain rate, not a prediction of how long a person would be supported.

## 8. What was verified, and how

The binding limit among the release law's four terms was identified by
measurement, not assumed: mobilisation capacity times liver function (0.105
mmol/min), against the glycogen term (19,786) and hepatic delivery (0.539, at a
measured 539 mL/min hepatic outflow). The circulation delivers five times more
hepatic flow than export can use.

A 60 s coupled run at the optimum sustained **zero unmet demand in all three
organ groups** — active, preserved and the auxiliary myocardium — over the full
minute, with 0 of 57 sustained pressure windows violated. That run supplies the
measured export rate and the measured drain above. Simulating the full 34 min was never
attempted: **verifying the drain rate over a minute is what licenses the
extrapolation**, and the extrapolation remains conditional on that rate holding
as the pool falls.

---

## 9. The sweep currently running

Two searches, both over **design controls only**. Scenario factors, compression
kappa, activation latency and the glucose floor are never searched: a search
that could choose them would lower the requirement instead of meeting it.

**Phase A — pressure boundary re-location at the settled horizon.**

| | |
|---|---|
| Varies | cavity scale (all four scales tied) |
| At | 60, 80 and 100 beats/min |
| Fixed | kappa 0.5, oracle release 0.5 s, reference physiology, nominal conductance |
| Method | Brent root-find on the continuous pressure margin, bracket [0.30, 0.70], tolerance 0.01 in scale |
| Horizon | **40 s**, with a plateau check on the final five cycle means (0.05 mmHg) |
| Cost | ~7 evaluations per rate at ~10 min each |

It answers: *what is the least cavity scale whose settled cycle-mean aortic
pressure meets the screen, at each rate?*

**Phase B — joint scale and coronary conductance.**

| | |
|---|---|
| Varies | cavity scale and coronary conductance together |
| Ladder | scale from the Phase A boundary upward: +0, +0.03, +0.06, +0.10, +0.15 |
| Inner search | bisect conductance over [0.5, 2.0] to 0.05 for the least value meeting auxiliary self-supply, re-checking pressure at every point |
| Reported | the cheapest feasible pair by auxiliary oxygen demand, plus an epsilon-constraint frontier over conductance ceilings |

It answers: *what is the cheapest design that satisfies pressure **and** feeds
its own myocardium, and how much coronary enlargement does that require?*

The ladder exists because the two constraints do not separate at the boundary:
raising conductance diverts flow from the aorta and costs 1–2 mmHg, which a
design sitting exactly on the pressure boundary cannot afford.
