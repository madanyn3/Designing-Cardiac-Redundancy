# T8.0 — Joint design optimization: results

2026-10-17. Protocol version 2. Code commit recorded in
[`configs/m8-protocol.json`](../configs/m8-protocol.json). Outputs under
`reports/m8-optimize/`.

**Every number below is conditional on four declared settings**, and none of the
conclusions survive their removal without being re-established:

| Condition | Value | Why it matters |
|---|---|---|
| Compression kappa | **0.5** | the most influential quantity in the self-supply constraint; at 1.0 the design space tested so far is infeasible |
| Physiology scenario | **reference** | protocol v2 requires judging a candidate on its worst declared scenario; adverse has not been run |
| Activation | **oracle, 0.5 s** | the M6 controller needs 2.5 s, at which nothing passes at any cavity scale |
| Horizon | **10 s** | a pass means "no modeled violation observed through 10 s", never support duration |

The nine unresolved requirements of `configs/m4-constraint-policy.json` remain
attached to every verdict, so the formal status of even a clean run is
`indeterminate`, not `satisfied`.

---

## 1. What was optimized

```
minimize    E(x) = total auxiliary oxygen demand (mL O2/min)
over        x = (cavity scale, period, coronary conductance)
subject to  the pressure triple (MAP 65, SBP 90, DBP 50) in every sustained window
            auxiliary self-supply: no sustained unmet demand in the pump group
```

Kappa, extraction, specific mass, work cost, suppression depth and substrate
capability are **not** in `x`. A search able to choose them would lower the
requirement rather than meet it.

Method: boundary tracing, justified by three measured properties — the objective
rises monotonically in cavity scale and rate, so the optimum lies on a constraint
boundary; conductance is nearly free energetically; and the constraints were
expected to nearly separate. **The third property turned out to be false at the
boundary** (§4), which the method handled by searching a two-dimensional ladder.

## 2. Phase A — the pressure boundary

Least cavity scale satisfying the pressure triple at each rate, with the
objective read from the same runs (55 evaluations, 18.9 min of compute):

| Rate | Least scale | Aux demand (mL O2/min) | Work (J) | MAP | Margin at boundary |
|---|---|---|---|---|---|
| 150 /min | 0.417 | 11.967 | 8.186 | 84.54 | +0.484 mmHg |
| 120 /min | 0.438 | 10.971 | 7.636 | 83.12 | +0.108 |
| 100 /min | 0.473 | 10.593 | 7.420 | 83.00 | +0.807 |
| 80 /min | 0.522 | 10.067 | 6.998 | 80.30 | +0.415 |
| **60 /min** | **0.621** | **9.769** | 6.550 | 77.56 | +0.343 |
| 48 /min | 0.803 | 10.630 | 6.752 | 78.33 | +0.243 |
| 40 /min | **none in bounds** | — | — | — | — |

Every boundary was verified with the failing point 0.007 below it, so the
monotonicity premise bisection relies on held at all six feasible rates. Each
boundary was then re-read through an independent continuous margin function and
found within 0.11–0.81 mmHg of the screen — two routes agreeing.

**MAP at the boundary is 77–85 mmHg, far above 65.** The binding metric is
**SBP ≥ 90**, not the mean the project has reported since M5.

**The objective is U-shaped in rate with an interior minimum.** Slower is not
monotonically better: 48 /min needs a 29% larger cavity and costs more oxygen
than 60 /min, and at 40 /min no scale in the declared bounds passes at all.
Mechanically, slower beating cuts the rate-proportional activation cost but
demands a larger cavity whose mass cost eventually dominates.

Refinement at the parabolic vertex (67 /min) and at 70 /min gave 9.881 and
9.757 mL O₂/min — a **1.3% spread across 60–70 /min**. The minimum is a shallow
basin, not a point, and is reported as such. Quoting an optimal rate more
precisely than "roughly 60–70 /min" would be false precision.

## 3. Phase B — the jointly feasible optimum

Adding the self-supply constraint and searching conductance jointly (48.3 min):

| Rate | Scale | Conductance | Aux demand | Work | MAP | Sustained unmet |
|---|---|---|---|---|---|---|
| **60 /min** | **0.651** | **0.969** | **10.257** | 6.893 | 79.77 | **0.000000** |
| 80 /min | 0.552 | 1.203 | 10.662 | — | 82.64 | 0.000000 |
| 60 /min | 0.681 | 0.969 | 10.743 | — | 81.85 | 0.000000 |
| 80 /min | 0.582 | 1.203 | 11.260 | — | 85.21 | 0.000000 |

Eight jointly feasible pairs were found. The best design is **cavity scale
0.651, 60 beats/min, coronary conductance 0.969**, costing **10.257 mL O₂/min**.

Against the first design that passed everything — scale 1.0 at 80 /min, 19.179
mL O₂/min — this is **47% less metabolic cost at a 35% smaller cavity**.

### The epsilon-constraint frontier is flat, and that is the important result

| Conductance ceiling | Cheapest feasible design |
|---|---|
| ≤ 1.00 | scale 0.651 at 60 /min, 10.257 mL O₂/min |
| ≤ 1.25 | identical |
| ≤ 1.50 | identical |
| ≤ 1.75 | identical |
| ≤ 2.00 | identical |

The minimum conductance the optimum needs is **0.969 — below nominal.** The
constraint is not binding, so **the result does not depend on the coronary
conductance assumption at all.** That assumption was the weakest link: the
0.5–2.0 bracket is an exploratory "A" scenario with no measured basis, and
`docs/uncertainty-plan.md` states that using it as a design control "requires
verified interface", which does not exist. The optimum needing no coronary
enlargement removes that dependency from the conclusion.

## 4. The constraints do not separate at the boundary

At the cheapest pressure-feasible scale — 0.621 at 60 /min and 0.522 at 80 /min —
**no conductance in the declared range makes the pump feed itself.** Both are
recorded as `self-supply unmet at the upper conductance bound`. The jointly
feasible design needs about +0.03 in scale above the pressure boundary.

This retires the near-separability premise that motivated the staged method.
Raising conductance diverts flow from the aorta and costs 1–2 mmHg (87.1 → 85.1
at scale 0.6), which a design sitting exactly on the pressure boundary cannot
afford. The two constraints bind in different places and the optimum lies on
neither boundary alone.

## 5. Compression kappa dominates everything

At scale 0.6, conductance 1.0, reference physiology:

| kappa | Pump sustained unmet | Status |
|---|---|---|
| 0 | 0.000000 mL | satisfied |
| 0.5 | 0.002351 mL | violated |
| 1 | 0.300505 mL | violated — 128× worse |

Kappa is the fraction of the auxiliary chamber's transmural pressure applied as
the closing pressure of its own coronary bed. Within one beat at scale 0.6 it
drives coronary flow from 1.33 mL/s in diastole down to 0.43 mL/s at peak
systole — a 3.3× swing, with the minimum occurring when aortic pressure is
near its *maximum*. With no tissue substrate buffer (L-M4-F01) demand must be
met at every instant, so that systolic dip is charged directly as unmet demand.

At kappa 0 even a 0.5-scale pump feeds itself completely. At kappa 1 a
0.6-scale pump misses 18% of its own request. **The feasibility conclusion is
governed by an untested model-form choice**: `configs/m2-frozen-overlay.yaml`
declared the cases 0/0.5/1 in M2 and `docs/uncertainty-plan.md` assigned them to
T7.4, but no runner exercised them before T8.2. This belongs at the top of the
T7.5 evidence priorities.

## 6. What this changes about the project's earlier conclusion

M5 and M7 reported that no tested configuration reached the 65 mmHg screen, with
0 of 84 global cases passing. That conclusion was **true of the design space
those milestones sampled**, which held cavity scale fixed at the 0.25 reference
fixture and varied conductance and contractility instead. Scale was always a
legal parameter in `BiologicalParameters`; it had simply never been varied,
because M5 froze it as a test fixture and M7 inherited that choice.

Three further quantities were frozen at the same fixture and had to be coupled
before optimization was meaningful: the auxiliary mass ratio, the activation
rate (hardcoded at 80 /min), and the auxiliary oxygen store. See
[`docs/m8-feasibility-criteria.md`](m8-feasibility-criteria.md).

## 7. Operating envelope: compression kappa

Phase B repeated at kappa 1 across both rates and all five ladder rungs:
**no design is feasible at any conductance in the declared range.** The
epsilon frontier is empty at every ceiling from 1.00 to 2.00, and every rung
is recorded as `self-supply unmet at the upper conductance bound`.

| kappa | Outcome in the region searched |
|---|---|
| 0 | self-supply satisfied even at cavity scale 0.5 |
| 0.5 | feasible at scale 0.651 and above, conductance 0.969 |
| 1 | no feasible design at any scale or conductance tested |

So **kappa <= 0.5 is a boundary of the operating envelope**, stated as a
condition on the result rather than as a failure of the design. Kappa is the
fraction of the auxiliary chamber's transmural pressure applied as the closing
pressure of its own coronary bed; within one beat at scale 0.6 it drives
coronary flow from 1.33 mL/s in diastole to 0.43 mL/s at peak systole, with the
minimum coinciding with near-peak aortic pressure. With no tissue substrate
buffer (L-M4-F01) that systolic dip is charged directly as unmet demand.

## 8. Operating envelope: physiology scenarios

Designs held fixed, coupled replay only, since no scenario factor enters the
hydraulic model (16 replays):

| Design | Reference | Adverse | Favourable |
|---|---|---|---|
| 0.651 @ 60 /min | 10.257 feasible | 15.843, violated (active + pump) | 7.921 feasible |
| 0.552 @ 80 /min | 10.662 feasible | 16.452, violated (pump) | 8.226 feasible |
| 0.681 @ 60 /min | 10.743 feasible | 16.592, violated (active + pump) | 8.296 feasible |
| 0.672 @ 80 /min | 13.038 feasible | 20.116, violated (pump) | 10.058 feasible |

Demand (mL O2/min) in each cell. Feasible in **2 of the 3 declared scenarios**,
reported as a fraction of tested scenarios and never as a probability. The
auxiliary heart's own cost swings 54-62% above reference under adverse
physiology and 23% below under favourable -- a wider spread than any design
choice produced.

At 60 /min the **active** group fails under adverse physiology, not only the
pump: brain, respiratory muscle and lungs stop being adequately supplied. The
80 /min designs fail on self-supply alone. The cheapest designs at reference
are therefore the most fragile.

## 9. Activation latency: the relation

Prescribed release stands for detection plus recruitment. Settling is measured
from release to the first sustained window satisfying the full triple:

| Design | Slope (s of adequacy per s of latency) | Max tolerable release |
|---|---|---|
| 0.651 @ 60 /min | 3.26 | 0.5 s |
| 0.651 @ 80 /min | 2.00 | 1.0 s |
| 0.800 @ 60 /min | 1.60 | 1.5 s |

A second of latency costs 1.6-3.3 s before the circulation is adequate: the
delay plus recovery from the arterial volume debt it accrued. **Capacity buys
latency tolerance** -- scale 0.651 to 0.800 at 60 /min cuts the slope from
3.26 to 1.60 and triples the tolerable delay. The cheapest design is the least
latency-tolerant, which is why the optimum requires an oracle detector.

**SBP is the last metric to satisfy in every case.** At 0.651 @ 60 /min with a
2 s release, MAP 65 is met from 3.0 s while SBP 90 is met only from 7.0 s. The
entire latency budget, and the binding constraint on cavity size, are set by
the systolic threshold.

## 10. Supported duration: the primary endpoint

The project's result is the interval over which vital-organ perfusion can be
sustained, not whether one instant satisfies a screen. The closed-form bound is
in [`src/cardiac_redundancy/duration.py`](../src/cardiac_redundancy/duration.py)
(20 tests).

Oxidative demand must be met by glucose; hepatic export is capped by a **rate**,
not by the store; the difference is drawn from the circulating pool.

| Scenario | Total demand | Net drain | **Supported duration** |
|---|---|---|---|
| Favourable | 87.2 mL O2/min | 0.248 mmol/min | **95.3 min** |
| Reference | 108.2 | 0.594 | **39.8 min** |
| Adverse | 144.8 | 0.977 | **24.2 min** |

The design trade runs as expected -- a larger pump costs duration:

| Design | Demand | Duration (reference) |
|---|---|---|
| scale 0.4 | 105.6 | 41.1 min |
| scale 0.651 | 109.5 | 39.1 min |
| scale 1.0 | 116.8 | 35.9 min |

**The export rate binds, not the store.** The liver holds 329.8 mmol, which at
the 0.21 mmol/min cap would last 1,570 min; export covers only 28.2 mL O2-eq/min
against 108 demanded, so the 23.7 mmol circulating pool sets the clock.

Measured elasticities: **demand -1.19, export +0.36.** Demand dominates because
the glucose draw is roughly four times the export cap, so the pool is emptied by
demand rather than refilled by export. The pump is 10-16% of total demand, so
suppression depth governs duration while the pump's share moves it
proportionally.

This is a **first-order bound** -- constant demand, export at its cap, no
declared glucose floor -- and is not yet a result. Transport lags
(`release_time_s`, `blood_target_time_s` of 60 s) and startup transients all
shorten it. It must be verified against a coupled long-horizon replay.

## 11. What is not established

- **The envelope is bounded, not open.** Feasibility holds at reference and
  favourable physiology with kappa <= 0.5; adverse physiology and kappa 1 lie
  outside it. These are conditions on the claim, established in sections 7-8.
- **The duration figures are an analytic bound**, not simulated outcomes
  (section 10).
- **Oracle activation does the heavy lifting.** At the M6 controller's 2.5 s
  detection latency nothing passes at any scale. The real requirement handed
  back to M6 is a detection-plus-recruitment budget, and the feasible region in
  (capacity, latency) has not been mapped.
- **Three design controls remain pinned**: left/right asymmetry, contractility
  and preservation conductance.
- **No verification search.** The protocol requires a matched-budget
  differential-evolution and space-filling comparison in the neighbourhood of
  the located optimum, since a structural argument cannot verify itself.
- **The pressure screen's three thresholds are not equally supported.** AHA
  2025 adult post-arrest guidance is expressed as a minimum MAP of 65 mmHg, with
  four randomized trials finding no benefit from higher MAP targets, and no adult
  systolic target is stated; systolic targets appear only in the pediatric
  guidance. SBP 90 and DBP 50 come from one perioperative observational cohort
  with radial measurement and a five-minute exposure definition. Since the
  systolic threshold binds both the cavity requirement (section 2) and the
  latency budget (section 9), a declared per-condition policy -- MAP enforced,
  systolic and diastolic reported with provenance -- is a pending amendment.
  An earlier draft of this work stated that AHA 2025 removed a previous systolic
  target; that could not be confirmed from the sources and is withdrawn.
- **Cavity scale is a capacity proxy**, not a determined anatomical mass. Scale
  0.651 is not a claim about implantable tissue volume.
- **The pending cycle-alignment amendment is unresolved**, and now matters more:
  `docs/m2-pressure-requirements.md` asks for per-activation-cycle evaluation,
  while the frozen policy uses 1 s windows, which mix 1.33 cycles at 80 /min and
  2.5 at 150 /min. Rate is now a live design variable, so this inconsistency
  touches every comparison across rates in §2.

## 12. Cost model correction

The protocol's measured single-evaluation cost (277 s) was taken at the
reference design under observed activation. Phase A measured **2.88 min per
hydraulic evaluation** under oracle activation with five parallel workers — 1.6×
higher, because releasing at 0.5 s instead of 2.5 s adds two seconds of
auxiliary pumping and the workers contend.

Phase A's 465-minute wall time was **not** compute: two 436-minute gaps in the
per-evaluation timestamps show the machine suspended overnight. Actual compute
was 18.9 min for 55 evaluations. Wall-clock estimates on this machine need to
account for suspension; compute estimates were accurate.

The Brent root-find on the continuous margin took **7 evaluations per line**
against bisection's 9 — a real but smaller saving than the ~5 predicted,
because the bracket used was wider than necessary.
