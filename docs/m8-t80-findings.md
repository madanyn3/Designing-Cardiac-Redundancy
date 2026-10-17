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

| Rate | Least scale | Aux demand (mL O₂/min) | Work (J) | MAP | Margin at boundary |
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

## 7. What is not established

- **No robustness.** Protocol v2 requires the worst declared scenario; only
  reference physiology has been run. The adverse triple may remove the optimum.
- **No kappa robustness.** The optimum is conditional on kappa 0.5 and the
  kappa 1 case was infeasible by a wide margin in the region tested.
- **Oracle activation does the heavy lifting.** At the M6 controller's 2.5 s
  detection latency nothing passes at any scale. The real requirement handed
  back to M6 is a detection-plus-recruitment budget, and the feasible region in
  (capacity, latency) has not been mapped.
- **Three design controls remain pinned**: left/right asymmetry, contractility
  and preservation conductance.
- **No verification search.** The protocol requires a matched-budget
  differential-evolution and space-filling comparison in the neighbourhood of
  the located optimum, since a structural argument cannot verify itself.
- **Cavity scale is a capacity proxy**, not a determined anatomical mass. Scale
  0.651 is not a claim about implantable tissue volume.
- **The pending cycle-alignment amendment is unresolved**, and now matters more:
  `docs/m2-pressure-requirements.md` asks for per-activation-cycle evaluation,
  while the frozen policy uses 1 s windows, which mix 1.33 cycles at 80 /min and
  2.5 at 150 /min. Rate is now a live design variable, so this inconsistency
  touches every comparison across rates in §2.

## 8. Cost model correction

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
