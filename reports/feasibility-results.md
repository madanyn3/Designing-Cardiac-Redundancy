# M8 feasibility results

2026-10-23. Protocol version 5. Arm A, the glucose-only primary analysis.
Figures in `reports/m8-figures/`; searches in `reports/m8-optimize/`,
`reports/m8-joint/`, `reports/m8-verification/`, `reports/m8-matrix/`.

**Every result below is conditional on five declared settings**, and none
survives their removal without being re-established:

| Condition | Value |
|---|---|
| Compression kappa | 0.5 (declared cases 0, 0.5, 1) |
| Physiology scenario | reference (of reference, adverse, favourable) |
| Activation | oracle, 0.5 s prescribed release |
| Screen horizon | 40 s, plateau confirmed |
| Fuel | all oxidative demand charged to glucose |
| Pressure margin | at least 3.0 mmHg above the enforced threshold (protocol v5) |

The nine unresolved requirements of `configs/m4-constraint-policy.json` remain
attached, so the formal status of even a clean run is *indeterminate*, not
*satisfied*. No claim of survival, neurological outcome or clinical
applicability is made anywhere in this report.

---

## 1. A feasible design exists, and it is not unique

Three independent searches agree. Boundary tracing located an optimum; a
matched-budget sample beat it; a joint search over all three design axes
reproduced the sample's value from a different direction.

| Method | Scale | Rate | Conductance | Auxiliary O₂ demand |
|---|---|---|---|---|
| Boundary tracing (discrete rates) | 0.496 | 60.0 /min | 1.203 | 7.812 mL/min |
| Matched-budget Sobol sample | 0.566 | 49.2 /min | **1.027** | **7.600** |
| Joint search (Sobol + compass) | 0.536 | 53.4 /min | 1.808 | **7.601** |

The two unconstrained searches land within **0.001 mL O₂/min** of each other
while differing materially in every design variable. The optimum is a **shallow
valley**, not a point: among the joint search's best eight designs, conductance
ranges 1.388–1.808 while demand varies only 7.601–7.614, a spread of 0.17 per
cent. The compass refinement improved on its own seed by 0.1 per cent across
four iterations.

### The unconstrained optimum is not a reportable design

**It holds +0.005 mmHg of pressure margin.** Changing its activation period by
0.0005 s, four hundredths of one per cent, flips it to violated with two bad
windows. The runner-up at 7.601 holds +0.040 mmHg and behaves the same way.

This is a **specification error, not a search error**. Minimizing a flat
objective against a hard threshold presses the answer against the boundary, and
the boundary is where the verdict is decided by window alignment rather than by
the design. Nine of 152 feasible designs hold under 0.5 mmHg and the cheapest
are all among them: the search concentrated on fragility by construction, and
it did so correctly, because nothing in the objective asked for margin.

Protocol v5 therefore declares a minimum margin of 3.0 mmHg, before the re-run
rather than after seeing which designs survive it. **Robustness is cheap:**

| Requirement | Cheapest design meeting it | Scale | Rate | Conductance | Margin |
|---|---|---|---|---|---|
| Feasible only | 7.600 mL O₂/min | 0.566 | 49.2 /min | 1.027 | **+0.005** |
| **Margin ≥ 3.0 mmHg** | **8.137** | 0.536 | 57.5 /min | 1.748 | +3.10 |
| Margin ≥ 3.0 and conductance ≤ 1.00 | 9.596 | 0.749 | 46.0 /min | 0.953 | **+8.97** |

Margin costs **7.1 per cent** more auxiliary oxygen. The reported design is
therefore **scale 0.536, 57.5 beats/min, conductance 1.748, at 8.137 mL O₂/min
with +3.10 mmHg of margin** — or the nominal-conductance design below, which
costs 17.9 per cent more and depends on one assumption fewer.

Thirty-six designs hold the margin and feed their own myocardium. The cheapest
few margin-holders fail self-supply, because margin wants a larger or faster
pump while self-supply wants a larger coronary bed, so the two cannot both be
bought at the cheap end.

Boundary tracing's 7.812 sits 2.7 per cent above the valley. The cause is
recorded rather than smoothed over: it traced boundaries only at 60, 80 and
100 beats/min, and the optimum lies near 49–53, between those grid lines. The
structural argument was sound; the rate axis was sampled too coarsely. **This
is what the protocol's verification requirement is for**, and it is reported as
a finding about the method, not hidden.

Where two designs cost the same, the one depending least on an unjustified
assumption is preferred (§2).

## 2. Feasibility does not depend on the coronary conductance assumption

The 0.5–2.0 conductance bracket is an exploratory assumption with no measured
basis, and `docs/uncertainty-plan.md` states that using it as a design control
requires a verified interface that does not exist. Pooling both searches:

| Conductance band | Cheapest margin-holding feasible design | Scale | Rate | Margin |
|---|---|---|---|---|
| **0.50–1.00 (at or below nominal)** | **9.596 mL O₂/min** | 0.749 | 46.0 /min | **+8.97** |
| 1.00–1.25 | 8.373 | 0.516 | 62.6 /min | +5.72 |
| 1.25–1.50 | 8.396 | 0.550 | 58.0 /min | +4.94 |
| 1.50–2.00 | **8.137** | 0.536 | 57.5 /min | +3.10 |

**A feasible design exists at or below nominal coronary conductance, and it is
the most robust design found.** It costs 17.9 per cent more oxygen than the
margin-constrained optimum and needs a 40 per cent larger cavity, but it holds
nearly three times the required pressure margin. The assumption buys
optimality, not feasibility, which is the materially weaker dependency and
should be stated that way.

Every band is populated, so there is no conductance value at which the design
problem becomes infeasible within the bracket. Margin *falls* monotonically as
conductance rises, which is the trade the bracket actually buys: a larger
coronary bed lets a smaller cavity feed itself, and a smaller cavity holds less
pressure margin.

## 3. The binding constraint is detection, not the pump

The core matrix held the design fixed across nine declared operating cases.
**One satisfied every screen: oracle activation.**

| Case | Pressure | Organs | MAP | Release |
|---|---|---|---|---|
| **oracle_activation** | **satisfied** | **satisfied** | 67.70 | 0.50 s |
| support_with_suppression | violated | pump | 67.63 | 2.55 s |
| delayed_activation | violated | pump | 67.63 | 2.55 s |
| ventilation_degraded | violated | pump | 67.63 | 2.55 s |
| partial_recovery | violated | pump | 105.08 | 2.55 s |
| support_no_suppression | violated | active, pump | 46.06 | 1.95 s |
| shared_failure | violated | active | 51.70 | 2.55 s |
| low_preload | violated | active, pump | 65.06 | 2.60 s |

`support_with_suppression` reaches MAP 67.63, **above** the enforced threshold,
and still fails: the observed controller releases at 2.55 s and the window
containing release lies outside the declared 2 s startup allowance. The design
is not what fails these cases.

The latency relation quantifies the requirement. Time to adequacy is linear in
latency with a design-dependent slope, because a delayed release both loses
windows and accrues an arterial volume debt to repay:

| Design | Slope (s of adequacy per s of latency) | Tolerable release |
|---|---|---|
| scale 0.387 at 80 /min | 2.00 | 1.0 s |
| scale 0.475 at 60 /min | 3.26 | 1.0 s |
| scale 0.650 at 60 /min | 1.60 | **1.5 s** |

**Every tolerable value lies below the M6 controller's 2.55 s.** Capacity buys
latency tolerance, so the requirement handed back to detection is a
detection-plus-recruitment budget under roughly 1.5 s, at a cavity above the
cost-optimal one.

## 4. Supported duration, and what sets it

| Scenario | Demand | Effective export | No floor | At 3.0 mmol/L floor |
|---|---|---|---|---|
| Favourable | 87.2 mL O₂/min | 0.400 mmol/min | 95.2 min | 38.1 min |
| **Reference** | 108.2 | 0.105 | **33.8 min** | **13.5 min** |
| Adverse | 144.8 | 0.025 | 22.5 min | 9.0 min |

A 60 s coupled replay sustained **zero unmet demand in all three organ groups**
and measured an export of 0.105 mmol/min, exactly the effective cap, giving a
drain consistent with the analytic 33.8 min. Analytic and simulated agree.

Three facts about this endpoint:

**The export rate binds, not the store.** Of the release law's four terms,
mobilisation capacity times liver function (0.105 mmol/min) binds; the glycogen
store term is 19,786 and hepatic delivery 0.539 at a measured 539 mL/min
hepatic outflow. The circulation delivers five times more hepatic flow than
export can use.

**The design barely moves it.** Elasticities are demand −1.03 and export +0.15,
and the pump is 7–16 per cent of total demand, so scale 0.4 to 1.0 spans only
34.8 to 31.0 min. Duration is a physiology-and-reporting result, not a design
result.

**The floor is a declared axis, not a value.** With none, duration means time to
complete exhaustion of the circulating pool, which is not a physiological
endpoint. 3.0 mmol/L is anchored to the trial-reporting level, measured symptom
thresholds, and where cerebral glucose transport begins to limit utilisation;
transport half-saturation estimates span 1.5–6.2 mmol/L, which is why a curve
is reported.

## 5. Where the envelope ends

| Boundary | Finding |
|---|---|
| Compression kappa | **Both** alternatives fail, for opposite reasons. At kappa 1 self-supply fails; at kappa 0 the pressure screen fails, because less compression means more coronary runoff and MAP falls to 63.44 against the 65 threshold. Feasibility holds only in a narrow band around 0.5. |
| Physiology | Feasible in **2 of 3** declared scenarios. Adverse raises auxiliary demand 54-62 per cent and fails the **active** group as well as the pump, in all four leading designs. |

T8.4 tested the four leading designs against five stresses: 4 of 4 survive
reference and favourable, **0 of 4 survive adverse, kappa 0 or kappa 1**. It
also found a **ranking reversal under both adverse and favourable physiology**,
in which the design cheapest at reference drops to third. Selecting on
reference alone was therefore the wrong criterion, and the reversal is reported
rather than resolved: resolving it would mean optimizing against one declared
scenario over the others, and none has a stronger claim.
| Detection latency | Tolerable 1.0–1.5 s against the controller's 2.55 s (§3). |

Kappa is the fraction of auxiliary chamber transmural pressure applied as its
own coronary closing pressure. Within one beat it drives coronary flow from
1.33 mL/s in diastole to 0.41 mL/s at peak systole, with the minimum coinciding
with near-peak aortic pressure. It had never been exercised by any runner
before T8.2, despite M2 declaring the cases and the uncertainty registry
assigning them to T7.4. **It is the single most influential quantity in the
self-supply constraint and rests on no measurement.**

## 6. What the earlier milestones concluded, and why it changed

M5 and M7 reported that no tested configuration reached the pressure screen,
with 0 of 84 global cases passing. That was true of the design space those
milestones sampled, which held cavity scale fixed at the 0.25 reference fixture
and varied conductance and contractility instead. Scale was always a legal
parameter; it had simply never been varied, because M5 froze it as a test
fixture and M7 inherited the choice. Inventory never constrained it: every
scale in 0.1–1.0 is admissible on the declared initial state.

Three further quantities were frozen at the same fixture and had to be coupled
before optimization was meaningful: the auxiliary mass ratio, the activation
rate charged in demand, and the auxiliary oxygen store.

Two screen definitions also changed, on evidence rather than convenience. The
enforced criterion became mean arterial pressure alone, because adult
post-arrest guidance states a MAP target and no systolic one while the systolic
and diastolic thresholds derive from one perioperative cohort; the strict
triple is still evaluated and reported for every run. And pressure windows
became cycle-aligned, because a fixed 1 s grid mixes 1.33 cycles at 80 /min and
2.5 at 150, making cross-rate comparisons incommensurable.

## 7. What is not established

- **Arm B is not included here.** All oxidative demand is charged to glucose,
  while resting metabolism draws most oxidative ATP from fat. This is the
  largest single bias in the duration figures and is addressed separately.
- Only reference physiology and kappa 0.5 are satisfied; the envelope ends
  where §5 says it does.
- Oracle activation does the heavy lifting; at the real controller latency
  nothing passes.
- Three design controls remain pinned: left/right asymmetry, contractility and
  preservation conductance.
- Cavity scale is a capacity proxy, not a determined anatomical mass.
- Regional pressure minima, CO₂ and pH, ischaemic and reperfusion injury,
  thrombosis and blood compatibility are not modelled.
- All results are *best found within the declared bounds and budget*. No global
  optimum is claimed.
