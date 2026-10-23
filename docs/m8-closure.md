# M8 closure — frozen feasibility experiments

2026-10-22. **M8 is closed at its declared experiment gate.** This is not a
demonstrated auxiliary heart, a validated repair window, or a survival claim.
M9 is the next milestone; preprint drafting has begun with the Methods section
and the Results material assembled here.

| Task | Completed evidence | Qualification |
|---|---|---|
| T8.1 | [Frozen protocol](m8-t81-protocol.md), [`configs/m8-protocol.json`](../configs/m8-protocol.json) | Version 4 after four recorded amendments, each with its reason, effect size and required re-runs |
| T8.0 | [Boundary trace, joint search and matched-budget verification](../reports/feasibility-results.md) | Best found within declared bounds and budget; no global optimum claimed |
| T8.2 | [Core matrix](../reports/m8-matrix/matrix-kappa0.5.json), nine operating cases | One case of nine satisfies every screen |
| T8.3 | [Arm A results](../reports/feasibility-results.md), [Arm B results](../reports/feasibility-results-armb.md), seven figures | Figures regenerate from saved output |
| T8.4 | [Arm A disconfirmation](../reports/m8-disconfirm/disconfirmation.json), [Arm B disconfirmation](../reports/m8b-disconfirm/disconfirmation.json) | Fractions of tested cases, never probabilities |
| T8.5 | [Answers to the research questions](m8-t85-answers.md) | Conditional on the declared envelope throughout |

460 regression tests pass.

## What the milestone establishes

**A feasible configuration exists, inside a narrow envelope.** A design
satisfying the pressure screen and auxiliary self-supply exists at cavity scale
0.566, 49.2 beats/min and essentially nominal coronary conductance, costing
7.600 mL O₂/min. Three independent searches agree, two of them landing within
0.001 mL O₂/min of each other from different directions. The optimum is a
shallow valley rather than a point: designs differing in all three variables
cost within 0.2 per cent of one another.

**The earlier negative result was an artefact of an unexamined constant.** M5
and M7 reported 0 of 84 cases reaching the pressure screen while holding cavity
scale at the 0.25 reference fixture. Scale was always a legal parameter, and
inventory never constrained it: every value in 0.1-1.0 is admissible on the
declared initial state. Three further quantities were frozen at the same
fixture and had to be coupled before optimization was meaningful — the
auxiliary mass ratio, the activation rate charged in demand, and the auxiliary
oxygen store.

**Detection, not the pump, is the binding constraint.** Eight of nine operating
cases fail, and the clearest reaches MAP 67.63, above the enforced threshold,
failing only because the controller releases at 2.55 s against a declared 2 s
allowance. Tolerable release is 1.0 to 1.5 s depending on cavity size. The
study converts a modelling assumption into a design requirement and hands it
back to detection; it cannot say whether that requirement is achievable.

**Feasibility does not depend on the coronary conductance assumption.** A
feasible design exists at or below nominal conductance, at 9.596 mL O₂/min and
cavity scale 0.749, against 7.600 at scale 0.566 just above nominal. The
unjustified 0.5-2.0 bracket buys optimality, not feasibility.

**Supported duration is a physiology result, not a design result.** The
hepatic export rate binds, not the glycogen store: the store would last
3,141 min at the effective cap against 33.8 for the circulating pool.
Elasticities are demand −1.03 and export +0.15, and the pump is 7 to 16 per
cent of total demand, so the whole admissible scale range spans 34.8 to 31.0
min. The reported figure moves by a factor of three on the declared glucose
floor, which the model does not set.

**Fat changes the fuel budget, not the machine.** Arm B raises supported
duration from 33.7 to 47.1 minutes by lowering the glucose fraction from 1.000
to 0.751. It does **not** admit a smaller pump: Arm B's optimum design is
feasible under Arm A's glucose-only metabolism as well, at 7.627 mL O₂/min
against 7.600, so the apparent size difference between the arms is an artefact
of two searches landing in different parts of a flat valley. Adipose perfusion,
not lipolytic capacity, limits the fuel at every flow examined, and feasibility
does not rest on the adipose policy: all three leading Arm B designs survive
with the depot suppressed.

**The envelope is narrow and was mapped by deliberate disconfirmation.** Of
four leading Arm A designs under five stresses, four survive reference and
favourable physiology and none survives adverse, kappa 0 or kappa 1. Both
compression alternatives fail for opposite reasons: at kappa 1 self-supply
fails, while at kappa 0 the pressure screen fails, because less compression
means more coronary runoff and MAP falls to 63.44. Feasibility holds only in a
narrow band around kappa 0.5, which rests on no measurement.

**The selection criterion did not survive its own test.** The ranking by
objective at reference inverts under both adverse and favourable physiology:
the design chosen as cheapest drops to third in each. The reversal is reported
rather than resolved, because resolving it would mean optimizing against one
declared scenario over the others and none has a stronger claim.

## Protocol amendments, and why each was made

| Version | Change | Reason |
|---|---|---|
| 2 | Auxiliary oxygen store sized by the shared five-second rule | The frozen pair was that rule evaluated once at the reference fixture; held fixed it aborted larger designs |
| 2 | Unresolved requirements attached to every policy | A clean numeric run would otherwise read as satisfied rather than indeterminate |
| 2 | Compression kappa declared as a model-form axis | Declared in M2, assigned to T7.4, never exercised by any runner |
| 3 | Enforced criterion is mean arterial pressure alone | Adult post-arrest guidance states a MAP target and no systolic one; the strict triple is still reported for every run |
| 3 | Cycle-aligned pressure windows | A fixed grid mixes 1.33 cycles at 80 beats/min and 2.5 at 150 |
| 4 | Horizon raised to 40 s with a plateau check | Cycle-mean pressure declines toward steady state, so a short horizon credits a value the design does not maintain |

Each amendment records its effect size. The horizon correction, for instance,
moves the located boundary by about one per cent of cavity scale: a
verdict-labelling correction at the margin, not a change of conclusion.

## Conditions carried into M9

1. Report the envelope, not the point. Every headline number is conditional on
   kappa 0.5, reference physiology, sub-second activation and the 40 s horizon,
   and the design fails outside them.
2. Report duration against the declared glucose floor as an axis. The model's
   native endpoint is exhaustion, which is not physiological, and the choice
   moves the answer by a factor of three.
3. Keep both arms visible. Arm A is the primary analysis; Arm B is a declared
   structural alternative. They agree on the machine and differ on the fuel
   budget: do not attribute design differences between them to the fuel model
   without a direct cross-check, which refuted exactly that claim once.
4. State that a feasible design exists at nominal coronary conductance, so
   feasibility does not rest on that bracket.
5. State the latency requirement as a requirement handed to detection, never as
   an achieved capability.
6. Do not describe any deficit as an injury dose, any duration as a survival
   time, or any pressure pass as brain protection. Regional pressure minima
   remain unresolved.
7. Carry forward R-M3-01, L-M4-F01, prescribed contraction, the absent
   fatty-acid and ketone pathways in Arm A, the assumed arterio-venous
   difference and depot mass in Arm B, and the unmodelled CO₂, pH, injury,
   reperfusion and thrombosis.

## What M8 could not determine

Whether any of this is survivable; whether the organ is constructible;
whether detection within 1.5 s is achievable; regional adequacy; the value of
compression kappa, which bounds the feasible region from both sides and is the
single most valuable missing datum in the study; and anything outside the
sampled domain, since every boundary reported is a boundary of what was tested.
