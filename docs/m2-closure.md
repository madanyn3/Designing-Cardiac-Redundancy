# M2 closure — circuit and specification gate

**2026-10-05 amendment:** Regional-flow parameterization is reopened following the user's evidence requirement. The frozen branch fractions lack source-to-value derivations. See the [provenance correction and published replacement candidates](m3-regional-flow-evidence.md). The historical specification freeze remains reproducible, but its allocations are not an evidence-qualified baseline for sizing.

2026-10-03. **M2 is closed at the roadmap's specification scope.** The user authorized progressing with the findings and prior inputs; remaining implementation choices were resolved as explicit exploratory assumptions, comparators or claim limitations. This does not certify the model against the full biological objective. No ODE simulation, pressure-feasibility result, organ-recovery prediction or optimized auxiliary size is claimed.

## Deliverables and task completion

| Task | Result and reviewable artifact |
|---|---|
| T2.1 Select reproducible baseline | User-accepted Smith structure with pinned source/provenance and unit/phase corrections; [source review](m2-source-review.md), [baseline configuration](../configs/baseline.yaml) |
| T2.2 Specify complete circuits | Native, dual support, selective support, oxygenator/unloading and separate auxiliary LV/RV realization; complete returns, preserved pools, kidney/coronary splits; [circuits](m2-circuits.md), [frozen overlay](../configs/m2-frozen-overlay.yaml) |
| T2.3 Mathematical specification | Source equations plus explicit project closure laws, domain limits, initialization and event semantics; [main specification](model-specification.md), [authoritative v1.0 supplement](m2-frozen-specification.md), [materialized hydraulic fixture](../data/m2-hydraulic-fixture.json) |
| T2.4 Failure scenarios | Full native LV/RV/septal activation loss, delayed/imbalanced support, limb restriction/drainage, no exchange, filling-domain failure and restoration interface; existing scenario list plus v1.0 event ordering |
| T2.5 Verification plan | Source comparisons, numerical tolerances, conservation, valve/pressure and closure-law tests; [verification plan](verification-plan.md) |

All five tasks are complete **as specifications/plans**, not implemented or validated models. Named support fixtures supply concrete initial volumes and actuator parameters for later interface tests; their values are not design conclusions. Gas-exchange constitutive laws and transport implementation stay in M4 as the roadmap states; M2 fixes their conservation interfaces and forbids implicit unlimited oxygen supply.

## Resolved choices

Adopt 5.25 L inventory, declared regional allocation/limb storage, open limb drainage, one-kidney split and lung-tissue return as A constructions. Add an ICP-plus-tone cerebral pressure-flow sensitivity, simple coronary compression sensitivity, fixed-R comparators and explicit conservative inlet-domain limits. Allocate auxiliary cavity blood from existing reservoirs before warm-up; never add blood at arrest. Freeze the published-evidence-informed pressure policy without calling it safe. [P01–P18 dispositions](m2-review-decisions.md) distinguish user direction from agent-selected numerical choices.

The aggregate brain model assumes regional distribution; it cannot prove that every region is adequately supplied. The simplified coronary closure omits capacitance/compression-driven displacement. Single-kidney filtration, hepatic/gut function, autonomous hypometabolic induction, preserved breathing drive, congestion injury limits, physical muscle geometry and recoverable duration are unresolved biology/model limitations or later work. They have not been “closed” as physiological facts.

## Verification evidence

Run from the project root:

```
python scripts/audit_m2_specification.py
python scripts/audit_m2_design.py
python scripts/audit_m2_freeze.py
```

The audits check source hashes, source algebra, volume inventory, graph continuity, resistance arithmetic, signed waterfall limiting cases, split conductances, auxiliary blood transfer, driver normalization and frozen-file hashes. [Freeze results](../data/m2-freeze-audit.json) and [snapshot manifest](../data/m2-freeze-manifest.json) record the static gate. A change to a frozen specification requires an explicit reviewed refreeze (`--freeze`); do not silently overwrite the manifest to hide a change.

## Handoff

### Note — progressive organ support and duration (D024)

The protected organ set is an experiment variable, rather than a fixed minimum set. Keep whole-brain coverage and the dependencies needed for the selected oxygenation method in every case; progressively add other organs and derive how the auxiliary-heart and circulatory-system requirements change at each support horizon. Also compare individual additions against the same core so the selected order does not determine the result. See the [organ-support study](organ-support-ladder.md) for the ladder, initial demand subtotals and evidence limits.

For each case, distinguish how long supplied organs meet modeled requirements, whether excluded organs remain recoverable, and whether their exclusion compromises the rest of the system. All tissues remain in blood-volume and oxygen-demand accounting. Preserving the possibility of return to the native heart requires its own maintenance/recovery assessment. More protected tissue increases the accounted protected demand, but does not establish a proportional increase in MAP or auxiliary-heart size.

This note supplements the M2 handoff without changing its frozen reference specification or closure status. M3 verifies added routing cases; M4 supplies oxygen/respiratory and organ constraints; M5 derives pump/circuit requirements; M7/M8 compare organ sets across duration and uncertainty. Current demand subtotals are not sizing or survival results.

**Next: M3.** Implement and reproduce the unchanged source configuration first. Then implement the regional/closure variants and diagnose their periodic initialization, conservation and model-domain behavior. Do not tune the project to hide a source reproduction failure. M4 closes gas/oxygen transport, extraction, respiratory work and finite demand transitions; M5 derives feasible pressure–flow/pump ranges including auxiliary self-demand. M6 handles detection/handover control; M7/M8 test uncertainty and conditional feasibility.

Reopen M2 if equations/initial conditions prove inconsistent, a required output cannot be represented at this fidelity, or sensitivities reverse a claimed conclusion. Numerical warm-up failure is a valid M3 finding and does not retroactively become evidence of physiological infeasibility. M0/M1 retain their evidence/assumption closures.
