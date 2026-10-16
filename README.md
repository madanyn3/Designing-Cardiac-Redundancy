# Designing Cardiac Redundancy in Humans

**Latest status (2026-10-16):** [M7 is closed at bounded assessment scope](docs/m7-closure.md). M8 has begun: the [T8.1 experiment protocol](docs/m8-t81-protocol.md) is frozen at version 1 with a two-column design/scenario registry, a tested endpoint evaluator and a measured search budget, at **347 passing regression tests**. [T7.4/T7.5 closure](docs/m7-t74-t75-closure.md) adds conservative regulation comparisons and finite-reserve diagnostics; [evidence priorities](docs/m7-evidence-priorities.md) are finalized. No viable biological pump or repair window is established; supply-limited contraction and full tissue buffering remain limitations. Preprint writing remains parked. Earlier checkpoints below are historical.

**M5 interpretation correction (2026-10-14):** The recorded auxiliary unmet oxygen-equivalent demand is glucose-limited in these runs, not demonstrated oxygen-delivery failure. Pressure failure remains unchanged. [Diagnosis and corrective experiments](docs/m5-failure-diagnosis.md).

## Current status — 2026-10-14

M5 is closed at qualified exploratory scope: **235 tests pass**. [M5 closure](docs/m5-closure.md), [final comparison and figure](reports/m5-final-comparison/report.md), [implementation guide](docs/m5-bilateral-implementation.md). Mandatory-organ routing, separate supported/preserved lungs, a matched mechanical reference and BIO startup delays are implemented. All final pressure screens fail and BIO self-supply remains inadequate; no viable design or repair window is established. M6 is next. Earlier checkpoint entries below are historical.

A computational feasibility study of hypothetical auxiliary circulation after loss of native cardiac pumping, motivated by a biological thought experiment.

The central question is: **Under what physiological assumptions could an auxiliary pump maintain organ oxygen delivery, and which constraints make this impossible?**

## Start here

- [Detailed M5–M7 report](docs/m5-m7-report.md): auxiliary-heart design and energy, selected-organ support, glycogen accounting, activation/handover, sensitivity, identifiability, structural uncertainty and the M8 handoff. Companion to the M0–M4 report; distinguishes verified implementation from unresolved biological feasibility.

- [Detailed M0–M4 report](docs/m0-m4-report.md): problem definition, evidence, assumptions, equations, milestone findings and remaining work. [Combined respiratory verification](reports/m4-respiratory/report.md), [normal oxygen reference](reports/m4-baseline/report.md).

- [M4 plan and progress](docs/m4-plan.md), [oxygen specification](docs/m4-oxygen-specification.md), and [first verification report](reports/m4-oxygen-foundation/report.md): oxygen transport, finite lung stores and [organ uptake](docs/m4-organ-uptake.md) implemented; [finite demand transitions](docs/m4-demand-transitions.md) also implemented; [pressure/oxygen assessment](docs/m4-constraints.md) also implemented; respiratory coupling and baseline verification are complete at exploratory scope; finite fuel accounting is verified at explicit restrictive scope. [Qualified M4 closure](docs/m4-closure.md), [fuel report](reports/m4-fuel/report.md).

- [M3 final report and closure decision](docs/m3-closure.md), [consolidated verification](reports/baseline-reproduction.md), and [anatomical assumptions](docs/m3-anatomical-assumptions.md). M3 is closed with user-accepted source-reproduction limitation R-M3-01.

- [Combined handover, stress and circuit findings](docs/m3-combined-circuits.md), [27-run verification report](reports/m3-circuits/report.md): prescribed restart/withdrawal, delayed and unequal support, domain limits, E1/E2 and a scaled auxiliary-heart fixture.

- [Configurable backup pump](docs/m3-backup-pump.md): independent right/left resistance and capacity; defaults match the prior fixture. [Hydraulic energy accounting](docs/m3-pump-energy.md) is implemented; [Handover and circuit verification](docs/m3-combined-circuits.md) is also complete for the declared protocol.

- [Current M3 status and next task](docs/m3-current-status.md), [whole-network controls](docs/m3-whole-network-controls.md), and [15-run verification report](reports/m3-whole-network/report.md): independent organ controls now extend beyond limbs.

- [Repair-window review of the selective pump](docs/m3-selective-implementation.md): what the verified fixture does, where it differs from the four-group support objective, and the next implementation contract.

- [Selective routing and backup-pump tests](reports/m3-selective/report.md), [implementation guide](docs/m3-selective-implementation.md): configurable regional storage, restriction/trapping/leakage, reopening and conserved blood-origin transport. Numerical fixtures; no recovery validation.

- [Selective perfusion and reperfusion research](docs/m3-selective-perfusion-reperfusion.md): arterial/venous controls, retained blood, blood-cell and anticoagulation dependencies, and decisions before selective backup support.

- [Scheduled native-arrest report](reports/m3-arrest/report.md) and [code guide](docs/m3-arrest-code-guide.md): no-support F0 transients, continuous-state preservation and 15 verified phase/solver runs.

- [Pressure-aware regional report](reports/m3-pressure-regional/report.md), [sensitivity results](reports/m3-pressure-sensitivity/report.md), and [code guide](docs/m3-pressure-closure-guide.md): cerebral closing pressure and coronary compression, verified separately from the fixed-resistance comparator.

- [Evidence-amended regional baseline](reports/m3-evidence-regional/report.md): separate digestive branches and referenced muscle/skin/adipose/respiratory estimates. Bone and limb-other are the two remaining provisional flow allocations; numerical run status is recorded in the report.

- [M3 code guide](docs/m3-code-guide.md): detailed walkthrough of equations, state arrays, valve events, numerical integration, outputs and tests; source files also contain explanatory comments/docstrings.
- [Numerical methods and model knobs](docs/m3-model-knobs.md): Radau explained, all current tuning controls, and commands for your own experiments.
- [M3 sweep report](reports/m3-sweeps/report.md): periodic baseline, 19 parameter/numerical cases, figures, archived failures and the unresolved thoracic-pressure reproduction discrepancy.

- [Project roadmap](docs/project-roadmap.md): ordered milestones M0–M10, detailed tasks, objectives, reading references, deliverables, and completion gates.
- [M0 research protocol](docs/research-protocol.md): problem statement, research questions, organ scope, endpoint definitions and completion audit.
- [M0 evidence findings](docs/m0-research-findings.md): six issues, published sources, inspection limits and downstream research gaps.
- [Decision log](docs/decision-log.md): scope choices, alternatives and conditions for revision.
- [M1 progress](docs/m1-progress.md): active evidence audit, claim ledger and next readings.
- [Working primer](docs/designing-cardiac-redundancy-in-humans-primer.md): revised research narrative, with original claims preserved in the archive.
- [M1 review packet](docs/m1-review-decisions.md): unresolved issues, evidence and proposed treatments for review.

## Current status

**T5.4/T5.5 started:** [BIO versus mechanical energy benchmark](reports/m5-architecture-foundation/report.md) and [biological trigger/activation evidence](docs/m5-t55-activation.md). 228 tests pass. The revised lung representation, matched closed-circuit runs and full startup sweeps remain open.

**T5.3 closed:** [Historical sizing audit and closure](docs/m5-t53-closure.md), [reproducible results](reports/m5-primer-audit/report.md). 217 tests pass. The new [mandatory organ set](docs/m5-supported-organ-set.md) is recorded; one-lung routing and updated mandatory-organ experiments remain to be implemented. M5 T5.4–T5.5 remain open.

**Joint exploration:** [Four-case BIO/conductance report](reports/m5-joint-preservation/report.md) verifies new hydraulics, recalculated self-demand and paired metabolic replays. Three cases pass the final-window pressure screen; all retain auxiliary coronary deficits. Coronary-specific delivery is the next design bottleneck.

**Selective-support amendment:** [Tunable organ demand profiles](docs/m5-organ-preservation.md) and [fixed-hydraulic sweep](reports/m5-organ-preservation/report.md) address the cardiac-only comparator gap. 206 tests pass. Preservation targets are explicit assumptions; joint pump/conductance optimization and recovery validation remain open.

**2026-10-14:** T5.2 closed at qualified energy-accounting scope; 198 tests pass. [Closure and next task](docs/m5-t52-closure.md), [native-preservation findings](reports/m5-native-preservation/report.md). M5 remains open; no feasible support duration established. The entries below retain earlier milestone checkpoints.

**M0-M2 are complete at their scope/evidence/specification gates. M3 is closed at verified-reconstruction scope with accepted limitation R-M3-01; the original reproduction screen remains failed.** The M3 closure suite had 102 passing tests. All declared project comparisons, including nine final storage-sensitivity runs, pass; all nine frozen M2 files are unchanged. The source baseline and independent reconstruction agree closely, but the thoracic-pressure output reduction remains 19.078% versus the original 9 +/- 5-point target. [Closure report](docs/m3-closure.md), [final audit](reports/m3-final/report.md). The user accepted the exception on 2026-10-07. M4 is closed at qualified exploratory accounting scope: 162 tests at the fuel checkpoint, 17 respiratory scenarios and 16 fuel scenarios with paired verification. Oxygen, glucose, glycogen, export and actual substrate-dependent oxygen use are explicitly accounted for. L-M4-F01 retains the missing tissue-glucose-buffer and joint physiological baseline limitation. [M4 closure](docs/m4-closure.md), [fuel report](reports/m4-fuel/report.md). The current full suite passes **191 tests**. M5 T5.1 is closed at bounded-model characterization scope: 1,500 static points, ten dynamic stage scenarios, four independent BIO scenarios and six joint oxygen/glucose screens are verified. All tested joint cases retain coronary demand deficits. [T5.1 results](reports/m5-t51/report.md), [closure and equations](docs/m5-t51-closure.md), [M5 plan](docs/m5-plan.md). T5.2 energy and self-demand has started, with the Smith-like two-ventricle BIO as the selected design. [Energy plan](docs/m5-t52-energy.md), [mechanical first increment](reports/m5-biological-energy/report.md), [normal-heart comparison and demand sensitivities](docs/m5-t52-normal-comparison.md). M5 remains open. Physiological recovery and an optimized fuel-limited repair window remain unestablished.

Latest study amendment: [organ set versus support duration](docs/organ-support-ladder.md). Start from brain and oxygenation dependencies, add protected organs progressively, and derive the changing pump/circuit requirements. The document distinguishes native-tissue demand subtotals from uncomputed sizing and recovery results.

Read the [M1 closure and handoff](docs/m1-closure.md) and [numerical assumption budget](docs/m1-assumption-budget.md). A synthetic 75 kg male provides the reference; demand reduction and oxygen-store alternatives are explicit assumptions. Unknown hydraulic/design inputs and unpredicted recovery endpoints remain visible. Pump requirements and support duration are later outputs.


Limb refinement: [adopted demand/perfusion amendment](docs/m1-limb-amendment.md). Nonmuscle limb pools are explicit; demand reduction and reduced-flow targets are independent. M1 remains closed at its evidence/assumption scope.


M2: [closure report](docs/m2-closure.md), [frozen specification](docs/m2-frozen-specification.md), [decision dispositions](docs/m2-review-decisions.md). Static audits pass; source baseline reproduction and hydraulic simulation remain M3 work, oxygen/respiratory transport M4, and auxiliary sizing M5.

- [M3 physical inventory verification](reports/m3-physical-inventory/report.md): first project-extension step, 5,250 mL accounting, source-coordinate identity and six new tests.

- [Regional circulation code guide](docs/m3-regional-code-guide.md): state layout, branch flows, hepatic balance and commands.
