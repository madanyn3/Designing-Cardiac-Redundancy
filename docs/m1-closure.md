# M1 completion and handoff

2026-10-02. **M1 is complete as an evidence/assumption audit.** The user explicitly authorized literature-informed assumptions and model limitations for this exploratory mathematical study. Completion does not imply validated physiology, a functioning simulator or demonstrated recovery.

## Completion gate

| Task | Final disposition and artifact |
|---|---|
| T1.1 Claim inventory | Complete: 62 stable original claims in the ledger; archive unchanged. |
| T1.2 Evidence audit | Complete at stated inspection depth: every claim retained with qualification, replaced or removed. `audit_resolution` makes final status explicit. Unsupported original assertions remain unclassified/inactive, rather than falsely relabeled E or A; replacement assumptions are labeled separately in the current documents. No uninspected full-text result is admitted as established. |
| T1.3 Parameter inventory | Complete: 210 rows with source/assumption/context or named deferred treatment. A reproducible native-tissue budget now supplies numerical starting scenarios. This inventory is not yet an executable circulation. |
| T1.4 Existing-support comparison | Complete at functional scope; no unsourced numerical device envelope is adopted. |
| T1.5 Presentation/reference repair | Complete through repair or explicit exclusion. Ambiguous historical author/year fragments are not active citations. Exact recovery of unused historical references is not a gate condition. |
| T1.6 Bibliography/search log | Complete for admitted M1 evidence at documented access depth. 35 records, with contextual/candidate records distinguished. Roadmap reading leads and new implementation sources are verified/imported when used in later milestones. |

The gate is met because every active numerical choice is traced to evidence or declared A, every original claim has a disposition, and every unresolved parameter has an owner and blocks a named dependent output. We have not claimed that all 210 quantities are measured or that every cited full text was accessible. This is a targeted evidence audit, not a systematic review.

## Numerical handoff

Read [assumption budget](m1-assumption-budget.md), [configuration](../data/m1-exploratory-assumptions.json), [derived results](../data/m1-demand-budget.json) and [per-parameter handoff](m1-parameter-handoff.csv).

The 75 kg synthetic male has a 255.85 mL O2/min resting reference budget. Total native-tissue demand under the two reduced-demand assumptions is 212.31 or 178.44 mL/min. Protected-set subtotals excluding modeled limb tissues are 197.43, 157.32 and 123.46 mL/min. These figures are conditional arithmetic, not pump-flow requirements, tissue viability bounds or demonstrated metabolic states. Biological pump self-demand is additional.

The registry contains 72 adopted assumptions, 22 conditional derived inputs, 16 contextual evidence candidates and 21 quarantined historical values. The other 79 entries are explicitly handed off: 36 model parameters, 16 design parameters, 25 unpredicted functional/recovery endpoints and 2 deferred limb outcomes. The first two groups must be populated before their dependent numerical experiments; the latter groups constrain what conclusions can be reported.

## Limitations to carry into the final report

1. The adult and organ partition are synthetic. Physiological compatibility is to be tested, not claimed from assembling plausible numbers.
2. Resting energy coefficients are literature-model transfers. Gut/lung allocation, caloric conversion and respiratory budget include assumptions.
3. Demand-reduction factors, timing and reversibility are hypothetical; maintenance, regional neurological function and organ-specific biochemical function are not independently identified.
4. The quiescent-heart case transfers animal arrest evidence to a cause-unspecified human pump-loss thought experiment. Other electrical/metabolic states can require more oxygen.
5. Diaphragm-supported breathing assumes preserved drive, airway and phrenic pathways. A simple ventilation budget does not establish sustained spontaneous breathing.
6. Tissue-store sensitivities include assumed saturation parameters and cross-muscle transfer. No-credit reference conservatively omits tissue buffering but retains blood/lung accounting; this does not prove a conservative injury prediction.
7. Oxygen deficit and pressure exposure do not predict recovery or safe arrest duration. Kidney filtration, liver biochemical output and gut barrier function are outside the first oxygen-balance claims.
8. Limb outcome exclusion is not physiological removal. Nonmuscle limb pools are now explicit assumed partitions; selective routing and reperfusion consequences require later work.

## Next milestone: M2

Select an accessible published LV/RV closed-loop model and inspect all equations/parameters. Specify complete whole-circulation and selective-routing candidates, native-lung and external-oxygenation alternatives, passive native paths, blood-volume/oxygen balances and failure intervention. Set baseline verification targets before implementation. Reconcile the M1 synthetic phenotype with that model explicitly. M3 implements and reproduces hydraulics; M4 adds demand/oxygen transport; M5 derives pump requirements. No implementation is claimed complete here.

## Verification

The budget script checks total mass, positive disjoint demand allocation and exact scenario coverage. Closure verification also checks parameter/claim ID uniqueness, source-key resolution, complete individual handoffs, local document links and the original archive hash. Numerical demand sensitivity checks are arithmetic checks only. Results are recorded in [progress](m1-progress.md).


Limb refinement: [adopted demand/perfusion amendment](m1-limb-amendment.md). Nonmuscle limb pools are explicit; demand reduction and reduced-flow targets are independent. M1 remains closed at its evidence/assumption scope.
