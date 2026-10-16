# Designing Cardiac Redundancy in Humans: Project Roadmap

**Active sequence (2026-10-15):** [M6 is closed at qualified control-scenario scope](m6-closure.md): 13 paired hydraulic scenarios, 12 paired metabolic replays and 254 tests. M7 uncertainty/sensitivity is underway, followed by M8 optimization. Preprint writing remains parked.

**M5 interpretation correction (2026-10-14):** The recorded auxiliary unmet oxygen-equivalent demand is glucose-limited in these runs, not demonstrated oxygen-delivery failure. Pressure failure remains unchanged. [Diagnosis and corrective experiments](m5-failure-diagnosis.md).

Created: 2026-09-30  
Status (2026-10-16): M0–M5 closed at their declared qualified gates. M3 retains R-M3-01; M4 retains L-M4-F01. M5 T5.1–T5.5 are complete at exploratory characterization/accounting scope: 235 tests, matched BIO/mechanical circulation and startup experiments. No viable auxiliary-heart design or repair window established. [M5 closure](m5-closure.md), [final comparison](../reports/m5-final-comparison/report.md). M6 is closed at qualified control-scenario scope; [closure](m6-closure.md), [results](../reports/m6-coupled/report.md). M7 is closed at bounded uncertainty/identifiability/structural-assessment scope ([closure](m7-closure.md), 291 tests). M8 is underway: [T8.1 protocol](m8-t81-protocol.md) is frozen at version 1 (347 tests); preprint writing remains parked.

## Research objective

**2026-10-03 study amendment:** [organ-support ladder](organ-support-ladder.md) makes protected organs and duration joint experiment variables. Carry it into M3 routing verification, M4 oxygen/function accounting, M5 pump/circuit sizing and M7/M8 uncertainty/experiments. M2's frozen reference remains unchanged; additional isolation/routing cases require explicit specification and verification. No stage demand subtotal is a pump specification or validated support time.

Determine under what explicitly stated physiological assumptions an auxiliary pump could maintain organ oxygen delivery after loss of native cardiac pumping, and identify the constraints that prevent it from doing so.

The project is a computational feasibility study motivated by a hypothetical biological redesign. Its success does not require a feasible backup heart. A reproducible negative result, a conditional feasibility boundary, or evidence that the question is underdetermined can each be a useful scientific result. A simulation meeting chosen constraints is not evidence of clinical safety or predicted survival.

The existing [primer](designing-cardiac-redundancy-in-humans-primer.md) is the starting hypothesis and eventual public-facing narrative. Its numerical design point is not a validated target.

## How to follow this plan

Work through M0–M10 in order. Read the listed material in order within each milestone; each entry says what to extract and which tasks it supports. Repeated links are intentional so the current milestone is self-contained. Complete the milestone's gate before relying on its outputs downstream. Record decisions and deviations instead of silently changing assumptions to obtain favorable results.

The references below are a curated starting syllabus, not a completed systematic review. Bibliographic records, abstracts, selected source text, or official documentation were inspected during planning; full-text evidence extraction remains a task. A listed source is not automatically verified support for every proposed use. Retrieve accessible full text through publisher, PubMed Central, library, or author routes. If unavailable, mark the evidence gap and do not invent page numbers, methods, parameter values, or conclusions. Follow references and newer citing studies when a model input requires evidence beyond this syllabus.

Use four labels throughout: **E** = source-supported evidence within its stated population and conditions; **A** = modeling assumption; **D** = result derived from equations or simulations; **S** = speculative biological capability. A paper about rodents, an isolated organ, or post-resuscitation care must retain that context.

## Milestone map

| Milestone | Main output | Depends on | Status |
|---|---|---|---|
| M0. Scope and research protocol | Questions, endpoints, exclusions, decision log | — | Complete: protocol gate, 2026-09-30 |
| M1. Evidence and primer audit | Claim ledger, source notes, parameter inventory | M0 | Complete: evidence/assumption gate, 2026-10-02 |
| M2. Circuit architectures and model specification | Complete circuit diagrams and equations | M1 | Complete: specification gate, 2026-10-03 |
| M3. Reproduce baseline haemodynamics | Verified baseline simulator | M2 | Complete with accepted source limitation R-M3-01, 2026-10-07 |
| M4. Oxygen transport and metabolic demand | Coupled supply–demand model | M3 | Closed at qualified exploratory scope, including fuel accounting and L-M4-F01, 2026-10-09 |
| M5. Auxiliary pump and energy budget | Architecture comparison and sizing bounds | M4 | Closed at qualified characterization/accounting scope; 235 tests, 2026-10-14 |
| M6. Detection, recovery, and failure modes | Controller and fault experiments | M5 | Closed at qualified control-scenario scope; 254 tests, 2026-10-15 |
| M7. Uncertainty and identifiability | Sensitivity rankings and uncertainty bounds | M6 | Closed at bounded assessment scope; [closure](m7-closure.md) |
| M8. Feasibility experiments | Frozen experiment suite and results | M7 | T8.1 frozen at version 1; T8.0 next |
| M9. Scientific synthesis | Technical report and completed primer | M8 | Not started |
| M10. Reproducible showcase | Release, figures, walkthrough, demo | M9 | Not started |

M0–M5 establish the core scientific model. M6–M8 make the research contribution. M9–M10 make it reviewable and presentable. Optional extensions appear after the core roadmap and must not delay a coherent first release.

## M0 — Define scope and a research protocol

**Objective:** Turn the thought experiment into falsifiable questions before writing a simulator.

### Reading, in order

1. [Guidelines for mechanistic modeling and analysis in cardiovascular research (2024)](https://pubmed.ncbi.nlm.nih.gov/38904851/). For T0.1–T0.4: extract the distinction between model purpose, assumptions, calibration, and model analysis.
2. [2025 AHA Guidelines, Part 11: Post–Cardiac Arrest Care](https://cpr.heart.org/en/resuscitation-science/cpr-and-ecc-guidelines/post-cardiac-arrest-care). For T0.2–T0.3: read blood pressure, oxygenation, ventilation, temperature control, and neuroprognostication sections. Record that these concern care after return of spontaneous circulation, not validation of an auxiliary-pump design.

### Tasks

- [x] **T0.1 — Write the research questions.** Primary: map the combinations of pump capacity, activation delay, oxygenation, and metabolic demand that satisfy explicit organ constraints. Secondary: determine dominant uncertainties, architecture tradeoffs, and failure cases.
- [x] **T0.2 — Define the modeled subject and event.** One idealized adult reference circulation with complete abrupt loss of both ventricles' active pumping; cause unspecified and no additional blood loss or vascular damage imposed. Explicit ventilation/gas-exchange dependency scenarios are defined in the protocol. Residual pumping and additional injury require separately labeled extensions; no claim to represent all clinical causes is made.
- [x] **T0.3 — Define endpoints.** Include organ flow, pressure, oxygen delivery/consumption, cumulative unmet modeled oxygen demand, and time to a declared constraint violation. Define each formula, units, threshold provenance, observation horizon, and handling of runs with no violation. Do not call these survival or brain-death predictions.
- [x] **T0.4 — Bound the first release.** Use a zero-dimensional closed-loop model; distinguish ideal pump performance from biological realizability. Exclude detailed electrophysiology, molecular reperfusion pathways, tissue engineering, and patient-specific clinical prediction from the core model.
- [x] **T0.5 — Start a decision log.** Record alternatives, reasons, evidence gaps, and what evidence would overturn a decision. Choose effort estimates only after a baseline reproduction attempt reveals the actual workload.

**Deliverables:** `docs/research-protocol.md`, `docs/decision-log.md`.

**Completion gate:** Every research question has measurable outputs; every endpoint has units and an interpretation; scope and exclusions are explicit; no success criterion assumes a particular pump size or hours of survival.

**Gate record (2026-09-30):** Completed in [research protocol](research-protocol.md) and [decision log](decision-log.md), supported by [issue-by-issue research findings](m0-research-findings.md) and [search record](m0-search-log.md). Numerical organ requirements and recoverability limits remain explicit M1–M4 evidence dependencies, not completed results. The protocol defines how missing thresholds produce an indeterminate judgment rather than an assumed pass.

## M1 — Audit the primer and build the evidence base

**Objective:** Separate supported physiology from extrapolation, unsupported claims, and proposed capabilities.

**Progress (2026-10-02):** M1 complete as an evidence/assumption audit. 62 claim dispositions, 210 parameter records, 35 bibliography records and a reproducible exploratory demand budget. Unknown quantities have explicit downstream owners or limit claims. See [closure](m1-closure.md) and [assumption budget](m1-assumption-budget.md).

### Reading, in order

1. [Kirino (1982), Delayed neuronal death in the gerbil hippocampus following ischemia](https://pubmed.ncbi.nlm.nih.gov/7093691/). For T1.2: extract species, ischemia method, observation times, and measured injury; assess exactly what it can support in the introduction.
2. [Sundgreen et al., Autoregulation of cerebral blood flow in patients resuscitated from cardiac arrest](https://pubmed.ncbi.nlm.nih.gov/11136927/). For T1.2–T1.3: extract measurement method and between-patient variation; avoid converting a pressure into a universal minimum cerebral flow.
3. [Richardson et al. (2021), ELSO interim ECPR guideline](https://pubmed.ncbi.nlm.nih.gov/33627592/). For T1.4: extract circulation/oxygenation functions, timing considerations, and practical limitations. This is a dated consensus document; check for updates before final synthesis.
4. [Yannopoulos et al. (2020), ARREST trial](https://pubmed.ncbi.nlm.nih.gov/33197396/), then [Suverein et al. (2023), INCEPTION trial](https://pubmed.ncbi.nlm.nih.gov/36720132/). For T1.4: compare eligibility, setting, intervention bundle, timing, endpoints, uncertainty, and limitations rather than pooling headline percentages.
5. [Induction of a torpor-like hypothermic and hypometabolic state in rodents by ultrasound (2023)](https://www.nature.com/articles/s42255-023-00804-z). For T1.2–T1.3: extract species, effect magnitude, induction time, duration, and reversibility; retain the rodent-to-human translation gap.

### Tasks

- [x] **T1.1 — Inventory every substantive claim.** Assign stable claim IDs to numerical, causal, comparative, and historical assertions in the primer. Include anatomy, timing, metabolism, evolution, pump geometry, and recovery claims—not only the sentences already carrying citations.
- [x] **T1.2 — Audit support.** Record exact supporting section/page/figure, population, experimental conditions, and limitations. Classify each claim E/A/D/S and mark retain, qualify, replace, remove, or unresolved. Correct the conflation of coronary autoregulation with pacemaking. Audit the universal injury timeline, evolutionary argument, spontaneous-recovery narrative, and sub-second trigger requirement.
- [x] **T1.3 — Build a parameter inventory.** Record symbol, units, nominal value if justified, range, context, source, source location, uncertainty type, and correlation constraints. Leave unknown values unknown; do not insert convenient defaults without an A label.
- [x] **T1.4 — Compare existing support.** Build a function table for CPR, VA-ECMO/ECPR, ventricular assistance, and the proposed biological pump. Distinguish native-lung dependence, oxygenation, flow support, activation assumptions, and shared failure modes. Obtain device-specific primary documentation if making quantitative device claims.
- [x] **T1.5 — Repair presentation separately from evidence.** Catalogue duplicated Layer 2, broken equations, spelling, incomplete citations, and inconsistent flow values. Preserve the original draft before substantive revision. Resolve incomplete historical references by title/DOI; do not infer an exact paper solely from an author-year pair.
- [x] **T1.6 — Build the bibliography and search log.** Record search terms, dates, inclusion reasons, contradictory evidence, inaccessible sources, and source versions. Add original studies for any parameter presently supported only by a secondary source.

**Deliverables:** `docs/claim-evidence-ledger.csv`, `docs/evidence-notes.md`, `docs/existing-support-comparison.md`, `data/parameters.csv`, `references/references.bib`, `docs/search-log.md`.

**Completion gate:** Every substantive primer claim has a disposition; every proposed model input has provenance or an explicit assumption/gap; no unsupported numerical claim is presented as established. Unresolved core inputs are bounded scenarios or block their dependent model component.

## M2 — Specify complete circuits and governing equations

**Objective:** Define systems that close blood-volume and oxygen balances before estimating their size.

**Progress (2026-10-03):** T2.1–T2.5 complete as specification/planning tasks. [Frozen v1.0](m2-frozen-specification.md) resolves implementation choices as explicit A assumptions and limitations. [M2 closure](m2-closure.md) records the static gate and M3 handoff; no dynamic simulation or organ-support validation has been run.

### Reading, in order

1. [Smith et al. (2004), Minimal haemodynamic system model including ventricular interaction and valve dynamics](https://pubmed.ncbi.nlm.nih.gov/15036180/). For T2.1–T2.3: obtain full equations, initial conditions, parameter definitions, and reported baseline outputs.
2. [Smith et al. (2006), Simulating transient ventricular interaction using a minimal cardiovascular system model](https://pubmed.ncbi.nlm.nih.gov/16400203/). For T2.2–T2.3: examine closed-loop coupling and ventricular interaction; decide what complexity the research question needs.
3. [Lorusso et al. (2021), ELSO interim VA-ECMO guideline](https://pubmed.ncbi.nlm.nih.gov/34339398/). For T2.2/T2.4: identify circulation and gas-exchange distinctions relevant to an extracorporeal comparator; check subsequent updates before using clinical statements in the final article.

### Tasks

- [x] **T2.1 — Choose a reproducible baseline.** Select a published model with sufficiently accessible equations and parameters. List any deviations and their consequences; if reproduction is impossible, select another documented model before implementation.
- [x] **T2.2 — Draw candidate circuits.** Specify systemic and pulmonary compartments, brain, coronary bed, remaining tissues, pump inlets/outlets, valves, and venous return. Compare baseline native circulation; a hypothetical auxiliary biventricular/native-lung circuit; and a separately labeled idealized pump-plus-oxygenator comparator. Add selective-organ perfusion only if its routing is fully specified.
- [x] **T2.3 — Write the mathematical specification.** Define states, flow laws, compliance/elastance relations, pressure references, valve behavior, initial conditions, blood-volume conservation, and event equations. Tabulate every unit conversion and map every parameter to M1.
- [x] **T2.4 — Specify failure scenarios.** State which paths remain open when native pumping stops. Check pulmonary transport, regurgitation/shunts, inlet filling, and volume accumulation. Treat obstruction and hemorrhage as distinct stress scenarios, not automatically solvable by another pump.
- [x] **T2.5 — Predefine verification targets.** Select source outputs for reproduction and tolerance criteria before tuning. Distinguish published numerical targets from broad physiological plausibility checks.

**Deliverables:** `docs/model-specification.md`, circuit figures, `configs/baseline.yaml`, `docs/verification-plan.md`.

**2026-10-02 design amendment:** Smith accepted by the user. [Auxiliary-heart design](m2-auxiliary-heart-design.md) adds pressure evidence, whole-brain/partial-organ objectives, limb-network resistance arithmetic and a chamber-based inverse-design interface. P13/P15 numerical assumptions are adopted as exploratory A choices in the v1.0 freeze, not established physiological minima. Actual output/chamber sizing remains M5 after M3/M4; no minimum viable SBP/DBP or physical septal mass is established by M2.

**Completion gate:** Each circuit traces a complete return path; no oxygenation or flow source is implicit; equations, units, and initial conditions are specified; baseline reproduction targets and tolerances are recorded.

## M3 — Implement and reproduce baseline haemodynamics

**Closed 2026-10-07:** user accepted R-M3-01, verified reconstruction with an explicit journal-reproduction limitation. [Closure report](m3-closure.md). The original failed screen remains in the audit; this acceptance supersedes pending-status notes below.

**Final review (2026-10-07):** technical work, anatomical-assumption disposition and reporting are complete. [Closure packet](m3-closure.md); 102 tests and the nine final storage-sensitivity runs pass. T3.3 retains the failed original thoracic-pressure screen. R-M3-01 is explicitly accepted by the user; M3 is closed with that limitation.

**User-directed sequence (2026-10-06):** the combined prescribed handover, delay/imbalance/depletion, E1/E2 and Smith-shaped auxiliary [verification task](m3-combined-circuits.md) is complete for its declared protocol: 102 tests and all 27 run comparisons pass (six runs correctly terminate at domain limits). That final documentation and audit task is now complete; the source discrepancy still requires the R-M3-01 closure decision. Whole M3 remains open until its completion gate is met.

**Latest progress (2026-10-06):** [Whole-network controls](m3-whole-network-controls.md) extend beyond limbs; 75 automated tests and all 15 new simulations pass. [Current status and next task](m3-current-status.md): pump resistance/capacity configuration is implemented ([guide](m3-backup-pump.md)); [work accounting](m3-pump-energy.md) is also implemented, with its 94-test checkpoint. The combined circuit verification now passes 102 tests and all 27 run comparisons. M3 remains open.

Execution plan: [M3 goals, sequence and M2 handoff](m3-plan.md). T3.1, T3.2 and T3.4 are complete for the declared implementation/verification scope; T3.3 remains partial because of the source-reproduction exception. Source numerics, physical inventory and the [regional resistor comparator](../reports/m3-regional/report.md) are verified. The evidence-amended regional split is closed by the user; [nominal pressure-aware closures](../reports/m3-pressure-regional/report.md) also pass numerical verification. Scheduled no-support arrest [F0 verification](../reports/m3-arrest/report.md) also passes. The basic B1/selective fixture also passes its [declared tests](../reports/m3-selective/report.md). The combined support/event/circuit task is verified; anatomical-assumption disposition and source-paper reconciliation/final reporting remain open.

**Objective:** Demonstrate that the numerical implementation solves the intended baseline model reliably.

### Reading, in order

1. [Smith et al. (2004)](https://pubmed.ncbi.nlm.nih.gov/15036180/). For T3.1/T3.3: work directly from the full model description selected in M2.
2. [SciPy `solve_ivp` documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html). For T3.1/T3.2: read methods, tolerances, event handling, and failure diagnostics; pin and record the version actually used.
3. [Ensuring reliability of safety-critical clinical applications of computational cardiac models](https://pubmed.ncbi.nlm.nih.gov/24376423/). For T3.2–T3.4: extract verification versus validation distinctions. Apply the credibility principles without implying clinical qualification.

### Tasks

- [x] **T3.1 — Implement modular equations.** Separate model equations, parameter loading, solver configuration, and plotting. Make the baseline runnable from a documented command; keep research logic out of notebook-only state.
- [x] **T3.2 — Verify conservation and limiting cases.** Check total blood volume, units, valve directionality, admissible compartment states, simple resistance/compliance cases with known solutions, and zero-pump behavior. Record failures; do not silently clip nonphysical values.
- [x] **T3.3 — Reproduce baseline outputs (closed with accepted limitation R-M3-01).** Generate relevant pressures, flows, chamber volumes, and pressure–volume loops. Quantify error against M2 targets; explain discrepancies instead of only judging visual similarity.
- [x] **T3.4 — Demonstrate numerical convergence.** Tighten tolerances and maximum step size, compare suitable solvers, and check event timing. Document warm-up to a periodic state and sensitivity to initial conditions.

**Deliverables:** `src/`, meaningful numerical tests in `tests/`, baseline configuration, `reports/baseline-reproduction.md`, generated baseline figures.

**Completion gate:** Predeclared reproduction and conservation tolerances pass; convergence changes key outputs less than the recorded acceptance tolerance; another run from the same configuration reproduces the result. Reproduction is labeled reproduction, not independent physiological validation.

## M4 — Add oxygen transport and finite metabolic transitions

**Started 2026-10-07:** [execution plan](m4-plan.md), [oxygen specification](m4-oxygen-specification.md), [foundation verification](../reports/m4-oxygen-foundation/report.md). T4.1 accounting is implemented at exploratory extraction/partition-assumption scope ([uptake](m4-organ-uptake.md)); its physiological parameter uncertainties remain. T4.2 finite demand transitions are implemented with the M1 hypothetical timing assumptions ([report](../reports/m4-demand-transitions/report.md)); T4.3 joint research-screen assessment is implemented ([report](../reports/m4-constraints/report.md)); regional/recovery thresholds remain unassessed. T4.4–T4.5 are now verified together at explicit assumption scope ([combined report](../reports/m4-respiratory/report.md)); finite liver glycogen/glucose accounting is now also verified ([fuel report](../reports/m4-fuel/report.md)). The [qualified closure](m4-closure.md) retains L-M4-F01 and the other physiological limitations. Pump efficiency/size optimization remains M5.

**Objective:** Determine whether modeled flow delivers adequate oxygen under explicit demand and gas-exchange assumptions.

### Reading, in order

1. [Identification of the critical oxygen delivery for anaerobic metabolism in critically ill septic and nonseptic humans](https://pubmed.ncbi.nlm.nih.gov/8411504/). For T4.1/T4.3: examine how delivery, consumption, and extraction were measured and why this population cannot supply universal organ thresholds.
2. [Absence of supply dependence of oxygen consumption in patients with septic shock](https://pubmed.ncbi.nlm.nih.gov/8305957/). For T4.3: inspect contrasting methods and potential coupling of calculated quantities; document why inferred thresholds may depend on measurement design.
3. [Rodent ultrasound-induced torpor study (2023)](https://www.nature.com/articles/s42255-023-00804-z). For T4.2: extract time courses and distinguish reduced demand from suppressed physiological function.
4. [Dankiewicz et al. (2021), TTM2 trial](https://pubmed.ncbi.nlm.nih.gov/34133859/). For T4.2/T4.4: distinguish temperature management after resuscitation from rapid metabolic suspension during arrest.
5. [Sundgreen et al., cerebral autoregulation after arrest](https://pubmed.ncbi.nlm.nih.gov/11136927/). For T4.3: revisit pressure/flow assumptions and their patient dependence.

### Tasks

- [x] **T4.1 — Specify oxygen accounting (implemented with explicit extraction/partition assumptions).** Implement arterial/venous content, mixing, compartment storage if retained, pulmonary replenishment, and organ uptake with consistent volume units. Source content-equation constants before use. Define extraction limits and the response when supply is insufficient.
- [x] **T4.2 — Separate organ demand functions (implemented as prescribed hypothetical scenarios).** Include brain, native myocardium, auxiliary myocardium when biological, and remaining tissues; resolve renal demand separately if kidney claims remain. Implement finite induction and recovery times. Include no suppression and several labeled hypothetical suppression depths, including the primer's 5% residual-demand scenario.
- [x] **T4.3 — Define constraints and deficit proxies (research screens; unknown biological limits unassessed).** Distinguish arterial pressure, perfusion pressure, delivered oxygen, and extractable oxygen. Define accumulated unmet modeled demand mathematically. Source any physiological thresholds; otherwise present threshold sensitivity rather than a claimed viability cutoff.
- [x] **T4.4 — Test ventilation scenarios (implemented/verified with explicit gas/respiratory assumptions).** Compare maintained gas exchange, reduced gas exchange, and absent replenishment. Avoid an unlimited arterial oxygen source when simulating apnea. Explain which CO2, acid–base, temperature, substrate, and microvascular processes are omitted and how they limit duration claims.
- [x] **T4.5 — Verify the coupling (numerical verification complete; independent matched physiological validation unavailable).** Test oxygen conservation with sources/sinks accounted for, mixing cases, zero uptake, absent gas exchange, and normal baseline behavior. Seek independent data for baseline oxygen outputs; mark unavailable validation explicitly.

**Deliverables:** Oxygen-model specification, source-backed parameter updates, coupled simulator, oxygen verification report, example supply–demand plots.

**Completion gate:** Blood and oxygen balances pass; finite oxygen stores behave consistently without replenishment; no suppression case works as a baseline; all deficit metrics are labeled proxies; no human 95% suppression capability is asserted.

## M5 — Implement the auxiliary pump and bound its requirements

**Current user scope:** [mandatory organ set](m5-supported-organ-set.md) includes native heart, whole brain, liver, respiratory muscles, one lung and one kidney, plus auxiliary self-supply. Earlier aggregate-lung/four-group runs do not implement this amendment. [T5.3 closure](m5-t53-closure.md) audits historical sizing only; anatomical/pump design remains unresolved.

**2026-10-14 objective-alignment amendment:** [Tunable organ preservation](m5-organ-preservation.md) carries the active-support and reduced preservation demands together. First verify a fixed-hydraulic demand sweep; then jointly explore conductance, pump workload and delivery. The earlier 3.8% whole-body change describes a cardiac-only comparator, not the selective design. Final sizing/architecture conclusions require this alignment.

**Started 2026-10-09:** [execution plan](m5-plan.md), [first bounded-stage characterization](../reports/m5-bounded-pump/report.md), [working requirements](../reports/pump-requirements.md). [Dynamic paired-stage tests](../reports/m5-dynamic/report.md), independent BIO controls and joint metabolic screens are verified. [T5.1 closure](m5-t51-closure.md), [results](../reports/m5-t51/report.md). T5.2 is closed at qualified accounting scope with BIO as the primary design: [closure](m5-t52-closure.md).

**Objective:** Estimate conditional hydraulic and energetic requirements without treating an ideal source as a realizable organ.

### Reading, in order

1. [Smith et al. (2004)](https://pubmed.ncbi.nlm.nih.gov/15036180/). For T5.1: use pressure–volume and valve relationships consistently with the native model.
2. [Suga et al. (1981), Regression of cardiac oxygen consumption on ventricular pressure-volume area in dog](https://pubmed.ncbi.nlm.nih.gov/7212076/). For T5.2: examine energetic relationships and experimental conditions; do not transfer canine coefficients directly to a hypothetical human organ.
3. [ELSO VA-ECMO guideline (2021)](https://pubmed.ncbi.nlm.nih.gov/34339398/). For T5.1/T5.4: revisit inlet/outlet and native-circulation interactions for the mechanical comparator.

### Tasks

- [x] **T5.1 — Add bounded pump models.** Start with an explicitly idealized capacity-limited pump, then a pressure/flow or elastance model. Include inlet/preload dependence and right/left circuit balance. Identify where an idealized source grants capabilities a biological pump does not have.
- [x] **T5.2 — Compute energy and self-demand.** Integrate hydraulic power across actual pressure differences for both circuits. Separate hydraulic output from metabolic/input energy. Represent uncertain efficiency, basal maintenance, and auxiliary coronary demand transparently; use sensitivity ranges when evidence is absent.
- [x] **T5.3 — Reproduce and audit the primer calculation.** Recalculate the 740 mL/min, 46 g, and 0.11 W estimates with clean equations and units. Compare thin-wall and exact spherical-shell volumes, assess the wall-thickness/radius ratio, test fixed-point convergence, and list excluded anatomical mass. Do not claim a corrected organ mass from geometry alone.
- [x] **T5.4 — Compare BIO with an abstract mechanical reference.** User-narrowed scope: BIO is primary; no materials/device-design search. Use the same native failure, documented initial inventory, organ demands, pressure/flow/pulsatility endpoints and oxygenation capability. Keep external mechanical energy separate from biological oxygen/glucose demand. A boundary-accounting benchmark alone does not close the matched closed-circuit comparison.
- [x] **T5.5 — Biological trigger, activation and readiness costs.** Use primary biosystem evidence to motivate separate detection/confirmation, signal, recruitment and activation-ramp parameters; test startup delay rather than assuming immediate full support. Retain dormancy, maintenance, perfusion while inactive, activation energy, stasis/thrombosis and anatomy as explicit costs/limitations. Measured reflex latency is not automatically a dormant-heart startup time.

**Deliverables:** Pump configurations, `reports/pump-requirements.md`, reproducible sizing calculation, architecture comparison figures.

**Completion gate:** Self-demand enters the budget when applicable; hydraulic and metabolic power are distinct; no unlimited-flow source is presented as physical; the original design point is reproducible and its limitations quantified or explicitly unresolved.

## M6 — Model detection, support, recovery, and shared failures

**Objective:** Test whether a controller improves modeled support under delays, noisy observations, and partial recovery. User scope: biological mechanisms motivate the policy; a biologically realized controller is not required. [Specification](m6-controller-specification.md), [synthetic benchmark](../reports/m6-policy-benchmark/report.md).

### Reading, in order

1. [Ã…ström and Murray, Feedback Systems, author-hosted second-edition resource](https://www.cds.caltech.edu/~murray/FBS/Second_Edition.html). For T6.1–T6.3: read the system modeling, stability, and feedback material relevant to the chosen simple controller; follow the authors' current-site link as needed.
2. [2025 AHA post–cardiac arrest guideline](https://cpr.heart.org/en/resuscitation-science/cpr-and-ecc-guidelines/post-cardiac-arrest-care). For T6.3: identify recovery concerns and limits of transferring post-resuscitation targets into the hypothetical circuit.
3. [ELSO ECPR guideline (2021)](https://pubmed.ncbi.nlm.nih.gov/33627592/). For T6.4: extract dependencies and complications to motivate fault cases, not to assign unsupported fault probabilities.

### Tasks

- [x] **T6.1 — Specify controller states.** Define monitoring, activation, support, recovery, and fault states with observables, thresholds, hysteresis, delays, and transition rules. Separate sensing latency from pump ramp-up and metabolic induction.
- [x] **T6.2 — Compare activation policies.** Use a simple threshold/delay policy as the baseline. Test noisy pressure/flow observations, false triggers, missed triggers, and loss of sensor input. Use identical scenarios for comparisons.
- [x] **T6.3 — Implement recovery as a scenario.** Ramp native pumping and metabolic demand independently, test re-arrest and pump competition, and compare abrupt versus graded handover. Do not infer prevention of reperfusion injury unless an independently supported injury model is added.
- [x] **T6.4 — Construct a failure-mode table.** Include absent ventilation, inadequate venous return, obstruction, depleted volume, loss of controller function, and common dependence on oxygen supply. Mark each modeled, bounded, excluded, or requiring a different architecture. Do not assign reliability percentages without data.
- [x] **T6.5 — Test controller invariants.** Check transition determinism, no rapid state chattering, actuator bounds, and declared responses to invalid states. Report both successful and adverse cases.

**Deliverables:** Controller specification and code, `docs/failure-modes.md`, recovery and fault experiment figures.

**Completion gate:** Baseline and controlled cases are directly comparable; all transitions and failure responses are specified; common-cause failures appear in the results; recovery claims remain within modeled physiology.

## M7 — Quantify uncertainty, sensitivity, and identifiability

**Objective:** Establish whether conclusions are robust or artifacts of selected parameters and model structure.

### Reading, in order

1. [Guidelines for mechanistic modeling and analysis in cardiovascular research (2024)](https://pubmed.ncbi.nlm.nih.gov/38904851/). For T7.1–T7.3: focus on calibration, parameter identifiability, and uncertainty.
2. [Herman and Usher (2017), SALib](https://joss.theoj.org/papers/10.21105/joss.00097). For T7.2: follow the linked official software documentation for the selected screening/global sensitivity method.
3. [Practical Identifiability and Uncertainty Quantification of a Pulsatile Cardiovascular Model](https://arxiv.org/abs/1710.07989). For T7.3: use the published [2018 journal version](https://doi.org/10.1016/j.mbs.2018.07.001) as a methodological example; its identifiable subset is not transferable to this model.
4. [Ensuring reliability of computational cardiac models](https://pubmed.ncbi.nlm.nih.gov/24376423/). For T7.4: revisit the boundary between numerical verification and predictive validity.

### Tasks

- [x] **T7.1 — Classify uncertainty.** See [uncertainty registry and protocol](uncertainty-plan.md). Separate variability between hypothetical subjects, uncertain physiology, numerical error, model-form choices, and speculative capabilities. Justify bounds and distributions. Use bounded scenarios where probability distributions lack evidence.
- [x] **T7.2 — Run sensitivity analysis.** [Qualified closure](m7-t72-closure.md): initial screening, factorial interaction checks and [84 joint global cases](../reports/m7-global/report.md), 274 tests, numerical and sample-size stability checks, bootstrap rank uncertainty. Closed at bounded six-factor/ten-second Morris screening scope; no continuous-domain Sobol or physiological convergence claim. Start with inexpensive screening, then apply global analysis to influential inputs. Track sample size, seeds, failed simulations, estimator convergence, and uncertainty in rankings. Do not apply independent-input assumptions to correlated parameters without justification.
- [x] **T7.3 — Assess identifiability.** Identify which inputs baseline observations can constrain. Use synthetic recovery experiments or profile-based analysis where appropriate. Do not fit many parameters to a few baseline numbers and call them uniquely estimated. Closed at conditional scope: [closure](m7-t73-closure.md), 26 local cases and seven nonlinear recovery/noise/compensation experiments.
- [x] **T7.4 — Test structural alternatives.** Compare plausible autoregulation, oxygen-extraction, demand-transition, and pump formulations. Separate a robust result from one that depends on a chosen equation form. Closed at bounded assessment scope: [structural closure](m7-t74-t75-closure.md), six paired extraction/transition replays, three paired regulation circulation/replays, nine finite-reserve diagnostics and prior verified pump comparison. Full tissue buffering and supply-limited contraction remain explicit unimplemented limitations.
- [x] **T7.5 — Set the next evidence priority.** Rank missing measurements by how much they alter the feasibility conclusion. Explain which new data would most improve the project. [Final evidence priorities](m7-evidence-priorities.md) identify decision impact, required observations, identifiability and horizon dependence.

**Deliverables:** `docs/uncertainty-plan.md`, sensitivity and identifiability reports, uncertainty figures, reproducible sampling configurations.

**Completion gate:** Dominant assumptions and unstable rankings are reported; sampling adequacy is assessed; fitted parameters have identifiability caveats; scenario fractions are not labeled probabilities of patient survival.

## M8 — Run the frozen feasibility experiments

**Objective:** Produce the central research results with traceable comparisons and a declared analysis plan.

### Reading, in order

1. Revisit [ARREST](https://pubmed.ncbi.nlm.nih.gov/33197396/) and [INCEPTION](https://pubmed.ncbi.nlm.nih.gov/36720132/). For T8.1/T8.4: use their population/endpoint differences to discipline interpretation; do not calibrate a hypothetical circulation to trial survival percentages.
2. Revisit [mechanistic modeling guidelines](https://pubmed.ncbi.nlm.nih.gov/38904851/). For T8.1–T8.4: check reporting of assumptions, calibration, and analysis against the planned study.

### Tasks

- [x] **T8.1 — Freeze the experiment protocol.** Save the code version, parameter provenance, scenario matrix, endpoint definitions, numerical settings, and analysis choices before the main runs. Record amendments and rerun affected analyses. Frozen at version 1: [protocol note](m8-t81-protocol.md), [`configs/m8-protocol.json`](../configs/m8-protocol.json), design/scenario registry and tested endpoint evaluator; 347 tests.
- [ ] **T8.0 — Optimize the bounded design.** Use M6 policies and M7 uncertainty/sensitivity to search pump/circulation settings across size/energy budgets. Compare analytical screening and a reproducible constrained stochastic search against a matched-budget sampling baseline; report best-found trade-offs, seed variability, boundary dependence and refined full-model verification. Do not optimize uncertain physiological demands freely or claim a global/clinical optimum. Freeze this search protocol under T8.1 before main runs.
- [ ] **T8.2 — Run the core matrix.** Include no auxiliary support; auxiliary support without suppression; support with finite hypothetical suppression; ventilation/gas-exchange degradation; delayed activation; partial native recovery; and selected shared failures. Compare architectures only under explicit capability differences.
- [ ] **T8.3 — Generate primary results.** Plot capacity versus activation-delay feasibility maps, organ oxygen-deficit time series, hydraulic/energy tradeoffs, uncertainty bands, and sensitivity rankings. Define exactly what a colored region means. Report solver failures separately from physiological constraint violations.
- [ ] **T8.4 — Seek disconfirmation.** Search deliberately for cases where the proposed design fails and where architecture rankings reverse. Report boundaries outside the sampled domain as unknown. Separate interpolation from extrapolation.
- [ ] **T8.5 — Answer the research questions.** State what is feasible within the model, what depends on speculative assumptions, and what cannot be decided. Report horizons with no violation as such, not as indefinite support.

**Deliverables:** Frozen run manifests, result tables, final scientific figures, `reports/feasibility-results.md`.

**Completion gate:** Every main conclusion maps to a saved experiment; controls and adverse results are included; all figures regenerate; the central claim remains conditional on assumptions and validation limits.

## M9 — Complete the technical report and the five-layer primer

**Objective:** Turn the evidence and computation into two coherent accounts for different readers.

### Reading, in order

1. [2025 AHA post–cardiac arrest guideline](https://cpr.heart.org/en/resuscitation-science/cpr-and-ecc-guidelines/post-cardiac-arrest-care), [Kirino](https://pubmed.ncbi.nlm.nih.gov/7093691/), and [TTM2](https://pubmed.ncbi.nlm.nih.gov/34133859/). For T9.2/T9.3: recheck clinical, injury-timing, and temperature claims against their actual scope.
2. [Rodent torpor study](https://www.nature.com/articles/s42255-023-00804-z), [ARREST](https://pubmed.ncbi.nlm.nih.gov/33197396/), and [INCEPTION](https://pubmed.ncbi.nlm.nih.gov/36720132/). For T9.2/T9.3: recheck extrapolation and balance in the speculative and comparative sections.
3. [Mechanistic modeling guidelines](https://pubmed.ncbi.nlm.nih.gov/38904851/). For T9.1/T9.4: use as a reporting review alongside the project's own verification plan.

### Tasks

- [ ] **T9.1 — Write the technical report.** Include research questions, existing approaches, model equations, provenance, numerical verification, available validation, experiments, uncertainty, results, and limitations. Separate findings from proposed future capabilities.
- [ ] **T9.2 — Finish the primer's five layers.** Layer 1: detection and activation tradeoffs. Layer 2: evidence versus hypothetical metabolic demand reduction. Layer 3: complete circuit and pump requirements. Layer 4: whole-brain protection and limits of deficit proxies. Layer 5: modeled handover versus unresolved biochemical reperfusion injury. Add methods/results context before conclusions.
- [ ] **T9.3 — Audit every final claim again.** Link claims to source locations or reproducible derivations; qualify species, conditions, and uncertainty. Remove unsupported evolutionary narratives unless a separate primary-literature review justifies keeping them. Refresh the literature search and record its date.
- [ ] **T9.4 — Perform an adversarial review.** Check omitted organ demand, impossible oxygen sources, double-counted flow, unsupported thresholds, geometry approximations, model dependence, and overstatement. Seek physiology and numerical-modeling review if accessible; record whether external review actually occurred.
- [ ] **T9.5 — Resolve presentation defects.** Repair equations, references, duplicated sections, terminology, captions, and figure accessibility. Use appropriate precision and show assumptions beside numerical results.

**Deliverables:** `docs/technical-report.md`, completed primer, updated bibliography/ledger, `docs/review-log.md`.

**Completion gate:** All central claims have evidence or derivations with limits; numerical claims match saved outputs; both documents agree; review issues are resolved or explicitly recorded as limitations. External review is not claimed if unavailable.

## M10 — Package a reproducible computational science showcase

**Objective:** Make the contribution easy to inspect, reproduce, and explain.

### Reading, in order

1. [Mechanistic modeling guidelines](https://pubmed.ncbi.nlm.nih.gov/38904851/). For T10.1/T10.2: revisit reproducibility and reporting expectations.
2. [SciPy solver documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html) and the [SALib paper/software links](https://joss.theoj.org/papers/10.21105/joss.00097). For T10.1: ensure numerical and sensitivity settings are captured for the actual dependency versions.

### Tasks

- [ ] **T10.1 — Build the reproduction path.** Pin dependencies, document installation, provide a quick baseline command and a full figure-generation command, capture seeds/configurations/software versions, and report expected runtime from measurement.
- [ ] **T10.2 — Reproduce from a clean environment.** Follow only the written instructions; regenerate headline results and compare against declared tolerances. Separate a smoke run from the full study. Record unavailable data or platform restrictions.
- [ ] **T10.3 — Prepare the portfolio narrative.** Explain the question, personal contribution, mathematical choices, verification, one principal result, strongest limitation, and what would change the conclusion. Include a short walkthrough and a one-page project summary.
- [ ] **T10.4 — Add an optional interactive explorer.** Let readers vary pump capacity, delay, and suppression within the studied domain. Display assumptions and model validity limits near outputs. Identify interpolation and prohibit presenting the explorer as a patient calculator.
- [ ] **T10.5 — Prepare a release.** Review code/data licenses, source acknowledgments, repository structure, and figure provenance. Create a release checklist and archive manifest. Public hosting or repository publication is a separate action, not implied by completing local documentation.

**Deliverables:** Root `README.md`, environment specification, runnable experiments, reproduction report, portfolio summary, release checklist, optional demo.

**Completion gate:** A fresh environment reproduces headline outputs; documentation explains both contribution and limitations; no hidden notebook state is required; links and figures work; the final artifact does not imply demonstrated biological feasibility or clinical effectiveness.

## Optional extensions after the core study

| Extension | Research purpose and prerequisite | First reading | Completion criterion |
|---|---|---|---|
| Bayesian calibration | Estimate uncertainty from suitable independent observations; requires identifiable parameters and licensed data after M7 | [Fast and robust parameter estimation with uncertainty quantification for cardiac function](https://pubmed.ncbi.nlm.nih.gov/36773593/) | Prior sensitivity, calibration diagnostics, held-out assessment, and limits of patient generalization reported |
| Multi-fidelity vascular model | Test whether spatial flow detail changes a specific conclusion; requires a demonstrated limitation of the zero-dimensional model | [Multilevel and multifidelity uncertainty quantification for cardiovascular hemodynamics](https://pubmed.ncbi.nlm.nih.gov/32336811/) | Documented comparison with the baseline and measured accuracy/cost tradeoff |
| Data-informed control | Optimize a declared support objective under actuator and observation constraints; requires M6–M8 | [Ã…ström and Murray, Feedback Systems](https://www.cds.caltech.edu/~murray/FBS/Second_Edition.html) | Beats the simple controller on held-out scenarios without hiding adverse cases |
| Model-based experimental prioritization | Determine which uncertain physiological measurement most changes the feasibility boundary; requires M7 | [Mechanistic modeling guidelines](https://pubmed.ncbi.nlm.nih.gov/38904851/) | Ranked evidence needs supported by quantitative changes in conclusions |

Do not add molecular models, CFD, machine learning, or detailed anatomy solely for sophistication. Each extension must answer a question the core model cannot and have a separate verification plan.

## Evidence and parameter record templates

**Claim ledger fields:** `claim_id`, `document_location`, `claim_text`, `classification_E_A_D_S`, `source_id`, `source_url`, `source_location`, `population_or_species`, `conditions`, `support_summary`, `limitations`, `disposition`, `verification_status`, `last_checked`.

**Parameter fields:** `parameter_id`, `symbol`, `description`, `units`, `nominal_value`, `lower_bound`, `upper_bound`, `distribution_if_justified`, `correlations`, `source_id`, `source_location`, `population_and_state`, `evidence_or_assumption`, `used_by`, `uncertainty_notes`.

**Reading note fields:** citation/version; question being answered; study design; population; measured versus inferred quantities; relevant methods and results; exact source locations; limitations; project claims/parameters affected; full-text verification status.

**Experiment manifest fields:** experiment ID; question; code version; configuration; parameter-table version; solver settings; seeds; endpoints; run status; output paths; numerical checks; interpretation limits.

## Project-wide definition of done

- The research questions are answered conditionally, including negative and unresolved findings.
- Core numerical results have verification evidence, provenance, and quantified uncertainty where feasible.
- All final factual claims have appropriately scoped references; assumptions and speculative capabilities remain visible.
- No model proxy is silently converted into clinical survival, neurological recovery, or brain-death prediction.
- The article, technical report, code, and figures tell the same story and can be reproduced from documented inputs.
- Deferred questions and optional extensions are explicit, so completion does not depend on solving every aspect of hypothetical human redesign.


Limb refinement: [adopted demand/perfusion amendment](m1-limb-amendment.md). Nonmuscle limb pools are explicit; demand reduction and reduced-flow targets are independent. M1 remains closed at its evidence/assumption scope.
