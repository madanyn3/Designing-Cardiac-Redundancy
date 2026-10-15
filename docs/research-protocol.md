# M0 research protocol

Version 1.4 — 2026-10-03. User-directed amendments: [R01–R12 outcomes](m1-review-outcomes.md), D019–D020 below, [problem/model-fit review](m2-problem-model-fit-review.md), and [D024 organ-set study](organ-support-ladder.md). Status: M0 scope/protocol gate complete; M1 evidence/assumption audit is complete; M2 specification gate closed 2026-10-03; model suitability remains conditional on M3/M4 verification and the declared limits. This is a computational feasibility protocol, not a clinical protocol.

Companion documents: [research findings and sources](m0-research-findings.md), [decision log](decision-log.md), [search record](m0-search-log.md), [roadmap](project-roadmap.md).

Labels: **E** published evidence within the stated experimental context; **A** study choice or modeling assumption; **D** mathematical/model-derived result; **S** speculative capability. A referenced paper does not turn a proposed study choice into evidence.

## 1. Problem statement

**D024 scope amendment:** the protected organ set is now an experiment variable crossed with support duration. Whole-brain coverage and the dependencies of the selected oxygenation method remain compulsory; add other organs individually and in nested sets. Distinguish protected tissues, actual perfusion routes and whole-body accounting. The previous protected set remains a reference, not a proven minimum. A brain/oxygenation-only comparator does not establish native-heart or whole-body recovery. The [organ-support study](organ-support-ladder.md) defines the ladder, duration endpoints, evidence limits and downstream work; its core demand subtotal is a proxy pending complete respiratory-dependency accounting.

Determine the feasible operating envelope of an abstract backup circulatory system, including selectively perfused configurations, following complete abrupt loss of left- and right-ventricular pumping in an adult initially at rest. Derive the flow, pressure, oxygenation, activation-time, and energy requirements needed to satisfy explicit organ-level support constraints over declared observation horizons. Evaluate how these requirements change with circuit routing, metabolic demand, and physiological uncertainty. Preserve near-complete neurological and other organ recovery as the intended biological objective, while reporting clearly when the model cannot establish recoverability.

**Current pressure and brain requirements:** explicitly derive and evaluate systemic MAP, systolic and diastolic pressure as summaries of the same central arterial waveform, alongside pulmonary and regional perfusion pressures. Seek the lowest feasible modeled systemic pressure for each declared demand/routing scenario and report its associated flow and energy tradeoffs. All brain regions remain protected recipients of perfusion; reduced demand and lower-than-resting flow are allowed subject to supply constraints. The aggregate brain compartment assumes distribution and cannot establish every region's adequacy. The [consolidated statement and requirement matrix](m2-problem-model-fit-review.md) govern the current M2 review; earlier full-resting-flow interpretations are superseded.

**A:** Neither pump technology, pump count, anatomical size, a target flow, nor a number of hours of successful support is stipulated. Failure to find a feasible region, or inability to identify a recovery boundary from available evidence, is a valid result.

## 2. Falsifiable research questions

| ID | Question | Measurable outputs / test |
|---|---|---|
| RQ1 | Which combinations of pump capacity, pressure head, activation delay, oxygenation and demand satisfy the declared constraints? | Feasibility maps, first violation time, limiting organ and constraint, organ pressure/flow/oxygen time series |
| RQ2 | Can selective perfusion reduce required pumping without violating constraints in protected or less-perfused tissues? | Matched whole-circulation versus routed-circulation comparisons; total flow, hydraulic energy, every compartment's oxygen deficit |
| RQ3 | How much does hypothesized demand reduction alter requirements, including induction and reversal? | No-suppression control; demand multiplier and transition-time sensitivity; maximum delay and conditional support horizon |
| RQ4 | What prevents continued support or handover to the native ventricles? | Oxygen-source exhaustion, inadequate head/flow, inlet limitations, volume pooling, constraint crossings during recovery |
| RQ5 | Which evidence gaps or model choices dominate the answer? | Sensitivity rankings and changes in feasibility boundaries across parameter and model-form alternatives |
| RQ6 | What can the available injury evidence establish about recovery? | Source-to-endpoint compatibility matrix; validated recovery endpoint if available, otherwise an explicit unidentifiable recovery boundary |

Model purpose, independent validation and uncertainty must remain distinct; this protocol adopts the workflow in [Colebank et al. (2024)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11442102/) (E). The simulator is not authorized to answer RQ6 through an oxygen-deficit cutoff invented for convenience.

## 3. Reference subject and event

**A — reference case:** one idealized healthy male in his twenties of athletic build at rest, initially supine, with a stable baseline circulation. Body/organ masses, blood volume, hemoglobin, baseline output and vascular parameters will come from a compatible source-backed reference configuration in M1–M3. No arbitrary mixing of unrelated population means into a purported patient.

At t = 0, both ventricles lose active pumping. Cause is unspecified. No additional hemorrhage, vascular obstruction or structural circulatory injury is imposed. This is a scope boundary, not a claim that all clinical arrests have these properties. Consequences of absent pumping are retained.

**A — mechanical interpretation:** remove active ventricular pressure generation; retain passive pressure–volume behavior, valves, vascular compliance and stored blood. Do not force every vessel flow to zero at t = 0 or erase ventricular volume. M2 must define the intervention in the selected equations. Arrhythmia and residual contractility are outside the reference event.

**A — respiratory scenarios:**

1. Required design scenario: diaphragm-supported native-lung exchange, with adequate respiratory drive, phrenic function, airway patency, ventilation/CO2 clearance and breathing work explicitly accounted for. Fractional diaphragm capacity is acceptable only if these requirements are met. Externally maintained ventilation remains a comparison condition; it does not demonstrate spontaneous breathing.
2. Reduced exchange, and absent oxygen replenishment with finite initial stores: dependency tests, not additional causal diagnoses.
3. External gas exchange: a separate architecture capability, with its energy and capacity requirements reported.

Baseline temperature is held at 37 °C as an A reference setting; suppression/temperature alternatives are separate scenarios. Initial posture is supine; height-dependent variants are optional sensitivity cases. Source-backed ranges and numerical respiratory parameters must be frozen before relevant experiments.

## 4. Model scope and organ policy

**A:** use a closed-loop 0D pressure–volume network with explicit LV and RV, systemic and pulmonary vessels, valves, oxygen transport and organ branches. [Smith et al. (2004)](https://pubmed.ncbi.nlm.nih.gov/15036180/) is the user-accepted hydraulic foundation, not yet numerically reproduced here. Oxygen transport, regional perfusion and respiratory coupling require extensions. Atria are lumped into inflow reservoirs in the selected baseline. [Model-fit review](m2-problem-model-fit-review.md) distinguishes supported structure from missing capabilities.

| Compartment | Intended preservation objective (A) | First-release measurable representation | Additional evidence needed before functional claims |
|---|---|---|---|
| Whole brain | Preserve capacity for near-complete neurological recovery; consciousness during support is not required | Cerebral flow, effective perfusion pressure, demand/uptake, storage and deficit | Regional vulnerability, suppression mechanism, delayed outcomes, validated linkage to cognition and function |
| Native myocardium | Preserve potential for restoration of pumping | Coronary supply and demand, including maintenance when mechanical work is absent | Human nonworking-heart demand and recovery relationship |
| Kidneys | Minimize required ongoing function; one functioning kidney is an exploratory candidate, not established sufficiency. Retain tissue maintenance/accounting | Renal pressure, flow, extraction and oxygen balance | Filtration/transport workload, regional oxygenation and injury evidence |
| Liver | Preserve recoverable metabolic/synthetic function | Separate arterial and portal inflow oxygen accounting | Compatible demand and pressure data; functional recovery evidence |
| Gut/splanchnic bed | Preserve tissue/barrier recoverability | Flow, oxygen uptake/deficit and portal outflow | Duration/severity-dependent barrier injury and recovery |
| Lungs | Account for gas exchange and tissue demand | Pulmonary flow, oxygen source, exchange capacity | Tissue oxygen budget, ventilation/CO2 dependence and regional limits |
| Remaining body | No tissue silently discarded | Residual volume, demand, flow and deficit; split further when necessary | Consequences of reduced perfusion and reperfusion; respiratory muscles if spontaneous ventilation is claimed |
| Backup tissue, if biological | Account for self-maintenance and work | Additional oxygen-demand term | Efficiency, basal demand and activation cost |

This is an accounting set, not an evidence-based ranking of expendable organs. No tissue has yet been assigned a validated zero-perfusion duration in this project; this does not mean every tissue requires uninterrupted flow. Temporary zero or reduced flow in resting limb skeletal muscle is an explicit research candidate, subject to duration, temperature, metabolic reserves, reperfusion and delayed functional outcomes. Human tourniquet metabolic observations motivate this candidate but do not validate whole-body exclusion during arrest [Haljamäe and Enger (1975)](https://pmc.ncbi.nlm.nih.gov/articles/PMC1343870/). See [skeletal-muscle follow-up](skeletal-muscle-perfusion.md). Selective circulation is an experimental routing variable. Hypothetical isolation/control capabilities are S; conservation and their hydraulic consequences are D.

**A:** liver and gut cannot be represented as independent parallel beds when portal supply is retained. Downstream oxygen content must reflect upstream extraction; the same blood must not be counted twice as new systemic flow. Other serial pathways require the same treatment.

## 5. Pressure and spatial detail

**D — initial network relation:** for a passive resistive path, Q = (P_in − P_out)/R. With elevation, use hydraulic potential P + rho*g*z in consistent SI units. Compliance adds storage; inertia is retained when required by the baseline. This provides location-dependent pressures without a 3D mesh.

**A:** record both upstream and downstream pressures. Brain pressure requires an explicit external-pressure/venous-outflow closure; P_arterial − max(P_venous, ICP) may be tested as a simplified effective gradient, not treated as a universal cerebral-flow law. Autoregulation, CO2 and venous pressure remain uncertainty axes. Human post-arrest variation motivates this policy [Sundgreen et al. (2001)](https://pubmed.ncbi.nlm.nih.gov/11136927/) (E).

Upgrade to more compartments or 1D modeling if regional pressure, wave transmission, collapsible vessels or distribution materially changes a conclusion and the 0D closure cannot be validated. 3D modeling is reserved for local geometry-dependent questions. M2 must document a comparison criterion before an upgrade; complexity alone is not a deliverable.

## 6. Endpoint dictionary

Time is in minutes for oxygen integration; simulator seconds must be converted explicitly. Blood flow Q is L blood/min; content C is mL O2/L blood; demand M and uptake U are mL O2/min. All D equations below are definitions/conservation statements or stated model closures, not clinical thresholds.

| Endpoint | Formula / definition | Units and interpretation |
|---|---|---|
| Systemic and pulmonary throughput | Flow across separately named circuit cross-sections | L/min; never add serial pulmonary and systemic flow as cardiac output |
| Organ inflow | Q_i = sum of incoming branch flows | L/min; also report signed path flows and reversals |
| Local pressures | P_in,i, P_out,i, P_external,i; chosen effective gradient | mmHg; distinguish arterial pressure from perfusion gradient |
| Central MAP/SBP/DBP | Cycle integral of P_AO divided by cycle duration; cycle maximum; cycle minimum | mmHg; three summaries of one waveform, not independently free inputs. For irregular/nonpulsatile runs use declared 1 s windows (A) and label window extrema/mean. Report arrest transients separately. |
| Oxygen delivery | DO2_i = sum_in Q*C | mL O2/min arriving at organ; multiple inlet contents handled separately |
| Actual oxygen uptake | U_i from coupled blood/tissue balance | mL O2/min consumed, not automatically equal to demanded metabolism |
| Tissue oxygen store | dS_i/dt = J_i − U_i, J_i = net blood-to-tissue transfer | S_i in mL O2; optional closure requires evidence-backed capacity |
| Blood-compartment oxygen balance | dB_i/dt = sum_in Q*C − sum_out Q*C − J_i | mL O2/min; at steady storage U_i equals inflow minus outflow oxygen |
| Extraction fraction | J_i/DO2_i where DO2_i > 0 | Dimensionless; undefined at zero delivery; do not silently set to zero |
| Instantaneous unmet demand | delta_i = max(M_i − U_i, 0) | mL O2/min; M_i remains distinct from supply-limited consumption |
| Cumulative unmet demand | D_i(t) = integral from 0 to t of delta_i dt | mL O2; a monotone exposure proxy, not a measured injury dose |
| Demand-normalized deficit | D_i / M_i,baseline | min of baseline-demand equivalent; not ischemic tolerance time |
| Pressure exposure | Integral max(P_min,i − P_eff,i, 0) dt | mmHg*min; exploratory exposure proxy |
| Hydraulic power | sum_pumps (P_out − P_in)*Q in Pa and m³/s | W; account for every powered stage |
| Input energy | Integral (P_hyd/eta + P_aux) dt in seconds | J; efficiency/auxiliary costs are A unless sourced; eta > 0 |
| First constraint violation | tau = inf{t in [0,H]: any declared g_j(t) < 0} | min; report organ, constraint and initial violations |

Reverse flow must use upwind oxygen content and direction-aware transfer; the simple positive-inflow notation is not permission to subtract oxygen using the wrong donor content. Storage-free limits are allowed only as documented approximations. Deficits do not disappear when perfusion returns; any separate repayment variable must not erase prior exposure.

### Threshold provenance and run classification

Every constraint record must contain: ID, organ, quantity, formula, units, comparator, numeric value or range, E/A status, source and source location, population/state, applicability limitation, instantaneous versus averaged evaluation, averaging window, allowed excursion duration, and uncertainty treatment.

**A reference screen:** systemic MAP >= 65 mmHg may be used as a clinical-context comparator only. The source recommendation is for adults after ROSC, not for an auxiliary circulation or a universal viability threshold [AHA 2025, section 5](https://cpr.heart.org/en/resuscitation-science/cpr-and-ecc-guidelines/post-cardiac-arrest-care) (E). Never label this comparator sufficient organ protection.

**2026-10-03 proposed pressure contract:** [source extraction and evaluation rules](m2-pressure-requirements.md) add A exploration screens SBP>=90 and DBP>=50 mmHg alongside MAP>=65, with stricter alternatives. These transfers include surgical radial-pressure evidence and cannot certify central or regional viability. Nominal MAP 80/100–70 waveform seed remains A. Screens operate on declared cycles/windows with no allowed violating window in the strict research classification; source exposure durations are not injury allowances. User review of the numerical policy remains open.

**A strict supply screens:** report M_i − U_i <= epsilon_i for each declared protected compartment, and separately report the same diagnostic for every compartment. Limb recovery constraints are deferred in the first feasibility claim; limb supply deficits remain reported. For both screens, epsilon_i is a declared numerical tolerance justified by convergence tests, not physiological tolerance. Report storage depletion and all pressure constraints separately. This conservative screen may reject delayed support; do not relax it silently.

Organ pressure minima, extraction ceilings and recovery-linked deficit limits remain **unresolved numerical inputs**. Explicit assumed store-accounting floors and demand scenarios are now adopted in the M1 assumption budget; they are not injury thresholds. Before M4/M8, either obtain compatible evidence, declare and sweep explicit A ranges, or mark the dependent feasibility judgment indeterminate. An unspecified threshold cannot be counted as a pass. No numerical deficit limit is adopted as a human recovery limit in M0.

For a complete declared constraint set: report pass within model, violation, numerical failure, or indeterminate. Record sensitivity to any finite excursion allowance. If no violation occurs, report tau > H (right-censored), not infinite support. If the simulation fails, do not substitute its last time for physiological failure.

### Observation and recovery windows

**A computational horizons:** H = 1, 5, 15, 30 and 60 minutes from pump loss; H = 120 minutes is an explicitly exploratory extension. These are sampling windows, not expected survivable durations. Extensions require reassessing omitted substrate, CO2, acid–base, thermal and injury processes. Record onset, activation and cumulative exposure before activation.

Native recovery is a prescribed scenario at t_restore, not an assumed spontaneous event. A 30-minute post-handover simulation window is an A hydraulic-stability check, not a recovery assessment. Continue reporting pre-handover deficits.

**A biological objective:** eventual near-baseline cognition, memory, motor/sensory and autonomic function, with preservation of the declared protected organs as the intended objective. Limb recovery is excluded from the first feasibility claim by the 2026-10-02 user decision. This does not establish whole-body recovery; systemic effects of limb underperfusion and reperfusion remain limitations or explicit model dependencies. Candidate clinical evidence windows are 90 days and later follow-up; the model does not predict these scores. The BOX trial measured neurological outcomes including cognition at three months, providing an endpoint example rather than a validation bridge [Kjaergaard et al. (2022)](https://pubmed.ncbi.nlm.nih.gov/36027564/) (E). A numerical margin for “near-complete” recovery remains a separate evidence/measurement-design task; no proxy substitutes for it.

## 7. Abstract pump and reduced circulation experiments

**A architecture families:** whole-circulation pumping through native lungs; selective routing through native lungs; circulation plus external gas exchange. Device count remains open. Specify inflow/outflow sites, blood return, passive native paths, all gas exchange and every added control capability before comparing configurations.

Characterize capacity by a pressure–flow envelope, inlet/preload dependence, activation delay, ramp time, controllability, efficiency and energy availability. A prescribed unlimited-flow source is only a mathematical control. Routing must respect blood volume, venous return and oxygen conservation. Passive pulmonary transit is neither guaranteed nor excluded: its required pressure balance must be solved.

For fixed architecture and scenario, map feasible combinations rather than report independently minimized head, flow and energy as if one device could realize all minima simultaneously. Present Pareto tradeoffs and limiting constraints. A necessary steady-state oxygen bound is Q_i >= M_i/(C_in,i*E_max,i) for a single-inlet, storage-free bed; it is not sufficient when pressure, heterogeneous supply, multiple inlets or evolving stores matter (D).

Core controls: normal baseline, no backup, immediate ideal support, delayed bounded support, no suppression, finite hypothetical suppression, reduced/absent oxygen replenishment and prescribed native recovery. Every selective-routing case must report losses in less-perfused compartments. Hemorrhage and obstruction are outside this event; add only as separately labeled future scope.

## 8. Evidence, validation and exclusions

Prefer direct human organ measurements and primary experiments for numeric inputs; retain animal species, anesthesia, temperature, surgical setting and reperfusion follow-up. Abstract-only sources may motivate scope but cannot supply uninspected equations or precise source locations. Record contradictory results. Do not convert resting energy expenditure into minimum viable demand.

M1 establishes parameter provenance and remaining gaps. M2 chooses accessible equations. M3 verifies baseline mechanics and numerical conservation. M4 verifies oxygen accounting against independent compatible data where available. M7 tests uncertainty and identifiability. M8 freezes thresholds, scenario ranges, code and horizon before the main experiments. Successful baseline reproduction does not validate arrest or injury prediction.

Core exclusions: detailed electrophysiology, cause-specific arrest treatment, molecular reperfusion simulation, tissue engineering, actual regional brain shutdown, device implantation design, patient-specific prediction and numerical human survival/recovery prediction without independent validation. Reduced demand is a scenario assumption; biological induction is S unless demonstrated under compatible conditions.

## 9. M0 completion audit and handoff

| Roadmap item | Completed artifact |
|---|---|
| T0.1 questions | Section 2 maps every question to measurable outputs |
| T0.2 subject/event | Section 3 defines reference event, passive mechanics and respiratory alternatives |
| T0.3 endpoints | Section 6 specifies equations, units, threshold provenance policy, horizons and censoring |
| T0.4 scope | Sections 4–8 define organs, fidelity, pump abstraction and exclusions |
| T0.5 decisions | [Decision log](decision-log.md) includes alternatives and reopening criteria |

The M0 gate is met as a scope/research-protocol milestone. This does not close organ physiology or numerical feasibility. Required downstream work is enumerated issue-by-issue in the findings document. No success criterion assumes a pump size or support duration. Effort estimates remain deferred until baseline reproduction.


## 2026-10-02 scope amendment

[R01–R12 outcomes](m1-review-outcomes.md) govern subsequent work. Whole-brain low-activity scenarios retain maintenance and essential respiratory-control requirements. Native myocardium has zero required pumping work but an unresolved quiescent maintenance demand. Kidney/liver functional fractions are exploratory; gut barrier and residual tissue maintenance cannot be silently deleted. Limb recovery is deferred, while limb volume, routing, oxygen exposure and systemic consequences remain accounted for. A scaled native heart is an initial pump candidate only: technology, size, count and scale remain open. Resting demand is a comparison baseline, not a proven physiological upper bound. Numerical derating requires a recorded source-transfer or assumption decision.


## Version 1.2: exploratory numerical assumptions

The user authorized educated assumptions to make the thought experiment quantitatively testable. [M1 assumption budget](m1-assumption-budget.md) is authoritative for the synthetic reference adult, resting demands, organ-specific total-demand multipliers, finite transitions, respiratory budget and store sensitivities. These supersede prior no-numeric-adoption wording. [M1 closure](m1-closure.md) lists limitations and deferred inputs; no physiological minima or human recovery times are certified. Scenario bounds are hypotheses, not statistical confidence limits.


## 2026-10-02 limb amendment (D017)

The user approved separate modest nonmuscle limb-demand reduction and reduced-perfusion scenarios. [Limb amendment](m1-limb-amendment.md) specifies assumed skin/bone/other/adipose pools and independent demand/flow factors. These pools join limb muscle outside the first protected recovery set; circulating blood remains conserved. Peripheral limb recovery is not inferred from whole-brain preservation. Earlier muscle-only numerical subtotals are superseded. Unmet demand and systemic/reperfusion limitations remain explicit.

## 2026-10-02 auxiliary-design amendment (D019)

**Historical interpretation below superseded by D020, 2026-10-03:** full brain perfusion means all regions, not full resting flow or mandatory unsuppressed demand.

Smith is accepted. The main inverse-design scenario now targets full baseline whole-brain flow with unsuppressed modeled brain demand, while exploring one protected kidney and fractional liver, nutritive lung and respiratory-muscle perfusion. Earlier reduced-brain-demand scenarios remain alternatives. Native heart, gut, residual body, the second kidney and limbs remain in mass/oxygen accounting; fractional flow is not automatic demand reduction. The [design amendment](m2-auxiliary-heart-design.md) proposes research pressure settings and documents why no universal viable systolic/diastolic pair can be claimed. Pump sizing remains an output of the coupled model.

## 2026-10-03 clarification (D020)

The user requires perfusion of every brain region, while seeking reduced overall brain demand. Remove the full-resting-flow constraint and restore M1 1.0/0.8/0.6 demand scenarios; 0.8 is the first exploratory reduced case, not a proven minimum. All regions remain protected and no regional zero-flow strategy is introduced. Regional distribution is assumed, not verified, by the aggregate brain compartment. Finite induction, reversibility and preserved respiratory drive are explicit dependencies. Jointly minimize feasible systemic MAP with pulmonary exchange, ventilation, diaphragm work, CO2 and organ supply constraints; increased lung/diaphragm flow does not imply proportional oxygenation. [Evidence and mathematical clarification](m2-brain-respiratory-coupling.md).


## 2026-10-03 M2 specification handoff (D023)

[M2 closure](m2-closure.md) adopts the frozen circuit/initialization/pressure-flow specification as exploratory A choices. This closes a specification gate, not the full biological objective. The [v1.0 supplement](m2-frozen-specification.md) governs current cerebral/coronary sensitivity laws, kidney split, auxiliary allocation and domain limits; older missing-closure wording is historical. No regional recovery or safe pressure minimum is established.
