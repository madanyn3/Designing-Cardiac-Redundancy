# Results

Draft for the preprint, written against saved runs. Every number carries the
file it came from, in `docs/claim-evidence-ledger.csv` under the M8 claim
range. Conditions common to all results: compression kappa 0.5, reference
physiology, prescribed 0.5 s activation, a 40 s settled horizon with a plateau
check, and the declared 3.0 mmHg pressure margin. The nine unresolved
requirements of §5.2 remain attached, so a numerically clean run is
*indeterminate*, not *satisfied*. No statement here concerns survival,
neurological outcome or clinical applicability.

## 1. A feasible design exists, and the earlier negative result was an artefact

Three independent searches — boundary tracing with a root-find on the signed
pressure margin, a matched-budget Sobol sample, and a joint compass search over
all three design controls — locate the same region of the design space
(Table 1). The two unconstrained searches agree to within 0.001 mL O₂/min while
differing materially in every control variable, which identifies the optimum as
a **shallow valley rather than a point**: among the joint search's eight best
designs, coronary conductance spans 1.388 to 1.808 while the objective varies
by 0.17 per cent.

Boundary tracing returned 7.812 mL O₂/min, 2.7 per cent above the valley,
because it traced boundaries only at 60, 80 and 100 beats/min while the optimum
lies near 49 to 57. The structural argument was sound and the rate axis was
sampled too coarsely. We report this as a finding about the method.

This supersedes an earlier negative result from the same codebase, which found
no configuration of 84 reaching the pressure screen. That earlier sweep held
cavity scale at a quarter-scale test fixture and varied conductance and
contractility instead. Cavity scale was always a legal parameter and inventory
never constrained it: every value in 0.1 to 1.0 is admissible on the declared
initial state. Three further quantities were frozen at the same fixture and had
to be coupled to the design before optimization was meaningful — the auxiliary
mass ratio, the activation rate charged in demand, and the auxiliary oxygen
store. **The negative result characterised an unexamined constant, not the
design space.**

## 2. Two designs are reportable, and they differ in what they assume

The unconstrained optimum is not reportable: it holds 0.005 mmHg of pressure
margin (§5.2). Under the declared 3.0 mmHg minimum, 36 designs satisfy both the
pressure screen and auxiliary self-supply. Robustness is inexpensive — the
margin costs **7.1 per cent** more auxiliary oxygen — and the choice that
remains is not about cost but about assumption (Table 2).

The cheaper design requires coronary conductance 1.748, inside an exploratory
0.5 to 2.0 bracket that rests on no measurement. The alternative sits at
conductance 0.953, at or below nominal, and therefore depends on that bracket
not at all; it costs 17.9 per cent more and needs a 40 per cent larger cavity,
and it holds nearly three times the required margin. Every conductance band in
the bracket contains a feasible margin-holding design, so **the bracket buys
optimality, not feasibility** — the materially weaker dependency.

Margin falls monotonically as conductance rises, which identifies the trade the
bracket actually represents: a larger coronary bed lets a smaller cavity feed
its own myocardium, and a smaller cavity holds less pressure margin. The two
reportable designs are the two ends of that trade.

## 3. Detection, not the pump, is the binding constraint

Holding the design fixed across nine declared operating cases, one satisfies
every screen: prescribed activation. **Five of nine satisfy every organ screen
including auxiliary self-supply**, and in four of those the auxiliary feeds its
own myocardium, every organ group is satisfied, and mean pressure sits well
above the enforced threshold — 68.69, 68.69, 68.69 and 104.35 mmHg — yet the
case fails, because release at 2.50 s puts five windows outside the declared
two-second startup allowance.

**In those four cases nothing about the design fails.** What fails is the
arrival time of support. An earlier version of this matrix was run at a design
holding +2.26 mmHg, below the declared margin, where those same cases also
failed the pump group and left room to read the design as simply inadequate.
That reading does not survive the re-run.

Measured against release time, the cost-optimal design's pressure margin runs
+3.10, +2.94, −0.67, −9.95, −15.95 and −21.00 mmHg at 0.5, 1.0, 1.5, 2.0, 2.5
and 3.0 s. Beyond 1.5 s the margin falls at roughly 12 mmHg per second, as a
late release both loses windows and accrues an arterial volume debt that must
be repaid.

Interpolating between solved points gives two thresholds per design — where the
screen fails, and where the declared margin is lost — and they separate the two
reportable designs sharply:

| Design | Margin at 0.5 s | Screen satisfied to | Margin held to |
|---|---|---|---|
| scale 0.536, conductance 1.748 | +3.10 | 1.41 s | **0.81 s** |
| scale 0.749, conductance 0.953 | +8.93 | 1.71 s | **1.57 s** |
| scale 0.650 (comparator) | +14.64 | 1.70 s | 1.58 s |

**Capacity buys latency tolerance, and it nearly doubles it.** The 18 per cent
more expensive design that escapes the conductance assumption (§2) also holds
its margin to 1.57 s rather than 0.81 s. The two advantages coincide, which
makes the trade in §2 a better bargain than cost alone suggests.

The requirement handed back to detection is therefore a detection-plus-
recruitment budget of **0.8 to 1.6 seconds depending on which design is
built**, against an implemented controller at 2.55 s. Every value lies below
it. The study converts a modelling assumption into a design requirement; it
cannot say whether that requirement is achievable.

## 4. Supported duration is set by physiology, not by the design

The primary endpoint is the time for which vital-organ perfusion is sustained
(Table 3). At reference physiology, charging all oxidative demand to glucose,
it is **33.8 min**; across the declared scenario range it spans 22.5 to
95.2 min. A 60 s coupled replay sustained zero unmet demand in all three organ
groups and measured hepatic export at exactly the effective cap, giving a drain
consistent with the analytic figure.

Three properties of this endpoint matter more than its value.

**The export rate binds, not the store.** Of the release law's four terms,
mobilisation capacity times liver function binds at 0.105 mmol/min. The
glycogen store term is 19,786 and hepatic delivery 0.539 at a measured
539 mL/min of hepatic outflow: the circulation delivers five times more hepatic
flow than export can use, and the store would outlast the circulating pool by
two orders of magnitude.

**The design barely moves it.** Elasticities are −1.03 to demand and +0.15 to
export, and the pump is 7 to 16 per cent of total demand, so the entire
admissible cavity range spans 34.8 to 31.0 min. Duration is a physiology and
reporting result rather than a design result.

**The reported figure depends on a declared endpoint.** With no glucose floor,
duration means time to complete exhaustion of the circulating pool, which is
not a physiological endpoint. At a 3.0 mmol/L floor the same configuration
gives 13.5 min. The floor moves the answer by a factor of three and the model
does not set it, so we report a curve rather than a value.

## 5. Fat changes the fuel budget, not the machine

A structural alternative permits fatty-acid oxidation, with visceral adipose
disaggregated out of the residual branch rather than added to the network, so
total cardiac output, the baseline fraction sum and total oxygen demand are
preserved by construction; both variants solve to the same systemic
conductance. Long-chain fatty acid does not cross the blood-brain barrier in
useful quantity, and ketogenesis requires days of fasting against a horizon of
minutes, so cerebral demand stays wholly on glucose and sets a floor on the
glucose fraction.

Fat raises supported duration from 33.7 to **47.1 min**, 13.5 to 18.9 against a
3.0 mmol/L floor, by lowering the glucose fraction from 1.000 to 0.751. It does
**not** admit a smaller pump, established two ways. The fat-capable optimum is
feasible under glucose-only metabolism as well, at 7.627 against 7.600 mL
O₂/min, so it needs no fat. And under the declared margin the two arms cost the
same to within 0.3 per cent — 8.137 against 8.164 — while the fuel rows are
unchanged. An earlier draft read the arms' unconstrained optima as fat buying a
22 per cent smaller cavity; that was two searches landing in different parts of
a flat valley, and we report the correction rather than the claim.

What is genuinely different is the limiting resource. Glucose-only, the hepatic
export rate binds. Fat-capable, **adipose perfusion binds at every flow
examined and lipolytic capacity never does.** The depot receives 52.5 mL/min,
delivers 27.0 mL O₂-equivalent/min, and covers 12.6 per cent of non-obligate
demand against a ceiling — fat covering all of it — of a 0.266 glucose fraction
and roughly 216 min. The entire gap between 12.6 per cent and that ceiling is
circulatory. **Fat availability under arrest is a perfusion question, not a
metabolic one.**

Opening the depot is a capability claim about the support circulation rather
than a design choice: the depot's vascular resistance is the subject's anatomy
and is identical under both policies, and the choice is categorical, because
any value above full conductance would mean dilating a vessel past its
physiological maximum. The policy is therefore not searched. It costs 0.91 mmHg
of mean pressure, consistent with the 2.1 per cent rise in effective systemic
conductance it produces.

## 6. Where the envelope ends

Four designs spanning the margin-holding feasible set — cavity scale 0.516 to
0.697, conductance 0.896 to 1.748 — were re-evaluated under five declared
stresses (Table 4).

**Self-supply fails at compression kappa 1 in four of four designs.** The
auxiliary myocardium cannot feed itself when its own chamber transmural
pressure is applied in full as coronary closing pressure, at any cavity size
tested. Within a single beat this mechanism drives coronary flow from 1.33 to
0.41 mL/s, with the minimum coinciding with near-peak aortic pressure. At
kappa 0 all four designs survive, at MAP 68.08 to 81.38. **The wall is
one-sided**: feasibility holds for kappa up to about 0.5 and fails above it.
Kappa rests on no measurement and decides whether the design is feasible at
all, which makes it the single most valuable missing datum in the study.

**Adverse physiology fails four of four**, raising auxiliary demand 55 to
60 per cent and failing the active organ group as well as the pump. Feasibility
holds in two of three declared scenarios.

**The ranking is stable.** Ordered by objective, the four designs come out
identically under reference, adverse and favourable physiology, so selecting on
reference is a defensible criterion among margin-holding designs.

The fat-capable arm carries one stress the glucose-only arm cannot, and passes
it: all three of its leading designs remain feasible with depot perfusion
withdrawn. Feasibility does not rest on the adipose policy; only duration does.

## 7. Four conclusions were corrected, three from one cause

We report this because it bears on how the rest should be read.

| Conclusion as first obtained | After cross-checking |
|---|---|
| Optimum at 7.600 mL O₂/min | Holds 0.005 mmHg margin; not reportable |
| No design survives kappa 0 | Four of four survive, MAP 68.08–81.38 |
| Objective ranking inverts under stress; criterion unsound | Identical ranking in all three scenarios |
| Fat admits a 22 per cent smaller cavity | Same machine to within 0.3 per cent |
| Tolerable release 1.0–1.5 s | Under 1 s once margin is required |

The first four share a mechanism. A flat objective minimised against a hard
threshold returns the most fragile design that technically qualifies; stressing
*that* design then measures the specification rather than the physiology. A
design 0.005 mmHg above a threshold moves below it under any perturbation, and
designs within 0.2 per cent of one another reorder under any change. The fifth
correction runs the other way: re-measuring latency on robust designs made the
requirement stricter.

Each correction came from an independent cross-check rather than from the run
that produced the original claim. The protocol now declares a minimum margin
and requires stressed samples to span the feasible set, and a regression test
fails if the selection reverts to the cheapest corner. **That gives the method
a guard against this failure; it is not evidence that no instance remains.**
