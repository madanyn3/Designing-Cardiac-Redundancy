# T8.5 — Answers to the research questions

2026-10-23. Arm A, the glucose-only analysis, protocol version 5. Every answer
is traced to a saved run; the supporting detail is in
[`reports/feasibility-results.md`](../reports/feasibility-results.md).

Standing qualification on all seven answers: nothing here is a statement about
survival, neurological outcome, injury avoidance or clinical applicability. The
nine unresolved requirements of `configs/m4-constraint-policy.json` remain
attached, so a numerically clean run is reported as *indeterminate*, not
*satisfied*.

---

## 1. Did we find any configuration satisfying all declared modeled constraints?

**Yes, within a bounded envelope, and only with an activation latency no
demonstrated detector achieves.**

A design satisfying the pressure screen and auxiliary self-supply exists at
cavity scale 0.536, 57.5 beats/min and coronary conductance 1.748, costing
8.137 mL O₂/min and holding 3.10 mmHg of pressure margin. Thirty-six designs
hold the declared margin and feed their own myocardium.

The cheapest design the searches found costs 7.600 mL O₂/min but holds
**+0.005 mmHg** of margin, and a four-hundredths-of-a-per-cent change in its
activation period flips it to violated. It is not a reportable design, and
protocol v5 declares a 3.0 mmHg margin for that reason. Margin costs 7.1 per
cent more oxygen. Three independent searches agree on the unconstrained valley,
two of them landing within 0.001 mL O₂/min of each other from different
directions, which is why the fragility is a property of the specification
rather than of any one search.

The envelope has three walls, and the design is inside all three only at
reference physiology, compression kappa 0.5, and a prescribed 0.5 s activation.
At the M6 controller's observed 2.55 s, **nothing passes at any cavity scale.**

T8.4 tested the four leading designs against five stresses. **None survived any
stress other than favourable physiology**: 0 of 4 at adverse, 0 of 4 at kappa 0,
0 of 4 at kappa 1. The feasible region is narrower than the single reference
result suggests, and section 6 says in which direction each wall lies.

## 2. What pump capacity, operating settings and circulation allocation did it require?

Cavity scale 0.536 of the native reference, which is a **capacity proxy and not
a determined anatomical mass**; 57.5 beats/min; coronary conductance 1.748; and
the frozen preservation allocation, with the mandatory set unchanged.

There is a second reportable design, and the choice between them is the clearest
engineering result in the milestone: **scale 0.749 at 46.0 beats/min with
coronary conductance 0.953**, costing 9.596 mL O₂/min and holding +8.97 mmHg.
It is 17.9 per cent more expensive and needs a 40 per cent larger cavity, but
it sits at or below nominal conductance, so it depends on the unjustified
0.5–2.0 bracket not at all, and it holds nearly three times the required
margin. Cheapness and independence from assumption point to opposite corners of
the feasible set.

The optimum is a **shallow valley rather than a point**. Among the leading
designs, conductance spans 1.388 to 1.808 while the objective varies 0.17 per
cent, and designs differing in all three variables cost within 0.2 per cent of
one another. Reporting a single parameter triple would overstate what the
search resolved.

## 3. How does it compare with the normal-heart reference?

The auxiliary heart consumes **7.600 mL O₂/min, about 7 per cent of the 108.6
mL O₂/min total demand** in the supported state. Across the explored range the
pump is 7 to 16 per cent of total demand.

The comparison that matters is not oxygen cost but what the circulation
delivers. The supported state meets a suppressed demand of 108.6 mL O₂/min
against a resting whole-body 255.8, so it is not a normal circulation operating
at lower output; it is a different circulation serving a reduced and selectively
allocated demand. Framing it as a fraction of normal cardiac output would
mislead.

## 4. Which resource or constraint limits it first?

**Detection latency.** Eight of nine declared operating cases fail, and the
most revealing is `support_with_suppression`: it reaches MAP 67.63, **above the
enforced threshold**, and fails only because the controller releases at 2.55 s
while the declared startup allowance ends at 2.0 s. The design is not what
fails these cases.

The ordering of binding constraints:

1. **Detection latency** — tolerable release is 1.0 to 1.5 s depending on
   cavity size; the controller delivers 2.55 s.
2. **Auxiliary self-supply** — which sets the minimum cavity, because the
   coronary bed scales with cavity size, so a smaller pump is worse at feeding
   itself. All four pressure-boundary designs are infeasible on self-supply
   even at the upper conductance bound.
3. **Substrate** — the softest of the three. It sets duration, not feasibility.

## 5. How do preservation obligations change the design?

Preservation is what makes the problem tractable: the supported state requests
108.6 mL O₂/min against a resting 255.8, and the mandatory set holds brain,
both coronaries, protected kidney, hepatic arterial and respiratory muscle at
full conductance while the rest is suppressed.

Two consequences are measurable. Removing suppression
(`support_no_suppression`) drops MAP to 46.06 and fails both the active and
pump groups, so the demand reduction is load-bearing rather than incidental.
And under adverse physiology the designs **cheapest at reference are the most
fragile**: at 60 beats/min the *active* group fails, not merely the pump, while
costlier designs fail on self-supply alone.

## 6. Which conclusions survive uncertainty and structural alternatives?

**Survives:**

- That cavity scale was the unexamined variable. M5 and M7 reported 0 of 84
  cases reaching the pressure screen while holding scale at the 0.25 fixture;
  inventory never constrained it, and every scale in 0.1–1.0 is admissible.
- That detection latency binds before the pump does.
- That the coronary bed scaling with cavity size makes self-supply easier for
  larger pumps, not harder.
- That feasibility does not depend on the coronary conductance assumption: a
  feasible design exists at or below nominal, at 9.596 mL O₂/min and scale
  0.749. The unjustified 0.5–2.0 bracket buys optimality, not feasibility.
- That the hepatic export **rate** limits duration, not the glycogen store,
  which would last 3,141 min at the effective cap against 33.8 for the pool.

**Does not survive:**

- Feasibility at **either** compression alternative, and the two fail for
  opposite reasons. At kappa 1 self-supply fails, as expected. At kappa 0 the
  **pressure** screen fails: less coronary compression means more coronary flow,
  which is more runoff from the aorta, and MAP falls to 63.44 against the
  65 mmHg threshold. An earlier statement that kappa 0 was benign came from a
  larger-cavity design with margin to spare; the optimized designs sit near MAP
  66 to 67 and a fraction of a millimetre of mercury flips them. Feasibility
  holds only in a **narrow band around kappa 0.5**, bounded by pressure below
  and self-supply above.
- Feasibility under adverse physiology, which raises auxiliary demand 54 to 62
  per cent and fails both the active and pump groups in all four designs.
- **The selection criterion itself.** The order by objective at reference
  inverts under *both* adverse and favourable physiology: the design chosen as
  cheapest drops to third in each. Selecting on reference alone was therefore
  the wrong criterion, which is precisely what T8.4 exists to detect. The
  reversal is reported rather than resolved, because resolving it would mean
  choosing a different scenario to optimize against, and no declared scenario
  has a stronger claim than the others.
- The duration figures under any fuel assumption other than glucose-only. This
  is the largest single bias in the analysis and is addressed in Arm B.

**Changed during the work, on evidence rather than convenience:** the enforced
pressure criterion became mean pressure alone, because adult post-arrest
guidance states a MAP target and no systolic one; and pressure windows became
cycle-aligned, because a fixed grid mixes different numbers of cycles across
rates. Both are recorded as protocol amendments with their reasons and the
analyses they required re-running.

## 7. What remains impossible to determine with this model?

- **Whether any of this is survivable.** Deficits are exposure proxies, not
  injury doses; no recovery-linked threshold exists for any tissue here.
- **Whether the organ is constructible.** Cavity scale is a hydraulic
  similarity parameter. Nothing in the model speaks to myocardial mass,
  implantation, immune compatibility or how such a pump would be controlled
  biologically.
- **Whether detection within 1.5 s is achievable.** The study turns a modelling
  assumption into a design requirement and hands it back; it cannot say whether
  the requirement can be met.
- **Regional adequacy.** The pressure screen is central-aortic. Regional organ
  pressure minima are unresolved, so no claim of brain protection follows from
  passing it.
- **The value of compression kappa.** It bounds the feasible region from both
  sides, pressure below and self-supply above, and rests on no measurement.
  This is the single most valuable missing datum in the study, and T8.4 raises
  rather than lowers its importance.
- **CO₂ and pH, ischaemic and reperfusion injury, thrombosis and blood
  compatibility**, none of which are modelled at all.
- **Anything outside the sampled domain.** Boundaries reported here are
  boundaries of what was tested, not of what is possible.

---

## 8. Arm B — where the answers change, and where they do not

2026-10-23. The fat-capable structural alternative, carried through the same
seven questions. Detail in
[`reports/feasibility-results-armb.md`](../reports/feasibility-results-armb.md).
Only the differences are stated; where Arm B is silent, Arm A's answer stands
unchanged.

**1. Is there a satisfying configuration?** Yes, and the margin requirement
applies equally: Arm B's search optimum holds +1.011 mmHg and is not
reportable. Sixteen of its 23 feasible designs clear 3.0 mmHg, the cheapest at
scale 0.443, 74.5 beats/min and conductance 1.471, costing 8.164 mL O₂/min.

**2. What capacity did it require?** **The same machine as Arm A, to within
0.3 per cent.** Under the declared margin, Arm A costs 8.137 mL O₂/min and
Arm B 8.164. An earlier draft read the arms' unconstrained optima as fat
buying a 22 per cent smaller cavity; a direct cross-check refuted that — Arm
B's optimum is feasible under Arm A's glucose-only metabolism at 7.627 against
7.600 — and the margin-constrained comparison now confirms it independently.
**Fat does not change the pump.**

**3. Against the normal-heart reference?** Unchanged. The fuel model does not
enter the hydraulics, and the circulation is the same circulation with one
territory named: total output, the baseline fraction sum and total demand are
preserved by construction, and both arms solve to 0.98644 mL/s/mmHg.

**4. What limits it first?** **The limiting resource changes, and this is Arm
B's substantive result.** In Arm A the hepatic glucose export *rate* binds. In
Arm B fat relieves part of that demand — glucose fraction 1.000 to 0.751,
supported duration 33.7 to 47.1 minutes, 13.5 to 18.9 against a 3.0 mmol/L
floor — but a new limit replaces it: **adipose perfusion, not lipolytic
capacity, binds at every flow examined.** The depot receives 52.5 mL/min and
delivers 27.0 mL O₂-eq/min, covering 12.6 per cent of non-obligate demand. Had
fat covered all of it the glucose fraction would fall to the brain's 0.266 and
duration reach roughly 216 min, so **the entire gap between 12.6 per cent and
that ceiling is circulatory.** Fat availability in arrest is a perfusion
question, not a metabolic one.

**5. How do preservation obligations change the design?** Arm B adds one
obligation Arm A does not have, and it is a **capability claim rather than a
design choice**: whether the support circulation opens the depot's existing
conductance. The depot's vascular resistance is the subject's anatomy and is
identical under both policies. The choice is categorical, because any value
above full conductance would mean dilating a vessel past its physiological
maximum, so the policy is **not searched** — the optimizer never selects the
assumption that flatters the result. Perfusing the depot costs 0.91 mmHg of
mean pressure at a fixed design, consistent with the 2.1 per cent rise in
effective systemic conductance it produces. That is the price of the fuel, and
it is modest.

**6. Which conclusions survive?** Arm B carries one stress Arm A cannot, and
passes it: **all three leading designs survive with the depot suppressed**, so
feasibility does not rest on the adipose policy — only the duration does. Of
the compression cases, 3 of 3 survive kappa 0 and **0 of 3 survive kappa 1**,
self-supply failing in every instance. The kappa 1 wall is the one boundary
both arms report identically, and it is the firmest negative result in the
milestone.

**7. What remains undetermined?** Arm A's list in full, plus three specific to
Arm B. The **arterio-venous fatty-acid difference is assumed, not sourced**, at
1.0 mmol/L against reported plasma concentrations of 0.5–1.0 fasting and 1.5–2
prolonged; fat supply scales linearly with it, so it deserves the evidence
treatment the pressure criterion and the glucose floor received. The **visceral
depot mass is assumed** at 2.0 kg. And **ketone contribution is zero**, correct
for a horizon of minutes and wrong for a multi-day one, since ketogenesis needs
days of fasting — which is why the brain's demand stays wholly on glucose here
and sets the floor on the glucose fraction.
