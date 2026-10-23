# Methods

*Draft for the arXiv preprint. Written in reader order rather than build order.
Every quantity here is traceable to a frozen configuration or a saved run; no
result appears in this section.*

---

## 1. Question and scope

We ask whether a hypothetical auxiliary biological pump, placed in parallel with
an arrested native heart, could maintain organ oxygen delivery for a useful
interval, and what such a pump would have to be. The study is entirely
computational. It establishes what a declared model permits, not what a person
would survive: no claim of survival, neurological outcome, injury avoidance or
clinical applicability is made anywhere in this work.

The reference subject is a resting adult with a physical blood inventory of
5250 mL and a whole-body resting oxygen demand of 255.848 mL O₂/min,
apportioned across organs (brain 48.11, liver 45.82, limb muscle 44.67, heart
20.79, kidneys 19.53 mL O₂/min, and the remainder). Native ventricular pumping
is lost abruptly at t = 0 and does not return.

Throughout, we distinguish four kinds of statement, following the project's
internal convention: **E**, an observation in its original published context;
**A**, that observation transferred into this model as an assumption; **D**, a
quantity produced by the simulation; and explicitly unresolved items, which are
carried as such rather than silently defaulted.

## 2. Circulation model

### 2.1 Native haemodynamics

The native circulation is the minimal closed-loop cardiovascular model of Smith
et al. (2004), which represents both ventricles with time-varying elastance,
ventricular interaction through a shared septum, a passive pericardium, and
inertial valve dynamics. We reconstructed it from the published equations and
verified the reconstruction against an independent implementation before any
extension (§6.1).

### 2.2 Regional extension

The single systemic resistance of the source model is replaced by **25 regional
branches** — brain, left and right coronary, two kidneys, hepatic arterial and
five splanchnic beds, lungs, respiratory muscle, four limb tissue compartments
and their returns, and a residual compartment — each with its own conductance
and each conserving volume through a shared incidence matrix. Branch
conductances were set so that baseline regional flow fractions match the
project's evidence-derived allocation; they are not independently fitted.

Two closures are pressure-dependent rather than purely resistive. Cerebral flow
uses a waterfall (Starling-resistor) law against an assumed intracranial
pressure, and coronary flow uses the same law against a closing pressure that
rises with the enclosing chamber's transmural pressure:

> p_close = max(p_pericardium, 0) + κ · max(p_chamber − p_pericardium, 0)

with **κ = 0.5** in the reference configuration and declared alternatives of 0
and 1. κ represents the fraction of a contracting chamber's pressure that is
transmitted to the vessels buried in its own wall. It is a model-form choice,
treated as an axis rather than a fitted parameter (§5.3).

The extended state vector carries **18 compartment volumes** (two native
ventricles, aorta, vena cava, pulmonary artery and vein, ten limb compartments,
and the two auxiliary chambers) together with eight valve flows, 26 states in
all. Total physical inventory is conserved to within 10⁻⁵ mL over every
accepted run.

### 2.3 Auxiliary pump

The auxiliary pump is a second Smith-shaped chamber pair with its own septum and
pericardium, connected in parallel through four additional inertial valves
(PU→aLV→AO on the left, VC→aRV→PA on the right) and supplied by its own coronary
branches. Its walls are obtained by similarity scaling of the native walls by a
**cavity scale** s: elastance and stiffness divide by s while dead and reference
volumes multiply by it, so that pressures are preserved under volume scaling.
Valve resistance and inertance, and the auxiliary coronary resistance, likewise
divide by s — so a larger pump is given a proportionally larger coronary bed.

The chambers are filled once, before arrest, by a declared physical transfer
from the pulmonary-vein and vena-cava pools; no blood is created. The pump is
inert until released, after which its activation follows a Gaussian driver of
period T with a one-second amplitude ramp.

**Cavity scale is a capacity proxy, not an anatomical mass.** Scaling wall
coordinates does not establish myocardial volume, implantability or
constructibility, and we make no such claim.

## 3. Oxygen and substrate transport

Metabolism is solved as a second stage over the saved circulation (§4.1). Each
organ has a prescribed oxygen demand, an extraction limit, and a small local
store sized at five seconds of its own reference demand. Oxygen enters at the
lungs through a ventilation model and is carried by blood whose oxygen content
follows the M4 binding relation.

Demand is expressed in oxygen units, but meeting it requires substrate, and the
model charges **all** oxidative demand to glucose at the stoichiometry of
complete oxidation: one mmol of glucose consumes six mmol of oxygen, or
**134.48 mL O₂ STPD per mmol**. Glucose is drawn from a circulating pool
(plasma 5.0 mmol/L, whole-blood concentration obtained by dividing by a
plasma-to-whole-blood ratio of 1.11, over 5250 mL) and replenished by hepatic
export from a glycogen store of 229 mmol/L over a 1.44 L liver.

Hepatic export is the minimum of four terms: requested replenishment; a
mobilisation capacity of 0.21 mmol/min multiplied by a liver-function fraction
of 0.5; the remaining accessible glycogen divided by a release timescale; and
forward hepatic flow multiplied by an allowed concentration increment. The
0.21 mmol/min figure derives from measured **fasting net glycogenolysis**, not
maximal stimulated output, and excludes gluconeogenesis.

Two restrictions of this substrate model are load-bearing and are reported as
axes rather than buried (§5.3): there is **no tissue substrate buffer**, so
demand must be met at every instant; and there is **no fatty-acid or ketone
pathway**, so a fuel that supplies the majority of resting oxidative ATP in a
real subject is absent. The primary analysis is run with the model as stated and
the consequences of these restrictions are reported alongside it.

## 4. Numerical methods

### 4.1 Two-stage solution

Circulation and metabolism are solved in sequence, not simultaneously. The
hydraulic state is integrated first; the saved trajectory of volumes and flows
is then replayed as a prescribed input to the metabolic model. **Metabolic state
never feeds back into the circulation.** Contraction therefore continues
unchanged under unmet demand, which is a declared limitation: the pump in this
model cannot starve, so a duration result is conditional on supply rather than
being terminated by it.

This one-way structure is also what makes the search affordable, since the
pressure screen can be decided from the hydraulic stage alone (§5.5).

### 4.2 Tolerances and paired verification

Hydraulics are integrated with an implicit Radau method with terminal valve
events, sampled at 5 ms. Metabolism is integrated with RK45 over the saved
trajectory. **Every accepted evaluation is run at two numerical settings:**

| Stage | Nominal | Refined |
|---|---|---|
| Hydraulic | rtol 10⁻⁷, max step 10 ms | rtol 10⁻⁹, max step 5 ms |
| Metabolic | rtol 10⁻⁹, max step 5 ms | rtol 10⁻¹¹, max step 2.5 ms |

A run is accepted only if the paired trajectories agree within a scaled error of
0.005, blood volume is conserved to 10⁻⁵ mL, septal root residuals stay below
10⁻⁸ mmHg, chamber work agrees with independently computed pressure–volume
quadrature within 0.05 J, and the oxygen and glucose ledgers close to 10⁻⁶ mL
and 10⁻⁸ mmol. Failures of these checks are recorded with their cause and
excluded from optimality claims; they are **never counted as passes**, and an
evaluation that fails to integrate can never be credited with satisfying a
screen.

Tolerances are frozen in the protocol and are not loosened to afford more
evaluations.

## 5. Experimental protocol

The protocol was frozen before the main runs and amended four times, each
amendment recorded with its reason, its effect size, and which analyses it
required re-running.

### 5.1 Design controls versus scenario factors

The single most important guardrail is a two-column split, enforced in code. A
search may vary only **design controls**: the four cavity scales, the activation
period, three contractilities, coronary conductance and the regional
preservation conductance. It may not vary **scenario factors** — oxygen
extraction limit, auxiliary myocardial mass coupling, work cost, demand
suppression depth, glycogen content, hepatic release capacity or liver function.
A search free to choose those would lower the requirement rather than meet it.
Scenario factors are instead evaluated as declared triples (reference, adverse,
favourable), and a candidate is judged on its **worst** declared scenario.

Size is not free. Because cavity scale sets pump capacity, it also sets
myocardial mass and hence oxygen cost: the auxiliary mass ratio is defined as
cavity scale multiplied by an uncertain specific-mass coupling, and the
activation component of myocardial demand is charged at the chosen rate 60/T
rather than at a fixed rate. The auxiliary oxygen store is sized by the same
five-second rule applied to every other organ.

### 5.2 What counts as adequate

Pressure is screened at the central aorta over complete windows. Three
thresholds are evaluated and reported for every run — mean ≥ 65, maximum ≥ 90
and minimum ≥ 50 mmHg — but **only the mean is enforced**. Adult post-arrest
haemodynamic guidance states a minimum mean arterial pressure and no systolic
target, and randomised trials comparing higher and lower mean targets found no
benefit from the higher; the systolic and diastolic thresholds derive from a
single perioperative observational cohort using radial measurement and a
five-minute exposure definition. Both criteria are reported for every design;
selecting one after seeing which passes is forbidden by the protocol.

Windows are **aligned to the activation cycle** once the pump is released, one
window per beat, with a fixed one-second window retained for the nonpulsatile
pre-release interval. A single grid would otherwise mix 1.33 cycles at 80
beats/min and 2.5 at 150, making comparisons across rates incommensurable.

A **startup allowance of two seconds measured from arrest** is declared in
advance. Windows beginning inside it are reported with their shortfall exposure
but do not decide the verdict. The allowance is measured from arrest rather than
from each design's own release, so a design that detects late pays for the
exposure rather than receiving a longer grace period.

Oxygen adequacy is screened per organ role — active, preserved, and the
auxiliary myocardium itself. The pump must feed its own muscle; this is a
constraint, not a diagnostic.

The screen horizon is **40 seconds**, with a plateau check on the final five
cycle means. Cycle-mean pressure declines toward a steady state after
activation, so the binding window is the last one evaluated and a horizon ending
mid-transient credits a value the design does not maintain. Any run whose final
cycles have not flattened to within 0.05 mmHg is flagged as horizon-limited.

A design must also hold a **minimum pressure margin of 3.0 mmHg** above the
enforced threshold to be reported as feasible. This was added to the protocol
after the first search, and the reason is a general one rather than a detail of
this model. The objective surface is flat — designs differing in all three
control variables cost within 0.2 per cent of one another — and the constraint
is a hard threshold, so minimizing cost drives the solution onto the boundary.
The design the search located satisfied every screen while holding
**0.005 mmHg** of margin, and changing its activation period by 0.0005 s, four
hundredths of one per cent, moved it to violated with two failing windows.
Nine of 152 feasible designs held under 0.5 mmHg and the cheapest were all
among them: the search concentrated on fragility by construction, correctly,
because nothing in the objective asked for margin. We therefore treat this as a
specification error rather than a search error. The margin was declared before
the affected analyses were re-run, not chosen after seeing which designs
survive it, and its cost is reported alongside it (§ Results).

The same reasoning governs how robustness is tested. Re-evaluating only the
cheapest designs under stress measures the fragility the objective selected, so
the protocol requires the stressed sample to **span the feasible set** rather
than its cheapest corner. Two conclusions drawn before this requirement were
reversed by it, and both are reported with their corrections rather than
silently replaced.

Nine unresolved requirements — regional organ pressure minima, recovery-linked
deficit thresholds, CO₂ clearance, reperfusion injury and others — are attached
to every verdict, so a numerically clean run is reported as *indeterminate*
rather than *satisfied*.

### 5.3 Declared axes

Four quantities are reported as axes rather than fixed, because each moves the
answer materially and none is settled by evidence:

| Axis | Values | What it represents |
|---|---|---|
| Compression κ | 0, 0.5, 1 | coronary self-compression by the contracting wall |
| Activation latency | 0.5–3.0 s | detection plus recruitment delay |
| Glucose floor | 0–3.5 mmol/L | plasma level below which the pool is unusable |
| Glucose fraction | 0.13–1.0 | share of oxidative demand charged to glucose |

The glucose floor is anchored at 3.0 mmol/L, the level at which trial reporting
is recommended, at which measured symptom thresholds cluster, and near which
cerebral glucose transport begins to limit utilisation; transport half-saturation
estimates span roughly 1.5–6.2 mmol/L, which is why a curve is reported rather
than a single value.

### 5.4 Endpoints

The primary endpoint is **supported duration**: the interval over which the
declared delivery screens can be met. Its substrate limit has a closed form —
the usable circulating pool divided by the net drain, where net drain is glucose
demand minus effective hepatic export — and this is verified against a coupled
replay of the drain rate rather than by simulating the full interval.

Secondary endpoints: first constraint violation bracket, startup shortfall
exposure, sustained status per constraint, chamber and positive transmural work,
auxiliary oxygen demand and unmet ledger, remaining glycogen and glucose, and
cavity scale as a capacity proxy.

A duration is always reported as *no modeled violation observed through H*,
never as indefinite support, a maximum recoverable interval, or a survival time.

### 5.5 Search

The objective is to minimise the auxiliary heart's own oxygen demand subject to
the pressure and self-supply screens, across the declared cases.

Three measured structural properties justify the method. The objective rises
monotonically with cavity scale and rate, so the optimum lies on a constraint
boundary rather than inside the feasible region; coronary conductance is nearly
free energetically, so its optimal value is the least that satisfies
self-supply; and the pressure screen depends on the hydraulic stage alone.

The search therefore **traces boundaries** rather than exploring a volume. For
each rate, the least feasible cavity scale is located by a safeguarded Brent
root-find on a continuous signed pressure margin, whose sign is tested to agree
with the boolean verdict. Bisection on a boolean would discard the magnitude
information; an unsafeguarded secant would be unsafe, because the margin is only
piecewise smooth where the binding window or metric changes. Monotonicity in
scale is a declared premise and is **checked** at every located boundary by
confirming that the point immediately below it fails.

The two constraints do not separate at the boundary: raising coronary
conductance diverts flow from the aorta and costs 1–2 mmHg, which a design
sitting exactly on the pressure boundary cannot afford. A second stage therefore
searches cavity scale and conductance jointly, stepping scale upward from each
pressure boundary and, at each step, bisecting conductance for the least value
meeting self-supply while re-checking pressure. The trade between cost and
coronary enlargement is reported as an ε-constraint frontier rather than
collapsed into a weighted sum.

Because a structural argument cannot verify itself, the located optimum is
checked against a **matched-budget scrambled Sobol sample** over a neighbourhood
of it. If boundary tracing cannot beat uniform sampling at equal cost, it has
demonstrated nothing.

All results are claimed as *best found within the tested bounds and budget*.
No global optimum is claimed.

### 5.6 Reproducibility

Every evaluation is addressed by a key carrying the design, the model-form case,
the activation latency and the screen horizon, so a cached trajectory cannot
answer a different question. The protocol file records the code commit, library
versions, and SHA-256 digests of every configuration and reference input.

## 6. Verification and limitations

### 6.1 Verification

The native reconstruction was checked against an independent implementation of
the source model before extension. Beyond the paired-tolerance checks of §4.2,
each run verifies blood-volume conservation, hepatic junction balance, oxygen
and glucose ledger closure, demand-integral identities, oxidation stoichiometry,
and that no organ store exceeds its capacity or falls below zero.

One reproduction discrepancy is carried openly rather than tuned away: the
source model's cardiac-output response to a thoracic-pressure change is 19.078%
against a 9 ± 5% target. It is stated wherever responses to thoracic pressure
are relevant.

### 6.2 Limitations carried into the results

- Metabolism does not feed back into the circulation; contraction is prescribed,
  so unmet demand never slows the pump.
- No tissue substrate buffer: demand must be met at every instant, so within-beat
  supply interruptions are charged directly as deficits.
- No fatty-acid or ketone oxidation; all oxidative demand is charged to glucose.
- Hepatic export represents only the glycogenolytic component, at fasting rather
  than stimulated rates, with no gluconeogenesis.
- The activation latencies examined include values no demonstrated detector
  achieves; latency is reported as a requirement rather than assumed.
- Cavity scale is a capacity proxy, not a determined anatomical mass.
- Regional pressure minima, CO₂ and pH, ischaemic and reperfusion injury,
  thrombosis and blood compatibility are not modelled at all.
- Coronary self-compression, demand suppression depth, extraction limits and
  liver function are declared assumptions, not measurements.

*Word count target for the submitted version: approximately 3,000 words.
Figures referenced from the Results section are not listed here.*
