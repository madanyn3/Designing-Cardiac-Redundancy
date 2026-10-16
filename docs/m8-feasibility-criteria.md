# What counts as adequate, and why — review before T8.0

2026-10-16. Working note for the four questions raised at the start of T8.0: how hard the
pressure thresholds really are, why the circulation takes seconds to reach them, whether a
prioritized territory set changes the answer, and where the glucose budget enters.

Nothing here is a protocol change yet. Changes to the frozen T8.1 protocol need a version
bump and a recorded amendment.

---

## 1. The three thresholds do not have equal standing

The screen is **MAP ≥ 65, SBP ≥ 90, DBP ≥ 50 mmHg** jointly, over complete windows. Its
provenance is [`docs/m2-pressure-requirements.md`](m2-pressure-requirements.md), and it is
worth being precise about what each number is, because **the component that currently
binds our minimum cavity size is the one with the weakest support.**

| Threshold | Source | Standing |
|---|---|---|
| MAP ≥ 65 | AHA 2025 post-ROSC recommendation (BP01) | Strongest. A clinical recommendation in a post-ROSC population. |
| DBP ≥ 50 | Ahuja 2020, 23,140 noncardiac surgical patients (BP05) | Observational association; radial measurement; defined on pressures sustained ≥ 5 min. |
| SBP ≥ 90 | Ahuja 2020 (BP05), same cohort | Same limitations **and** AHA 2025 *removed* its previous systolic target because the trials evaluated MAP. |

Two further points from the contract itself:

- **65 is not a cliff.** The BOX trial (BP02) randomised post-arrest MAP targets of 63 versus
  77 and found no significant difference in death, severe disability or coma. The contract
  draws the conclusion explicitly: this "does not establish equivalence or a safe
  physiological floor" — it cuts both ways, and the honest reading is uncertainty, not a
  sharp boundary at 65.
- **Pulse pressure was deliberately not adopted.** The contract records the 35 mmHg PP
  association as "an observational comparator, not a mandatory pump constraint," and states
  that the PP association "is not a demonstrated need for that pulsatility in an auxiliary
  pump."

So an auxiliary design rejected *only* because SBP sits at 85 instead of 90 is being rejected
by a transferred systolic threshold that the relevant clinical body has itself retired. That
needs saying wherever the minimum feasible capacity is reported.

### What the model can and cannot say about the brain specifically

The contract is explicit: **"Do not declare an independent cerebral-pressure minimum."**
Hori 2017 (BP03) gives a population autoregulation lower limit of 65 ± 12 mmHg with
individual limits crossing the population optimum in substantial subsets, and Sundgreen 2001
(BP04) found autoregulation impaired in 8 of 18 post-arrest patients with markedly varying
limits. The contract's conclusion is that "a pressure that passes an aggregate screen may
still fail an individual's cerebral supply" — i.e. the evidence argues for *heterogeneity*,
not for a number.

Mechanically, the model has a **diagnostic** rather than a constraint:
`P_AO − max(P_VC, ICP)` with ICP 10 mmHg (5/10/15 sensitivity), and
`docs/m2-pressure-requirements.md` §3 states that "the current brain resistor does not
incorporate ICP into its flow law." So cerebral perfusion pressure is reported, not enforced,
and brain flow in the model follows from conductance and pressure like any other bed.

**Consequence for the claim we can make:** the pressure screen is a central-aorta research
screen. It is not an organ-wise safety composite, and no result from it can be stated as
brain protection.

## 2. The premise about how the target was built

The screen was **not** constructed as the maximum over individual organs' safe minima. It is
a transfer of a central-aorta clinical reference, and §3 of the contract lists every organ
row as *unresolved*: no established minimum for native coronary maintenance, none for the
contracting auxiliary heart, renal MAP 75 is a sensitivity setting rather than a requirement,
liver anchors are labelled assumptions, and the pulmonary row is a reproduction benchmark.

That matters for your prioritization question: there is no per-organ stack to re-prioritize,
because the per-organ minima were never established. What the model *does* have is per-organ
**oxygen delivery and unmet-demand ledgers**, which are computed and conserved. Prioritization
therefore has to be expressed in delivery terms, not pressure terms.

## 3. Prioritization is already partly implemented — and it has the right sign

`BilateralSupport.apply_support_controls`
([bilateral_support.py#L51](../src/cardiac_redundancy/bilateral_support.py)) already keeps a
mandatory set at full conductance and suppresses the rest:

```
mandatory = {brain, coronary_left, coronary_right, kidney_protected,
             liver_arterial, respiratory_muscle}      → factor 1.0
everything else supplied from AO                      → factor = preservation_conductance
lungs / pulmonary_exchange                            → from the lung allocation
```

So brain and respiratory muscle are already protected, and `preservation_conductance` (a
design control, bounds 0.1–0.5, currently fixed at 0.25 in every run to date) is the
prioritization dial. **It has not been swept once.**

Two reasons it is the most promising untested lever for lowering cavity size:

1. **It raises pressure.** Suppressing non-priority beds lowers total systemic conductance,
   so the same auxiliary output produces a higher mean pressure.
2. **It shortens the transient** (see §4): less runoff means the arterial compartment
   refills faster.

Narrowing the mandatory set further — your suggestion, e.g. dropping `kidney_protected` and
`liver_arterial` — pushes in the same direction and is a declared scenario choice rather than
a model change. The cost is explicit and must be reported, not hidden: those territories then
accumulate oxygen deficits, and the model cannot say what that deficit means for organ
function, because recovery-linked deficit thresholds are in the unresolved list.

**Proposed amendment (not yet applied):** report a second, clearly labelled criterion
alongside the central pressure screen —

> **Prioritized-territory adequacy:** no sustained unmet oxygen demand in
> {brain, respiratory muscle, lungs} ∪ {auxiliary myocardium}, with the deficits accumulated
> by de-prioritized territories reported in full.

This is *better supported by this model* than the pressure triple, because per-organ delivery
is computed and conserved whereas regional pressure minima are unresolved. It must be
reported **in addition to**, never instead of, the pressure screen — otherwise it is
relabelling a failure as a pass, which the protocol forbids.

## 4. Why it takes seconds — the mechanism, measured

The settling time is **arterial volume re-accumulation**, and it is quantifiable. Aortic
compartment volume and pressure, scale 0.6, from the saved trajectories:

| t (s) | Oracle release at 0.5 s | Observed release at 2.55 s |
|---|---|---|
| 0.0 | 308.6 mL / 109.7 mmHg | 308.6 mL / 109.7 mmHg |
| 0.5 | 280.9 mL / 90.5 mmHg | 280.9 mL / 90.5 mmHg |
| 1.0 | 271.1 mL / 83.7 mmHg | 258.8 mL / 75.2 mmHg |
| 2.0 | 268.0 mL / 81.6 mmHg | 226.4 mL / 52.8 mmHg |
| 2.6 | 276.8 mL / 87.6 mmHg | **212.3 mL / 43.1 mmHg** |
| 4.0 | 283.7 mL / 92.4 mmHg | 243.4 mL / 64.6 mmHg |
| 6.0 | 259.9 mL / 76.0 mmHg | 275.9 mL / 87.0 mmHg |

Three things follow.

**The arterial compartment is stiff.** 212 → 284 mL spans 43 → 93 mmHg, so about
0.7 mmHg/mL, i.e. a lumped arterial compliance near 1.4 mL/mmHg. Small volume changes move
pressure a lot — which is why the recovery is steep once net inflow is positive.

**Delay costs twice.** Drainage runs at roughly 55 mL/s in the first half second. By 0.5 s the
deficit is ~28 mL; by 2.55 s it is ~96 mL. So a late release loses the early windows *and*
hands the pump three times the volume debt to repay. Observed release refills 212 → 276 mL in
3.4 s, a net ~19 mL/s — which puts the settling time at roughly deficit ÷ net inflow, and
explains the ~6 s directly.

**Factors beyond activation delay and cavity size.** From the mechanism, the settling time
depends on:

| Factor | In the model | Status |
|---|---|---|
| Arterial volume deficit at release | drainage during pre-activation | set by delay; the dominant term |
| Net inflow = output − runoff | `preservation_conductance` sets runoff | **design control, never swept** |
| Activation ramp | `auxiliary_ramp_s`, fixed at 1 s | declared capability, could be a scenario |
| Inlet resistance / inertance | scale as 1/cavity scale in `BiologicalPumpModel` | moves with scale already |
| Contractility | `left/right/septal_contractility`, 0.5–2.0 | **design control, never swept** |
| Stroke volume per beat | `period_s`, 0.4–1.5 s | **design control, never swept** |
| Arterial compliance | from `baseline.yaml` | not a design control |
| Donor pool filling (PU, VC) | preload available to the auxiliary inlets | scenario-dependent |

So: at least three unswept design controls, one of which (`preservation_conductance`) acts on
pressure *and* settling time *and* total oxygen demand in the same favourable direction.

## 5. Where the glucose budget binds — measured

Stage 2 of the evaluator ([`scripts/run_m8_metabolism.py`](../scripts/run_m8_metabolism.py))
charges the auxiliary heart for its own size and rate, as T8.1 requires. On the scale 0.6
design that passes the pressure screen under oracle activation, reference physiology:

| Organ role | Requested O₂ | Unmet | Sustained status |
|---|---|---|---|
| Active (brain, respiratory muscle, lungs, protected kidney, liver, coronary) | 11.159 mL | 0.000000 | satisfied |
| Preserved (suppressed territories) | 5.203 mL | 0.000000 | satisfied |
| **Auxiliary myocardium (the pump itself)** | 1.885 mL | **0.003042 mL** | **violated** |

Auxiliary demand is 11.624 mL O₂/min at a coupled mass ratio of 0.600 — 2.4× the reference
fixture's charge, exactly as the size coupling intends. Liver glycogen is essentially
untouched (329.74 of 329.76 mmol remaining).

**The binding constraint is the pump's own instantaneous self-supply, and nothing else.**
Every other territory is satisfied. The shortfall is 0.16% of the pump's own request, it is
glucose-limited rather than oxygen-limited (M7: the oxygen-limited term is zero at every
saved sample), and its two known causes are both declared model limitations: within-beat
coronary compression, and **no tissue glucose buffer at all** (L-M4-F01) so the myocardium
must be supplied instantaneously at every sample. T7.4's finite-reserve diagnostic already
showed that a modest reserve covers much of a short deficit.

This is the tradeoff you asked for, and it is sharper than expected: **capacity buys pressure
and costs self-supply**, and the two requirements are currently satisfied in disjoint regions
of cavity scale — but only because the self-supply side is evaluated with no buffer. Whether
the gap is real or an artefact of L-M4-F01 is now the single highest-value question in M8, and
it is a structural question, not a search question.

## 6. What I propose doing next, in order

1. **Sweep the three unswept design controls** at low cavity scale (`preservation_conductance`,
   `period_s`, contractility) under oracle activation. If a smaller, slower, more strongly
   suppressed design clears the screen, the minimum capacity — and therefore the self-supply
   burden — drops.
2. **Map latency as its own axis.** 0.5 s is an oracle; 2.55 s is the current observed policy.
   The feasible region in (capacity, latency) is a more honest result than either endpoint and
   is a direct M6 design requirement.
3. **Decide the buffer question before the main search**, because it determines whether the
   objective is "satisfy self-supply" or "report a bounded artefact". This is a protocol
   decision, not a result.
4. **Then** run the search, with glucose cost in the objective rather than as a post-hoc check.
5. Report exposure — time to first satisfied window, integral shortfall — as primary outputs
   alongside the strict pass/fail, which is what the M2 contract asks for and what makes a
   transient legible.

## 7. One protocol detail to fix

The M2 contract says to evaluate "integral mean, maximum and minimum over each auxiliary
activation cycle" and to use the 1 s window only "for irregular/nonpulsatile runs". The frozen
T8.1 policy uses 1 s windows throughout, which at a 0.75 s period mixes 1.33 cycles and will
mix differently once `period_s` varies. Cycle-aligned evaluation should be added before the
main runs, as a version 2 amendment.
