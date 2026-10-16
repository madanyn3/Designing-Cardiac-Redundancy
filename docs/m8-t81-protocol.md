# T8.1 — Frozen M8 experiment protocol

2026-10-16. **Status: frozen at version 1.** Machine-readable form:
[`configs/m8-protocol.json`](../configs/m8-protocol.json), generated from the code by
[`scripts/freeze_m8_protocol.py`](../scripts/freeze_m8_protocol.py) at commit `f198c8d`.
Amendments must bump `version`, say what changed, and rerun the affected analyses.

This note explains what was frozen and why, with pointers into the source. It adds
no results. Numbers here are either measured in this task or copied from an M5/M7
report; measurements made here are labelled **measured**.

---

## 1. What M8 is searching over, and what it is not

The single most important guardrail is a two-column split, enforced in code by
[`src/cardiac_redundancy/design.py`](../src/cardiac_redundancy/design.py).

**Column A — design controls (`DesignPoint`).** Things an engineer of the auxiliary
heart could choose. A search may move these freely inside their declared bounds.

| Control | Bounds | Where it acts |
|---|---|---|
| `left_scale`, `right_scale`, `septal_scale`, `pericardial_scale` | 0.1 – 1.0 | cavity allocation, valve R/L and wall laws in `BiologicalPumpModel.__init__` ([biological_pump.py#L35-L50](../src/cardiac_redundancy/biological_pump.py)) |
| `period_s` | 0.4 – 1.5 s (40–150 /min) | auxiliary driver period in `ControlledCirculation.pressures` |
| `left_contractility`, `right_contractility`, `septal_contractility` | 0.5 – 2.0 | wall elastance via `scaled_wall` |
| `coronary_conductance` | 0.5 – 2.0 | auxiliary coronary resistance in `SensitivityCirculation.flows` |
| `preservation_conductance` | 0.1 – 0.5 | regional conductance of the non-mandatory beds, `apply_support_controls` ([bilateral_support.py#L51](../src/cardiac_redundancy/bilateral_support.py)) |

**Column B — scenario factors (`ScenarioPoint`).** Uncertain or speculative
physiology: extraction limit, specific myocardial mass, work cost, suppressed-demand
depth, liver glycogen, release cap and hepatic export fraction. These are **evaluated
as declared scenarios, never chosen by a search.** A search allowed to pick them would
lower the requirement rather than meet it — it would discover that a pump looks
excellent once you assume tissues need less oxygen.

Three scenario triples are declared: `reference` (the M5/M7 values), `adverse` and
`favourable`, each taking every factor to the same end of its interval. These are
deliberate corners, not probabilities.

`registry()` asserts that the two columns partition the declared bounds exactly, and
`tests/test_design.py` fails if a factor ever moves columns.

## 2. Two couplings that close free lunches

A search finds whichever bookkeeping gap is cheapest. Two were open.

**Size was free.** M7 varied `mass_ratio` — the auxiliary myocardial mass that sets its
oxygen cost — as an independent demand parameter, which was correct there, because
cavity scale was never varied (`run_m7_global_hydraulics.py` leaves the scales at the
0.25 reference and changes only contractility). In M8 cavity scale *is* a control, so an
uncoupled `mass_ratio` would let the optimizer buy a four-times-larger pump without
paying for its muscle. The protocol therefore declares

> `mass_ratio = cavity_scale × mass_coupling`, with `mass_coupling ∈ [0.8, 1.6]` carrying
> the uncertainty in specific mass per unit cavity volume.

At the reference point this returns exactly 0.25, reproducing the M5/M7 demand parameter —
a continuity check in `test_reference_point_reproduces_the_verified_fixture`.

**Rate was free.** `side_demand` charges an activation cost proportional to
`rate_bpm/80`, but [`scripts/run_m7_global_metabolism.py:202`](../scripts/run_m7_global_metabolism.py)
passes the literal `80`. That was right for M7, where the period was always 0.75 s =
80 /min. With `period_s` as a control it becomes a free lunch: faster beating, more output,
no extra activation cost. The M8 evaluator must pass `design.rate_bpm = 60/period_s`.

Coronary capacity is *already* coupled: the auxiliary coronary resistance divides by the
cavity scale in `BiologicalPumpModel.flows`, so a larger pump gets a proportionally larger
coronary bed, and `coronary_conductance` multiplies that as a separate design hypothesis.

## 3. What counts as success

Implemented in [`src/cardiac_redundancy/feasibility.py`](../src/cardiac_redundancy/feasibility.py),
which reuses the screens from `constraints.py` so an M8 verdict is the same arithmetic the
earlier milestones reported.

**The pressure screen is the declared triple, not the mean alone.** From
`configs/auxiliary-heart-design.yaml`: MAP ≥ 65, SBP ≥ 90, DBP ≥ 50 mmHg, jointly, over
complete 1 s windows, with zero violating windows. M5 and M7 reported the mean-65 screen;
the frozen M8 target is all three. A design can satisfy the mean and still fail.

**The startup allowance is fixed at 2 s from arrest**, and the same for every design —
1.5 s of activation latency plus ramp (`ActivationTiming` defaults: 0.25 + 0.05 + 0.2 + 1.0),
rounded up to the next whole window. Windows beginning before it are reported with their
shortfall exposure but do not decide the verdict. Measuring the allowance from arrest rather
than from each design's own release is deliberate: a design that detects late pays for the
exposure instead of receiving a longer grace period
(`test_late_release_does_not_earn_a_longer_allowance`).

**Oxygen screens run per group** — active, preserved, and the pump itself. The auxiliary
heart must supply its own myocardium; that is a constraint, not a diagnostic.

**Numerical failure and model-domain exit stay separate from physiological violation**, so a
design can never be credited with a pass it obtained by failing to integrate
(`test_a_domain_exit_cannot_be_credited_as_a_pass`).

**A pass means "no violation of the declared screens through the tested horizon."** Not
support duration, not survival, not recovery. The unresolved requirements from
`configs/m4-constraint-policy.json` — regional pressure minima, CO₂/pH, injury, reperfusion,
thrombosis and the rest — remain attached to every verdict, which is why a clean run can
report `indeterminate` rather than a pass when they are declared.

## 4. What one evaluation costs, and the budget that follows

**Measured** by [`scripts/profile_m8_evaluation.py`](../scripts/profile_m8_evaluation.py)
on the unchanged M7 pipeline, both paired numerical settings included:

| Stage | Wall time |
|---|---|
| Hydraulic (Radau, nominal + refined) | **110.7 s** |
| Metabolic replay (RK45, nominal + refined) | **166.3 s** |
| One complete evaluation | **277.1 s** |

Recorded in `reports/m8-protocol/timing/single-evaluation.json`. Machine-specific, not a
portable claim.

This makes the search affordable only because of a structural fact: **the pipeline is
one-way.** Hydraulics are solved first and metabolism is replayed over the saved
trajectory, so metabolic state never feeds back into pressures or flows. The pressure
screen therefore needs only the hydraulic stage. The protocol stages the search:

- **Stage 1** — pressure feasibility on hydraulics alone: 3000 evaluations ≈ **7.7 h** on 12 workers.
- **Stage 2** — coupled evaluation of candidates that already pass pressure, across scenarios: 600 evaluations ≈ **3.9 h**.

If supply-limited contraction is ever implemented, the coupling becomes two-way and the
stages must be merged again. That dependency is recorded in the protocol.

Tolerances are frozen and are **not** to be loosened to afford more evaluations; that was
tried and rejected in the M7 pilot.

## 5. One bound established before searching

A larger auxiliary cavity is filled by a one-time transfer from the pulmonary-vein and
vena-cava pools, and total blood volume is fixed at 5250 mL. So the first question is
whether inventory limits the pump's size.

**It does not.** [`scripts/screen_m8_domain_bounds.py`](../scripts/screen_m8_domain_bounds.py)
builds the model at every scale in the declared range and checks the initial state
(`reports/m8-protocol/domain-bounds.json`):

| Cavity scale | aLV + aRV | PU margin | VC margin |
|---|---|---|---|
| 0.10 | 9.5 + 9.1 mL | 771.3 mL | 259.5 mL |
| **0.25 (reference)** | **23.7 + 22.7 mL** | 757.1 mL | 245.9 mL |
| 0.50 | 47.3 + 45.4 mL | 733.4 mL | 223.2 mL |
| 1.00 | 94.7 + 90.7 mL | 686.0 mL | 177.8 mL |

Every scale in 0.1–1.0 is admissible; the binding margin is an unrelated limb compartment
at 1.5 mL. The reference auxiliary heart holds about 24 mL per cavity, and nothing in the
inventory stopped it being four times larger.

**This matters for how M8 is framed.** The quarter scale was a declared test fixture —
`configs/m2-frozen-overlay.yaml` calls it a "non-optimized hydraulic test fixture, not
proposed anatomical size" — and no coupled run has ever been done above it. The M5/M7
pressure failure (reference MAP **measured** at 49.01 mmHg, auxiliary left output
1.8051 L/min, release at 2.6 s) is therefore not yet evidence that no admissible design
reaches the screen. It is evidence that *this fixture* does not. Cavity scale is the
largest unexplored axis in the study, and T8.0 explores it.

The caveat stays attached: scale is a capacity proxy. An admissible scale is not a
demonstrated anatomical mass, and admissibility is not adequacy.

## 6. Deliverables and gate

| T8.1 requirement | Where |
|---|---|
| Frozen protocol | `configs/m8-protocol.json` (version 1, commit recorded, 20 input digests) |
| Parameter / constraint registry | `src/cardiac_redundancy/design.py`, `registry()` |
| Tested endpoint evaluator | `src/cardiac_redundancy/feasibility.py` |
| Code version, provenance, numerical settings, analysis choices | recorded in the protocol file |
| Tests | `tests/test_design.py` (23), `tests/test_feasibility.py` (21), `tests/test_m8_protocol.py` (12) |

`tests/test_m8_protocol.py` fails if the frozen file and the validators ever disagree, so
protocol drift is a test failure rather than a silent inconsistency.

## 7. Open decisions carried into T8.0

1. **Post-shortfall behaviour.** Contraction is prescribed: unmet demand never reduces
   pumping (M7 condition 2, L-M5-E03). A conservative reading ends the supported-performance
   claim at auxiliary metabolic inadequacy. The alternative — supply-limited contraction —
   makes the coupling two-way, invalidates the staged budget above, and is a larger piece of
   work than M8.
2. **Horizon.** The protocol declares 10 s, inherited from M6/M7. A longer horizon changes
   what the finite stores mean and multiplies every cost in §4.
3. **Objective shape.** Currently: pressure margin first, then capacity and energy among
   passing candidates. If Stage 1 finds nothing, the result is a bounded negative one and the
   objective becomes "how far from the screen, and why" — worth planning for before the runs,
   not after.
