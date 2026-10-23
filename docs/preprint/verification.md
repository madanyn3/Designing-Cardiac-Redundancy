# Verification

Draft for the preprint. This section states what was checked, what failed, and
what is carried unresolved. Two named limitations — **R-M3-01** and
**L-M4-F01** — are stated here rather than in a footnote, because both bound
how the results may be read.

## 1. Reproduction of the source model, and the discrepancy that failed

The native circulation is a reconstruction of a published minimal
cardiovascular model. Before any extension, it was reproduced twice: once
directly, and once as an independent implementation written from the paper's
equations. The two agree to within **0.0032 per cent** on every compared
metric, which establishes that the reconstruction is faithful to its own
specification.

It is not faithful to one published response. The source model's reported
cardiac-output decrease under a thoracic-pressure change is 9 ± 5 per cent; the
reconstruction gives **19.078 per cent**, outside that tolerance.

This is recorded as **R-M3-01** and carried as an accepted, explicit
limitation. No source parameter and no tolerance was adjusted to close it. The
alternative — tuning a parameter until the screen passes — would have produced
a model that matched one published number while losing the property that makes
reproduction meaningful. The failure is stated wherever responses to thoracic
pressure are relevant, and no result in this paper depends on a
thoracic-pressure manoeuvre.

Six coarse baseline screens pass. The original milestone gate, which required
the thoracic screen among them, **remains failed**, and the reconstruction was
closed under an explicitly recorded exception rather than by redefining the
gate.

## 2. The substrate model is verified for oxygen, not for whole metabolism

The oxygen and glucose subsystem closes its ledgers exactly (§3). It does not
follow that the substrate model reproduces normal physiology, and the
no-suppression baseline gate is satisfied **conditionally for the verified
oxygen subsystem only**. This is **L-M4-F01**.

Specifically: there is no intracellular or interstitial glucose buffer, and in
the primary analysis no alternative substrate. Demand must be met at every
instant from concurrent delivery, so a within-beat shortfall is recorded as a
deficit where a real tissue would draw on a local pool. This makes the
glucose-only arm **conservative on short timescales and biased on long ones** —
the bias the fat-capable arm exists to bound. Neither buffering nor alternative
substrates is a reason to erase a recorded deficit or to grant free energy, and
neither is treated as one.

## 3. What every run checks

Each evaluation verifies, and aborts on violation rather than reporting a
number:

- **Blood-volume conservation** across all compartments.
- **Hepatic junction balance**, since the splanchnic split introduces a
  confluence that an incorrect conductance assembly would silently unbalance.
- **Oxygen and glucose ledger closure**, including that cumulative deficit is
  non-decreasing: a backward step would mean a deficit had been repaid, which
  the model does not permit.
- **Demand-integral identities** — delivered plus unmet equals requested,
  per organ and in total.
- **Oxidation stoichiometry**, generalised to a per-organ glucose share so the
  identity holds in both fuel arms rather than only where all demand is
  glucose.
- **Store bounds**: no organ store exceeds capacity or falls below zero.

Ledger tolerances are 10⁻⁷ mL, far below any reported quantity, so these are
structural checks rather than accuracy checks.

## 4. Numerical settings are verified in pairs

Every hydraulic result is computed at a nominal and a refined tolerance setting
and accepted only if the two agree within a scaled error bound (Methods §4.2).
A run that fails the pair is retried once and, if it fails again, is **recorded
as a numerical failure rather than as a physiological result**. Across the
reported searches, numerical failures number zero and retries are recorded
where they occurred.

The screen horizon is itself verified: cycle-mean pressure declines toward
steady state after activation, so a run whose final five cycle means have not
flattened to within 0.05 mmHg is flagged as horizon-limited rather than
credited with its last value.

## 5. Searches are verified against a baseline that could beat them

A structural search argument cannot verify itself. The optimum located by
boundary tracing was therefore checked against a **matched-budget scrambled
Sobol sample** at equal evaluation cost, and a joint compass search was run
from a different direction. The sample **beat** the traced result by 2.7 per
cent, and the cause was identified as a rate axis traced at 60, 80 and 100
beats/min while the optimum lies near 49 to 57. We report this as a finding
about the method rather than adjusting the account of it.

Caching is verified too, because a cache that answers a different question is
indistinguishable from a result. Every evaluation is addressed by a key
carrying the design, model-form case, activation latency and screen horizon. A
defect in which the horizon was absent from that key — pairing 40-second
pressure verdicts with 10-second self-supply verdicts across 215 stored runs —
was found and is covered by a regression test.

## 6. Claims are audited against their evidence

Every numerical claim in the Results section maps to a saved run or a citation,
recorded with its kind, verification status and qualification, and the audit is
**executable**: it fails when a claim loses its evidence. The first pass over
48 claims found three citing paths that did not exist, a rate of one in
sixteen, and one set of figures that was correct but rested on an uncommitted
run cache and so was not reproducible from a clean clone.

The audit also enforces the discipline this subject requires: no modelled
proxy is restated as survival, neurological outcome or an injury dose; scenario
results are fractions of tested cases and never probabilities; and a model
output is never presented as a measurement in a subject.

## 7. The regression suite

475 automated tests cover the model, the protocol and the reported claims.
Tests were added specifically for defects found during the work, so that each
failure mode is caught rather than remembered: horizon-free cache keys, a
silently ignored activation-latency argument, a cycle-remainder rule that
voided valid runs, the omission of a liver-function factor from the duration
calculation, seam defects between the two fuel arms, selection of designs
below the declared pressure margin, and claims whose cited evidence has moved.

## 8. What verification does not establish

Verification here means internal consistency, reproducibility and traceability.
It does **not** establish that the model is an adequate representation of a
human in cardiac arrest. There is no validation against any measurement from
an arrested subject, because the configuration the study examines does not
exist. Every result is therefore a statement about a declared model under
declared assumptions, and the two named limitations above bound it further.
