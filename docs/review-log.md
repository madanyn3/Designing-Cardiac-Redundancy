# Review log

Records audit passes over preprint text. An entry states what was checked, what
was found, and what changed — not that a check happened.

---

## 2026-10-23 — Claim audit of the Results section (M8-7.2 gate, first pass)

**Scope.** Every numerical claim in [`preprint/results.md`](preprint/results.md),
48 claims, recorded in
[`preprint/results-claims.csv`](preprint/results-claims.csv) with its source
file or citation, kind, verification status and qualification. The audit is kept
executable as `tests/test_results_claims.py`, so a claim that loses its evidence
fails the suite rather than surviving in the text.

**Composition.** 26 model outputs, 17 derived quantities, 4 cited to literature,
1 declared assumption. No claim is unsourced.

**Found and fixed.**

1. **Three claims cited paths that do not exist** — the boundary trace, the
   earlier 84-case result and the cross-arm feasibility check. All three
   existed under different names
   (`phase-a-40s-kappa0.5.json`, `m7-global/summary-12.json`, and a folder whose
   period is written to six decimal places). The claims were correct; their
   citations were not. This is the defect the audit exists to catch, and it
   appeared in the first pass at a rate of 1 in 16.

2. **The latency numbers rested on an uncommitted run cache.** The re-measured
   relation was read directly from solve folders under `reports/m8-optimize/`,
   which are not committed, so the reported figures were correct but not
   reproducible from a clean clone. Collected into
   `reports/m8-latency/latency-margin.json` by
   `scripts/collect_m8_latency_margin.py`, which solves nothing and reports a
   missing solve as missing. Three claims whose only evidence is still a cache
   path now name the script that regenerates it, and a test enforces that.

3. **The latency conclusion was wrong as stated, and the artifact corrected
   it.** The text said the detection budget is "under one second". That is true
   of the cost-optimal design, which holds its margin to 0.81 s, but the
   nominal-conductance design holds to 1.57 s. Collecting all five designs into
   one artifact made the comparison visible. The claim is now stated per design,
   0.8 to 1.6 s, and it strengthened a separate result: the design that escapes
   the coronary conductance assumption also tolerates roughly twice the
   detection delay, so those two advantages coincide.

**Discipline checks, all passing as tests.** No modelled proxy is converted into
survival, neurological outcome or injury dose. Scenario results are stated as
fractions of tested cases and never as probabilities. Model outputs are not
presented as measurements in a subject — including the 539 mL/min hepatic
outflow and the 0.105 mmol/min export rate, both of which read as measurements
and are explicitly qualified. Species and conditions remain attached to
transferred coefficients.

**Not done in this pass.** The Methods and Verification sections have not been
audited. R-M3-01 and L-M4-F01 are carried in the limitations but their statement
in the Verification section is unverified. External review has not happened;
nothing in the text implies that it has.

---

## 2026-10-23 — T8.2 matrix re-run, and a cache key that omitted the design

**Why.** The nine-case operating matrix had been run at a design holding
+2.26 mmHg of pressure margin, below the 3.0 mmHg minimum declared in protocol
v5. Its conclusion concerns which window the release lands in rather than
sustained pressure, so it was not expected to move; it was re-run because a
reported matrix should stand on a reportable design.

**A defect found while doing it.** The matrix cache key was built from the case
name and the compression case alone, omitting the design. Re-pointing the
matrix at a different design therefore returned the previous design's runs:
every MAP came back byte-identical to four decimal places across nine cases,
and nothing reported that the requested design had never been solved. The run
wrote a summary labelled with the new design whose rows were the old design's
results.

This is the second instance of the same class. The first was a key missing the
screen horizon, which paired 40-second pressure verdicts with 10-second
self-supply verdicts across 215 stored runs. The pattern is that **a cache keyed
on less than the full question silently answers a different one**, and in both
cases the symptom was results that looked plausible rather than an error.

Recording the design in the artifact is not sufficient, since the mislabelled
summary did exactly that. `tests/test_m8_matrix_key.py` therefore checks four
things: that two declared designs cannot collide in the key, that the artifact
names its design, that the design holds the declared margin, and that each
row's reported MAP **matches the solve in the folder the design names**. The
last of these failed against the mislabelled artifact, which is how it was
confirmed rather than assumed.

**The re-run strengthened the conclusion.** At the margin-constrained design,
**five of nine cases satisfy every organ screen including auxiliary
self-supply**, against one of nine at the superseded design. In four of those
five the auxiliary feeds its own myocardium, every organ group is satisfied,
and mean pressure sits at 68.69 to 104.35 mmHg — and the case fails only
because release at 2.50 s puts five windows outside the declared two-second
startup allowance. In those four cases nothing about the design fails.

At the superseded design the same cases also failed the pump group, which left
room to read the design as inadequate. That reading is now closed off: what
fails is the arrival time of support, not the support. Two claims in the
Results ledger were updated accordingly, and the earlier "one of nine" organ
figure is recorded as superseded rather than deleted.
