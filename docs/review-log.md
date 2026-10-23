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
