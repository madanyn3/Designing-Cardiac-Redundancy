# M3 consolidated baseline and project verification report

Updated 2026-10-07. **Engineering verification is complete; the original journal-reproduction gate is not fully passed.** The user explicitly accepted R-M3-01 on 2026-10-07, closing M3 at verified-reconstruction scope while preserving the failed original screen.

Read the [M3 closure report and M4 handoff](../docs/m3-closure.md), [final evidence audit](m3-final/report.md), [source disposition](../docs/m3-source-final-disposition.md), and [anatomical-assumption register](../docs/m3-anatomical-assumptions.md). The [earlier implementation-progress report](m3-baseline-progress-archive.md) is retained as history.

## Source reproduction

| Quantity | Result | Original screen |
|---|---:|---|
| LV stroke volume | 65.9992 mL | Pass: 70 +/- 7 mL |
| Cardiac output | 5.27993 L/min | Pass: 5.6 +/- .56 L/min |
| Aortic systolic/diastolic | 116.502 / 79.044 mmHg | Pass: 120/80 +/- 10 mmHg |
| Pulmonary systolic/diastolic | 25.676 / 7.808 mmHg | Pass: 25/8 +/- 5 mmHg |
| Integral aortic MAP | 97.776 mmHg | Reported output, not a new adequacy threshold |
| Thoracic-pressure CO reduction | 19.078449% | **Fail:** 9 +/- 5 percentage points |
| Independent reconstruction | Maximum metric difference 0.003206% | Pass: <1% |

The [source sweeps](m3-sweeps/report.md) include nominal/refined Radau, BDF, repeat, output-grid and initial-state comparisons. All declared numerical comparisons pass. The [independent CVSX execution](m3-source-audit/report.md) confirms the reconstructed dynamics, including the disputed response. It does not identify the original Figure 8 experiment. No source parameter or acceptance tolerance was changed to force agreement.

## Project extensions

Verified components include physical/unstressed volume accounting; evidence-amended regional flows; cerebral and coronary closures; native arrest; independent regional conductance controls; conserved serial hepatic flow; reversible storage/routing/marker fixtures; configurable support; mechanical energy accounting; prescribed native restart and support withdrawal; delayed/unequal support; explicit model-domain exits; E1/E2 hydraulic topologies; and Smith-shaped auxiliary mechanics with independent periodic warm-up.

[The combined report](m3-circuits/report.md) covers 27 runs: 21 complete and six correctly localize domain stops. Its maximum scaled trajectory difference is 0.010049%, valve timing difference 0.002868 ms and domain-stop timing difference 0.041486 ms. The quarter-scale biological fixture does not maintain aortic pressure without native restart; numerical implementation success is not sufficient support.

[Final storage sensitivity](m3-storage-sensitivity/report.md) adds nine runs for 300/600/900 mL initial limb inventories with correspondingly scaled U/C, fixed total physical blood and independently warmed baselines. All comparisons pass. Maximum numerical difference is 0.008819%; retained limb blood changes substantially with assumed storage, while the tested fixed-head flow changes are modest. This is not an anatomical confidence interval.

## Verification and reproducibility

The fresh [automated suite](m3-final-tests.xml) has 102 passing tests. All three static M2 audits pass, with nine frozen hashes unchanged. The [machine-readable audit](m3-final/closure-audit.json) inventories separate pass/fail outcomes and evidence hashes. Historic ODE runs remain archived rather than being relabeled as newly executed. Scripts, effective configurations, snapshots, solver settings, events and trajectories are linked from each experiment report.

Regenerate the consolidated audit with `.venv/Scripts/python.exe scripts/audit_m3_closure.py`. The closure packet links the experiment-specific commands. Full fresh-machine environment reconstruction has not been performed.

## Scope and closure

Regional storage, anatomical isolation/collaterals, brain spatial supply, gas exchange, hemocompatibility, tissue recovery and finite metabolic fuel remain explicit limitations or later milestones. Mechanical work is not total input energy or liver glycogen consumption. No recoverable support duration or optimized pump is established.

The original reproduction gate remains failed. M3 is closed with the user-accepted R-M3-01 limitation described in [the closure report](../docs/m3-closure.md). M4 is next; its implementation remains unstarted.
