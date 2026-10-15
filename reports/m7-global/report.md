# T7.2 global sensitivity: 12 trajectories

84 joint cases verified, with paired numerical solves. No failed cases omitted. Main-study retries: 0. Four-level Morris design, seed7202; 5000 whole-trajectory bootstrap replicates, seed7203. 274 regression tests pass.

Mu-star is mean absolute elementary effect per full input range. Sigma mixes nonlinearity and interactions. Bootstrap intervals quantify finite-trajectory sampling uncertainty under this artificial scenario design; they do not describe a human population.

## Ten-second results

| Output | Factor | Mu-star | 95% bootstrap interval | Sigma | Rank interval |
|---|---|---:|---|---:|---|
| total_deficit_mL | extraction_limit | 0.678239 | [0.414377, 0.939541] | 0.483835 | [1, 1] |
| total_deficit_mL | preservation_conductance | 0.269863 | [0.117308, 0.450051] | 0.307576 | [2, 4] |
| total_deficit_mL | coronary_conductance | 0.260883 | [0.166424, 0.356357] | 0.173709 | [2, 3] |
| total_deficit_mL | mass_ratio | 0.144487 | [0.0651402, 0.236982] | 0.16108 | [3, 6] |
| total_deficit_mL | contractility | 0.0904701 | [0.0452694, 0.142693] | 0.129012 | [4, 6] |
| total_deficit_mL | work_cost_multiplier | 0.0897377 | [0.0525634, 0.128863] | 0.0710374 | [4, 6] |
| MAP_mmHg | preservation_conductance | 13.6327 | [12.5664, 14.8361] | 2.04402 | [1, 1] |
| MAP_mmHg | contractility | 2.38229 | [1.91749, 2.911] | 0.914532 | [2, 2] |
| MAP_mmHg | coronary_conductance | 0.721884 | [0.65025, 0.795446] | 0.134171 | [3, 3] |
| aux_output_L_min | preservation_conductance | 0.727924 | [0.628576, 0.811383] | 0.168331 | [1, 1] |
| aux_output_L_min | contractility | 0.123413 | [0.0795207, 0.168482] | 0.0810672 | [2, 3] |
| aux_output_L_min | coronary_conductance | 0.0792377 | [0.0235194, 0.143575] | 0.124654 | [2, 3] |

Zero or numerically unresolved effects are omitted from ranking only, never from the simulations or full JSON results. Separate bootstrap intervals do not constitute a simultaneous ranking confidence set.

## Sampling adequacy

| Output | Previous prefix rank | Current rank | Relative change in current top two | Primary stability rule |
|---|---|---|---|---|
| total_deficit_mL | extraction_limit > preservation_conductance > coronary_conductance > mass_ratio > contractility > work_cost_multiplier | extraction_limit > preservation_conductance > coronary_conductance > mass_ratio > contractility > work_cost_multiplier | {'extraction_limit': 0.2449754429479172, 'preservation_conductance': 0.07431688261996651} | True |
| MAP_mmHg | preservation_conductance > contractility > coronary_conductance | preservation_conductance > contractility > coronary_conductance | {'preservation_conductance': 0.00946395402443868, 'contractility': 0.0718428112630986} | True |
| aux_output_L_min | preservation_conductance > contractility > coronary_conductance | preservation_conductance > contractility > coronary_conductance | {'preservation_conductance': 0.02812277521391386, 'contractility': 0.07792786610708229} | True |

Unstable registered primary outputs: none. Next registered stage: none. The rule compares the previous/current top-two sets and allows at most 30% change in their mu-star magnitudes. It is an exploratory stopping rule, not proof of asymptotic convergence. At sixteen trajectories any unresolved ranking is retained explicitly.

## Window and numerical checks

| Output | 4 s | 7 s | 10 s | Numerical rank match at 10 s |
|---|---|---|---|---|
| total_deficit_mL | extraction_limit > preservation_conductance > coronary_conductance > contractility > mass_ratio > work_cost_multiplier | extraction_limit > preservation_conductance > coronary_conductance > mass_ratio > contractility > work_cost_multiplier | extraction_limit > preservation_conductance > coronary_conductance > mass_ratio > contractility > work_cost_multiplier | True |
| MAP_mmHg | preservation_conductance > contractility > coronary_conductance | preservation_conductance > contractility > coronary_conductance | preservation_conductance > contractility > coronary_conductance | True |
| aux_output_L_min | contractility > preservation_conductance > coronary_conductance | preservation_conductance > contractility > coronary_conductance | preservation_conductance > contractility > coronary_conductance | True |

## Disjoint-half diagnostic

| Output | First half mu-star ranking | Second half mu-star ranking |
|---|---|---|
| total_deficit_mL | extraction_limit > preservation_conductance > mass_ratio > coronary_conductance > contractility > work_cost_multiplier | extraction_limit > coronary_conductance > preservation_conductance > work_cost_multiplier > contractility > mass_ratio |
| MAP_mmHg | preservation_conductance > contractility > coronary_conductance | preservation_conductance > contractility > coronary_conductance |
| aux_output_L_min | preservation_conductance > contractility > coronary_conductance | preservation_conductance > contractility > coronary_conductance |

Disjoint halves reveal dependence on sampled backgrounds. Rank disagreements remain visible even if the nested-prefix stopping rule passes.


## Scenario outcomes

Final-window MAP spans 39.4997–57.8211 mmHg; 0 of 84 points reach 65 mmHg. Total deficit spans 0.000366526–0.954152 mL O2 equivalent. Counts are neither feasibility nor survival probabilities.

Maximum additional oxygen limitation at saved samples is 0 mL/s. Cumulative ledgers remain primary; sampled attribution and quadrature errors are saved separately.

## Scope and limitations

These runs directly evaluate joint interior grid points and avoid the invalid corner surrogate. Six influential coordinates are varied together; anatomy groups remain linked and all other physiology stays fixed. The design does not estimate continuous-domain Sobol indices, optimized pump parameters, injury, reversible recovery or long-duration support. Variability among elementary effects is a reason to retain uncertainty, not to select favorable defaults. Wider physiology and alternative equations remain separate tasks. Contraction remains prescribed despite unmet metabolic demand: pressure/output sensitivities are conditional on that driver, not validated energy-limited pump performance.

The faster pilot was rejected after capacity checks failed; the main study uses original strict settings. All old frozen inputs must retain their hashes. Individual conservation, capacity, stoichiometry, event and refinement checks remain in case results.

Protocol: [joint global protocol](../../docs/m7-global-protocol.md). Methods: [Morris 1991](https://doi.org/10.1080/00401706.1991.10484804), [Campolongo et al. 2007](https://doi.org/10.1016/j.envsoft.2006.10.004), [SALib documentation](https://salib.readthedocs.io/en/latest/api/SALib.analyze.html#module-SALib.analyze.morris). Reproduce with scripts/run_m7_global.py followed by scripts/report_m7_global.py, passing the trajectory count.
