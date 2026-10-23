"""T8.1 freeze: write configs/m8-protocol.json from the validated code.

The registry, the evaluation policy and the declared bounds are read out of
src/cardiac_redundancy/design.py and feasibility.py rather than retyped, so a
protocol file cannot disagree with the validators the runs will use.

Also recorded: the code commit, library versions, the measured single-
evaluation cost that the search budget is derived from, and sha256 digests of
every input the protocol depends on. Amendments must bump `version` and state
what changed; reruns of affected analyses are required, not optional.
"""
from dataclasses import asdict, fields
from pathlib import Path
import hashlib, json, platform, subprocess, sys
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / 'src'))
import numpy, scipy
from cardiac_redundancy.design import BOUNDS, DesignPoint, ScenarioPoint, registry
from cardiac_redundancy.feasibility import CRITERIA, FeasibilityPolicy

#: The policy the main runs use. The dataclass defaults stay at the stricter
#: fixed-window strict-triple setting so code that does not read this file is
#: never silently relaxed; the declared values here are what T8.2 onward apply.
DECLARED_POLICY = FeasibilityPolicy(window_mode='cycle', criterion='map_primary', horizon_s=40.)

#: Minimum pressure margin a design must hold to be reported as feasible.
#: A hard threshold with a flat objective selects the most fragile point that
#: technically qualifies: the located optimum held +0.005 mmHg, and a 0.04 per
#: cent change in its activation period flipped it to violated. The margin is
#: declared here, before the re-run, rather than chosen after seeing which
#: designs survive.
MINIMUM_PRESSURE_MARGIN_MMHG = 3.0

OUT = ROOT / 'configs/m8-protocol.json'
POLICY_INPUTS = ('configs/baseline.yaml', 'configs/m3-evidence-regional.yaml',
    'configs/m2-frozen-overlay.yaml', 'configs/auxiliary-heart-design.yaml',
    'configs/m4-constraint-policy.json', 'configs/m4-oxygen-reference.json',
    'configs/m4-fuel-reference.json', 'configs/m4-organ-mapping.json',
    'configs/m5-myocardial-demand.json', 'configs/m5-final-comparison.json',
    'configs/m5-organ-preservation.json', 'configs/m6-controller.json',
    'configs/m6-coupled-protocol.json', 'data/m1-demand-budget.json',
    'src/cardiac_redundancy/design.py', 'src/cardiac_redundancy/feasibility.py',
    'reports/m5-joint-preservation/reference_g25/fine-refined.npz',
    'reports/m4-baseline/native-cycle.npz', 'reports/m4-baseline/store-5-refined.npz',
    'reports/m5-normal-comparison/verification.json')

#: Declared physiology triples. A search may not move between them; each is
#: evaluated as its own scenario. Adverse and favourable take every factor to
#: the same end of its declared interval, which is a deliberately conservative
#: and deliberately optimistic corner, not a probability statement.
SCENARIOS = dict(
    reference=dict(),
    adverse=dict(extraction_limit=.25, mass_coupling=1.6, work_cost_multiplier=1.5,
                 preservation_fraction=.5, glycogen_mmol_L_liver=150.,
                 release_cap_mmol_min=.1, liver_function_fraction=.25),
    favourable=dict(extraction_limit=.75, mass_coupling=.8, work_cost_multiplier=.75,
                    preservation_fraction=.1, glycogen_mmol_L_liver=400.,
                    release_cap_mmol_min=.4, liver_function_fraction=1.))

#: Model-form alternatives that cross every physiology scenario. These are not
#: design controls: a search may not choose them, because choosing them would
#: pick a favourable physiology instead of meeting the requirement. A candidate
#: must satisfy the screens across the declared cases, and where it cannot, the
#: result is reported as conditional on the case.
MODEL_FORM = dict(
    compression_kappa=dict(
        reference=.5, cases=[0., .5, 1.],
        acts_on='fraction of auxiliary chamber transmural pressure applied as the '
                'closing pressure of its own coronary bed (vascular waterfall)',
        provenance='configs/m2-frozen-overlay.yaml compression_kappa_cases; declared in M2 '
                   'and assigned to T7.4 by docs/uncertainty-plan.md, never exercised by any '
                   'runner before T8.2',
        measured_effect='at scale 0.6, conductance 1.0, reference physiology: auxiliary unmet '
                        '0.000000 mL at kappa 0, 0.002351 at kappa 0.5, 0.300505 at kappa 1'))

#: T8.2 operating cases, each a prescribed circumstance rather than a design
#: choice. Keys match the existing M6/M7 runner case flags.
OPERATING_CASES = dict(
    no_support=dict(policy='none', meaning='native pumping loss with no auxiliary release'),
    support_no_suppression=dict(policy='observed', preservation_override=1.,
        meaning='auxiliary support with organ demand left at baseline'),
    support_with_suppression=dict(policy='observed',
        meaning='auxiliary support with the declared suppressed-demand scenario'),
    oracle_activation=dict(policy='fixed', release_s=.5,
        meaning='prescribed release at the activation-timing latency; an upper reference'),
    delayed_activation=dict(policy='observed', meaning='observed detection latency',
        note='the observed policy released at 2.6 s in the reference run'),
    ventilation_degraded=dict(policy='observed', apnea=4.,
        meaning='gas-exchange loss during support'),
    partial_recovery=dict(policy='observed', restart=6., reopen=6.,
        meaning='native restart with regional reopening and handover'),
    shared_failure=dict(policy='observed', sensor_loss=3.,
        meaning='sensing loss during support; command is latched'),
    low_preload=dict(policy='observed', deplete_PU=True,
        meaning='reduced donor-pool filling'))

NUMERICAL = dict(
    hydraulic=dict(solver='Radau', pairs=[dict(label='nominal', rtol=1e-7, max_step_s=.01),
                                          dict(label='refined', rtol=1e-9, max_step_s=.005)],
                   sample_step_s=.005, paired_scaled_error_limit=.005,
                   conservation_limits=dict(blood_error_mL=1e-5, state_jump_mL=1e-10,
                                            root_residual_mmHg=1e-8, PV_quadrature_error_J=.05)),
    metabolic=dict(solver='RK45', pairs=[dict(label='nominal', rtol=1e-9, max_step_s=.005),
                                         dict(label='refined', rtol=1e-11, max_step_s=.0025)],
                   paired_state_limits=dict(oxygen_mL=.002, glucose_mmol=1e-5),
                   balance_limits=dict(oxygen_mL=1e-6, glucose_mmol=1e-8)),
    controller=dict(sample_period_s=.05, filter_window_s=.75, max_sample_gap_s=.1),
    failed_run_policy='a domain exit or failed paired check is recorded with its cause and '
                      'excluded from optimality claims; it is never counted as a pass, and '
                      'never silently resampled',
    tolerance_policy='tolerances are frozen here; they are not loosened to afford more '
                     'evaluations, which was tried and rejected in the M7 pilot')


def digest(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


def git(*args):
    try:
        return subprocess.run(['git', '-C', str(ROOT), *args], capture_output=True,
                              text=True, check=True).stdout.strip()
    except Exception:
        return None


def budget(evaluation_s, hydraulic_s, workers=12):
    """Search budget derived from the measured cost, not chosen first.

    The hydraulic stage alone decides the pressure screen, because metabolism
    is replayed over saved hydraulics and never feeds back into them. Pressure
    feasibility can therefore be screened at the hydraulic cost, with the
    coupled stage spent only on candidates that already pass.
    """
    stage1, stage2 = 3000, 600
    return dict(measured_single_evaluation_s=evaluation_s, measured_hydraulic_s=hydraulic_s,
        workers=workers,
        stage1=dict(what='pressure feasibility search on the hydraulic stage',
                    evaluations=stage1, wall_hours=round(stage1*hydraulic_s/workers/3600, 2)),
        stage2=dict(what='coupled evaluation of pressure-passing candidates across scenarios',
                    evaluations=stage2, wall_hours=round(stage2*evaluation_s/workers/3600, 2)),
        justification='staging is sound only because the pipeline is one-way; if supply-limited '
                      'contraction is ever implemented the stages must be merged again')


def main():
    timing = json.loads((ROOT / 'reports/m8-protocol/timing/single-evaluation.json').read_text(encoding='utf-8'))
    bounds_ok = {k: list(v[0]) for k, v in BOUNDS.items()}
    protocol = dict(
        version=5, task='T8.1', frozen_on='2026-10-22',
        amendments=[
            dict(version=2, change='auxiliary myocardial oxygen store sized by the shared '
                 'store_seconds=5 rule at the design demand, replacing the absolute pair in '
                 'configs/m5-myocardial-demand.json',
                 reason='the frozen pair is five seconds of a 7.0/3.0 mL/min demand, i.e. the '
                        'same rule evaluated once at the reference fixture. Held fixed, a larger '
                        'pump carries a smaller pump\'s store and the run aborts on "Oxygen '
                        'capacity exceeded" as soon as supply exceeds demand',
                 reruns='none required: scales 0.4 and 0.6 returned identical deficits, so the '
                        'constant only blocked evaluation. It also confirms the deficit is '
                        'glucose-limited, not oxygen-store-limited'),
            dict(version=2, change='the nine unresolved requirements from '
                 'configs/m4-constraint-policy.json are attached to every FeasibilityPolicy',
                 reason='the default was empty, so a clean numeric run would have been labelled '
                        'satisfied_within_declared_model instead of indeterminate',
                 reruns='none: no verdict had been issued before the fix'),
            dict(version=2, change='compression_kappa declared as a model-form scenario axis '
                 'with cases 0/0.5/1',
                 reason='it is the most influential quantity found in the self-supply constraint '
                        'and had never been varied by any runner, despite M2 declaring the cases '
                        'and the uncertainty registry assigning them to T7.4',
                 reruns='self-supply results obtained at kappa 0.5 are now conditional on that '
                        'case and are reported as such'),
            dict(version=2, change='search method: boundary tracing with bisection as primary, '
                 'differential evolution demoted to neighbourhood verification',
                 reason='the objective proved monotone in scale and rate, placing the optimum on '
                        'the constraint boundary; a global population search would spend '
                        'thousands of evaluations rediscovering that structure',
                 reruns='none: no search had been run'),
            dict(version=2, change='pressure_endpoint additionally returns per-window records, a '
                 'union violation count and the first fully satisfied window',
                 reason='a per-metric count was mistaken for the overall total when reporting',
                 reruns='none: reporting only, no screen arithmetic changed')],
        amendments_v3=[
            dict(version=3, change='pressure windows are cycle-aligned once the auxiliary pump '
                 'is released, one activation period per window, with the declared 1 s window '
                 'kept for the nonpulsatile pre-release interval',
                 reason='docs/m2-pressure-requirements.md asks for the integral mean, maximum '
                        'and minimum over each auxiliary activation cycle and reserves the fixed '
                        'window for irregular or nonpulsatile runs. A single 1 s grid mixes 1.33 '
                        'cycles at 80 beats/min and 2.5 at 150, so every comparison across rates '
                        'was measuring different things',
                 reruns='Phase A and Phase B must be re-run; the located optimum and the '
                        'U-shape in rate both depend on cross-rate comparisons'),
            dict(version=3, change='the enforced criterion is map_primary (mean AO pressure '
                 'against 65 mmHg); the strict triple is still evaluated and reported for every '
                 'run under status_by_criterion',
                 reason='AHA 2025 Part 11 states a minimum MAP of 65 mmHg for adults after '
                        'arrest and no adult systolic target, and four randomized trials found '
                        'no benefit from higher MAP targets. SBP 90 and DBP 50 derive from one '
                        'perioperative observational cohort using radial measurement and a '
                        'five-minute exposure definition, a different population and site. The '
                        'systolic threshold was binding both the minimum cavity scale and the '
                        'whole activation-latency budget, so which criterion is enforced had to '
                        'be settled on the evidence and declared before the runs, not after',
                 withdrawn_claim='an earlier draft stated that AHA 2025 removed a previous '
                                 'systolic target; that could not be confirmed from the sources '
                                 'and is withdrawn. What is supported is that the adult guidance '
                                 'is expressed as a MAP target and states no systolic one',
                 sources=['AHA 2025 Part 11 Post-Cardiac Arrest Care, doi 10.1161/CIR.0000000000001375',
                          'Adult Post-Cardiac Arrest Care Algorithm, AHA 2025',
                          'Ahuja 2020 PMID 31939844 (origin of SBP 90 / DBP 50)',
                          'Kjaergaard 2022 BOX trial PMID 36027564 (MAP 63 vs 77)'],
                 reruns='all pressure verdicts are recomputed; both criteria are reported so no '
                        'earlier result is discarded'),
            dict(version=3, change='supported duration is the primary endpoint',
                 reason='the research question is how long vital-organ perfusion can be '
                        'sustained, not whether one instant satisfies a screen. The screens '
                        'become the condition under which a duration is counted',
                 reruns='none: duration is computed from demands the existing runs already '
                        'report, and the analytic bound requires long-horizon verification')],
        minimum_pressure_margin_mmHg=MINIMUM_PRESSURE_MARGIN_MMHG,
        amendments_v5=[
            dict(version=5, change='a design must hold at least %.1f mmHg of pressure margin '
                 'to be reported as feasible' % MINIMUM_PRESSURE_MARGIN_MMHG,
                 reason='minimizing a flat objective against a hard threshold selects the most '
                        'fragile qualifying point. The located optimum held +0.005 mmHg, and a '
                        '0.0005 s change in its period, 0.04 per cent, flipped it to violated '
                        'with two bad windows. Nine of 152 feasible designs held under 0.5 mmHg '
                        'and the cheapest were all among them, so the search concentrated on '
                        'fragility by construction. This is a specification error rather than a '
                        'search error',
                 cost='a 2 mmHg margin costs 2.2 per cent more auxiliary oxygen demand and a '
                      '3 mmHg margin 3.1 per cent; the 3 mmHg design also runs at coronary '
                      'conductance 0.500, below nominal, so it depends on no conductance '
                      'assumption',
                 reruns='the reported optimum, and any disconfirmation based on it. The earlier '
                        'disconfirmation tested the four cheapest designs, all with margins '
                        'under 0.2 mmHg, so it measured the fragility the objective selected '
                        'rather than a property of the arm'),
            dict(version=5, change='disconfirmation samples a spread of the feasible set rather '
                 'than its cheapest members',
                 reason='two conclusions were corrected by cross-checks that found the opposite '
                        'outside the sampled corner: that fat admits a smaller pump, and that no '
                        'design survives compression kappa 0. Both followed from generalising '
                        'from the cheapest corner of a flat valley',
                 reruns='Arm A disconfirmation')],
        amendments_v4=[
            dict(version=4, change='screen horizon raised from 10 s to 40 s, with a plateau '
                 'check on the final five cycle means (0.05 mmHg spread) reported for every run',
                 reason='cycle-mean aortic pressure declines toward a steady state after '
                        'activation, so a horizon ending mid-transient credits a value the design '
                        'does not maintain. Measured at scale 0.375 and 80 beats/min: the tenth '
                        'cycle reads 65.42 mmHg and passes, while the design settles at 65.00 and '
                        'dips to 64.95 around cycles 21 to 60. Settling took 24 cycles, which is '
                        '36 s at the slowest declared rate, so 40 s covers the rate range',
                 effect_size='the pressure-scale slope is about 101 mmHg per unit cavity scale, '
                             'so a 0.45 mmHg steady-state shortfall moves the located boundary by '
                             'roughly 0.005 in scale, about 1 per cent. This is a verdict-labelling '
                             'correction at the margin, not a change of design conclusion',
                 cost='about 12 min per evaluation against 3 min at the ten-second horizon',
                 reruns='Phase A and Phase B boundary locations; the duration endpoint, the kappa '
                        'results, the coronary-bed scaling, the arterial volume-debt mechanism and '
                        'the latency slopes are all screen-independent and unaffected'),
            dict(version=4, change='the hydraulic cache key carries the horizon',
                 reason='without it a ten-second trajectory would silently answer a forty-second '
                        'request, the same class of defect as the earlier release-time key',
                 reruns='none')],
        pending_amendments=[
            dict(change='a declared lower limit on blood glucose, or duration reported as a '
                 'curve against that limit',
                 reason='the model declares no glucose floor, so duration means time to complete '
                        'exhaustion of the circulating pool: 33.6 min at the measured drain. '
                        'Against a 3.0 mmol/L plasma floor the same run gives 13.4 min. The '
                        'headline moves by a factor of three on a parameter never declared',
                 status='must be settled before any duration is reported as a result'),
            dict(change='long-horizon verification of the analytic duration bound',
                 reason='the 39.8 min reference figure assumes constant demand, hepatic export '
                        'already at its cap and no glucose floor; transport lags and startup '
                        'transients all shorten it',
                 status='required before any duration is reported as a result')],
        classification='frozen search and evaluation protocol; declared assumptions, not '
                       'measured limits, and no clinical or survival claim',
        code=dict(commit=git('rev-parse', 'HEAD'), branch=git('rev-parse', '--abbrev-ref', 'HEAD'),
                  dirty=bool(git('status', '--porcelain', '--', 'src', 'scripts', 'configs')),
                  python=platform.python_version(), numpy=numpy.__version__, scipy=scipy.__version__),
        objective=dict(
            primary='maximize the sustained-interval margin against the declared pressure screen, '
                    'then minimize auxiliary capacity and energy among candidates that satisfy it',
            aggregation='a candidate is reported on its worst declared scenario, not its mean',
            censoring='a run with no violation through the horizon is reported as '
                      '"no modeled violation observed through H", never as indefinite support',
            tie_breaking=['smaller cavity scale', 'lower chamber work', 'lower auxiliary demand'],
            forbidden=['choosing scenario factors to lower the requirement',
                       'crediting reduced consumption after unmet demand as efficiency',
                       'relabelling the screen to obtain a pass']),
        registry=registry(), declared_bounds=bounds_ok,
        reference_design=asdict(DesignPoint()), reference_scenario=asdict(ScenarioPoint()),
        scenarios=SCENARIOS, model_form=MODEL_FORM, operating_cases=OPERATING_CASES,
        evaluation_policy=asdict(DECLARED_POLICY),
        criteria={k: list(v) for k, v in CRITERIA.items()},
        constraints=dict(
            pressure='joint mean/maximum/minimum AO screen over complete 1 s windows',
            oxygen_active='no sustained unmet demand beyond numerical tolerance in the active group',
            oxygen_preservation='same screen on the preserved group',
            auxiliary_self_supply='same screen on the pump group; the auxiliary heart must '
                                  'supply its own myocardium',
            model_domain='positive cavity volumes and nonnegative vascular margins throughout',
            inventory='total blood volume conserved to the declared tolerance',
            unassessed=json.loads((ROOT / 'configs/m4-constraint-policy.json').read_text(
                encoding='utf-8'))['unresolved_requirements']),
        primary_endpoint='supported duration in minutes, with the declared screens as the '
                         'condition under which the interval is counted',
        endpoints=['supported duration (minutes)', 'first constraint violation bracket', 'startup shortfall exposure',
                   'sustained-interval status per constraint', 'chamber work and positive '
                   'transmural work', 'auxiliary oxygen demand and unmet ledger',
                   'glycogen and glucose remaining', 'cavity scale as a capacity proxy'],
        numerical=NUMERICAL,
        search=dict(
            objective='minimize total auxiliary oxygen demand (mL O2/min) computed from the '
                      'design and its own hydraulic work; chamber work reported alongside',
            formulation='robust constrained minimization: the design vector excludes every '
                        'scenario and model-form factor, and the screens must hold across the '
                        'declared cases, with the worst case reported',
            primary='boundary tracing with bisection, justified by three measured structural '
                    'properties: the objective rises monotonically in cavity scale and rate, so '
                    'the optimum lies on the pressure-feasibility boundary rather than inside '
                    'it; coronary conductance is nearly free energetically (11.624 to 11.589 '
                    'mL O2/min from 1.0 to 2.0) and only relaxes self-supply, so its optimal '
                    'value is the least that satisfies the constraint; and the constraints '
                    'nearly separate, since conductance and kappa move MAP by under 2 mmHg',
            monotonicity_premise='bisection assumes pressure feasibility is monotone in cavity '
                                 'scale at fixed rate. Observed at 80/min (0.25, 0.4 fail; 0.6, '
                                 '0.8, 1.0 pass) and 150/min (0.3, 0.4 fail; 0.5 pass). The '
                                 'premise is declared and must be checked by evaluating points '
                                 'either side of every located boundary, not assumed silently',
            secondary_objective='minimum coronary conductance, traced by the epsilon-constraint '
                                'method rather than folded into a weighted sum',
            verification='constrained differential evolution (scipy.optimize) in the '
                         'neighbourhood of the located optimum, plus a matched-budget '
                         'space-filling sample; required, not optional, since a structural '
                         'argument cannot be allowed to verify itself',
            seeds=[8101, 8102, 8103],
            fallbacks='mesh adaptive direct search, or Bayesian optimization with a separate '
                      'feasibility classifier, if the monotonicity premise fails',
            polish=False, polish_reason='objective is discontinuous at domain exits and '
                                        'screen boundaries',
            claim='best found within the tested bounds and budget; not a global optimum'),
        budget=budget(timing['evaluation_s'], timing['hydraulic_s']),
        inherited_limitations=['R-M3-01 source thoracic-pressure response',
            'L-M4-F01 no tissue glucose buffer; substrate use is instantaneous',
            'prescribed BIO contraction: unmet demand never reduces pumping',
            'supply-limited contraction and full tissue buffering unimplemented',
            'regional pressure minima, CO2/pH, injury, reperfusion and thrombosis unassessed',
            'cavity scale is a capacity proxy, not a determined anatomical mass'],
        inputs={p: digest(p) for p in POLICY_INPUTS},
        amendment_policy='bump version, record what changed and why, and rerun every analysis '
                         'that depended on the previous version')
    OUT.write_text(json.dumps(protocol, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    print(f'froze {OUT.relative_to(ROOT)} at commit {protocol["code"]["commit"]}')
    print(f'stage 1 {protocol["budget"]["stage1"]["wall_hours"]} h, '
          f'stage 2 {protocol["budget"]["stage2"]["wall_hours"]} h')
    return protocol


if __name__ == '__main__':
    main()
