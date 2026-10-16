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
from cardiac_redundancy.feasibility import FeasibilityPolicy

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
        version=2, task='T8.1', frozen_on='2026-10-16',
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
        pending_amendments=[
            dict(change='cycle-aligned pressure evaluation',
                 reason='docs/m2-pressure-requirements.md asks for the integral mean, maximum and '
                        'minimum over each auxiliary activation cycle, reserving the fixed '
                        '1 s window for irregular or nonpulsatile runs. The frozen policy uses '
                        '1 s windows throughout, which mixes 1.33 cycles at the reference period '
                        'and will mix differently once period_s varies',
                 status='must be settled before the main runs; it changes screen arithmetic and '
                        'would require rerunning every pressure verdict')],
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
        evaluation_policy=asdict(FeasibilityPolicy()),
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
        endpoints=['first constraint violation bracket', 'startup shortfall exposure',
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
