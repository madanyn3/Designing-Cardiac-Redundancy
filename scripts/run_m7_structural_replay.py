"""Paired conservative replay of extraction and demand equation alternatives.

Hydraulics are fixed. This experiment cannot change MAP or repair the prescribed
pump's lack of energetic feedback. Dependency substitutions are process-local.
"""
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
import run_m7_global_metabolism as replay
from cardiac_redundancy.structural_closures import flow_extraction, transition_fraction

OUT = ROOT/'reports/m7-structural-replay'
BASE = ROOT/'reports/m7-factorial/hydraulic_nominal'
OriginalRespiration = replay.BilateralRespiration
OriginalGlucose = replay.GlucoseTransport
OriginalPreservation = replay.OrganPreservation


def experiment(spec):
    case = json.loads((BASE/'result.json').read_text())['case']
    case.update(name=spec['name'], hydraulic_source=str(BASE.relative_to(ROOT)), structural_alternative=spec)
    cycle = np.load(BASE/'hydraulic-refined.npz')
    times = cycle['time_s']
    mask = times >= times[-1]-2
    qref = np.trapezoid(np.abs(cycle['flows_mL_s'][mask]), times[mask], axis=0)/(times[-1]-times[mask][0])

    class AlternativeRespiration(OriginalRespiration):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.indices = [[e[0] for e in self.edges].index(o.edge) for o in self.organs]
            self.reference_organs = self.organs
        def evaluate(self, volumes, state, flows, **kwargs):
            if not spec.get('oxygen_flow'):
                return super().evaluate(volumes, state, flows, **kwargs)
            fractions = flow_extraction(np.asarray(flows)[self.indices], qref[self.indices], .75)
            self.organs = tuple(replace(o, extraction_limit=float(e)) for o,e in zip(self.reference_organs, fractions))
            try:
                return super().evaluate(volumes, state, flows, **kwargs)
            finally:
                self.organs = self.reference_organs

    class AlternativeGlucose(OriginalGlucose):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.indices = [[e[0] for e in self.edges].index(o.edge) for o in self.organs]
            self.reference_extraction = self.extraction.copy()
        def evaluate(self, volumes, state, flows, *args, **kwargs):
            if spec.get('glucose_flow'):
                self.extraction = flow_extraction(np.asarray(flows)[self.indices], qref[self.indices], self.reference_extraction)
            return super().evaluate(volumes, state, flows, *args, **kwargs)

    class AlternativePreservation(OriginalPreservation):
        def demand(self, time_s):
            form = spec.get('transition', 'step')
            if form == 'step':
                return super().demand(time_s)
            return self.reference*np.array([transition_fraction(time_s, e['fraction'], 5., form) for e in self.entries])

    replay.OUT = OUT
    replay.BilateralRespiration = AlternativeRespiration
    replay.GlucoseTransport = AlternativeGlucose
    replay.OrganPreservation = AlternativePreservation
    result = replay.run(case)
    result['structural_metadata'] = dict(reference_flow_mL_s=qref.tolist(),
        reference_flow_definition='mean absolute edge flow over final two seconds of frozen nominal hydraulics',
        demand_form=spec.get('transition', 'step'), duration_or_time_constant_s=5. if spec.get('transition') else 0.,
        profile_interpretation='profile fractions/roles retained; structural demand form overrides stored profile ramp_s; auxiliary demand remains workload/activation based')
    replay.write(OUT/case['name']/'result.json', result)
    return result


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    specs = [dict(name='reference'), dict(name='oxygen_flow', oxygen_flow=True),
             dict(name='glucose_flow', glucose_flow=True),
             dict(name='linear_5s', transition='linear'),
             dict(name='exponential_5s', transition='exponential'),
             dict(name='combined', oxygen_flow=True, glucose_flow=True, transition='exponential')]
    paths = [Path(__file__), ROOT/'scripts/run_m7_global_metabolism.py',
             ROOT/'src/cardiac_redundancy/structural_closures.py',
             ROOT/'docs/m7-structural-replay-protocol.md',
             ROOT/'reports/m7-factorial/manifest.json', BASE/'result.json', BASE/'hydraulic-refined.npz']
    manifest = dict(specs=specs, hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    manifest_path = OUT/'input-manifest.json'
    if manifest_path.exists() and json.loads(manifest_path.read_text()) != manifest:
        raise RuntimeError('Structural replay inputs changed')
    replay.write(manifest_path, manifest)
    with ProcessPoolExecutor(max_workers=3) as pool:
        for result in pool.map(experiment, specs):
            print(result['case']['name'], result['status'], flush=True)


if __name__ == '__main__':
    main()
