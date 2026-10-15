"""Static M2 specification checks. Does not integrate a cardiovascular model."""
from pathlib import Path
import json
import math
import hashlib
import csv

ROOT = Path(__file__).resolve().parents[1]

def audit():
    b = json.loads((ROOT/'configs/baseline.yaml').read_text())
    p = json.loads((ROOT/'configs/project-extension.yaml').read_text())
    manifest = json.loads((ROOT/'references/m2/source-manifest.json').read_text())
    for path, digest in manifest['local_files'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == digest, path
    assert abs(sum(b['initial_volumes_mL'].values())-1500) < 1e-9
    assert abs(sum(p['initial_physical_volumes_mL'].values())-5250) < 1e-9
    assert abs(sum(p['unstressed_offsets_mL'].values())+sum(x['unstressed_mL'] for x in p['limb_pools'].values())-3750) < 1e-9
    assert abs(sum(x['volume_initial_mL'] for x in p['limb_pools'].values())-600) < 1e-9
    assert all(v > 0 for v in p['initial_physical_volumes_mL'].values())
    assert abs(sum(x['baseline_fraction'] or 0 for x in p['branches'])-1) < 1e-12
    assert len({x['name'] for x in p['branches']}) == len(p['branches'])
    pressures = {'AO':100,'VC':2,'PU':4,'H':12}
    pressures.update({k:5 for k in p['limb_pools']})
    junction_net = {'H':0, **{k:0 for k in p['limb_pools']}}
    for edge in p['branches']:
        assert edge['R_mmHg_s_mL'] > 0
        q=(pressures[edge['from']]-pressures[edge['to']])/edge['R_mmHg_s_mL']
        assert abs(q-edge['target_flow_mL_s']) < 1e-10
        if edge['from'] in junction_net: junction_net[edge['from']] -= q
        if edge['to'] in junction_net: junction_net[edge['to']] += q
    assert max(abs(x) for x in junction_net.values()) < 1e-10
    # Check closed-loop topology and incidence conservation independently of flows.
    native=[('PU','LV'),('LV','AO'),('VC','RV'),('RV','PA'),('PA','PU')]
    edges=native+[(e['from'],e['to']) for e in p['branches']]
    nodes=set(sum(([a,z] for a,z in edges),[]))
    for origin in nodes:
        reached={origin}
        while True:
            after=reached|{z for a,z in edges if a in reached}
            if after==reached: break
            reached=after
        assert reached==nodes, (origin,reached)
    balances=dict.fromkeys(nodes,0.0)
    for i,(a,z) in enumerate(edges):
        q=(-1)**i*(i+1)/7
        balances[a]-=q;balances[z]+=q
    assert abs(sum(balances.values())) < 1e-10
    # Passive and active septal residual is strictly increasing: sample root checks.
    pars=b['source_parameters']
    def wall(comp,suffix,v,e):
        d=pars[comp]
        return e*d['E_es_'+suffix]*(v-d['V_d_'+suffix])+(1-e)*d['P_0_'+suffix]*math.expm1(d['lambda_'+suffix]*(v-d['V_0_'+suffix]))
    roots=[]
    for e in [0,.5,1]:
        def residual(v):
            return wall('septum','spt',v,e)-wall('lvf_calculator','lvf',94.6812-v,e)+wall('rvf_calculator','rvf',90.7302+v,e)
        lo,hi=-100,100
        assert residual(lo)<0<residual(hi)
        for _ in range(100):
            mid=(lo+hi)/2
            if residual(mid)>0: hi=mid
            else: lo=mid
        v=(lo+hi)/2;assert abs(residual(v))<1e-8
        roots.append({'activation':e,'septal_coordinate_mL':v,'residual_mmHg':residual(v)})
    assert b['pressure_unit_override'].startswith('mmHg')
    assert b['source_parameters']['heart_parameters']['V_tot']==5.5
    # The phase repeats at long horizons instead of terminating after 14 beats.
    d=b['driver'];T=d['period_s']
    for t in [.01,.2,.5]:
        def drive(x): return d['A']*math.exp(-d['B_per_s2']*((x % T)-d['C_s'])**2)
        assert abs(drive(t)-drive(t+100*T)) < 1e-12
    rows=list(csv.DictReader((ROOT/'data/m2-source-parameters.csv').open(encoding='utf-8')))
    assert len(rows)==len({(r['component'],r['symbol']) for r in rows})
    return {'status':'passed static specification checks only','source_parameter_rows':len(rows),
            'physical_volume_mL':5250,'limb_pool_volume_mL':600,'source_coordinate_volume_mL':1500,
            'regional_fraction_sum':1,'septal_checks':roots,
            'checks':['pinned hashes','volume accounting','positive resistances','anchor flow reconstruction',
                      'portal/pool continuity','complete native return paths','incidence conservation',
                      'passive/active septal roots','periodic driver beyond 14 beats','unique source rows'],
            'not_done':['ODE integration','waveform reproduction','physiological validation','support sizing']}

if __name__=='__main__':
    result=audit()
    (ROOT/'data/m2-specification-audit.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))
