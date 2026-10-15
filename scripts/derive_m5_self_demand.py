"""Derive declared demand scenarios from fine BIO work and normal reference."""
from pathlib import Path
import json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from cardiac_redundancy.biological_pump import BiologicalPumpModel,BiologicalParameters
from cardiac_redundancy.myocardial_demand import DemandAssumptions,side_demand
from cardiac_redundancy.pump_energy import MMHG_ML_TO_J as K

def main():
    read=lambda p:json.loads((ROOT/p).read_text(encoding='utf-8'))
    cfg=read('configs/m5-myocardial-demand.json');normal=read('reports/m5-normal-comparison/verification.json')['runs'][-1]['metrics']['cycle']
    native=read('data/m1-demand-budget.json')['baseline_mL_min']['heart'];weights=read('configs/m4-organ-mapping.json')['mapping']['heart'];results=[]
    for name in ('reference','faster','stronger_left','asymmetric'):
        params=read(f'reports/m5-biological/{name}/result.json')['parameters'];a=np.load(ROOT/f'reports/m5-biological-energy/{name}-refined.npz');t=a['time_s'];y=a['states'];p=a['pressures_mmHg']
        m=BiologicalPumpModel(ROOT/'configs/baseline.yaml',ROOT/'configs/m3-evidence-regional.yaml',ROOT/'configs/m2-frozen-overlay.yaml',biological=BiologicalParameters(**params),arterial_fraction=.2,inlet_resistance_fraction=.1)
        pars=m.native.base.source.parameters['pericardium'];external=m.native.base.source.config['pressure_external_mmHg']['pericardium']
        peri=external+pars['P_0_pcd']*np.expm1(pars['lambda_pcd']*((y[16]+y[17])/params['pericardial_scale']-pars['V_0_pcd']))
        variants={'reference':cfg['reference']}
        if name=='reference':variants.update({k:cfg['reference']|v for k,v in cfg['sensitivities'].items()})
        for variant,assumptions in variants.items():
            sides={}
            for side,vi,qi in [('left',16,22),('right',17,24)]:
                positive=float(np.trapezoid(np.maximum(K*(p[:,vi]-peri)*(y[qi+1]-y[qi]),0),t)/(t[-1]-t[0]))
                c=max(params[side+'_contractility'],params['septal_contractility'])
                sides[side]=side_demand(native*weights['coronary_'+side],positive,normal['sides'][side]['positive_transmural_mean_W'],60/params['period_s'],c,DemandAssumptions(**assumptions),cfg['normal_rate_bpm'])
                sides[side]['positive_work_W']=positive
            results.append(dict(name='BIO_'+name+'_'+variant,architecture='BIO',trajectory=name,assumptions=assumptions,sides=sides,auxiliary_mL_min=sum(x['total_mL_min'] for x in sides.values()),auxiliary_fraction_normal=sum(x['total_mL_min'] for x in sides.values())/native))
    out=ROOT/'reports/m5-self-demand';out.mkdir(parents=True,exist_ok=True)
    (out/'scenarios.json').write_text(json.dumps(dict(classification=cfg['classification'],normal_heart_O2_mL_min=native,scenarios=results),indent=2)+'\n')
    for r in results:print(r['name'],r['auxiliary_mL_min'])
if __name__=='__main__':main()
