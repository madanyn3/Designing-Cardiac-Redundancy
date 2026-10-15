"""Read-only reanalysis of historical selective runs against the new objective."""
import hashlib,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'reports/m3-repair-window-review'


def main():
    rows=[];hashes={}
    for name in ('open-nominal','excluded-nominal','trapped-nominal'):
        p=ROOT/'reports/m3-selective'/name
        r=json.loads((p/'result.json').read_text());d=np.load(p/'segment-1.npz')
        t=d['time_s']-r['arrest_s'];mask=(t>=2)&(t<=4+1e-9)
        def average(x):return float(np.trapezoid(x[mask],t[mask])/(t[mask][-1]-t[mask][0]))
        q={n:d['flows'][:,i] for i,n in enumerate(r['flow_names'])}
        pr={n:d['pressures'][:,i] for i,n in enumerate(r['volume_names']+['H'])}
        rows.append(dict(case=name,window_s=[float(t[mask][0]),float(t[mask][-1])],
            right_support_L_min=average(q['support_right'])*.06,
            left_support_L_min=average(q['support_left'])*.06,
            total_net_port_hydraulic_power_W=average((pr['PA']-pr['VC'])*q['support_right']+(pr['AO']-pr['PU'])*q['support_left'])*.000133322,
            kidney_pair_mL_min=average(q['kidney_protected']+q['kidney_other'])*60))
        for f in (p/'result.json',p/'segment-1.npz'):
            hashes[str(f.relative_to(ROOT))]=hashlib.sha256(f.read_bytes()).hexdigest()
    OUT.mkdir(parents=True,exist_ok=True)
    result=dict(scope='Transient 2–4 s averages, not periodic support or design minima. Net hydraulic power at fluid ports; excludes losses, drive power and tissue metabolism.',
                conversion='1 mmHg*mL/s = 0.000133322 W',cases=rows,input_hashes=hashes)
    (OUT/'metrics.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
