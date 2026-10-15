"""T5.1 static capability characterization, not circuit optimization or sizing."""
from pathlib import Path
import csv,hashlib,itertools,json,sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from cardiac_redundancy.bounded_pump import PreloadGate,bounded_stage_flow
from cardiac_redundancy.backup_pump import PumpStage
from cardiac_redundancy.selective import finite_head_flow


def main():
    path=ROOT/'configs/m5-bounded-pump.json';cfg=json.loads(path.read_text());grid=cfg['grid'];rows=[]
    for side in ('right','left'):
        gate=PreloadGate(**cfg[side+'_gate']);resistance=cfg[side]['resistance_mmHg_s_mL']
        for pin,pout,head,capacity in itertools.product(grid['inlet_mmHg'],grid['outlet_mmHg'],grid['head_mmHg'],grid['capacity_mL_s']):
            stage=PumpStage(resistance,capacity);q=bounded_stage_flow(pin,pout,head,stage,gate)
            ideal=finite_head_flow(pin,pout,head,resistance,capacity)
            assert 0<=q<=capacity*gate.fraction(pin) and q<=ideal
            rows.append(dict(side=side,inlet_mmHg=pin,outlet_mmHg=pout,head_mmHg=head,capacity_mL_s=capacity,
                             preload_fraction=gate.fraction(pin),ideal_flow_mL_s=ideal,bounded_flow_mL_s=q))
    out=ROOT/'reports/m5-bounded-pump';out.mkdir(parents=True,exist_ok=True)
    with (out/'capabilities.csv').open('w',newline='',encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    changed=sum(r['bounded_flow_mL_s']<r['ideal_flow_mL_s'] for r in rows)
    result=dict(status='verified_static_characterization',points=len(rows),preload_limited_points=changed,
        dynamic_circulation='not yet characterized',physical_sizing='not inferred',
        hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                (path,Path(__file__),ROOT/'src/cardiac_redundancy/bounded_pump.py',ROOT/'src/cardiac_redundancy/backup_pump.py',ROOT/'src/cardiac_redundancy/selective.py')})
    (out/'result.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    (out/'report.md').write_text(f'''# T5.1 first increment: bounded inlet capability

2026-10-09. **T5.1 started, not closed.** {len(rows)} static stage operating points checked; {changed} have lower flow than the ungated comparator because of the assumed preload gate. [Raw grid](capabilities.csv), [configuration](../../configs/m5-bounded-pump.json), [provenance](result.json), [M5 plan](../../docs/m5-plan.md).

Let x = inlet pressure minus the configured external pressure. The assumed gate is g = clamp((x - p_closed)/(p_full - p_closed), 0, 1). Then Q = min(g × capacity, max(0, (P_in + head - P_out)/R)). Inlet availability limits capacity; it does not impose an output flow.

The default 0/4 mmHg thresholds are **A characterization fixtures**, not physiological suction or venous-collapse limits. Fully open operation exactly matches the M3 port law. At inlet 2 mmHg with 120 mL/s nominal capacity, the assumed available capacity is 60 mL/s; sufficient head still cannot exceed it. At/below closure pressure, additional head cannot produce flow. A filled enabled port may pass passive forward flow at zero head; isolation remains explicit.

The B1 adapter preserves the existing physical state, incidence matrix and regional controls. Tests verify exact fully-open RHS regression, independent stage limits and global blood conservation. These checks do not demonstrate balanced operation over time: pulmonary/systemic inventory shifts and the nutritive lung return require the next dynamic sweep.

This gate grants an ideal actuator and check valve, with prescribed head and no electrical/metabolic actuator dynamics. It does not describe a muscle wall, septum, chamber volume, maximum contractility, blood damage or efficiency. The M3 Smith-shaped fixture remains the separate biological starting point.

Use actual gated flow for energy accounting. The historical `pump_energy.pump_power` helper reconstructs ungated flow and must not be applied directly to this subclass; extending the energy interface is part of T5.2. No new input-energy or tissue-mass result is claimed.

Next: paired dynamic B1 stage sweeps for preload, unequal support and reservoir limits, followed by Smith-shaped chamber/rate controls. Fuel screening must carry M4's explicit L-M4-F01 limitation and every required preservation tissue.

Reproduce: `.venv\\Scripts\\python.exe scripts/characterize_m5_pump.py`.
''',encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
