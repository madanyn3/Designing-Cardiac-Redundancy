"""Run the first project-extension gate: physical inventory and coordinate identity."""
import hashlib
import json
from pathlib import Path
import platform
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
import numpy as np
import scipy
from cardiac_redundancy.inventory import PhysicalInventory, PhysicalSourceReference
from cardiac_redundancy.integration import integrate_source


def main():
    out=ROOT/'reports/m3-physical-inventory'
    out.mkdir(parents=True,exist_ok=True)
    source,project=ROOT/'configs/baseline.yaml',ROOT/'configs/project-extension.yaml'
    inventory=PhysicalInventory(source,project)
    model=PhysicalSourceReference(source,project)
    a=integrate_source(model.source,1.5)
    b=integrate_source(model,1.5)
    coordinate_error=np.max(abs(b.states-model.translation[:,None]-a.states),axis=1)
    event_labels_equal=[(e['valve'],e['action']) for e in a.events]==[(e['valve'],e['action']) for e in b.events]
    event_error=max(abs(x['time_s']-y['time_s']) for x,y in zip(a.events,b.events)) if event_labels_equal else None
    np.savez_compressed(out/'translation-trajectories.npz',time_s=a.time,
                        source_states=a.states,physical_states=b.states,translation=model.translation)
    frozen=json.loads((ROOT/'data/m2-freeze-manifest.json').read_text())['sha256']
    frozen_ok=all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in frozen.items())
    result={'scope':'Step 1 physical accounting and pre-split coordinate verification; no regional circulation yet',
        'total_physical_mL':float(inventory.initial_volumes.sum()),
        'total_unstressed_mL':float(inventory.unstressed_volumes.sum()),
        'central_coordinate_mL':float(inventory.coordinate_volumes(inventory.initial_volumes)[:6].sum()),
        'limb_physical_mL':float(inventory.initial_volumes[6:].sum()),
        'limb_coordinate_mL':float(inventory.coordinate_volumes(inventory.initial_volumes)[6:].sum()),
        'initial_vascular_pressure_mmHg':inventory.vascular_pressures(inventory.initial_volumes),
        'translation_error_by_state':coordinate_error.tolist(),
        'translation_absolute_limit_mL_or_mL_s':2e-4,
        'valve_sequence_equal':event_labels_equal,'maximum_valve_time_difference_s':event_error,
        'valve_time_limit_s':1e-5,
        'relative_physical_volume_drift':float(np.max(abs(b.states[:6].sum(axis=0)-5250))/5250),
        'frozen_nine_files_unchanged':frozen_ok}
    result['passed']=bool(max(coordinate_error)<2e-4 and event_labels_equal and event_error<1e-5
                          and result['relative_physical_volume_drift']<1e-8 and frozen_ok)
    paths=[source,project,ROOT/'configs/m2-frozen-overlay.yaml',Path(__file__),
           *sorted((ROOT/'src/cardiac_redundancy').glob('*.py'))]
    manifest={'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,
        'duration_s':1.5,'sample_step_s':.001,'solver':'Radau','settings':model.settings,
        'hashes':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    for p in paths:
        target=out/'implementation-snapshot'/p.relative_to(ROOT)
        target.parent.mkdir(parents=True,exist_ok=True); target.write_bytes(p.read_bytes())
    for name,data in [('result.json',result),('manifest.json',manifest)]:
        (out/name).write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))
    return 0 if result['passed'] else 1


if __name__=='__main__': raise SystemExit(main())
