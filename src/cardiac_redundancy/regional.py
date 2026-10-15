"""Step-2 regional hydraulic comparator: physical reservoirs and signed resistors.

Uses M2 topology and kidney/coronary splits. Brain and coronary beds are the
declared simple-resistor comparators; pressure-dependent closures are step 3.
No organ oxygen demand, support pump or physiological flow floor is imposed.
"""
from dataclasses import dataclass
import json
from pathlib import Path

import numpy as np

from .inventory import PhysicalInventory
from .source import SourceModel, EDGES, VALVES
from .mechanics import ModelDomainError


@dataclass(frozen=True)
class Branch:
    name: str
    donor: str
    recipient: str
    resistance: float


class RegionalModel:
    """State: eleven physical volumes, then four native valve flows (15 total).

    H is an algebraic hepatic junction, with zero storage. Gut and hepatic artery
    enter H; hepatic outflow enters VC. H is never counted as extra blood or as
    an independent additional contribution to systemic output.
    """
    def __init__(self, source_path, project_path, overlay_path):
        self.source=SourceModel(source_path)
        self.inventory=PhysicalInventory(source_path,project_path)
        self.config=self.source.config
        self.settings=self.source.settings
        self.project=json.loads(Path(project_path).read_text(encoding='utf-8-sig'))
        self.overlay=json.loads(Path(overlay_path).read_text(encoding='utf-8-sig'))
        self.volume_names=self.inventory.names
        self.volume_count=len(self.volume_names)
        self.initial_state=np.r_[self.inventory.initial_volumes,self.source.initial_state[6:]]
        branches=[]
        for row in self.project['branches']:
            splits=[(row['name'],1.)]
            if row['name']=='kidneys':
                splits=list(zip(self.overlay['kidneys']['replace_pair_with'],self.overlay['kidneys']['baseline_weights']))
            elif row['name']=='heart':
                c=self.overlay['coronary']
                splits=[('coronary_left',c['left_conductance_fraction']),('coronary_right',c['right_conductance_fraction'])]
            if not np.isclose(sum(w for _,w in splits),1) or any(w<=0 for _,w in splits):
                raise ValueError('Invalid split conductances')
            for name,weight in splits:
                r=row['R_mmHg_s_mL']/weight
                if not np.isfinite(r) or r<=0: raise ValueError('Invalid regional resistance')
                branches.append(Branch(name,row['from'],row['to'],r))
        self.branches=tuple(branches)
        # Every branch touching the zero-storage hepatic junction participates
        # in its pressure solve. This supports both the historical gut aggregate
        # and the evidence-amended set of independently perfused portal organs.
        self.hepatic_branches=tuple(b for b in branches if 'H' in (b.donor,b.recipient))
        if not self.hepatic_branches or any(
                (b.recipient if b.donor=='H' else b.donor) not in self.volume_names
                for b in self.hepatic_branches):
            raise ValueError('H must connect to physical pressure compartments')
        self.branch_by_name={b.name:b for b in branches}
        if len(self.branch_by_name)!=len(branches): raise ValueError('Duplicate regional branch names')
        self.flow_names=VALVES+('pulmonary_exchange',)+tuple(b.name for b in branches)
        endpoints=EDGES[:4]+(('PA','PU'),)+tuple((b.donor,b.recipient) for b in branches)
        nodes=self.volume_names+('H',)
        self.incidence=np.zeros((len(nodes),len(endpoints)))
        for k,(a,b) in enumerate(endpoints):
            self.incidence[nodes.index(a),k]=-1
            self.incidence[nodes.index(b),k]=1
        self.incidence.flags.writeable=False

    def domain_margins(self,state):
        return self.inventory.domain_margins(np.asarray(state)[:self.volume_count])

    def pressures(self,t,state,active=True):
        state=np.asarray(state)
        if state.shape!=(self.volume_count+4,) or not np.all(np.isfinite(state)):
            raise ModelDomainError('Invalid regional state')
        # Reuse source cardiac wall/septal/pericardial mechanics without changing
        # LV/RV physical cavity coordinates. The temporary central vascular
        # pressures are all replaced by the physical V-U laws immediately below.
        p,m=self.source.pressures(t,np.r_[state[:6],state[self.volume_count:]],active=active)
        p.update(self.inventory.vascular_pressures(state[:self.volume_count]))
        # Kirchhoff balance at zero-storage H, valid for signed/reversed flows.
        conductance=sum(1/b.resistance for b in self.hepatic_branches)
        p['H']=sum(p[b.recipient if b.donor=='H' else b.donor]/b.resistance
                   for b in self.hepatic_branches)/conductance
        return p,m

    def flows(self,t,state,open_valves,active=True,pressures=None):
        p=self.pressures(t,state,active=active)[0] if pressures is None else pressures
        modes=np.asarray(open_valves,dtype=bool)
        if modes.shape!=(4,): raise ValueError('Four valve modes required')
        valve=np.where(modes,np.asarray(state)[self.volume_count:],0.)
        qp=(p['PA']-p['PU'])/self.source.parameters['heart_parameters']['R_pul']
        q=np.r_[valve,qp,[(p[b.donor]-p[b.recipient])/b.resistance for b in self.branches]]
        return dict(zip(self.flow_names,q))

    def rhs(self,t,state,open_valves,active=True):
        p,_=self.pressures(t,state,active=active)
        flows=self.flows(t,state,open_valves,active=active,pressures=p)
        dv=self.incidence@np.array(list(flows.values()))
        if abs(dv[-1])>1e-9: raise RuntimeError('Hepatic algebraic flow imbalance')
        dq=np.zeros(4)
        for i,(a,b) in enumerate(EDGES[:4]):
            if open_valves[i]:
                dq[i]=(p[a]-p[b]-self.source.resistances[i]*flows[VALVES[i]])/self.source.inertances[i]
        return np.r_[dv[:-1],dq]
