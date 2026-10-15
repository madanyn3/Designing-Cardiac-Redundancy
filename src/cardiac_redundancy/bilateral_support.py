"""Two parallel lung paths and amended mandatory-organ controls.

PA/PU remain shared vascular reservoirs; no unmeasured regional blood storage
is invented. Separate gas pools, ventilation, nutritive routes and conductances
represent two lungs. Equal native shares are exploratory assumptions.
"""
from dataclasses import dataclass
import math
import numpy as np
from .biological_pump import BiologicalPumpModel
from .respiratory import RespiratoryOxygenTransport
from .oxygen import oxygen_content

@dataclass(frozen=True)
class LungAllocation:
    supported_share: float = .5
    supported_perfusion: float = 1.
    preserved_perfusion: float = .1
    supported_nutrition: float = 1.
    preserved_nutrition: float = .25
    supported_ventilation: float = 1.
    preserved_ventilation: float = 0.
    def __post_init__(self):
        if any(not math.isfinite(v) or not 0<=v<=1 for v in vars(self).values()) or not 0<self.supported_share<1:
            raise ValueError('Finite allocation fractions in range required')
        if self.total('perfusion')<=0 or self.total('nutrition')<=0:
            raise ValueError('At least one perfused/nourished lung required')
    def weights(self,kind):
        return np.array([self.supported_share*getattr(self,'supported_'+kind),
                         (1-self.supported_share)*getattr(self,'preserved_'+kind)])
    def total(self,kind):return float(self.weights(kind).sum())

class BilateralSupport(BiologicalPumpModel):
    """BIO or ideal externally actuated chambers with the same bounded law.

Mechanical mode retains identical cavities/valves/pressure-law parameters but
removes auxiliary myocardial coronary branches. It is a hydraulic actuator
reference, not a claim of realizable hardware or calibrated efficiency.
"""
    def __init__(self,*args,allocation=LungAllocation(),device='BIO',**kwargs):
        super().__init__(*args,**kwargs)
        if device not in ('BIO','mechanical'):raise ValueError('Unknown device')
        self.device=device;self.allocation=allocation
        self.edges=[('lung_supported',a,b) if name=='lungs' else (name,a,b) for name,a,b in self.edges]
        self.edges.extend([('lung_preserved','AO','PU'),('pulmonary_preserved','PA','PU')])
        if device=='mechanical':self.edges=[e for e in self.edges if not e[0].startswith('aux_coronary_')]
        self.flow_names=tuple(e[0] for e in self.edges)
        self.incidence=np.zeros((self.volume_count+1,len(self.edges)));nodes=self.volume_names+('H',)
        for i,(_,a,b) in enumerate(self.edges):self.incidence[nodes.index(a),i]=-1;self.incidence[nodes.index(b),i]=1

    def apply_support_controls(self,time,preservation=.25):
        if not math.isfinite(preservation) or not 0<preservation<=1:raise ValueError('Positive preservation conductance required')
        mandatory={'brain','coronary_left','coronary_right','kidney_protected','liver_arterial','respiratory_muscle'}
        supplies={b.name for b in self.native.base.branches if b.donor=='AO'}
        factors={name:(1. if name in mandatory else preservation) for name in supplies}
        factors['lungs']=self.allocation.total('nutrition')
        factors['pulmonary_exchange']=self.allocation.total('perfusion')
        self.native.set_regional_controls(time=time,factors=factors)

    def flows(self,t,state,open_valves,active=True,pressures=None):
        # The parent selects native flows using self.flow_names. Temporarily
        # recover the renamed nutritive route directly from the native model.
        p=self.pressures(t,state,active)[0] if pressures is None else pressures
        q=super().flows(t,state,open_valves,active,p)
        native=self.native.flows(t,self.native_state(state),open_valves[:4],active,pressures=p)
        w=self.allocation.weights('nutrition')/self.allocation.total('nutrition')
        q['lung_supported']=native['lungs']*w[0];q['lung_preserved']=native['lungs']*w[1]
        w=self.allocation.weights('perfusion')/self.allocation.total('perfusion')
        q['pulmonary_exchange']=native['pulmonary_exchange']*w[0];q['pulmonary_preserved']=native['pulmonary_exchange']*w[1]
        return {name:q[name] for name in self.flow_names}

class BilateralRespiration(RespiratoryOxygenTransport):
    """Keep original total alveolar O2 plus an independent preserved-lung O2.

The otherwise unused oxygenator-gas slot holds preserved-lung gas; supported
lung gas is total minus preserved. This preserves the fuel interface/state
size without double counting. Oxygenators are explicitly unavailable here.
"""
    def __init__(self,*args,allocation=LungAllocation(),**kwargs):
        super().__init__(*args,**kwargs)
        if self.oxygenator_index is not None:raise ValueError('Bilateral fixture has no oxygenator')
        self.allocation=allocation
        self.other_index=[e[0] for e in self.edges].index('pulmonary_preserved')

    def initial_state(self,*args,**kwargs):
        z=super().initial_state(*args,**kwargs)
        z[self.base_size+1]=z[len(self.names)]*(1-self.allocation.supported_share)
        return z

    def evaluate(self,volumes,state,flows,**kwargs):
        state=np.asarray(state);n=len(self.names);b=self.base_size;p=self.parameters
        caps=p.alveolar_volume_mL_STPD*np.array([self.allocation.supported_share,1-self.allocation.supported_share])
        gas=np.array([state[n]-state[b+1],state[b+1]])
        if np.any(gas < -1e-7) or np.any(gas>caps+1e-7):raise ValueError('Regional lung gas outside domain')
        compatible=state.copy();compatible[b+1]=0.
        dz,diag=super().evaluate(volumes,compatible,flows,**kwargs)
        actual_v=diag['ventilation_mL_s_BTPS'];fraction=gas/caps
        q0=flows[self.pulmonary_index];donor0,rec0=('PA','PU') if q0>=0 else ('PU','PA')
        content0=state[self.names.index(donor0)]/volumes[self.names.index(donor0)]
        total_fraction=state[n]/p.alveolar_volume_mL_STPD
        old_transfer=abs(q0)*p.exchange_fraction*(oxygen_content((p.barometric_mmHg-p.water_vapor_mmHg)*total_fraction,p)-content0)
        dz[self.names.index(rec0)]-=old_transfer
        transfers=[]
        for idx,frac in zip((self.pulmonary_index,self.other_index),fraction):
            q=flows[idx];donor,recipient=('PA','PU') if q>=0 else ('PU','PA')
            content=state[self.names.index(donor)]/volumes[self.names.index(donor)]
            x=abs(q)*p.exchange_fraction*(oxygen_content((p.barometric_mmHg-p.water_vapor_mmHg)*frac,p)-content)
            dz[self.names.index(recipient)]+=x;transfers.append(float(x))
        vents=actual_v*self.allocation.weights('ventilation')
        external=vents*p.btps_to_stpd*(p.inspired_fraction-fraction)
        dz[n]=external.sum()-sum(transfers);dz[n+1]=external.sum()
        dz[b+1]=external[1]-transfers[1];dz[b+2]=0.
        diag.update(lung_gas_O2_mL=gas,lung_fractions=fraction,lung_transfer_mL_s=transfers,
                    lung_ventilation_mL_s=vents,lung_external_mL_s=external)
        return dz,diag

    def balance(self,state):
        z=np.asarray(state)
        return super().balance(z)-z[...,self.base_size+1]
