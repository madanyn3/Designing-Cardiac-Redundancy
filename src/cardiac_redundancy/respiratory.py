"""M4 finite gas sources, assumed ventilatory feedback and tissue-store sensitivity.

These are conservative exploratory closures, not a respiratory mechanics,
CO2, injury or recovery model. See docs/m4-respiratory-verification.md.
"""
from dataclasses import dataclass
import math
import numpy as np
from .organ_oxygen import OrganOxygenTransport
from .oxygen import oxygen_content


@dataclass(frozen=True)
class RespiratoryParameters:
    reference_ventilation_mL_s_BTPS: float = 4690/60
    response_time_s: float = 5.
    oxygenator_gas_volume_mL_STPD: float = 100.
    oxygenator_initial_fraction: float = 1.
    oxygenator_inspired_fraction: float = 1.
    oxygenator_pressure_mmHg: float = 713.
    oxygenator_exchange_fraction: float = 1.
    oxygenator_transfer_limit_mL_s: float = 5.
    store_fill_time_s: float = 1.
    store_release_time_s: float = 1.

    def __post_init__(self):
        if any(not math.isfinite(x) for x in vars(self).values()): raise ValueError('Finite respiratory inputs required')
        if min(self.response_time_s,self.oxygenator_gas_volume_mL_STPD,self.oxygenator_pressure_mmHg,self.store_fill_time_s,self.store_release_time_s)<=0:
            raise ValueError('Positive volumes, pressures and response times required')
        if min(self.reference_ventilation_mL_s_BTPS,self.oxygenator_transfer_limit_mL_s)<0:raise ValueError('Nonnegative capacities required')
        if any(not 0<=x<=1 for x in (self.oxygenator_initial_fraction,self.oxygenator_inspired_fraction,self.oxygenator_exchange_fraction)):
            raise ValueError('Fractions must be in [0,1]')


class RespiratoryOxygenTransport(OrganOxygenTransport):
    """Adds ventilation, oxygenator O2/external ledger, and optional tissue stores.

    Tissue states hold USABLE oxygen equivalents, not total tissue oxygen.
    Capacity = baseline demand * assumed storage seconds / 60. Zero is M1's
    no-credit reference. Nonzero cases are sensitivity assumptions.
    """
    def __init__(self,*args,respiratory=RespiratoryParameters(),store_seconds=0.,**kwargs):
        super().__init__(*args,**kwargs)
        self.respiratory=respiratory
        self.base_size=len(self.names)+2+2*self.organ_count
        self.reference_demand=np.array([o.demand_mL_min for o in self.organs])
        seconds=np.broadcast_to(np.asarray(store_seconds,dtype=float),self.reference_demand.shape)
        if np.any(~np.isfinite(seconds)) or np.any(seconds<0):raise ValueError('Invalid store credit')
        self.store_capacity=self.reference_demand*seconds/60
        self.resp_indices=[j for j,o in enumerate(self.organs) if o.parent_pool=='respiratory_muscle']
        self.oxygenator_index=next((j for j,e in enumerate(self.edges) if e[0]=='oxygenator'),None)

    def initial_state(self,volumes,po2_mmHg=40.,alveolar_fraction=.14,store_fraction=1.):
        if not math.isfinite(store_fraction) or not 0<=store_fraction<=1:raise ValueError('Invalid initial store fraction')
        base=super().initial_state(volumes,po2_mmHg,alveolar_fraction)
        p=self.respiratory
        gas=p.oxygenator_gas_volume_mL_STPD*p.oxygenator_initial_fraction if self.oxygenator_index is not None else 0.
        return np.r_[base,p.reference_ventilation_mL_s_BTPS,gas,0.,self.store_capacity*store_fraction]

    def evaluate(self,volumes,state,flows,*,mode='maintained',ventilation_fraction=1.,drive_fraction=1.,
                 airway_fraction=1.,sweep_mL_s_STPD=0.,demand_mL_min=None):
        """Modes: maintained/reduced prescribed boundary, apnea, or assumed feedback.

        Feedback ventilation target scales with ACTUAL respiratory O2 uptake
        divided by the reference respiratory requirement, not by reduced demand.
        Drive/airway fractions are independent declared capabilities.
        """
        p=self.respiratory;b=self.base_size;state=np.asarray(state,dtype=float);flows=np.asarray(flows,dtype=float)
        if state.shape!=(b+3+self.organ_count,) or np.any(~np.isfinite(state)):raise ValueError('Invalid respiratory state')
        if mode not in ('maintained','reduced','apnea','feedback'):raise ValueError('Unknown respiratory mode')
        if any(not math.isfinite(x) or not 0<=x<=1 for x in (ventilation_fraction,drive_fraction,airway_fraction)) or not math.isfinite(sweep_mL_s_STPD) or sweep_mL_s_STPD<0:
            raise ValueError('Invalid respiratory boundary')
        v,gas,_=state[b:b+3];stores=state[b+3:]
        # Implicit Jacobian probes and explicit Runge-Kutta internal stages
        # can cross a boundary before error control rejects the step. Permit
        # a trial extension, but never clip state or inventories;
        # accepted trajectories are checked separately at 1e-7 mL tolerance.
        store_trial=1e-3*np.maximum(1.,self.store_capacity)
        if v < -1e-7 or np.any(stores < -store_trial) or np.any(stores>self.store_capacity+store_trial) or gas < -1e-7 or gas>p.oxygenator_gas_volume_mL_STPD*(1+1e-3):
            raise ValueError(f'Respiratory store outside model domain: ventilation={v}, gas={gas}, store_min={stores.min()}, store_excess={np.max(stores-self.store_capacity)}')
        d=np.array(self.reference_demand if demand_mL_min is None else demand_mL_min,dtype=float)
        if d.shape!=self.reference_demand.shape or np.any(~np.isfinite(d)) or np.any(d<0):raise ValueError('Invalid metabolic demand')
        actual_v=0. if mode=='apnea' else max(0.,v) if mode=='feedback' else p.reference_ventilation_mL_s_BTPS*ventilation_fraction
        # Charge stores only with oxygen actually extracted from passing blood.
        requested=d+60*np.maximum(self.store_capacity-stores,0)/p.store_fill_time_s
        base_rate,diag=super().evaluate(volumes,state[:b],flows,actual_v,requested)
        exchange=diag['consumed'].copy()
        uptake=np.minimum(d/60,exchange+np.maximum(stores,0)/p.store_release_time_s)
        store_rate=exchange-uptake
        n=len(self.names);m=self.organ_count
        base_rate[n+2:n+2+m]=uptake;base_rate[n+2+m:]=d/60-uptake
        resp_ref=self.reference_demand[self.resp_indices].sum()/60
        if mode=='feedback' and resp_ref<=0:raise ValueError('Feedback requires explicit respiratory-muscle reference demand')
        available=min(1.,float(uptake[self.resp_indices].sum()/resp_ref)) if resp_ref>0 else 0.
        target=p.reference_ventilation_mL_s_BTPS*ventilation_fraction*drive_fraction*airway_fraction*available if mode=='feedback' else actual_v
        dv=(target-v)/p.response_time_s
        transfer=external=0.
        if self.oxygenator_index is not None:
            q=flows[self.oxygenator_index];_,a,c=self.edges[self.oxygenator_index]
            donor,recipient=(a,c) if q>=0 else (c,a)
            incoming=state[self.names.index(donor)]/volumes[self.names.index(donor)]
            fraction=max(0.,gas)/p.oxygenator_gas_volume_mL_STPD
            raw=abs(q)*p.oxygenator_exchange_fraction*(oxygen_content(p.oxygenator_pressure_mmHg*fraction,self.parameters)-incoming)
            transfer=float(np.clip(raw,-p.oxygenator_transfer_limit_mL_s,p.oxygenator_transfer_limit_mL_s))
            base_rate[self.names.index(recipient)]+=transfer
            external=sweep_mL_s_STPD*(p.oxygenator_inspired_fraction-fraction)
        elif sweep_mL_s_STPD:raise ValueError('Sweep gas supplied without oxygenator')
        diag.update(demand=d/60,consumed=uptake,unmet=d/60-uptake,blood_extracted=exchange,
            store_rate=store_rate,ventilation_mL_s_BTPS=actual_v,ventilation_target_mL_s_BTPS=target,
            respiratory_availability=available,oxygenator_transfer_mL_s=transfer,oxygenator_external_mL_s=external)
        return np.r_[base_rate,dv,external-transfer,external,store_rate],diag

    def rhs(self,volumes,state,flows,**kwargs):return self.evaluate(volumes,state,flows,**kwargs)[0]

    def balance(self,state):
        state=np.asarray(state);b=self.base_size
        return super().balance(state[...,:b])+state[...,b+1]-state[...,b+2]+state[...,b+3:].sum(axis=-1)
