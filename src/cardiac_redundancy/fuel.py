"""Finite liver glycogen and transported glucose, with optional O2 coupling.

Amounts are mmol glucose equivalents, blood volumes mL and time seconds.
This restrictive model supplies no gluconeogenesis, fat, ketones or exogenous
glucose. A fractional glucose share leaves other substrates UNASSESSED.
See docs/m4-fuel-accounting.md for evidence and the allocation assumption.
"""
from dataclasses import dataclass
import math
import numpy as np


GLUCOSE_MG_PER_MMOL = 180.156
GLUCOSYL_MG_PER_MMOL = 162.141
O2_ML_STPD_PER_MMOL_GLUCOSE = 6 * 22.414


@dataclass(frozen=True)
class FuelParameters:
    liver_volume_L: float = 1.44
    glycogen_mmol_L_liver: float = 229.
    accessible_fraction: float = 1.
    initial_plasma_glucose_mmol_L: float = 5.
    plasma_to_whole_blood_ratio: float = 1.11
    release_cap_mmol_min: float = .21
    release_time_s: float = 1.
    export_increment_mmol_L: float = 1.
    blood_target_time_s: float = 60.
    nonoxidative_mmol_min: float = 0.
    nonoxidative_release_time_s: float = 1.

    def __post_init__(self):
        if any(not math.isfinite(v) for v in vars(self).values()):
            raise ValueError('Finite fuel parameters required')
        if min(self.liver_volume_L, self.plasma_to_whole_blood_ratio,
               self.release_time_s, self.blood_target_time_s,
               self.nonoxidative_release_time_s) <= 0:
            raise ValueError('Positive fuel volumes, conversion and times required')
        if min(self.glycogen_mmol_L_liver, self.initial_plasma_glucose_mmol_L,
               self.release_cap_mmol_min, self.export_increment_mmol_L,
               self.nonoxidative_mmol_min) < 0 or not 0 <= self.accessible_fraction <= 1:
            raise ValueError('Invalid fuel stock or capacity')

    @property
    def glycogen_mmol(self):
        return self.liver_volume_L*self.glycogen_mmol_L_liver*self.accessible_fraction

    @property
    def blood_concentration_mmol_mL(self):
        return self.initial_plasma_glucose_mmol_L/self.plasma_to_whole_blood_ratio/1000


class GlucoseTransport:
    """Conservative signed advection and serial uptake before/after hepatic mixing.

    State: blood[n], accessible glycogen, cumulative organ uptake[m], cumulative
    unmet request[m], cumulative release, nonoxidative uptake, nonoxidative unmet.
    The release ledger is diagnostic, NOT an additional physical inventory.
    Release enters the recipient of forward hepatic outflow; zero/reverse
    hepatic outflow blocks export in this deliberately restrictive closure.
    """
    def __init__(self, names, edges, organs, parameters=FuelParameters(),
                 extraction_limit=.5, hepatic_edge='hepatic_outflow'):
        self.names=tuple(names); self.edges=tuple(edges); self.organs=tuple(organs)
        self.parameters=parameters; self.n=len(names); self.m=len(organs)
        self.by_edge={o.edge:j for j,o in enumerate(organs)}
        if len(set(names))!=self.n or 'H' in names or len(self.by_edge)!=self.m:
            raise ValueError('Unique physical volumes and organ edges required')
        if len({o.name for o in organs})!=self.m or len({e[0] for e in edges})!=len(edges):
            raise ValueError('Unique organ and flow names required')
        if any(a not in (*names,'H') or b not in (*names,'H') or a==b for _,a,b in edges):
            raise ValueError('Unknown or self-connected transport node')
        if not set(self.by_edge)<=set(e[0] for e in edges):raise ValueError('Unknown uptake edge')
        self.hepatic_index=next((j for j,e in enumerate(edges) if e[0]==hepatic_edge),None)
        if self.hepatic_index is None or edges[self.hepatic_index][1]!='H' or edges[self.hepatic_index][2]=='H':
            raise ValueError('Forward hepatic outflow H to physical recipient required')
        self.extraction=np.broadcast_to(np.asarray(extraction_limit,dtype=float),(self.m,)).copy()
        if np.any(~np.isfinite(self.extraction)) or np.any((self.extraction<0)|(self.extraction>1)):
            raise ValueError('Glucose extraction fractions must be in [0,1]')
        self.size=self.n+1+2*self.m+3

    def initial_state(self, volumes):
        v=np.asarray(volumes,dtype=float)
        if v.shape!=(self.n,) or np.any(~np.isfinite(v)) or np.any(v<=0):
            raise ValueError('Positive physical volumes required')
        return np.r_[v*self.parameters.blood_concentration_mmol_mL,
                     self.parameters.glycogen_mmol,np.zeros(2*self.m+3)]

    def evaluate(self,volumes,state,flows,demand_mmol_s,*,liver_function_fraction=1.,actual_mmol_s=None,numerical_trial=False):
        """Reserve against requested demand, or withdraw an explicitly smaller use.

        `actual_mmol_s` supports conservative oxygen coupling: unused reserved
        substrate returns to blood, rather than disappearing or being burned
        without oxygen. The physiological release-cap transfer and the liver
        function fraction are independent, declared assumptions.
        """
        n=self.n;m=self.m;p=self.parameters
        v=np.asarray(volumes,dtype=float);y=np.asarray(state,dtype=float);q=np.asarray(flows,dtype=float)
        d=np.asarray(demand_mmol_s,dtype=float)
        if v.shape!=(n,) or np.any(~np.isfinite(v)) or np.any(v<=0):raise ValueError('Invalid volumes')
        # Implicit solver iterates can cross zero before step rejection. The
        # optional rate extension uses nonnegative available inventories, never
        # clips the state, and requires separate accepted-state domain checks.
        if y.shape!=(self.size,) or np.any(~np.isfinite(y)) or (not numerical_trial and np.min(y[:n+1]) < -1e-7):
            raise ValueError('Invalid glucose/glycogen state')
        if q.shape!=(len(self.edges),) or np.any(~np.isfinite(q)):raise ValueError('Invalid flows')
        if d.shape!=(m,) or np.any(~np.isfinite(d)) or np.any(d<0):raise ValueError('Invalid glucose demand')
        if not math.isfinite(liver_function_fraction) or not 0<=liver_function_fraction<=1:
            raise ValueError('Invalid liver function fraction')
        use=None if actual_mmol_s is None else np.asarray(actual_mmol_s,dtype=float)
        if use is not None and (use.shape!=(m,) or np.any(~np.isfinite(use)) or np.any(use<0) or np.any(use>d+1e-12)):
            raise ValueError('Invalid actual glucose use')
        c=dict(zip(self.names,np.maximum(y[:n],0)/v))
        change=dict.fromkeys(self.names,0.); delivered=np.zeros(m); uptake=np.zeros(m)
        oriented=[(name,a if f>=0 else b,b if f>=0 else a,abs(f)) for (name,a,b),f in zip(self.edges,q)]
        incoming=sum(f for _,_,b,f in oriented if b=='H')
        outgoing=sum(f for _,a,_,f in oriented if a=='H')
        if abs(incoming-outgoing)>1e-7:raise ValueError('Unbalanced hepatic junction')
        incoming_glucose=0.

        def move(edge):
            nonlocal incoming_glucose
            name,a,b,f=edge;flux=f*c[a];u=0.
            if name in self.by_edge:
                j=self.by_edge[name];delivered[j]=flux;capacity=self.extraction[j]*flux
                u=min(d[j],capacity) if use is None else use[j]
                if u>capacity+1e-10:raise ValueError('Actual glucose use exceeds delivery capacity')
                uptake[j]=u
            if a!='H':change[a]-=flux
            if b=='H':incoming_glucose+=flux-u
            else:change[b]+=flux-u

        for edge in oriented:
            if edge[1]!='H':move(edge)
        c['H']=incoming_glucose/incoming if incoming else 0.
        for edge in oriented:
            if edge[1]=='H':move(edge)
        # This controller is bookkeeping, not a hormone/biochemical model.
        # It never credits unmodeled gluconeogenesis or uses reverse hepatic flow.
        blood_gap=max(0.,p.blood_concentration_mmol_mL*v.sum()-y[:n].sum())
        requested_release=d.sum()+p.nonoxidative_mmol_min/60+blood_gap/p.blood_target_time_s
        release=min(requested_release,p.release_cap_mmol_min/60*liver_function_fraction,
                    max(0.,y[n])/p.release_time_s,
                    max(0.,q[self.hepatic_index])*p.export_increment_mmol_L/1000)
        recipient=self.edges[self.hepatic_index][2];change[recipient]+=release
        nonox_demand=p.nonoxidative_mmol_min/60*v/v.sum()
        nonox=np.minimum(nonox_demand,np.maximum(y[:n],0)/p.nonoxidative_release_time_s)
        blood_rate=np.array([change[name] for name in self.names])-nonox
        rate=np.r_[blood_rate,-release,uptake,d-uptake,release,nonox.sum(),(nonox_demand-nonox).sum()]
        return rate,dict(delivered_mmol_s=delivered,consumed_mmol_s=uptake,demand_mmol_s=d,
                         unmet_mmol_s=d-uptake,release_mmol_s=release,hepatic_mmol_mL=c['H'])

    def balance(self,state):
        y=np.asarray(state);n=self.n;m=self.m
        return y[...,:n+1].sum(axis=-1)+y[...,n+1:n+1+m].sum(axis=-1)+y[...,-2]


class FuelOxygenTransport:
    """Couple the existing respiratory O2 model to a conserved glucose supply.

    Glucose reservations follow the full requested demand in vascular order.
    They cap O2 demand BEFORE actual oxidation and ventilatory feedback. Actual
    oxidation then determines actual glucose withdrawal. Unused reservations
    remain in circulation; no instantaneous downstream reallocation is assumed.
    This conservative allocation can underuse supply when upstream O2 is low.
    Original demand, including any unmet amount, remains in both ledgers.
    """
    def __init__(self,oxygen,glucose,glucose_share=1.):
        if tuple(oxygen.names)!=glucose.names or tuple(oxygen.edges)!=glucose.edges or tuple(oxygen.organs)!=glucose.organs:
            raise ValueError('Oxygen and glucose must use identical ordered networks and demands')
        self.oxygen=oxygen;self.glucose=glucose
        self.oxygen_size=oxygen.base_size+3+oxygen.organ_count
        self.share=np.broadcast_to(np.asarray(glucose_share,dtype=float),(glucose.m,)).copy()
        if np.any(~np.isfinite(self.share)) or np.any((self.share<0)|(self.share>1)):
            raise ValueError('Glucose shares must be in [0,1]')

    def initial_state(self,volumes,**kwargs):
        return np.r_[self.oxygen.initial_state(volumes,**kwargs),self.glucose.initial_state(volumes)]

    def evaluate(self,volumes,state,flows,*,demand_mL_min=None,liver_function_fraction=1.,**respiratory):
        k=self.oxygen_size;g=self.glucose;o=self.oxygen
        state=np.asarray(state,dtype=float)
        if state.shape!=(k+g.size,):raise ValueError('Invalid coupled state')
        d=np.array(o.reference_demand if demand_mL_min is None else demand_mL_min,dtype=float)
        if d.shape!=(g.m,) or np.any(~np.isfinite(d)) or np.any(d<0):raise ValueError('Invalid oxygen demand')
        gd=d/60*self.share/O2_ML_STPD_PER_MMOL_GLUCOSE
        _,reservation=g.evaluate(volumes,state[k:],flows,gd,liver_function_fraction=liver_function_fraction)
        ceiling=d.copy();selected=self.share>0
        ceiling[selected]=np.minimum(d[selected],reservation['consumed_mmol_s'][selected]
                                     *O2_ML_STPD_PER_MMOL_GLUCOSE*60/self.share[selected])
        ro,diag=o.evaluate(volumes,state[:k],flows,demand_mL_min=ceiling,**respiratory)
        actual=diag['consumed']*self.share/O2_ML_STPD_PER_MMOL_GLUCOSE
        rg,gdg=g.evaluate(volumes,state[k:],flows,gd,liver_function_fraction=liver_function_fraction,actual_mmol_s=actual)
        n=g.n;m=g.m
        # O2's requested-demand ledger must retain the user's original target.
        ro[n+2+m:n+2+2*m]=d/60-diag['consumed']
        diag.update(demand=d/60,unmet=d/60-diag['consumed'],
                    glucose_ceiling_mL_s=ceiling/60,glucose=gdg,
                    other_substrates_unassessed=bool(np.any(self.share<1)))
        return np.r_[ro,rg],diag

    def rhs(self,volumes,state,flows,**kwargs):return self.evaluate(volumes,state,flows,**kwargs)[0]

    def balances(self,state):
        y=np.asarray(state);k=self.oxygen_size
        return self.oxygen.balance(y[...,:k]),self.glucose.balance(y[...,k:])
