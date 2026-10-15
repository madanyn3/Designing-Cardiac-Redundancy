"""Independent regional conductance controls on the verified selective fixture.

Roles describe the objective; they never prescribe flow or imply recoverability.
The hepatic zero-storage pressure is solved using the controlled conductances
before any flow is evaluated. Native valves and support ports are not organ gates.
"""
import math
from .selective import SelectiveModel


class WholeNetworkModel(SelectiveModel):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self.regional_factors={n:1. for n in self.flow_names if n not in
            ('MV','AV','TV','PV','support_right','support_left')}
        # Valve names are obtained from the parent as well, avoiding aliases.
        for n in self.base.flow_names[:4]:self.regional_factors.pop(n,None)
        self.regional_ramp=None

    def regional_at(self,time):
        if self.regional_ramp is None:return dict(self.regional_factors)
        start,duration,old,new=self.regional_ramp
        f=min(1.,max(0.,(time-start)/duration))
        return {n:old[n]+f*(new[n]-old[n]) for n in new}

    def set_regional_controls(self,*,time,factors,ramp_s=0.):
        """Call at a segment boundary; factors multiply conductance, not flow targets."""
        if not math.isfinite(time) or not math.isfinite(ramp_s) or ramp_s<0:raise ValueError('Invalid regional event')
        if set(factors)-set(self.regional_factors):raise ValueError('Unknown or non-regional edge')
        if any(not math.isfinite(f) or not 0<=f<=1 for f in factors.values()):raise ValueError('Regional factors must be 0..1')
        old=self.regional_at(time);new={**old,**factors}
        # A disconnected zero-storage node has no defined pressure. Reject it
        # explicitly rather than insert an arbitrary pressure or delete blood.
        if not any(new[b.name]>0 for b in self.base.hepatic_branches):raise ValueError('Isolated hepatic junction')
        self.regional_factors=new
        self.regional_ramp=(time,ramp_s,old,new) if ramp_s else None

    def pressures(self,t,state,active=True):
        p,meta=super().pressures(t,state,active=active)
        factors=self.regional_at(t)
        branches=self.base.hepatic_branches
        total=sum(factors[b.name]/b.resistance for b in branches)
        if total<=0:raise ValueError('Isolated hepatic junction')
        p['H']=sum(p[b.recipient if b.donor=='H' else b.donor]*factors[b.name]/b.resistance for b in branches)/total
        return p,meta

    def flows(self,t,state,open_valves,active=True,pressures=None):
        p=self.pressures(t,state,active=active)[0] if pressures is None else pressures
        q=super().flows(t,state,open_valves,active=active,pressures=p)
        factors=self.regional_at(t)
        for n,f in factors.items():q[n]*=f
        return q


def validate_roles(model,scenario):
    """Every tissue supply edge belongs to exactly one role; no silent omissions."""
    supplies={b.name for b in model.base.branches if b.donor=='AO'}
    active=set(scenario['active_supply_edges']);preservation=set(scenario['preservation_supply_edges'])
    if active & preservation or active|preservation!=supplies:raise ValueError('Supply roles must partition all tissue supply edges')
    if active!={'brain','lungs','respiratory_muscle','liver_arterial'}:raise ValueError('Active set disagrees with repair-window objective')
    if not {'coronary_left','coronary_right'}<=preservation:raise ValueError('Native heart preservation must remain explicit')
    if set(scenario['dependency_edges'])!={'pulmonary_exchange','hepatic_outflow'}:raise ValueError('Missing pulmonary/hepatic dependency')
    return True
