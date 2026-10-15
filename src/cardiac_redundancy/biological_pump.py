"""Bounded exploratory controls for the Smith-shaped auxiliary fixture.

Bounds are search assumptions, not measured biological limits. Scaling wall
coordinates does not establish myocardial mass or anatomical dimensions.
"""
from dataclasses import dataclass, fields
import math
import numpy as np
from .circuit import CircuitModel
from .mechanics import Wall, septal_root
from .pressure_regional import waterfall

@dataclass(frozen=True)
class BiologicalParameters:
    left_scale: float = .25
    right_scale: float = .25
    septal_scale: float = .25
    pericardial_scale: float = .25
    period_s: float = .75
    left_contractility: float = 1.
    right_contractility: float = 1.
    septal_contractility: float = 1.
    def __post_init__(self):
        for f in fields(self):
            x=getattr(self,f.name)
            lo,hi=(.1,1.) if 'scale' in f.name else ((.4,1.5) if f.name=='period_s' else (.5,2.))
            if isinstance(x,bool) or not isinstance(x,(int,float)) or not math.isfinite(x) or not lo<=x<=hi:
                raise ValueError(f'{f.name} must be within [{lo}, {hi}]')

def scaled_wall(w,scale,contractility):
    return Wall(w.elastance*contractility/scale,w.dead_volume*scale,
                w.pressure_scale,w.stiffness/scale,w.reference_volume*scale)

class BiologicalPumpModel(CircuitModel):
    def __init__(self,*args,biological=BiologicalParameters(),**kwargs):
        super().__init__(*args,architecture='BIO',**kwargs)
        self.biological=biological
        self.auxiliary_period_s=biological.period_s
        scales=np.array([biological.left_scale,biological.right_scale])
        cavities=scales*self.native.base.source.initial_state[:2]
        # Replace the legacy allocation with a conservative donor transfer.
        for j,donor in enumerate(('PU','VC')):
            self.initial_state[self.volume_names.index(donor)]+=self.initial_state[16+j]-cavities[j]
        self.initial_state[16:18]=cavities
        if min(self.domain_margins(self.initial_state).values())<=0:
            raise ValueError('Auxiliary allocation exceeds donor-domain inventory')
        self.valve_R[4:]=self.native.base.source.resistances/np.repeat(scales,2)
        self.valve_L[4:]=self.native.base.source.inertances/np.repeat(scales,2)
        self.auxiliary_walls=tuple(scaled_wall(w,s,c) for w,s,c in zip(self.native.base.source.walls,
            (*scales,biological.septal_scale),(biological.left_contractility,biological.right_contractility,biological.septal_contractility)))

    def pressures(self,t,state,active=True):
        p,meta=self.native.pressures(t,self.native_state(state),active=active)
        lv,rv=state[16:18]
        if min(lv,rv)<=0:raise ValueError('Nonpositive auxiliary cavity')
        elapsed=0. if self.trigger is None else t-self.trigger
        ramp=0. if self.trigger is None else min(1.,max(0.,elapsed/self.auxiliary_ramp_s))
        e=ramp*math.exp(-45*((elapsed%self.auxiliary_period_s)/self.auxiliary_period_s-.5)**2)
        if self.withdrawal is not None:e*=1-min(1.,max(0.,t-self.withdrawal))
        septum,residual=septal_root(lv,rv,e,*self.auxiliary_walls)
        source=self.native.base.source;pars=source.parameters['pericardium'];s=self.biological.pericardial_scale
        peri=source.config['pressure_external_mmHg']['pericardium']+pars['P_0_pcd']*math.expm1(pars['lambda_pcd']*((lv+rv)/s-pars['V_0_pcd']))
        p.update(aLV=self.auxiliary_walls[0].pressure(lv-septum,e)+peri,
                 aRV=self.auxiliary_walls[1].pressure(rv+septum,e)+peri,aPERI=peri)
        meta.update(aux_activation=e,aux_septal_residual_mmHg=residual,aux_septal_coordinate_mL=septum)
        return p,meta

    def flows(self,t,state,open_valves,active=True,pressures=None):
        p=self.pressures(t,state,active)[0] if pressures is None else pressures
        q=super().flows(t,state,open_valves,active,p)
        for side,node,scale in (('left','aLV',self.biological.left_scale),('right','aRV',self.biological.right_scale)):
            pc=max(p['aPERI'],0)+self.native.base.closure_settings['kappa']*max(p[node]-p['aPERI'],0)
            q['aux_coronary_'+side]=waterfall(p['AO'],p['VC'],pc,self.native.base.branch_by_name['coronary_'+side].resistance/scale)
        return q
