"""M3 circuit fixtures: B1 dual ports, E1/E2 ideal oxygenator, biological heart.

Oxygenator labels specify topology only; no oxygen source or gas law is present.
The auxiliary similarity fixture has two cavities, its own septum/pericardium,
four inertial valves and coronary loads. Its cavity allocation occurs once.
"""
import math
import numpy as np
from .backup_pump import BackupPumpModel, PumpStage
from .selective import finite_head_flow
from .mechanics import septal_root
from .pressure_regional import waterfall
from .source import EDGES, VALVES


class CircuitModel:
    def __init__(self, *args, architecture='B1', **kwargs):
        if architecture not in ('B1','E1','E2','BIO'):
            raise ValueError('Unknown circuit architecture')
        self.native=BackupPumpModel(*args,**kwargs)
        self.architecture=architecture
        self.settings=self.native.settings
        self.config=self.native.config
        self.volume_names=self.native.volume_names
        self.valve_names=VALVES
        self.valve_edges=EDGES[:4]
        self.alpha=self.native.base.overlay['auxiliary_fixture']['alpha']
        self.trigger=None
        self.withdrawal=None
        self.auxiliary_ramp_s=1.
        self.auxiliary_period_s=.75
        self.oxygenator_stage=PumpStage()
        self.edges=[e for e in self.native.edges if not e[0].startswith('support_')]
        y=self.native.initial_state.copy()
        if architecture=='BIO':
            self.volume_names+=('aLV','aRV')
            self.valve_names+=tuple('a'+v for v in VALVES)
            self.valve_edges+=(('PU','aLV'),('aLV','AO'),('VC','aRV'),('aRV','PA'))
            cavities=self.alpha*self.native.base.source.initial_state[:2]
            # One-time physical transfer before warm-up, not at activation.
            y[self.volume_names.index('PU')]-=cavities[0]
            y[self.volume_names.index('VC')]-=cavities[1]
            y=np.r_[y[:16],cavities,y[16:],np.zeros(4)]
            self.edges.extend((name,a,b) for name,(a,b) in zip(self.valve_names[4:],self.valve_edges[4:]))
            self.edges.extend([('aux_coronary_left','AO','VC'),('aux_coronary_right','AO','VC')])
        elif architecture=='B1':
            self.edges.extend([('support_right','VC','PA'),('support_left','PU','AO')])
        else:
            self.edges.append(('oxygenator','VC','AO'))
            if architecture=='E2':self.edges.append(('support_left','PU','AO'))
        self.volume_count=len(self.volume_names)
        self.flow_names=tuple(e[0] for e in self.edges)
        self.valve_enabled=[True]*4+([False]*4 if architecture=='BIO' else [])
        self.initial_state=y
        self.incidence=np.zeros((self.volume_count+1,len(self.edges)))
        nodes=self.volume_names+('H',)
        for i,(_,a,b) in enumerate(self.edges):
            self.incidence[nodes.index(a),i]=-1
            self.incidence[nodes.index(b),i]=1
        self.valve_R=np.tile(self.native.base.source.resistances,2 if architecture=='BIO' else 1)
        self.valve_L=np.tile(self.native.base.source.inertances,2 if architecture=='BIO' else 1)
        if architecture=='BIO':
            self.valve_R[4:]/=self.alpha;self.valve_L[4:]/=self.alpha

    def native_state(self,state):
        return np.r_[state[:16],state[self.volume_count:self.volume_count+4]]

    def release_auxiliary(self,time):
        if self.architecture!='BIO' or self.trigger is not None or not math.isfinite(time):
            raise ValueError('Auxiliary release requires an isolated biological fixture')
        self.trigger=time
        self.valve_enabled[4:]=[True]*4

    def domain_margins(self,state):
        margins=self.native.domain_margins(self.native_state(state))
        if self.architecture=='BIO':margins.update(aLV=state[16],aRV=state[17])
        return margins

    def pressures(self,t,state,active=True):
        p,meta=self.native.pressures(t,self.native_state(state),active=active)
        if self.architecture=='BIO':
            # Similarity scaling is evaluated in source coordinates V/alpha:
            # E and lambda /alpha; Vd,V0 *alpha gives identical pressures.
            lv,rv=np.asarray(state[16:18])/self.alpha
            if min(lv,rv)<=0:raise ValueError('Nonpositive auxiliary cavity')
            elapsed=0. if self.trigger is None else t-self.trigger
            ramp=0. if self.trigger is None else min(1.,max(0.,elapsed/self.auxiliary_ramp_s))
            e=ramp*math.exp(-45*((elapsed%self.auxiliary_period_s)/self.auxiliary_period_s-.5)**2)
            if self.withdrawal is not None:
                e*=1-min(1.,max(0.,t-self.withdrawal))
            source=self.native.base.source
            septum,residual=septal_root(lv,rv,e,*source.walls)
            pars=source.parameters['pericardium']
            peri=source.config['pressure_external_mmHg']['pericardium']+pars['P_0_pcd']*math.expm1(pars['lambda_pcd']*(lv+rv-pars['V_0_pcd']))
            p.update(aLV=source.walls[0].pressure(lv-septum,e)+peri,
                     aRV=source.walls[1].pressure(rv+septum,e)+peri,aPERI=peri)
            meta.update(aux_activation=e,aux_septal_residual_mmHg=residual,
                        aux_septal_coordinate_mL=septum*self.alpha)
        return p,meta

    def flows(self,t,state,open_valves,active=True,pressures=None):
        p=self.pressures(t,state,active)[0] if pressures is None else pressures
        base=self.native.flows(t,self.native_state(state),open_valves[:4],active,pressures=p)
        q={n:base[n] for n in self.flow_names if n in base}
        if self.architecture=='BIO':
            for i,n in enumerate(self.valve_names[4:],4):
                q[n]=state[self.volume_count+i] if self.valve_enabled[i] and open_valves[i] else 0.
            for side,node in (('left','aLV'),('right','aRV')):
                pc=max(p['aPERI'],0)+self.native.base.closure_settings['kappa']*max(p[node]-p['aPERI'],0)
                resistance=self.native.base.branch_by_name['coronary_'+side].resistance/self.alpha
                q['aux_coronary_'+side]=waterfall(p['AO'],p['VC'],pc,resistance)
        elif self.architecture in ('E1','E2'):
            c=self.native.control_at(t)
            # Use left-head ramp for the VC-to-AO fixture, with separately
            # configurable resistance/capacity; no extra oxygenator storage.
            q['oxygenator']=finite_head_flow(p['VC'],p['AO'],c['left_head'],
                self.oxygenator_stage.resistance_mmHg_s_mL,self.oxygenator_stage.capacity_mL_s) if c['pump_on'] else 0.
        return q

    def rhs(self,t,state,open_valves,active=True):
        p,_=self.pressures(t,state,active)
        q=self.flows(t,state,open_valves,active,p)
        dv=self.incidence@np.array([q[n] for n in self.flow_names])
        if abs(dv[-1])>1e-9:raise RuntimeError('Hepatic conservation failure')
        dq=np.zeros(len(self.valve_names))
        for i,(a,b) in enumerate(self.valve_edges):
            if self.valve_enabled[i] and open_valves[i]:
                dq[i]=(p[a]-p[b]-self.valve_R[i]*q[self.valve_names[i]])/self.valve_L[i]
        return np.r_[dv[:-1],dq]
