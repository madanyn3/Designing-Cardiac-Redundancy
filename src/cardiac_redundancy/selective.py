"""Exploratory arterial/venous storage, reversible routing and B1 support.

The five tissue territories remain pooled across limbs; they are not four
anatomical limbs. Storage shares and resistance shares are explicit numerical
fixtures, not measured anatomy. All control changes preserve continuous blood
states. No clotting, cell injury, oxygen uptake or clinical recovery is modeled.
"""
import math
import numpy as np
from .pressure_regional import PressureRegionalModel
from .source import EDGES, VALVES


def finite_head_flow(pin, pout, head, resistance, capacity):
    """Check-valved, finite-head/finite-capacity abstract pump [mL/s]."""
    if not all(math.isfinite(v) for v in (pin,pout,head,resistance,capacity)) or resistance<=0 or min(head,capacity)<0:
        raise ValueError('Invalid pump parameter')
    return min(capacity,max(0.,(pin+head-pout)/resistance))


class SelectiveModel:
    """16 physical volume states followed by the four native valve flows.

    Each old limb pool is partitioned into a proximal arterial/microvascular
    storage and a distal venous storage. Inflow gate precedes proximal storage;
    the intervening resistor remains open when both external gates close.
    Trapping therefore conserves combined regional volume, not each subvolume.
    """
    def __init__(self, source, project, overlay, *, arterial_fraction, inlet_resistance_fraction):
        if not 0<arterial_fraction<1 or not 0<inlet_resistance_fraction<1:
            raise ValueError('Storage and resistance fractions must lie strictly between zero and one')
        self.base=PressureRegionalModel(source,project,overlay)
        self.settings=self.base.settings; self.config=self.base.config
        self.fraction=arterial_fraction; self.inlet_fraction=inlet_resistance_fraction
        self.beds=tuple(b.name for b in self.base.branches if b.name.startswith('limb_') and not b.name.endswith('_return'))
        self.pool_names=tuple(self.base.branch_by_name[n].recipient for n in self.beds)
        self.arterial_names=tuple('A_'+n for n in self.beds)
        self.volume_names=self.base.volume_names+self.arterial_names
        self.volume_count=len(self.volume_names)
        self.U=np.r_[self.base.inventory.unstressed_volumes,np.zeros(len(self.beds))]
        self.C={r.name:r.compliance_mL_mmHg for r in self.base.inventory.reservoirs if r.compliance_mL_mmHg is not None}
        v=np.r_[self.base.initial_state[:11],np.zeros(len(self.beds))]
        self.edges=[]
        for name in VALVES:self.edges.append((name,*EDGES[VALVES.index(name)]))
        self.edges.append(('pulmonary_exchange','PA','PU'))
        for b in self.base.branches:
            if b.name in self.beds:
                a='A_'+b.name
                self.edges.extend([(b.name,b.donor,a),(b.name+'_exchange',a,b.recipient),(b.name+'_collateral',b.donor,a)])
            else:self.edges.append((b.name,b.donor,b.recipient))
        self.edges.extend([('support_right','VC','PA'),('support_left','PU','AO')])
        self.flow_names=tuple(e[0] for e in self.edges)
        self.incidence=np.zeros((self.volume_count+1,len(self.edges)))
        nodes=self.volume_names+('H',)
        for j,(_,a,b) in enumerate(self.edges):
            self.incidence[nodes.index(a),j]=-1; self.incidence[nodes.index(b),j]=1
        for bed,pool,a in zip(self.beds,self.pool_names,self.arterial_names):
            i=self.volume_names.index(pool); j=self.volume_names.index(a)
            oldU=self.U[i]; oldV=v[i]; oldC=self.C[pool]
            self.U[i]=(1-arterial_fraction)*oldU; self.U[j]=arterial_fraction*oldU
            v[i]=(1-arterial_fraction)*oldV; v[j]=arterial_fraction*oldV
            # Reference construction: AO=100, distal pool=5 mmHg. Divide the
            # original 95-mmHg inlet drop between gate and exchange resistors.
            pa=100-inlet_resistance_fraction*95
            self.C[a]=arterial_fraction*(oldV-oldU)/pa
            self.C[pool]=(1-arterial_fraction)*oldC
        self.initial_state=np.r_[v,self.base.initial_state[11:]]
        self.controls=dict(inlet=1.,outlet=1.,collateral=0.,right_head=0.,left_head=0.,pump_on=False)
        self.ramp=None

    def set_controls(self, *, time, ramp_s=0., **changes):
        """Configure a new segment; caller integrates exactly to time first."""
        if not math.isfinite(time) or not math.isfinite(ramp_s) or ramp_s<0:raise ValueError('Invalid event time/ramp')
        old=self.control_at(time)
        if set(changes)-set(old):raise ValueError('Unknown control')
        new={**old,**changes}
        for n in ('inlet','outlet','collateral'):
            if not math.isfinite(new[n]) or not 0<=new[n]<=1:raise ValueError('Conductance fraction must be 0..1')
        for n in ('right_head','left_head'):
            if not math.isfinite(new[n]) or new[n]<0:raise ValueError('Nonnegative finite head required')
        if not isinstance(new['pump_on'],bool):raise ValueError('pump_on must be boolean')
        self.controls=new
        self.ramp=(time,ramp_s,old,new) if ramp_s else None

    def control_at(self,t):
        if self.ramp is None:return dict(self.controls)
        start,duration,old,new=self.ramp
        f=min(1.,max(0.,(t-start)/duration))
        # Pump release is discrete; its head is ramped. Inlet/outlet controls
        # are continuous conductance ramps. An off port is always closed.
        return {k:new[k] if k=='pump_on' else old[k]+f*(new[k]-old[k]) for k in new}

    def domain_margins(self,state):
        return dict(zip(self.volume_names,map(float,np.asarray(state)[:self.volume_count]-self.U)))

    def pressures(self,t,state,active=True):
        y=np.asarray(state)
        if y.shape!=(self.volume_count+4,) or not np.all(np.isfinite(y)):raise ValueError('Invalid selective state')
        # Parent mechanics only need the six central volumes and four flows.
        old=np.r_[y[:11].copy(),y[self.volume_count:]]
        for pool,a in zip(self.pool_names,self.arterial_names):old[self.base.volume_names.index(pool)]+=y[self.volume_names.index(a)]
        p,meta=self.base.pressures(t,old,active=active)
        for name in self.pool_names+self.arterial_names:
            i=self.volume_names.index(name); p[name]=(y[i]-self.U[i])/self.C[name]
        return p,meta

    def flows(self,t,state,open_valves,active=True,pressures=None):
        p=self.pressures(t,state,active=active)[0] if pressures is None else pressures
        old=np.r_[np.asarray(state)[:11],np.asarray(state)[self.volume_count:]]
        base=self.base.flows(t,old,open_valves,active=active,pressures=p)
        c=self.control_at(t); q={}
        for name,donor,recipient in self.edges:
            if name in self.beds:
                r=self.base.branch_by_name[name].resistance*self.inlet_fraction
                q[name]=c['inlet']*(p[donor]-p[recipient])/r
            elif name.endswith('_exchange') and name!='pulmonary_exchange':
                bed=name.removesuffix('_exchange'); r=self.base.branch_by_name[bed].resistance*(1-self.inlet_fraction)
                q[name]=(p[donor]-p[recipient])/r
            elif name.endswith('_collateral'):
                bed=name.removesuffix('_collateral'); r=self.base.branch_by_name[bed].resistance*self.inlet_fraction
                q[name]=c['collateral']*(p[donor]-p[recipient])/r
            elif name.startswith('limb_') and name.endswith('_return'):
                q[name]=c['outlet']*base[name]
            elif name.startswith('support_'):
                head=c['right_head'] if name=='support_right' else c['left_head']
                q[name]=finite_head_flow(p[donor],p[recipient],head,.1,120.) if c['pump_on'] else 0.
            else:q[name]=base[name]
        return q

    def rhs(self,t,state,open_valves,active=True):
        p,_=self.pressures(t,state,active=active); q=self.flows(t,state,open_valves,active=active,pressures=p)
        dv=self.incidence@np.array([q[n] for n in self.flow_names])
        if abs(dv[-1])>1e-9:raise RuntimeError('Hepatic balance failed')
        dq=np.zeros(4)
        for i,(a,b) in enumerate(EDGES[:4]):
            if open_valves[i]:dq[i]=(p[a]-p[b]-self.base.source.resistances[i]*q[VALVES[i]])/self.base.source.inertances[i]
        return np.r_[dv[:-1],dq]


def transport_rhs(volumes,amounts,edges,flows):
    """Conservative passive-marker amount transfer; signed flow uses its donor.

    For verification only: no oxygen reactions, anticoagulant pharmacology,
    coagulation or cell-damage kinetics are inferred from this passive marker.
    Edges here must connect physical nodes (eliminate zero-storage H first).
    """
    volumes=np.asarray(volumes); amounts=np.asarray(amounts)
    if np.any(volumes<=0):raise ValueError('Positive mixing volumes required')
    derivative=np.zeros_like(amounts,dtype=float)
    for (a,b),flow in zip(edges,flows):
        donor=a if flow>=0 else b
        transfer=flow*amounts[donor]/volumes[donor]
        derivative[a]-=transfer; derivative[b]+=transfer
    return derivative
