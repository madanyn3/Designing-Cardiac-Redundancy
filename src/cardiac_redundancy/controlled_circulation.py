"""Continuous actuator envelopes on the frozen bilateral circulation.

Native restoration is a prescribed scenario, not a controller observation.
Changing an envelope does not change any volume or valve-flow state.
"""
import math
from .bilateral_support import BilateralSupport
from .mechanics import septal_root
from .source import driver


class Envelope:
    def __init__(self, initial=0.):
        self.initial = initial
        self.events = []

    def value(self, t):
        event = next((e for e in reversed(self.events) if e[0] <= t), None)
        if event is None: return self.initial
        start, duration, old, target = event
        f = 1. if duration == 0 else min(1., max(0., (t-start)/duration))
        return old+(target-old)*f

    def set(self, t, target, duration):
        if not all(math.isfinite(v) for v in (t,target,duration)) or duration < 0 or not 0 <= target <= 1:
            raise ValueError('Invalid actuator envelope')
        if self.events and t < self.events[-1][0]: raise ValueError('Unordered actuator event')
        self.events.append((t,duration,self.value(t),target))


class ControlledCirculation(BilateralSupport):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.aux_envelope = Envelope()
        self.native_envelope = Envelope()

    def pressures(self,t,state,active=True):
        # Reuse vascular/pericardial pressures and regional storage exactly.
        # Both free walls AND septum receive the same scaled driver.
        p, meta = super().pressures(t,state,active=False)
        source = self.native.base.source
        d = source.config['driver']
        e = self.native_envelope.value(t)*driver(t,d['period_s'],d['A'],d['B_per_s2'],d['C_s'])
        lv,rv = state[:2]
        septum,residual = septal_root(lv,rv,e,*source.walls)
        p['LV'] = source.walls[0].pressure(lv-septum,e)+p['PERI']
        p['RV'] = source.walls[1].pressure(rv+septum,e)+p['PERI']
        meta.update(activation=e,septal_residual_mmHg=residual,septal_coordinate_mL=septum)
        elapsed = 0. if self.trigger is None else max(0.,t-self.trigger)
        e = self.aux_envelope.value(t)*math.exp(-45*((elapsed%self.auxiliary_period_s)/self.auxiliary_period_s-.5)**2)
        lv,rv = state[16:18]
        septum,residual = septal_root(lv,rv,e,*self.auxiliary_walls)
        p['aLV'] = self.auxiliary_walls[0].pressure(lv-septum,e)+p['aPERI']
        p['aRV'] = self.auxiliary_walls[1].pressure(rv+septum,e)+p['aPERI']
        meta.update(aux_activation=e,aux_septal_residual_mmHg=residual,aux_septal_coordinate_mL=septum)
        return p,meta
