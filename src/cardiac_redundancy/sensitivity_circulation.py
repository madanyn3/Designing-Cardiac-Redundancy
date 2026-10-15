"""M7 explicit coronary-conductance sensitivity on the frozen M6 circuit.

This multiplier is an independent anatomical/design hypothesis, not a measured
remodeling law. Both auxiliary coronary flows remain in the incidence matrix,
so their volume transfers are conservative. Native coronary flow is unchanged
at a common prescribed pressure/state; coupled pressures may subsequently move.
"""
import math
from .controlled_circulation import ControlledCirculation

class SensitivityCirculation(ControlledCirculation):
    def __init__(self,*args,coronary_conductance=1.,**kwargs):
        if isinstance(coronary_conductance,bool) or not math.isfinite(coronary_conductance) or coronary_conductance<=0:
            raise ValueError('Positive finite coronary conductance multiplier required')
        self.coronary_conductance=coronary_conductance
        super().__init__(*args,**kwargs)

    def flows(self,t,state,open_valves,active=True,pressures=None):
        q=super().flows(t,state,open_valves,active,pressures)
        for side in ('left','right'):
            q['aux_coronary_'+side]*=self.coronary_conductance
        return q
