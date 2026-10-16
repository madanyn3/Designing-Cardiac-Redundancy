"""Ideal bounded cerebral flow regulation as a structural comparator.

No response delay, CO2 or ICP dynamics are claimed. Closing pressure and signed
passivity remain in the inherited flow law; the incidence matrix transports
the modified flow conservatively between the same physical volumes.
"""
import math
from .sensitivity_circulation import SensitivityCirculation


def regulated_flow(passive_flow, target, maximum_gain):
    if not all(math.isfinite(v) for v in (passive_flow, target, maximum_gain)) or target <= 0 or maximum_gain < 1:
        raise ValueError('Finite flow, positive target and maximum gain >=1 required')
    if passive_flow <= 0:
        return passive_flow
    gain = min(maximum_gain, max(.5, target/passive_flow))
    return passive_flow*gain


class RegulatedCirculation(SensitivityCirculation):
    def __init__(self, *args, regulation_gain=None, **kwargs):
        if regulation_gain is not None:
            regulated_flow(1., 1., regulation_gain)
        self.regulation_gain = regulation_gain
        super().__init__(*args, **kwargs)
        self.cerebral_target = next(r['target_flow_mL_s'] for r in self.native.base.project['branches'] if r['name']=='brain')

    def flows(self, t, state, open_valves, active=True, pressures=None):
        q = super().flows(t, state, open_valves, active, pressures)
        if self.regulation_gain is not None:
            q['brain'] = regulated_flow(q['brain'], self.cerebral_target, self.regulation_gain)
        return q
