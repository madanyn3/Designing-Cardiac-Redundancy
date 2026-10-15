"""Prescribed arrest-to-activation timing; not a validated biological sensor.

Detection/confirmation, signal transit, and recruitment delays are distinct.
A simulation must stop at release_time and call release there: enabling valves
in advance would incorrectly alter pre-trigger flow. Effective support remains
an output of pressure/flow/substrate checks, not the end of the ramp.
"""
from dataclasses import dataclass
import math

@dataclass(frozen=True)
class ActivationTiming:
    detection_s: float = .25
    signal_s: float = .05
    recruitment_s: float = .2
    ramp_s: float = 1.

    def __post_init__(self):
        if any(isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v) or v<0 for v in vars(self).values()) or self.ramp_s<=0:
            raise ValueError('Finite nonnegative delays and positive ramp required')

    @property
    def latency_s(self):return self.detection_s+self.signal_s+self.recruitment_s

    def release_time(self,arrest_s):
        if not math.isfinite(arrest_s):raise ValueError('Finite arrest time required')
        return arrest_s+self.latency_s

    def amplitude(self,elapsed_s):
        if not math.isfinite(elapsed_s):raise ValueError('Finite time required')
        return min(1.,max(0.,(elapsed_s-self.latency_s)/self.ramp_s))

    def release(self,model,arrest_s,current_time_s):
        if not math.isfinite(current_time_s) or not math.isclose(current_time_s,self.release_time(arrest_s),rel_tol=0,abs_tol=1e-10):
            raise ValueError('Integrate to the exact release event before opening valves')
        if model.architecture!='BIO' or model.trigger is not None:
            raise ValueError('Unreleased BIO model required')
        model.auxiliary_ramp_s=self.ramp_s
        model.release_auxiliary(current_time_s)
