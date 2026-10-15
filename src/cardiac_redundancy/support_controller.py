"""Sampled supervisory policy, not a biologically realized controller.

Inputs must be causal filtered observations, never future simulation truth.
Commands do not change blood states; a circulation adapter must integrate to
release events and apply continuous actuator envelopes separately.
"""
from dataclasses import dataclass
import math


@dataclass(frozen=True)
class ControllerParameters:
    trigger_MAP_mmHg: float = 60.
    reset_MAP_mmHg: float = 65.
    native_failure_L_min: float = .5
    recovery_MAP_mmHg: float = 70.
    recovery_native_L_min: float = 2.
    confirmation_s: float = .25
    signal_s: float = .05
    recruitment_s: float = .2
    ramp_s: float = 1.
    recovery_dwell_s: float = 2.
    withdrawal_s: float = 2.
    max_sample_gap_s: float = .1

    def __post_init__(self):
        if any(isinstance(v, bool) or not isinstance(v, (int, float))
               or not math.isfinite(v) or v < 0 for v in vars(self).values()):
            raise ValueError('Finite nonnegative controller parameters required')
        if not self.trigger_MAP_mmHg < self.reset_MAP_mmHg <= self.recovery_MAP_mmHg:
            raise ValueError('Ordered pressure thresholds required')
        if self.recovery_native_L_min <= self.native_failure_L_min:
            raise ValueError('Recovery flow must exceed failure threshold')
        if min(self.ramp_s, self.withdrawal_s, self.max_sample_gap_s) <= 0:
            raise ValueError('Positive ramp, withdrawal and sampling gap required')


class SupportController:
    """MONITOR -> CONFIRM -> ACTIVATE -> SUPPORT -> RECOVER -> MONITOR.

    FAULT is latched and holds the last command (an assumed policy, not proven
    safe). SUPPORT means commanded support; physiological adequacy is separate.
    Recovery requires a supplied adequacy flag AND native-flow/pressure checks.
    """
    def __init__(self, parameters=ControllerParameters()):
        self.p = parameters
        self.state = 'MONITOR'
        self.command = 0.
        self.last_time = None
        self.confirm_since = self.recovery_since = None
        self.release_s = self.ramp_start = self.withdraw_start = None
        self.ramp_from = 0.
        self.events = []
        self.fault_reason = None

    def transition(self, state, t, reason):
        if self.state != state:
            self.events.append(dict(time_s=t, before=self.state, after=state, reason=reason))
            self.state = state

    def update(self, t, mean_aortic_mmHg, native_output_L_min, *,
               sensor_valid=True, model_valid=True, adequate_for_recovery=False):
        if not math.isfinite(t) or t < 0 or (self.last_time is not None and t <= self.last_time):
            raise ValueError('Strictly increasing nonnegative sample times required')
        gap = self.last_time is not None and t-self.last_time > self.p.max_sample_gap_s + 1e-12
        self.last_time = t
        if self.state == 'FAULT':
            return self.command
        valid = (sensor_valid and mean_aortic_mmHg is not None and native_output_L_min is not None
                 and math.isfinite(mean_aortic_mmHg) and math.isfinite(native_output_L_min))
        if gap or not valid or not model_valid:
            self.fault_reason = 'sample_gap' if gap else 'invalid_sensor' if not valid else 'model_domain'
            self.transition('FAULT', t, self.fault_reason)
            return self.command
        p = self.p
        low = mean_aortic_mmHg < p.trigger_MAP_mmHg and native_output_L_min < p.native_failure_L_min
        reset = mean_aortic_mmHg >= p.reset_MAP_mmHg or native_output_L_min >= p.native_failure_L_min
        recovered = (mean_aortic_mmHg >= p.recovery_MAP_mmHg
                     and native_output_L_min >= p.recovery_native_L_min and adequate_for_recovery)
        if self.state == 'MONITOR' and low:
            self.confirm_since = t
            self.transition('CONFIRM', t, 'low_pressure_and_native_flow')
        if self.state == 'CONFIRM':
            if reset:
                self.confirm_since = None
                self.transition('MONITOR', t, 'confirmation_cancelled')
            elif t-self.confirm_since >= p.confirmation_s-1e-12:
                # Confirmation uses observed sample time, never backdates a
                # decision to the unobserved instant a timer expired.
                self.release_s = t+p.signal_s+p.recruitment_s
                self.ramp_start = self.release_s
                self.ramp_from = self.command
                self.transition('ACTIVATE', t, 'confirmed')
        if self.state == 'ACTIVATE':
            fraction = min(1., max(0., (t-self.ramp_start)/p.ramp_s))
            self.command = self.ramp_from+(1-self.ramp_from)*fraction
            if fraction >= 1-1e-12:
                self.command = 1.
                self.transition('SUPPORT', t, 'ramp_complete_not_adequacy')
        if self.state == 'SUPPORT':
            if not recovered:
                self.recovery_since = None
            elif self.recovery_since is None:
                self.recovery_since = t
            if self.recovery_since is not None and t-self.recovery_since >= p.recovery_dwell_s-1e-12:
                self.withdraw_start = t
                self.transition('RECOVER', t, 'sustained_native_recovery')
        if self.state == 'RECOVER':
            self.command = max(0., 1-(t-self.withdraw_start)/p.withdrawal_s)
            if not recovered:
                # Abort withdrawal on any loss of recovery qualification.
                # Ramp from the current command, avoiding an instantaneous jump.
                self.ramp_from = self.command
                self.ramp_start = self.release_s = t
                self.recovery_since = None
                self.transition('ACTIVATE', t, 'recovery_lost')
            elif self.command <= 0:
                self.recovery_since = self.confirm_since = None
                self.release_s = None
                self.transition('MONITOR', t, 'withdrawal_complete')
        return self.command
