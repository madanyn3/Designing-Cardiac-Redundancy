"""Configurable finite-head support ports; no pump-size or efficiency model.

This extension leaves historical fixture implementations unchanged. Resistance
and capacity are fixed for a run, while existing head ramps and release controls
remain available. A zero capacity independently isolates either support stage.
"""
from dataclasses import dataclass
import json
import math
from pathlib import Path

from .selective import finite_head_flow
from .whole_network import WholeNetworkModel


@dataclass(frozen=True)
class PumpStage:
    """Internal hydraulic resistance [mmHg s/mL] and capacity [mL/s]."""

    resistance_mmHg_s_mL: float = .1
    capacity_mL_s: float = 120.

    def __post_init__(self):
        for value in (self.resistance_mmHg_s_mL, self.capacity_mL_s):
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
                raise ValueError('Pump parameters must be finite numbers')
        if self.resistance_mmHg_s_mL <= 0 or self.capacity_mL_s < 0:
            raise ValueError('Pump resistance must be positive and capacity nonnegative')


@dataclass(frozen=True)
class PumpParameters:
    right: PumpStage = PumpStage()
    left: PumpStage = PumpStage()

    def __post_init__(self):
        if not isinstance(self.right, PumpStage) or not isinstance(self.left, PumpStage):
            raise ValueError('Both support stages require PumpStage parameters')

    @classmethod
    def from_json(cls, path):
        """Require explicit units and both stages; reject silent spelling errors."""
        data = json.loads(Path(path).read_text(encoding='utf-8'))
        fields = {'resistance_mmHg_s_mL', 'capacity_mL_s'}
        if not isinstance(data, dict) or set(data) != {'right', 'left'}:
            raise ValueError('Pump configuration must contain exactly right and left')
        if any(not isinstance(s, dict) or set(s) != fields for s in data.values()):
            raise ValueError('Each stage requires resistance_mmHg_s_mL and capacity_mL_s')
        return cls(**{name: PumpStage(**stage) for name, stage in data.items()})


class BackupPumpModel(WholeNetworkModel):
    """Whole-network hydraulics with independently parameterized support stages."""

    def __init__(self, *args, pump_parameters=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.pump_parameters = PumpParameters() if pump_parameters is None else pump_parameters
        if not isinstance(self.pump_parameters, PumpParameters):
            raise ValueError('Expected validated PumpParameters')

    def flows(self, t, state, open_valves, active=True, pressures=None):
        p = self.pressures(t, state, active=active)[0] if pressures is None else pressures
        q = super().flows(t, state, open_valves, active=active, pressures=p)
        controls = self.control_at(t)
        # Only replace support-edge constitutive laws. The inherited incidence
        # matrix still transfers every mL from VC to PA or from PU to AO.
        for side, donor, recipient in (('right', 'VC', 'PA'), ('left', 'PU', 'AO')):
            stage = getattr(self.pump_parameters, side)
            q['support_' + side] = finite_head_flow(
                p[donor], p[recipient], controls[side + '_head'],
                stage.resistance_mmHg_s_mL, stage.capacity_mL_s,
            ) if controls['pump_on'] else 0.
        return q
