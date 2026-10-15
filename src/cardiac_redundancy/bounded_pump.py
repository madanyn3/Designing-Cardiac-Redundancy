"""T5.1 assumed preload-limited extension of the verified abstract pump.

The gate is a configurable constitutive assumption, not a human suction or
venous-collapse threshold. It bounds available capacity without prescribing
flow or bypassing the circulation's physical volume balances.
"""
from dataclasses import dataclass
import math
from .backup_pump import BackupPumpModel, PumpStage
from .selective import finite_head_flow


@dataclass(frozen=True)
class PreloadGate:
    closed_below_mmHg: float = 0.
    full_above_mmHg: float = 4.
    external_pressure_mmHg: float = 0.

    def __post_init__(self):
        if any(isinstance(x,bool) or not isinstance(x,(int,float)) or not math.isfinite(x) for x in vars(self).values()):
            raise ValueError('Finite numeric preload parameters required')
        if self.full_above_mmHg<=self.closed_below_mmHg:
            raise ValueError('Full-capacity threshold must exceed closure threshold')

    def fraction(self,inlet_pressure_mmHg):
        if not math.isfinite(inlet_pressure_mmHg):raise ValueError('Finite inlet pressure required')
        preload=inlet_pressure_mmHg-self.external_pressure_mmHg
        return min(1.,max(0.,(preload-self.closed_below_mmHg)/(self.full_above_mmHg-self.closed_below_mmHg)))


def bounded_stage_flow(pin,pout,head,stage=PumpStage(),gate=PreloadGate(),enabled=True):
    """Q=min(gate*capacity, max(0,(Pin+head-Pout)/R)) [mL/s].

    External pressure affects inlet availability, not the port pressure rise.
    Zero command head still allows passive forward flow if enabled and filled.
    Disable/isolate explicitly. No efficiency or metabolic claim is made here.
    """
    if any(not math.isfinite(x) for x in (pin,pout,head)) or head<0:
        raise ValueError('Finite pressures and nonnegative pump head required')
    if not isinstance(stage,PumpStage) or not isinstance(gate,PreloadGate):
        raise ValueError('Validated stage and gate required')
    return finite_head_flow(pin,pout,head,stage.resistance_mmHg_s_mL,
                            stage.capacity_mL_s*gate.fraction(pin)) if enabled else 0.


class PreloadPumpModel(BackupPumpModel):
    """B1 adapter: inherit state, routing, events and regional controls unchanged.

    Right VC→PA and left PU→AO stages have independent gates and capacities.
    Use actual stage flow when computing energy; the historical pump_power
    helper reconstructs the ungated law and is not valid for this subclass.
    """
    def __init__(self,*args,right_gate=PreloadGate(),left_gate=PreloadGate(),**kwargs):
        super().__init__(*args,**kwargs)
        if not isinstance(right_gate,PreloadGate) or not isinstance(left_gate,PreloadGate):
            raise ValueError('Validated preload gates required')
        self.right_gate=right_gate;self.left_gate=left_gate

    def flows(self,t,state,open_valves,active=True,pressures=None):
        p=self.pressures(t,state,active=active)[0] if pressures is None else pressures
        q=super().flows(t,state,open_valves,active=active,pressures=p)
        controls=self.control_at(t)
        for side,donor,recipient in (('right','VC','PA'),('left','PU','AO')):
            q['support_'+side]=bounded_stage_flow(p[donor],p[recipient],controls[side+'_head'],
                getattr(self.pump_parameters,side),getattr(self,side+'_gate'),controls['pump_on'])
        return q
