"""Flow-limited organ uptake with conservative serial hepatic transport.

This is an exploratory extraction-cap closure, not a tissue viability model.
Demands are independent of hydraulic controls. See docs/m4-organ-uptake.md.
"""
from dataclasses import dataclass
import math
import numpy as np

from .oxygen import ConservativeOxygenTransport, OxygenParameters
from .transport import marker_rhs


@dataclass(frozen=True)
class OrganDemand:
    name: str
    edge: str
    demand_mL_min: float
    extraction_limit: float
    parent_pool: str

    def __post_init__(self):
        if not math.isfinite(self.demand_mL_min) or self.demand_mL_min < 0:
            raise ValueError('Demand must be finite and nonnegative')
        if not math.isfinite(self.extraction_limit) or not 0 <= self.extraction_limit <= 1:
            raise ValueError('Extraction limit must be in [0, 1]')


def map_demands(budget, mapping, flow_names, extraction_limit, auxiliary_total_mL_min=None):
    """Partition every M1 pool once; mapping weights are fixed, never live flow.

    BIO requires an explicit auxiliary demand input. Its hydraulic scale alone
    cannot determine its metabolism. Caller may supply an assumed demand but
    must report it as such. Native demand is retained independently.
    """
    if set(mapping) != set(budget):
        raise ValueError('Mapping must cover exactly the supplied metabolic pools')
    organs = []
    for pool, weights in mapping.items():
        if not weights or any(not math.isfinite(w) or w <= 0 for w in weights.values()):
            raise ValueError('Positive finite partition weights required')
        if not math.isclose(sum(weights.values()), 1., abs_tol=1e-12):
            raise ValueError('Each metabolic pool must partition to one')
        for edge, weight in weights.items():
            organs.append(OrganDemand(edge, edge, budget[pool]*weight, extraction_limit, pool))
    aux = {'aux_coronary_left', 'aux_coronary_right'} & set(flow_names)
    if aux:
        if len(aux) != 2 or auxiliary_total_mL_min is None:
            raise ValueError('Biological pump requires explicit auxiliary myocardial demand')
        for edge, weight in mapping['heart'].items():
            organs.append(OrganDemand('aux_'+edge, 'aux_'+edge,
                auxiliary_total_mL_min*weight, extraction_limit, 'auxiliary_heart'))
    elif auxiliary_total_mL_min is not None:
        raise ValueError('Auxiliary demand supplied without coronary routes')
    if len({o.edge for o in organs}) != len(organs):
        raise ValueError('Multiple metabolic assignments to one exchange edge')
    if any(o.edge not in flow_names for o in organs):
        raise ValueError('Metabolic route missing from hydraulic model')
    return tuple(organs)


class OrganOxygenTransport(ConservativeOxygenTransport):
    """State: foundation states, cumulative uptake per organ, then unmet demand.

    An uptake edge removes oxygen from the traversing stream, before that
    stream enters its recipient. For H, process all incoming edges first,
    mix their DEPLETED oxygen fluxes, and only then process outgoing edges.
    This same ordering works if the hepatic circulation reverses.
    """
    def __init__(self, names, edges, organs, parameters=OxygenParameters()):
        super().__init__(names, edges, parameters)
        self.organs = tuple(organs)
        self.organ_count = len(self.organs)
        self.by_edge = {o.edge: j for j, o in enumerate(self.organs)}
        if len(self.by_edge) != len(self.organs) or len({o.name for o in organs}) != len(organs):
            raise ValueError('Unique organ names and exchange edges required')
        allowed = {e[0] for e in self.edges} - {'pulmonary_exchange'}
        if not set(self.by_edge) <= allowed:
            raise ValueError('Unknown tissue route or uptake on pulmonary gas exchange')

    def initial_state(self, volumes, po2_mmHg=40., alveolar_fraction=.14):
        base = super().initial_state(volumes, po2_mmHg, alveolar_fraction)
        return np.r_[base, np.zeros(2*self.organ_count)]

    def evaluate(self, volumes, state, flows, ventilation_mL_s_BTPS=0., demand_mL_min=None):
        """Return derivative and instantaneous organ fluxes [mL O2/s].

        Optional demand vector is for explicit future transition schedules;
        it must never be inferred from instantaneous conductance or flow.
        """
        n = len(self.names); m = self.organ_count; base_size = n+2
        state = np.asarray(state, dtype=float); volumes = np.asarray(volumes, dtype=float)
        flows = np.asarray(flows, dtype=float)
        if state.shape != (base_size+2*m,) or np.any(~np.isfinite(state)):
            raise ValueError('Invalid organ oxygen state')
        base_rate = super().rhs(volumes, state[:base_size], flows, ventilation_mL_s_BTPS)
        demand = np.asarray([o.demand_mL_min for o in self.organs] if demand_mL_min is None else demand_mL_min, dtype=float)/60
        if demand.shape != (m,) or np.any(~np.isfinite(demand)) or np.any(demand < 0):
            raise ValueError('Invalid scheduled demand')
        concentrations = dict(zip(self.names, state[:n]/volumes))
        change = dict.fromkeys(self.names, 0.)
        delivered = np.zeros(m); capacity = np.zeros(m); uptake = np.zeros(m)
        oriented = [(name, a if q >= 0 else b, b if q >= 0 else a, abs(q))
                    for (name, a, b), q in zip(self.edges, flows)]
        incoming = sum(q for _, _, b, q in oriented if b == 'H')
        outgoing = sum(q for _, a, _, q in oriented if a == 'H')
        if abs(incoming-outgoing) > 1e-7:
            raise ValueError('Unbalanced hepatic junction')
        incoming_oxygen = 0.

        def move(edge):
            nonlocal incoming_oxygen
            name, donor, recipient, q = edge
            flux = q*concentrations[donor]
            consumed = 0.
            if name in self.by_edge:
                j = self.by_edge[name]
                delivered[j] = flux
                capacity[j] = self.organs[j].extraction_limit*flux
                consumed = uptake[j] = min(demand[j], capacity[j])
            if donor != 'H':
                change[donor] -= flux
            if recipient == 'H':
                incoming_oxygen += flux-consumed
            else:
                change[recipient] += flux-consumed

        # Physical-compartment donors are known; H donors must wait for mixing.
        for edge in oriented:
            if edge[1] != 'H':
                move(edge)
        concentrations['H'] = incoming_oxygen/incoming if incoming else 0.
        for edge in oriented:
            if edge[1] == 'H':
                move(edge)
        # Preserve the verified foundation's gas source; replace only advection.
        passive = marker_rhs(volumes, state[:n], self.edges, flows, self.names)
        base_rate[:n] += np.array([change[name] for name in self.names])-passive
        deficit = demand-uptake
        derivative = np.r_[base_rate, uptake, deficit]
        return derivative, dict(delivered=delivered, extractable=capacity,
            consumed=uptake, demand=demand, unmet=deficit, hepatic_content=concentrations['H'])

    def rhs(self, volumes, state, flows, ventilation_mL_s_BTPS=0., demand_mL_min=None):
        return self.evaluate(volumes, state, flows, ventilation_mL_s_BTPS, demand_mL_min)[0]

    def balance(self, state):
        """Blood + lung - net external input + cumulative ACTUAL consumption."""
        state = np.asarray(state); n = len(self.names); m = self.organ_count
        return state[..., :n+1].sum(axis=-1)-state[..., n+1]+state[..., n+2:n+2+m].sum(axis=-1)
