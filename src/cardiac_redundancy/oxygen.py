"""M4 foundation: conservative blood oxygen transport and a finite gas pool.

Amounts are mL O2 STPD, blood volumes mL, flows mL/s, and time seconds.
This stage deliberately has ZERO tissue consumption. It verifies transport,
not organ support. See docs/m4-oxygen-specification.md for sources and limits.
"""
from dataclasses import dataclass
import math

import numpy as np
from scipy.optimize import brentq

from .transport import marker_rhs


@dataclass(frozen=True)
class OxygenParameters:
    hemoglobin_g_dL: float = 15.0
    binding_mL_g: float = 1.34
    dissolved_mL_dL_mmHg: float = 0.0031
    p50_mmHg: float = 26.8
    hill_alpha: float = 2.82
    hill_beta: float = 1.20
    hill_gamma_mmHg: float = 29.25
    alveolar_volume_mL_BTPS: float = 2500.0
    barometric_mmHg: float = 760.0
    water_vapor_mmHg: float = 47.0
    temperature_K: float = 310.15
    inspired_fraction: float = 0.21
    exchange_fraction: float = 1.0

    def __post_init__(self):
        if not all(math.isfinite(v) for v in vars(self).values()):
            raise ValueError('Oxygen parameters must be finite')
        positive = (self.binding_mL_g, self.dissolved_mL_dL_mmHg,
                    self.p50_mmHg, self.hill_gamma_mmHg,
                    self.alveolar_volume_mL_BTPS, self.temperature_K)
        if min(positive) <= 0 or self.hemoglobin_g_dL < 0:
            raise ValueError('Invalid oxygen content/store parameters')
        if not 0 <= self.water_vapor_mmHg < self.barometric_mmHg:
            raise ValueError('Invalid gas pressures')
        if not 0 <= self.hill_beta < self.hill_alpha:
            raise ValueError('Invalid saturation curve')
        if not (0 <= self.inspired_fraction <= 1 and 0 <= self.exchange_fraction <= 1):
            raise ValueError('Gas fractions must be in [0, 1]')

    @property
    def btps_to_stpd(self):
        """Ideal-gas conversion to dry 0 C, 760 mmHg reference volume."""
        return (self.barometric_mmHg-self.water_vapor_mmHg)/760 * 273.15/self.temperature_K

    @property
    def alveolar_volume_mL_STPD(self):
        return self.alveolar_volume_mL_BTPS * self.btps_to_stpd


def saturation(po2_mmHg, p=OxygenParameters()):
    """Dash et al. 2015 variable-Hill curve at fixed standard conditions.

    Only the standard-condition O2 curve is implemented, not the paper's
    CO2/pH/temperature-dependent model. Zero pressure has exact zero saturation.
    """
    po2 = np.asarray(po2_mmHg, dtype=float)
    if np.any(~np.isfinite(po2)) or np.any(po2 < 0):
        raise ValueError('PO2 must be finite and nonnegative')
    exponent = p.hill_alpha-p.hill_beta*10**(-po2/p.hill_gamma_mmHg)
    ratio = (po2/p.p50_mmHg)**exponent
    return ratio/(1+ratio)


def oxygen_content(po2_mmHg, p=OxygenParameters()):
    """mL O2 STPD / mL blood; divide the published mL/dL formula by 100."""
    return (p.binding_mL_g*p.hemoglobin_g_dL*saturation(po2_mmHg, p)
            + p.dissolved_mL_dL_mmHg*np.asarray(po2_mmHg))/100


def oxygen_pressure(content_mL_mL, p=OxygenParameters()):
    """Invert content numerically, including dissolved oxygen (scalar input)."""
    c = float(content_mL_mL)
    if not math.isfinite(c) or c < 0:
        raise ValueError('Content must be finite and nonnegative')
    if c == 0:
        return 0.0
    # Dissolved oxygen alone reaches c at this upper bracket.
    upper = max(1.0, 100*c/p.dissolved_mL_dL_mmHg)
    return brentq(lambda pressure: float(oxygen_content(pressure, p))-c, 0, upper)


class ConservativeOxygenTransport:
    """One-way coupling to hydraulic volumes/flows; no organ uptake yet.

    State: oxygen amount per physical blood compartment, alveolar oxygen,
    cumulative NET external oxygen exchange. H is an algebraic mixer and
    receives no storage. Pump edges transport existing oxygen only.
    """
    def __init__(self, names, edges, parameters=OxygenParameters()):
        self.names = tuple(names)
        self.edges = tuple(edges)
        self.parameters = parameters
        if len(set(names)) != len(names) or 'H' in names or not names:
            raise ValueError('Unique physical compartment names required')
        edge_names = [edge[0] for edge in edges]
        if len(set(edge_names)) != len(edge_names):
            raise ValueError('Unique edge names required')
        for _, a, b in edges:
            if a not in (*names, 'H') or b not in (*names, 'H'):
                raise ValueError('Unknown hydraulic endpoint')
        if edge_names.count('pulmonary_exchange') != 1:
            raise ValueError('Exactly one pulmonary exchange edge required')
        self.pulmonary_index = edge_names.index('pulmonary_exchange')
        if tuple(edges[self.pulmonary_index][1:]) != ('PA', 'PU'):
            raise ValueError('Pulmonary exchange must connect PA to PU')

    def initial_state(self, volumes, po2_mmHg=40.0, alveolar_fraction=0.14):
        """Uniform PO2 default is a TEST initial condition, not a resting state."""
        if not math.isfinite(alveolar_fraction) or not 0 <= alveolar_fraction <= 1:
            raise ValueError('Invalid alveolar fraction')
        volumes = np.asarray(volumes, dtype=float)
        if volumes.shape != (len(self.names),) or np.any(~np.isfinite(volumes)) or np.any(volumes <= 0):
            raise ValueError('Positive finite physical volumes required')
        amounts = volumes*oxygen_content(po2_mmHg, self.parameters)
        if amounts.shape != volumes.shape:
            raise ValueError('Initial PO2 shape mismatch')
        return np.r_[amounts, self.parameters.alveolar_volume_mL_STPD*alveolar_fraction, 0.0]

    def rhs(self, volumes, state, flows, ventilation_mL_s_BTPS=0.0):
        n = len(self.names)
        volumes = np.asarray(volumes, dtype=float)
        state = np.asarray(state, dtype=float)
        flows = np.asarray(flows, dtype=float)
        if volumes.shape != (n,) or state.shape != (n+2,) or flows.shape != (len(self.edges),):
            raise ValueError('Hydraulic/oxygen state shape mismatch')
        if any(np.any(~np.isfinite(a)) for a in (volumes, state, flows)) or np.any(volumes <= 0):
            raise ValueError('Finite states and positive volumes required')
        if not math.isfinite(ventilation_mL_s_BTPS) or ventilation_mL_s_BTPS < 0:
            raise ValueError('Ventilation must be finite and nonnegative')
        p = self.parameters
        gas_fraction = state[n]/p.alveolar_volume_mL_STPD
        if np.any(state[:n] < 0) or not 0 <= gas_fraction <= 1:
            raise ValueError('Oxygen store outside model domain')
        # Existing signed upwind transport also solves instantaneous H mixing.
        blood_rate = marker_rhs(volumes, state[:n], self.edges, flows, self.names)
        q = flows[self.pulmonary_index]
        donor, recipient = ('PA', 'PU') if q >= 0 else ('PU', 'PA')
        donor_content = state[self.names.index(donor)]/volumes[self.names.index(donor)]
        alveolar_po2 = (p.barometric_mmHg-p.water_vapor_mmHg)*gas_fraction
        # Positive transfer is gas -> blood; negative allows blood -> gas.
        # Reverse-flow equilibration is an explicit mathematical assumption.
        transfer = abs(q)*p.exchange_fraction*(oxygen_content(alveolar_po2, p)-donor_content)
        blood_rate[self.names.index(recipient)] += transfer
        # Fixed-volume gas replacement is an O2-only closure, not a CO2 model.
        external = ventilation_mL_s_BTPS*p.btps_to_stpd*(p.inspired_fraction-gas_fraction)
        return np.r_[blood_rate, external-transfer, external]

    def balance(self, state):
        """Constant without consumption: blood + lung - net external input."""
        state = np.asarray(state)
        return state[..., :-1].sum(axis=-1)-state[..., -1]
