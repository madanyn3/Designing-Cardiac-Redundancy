"""Explicit organ demand schedules for exploratory preservation scenarios.

Fractions prescribe a hypothetical metabolic state; reduced flow never selects
that state automatically. Positive targets avoid treating excluded tissue as
free to discard. Auxiliary demand is calculated independently from pump work.
"""
import math
import numpy as np
from .preservation_energy import native_fraction, integrated_native_fraction

class OrganPreservation:
    def __init__(self, organs, entries):
        self.organs = tuple(organs)
        if set(entries) != {o.name for o in organs}:
            raise ValueError('Profile must cover every organ exactly once')
        self.entries = tuple(dict(entries[o.name]) for o in organs)
        for o, e in zip(organs, self.entries):
            if e['role'] not in ('active', 'preservation', 'pump'):
                raise ValueError('Unknown organ role')
            native_fraction(0, e['fraction'], e['ramp_s'])
            if (o.parent_pool == 'auxiliary_heart') != (e['role'] == 'pump'):
                raise ValueError('Auxiliary myocardium requires independent pump role')
            if e['role'] == 'pump' and (e['fraction'] != 1 or e['ramp_s'] != 0):
                raise ValueError('Do not discount the workload-derived pump demand')
        self.reference = np.array([o.demand_mL_min for o in organs])

    def demand(self, time_s):
        """Requested oxygen in mL/min, independent of delivered supply."""
        return self.reference * np.array([native_fraction(time_s,e['fraction'],e['ramp_s']) for e in self.entries])

    def integrated(self, time_s):
        """Analytic original-demand integral in mL O2, including transitions."""
        return self.reference / 60 * np.array([integrated_native_fraction(time_s,e['fraction'],e['ramp_s']) for e in self.entries])

    def groups(self):
        return {role: [i for i,e in enumerate(self.entries) if e['role']==role]
                for role in ('active','preservation','pump')}
