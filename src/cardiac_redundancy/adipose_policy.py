"""Arm B support policy: whether visceral adipose is perfused or suppressed.

This is a **policy** choice over the subject's existing vasculature, not a
design change to it. The depot's branch resistance is anatomy and stays fixed;
what the support policy elects is the fraction of that conductance the
circulation opens, exactly as it already does for brain, coronary, protected
kidney, hepatic arterial and respiratory muscle.

The choice is categorical rather than continuous: adipose either joins the
mandatory set at full conductance or stays with the preserved territories at
the declared preservation conductance. Treating it as a continuous knob would
invite exactly the error the two-column registry exists to prevent, since any
value above full conductance would mean dilating a vessel beyond its
physiological maximum.

The hydraulic and metabolic sides must agree. A depot perfused at full
conductance but still carrying a suppressed metabolic demand, or the reverse,
would be an inconsistency no single result would reveal, so both are derived
here from one flag.
"""
from dataclasses import dataclass

ADIPOSE_BRANCH = 'visceral_adipose'


@dataclass(frozen=True)
class AdiposePolicy:
    """One flag, applied consistently to circulation and metabolism."""
    perfused: bool = False

    def __post_init__(self):
        if not isinstance(self.perfused, bool):
            raise ValueError('Adipose policy requires an explicit boolean')

    @property
    def additional_mandatory(self):
        """Branches to add to the mandatory set in apply_support_controls."""
        return (ADIPOSE_BRANCH,) if self.perfused else ()

    def organ_entry(self, preservation_fraction):
        """The OrganPreservation entry for the depot under this policy.

        Perfused, the depot carries its full demand and is counted as an active
        territory. Suppressed, it is reduced like any other preserved tissue.
        """
        if self.perfused:
            return dict(role='active', fraction=1., ramp_s=0.)
        return dict(role='preservation', fraction=float(preservation_fraction), ramp_s=0.)

    def conductance_factor(self, preservation_conductance):
        """The fraction of the depot's own conductance that is opened."""
        return 1. if self.perfused else float(preservation_conductance)

    def describe(self):
        return dict(perfused=self.perfused, branch=ADIPOSE_BRANCH,
                    meaning=('visceral adipose joins the mandatory set at full conductance'
                             if self.perfused else
                             'visceral adipose stays with the preserved territories'),
                    scope='a support-policy choice over existing vasculature; the depot '
                          'resistance is unchanged in either case')


def apply(model, policy, time=0., preservation_conductance=.25):
    """Set regional controls for this policy, returning the resulting factors."""
    if not isinstance(policy, AdiposePolicy):
        raise ValueError('Validated adipose policy required')
    model.apply_support_controls(time, preservation_conductance,
                                 additional_mandatory=policy.additional_mandatory)
    return dict(model.native.regional_factors)


def organ_entries(entries, policy, preservation_fraction):
    """Add or update the depot's entry in a metabolic preservation profile."""
    if not isinstance(policy, AdiposePolicy):
        raise ValueError('Validated adipose policy required')
    out = {k: dict(v) for k, v in entries.items()}
    out[ADIPOSE_BRANCH] = policy.organ_entry(preservation_fraction)
    return out
