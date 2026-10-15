"""Smith wall and flow equations, per docs/model-specification.md sections 1–2.

Source attribution: Smith et al. (2004), DOI 10.1016/j.medengphy.2003.10.001;
PMR numerical source and license retained in references/m2. Implementation is
written from the documented equations, not copied from the CVSX reference code.
Units: seconds, mL, mmHg, mL/s throughout.
"""
from dataclasses import dataclass
import math

from scipy.optimize import brentq


class ModelDomainError(ValueError):
    """Invalid input or failed constitutive evaluation; never silently clipped."""


@dataclass(frozen=True)
class Wall:
    """Immutable parameters for one ventricular free wall or the septum.

    elastance: E, active pressure/volume slope [mmHg/mL].
    dead_volume: Vd, zero-pressure intercept of the active law [mL].
    pressure_scale: P0, scale of the passive exponential [mmHg].
    stiffness: lambda, inverse-volume exponential coefficient [1/mL].
    reference_volume: V0, zero-pressure reference of the passive law [mL].

    These parameters describe a lumped pressure law, not muscle anatomy. In
    particular, the wall coordinate supplied to pressure() is not tissue mass
    or a separately conserved blood volume.
    """
    elastance: float
    dead_volume: float
    pressure_scale: float
    stiffness: float
    reference_volume: float

    def __post_init__(self):
        """Reject material parameters that invalidate the monotone wall law."""
        values = (self.elastance, self.dead_volume, self.pressure_scale,
                  self.stiffness, self.reference_volume)
        if not all(math.isfinite(x) for x in values):
            raise ModelDomainError("Nonfinite wall parameter")
        if min(self.elastance, self.pressure_scale, self.stiffness) <= 0:
            raise ModelDomainError("Wall material parameters must be positive")

    def pressure(self, volume, activation):
        """Return wall pressure relative to the pericardium [mmHg].

        F(w,e) = e*E*(w-Vd) + (1-e)*P0*(exp(lambda*(w-V0))-1).
        `volume` is the geometric coordinate w [mL]; it may be negative.
        `activation` is dimensionless e in [0,1]. e=0 leaves passive elastic
        pressure; it does not imply zero chamber pressure or an empty chamber.
        """
        if not math.isfinite(volume) or not 0 <= activation <= 1:
            raise ModelDomainError("Invalid wall coordinate or activation")
        try:
            active = activation * self.elastance * (volume - self.dead_volume)
            # Skip the unused exponential at e=1: evaluating it first could
            # overflow even though its physical coefficient would be zero.
            # expm1(x) evaluates exp(x)-1 accurately when x is near zero.
            passive = 0.0 if activation == 1 else (
                (1 - activation) * self.pressure_scale *
                math.expm1(self.stiffness * (volume - self.reference_volume)))
            result = active + passive
        except OverflowError as exc:
            raise ModelDomainError("Wall exponential overflow") from exc
        if not math.isfinite(result):
            raise ModelDomainError("Nonfinite wall pressure")
        return result


def driver(t, period, amplitude, width, center):
    """Periodic dimensionless Gaussian activation, not electrical conduction.

    t, period and center are seconds; width is 1/s^2, amplitude dimensionless.
    Modulo restarts phase every period, including beyond the first few beats.
    The source has period=.75 s (80/min), center=.375 s and width=80/s^2.
    This function generates contraction timing; it does not integrate a state.
    """
    if not all(math.isfinite(v) for v in (t, period, amplitude, width, center)):
        raise ModelDomainError("Nonfinite driver input")
    if period <= 0 or width <= 0 or not 0 <= amplitude <= 1:
        raise ModelDomainError("Invalid driver parameters")
    return amplitude * math.exp(-width * (t % period - center) ** 2)


def septal_root(lv_volume, rv_volume, activation, left, right, septum,
                tolerance=1e-8, hint=0.0):
    """Return (septal coordinate [mL], pressure residual [mmHg]).

    Find v such that F_spt(v)-F_left(V_LV-v)+F_right(V_RV+v)=0.
    The septum changes both free-wall coordinates with opposite signs. Its
    position is solved algebraically at each pressure evaluation; it is not
    an eleventh state and contributes no additional blood to the inventory.

    The residual derivative is the sum of three positive wall derivatives,
    so there is at most one root for valid parameters. Expand a bracket about
    `hint`, then use Brent's safeguarded scalar solve. A hint changes search
    cost, not the equation. No continuation cache is maintained here; the
    current source RHS uses the default hint=0 on every evaluation.

    `tolerance` bounds the pressure residual, independently of Brent's volume
    coordinate tolerances. An overflowing search evaluation or absent bracket
    fails explicitly; it does not return a stale or clipped root.
    """
    if not math.isfinite(hint) or not math.isfinite(tolerance) or tolerance <= 0:
        raise ModelDomainError("Invalid root settings")

    def residual(v):
        return (septum.pressure(v, activation)
                - left.pressure(lv_volume - v, activation)
                + right.pressure(rv_volume + v, activation))

    half_width = 1.0  # Initial search interval is hint +/- 1 mL.
    for _ in range(20):
        lo, hi = hint - half_width, hint + half_width
        flo, fhi = residual(lo), residual(hi)
        if flo <= 0 <= fhi:
            root = brentq(residual, lo, hi, xtol=1e-13, rtol=1e-14)
            # Convergence in root position alone does not guarantee pressure
            # balance when the exponential wall law is steep.
            error = residual(root)
            if abs(error) > tolerance:
                raise ModelDomainError(f"Septal residual exceeds tolerance: {error}")
            return root, error
        half_width *= 2
    raise ModelDomainError("No finite septal bracket found")


def vascular_pressure(volume, elastance, external=0.0, unstressed=0.0):
    """P = P_external + E*(V-U) [mmHg] for a linear reservoir.

    V and U are mL; E is mmHg/mL and compliance is 1/E. The external pressure
    is added after the elastic term, not multiplied by E. This local relation
    does not impose the later project model's vascular collapse boundary.
    """
    if not all(math.isfinite(v) for v in (volume, elastance, external, unstressed)) or elastance <= 0:
        raise ModelDomainError("Invalid vascular relation")
    return external + elastance * (volume - unstressed)


def resistive_flow(upstream, downstream, resistance):
    """Signed Q=(P_up-P_down)/R [mL/s], with R in mmHg*s/mL.

    This is an ordinary bidirectional resistance, not a one-way valve. Its
    sign tells the network which way blood actually travels along an edge.
    """
    if not all(math.isfinite(v) for v in (upstream, downstream, resistance)) or resistance <= 0:
        raise ModelDomainError("Invalid resistive relation")
    return (upstream - downstream) / resistance


def inertial_derivative(gradient, flow, resistance, inertance):
    """Return dQ/dt=(delta_P-R*Q)/L [mL/s^2] for an open valve.

    gradient is upstream minus downstream pressure [mmHg]; inertance L is
    mmHg*s^2/mL. Positive flow can persist while gradient is negative because
    it takes time to decelerate. The hybrid integrator, not this equation,
    closes the valve when forward flow reaches zero.
    """
    if not all(math.isfinite(v) for v in (gradient, flow, resistance, inertance)) or min(resistance, inertance) <= 0:
        raise ModelDomainError("Invalid inertial relation")
    return (gradient - resistance * flow) / inertance
