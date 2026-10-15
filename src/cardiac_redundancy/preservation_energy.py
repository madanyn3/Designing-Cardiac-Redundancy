"""Declared native-heart demand transitions and non-credit energy reporting.

A decline in prescribed demand is an assumption, not inferred from pressure or
lack of contraction. The oxygen-energy equivalent is a reporting convention,
not ATP yield, efficiency calibration, or an added source of energy.
"""
import math
import numpy as np

def native_fraction(time_s,target,ramp_s):
    if not all(math.isfinite(v) for v in (time_s,target,ramp_s)) or time_s<0 or not 0<target<=1 or ramp_s<0:raise ValueError('Invalid preservation transition')
    if ramp_s==0:return target
    return 1-(1-target)*min(time_s/ramp_s,1.)

def integrated_native_fraction(time_s,target,ramp_s):
    native_fraction(time_s,target,ramp_s)
    if ramp_s==0:return time_s*target
    t=min(time_s,ramp_s)
    return t-(1-target)*t*t/(2*ramp_s)+max(0,time_s-ramp_s)*target

def oxygen_energy(oxygen_mL,kcal_per_L=4.85):
    a=np.asarray(oxygen_mL,dtype=float)
    if np.any(~np.isfinite(a)) or np.any(a< -1e-10) or not math.isfinite(kcal_per_L) or kcal_per_L<=0:raise ValueError('Invalid oxygen-energy reporting input')
    return a*kcal_per_L*4.184
