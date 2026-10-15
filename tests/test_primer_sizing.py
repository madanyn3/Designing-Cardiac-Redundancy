from dataclasses import replace
import math
import pytest
from cardiac_redundancy.primer_sizing import PrimerInputs,coefficients,geometry,fixed_point,iterate,steady_power_W

@pytest.mark.parametrize('exact',[False,True])
def test_independent_volume_and_feedback_identities(exact):
    p=PrimerInputs();r=fixed_point(p,exact);a=r['radius_cm'];h=r['wall_cm']
    volume=2*4*math.pi*((a+h)**3-a**3)/3 if exact else 2*4*math.pi*a*a*h
    assert r['two_shell_mass_g']==pytest.approx(volume*p.density_g_mL)
    assert r['flow_mL_min']==pytest.approx(p.external_mL_min+p.self_flow_mL_min_g*r['two_shell_mass_g'])

@pytest.mark.parametrize('seed',[0,700,1500])
@pytest.mark.parametrize('exact',[False,True])
def test_iteration_reaches_analytic_solution(seed,exact):
    p=PrimerInputs();r=iterate(p,seed,1e-8,exact)
    assert abs(r['values'][-1]-fixed_point(p,exact)['flow_mL_min'])<=1e-8
    assert r['error_bound_mL_min']<=1e-8

def test_no_positive_solution_at_or_above_unit_gain():
    p=PrimerInputs();k=coefficients(p)[2]
    for gain in [1,1.1]:
        bad=replace(p,self_flow_mL_min_g=gain/k)
        with pytest.raises(ValueError):fixed_point(bad)
        with pytest.raises(ValueError):iterate(bad,700)

def test_si_power_and_series_heads():
    assert steady_power_W(1000,1)==pytest.approx(133.322387415e-3/60)
    assert steady_power_W(740,80)==pytest.approx(steady_power_W(740,70)+steady_power_W(740,10))

def test_invalid_inputs():
    with pytest.raises(ValueError):PrimerInputs(ejection_fraction=1.1)
    with pytest.raises(ValueError):PrimerInputs(stress_kPa=0)
    with pytest.raises(ValueError):iterate(PrimerInputs(),float('nan'))
