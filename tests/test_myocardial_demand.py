from dataclasses import replace
import pytest
from cardiac_redundancy.myocardial_demand import DemandAssumptions,side_demand

def test_normal_reference_is_identity_not_validation():
    a=DemandAssumptions(mass_ratio=1)
    assert side_demand(20,2,2,80,1,a)['total_mL_min']==pytest.approx(20)

def test_maintenance_survives_no_work_and_no_activation():
    d=side_demand(20,0,2,0,1)
    assert d['total_mL_min']==pytest.approx(.75)
    assert d['activation_mL_min']==d['work_mL_min']==0

def test_independent_work_rate_contractility_and_mass_costs():
    d=side_demand(20,1,2,80,1)
    fast=side_demand(20,1,2,120,1)
    assert fast['activation_mL_min']==pytest.approx(d['activation_mL_min']*1.5)
    assert fast['work_mL_min']==d['work_mL_min']
    strong=side_demand(20,1,2,80,1.5)
    assert strong['activation_mL_min']==pytest.approx(d['activation_mL_min']*1.5)
    heavier=side_demand(20,1,2,80,1,DemandAssumptions(mass_ratio=.5))
    assert heavier['maintenance_mL_min']==pytest.approx(2*d['maintenance_mL_min'])
    assert heavier['work_mL_min']==d['work_mL_min']

@pytest.mark.parametrize('kwargs',[{'mass_ratio':0},{'work_fraction':.5},{'mass_ratio':float('nan')},{'activation_fraction':True}])
def test_invalid_assumptions(kwargs):
    with pytest.raises(ValueError):DemandAssumptions(**kwargs)

def test_no_negative_work_energy_credit():
    with pytest.raises(ValueError):side_demand(20,-1,2,80,1)
