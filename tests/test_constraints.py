import numpy as np
import pytest
from cardiac_redundancy.constraints import pressure_windows, oxygen_screen, positive_exposure, combined_status


def test_pressure_uses_integral_and_joint_screens_not_map_formula():
    r = pressure_windows([0,.2,1], [50,100,50], dict(mean=70,maximum=90,minimum=50))
    assert r['windows'][0]['metrics_mmHg']['mean'] == pytest.approx(75)
    assert r['windows'][0]['window_pressure_range_mmHg'] == 50
    assert r['status'] == 'satisfied'
    r = pressure_windows([0,1], [80,80], dict(mean=65,maximum=90,minimum=50))
    assert r['status']=='violated'  # adequate mean alone cannot pass


def test_incomplete_window_is_not_silently_passed():
    r = pressure_windows([0,1,1.5], [100]*3, dict(mean=65,maximum=90,minimum=50))
    assert r['status']=='unassessed' and r['incomplete_tail_s']==[1,1.5]


def test_exact_exposure_triangle_and_first_violation():
    assert positive_exposure([0,2], [-1,1]) == dict(duration_s=1.,integral=.5)
    r=oxygen_screen([0,1,2],[1]*3,[1,.5,0],[2,1,0],[0,.25,1],60,.01)
    assert r['first_violation']['bracket_s']==[0,1]
    assert r['extraction_undefined_samples']==1
    assert r['baseline_demand_equivalent_min']==pytest.approx(1/60)


def test_unknown_recovery_blocks_overall_success():
    assert combined_status('valid','satisfied',['satisfied'],['recovery'])['status']=='indeterminate'
    r=combined_status('valid','violated',['satisfied'],['recovery'])
    assert r['status']=='violated_research_screen' and r['unresolved_requirements']==['recovery']
    assert combined_status('numerical_failure','violated',[],['recovery'])['status']=='numerical_failure'
    assert combined_status('model_domain_limit','satisfied',[],[])['status']=='model_domain_limit'


def test_zero_delivery_and_initial_failure():
    r=oxygen_screen([0,1],[1,1],[0,0],[0,0],[0,1],60,1e-6)
    assert r['extraction_fraction_range'] is None and r['first_violation']['initial']
    with pytest.raises(ValueError): pressure_windows([0,0],[80,80],dict(mean=65,maximum=90,minimum=50))
    with pytest.raises(ValueError): oxygen_screen([0,1],[1,1],[2,2],[2,2],[0,0],60,0)


def test_tissue_storage_separates_consumption_from_blood_extraction():
    r=oxygen_screen([0,1],[1,1],[1,1],[0,0],[0,0],60,1e-8,blood_extraction=[0,0])
    assert r['status']=='satisfied' and r['extraction_fraction_range'] is None
