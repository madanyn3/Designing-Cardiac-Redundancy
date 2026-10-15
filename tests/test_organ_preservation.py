import numpy as np
import pytest
from cardiac_redundancy.organ_oxygen import OrganDemand
from cardiac_redundancy.organ_preservation import OrganPreservation

def fixture():
    organs=[OrganDemand('brain','brain',48,0.75,'brain'),
            OrganDemand('kidney','kidney',20,0.75,'kidneys'),
            OrganDemand('aux','aux',5,0.75,'auxiliary_heart')]
    entries={'brain':dict(role='active',fraction=.6,ramp_s=5),
             'kidney':dict(role='preservation',fraction=.25,ramp_s=0),
             'aux':dict(role='pump',fraction=1,ramp_s=0)}
    return organs,entries

def test_independent_targets_and_transition_accounting():
    organs,entries=fixture();p=OrganPreservation(organs,entries)
    np.testing.assert_allclose(p.demand(0),[48,5,5])
    np.testing.assert_allclose(p.demand(8),[28.8,5,5])
    t=np.linspace(0,8,8001)
    np.testing.assert_allclose(np.trapezoid([p.demand(x)/60 for x in t],t,axis=0),p.integrated(8))
    assert organs[1].demand_mL_min==20
    assert p.groups()=={'active':[0],'preservation':[1],'pump':[2]}

@pytest.mark.parametrize('change',['missing','unknown','zero','negative','nan','pump_discount','wrong_role'])
def test_invalid_profiles_fail_closed(change):
    organs,e=fixture()
    if change=='missing':e.pop('kidney')
    elif change=='unknown':e['typo']=dict(e['kidney'])
    elif change=='zero':e['kidney']['fraction']=0
    elif change=='negative':e['kidney']['ramp_s']=-1
    elif change=='nan':e['brain']['fraction']=float('nan')
    elif change=='pump_discount':e['aux']['fraction']=.1
    else:e['aux']['role']='preservation'
    with pytest.raises(ValueError):OrganPreservation(organs,e)
