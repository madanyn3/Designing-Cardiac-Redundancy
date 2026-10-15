import json
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from cardiac_redundancy.inventory import PhysicalInventory, PhysicalSourceReference
from cardiac_redundancy.integration import integrate_source
from cardiac_redundancy.mechanics import ModelDomainError

ROOT=Path(__file__).resolve().parents[1]


def test_presplit_physical_ode_is_source_coordinate_translation():
    model=PhysicalSourceReference(ROOT/'configs/baseline.yaml',ROOT/'configs/project-extension.yaml')
    assert model.initial_state[:6].sum()==pytest.approx(5250)
    assert model.translation[3]==3000
    native=integrate_source(model.source,.75)
    physical=integrate_source(model,.75)
    np.testing.assert_allclose(physical.states-model.translation[:,None],native.states,rtol=0,atol=2e-4)
    np.testing.assert_allclose(physical.states[:6].sum(axis=0),5250,rtol=0,atol=1e-6)
    assert [(e['valve'],e['action']) for e in physical.events]==[(e['valve'],e['action']) for e in native.events]
    assert max(abs(a['time_s']-b['time_s']) for a,b in zip(physical.events,native.events))<1e-5


@pytest.fixture
def inventory():
    return PhysicalInventory(ROOT/"configs/baseline.yaml",ROOT/"configs/project-extension.yaml")


def test_frozen_physical_partition_and_pressure(inventory):
    v=inventory.initial_volumes
    u=inventory.unstressed_volumes
    assert v.sum()==pytest.approx(5250)
    assert u.sum()==pytest.approx(3750)
    assert (v-u)[:6].sum()==pytest.approx(1350)
    assert v[6:].sum()==pytest.approx(600)
    assert (v-u)[6:].sum()==pytest.approx(150)
    pressure=inventory.vascular_pressures(v)
    assert pressure['VC']==pytest.approx(179.7803*.0059)
    for name in inventory.names[6:]:
        assert pressure[name]==pytest.approx(5)
    # Independent source pressure comparison: all central vascular pressures
    # except VC match after offsets; VC gives up 150 stressed mL to limb pools.
    source=json.loads((ROOT/'configs/baseline.yaml').read_text())
    for name,component in [('AO','aorta'),('PA','pulmonary_artery'),('PU','pulmonary_vein')]:
        expected=source['initial_volumes_mL'][name]*source['source_parameters'][component]['E_es_'+name.lower()]+source['pressure_external_mmHg'][name]
        assert pressure[name]==pytest.approx(expected)
    v[:]=0
    assert inventory.initial_volumes.sum()==pytest.approx(5250)


def test_transfers_preserve_physical_blood_with_reverse_and_closed_edges(inventory):
    edges=[('AO','B_limb_muscle'),('B_limb_muscle','VC'),('PU','LV')]
    a=inventory.incidence(edges)
    np.testing.assert_array_equal(a.sum(axis=0),0)
    rates=a@np.array([0.,-2.3,7.1])
    assert rates.sum()==pytest.approx(0,abs=1e-14)
    assert rates[inventory.names.index('B_limb_muscle')]==pytest.approx(2.3)


def test_physical_reservoir_drainage_matches_analytic_solution(inventory):
    # Two limb reservoirs exchange through a resistor, with all other routes
    # closed. Analytic exponential pressure equilibration checks U handling,
    # conservation and signed edge assembly independently of heart mechanics.
    donor,recipient='B_limb_muscle','B_limb_skin'
    i,j=map(inventory.names.index,(donor,recipient))
    c1,c2=[inventory.reservoirs[k].compliance_mL_mmHg for k in (i,j)]
    r=2.
    v=inventory.initial_volumes
    v[i]+=10; v[j]-=10
    a=inventory.incidence([(donor,recipient)])
    def rhs(t,y):
        p=inventory.vascular_pressures(y)
        return a[:,0]*(p[donor]-p[recipient])/r
    sol=solve_ivp(rhs,(0,20),v,rtol=1e-10,atol=1e-10,dense_output=True)
    assert sol.success
    times=np.linspace(0,20,101)
    ys=sol.sol(times)
    p0=inventory.vascular_pressures(v)
    expected=(p0[donor]-p0[recipient])*np.exp(-times*(1/c1+1/c2)/r)
    observed=np.array([inventory.vascular_pressures(y)[donor]-inventory.vascular_pressures(y)[recipient] for y in ys.T])
    np.testing.assert_allclose(observed,expected,atol=1e-8,rtol=1e-8)
    np.testing.assert_allclose(ys.sum(axis=0),5250,atol=1e-9,rtol=0)
    for y in ys.T: inventory.validate(y)


def test_domain_is_physical_and_does_not_clip(inventory):
    v=inventory.initial_volumes
    vc=inventory.names.index('VC')
    v[vc]=inventory.unstressed_volumes[vc]
    inventory.validate(v)  # positive physical V at zero transmural pressure
    v[vc]-=.1
    assert inventory.domain_margins(v)['VC']<0
    assert inventory.vascular_pressures(v)['VC']<0
    before=v.copy()
    with pytest.raises(ModelDomainError): inventory.validate(v)
    np.testing.assert_array_equal(v,before)


def test_inconsistent_inventory_is_rejected(tmp_path):
    c=json.loads((ROOT/'configs/project-extension.yaml').read_text())
    c['initial_physical_volumes_mL']['VC']+=1
    p=tmp_path/'bad.json'; p.write_text(json.dumps(c))
    with pytest.raises(ValueError,match='sum to target'):
        PhysicalInventory(ROOT/'configs/baseline.yaml',p)
