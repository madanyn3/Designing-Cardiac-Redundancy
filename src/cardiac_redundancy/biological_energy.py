"""BIO hydraulic energy accounting, independent of metabolic calibration.

Pressures are gauge mmHg, flows mL/s. Signed boundary work is retained.
Pericardial exchange is separate from transmural cavity work; neither is PVA
or a uniquely assigned septal muscle energy. No oxygen credit is inferred.
"""
import numpy as np
from .pump_energy import MMHG_ML_TO_J

def chamber_power(pin,pc,pout,peri,qin,qout,rin,rout,lin,lout):
    vals=np.array([pin,pc,pout,peri,qin,qout,rin,rout,lin,lout],dtype=float)
    if not np.all(np.isfinite(vals)) or min(rin,rout,lin,lout)<=0 or min(qin,qout)<-1e-7:
        raise ValueError('Finite pressures, positive R/L and forward valve flow required')
    k=MMHG_ML_TO_J
    # Multiply the inertial valve equations by Q; at a closed valve Q=0,
    # so its contribution is zero regardless of the pressure difference.
    inlet_inertia=qin*(pin-pc-rin*qin)
    outlet_inertia=qout*(pc-pout-rout*qout)
    chamber=k*pc*(qout-qin)
    external=k*peri*(qout-qin)
    port=k*(pout*qout-pin*qin)
    loss=k*(rin*qin**2+rout*qout**2)
    inertia=k*(inlet_inertia+outlet_inertia)
    return dict(chamber_W=chamber,transmural_W=chamber-external,
        pericardial_exchange_W=external,port_W=port,valve_loss_W=loss,
        valve_storage_rate_W=inertia,valve_stored_J=.5*k*(lin*qin**2+lout*qout**2),
        balance_residual_W=chamber-port-loss-inertia)

def biological_power(model,time,state):
    if model.architecture!='BIO':raise ValueError('BIO model required')
    p,_=model.pressures(time,state,active=False);n=model.volume_count
    result={}
    for side,node,inlet,outlet,j in (('left','aLV','PU','AO',4),('right','aRV','VC','PA',6)):
        result[side]=chamber_power(p[inlet],p[node],p[outlet],p['aPERI'],
            state[n+j],state[n+j+1],model.valve_R[j],model.valve_R[j+1],model.valve_L[j],model.valve_L[j+1])
    return result
