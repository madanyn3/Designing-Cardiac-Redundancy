"""Passive blood-marker transport on a saved hydraulic network.

This tracks origin/washout of blood, not oxygen chemistry, viable RBC counts,
anticoagulant effectiveness or clot formation. H is a zero-storage mixing
junction; its concentration is solved from its instantaneous inflows.
"""
import numpy as np


def marker_rhs(volumes,amounts,edges,flows,names):
    volumes=np.asarray(volumes);amounts=np.asarray(amounts)
    if np.any(volumes<=0):raise ValueError('Positive transport volumes required')
    concentration=dict(zip(names,amounts/volumes))
    incoming=0.;weighted=0.;net=0.
    for (_,a,b),q in zip(edges,flows):
        if a!='H' and b!='H':continue
        into=q if b=='H' else -q;net+=into
        if into>0:
            donor=a if b=='H' else b
            incoming+=into;weighted+=into*concentration[donor]
    if abs(net)>1e-7:raise ValueError('Unbalanced zero-storage junction')
    concentration['H']=weighted/incoming if incoming else 0.
    change=dict.fromkeys(names,0.)
    for (_,a,b),q in zip(edges,flows):
        donor=a if q>=0 else b
        transfer=q*concentration[donor]
        if a!='H':change[a]-=transfer
        if b!='H':change[b]+=transfer
    return np.array([change[n] for n in names])
