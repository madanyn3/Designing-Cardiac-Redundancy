"""Exact contrasts on an equally weighted, complete two-level factorial.

Variance shares describe this discrete design measure only, not a population
or a continuous-domain Sobol index. Numerical pairs are not replicates.
"""
import itertools
import numpy as np

def decompose(signs,values,names):
    x=np.asarray(signs,dtype=float);y=np.asarray(values,dtype=float);d=len(names)
    expected=set(itertools.product((-1.,1.),repeat=d))
    if x.shape!=(2**d,d) or set(map(tuple,x))!=expected or y.shape!=(2**d,) or not np.isfinite(y).all():
        raise ValueError('Complete unique finite two-level design required')
    terms={};reconstruction=np.full(len(y),y.mean())
    for count in range(1,d+1):
        for idx in itertools.combinations(range(d),count):
            basis=np.prod(x[:,idx],axis=1);coefficient=float(np.mean(y*basis))
            key=' * '.join(names[i] for i in idx)
            terms[key]=dict(coefficient=coefficient,effect=2*coefficient,order=count)
            reconstruction+=coefficient*basis
    variance=float(np.mean((y-y.mean())**2))
    for term in terms.values():term['variance_share']=term['coefficient']**2/variance if variance>1e-24 else None
    totals={name:sum(t['coefficient']**2 for key,t in terms.items() if name in key.split(' * '))/variance if variance>1e-24 else None for name in names}
    return dict(mean=float(y.mean()),variance=variance,terms=terms,total_shares=totals,reconstruction_error=float(np.max(abs(y-reconstruction))))

def interpolate(decomposition,values,bounds):
    coded={name:(2*values[name]-lo-hi)/(hi-lo) for name,(lo,hi) in bounds.items()}
    return decomposition['mean']+sum(t['coefficient']*np.prod([coded[n] for n in key.split(' * ')]) for key,t in decomposition['terms'].items())
