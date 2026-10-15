"""Seeded Morris trajectories and whole-trajectory bootstrap diagnostics.

Inputs are normalized to [0,1]. Elementary effects are divided by their signed
normalized step, so their units are output units per full parameter range.
This is a grid-based global screening measure, not a Sobol variance fraction.
"""
import numpy as np

def trajectories(count,dimension,seed,levels=4):
    if count<1 or dimension<1 or levels<4 or levels%2:raise ValueError('Positive dimensions and even levels >=4 required')
    rng=np.random.default_rng(seed);delta=levels/(2*(levels-1));out=[]
    for _ in range(count):
        base=rng.integers(0,levels//2,size=dimension)/(levels-1)
        direction=rng.choice([-1,1],dimension);x=base+np.where(direction<0,delta,0)
        rows=[x.copy()]
        for j in rng.permutation(dimension):
            x=x.copy();x[j]+=direction[j]*delta;rows.append(x)
        out.append(np.array(rows))
    return np.array(out)

def elementary_effects(x,y):
    x=np.asarray(x,float);y=np.asarray(y,float)
    if x.ndim!=3 or y.shape!=x.shape[:2] or x.shape[1]!=x.shape[2]+1 or not np.isfinite(x).all() or not np.isfinite(y).all():raise ValueError('Complete finite trajectory arrays required')
    if x.min() < -1e-10 or x.max()>1+1e-10:raise ValueError('Normalized inputs outside bounds')
    effects=[]
    for xx,yy in zip(x,y):
        row=np.full(x.shape[2],np.nan)
        for dx,dy in zip(np.diff(xx,axis=0),np.diff(yy)):
            changed=np.flatnonzero(abs(dx)>1e-10)
            if len(changed)!=1 or np.isfinite(row[changed[0]]):raise ValueError('Each coordinate must change exactly once')
            j=changed[0];row[j]=dy/dx[j]
        if not np.isfinite(row).all():raise ValueError('Missing coordinate')
        effects.append(row)
    return np.array(effects)

def summarize(effects,seed=7203,resamples=5000):
    e=np.asarray(effects,float)
    if e.ndim!=2 or len(e)<2 or not np.isfinite(e).all():raise ValueError('At least two finite trajectories required')
    rng=np.random.default_rng(seed);idx=rng.integers(len(e),size=(resamples,len(e)))
    boot=np.abs(e[idx]).mean(axis=1);mu=np.abs(e).mean(axis=0)
    # Whole trajectories resample related effects together. Rank intervals are
    # conditional finite-sample diagnostics, not probabilities of physiology.
    order=np.argsort(-boot,axis=1);ranks=np.empty_like(order)
    np.put_along_axis(ranks,order,np.arange(1,e.shape[1]+1)[None,:],axis=1)
    low,high=np.quantile(boot,[.025,.975],axis=0)
    ranklow,rankhigh=np.quantile(ranks,[.025,.975],axis=0)
    return dict(mu=e.mean(axis=0).tolist(),mu_star=mu.tolist(),sigma=e.std(axis=0,ddof=1).tolist(),mu_star_interval=[low.tolist(),high.tolist()],rank_interval=[ranklow.tolist(),rankhigh.tolist()],top_two_frequency=(ranks<=2).mean(axis=0).tolist(),resamples=resamples,seed=seed)
