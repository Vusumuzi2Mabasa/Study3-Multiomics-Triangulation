"""Wholly synthetic paired genotype and RNA; truth never enters model fitting."""
import numpy as np
from scipy.special import expit, ndtr

SCENARIOS=['null','regulation_only','abundance_only','convergent','interaction_only']


def generate(n,seed,scenario):
    if scenario not in SCENARIOS: raise ValueError(scenario)
    rng=np.random.default_rng(seed)
    ancestry=rng.choice(3,n,p=[.5,.3,.2]);batch=rng.integers(0,3,n)
    age=rng.normal(size=n);cell=rng.normal(size=n)
    c=np.c_[age,cell,ancestry==1,ancestry==2,batch==1,batch==2].astype(float)
    g=np.empty((n,48));ld=.65**np.abs(np.arange(4)[:,None]-np.arange(4)[None,:])
    for locus in range(12):
        u=ndtr(rng.multivariate_normal(np.zeros(4),ld,size=n))
        p=np.clip(np.array([.18,.25,.33,.4])+(ancestry[:,None]-.7)*.025,.08,.46)
        g[:,4*locus:4*locus+4]=(u>(1-p)**2).astype(int)+(u>1-p*p).astype(int)
    y=rng.binomial(1,expit(-.6+.25*age+.2*cell+.15*(ancestry==2)))
    h=rng.normal(size=(n,4));beta=np.zeros(24);delta=np.zeros(24);gamma=np.zeros(24)
    if scenario in ['regulation_only','convergent']:beta[:6]=.5
    if scenario in ['abundance_only','convergent']:gamma[:6]=.55
    if scenario=='interaction_only':delta[:6]=np.array([1,-1,1,-1,1,-1])*.55
    latent=np.empty((n,24))
    for j in range(24):
        variant=4*(j//2)+(0 if j%2==0 else 2)
        centered=g[:,variant]-2*np.array([.18,.25,.33,.4])[variant%4]
        latent[:,j]=5+.18*age+.18*cell+.1*c[:,2]-.1*c[:,3]+.12*c[:,4]-.12*c[:,5]+.45*h[:,j//6]+beta[j]*centered+gamma[j]*y+delta[j]*centered*y+rng.normal(0,.5,n)
    library=rng.lognormal(np.log(20e6),.2,n);mu=library[:,None]/1e6*2**latent
    counts=rng.negative_binomial(30,30/(30+mu));rna=np.log2(counts/library[:,None]*1e6+1)
    return dict(genotype=g,RNA=rna,counts=counts,library=library,covariates=c,y=y,h=h,beta=beta,delta=delta,gamma=gamma,seed=seed,scenario=scenario)


def qc(d):
    n=len(d['y'])
    assert d['genotype'].shape==(n,48) and d['RNA'].shape==(n,24)
    assert np.isin(d['genotype'],[0,1,2]).all()
    assert np.isfinite(d['RNA']).all() and np.isfinite(d['covariates']).all()
    assert (d['counts']>=0).all() and np.all(d['counts']==d['counts'].astype(int))
    assert np.all(d['library']>0) and len(np.unique(d['y']))==2
    assert np.min(np.bincount(d['y']))>=20
    return dict(status='passed',n=n,TNBC=int(d['y'].sum()),minimum_MAF=float(np.minimum(d['genotype'].mean(0)/2,1-d['genotype'].mean(0)/2).min()))
