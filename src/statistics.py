"""Association, regional testing and restricted HC2 wild bootstrap."""
import numpy as np
from scipy.special import expit
from scipy.stats import chi2, f, hypergeom
from scipy.optimize import minimize


def bh(p):
    p=np.asarray(p);ix=np.argsort(p);out=np.empty(len(p))
    out[ix]=np.minimum(1,np.minimum.accumulate((p[ix]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1])
    return out


def logistic(x,y):
    def objective(b):
        eta=x@b
        return np.logaddexp(0,eta).sum()-y@eta,x.T@(expit(eta)-y)
    fit=minimize(objective,np.zeros(x.shape[1]),jac=True,method='BFGS',options={'gtol':1e-7,'maxiter':500})
    assert np.max(np.abs(objective(fit.x)[1]))<1e-3,'Logistic fit gradient too large'
    return float(-fit.fun),fit.x


def factor_associations(z,c,y):
    z=(z-z.mean(0))/np.maximum(z.std(0),1e-8);x=np.c_[np.ones(len(y)),c]
    null,_=logistic(x,y);p=[];effect=[]
    for j in range(z.shape[1]):
        value,b=logistic(np.c_[x,z[:,j]],y)
        p.append(float(chi2.sf(max(0,2*(value-null)),1)));effect.append(float(b[-1]))
    return dict(p=p,q=bh(p).tolist(),effect=effect)


def eqtl(d,draws=0,seed=0,fixture=False):
    n=len(d['y']);base=np.c_[np.ones(n),d['covariates'],d['y']]
    p=np.empty((24,2));pw=np.empty((24,2)) if draws else None
    rng=np.random.default_rng(seed);multipliers=rng.choice([-1.,1.],size=(n,draws)) if draws else None
    saved=None
    for locus in range(12):
        g=d['genotype'][:,4*locus:4*locus+4];x=np.c_[base,g,g*d['y'][:,None]]
        assert np.linalg.matrix_rank(x)==x.shape[1]
        q,_=np.linalg.qr(x,mode='reduced');y=d['RNA'][:,2*locus:2*locus+2]
        theta=q.T@y;rss=np.sum((y-q@theta)**2,axis=0);df=n-x.shape[1]
        for family,p0 in enumerate([base.shape[1],base.shape[1]+4]):
            obs=np.sum(theta[p0:]**2,axis=0)/(x.shape[1]-p0)/(rss/df)
            p[2*locus:2*locus+2,family]=f.sf(obs,x.shape[1]-p0,df)
            if draws:
                q0=q[:,:p0];r=y-q0@(q0.T@y);adjusted=r/np.sqrt(1-np.sum(q0*q0,axis=1))[:,None]
                for j in range(2):
                    projection=q.T@(adjusted[:,j,None]*multipliers)
                    den=np.maximum(np.sum(adjusted[:,j]**2)-np.sum(projection**2,axis=0),1e-12)/df
                    stars=np.sum(projection[p0:]**2,axis=0)/(x.shape[1]-p0)/den
                    pw[2*locus+j,family]=(1+np.sum(stars>=obs[j]))/(draws+1)
                    if fixture and locus==0 and family==0 and j==0:
                        saved=dict(x=x.tolist(),y=y[:,j].tolist(),p0=p0,weights=multipliers[:,:17].tolist(),F=float(obs[j]),stars=stars[:17].tolist())
    result=dict(p=p.tolist(),q=np.column_stack([bh(p[:,i]) for i in range(2)]).tolist())
    if draws:result.update(wild_p=pw.tolist(),wild_q=np.column_stack([bh(pw[:,i]) for i in range(2)]).tolist())
    if fixture:result['fixture']=saved
    return result


def levels(genes):
    genes=set(map(int,genes))
    return dict(gene=genes,locus={j//2 for j in genes},module={j//6 for j in genes})


def overlap(a,b):
    a=set(a);b=set(b)
    return dict(intersection=len(a&b),union=len(a|b),jaccard=len(a&b)/len(a|b) if a|b else None)


def module_enrichment(genes):
    genes=set(map(int,genes));p=[float(hypergeom.sf(len(genes & set(range(6*m,6*m+6)))-1,24,6,len(genes))) for m in range(4)]
    return dict(p=p,q=bh(p).tolist(),modules=np.flatnonzero(bh(p)<.05).tolist())
