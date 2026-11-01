"""Official two-view MOFA+ and honest test-participant prediction."""
import contextlib,io,time,hashlib
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score,balanced_accuracy_score,brier_score_loss


def prep_fit(views,c):
    # Ancestry and batch adjustment, no subtype, age or cellular composition removal.
    x=np.c_[np.ones(len(c)),c[:,2:]];out=[];params=[]
    for a in views:
        coef=np.column_stack([np.linalg.lstsq(x[np.isfinite(a[:,j])],a[np.isfinite(a[:,j]),j],rcond=None)[0] for j in range(a.shape[1])])
        z=a-x@coef;scale=np.maximum(np.nanstd(z,axis=0),.1)
        out.append(z/scale);params.append((coef,scale))
    return out,params


def prep_apply(views,c,params):
    x=np.c_[np.ones(len(c)),c[:,2:]]
    return [(a-x@b)/s for a,(b,s) in zip(views,params)]


def mofa(views,k,seed,iterations,perturb=False):
    from mofapy2.run.entry_point import entry_point
    log=io.StringIO()
    with contextlib.redirect_stdout(log),contextlib.redirect_stderr(log):
        e=entry_point();e.set_data_options(scale_views=False,center_groups=False)
        e.set_data_matrix([[v.copy()] for v in views],likelihoods=['gaussian','gaussian'],views_names=['genotype','RNA'])
        e.set_model_options(factors=k,ard_weights=True,spikeslab_weights=True,ard_factors=False,spikeslab_factors=False)
        e.set_train_options(iter=iterations,convergence_mode='medium',seed=seed,quiet=True,verbose=False,dropR2=None)
        e.build();node=e.model.nodes['Z']
        if perturb:
            node.Q.params['mean']+=np.random.default_rng(seed+1729).normal(0,.25,size=node.Q.params['mean'].shape);node.updateExpectations()
        ih=hashlib.sha256(node.getExpectation().tobytes()).hexdigest();e.run()
        w=e.model.nodes['W'].getExpectation();z=e.model.nodes['Z'].getExpectation()
        tau=[np.mean(t,axis=0) if np.ndim(t)==2 else np.asarray(t) for t in e.model.nodes['Tau'].getExpectation()]
        elbo=np.asarray(e.model.getTrainingStats()['elbo']);elbo=elbo[np.isfinite(elbo)]
    return dict(weights=w,scores=z,precision=tau,elbo=float(elbo[-1]),converged='Converged!' in log.getvalue(),initialization_hash=ih)


def project(model,views):
    w=np.concatenate(model['weights']);tau=np.concatenate(model['precision']);x=np.concatenate(views,axis=1)
    cov=np.linalg.inv(np.eye(w.shape[1])+(w.T*tau)@w)
    return (x*tau)@w@cov


def tune(d,seed,config):
    raw=[d['genotype'],d['RNA']];rng=np.random.default_rng(seed+79)
    masks=[rng.uniform(size=v.shape)<config['mask_fraction'] for v in raw]
    masked=[np.where(m,np.nan,v) for v,m in zip(raw,masks)]
    a,params=prep_fit(masked,d['covariates']);target=prep_apply(raw,d['covariates'],params);candidates=[]
    def starts(v,k):
        models=[mofa(v,k,seed+s,config['mofa_iterations'],perturb=(i==1)) for i,s in enumerate(config['mofa_starts'])]
        assert len({m['initialization_hash'] for m in models})==2
        return max(models,key=lambda m:m['elbo']),[dict(hash=m['initialization_hash'],elbo=m['elbo'],converged=m['converged']) for m in models]
    for k in config['factor_grid']:
        model,initial=starts(a,k)
        errors=[float(np.mean((target[j][masks[j]]-(model['scores']@model['weights'][j].T)[masks[j]])**2)) for j in range(2)]
        candidates.append(dict(k=k,validation_mse=float(np.mean(errors)),view_errors=errors,initializations=initial))
    k=min(candidates,key=lambda r:(r['validation_mse'],r['k']))['k']
    a,params=prep_fit(raw,d['covariates']);model,initial=starts(a,k)
    return model,params,k,candidates,initial


def classifier(a,b,train,test):
    scaler=StandardScaler().fit(a);fit=LogisticRegression(C=1,max_iter=1000).fit(scaler.transform(a),train['y'])
    p=fit.predict_proba(scaler.transform(b))[:,1]
    return dict(auc=float(roc_auc_score(test['y'],p)),balanced_accuracy=float(balanced_accuracy_score(test['y'],p>=.5)),brier=float(brier_score_loss(test['y'],p))),p


def run(train,test,seed,config):
    started=time.perf_counter();model,params,k,candidates,initial=tune(train,seed,config)
    a=prep_apply([train['genotype'],train['RNA']],train['covariates'],params);b=prep_apply([test['genotype'],test['RNA']],test['covariates'],params)
    za=project(model,a);zb=project(model,b);records=[];predictions={}
    pca=PCA(n_components=k,svd_solver='full').fit(np.concatenate(a,axis=1))
    methods={
        'MOFA+':(za,zb),
        'PCA':(pca.transform(np.concatenate(a,axis=1)),pca.transform(np.concatenate(b,axis=1))),
        'Genotype ridge':(a[0],b[0]),
        'RNA ridge':(a[1],b[1]),
        'Combined ridge':(np.concatenate(a,axis=1),np.concatenate(b,axis=1)),
        'Clinical only':(np.empty((len(train['y']),0)),np.empty((len(test['y']),0)))
    }
    for name,(x,xt) in methods.items():
        metrics,p=classifier(np.c_[train['covariates'],x],np.c_[test['covariates'],xt],train,test)
        records.append(dict(method=name,status='ok',**metrics,k=k if name in ['MOFA+','PCA'] else None,diagnostic_pass=bool(model['converged']) if name=='MOFA+' else True));predictions[name]=p.tolist()
    return dict(records=records,predictions=predictions,candidates=candidates,selected_k=k,final_initializations=initial,seconds=time.perf_counter()-started),dict(weights_genotype=model['weights'][0],weights_RNA=model['weights'][1],scores=za,precision_genotype=model['precision'][0],precision_RNA=model['precision'][1],test_scores=zb)
