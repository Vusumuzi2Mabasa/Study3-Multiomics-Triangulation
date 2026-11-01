"""Documented convergence-budget extension; original 150-iteration fits retained."""
from pathlib import Path
import json,time,hashlib
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from scipy.optimize import linear_sum_assignment
from sklearn.decomposition import PCA
from simulation import generate
from statistics import factor_associations
from mofa_engine import tune,prep_apply,prep_fit,project,classifier,mofa
from run import nominations,save,lock

ROOT=Path(__file__).resolve().parents[1]


def refine(path):
    out=ROOT/'results/refined'/path.name
    if out.exists():return
    row=json.loads(path.read_text());config=json.loads((ROOT/'protocol/design.json').read_text())
    settings=dict(config['mofa'],mofa_iterations=1000)
    seed=row['seed'];n=row['n'];scenario=row['scenario'];began=time.perf_counter()
    train=generate(n,seed,scenario);test=generate(config['test_size'],row['test_seed'],scenario)
    model,params,k,candidates,initial=tune(train,seed,settings)
    a=prep_apply([train['genotype'],train['RNA']],train['covariates'],params);b=prep_apply([test['genotype'],test['RNA']],test['covariates'],params)
    za=project(model,a);zb=project(model,b)
    assoc=factor_associations(za,train['covariates'],train['y']);validation=factor_associations(zb,test['covariates'],test['y'])
    row.update(selected_k=k,candidates=candidates,initializations=initial,mofa_converged=bool(model['converged']),mofa_association=assoc,mofa_validation=validation,mofa_genes=nominations(model,assoc),mofa_confirmed_genes=nominations(model,dict(q=np.maximum(assoc['q'],validation['q']).tolist())))
    pca=PCA(n_components=k,svd_solver='full').fit(np.concatenate(a,axis=1))
    for name,(x,xt) in {'MOFA factors':(za,zb),'PCA scores':(pca.transform(np.concatenate(a,axis=1)),pca.transform(np.concatenate(b,axis=1)))}.items():
        metrics,p=classifier(np.c_[train['covariates'],x],np.c_[test['covariates'],xt],train,test)
        row['prediction'][name]=dict(**metrics,probability=p.tolist())
    row['reconstruction_r2']=[float(1-np.sum((a[v]-za@model['weights'][v].T)**2)/np.sum((a[v]-a[v].mean(0))**2)) for v in range(2)]
    row['stability']=[]
    if row['replicate']==0:
        rng=np.random.default_rng(seed+69000)
        for boot in range(config['stability_bootstraps_first_replicate']):
            idx=rng.integers(0,n,n);aa,_=prep_fit([train['genotype'][idx],train['RNA'][idx]],train['covariates'][idx])
            bm=mofa(aa,k,seed+80000+boot,1000)
            w=np.concatenate(model['weights']);wb=np.concatenate(bm['weights'])
            sim=np.abs(w.T@wb)/(np.linalg.norm(w,axis=0)[:,None]*np.linalg.norm(wb,axis=0)[None,:]+1e-12);ii,jj=linear_sum_assignment(-sim)
            ba=factor_associations(project(bm,aa),train['covariates'][idx],train['y'][idx])
            row['stability'].append(dict(matched_weight_cosine=sim[ii,jj].tolist(),genes=nominations(bm,ba),converged=bool(bm['converged'])))
    folder=ROOT/'results/refined_models';folder.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(folder/(path.stem+'.npz'),weights_G=model['weights'][0],weights_RNA=model['weights'][1],scores=za,test_scores=zb,train_y=train['y'],test_y=test['y'],train_cov=train['covariates'],test_cov=test['covariates'])
    row['iteration_cap']=1000;row['refinement_seconds']=time.perf_counter()-began
    save(out,row)


def main():
    lock();path=ROOT/'protocol/convergence_amendment.json'
    value=dict(reason='Only 52 of 200 final MOFA fits met the medium ELBO stopping criterion at the initial 150-iteration cap (7/100 at N240, 45/100 at N480). Extend the cap for every cohort, including masked tuning candidates, starts, final fits and fixed-fixture bootstrap fits. Retain every original result. No effect size, seed, threshold, factor grid, start or DIABLO/eQTL setting changes.',original_cap=150,new_cap=1000,primary_reporting='Use the uniform extended-budget MOFA results and retain original-budget sensitivity results. Report all remaining convergence failures.',source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    if path.exists():assert json.loads(path.read_text())==value
    else:save(path,value)
    paths=sorted((ROOT/'results/cohorts').glob('*.json'))
    assert len(paths)==200
    with ProcessPoolExecutor(max_workers=2) as pool:
        fs=[pool.submit(refine,p) for p in paths]
        for i,f in enumerate(as_completed(fs),1):f.result();print('Extended MOFA',i,'/',len(paths),flush=True)


if __name__=='__main__':main()
