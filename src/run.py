"""Frozen, resumable repeated cohorts; training never sees test outcomes."""
from pathlib import Path
import argparse,json,hashlib,os,time
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from scipy.optimize import linear_sum_assignment
from simulation import generate,qc,SCENARIOS
from statistics import eqtl,factor_associations
from mofa_engine import tune,prep_apply,project,classifier,mofa
from sklearn.decomposition import PCA

ROOT=Path(__file__).resolve().parents[1]


def save(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_suffix('.tmp');temporary.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n');temporary.replace(path)


def lock():
    names=['protocol/design.json','src/simulation.py','src/statistics.py','src/run.py','src/mofa_engine.py','src/diablo.R']
    value={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names}
    path=ROOT/'protocol/source_lock.json'
    if path.exists():assert json.loads(path.read_text())==value,'Frozen analysis source changed'
    else:save(path,value)


def nominations(model,association):
    selected=np.flatnonzero(np.array(association['q'])<.05)
    genes=set()
    for j in selected:genes.update(np.argsort(-np.abs(model['weights'][1][:,j]),kind='stable')[:5].tolist())
    return sorted(genes)


def main_cohort(task):
    scenario,n,rep=task;out=ROOT/f'results/cohorts/{scenario}_{n}_{rep:02}.json'
    if out.exists():return
    config=json.loads((ROOT/'protocol/design.json').read_text());seed=config['seed_base']+SCENARIOS.index(scenario)*1000000+n*1000+rep
    train=generate(n,seed,scenario);test=generate(config['test_size'],seed+500000000,scenario);began=time.perf_counter()
    check=[qc(train),qc(test)]
    model,params,k,candidates,initial=tune(train,seed,config['mofa'])
    a=prep_apply([train['genotype'],train['RNA']],train['covariates'],params);b=prep_apply([test['genotype'],test['RNA']],test['covariates'],params)
    za=project(model,a);zb=project(model,b)
    assoc=factor_associations(za,train['covariates'],train['y']);validation=factor_associations(zb,test['covariates'],test['y'])
    genes=nominations(model,assoc);confirmed=nominations(model,dict(q=np.maximum(assoc['q'],validation['q']).tolist()))
    pca=PCA(n_components=k,svd_solver='full').fit(np.concatenate(a,axis=1))
    inputs={'MOFA factors':(za,zb),'PCA scores':(pca.transform(np.concatenate(a,axis=1)),pca.transform(np.concatenate(b,axis=1))),'RNA ridge':(a[1],b[1]),'Combined ridge':(np.concatenate(a,axis=1),np.concatenate(b,axis=1)),'Clinical only':(np.empty((n,0)),np.empty((len(test['y']),0)))}
    prediction={}
    for name,(x,xt) in inputs.items():
        metrics,prob=classifier(np.c_[train['covariates'],x],np.c_[test['covariates'],xt],train,test)
        prediction[name]=dict(**metrics,probability=prob.tolist())
    r2=[]
    for view in range(2):
        observed=a[view];reconstructed=za@model['weights'][view].T
        r2.append(float(1-np.sum((observed-reconstructed)**2)/np.sum((observed-observed.mean(0))**2)))
    stability=[]
    if rep==0:
        rng=np.random.default_rng(seed+69000)
        for boot in range(config['stability_bootstraps_first_replicate']):
            idx=rng.integers(0,n,n);resampled={key:(value[idx] if isinstance(value,np.ndarray) and value.shape[0]==n else value) for key,value in train.items()}
            from mofa_engine import prep_fit
            aa,_=prep_fit([resampled['genotype'],resampled['RNA']],resampled['covariates'])
            bootmodel=mofa(aa,k,seed+80000+boot,config['mofa']['mofa_iterations'])
            w=np.concatenate(model['weights']);wb=np.concatenate(bootmodel['weights'])
            similarity=np.abs(w.T@wb)/(np.linalg.norm(w,axis=0)[:,None]*np.linalg.norm(wb,axis=0)[None,:]+1e-12)
            rows,cols=linear_sum_assignment(-similarity)
            assocboot=factor_associations(project(bootmodel,aa),resampled['covariates'],resampled['y'])
            stability.append(dict(matched_weight_cosine=similarity[rows,cols].tolist(),genes=nominations(bootmodel,assocboot),converged=bool(bootmodel['converged'])))
    qtl=eqtl(train,config['wild_draws'],seed+90000,fixture=(rep==0 and scenario=='null' and n==240))
    qtltest=eqtl(test)
    # Inputs are released so the official R implementation is independently reproducible.
    folder=ROOT/f'data/cohorts/{scenario}_{n}_{rep:02}';folder.mkdir(parents=True,exist_ok=True)
    for label,d in [('train',train),('test',test)]:
        for name in ['genotype','RNA','covariates','y']:np.savetxt(folder/f'{label}_{name}.csv',d[name],delimiter=',',fmt='%.12g')
    fixturepath=ROOT/f'results/models/{scenario}_{n}_{rep:02}.npz';fixturepath.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(fixturepath,weights_G=model['weights'][0],weights_RNA=model['weights'][1],scores=za,test_scores=zb,train_y=train['y'],test_y=test['y'],train_cov=train['covariates'],test_cov=test['covariates'])
    value=dict(status='ok',scenario=scenario,n=n,replicate=rep,seed=seed,test_seed=seed+500000000,qc=check,selected_k=k,candidates=candidates,initializations=initial,mofa_converged=bool(model['converged']),mofa_association=assoc,mofa_validation=validation,mofa_genes=genes,mofa_confirmed_genes=confirmed,stability=stability,reconstruction_r2=r2,prediction=prediction,eqtl=qtl,eqtl_validation=qtltest,truth=dict(genetic=(train['beta']!=0)|(train['delta']!=0),subtype=(train['gamma']!=0)|(train['delta']!=0),interaction=train['delta']!=0),test_y=test['y'].tolist(),seconds=time.perf_counter()-began)
    value['truth']={key:v.tolist() for key,v in value['truth'].items()}
    save(out,value)


def calibration(task):
    scenario,n,rep=task;out=ROOT/f'results/calibration/{scenario}_{n}_{rep:04}.json'
    if out.exists():return
    config=json.loads((ROOT/'protocol/design.json').read_text());seed=config['calibration_seed_base']+SCENARIOS.index(scenario)*1000000+n*1000+rep
    d=generate(n,seed,scenario);qc(d);result=eqtl(d,config['wild_draws'],seed+90000)
    save(out,dict(status='ok',scenario=scenario,n=n,replicate=rep,seed=seed,**result))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--stage',choices=['pilot','main','calibration'],required=True);ap.add_argument('--workers',type=int,default=2);args=ap.parse_args()
    config=json.loads((ROOT/'protocol/design.json').read_text())
    if args.stage=='pilot':
        main_cohort(('null',240,0));main_cohort(('convergent',240,0));return
    lock();scenarios=config['scenarios'] if args.stage=='main' else config['calibration_scenarios'];reps=config['replicates'] if args.stage=='main' else config['calibration_replicates']
    tasks=[(s,n,r) for s in scenarios for n in config['sample_sizes'] for r in range(reps)]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        fs=[pool.submit(main_cohort if args.stage=='main' else calibration,t) for t in tasks]
        for i,future in enumerate(as_completed(fs),1):future.result();print(args.stage,i,'/',len(tasks),flush=True)


if __name__=='__main__':main()
