"""Objective-specific recovery and descriptive cross-method triangulation."""
from pathlib import Path
import itertools,json,math
import numpy as np
import pandas as pd
from scipy.stats import t
from statistics import levels,overlap,module_enrichment

ROOT=Path(__file__).resolve().parents[1]
METHODS=['MOFA','DIABLO','eQTL']


def wilson(k,n):
    p=k/n;z=1.95996398454;d=1+z*z/n
    center=(p+z*z/(2*n))/d;half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return float(center-half),float(center+half)


def summarize(frame,groups,value):
    rows=[]
    for keys,cell in frame.groupby(groups,dropna=False,sort=False):
        if not isinstance(keys,tuple):keys=(keys,)
        values=cell[value].dropna().to_numpy(float);n=len(values)
        mean=float(np.mean(values)) if n else None
        se=float(np.std(values,ddof=1)/np.sqrt(n)) if n>1 else None
        delta=float(t.ppf(.975,n-1)*se) if n>1 else None
        rows.append(dict(zip(groups,keys),metric=value,replicates=n,mean=mean,mcse=se,lower=mean-delta if delta is not None else None,upper=mean+delta if delta is not None else None))
    return pd.DataFrame(rows)


def main():
    cfg=json.loads((ROOT/'protocol/design.json').read_text());rows=[];prediction=[];factor=[];tri=[];pair=[];markers=[];stability=[];budgets=[];interactions=[]
    files=sorted((ROOT/'results/final').glob('*.json'));assert len(files)==200
    for path in files:
        d=json.loads(path.read_text());r=json.loads((ROOT/'results/diablo'/path.name).read_text());assert d['status']==r['status']=='ok',path.name
        meta={key:d[key] for key in ['scenario','n','replicate']}
        truth_g=np.array(d['truth']['genetic']);truth_s=np.array(d['truth']['subtype']);truth_i=np.array(d['truth']['interaction'])
        for method,field in [('Ordinary F','q'),('Restricted wild F','wild_q')]:
            selected_i=np.array(d['eqtl'][field])[:,1]<.05
            tp=int(np.sum(selected_i&truth_i));fp=int(np.sum(selected_i&~truth_i))
            interactions.append(dict(**meta,method=method,selected=int(selected_i.sum()),FDP=fp/max(int(selected_i.sum()),1),sensitivity=tp/int(truth_i.sum()) if truth_i.sum() else None))
        sets={'MOFA':set(d['mofa_genes']),'DIABLO':set(np.atleast_1d(r['genes']).astype(int)),'eQTL':set(np.flatnonzero(np.array(d['eqtl']['q'])[:,0]<.05)),'Wild eQTL':set(np.flatnonzero(np.array(d['eqtl']['wild_q'])[:,0]<.05)),'MOFA confirmed':set(d['mofa_confirmed_genes']),'eQTL replicated':set(np.flatnonzero((np.array(d['eqtl']['q'])[:,0]<.05)&(np.array(d['eqtl_validation']['q'])[:,0]<.05)))}
        for method,selected in sets.items():
            truth=truth_g if 'eQTL' in method else truth_s
            true=set(np.flatnonzero(truth));tp=len(selected&true);fp=len(selected-true)
            rows.append(dict(**meta,method=method,selected=len(selected),true_selected=tp,false_selected=fp,FDP=fp/max(1,len(selected)),sensitivity=tp/len(true) if len(true) else None))
            for gene in range(24):markers.append(dict(**meta,method=method,gene=f'SIM_GENE_{gene+1:02}',gene_index=gene,locus=f'SIM_LOCUS_{gene//2+1:02}',module=f'SIM_MODULE_{gene//6+1:02}',selected=gene in selected,true_genetic=bool(truth_g[gene]),true_subtype=bool(truth_s[gene]),true_interaction=bool(truth_i[gene]),omnibus_q=d['eqtl']['q'][gene][0],interaction_q=d['eqtl']['q'][gene][1],wild_omnibus_q=d['eqtl']['wild_q'][gene][0]))
        mapped={method:dict(gene=selected,locus=levels(selected)['locus'],module=set(module_enrichment(selected)['modules'])) for method,selected in sets.items()}
        for level in ['gene','locus','module']:
            for a,b in itertools.combinations(METHODS,2):pair.append(dict(**meta,level=level,comparison=a+' / '+b,**overlap(mapped[a][level],mapped[b][level])))
            for qtl in ['eQTL','Wild eQTL']:
                common=mapped['MOFA'][level]&mapped['DIABLO'][level]&mapped[qtl][level]
                truegenes=set(np.flatnonzero(truth_g&truth_s));true=levels(truegenes)[level]
                tri.append(dict(**meta,level=level,qtl_reference=qtl,triple_count=len(common),any_triple=bool(common),true_triple=len(common&true),false_triple=len(common-true),FDP=len(common-true)/max(len(common),1),sensitivity=len(common&true)/len(true) if true else None))
        for name,v in d['prediction'].items():prediction.append(dict(**meta,method=name,balanced_accuracy=v['balanced_accuracy'],AUC=v['auc'],Brier=v['brier']))
        prediction.append(dict(**meta,method='DIABLO',balanced_accuracy=r['balanced_accuracy'],AUC=None,Brier=None))
        for j,(p,q,qt,e,et) in enumerate(zip(d['mofa_association']['p'],d['mofa_association']['q'],d['mofa_validation']['q'],d['mofa_association']['effect'],d['mofa_validation']['effect'])):
            factor.append(dict(**meta,factor=j+1,discovery_p=p,discovery_q=q,validation_q=qt,discovery_log_OR=e,validation_log_OR=et,replicated=q<.05 and qt<.05,direction_consistent=e*et>0))
        budgets.append(dict(**meta,k=d['selected_k'],active_factors=sum(d['active_factors']),converged=d['mofa_converged'],genotype_R2=d['reconstruction_r2'][0],RNA_R2=d['reconstruction_r2'][1],diablo_converged=r['converged'],initial_converged=json.loads((ROOT/'results/cohorts'/path.name).read_text())['mofa_converged']))
        for boot,v in enumerate(d['stability']):
            stability.append(dict(**meta,method='MOFA',bootstrap=boot+1,matched_cosine=float(np.mean(v['matched_weight_cosine'])),marker_jaccard=overlap(sets['MOFA'],v['genes'])['jaccard'],converged=v['converged']))
        for boot,v in enumerate(r['stability']):stability.append(dict(**meta,method='DIABLO',bootstrap=boot+1,matched_cosine=None,marker_jaccard=overlap(sets['DIABLO'],v['genes'])['jaccard'],converged=None))
    frames={name:pd.DataFrame(value) for name,value in [('marker_recovery',rows),('prediction',prediction),('factor_associations',factor),('triangulation',tri),('pairwise_overlap',pair),('marker_selection',markers),('bootstrap_stability',stability),('mofa_diagnostics',budgets),('interaction_recovery',interactions)]}
    for name,frame in frames.items():frame.to_csv(ROOT/f'results/{name}.csv',index=False)
    summary=[]
    for name,columns,values in [('prediction',['scenario','n','method'],['balanced_accuracy','AUC','Brier']),('marker_recovery',['scenario','n','method'],['selected','FDP','sensitivity']),('interaction_recovery',['scenario','n','method'],['selected','FDP','sensitivity']),('triangulation',['scenario','n','level','qtl_reference'],['triple_count','any_triple','FDP','sensitivity']),('pairwise_overlap',['scenario','n','level','comparison'],['jaccard'])]:
        combined=pd.concat([summarize(frames[name],columns,v) for v in values],ignore_index=True);combined.to_csv(ROOT/f'results/{name}_summary.csv',index=False)
    frequency=frames['marker_selection'].groupby(['scenario','n','method','gene','gene_index','locus','module'],as_index=False).agg(selection_frequency=('selected','mean'),replicates=('selected','size'),true_genetic=('true_genetic','first'),true_subtype=('true_subtype','first'))
    assert (frequency.replicates==20).all();frequency.to_csv(ROOT/'results/selection_frequency.csv',index=False)
    cal=[]
    for path in sorted((ROOT/'results/calibration').glob('*.json')):
        d=json.loads(path.read_text());assert d['status']=='ok'
        for method,field in [('Ordinary F','q'),('Restricted wild F','wild_q')]:
            for family,j in [('omnibus',0),('interaction',1)]:cal.append(dict(scenario=d['scenario'],n=d['n'],replicate=d['replicate'],method=method,family=family,any_false=bool(np.any(np.array(d[field])[:,j]<.05)),false_genes=int(np.sum(np.array(d[field])[:,j]<.05))))
    cal=pd.DataFrame(cal);assert len(cal)==8000
    calibration=[]
    for key,cell in cal.groupby(['scenario','n','method','family']):
        n=len(cell);k=int(cell.any_false.sum());lo,hi=wilson(k,n)
        calibration.append(dict(zip(['scenario','n','method','family'],key),cohorts=n,rejections=k,rate=k/n,lower=lo,upper=hi))
    cal.to_csv(ROOT/'results/calibration_cohort_metrics.csv',index=False);pd.DataFrame(calibration).to_csv(ROOT/'results/calibration_summary.csv',index=False)
    example=frames['marker_selection'].query("scenario=='convergent' and n==480 and replicate==0").copy();example.to_csv(ROOT/'results/example_markers.csv',index=False)
    print('Summarised 200 independent pairs, 2000 calibration cohorts and every candidate gene')


if __name__=='__main__':main()
