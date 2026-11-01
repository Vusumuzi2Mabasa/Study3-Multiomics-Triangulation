"""Suppress numerically dormant factors before downstream inference/prediction."""
from pathlib import Path
import json,hashlib
import numpy as np
from statistics import factor_associations
from mofa_engine import classifier
from run import save,nominations

ROOT=Path(__file__).resolve().parents[1]


def main():
    amendment=dict(reason='Cross-language checks found ARD-pruned factor scores around 1e-130. A generic standardizer can amplify these numerical remnants into predictors. Uniformly zero scores with discovery SD below 1e-8 before factor associations and downstream classification; use discovery alone to define active factors. Preserve all raw fits and predictions.',threshold=1e-8,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    manifest=ROOT/'protocol/dormant_factor_amendment.json'
    if manifest.exists():assert json.loads(manifest.read_text())==amendment
    else:save(manifest,amendment)
    files=sorted((ROOT/'results/refined').glob('*.json'));assert len(files)==200
    for path in files:
        row=json.loads(path.read_text());source=ROOT/'results/refined_models'/(path.stem+'.npz')
        arrays=dict(np.load(source));za=arrays['scores'].copy();zb=arrays['test_scores'].copy();active=np.std(za,axis=0)>=1e-8
        za[:,~active]=0;zb[:,~active]=0
        # Constant columns have no tested coefficient; retain nominal-budget q=1 for reporting.
        assoc=factor_associations(za,arrays['train_cov'],arrays['train_y']);validation=factor_associations(zb,arrays['test_cov'],arrays['test_y'])
        for result in [assoc,validation]:
            for j in np.flatnonzero(~active):result['p'][j]=1.;result['effect'][j]=0.
            from statistics import bh
            result['q']=bh(result['p']).tolist()
        model={'weights':[arrays['weights_G'],arrays['weights_RNA']]}
        row.update(active_factors=active.tolist(),mofa_association=assoc,mofa_validation=validation,mofa_genes=nominations(model,assoc),mofa_confirmed_genes=nominations(model,dict(q=np.maximum(assoc['q'],validation['q']).tolist())))
        metrics,p=classifier(np.c_[arrays['train_cov'],za],np.c_[arrays['test_cov'],zb],dict(y=arrays['train_y']),dict(y=arrays['test_y']))
        row['prediction']['MOFA factors']=dict(**metrics,probability=p.tolist())
        save(ROOT/'results/final'/path.name,row)
        out=ROOT/'results/final_models'/source.name;out.parent.mkdir(parents=True,exist_ok=True)
        arrays['scores']=za;arrays['test_scores']=zb;np.savez_compressed(out,**arrays)
    print('Finalised 200 pairs using a training-only active-factor rule')


if __name__=='__main__':main()
