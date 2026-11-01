"""Release integrity and cross-language fixtures, independent of marker overlap."""
from pathlib import Path
import argparse,hashlib,json
import numpy as np
from simulation import generate,SCENARIOS,qc
from statistics import bh
from mofa_engine import classifier

ROOT=Path(__file__).resolve().parents[1]


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--fixture-only',action='store_true');a=ap.parse_args()
    phase='final' if (ROOT/'results/final/null_240_00.json').exists() else 'cohorts'
    modelphase='final_models' if phase=='final' else 'models'
    baseline=json.loads((ROOT/f'results/{phase}/null_240_00.json').read_text())
    model=np.load(ROOT/f'results/{modelphase}/null_240_00.npz')
    fixture=baseline['eqtl']['fixture'];fixture.update(z=model['scores'].tolist(),covariates=model['train_cov'].tolist(),outcome=model['train_y'].tolist(),factor_p=baseline['mofa_association']['p'])
    (ROOT/'results/verification_fixture.json').write_text(json.dumps(fixture)+'\n',encoding='utf-8',newline='\n')
    if a.fixture_only:return
    lock=json.loads((ROOT/'protocol/source_lock.json').read_text())
    for name,digest in lock.items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    amendment=json.loads((ROOT/'protocol/convergence_amendment.json').read_text())
    assert hashlib.sha256((ROOT/'src/refine_mofa.py').read_bytes()).hexdigest()==amendment['source_sha256']
    dormant=json.loads((ROOT/'protocol/dormant_factor_amendment.json').read_text())
    assert hashlib.sha256((ROOT/'src/finalize.py').read_bytes()).hexdigest()==dormant['source_sha256']
    assert hashlib.sha256((ROOT/'src/mofa_engine.py').read_bytes()).hexdigest()==hashlib.sha256((ROOT/'provenance/mofa_engine_original.py').read_bytes()).hexdigest()
    seeds=set();converged=0;diablo_converged=0;replications=0
    for path in sorted((ROOT/'results/final').glob('*.json')):
        row=json.loads(path.read_text());original=json.loads((ROOT/'results/cohorts'/path.name).read_text());r=json.loads((ROOT/'results/diablo'/path.name).read_text())
        assert row['status']==r['status']=='ok'
        assert row['eqtl']==original['eqtl'] and row['eqtl_validation']==original['eqtl_validation']
        assert row['seed'] not in seeds and row['test_seed'] not in seeds;seeds.update([row['seed'],row['test_seed']])
        for entry in [row['eqtl'],row['eqtl_validation']]:
            for pfield,qfield in [('p','q')]+([('wild_p','wild_q')] if 'wild_p' in entry else []):
                p=np.array(entry[pfield]);q=np.array(entry[qfield]);assert p.shape==q.shape==(24,2)
                assert np.isfinite(p).all() and np.all((p>=0)&(p<=1))
                np.testing.assert_allclose(q,np.column_stack([bh(p[:,j]) for j in range(2)]),atol=1e-15)
        z=np.load(ROOT/'results/final_models'/(path.stem+'.npz'))
        train=dict(y=z['train_y']);test=dict(y=z['test_y'])
        metrics,prob=classifier(np.c_[z['train_cov'],z['scores']],np.c_[z['test_cov'],z['test_scores']],train,test)
        np.testing.assert_allclose(prob,row['prediction']['MOFA factors']['probability'],atol=1e-12)
        assert len(r['test_class'])==360 and set(r['test_class'])<={0,1}
        y=np.array(row['test_y']);pred=np.array(r['test_class']);ba=.5*(np.mean(pred[y==0]==0)+np.mean(pred[y==1]==1))
        assert abs(ba-r['balanced_accuracy'])<1e-12
        converged+=row['mofa_converged'];diablo_converged+=r['converged'];replications+=1
    assert replications==200
    for path in (ROOT/'results/calibration').glob('*.json'):
        row=json.loads(path.read_text());assert row['seed'] not in seeds;seeds.add(row['seed'])
    assert len(seeds)==2400
    assert json.loads((ROOT/'results/R_validation.json').read_text())['status']=='passed'
    report=dict(status='passed',training_test_pairs=200,independent_test_participants_per_pair=360,fresh_calibration_cohorts=2000,distinct_cohort_seeds=2400,MOFA_final_converged=converged,DIABLO_converged=diablo_converged,original_iteration_budget_records_preserved=200,scientific_source_locks_verified=True,probability_replays=200,DIABLO_balanced_accuracy_replays=200,disclosure='Every measurement, gene, coordinate and module is synthetic; no causal or clinical validation claim.')
    (ROOT/'results/audit.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps(report,indent=2))


if __name__=='__main__':main()
