import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from simulation import generate,qc
from statistics import bh,eqtl,levels,overlap,logistic
from mofa_engine import classifier


class Science(unittest.TestCase):
    def test_seed_and_truth(self):
        a=generate(240,991,'convergent');b=generate(240,991,'convergent')
        np.testing.assert_array_equal(a['RNA'],b['RNA']);qc(a)
        self.assertEqual(np.count_nonzero(a['beta']),6);self.assertEqual(np.count_nonzero(a['gamma']),6)
        self.assertEqual(np.count_nonzero(generate(240,991,'regulation_only')['gamma']),0)
        self.assertEqual(np.count_nonzero(generate(240,991,'interaction_only')['delta']),6)

    def test_bh(self):
        np.testing.assert_allclose(bh([.01,.04,.03,.8]),[.04,.053333333333,.053333333333,.8])

    def test_bootstrap_refits(self):
        d=generate(240,995,'null');r=eqtl(d,99,996,True);a=r['fixture']
        x=np.array(a['x']);y=np.array(a['y']);x0=x[:,:a['p0']]
        def statistic(y):
            e0=y-x0@np.linalg.lstsq(x0,y,rcond=None)[0];e=y-x@np.linalg.lstsq(x,y,rcond=None)[0]
            return ((e0@e0-e@e)/(x.shape[1]-x0.shape[1]))/((e@e)/(len(y)-x.shape[1]))
        self.assertAlmostEqual(statistic(y),a['F'],11)
        q0=np.linalg.qr(x0,mode='reduced')[0];mu=q0@(q0.T@y);res=(y-mu)/np.sqrt(1-np.sum(q0*q0,axis=1))
        for w,value in zip(np.array(a['weights']).T,a['stars']):self.assertAlmostEqual(statistic(mu+res*w),value,10)
        self.assertTrue(np.all(np.asarray(r['wild_p'])>=.01))

    def test_mapping(self):
        self.assertEqual(levels([0,1,6])['locus'],{0,3});self.assertEqual(levels([0,1,6])['module'],{0,1})
        self.assertIsNone(overlap([],[])['jaccard']);self.assertEqual(overlap([0,1],[1,2])['jaccard'],1/3)

    def test_logistic_gradient(self):
        d=generate(240,997,'null');ll,b=logistic(np.c_[np.ones(240),d['covariates']],d['y'])
        self.assertTrue(np.isfinite(ll));self.assertEqual(len(b),7)

    def test_test_labels_do_not_affect_predictions(self):
        rng=np.random.default_rng(998);x=rng.normal(size=(120,3));xt=rng.normal(size=(80,3))
        y=rng.binomial(1,.45,120);yt=rng.binomial(1,.45,80)
        _,a=classifier(x,xt,dict(y=y),dict(y=yt));_,b=classifier(x,xt,dict(y=y),dict(y=1-yt))
        np.testing.assert_array_equal(a,b)


if __name__=='__main__':unittest.main()
