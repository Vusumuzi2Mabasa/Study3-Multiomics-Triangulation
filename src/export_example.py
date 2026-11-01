"""Export a fixed first-replicate cohort with explicit synthetic identifiers."""
from pathlib import Path
import json
import pandas as pd
from simulation import generate,qc,SCENARIOS

ROOT=Path(__file__).resolve().parents[1]


def main():
    cfg=json.loads((ROOT/'protocol/design.json').read_text());seed=cfg['seed_base']+SCENARIOS.index('convergent')*1000000+480000
    d=generate(480,seed,'convergent');out=ROOT/'data/example';out.mkdir(parents=True,exist_ok=True)
    ids=[f'SIM_PARTICIPANT_{i+1:03}' for i in range(480)]
    for name in ['genotype','RNA','counts']:
        columns=[f'SIM_VARIANT_{j+1:02}' for j in range(48)] if name=='genotype' else [f'SIM_GENE_{j+1:02}' for j in range(24)]
        frame=pd.DataFrame(d[name],columns=columns);frame.insert(0,'participant',ids);frame.to_csv(out/f'{name}.csv',index=False)
    frame=pd.DataFrame(d['covariates'],columns=['age_z','cell_covariate_z','ancestry_B','ancestry_C','batch_B','batch_C']);frame.insert(0,'participant',ids);frame['TNBC']=d['y'];frame['library_size']=d['library'];frame.to_csv(out/'metadata.csv',index=False)
    rows=[]
    for j in range(24):rows.append(dict(gene=f'SIM_GENE_{j+1:02}',locus=f'SIM_LOCUS_{j//2+1:02}',module=f'SIM_MODULE_{j//6+1:02}',synthetic_TSS=115000+(j%2)*10000,planted_variant=f'SIM_VARIANT_{4*(j//2)+(0 if j%2==0 else 2)+1:02}',beta=d['beta'][j],delta=d['delta'][j],gamma=d['gamma'][j]))
    pd.DataFrame(rows).to_csv(out/'gene_truth.csv',index=False)
    pd.DataFrame([dict(variant=f'SIM_VARIANT_{j+1:02}',locus=f'SIM_LOCUS_{j//4+1:02}',synthetic_position=100000+(j%4)*10000) for j in range(48)]).to_csv(out/'variant_map.csv',index=False)
    manifest=dict(disclosure='All measurements, participants, coordinates, genes and modules are synthetic.',scenario='convergent',seed=seed,replicate=0,**qc(d))
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8',newline='\n')


if __name__=='__main__':main()
