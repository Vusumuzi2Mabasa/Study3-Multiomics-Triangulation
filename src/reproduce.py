"""Reproduce all scientific stages in an explicitly NEW output directory."""
from pathlib import Path
import argparse,os,shutil,subprocess,sys

ROOT=Path(__file__).resolve().parents[1]


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);ap.add_argument('--rscript',default=shutil.which('Rscript'));ap.add_argument('--quarto',default=shutil.which('quarto'));a=ap.parse_args()
    if not a.rscript or not a.quarto:ap.error('Rscript and Quarto are required')
    target=a.output.resolve()
    if target.exists():ap.error('Choose a new directory; existing results are never deleted')
    target.mkdir(parents=True)
    for name in ['src','tests','protocol','provenance']:shutil.copytree(ROOT/name,target/name,ignore=shutil.ignore_patterns('__pycache__'))
    for name in ['requirements.txt','LICENSE']:shutil.copyfile(ROOT/name,target/name)
    for name in ['work','results','figures','report','data/example']:(target/name).mkdir(parents=True,exist_ok=True)
    env=os.environ.copy()
    for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:env[key]='1'
    def run(command):subprocess.run(command,cwd=target,env=env,check=True)
    run([sys.executable,'-m','unittest','discover','-s','tests','-v'])
    run([sys.executable,'src/run.py','--stage','main','--workers','2'])
    run([a.rscript,'--vanilla','src/diablo.R','.'])
    run([sys.executable,'src/run.py','--stage','calibration','--workers','2'])
    for name in ['refine_mofa','finalize','analysis','figures','export_example']:run([sys.executable,f'src/{name}.py'])
    run([sys.executable,'src/verify.py','--fixture-only']);run([a.rscript,'--vanilla','src/check_R.R','.'])
    run([sys.executable,'src/verify.py']);run([sys.executable,'src/paper.py']);run([sys.executable,'src/build_word.py','--quarto',a.quarto])
    print('Scientific reproduction finished. Render and visually review the Word report before generating publication manifests.',target)


if __name__=='__main__':main()
