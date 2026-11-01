"""Publication figures, with synthetic modules explicitly distinguished from biology."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from statistics import module_enrichment
from analysis import wilson
from mofa_engine import prep_fit

ROOT=Path(__file__).resolve().parents[1]
SC=['null','regulation_only','abundance_only','convergent','interaction_only']
LABEL=['Null','Regulation only','Abundance only','Convergent','Interaction only']
METHODS=['MOFA','DIABLO','eQTL'];COLORS=['#2166ac','#b2182b','#1b7837']


def read(name):return pd.read_csv(ROOT/f'results/{name}.csv',keep_default_na=False,na_values=[''])


def save(fig,name):
    fig.savefig(ROOT/f'figures/{name}.png',dpi=180,bbox_inches='tight',facecolor='white')
    p=ROOT/f'figures/{name}.svg';fig.savefig(p,bbox_inches='tight',facecolor='white');plt.close(fig)
    p.write_text('\n'.join(line.rstrip() for line in p.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8',newline='\n')


def main():
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(figsize=(8,4));ax.axis('off')
    ax.text(.5,.95,'Matched genotype and RNA measurements',ha='center',fontsize=15,weight='bold')
    ax.text(.5,.86,'All participants, genes, LD and modules are synthetic',ha='center',fontsize=10)
    for x,title,body,color in zip([.16,.5,.84],['MOFA+','DIABLO','cis-eQTL'],['Unsupervised factors\nAdjusted subtype associations\nHeld-out association checks','Training-only supervised tuning\nIndependent prediction\nFeature selection stability','Regional genetic associations\nSubtype-dependent slopes\nIndependent null calibration'],COLORS):
        ax.annotate('',xy=(x,.68),xytext=(.5,.81),arrowprops=dict(arrowstyle='->',color='grey'))
        ax.text(x,.61,title,ha='center',fontsize=14,color=color,weight='bold');ax.text(x,.42,body,ha='center',va='center',fontsize=9)
        ax.annotate('',xy=(.5,.15),xytext=(x,.27),arrowprops=dict(arrowstyle='->',color='grey'))
    ax.text(.5,.06,'Triangulation at gene, locus and synthetic module levels',ha='center',fontsize=12,weight='bold');save(fig,'study_design')
    d=json.loads((ROOT/'results/final/convergent_480_00.json').read_text());a=np.load(ROOT/'results/final_models/convergent_480_00.npz')
    fig,axes=plt.subplots(3,1,figsize=(7.6,7.3),gridspec_kw={'height_ratios':[1,1,1.2]})
    x=np.arange(d['selected_k'])
    axes[0].bar(x-.17,-np.log10(np.maximum(d['mofa_association']['q'],1e-15)),.34,label='Discovery',color=COLORS[0]);axes[0].bar(x+.17,-np.log10(np.maximum(d['mofa_validation']['q'],1e-15)),.34,label='Independent participants',color=COLORS[1]);axes[0].axhline(-np.log10(.05),ls=':',color='grey');axes[0].set(ylabel='−log10 factor q',xticks=x,xticklabels=x+1,title='Prespecified convergent fixture, N=480');axes[0].legend(fontsize=8)
    folder=ROOT/'data/cohorts/convergent_480_00';c=np.loadtxt(folder/'train_covariates.csv',delimiter=',');views=[np.loadtxt(folder/f'train_{v}.csv',delimiter=',') for v in ['genotype','RNA']];prepared,_=prep_fit(views,c)
    contributions=np.array([[1-np.sum((prepared[v]-np.outer(a['scores'][:,k],a[name][:,k]))**2)/np.sum(prepared[v]**2) for k in x] for v,name in enumerate(['weights_G','weights_RNA'])])
    pd.DataFrame(contributions.T,columns=['genotype_single_factor_R2','RNA_single_factor_R2']).assign(factor=x+1).to_csv(ROOT/'results/example_factor_reconstruction.csv',index=False)
    image=axes[1].imshow(contributions,aspect='auto',cmap='Blues',vmin=0,vmax=max(.2,contributions.max()));axes[1].set(xticks=x,xticklabels=x+1,yticks=[0,1],yticklabels=['Genotype','RNA'],xlabel='Factor',title='Single factor reconstruction R², not additive variance fractions');fig.colorbar(image,ax=axes[1],fraction=.025,pad=.02)
    strongest=int(np.argmin(d['mofa_association']['q']));w=a['weights_RNA'][:,strongest];ix=np.argsort(-np.abs(w))[:8]
    axes[2].bar(np.arange(len(ix)),w[ix],color=COLORS[0]);axes[2].set(xticks=np.arange(len(ix)),xticklabels=[f'G{j+1:02}' for j in ix],ylabel='RNA loading',title=f'Largest RNA weights in factor {strongest+1}; G labels are synthetic')
    fig.tight_layout();save(fig,'mofa_factors')
    predictions=read('prediction_summary');methods=['MOFA factors','PCA scores','DIABLO','RNA ridge','Combined ridge','Clinical only'];colors=['#2166ac','#762a83','#b2182b','#1b7837','#e08214','#777777']
    fig,axes=plt.subplots(2,1,figsize=(8,5.8),sharex=True)
    for ax,n in zip(axes,[240,480]):
        for j,(method,color) in enumerate(zip(methods,colors)):
            cell=predictions[(predictions.n==n)&(predictions.method==method)&(predictions.metric=='balanced_accuracy')].set_index('scenario').loc[SC]
            ax.errorbar(np.arange(5)+(j-2.5)*.10,cell['mean'],yerr=[cell['mean']-cell.lower,cell.upper-cell['mean']],fmt='o',ms=3,capsize=2,label=method,color=color)
        ax.axhline(.5,color='grey',ls=':');ax.set(title=f'Independent prediction, training N={n}',ylabel='Balanced accuracy',ylim=(.4,.9))
    axes[-1].set(xticks=np.arange(5),xticklabels=LABEL);h,l=axes[0].get_legend_handles_labels();fig.legend(h,l,loc='upper center',ncol=3,bbox_to_anchor=(.5,1.06),fontsize=8);fig.tight_layout();save(fig,'prediction')
    frequency=read('selection_frequency');fig,axes=plt.subplots(1,4,figsize=(9,5.3),sharey=True)
    for ax,scenario in zip(axes,SC[1:]):
        matrix=frequency[(frequency.n==480)&(frequency.scenario==scenario)&frequency.method.isin(METHODS)].pivot(index='gene_index',columns='method',values='selection_frequency')[METHODS]
        image=ax.imshow(matrix,cmap='Blues',vmin=0,vmax=1,aspect='auto');ax.set(title=LABEL[SC.index(scenario)],xticks=[0,1,2],xticklabels=['MOFA','DIABLO','eQTL'],yticks=np.arange(24),yticklabels=[f'G{j+1:02}' for j in range(24)]);ax.tick_params(axis='x',rotation=50)
        ax.axhline(5.5,color='black',ls=':',lw=.8)
    fig.colorbar(image,ax=axes.tolist(),fraction=.025,pad=.02,label='Selection frequency across 20 independent cohorts');save(fig,'marker_frequency')
    tri=read('triangulation');fig,axes=plt.subplots(2,3,figsize=(9,5.5),sharey=True)
    for row,n in enumerate([240,480]):
        for col,level in enumerate(['gene','locus','module']):
            ax=axes[row,col]
            for j,(method,color) in enumerate(zip(['eQTL','Wild eQTL'],['#2166ac','#b2182b'])):
                rates=[];lo=[];hi=[]
                for scenario in SC:
                    values=tri[(tri.n==n)&(tri.scenario==scenario)&(tri.level==level)&(tri.qtl_reference==method)].any_triple.astype(int)
                    rate=float(values.mean());l,h=wilson(int(values.sum()),len(values));rates.append(rate);lo.append(rate-l);hi.append(h-rate)
                ax.errorbar(np.arange(5)+(j-.5)*.15,rates,yerr=np.maximum([lo,hi],0),fmt='o',ms=3,capsize=2,label=method,color=color)
            ax.set(title=f'{level.capitalize()}, N={n}',ylim=(-.03,1.04),xticks=np.arange(5),xticklabels=['Null','Reg','Abund','Conv','Inter']);ax.tick_params(axis='x',rotation=35)
            if col==0:ax.set_ylabel('Probability of any triple overlap')
    h,l=axes[0,0].get_legend_handles_labels();fig.legend(h,l,loc='upper center',ncol=2,bbox_to_anchor=(.5,1.04));fig.tight_layout();save(fig,'triangulation')
    cal=read('calibration_summary');fig,axes=plt.subplots(2,2,figsize=(7,5.5),sharey=True)
    for row,n in enumerate([240,480]):
        for col,family in enumerate(['omnibus','interaction']):
            ax=axes[row,col]
            for j,(method,color) in enumerate(zip(['Ordinary F','Restricted wild F'],['#2166ac','#b2182b'])):
                cell=cal[(cal.n==n)&(cal.family==family)&(cal.method==method)].set_index('scenario').loc[['null','abundance_only']]
                ax.errorbar(np.arange(2)+(j-.5)*.14,cell.rate,yerr=[cell.rate-cell.lower,cell.upper-cell.rate],fmt='o',capsize=2,ms=4,color=color,label=method)
            ax.axhline(.05,color='grey',ls=':');ax.set(title=f'{family.capitalize()}, N={n}',xticks=[0,1],xticklabels=['Null','Abundance only'],ylim=(0,max(.15,cal.upper.max()+.01)),ylabel='Any false gene after BH')
    h,l=axes[0,0].get_legend_handles_labels();fig.legend(h,l,loc='upper center',ncol=2,bbox_to_anchor=(.5,1.04));fig.tight_layout();save(fig,'calibration')
    stability=read('bootstrap_stability');fig,axes=plt.subplots(1,2,figsize=(7,3.6))
    for ax,method in zip(axes,['MOFA','DIABLO']):
        cell=stability[stability.method==method];groups=[cell[cell.n==n].marker_jaccard.dropna() for n in [240,480]]
        ax.boxplot(groups,tick_labels=['N240','N480'],showmeans=True);ax.set(title=f'{method} fixed-fixture bootstrap markers',ylabel='Jaccard with original signature',ylim=(-.03,1.03))
    fig.tight_layout();save(fig,'bootstrap_stability')
    example=read('example_markers');fig,ax=plt.subplots(figsize=(6,3.4));matrix=[]
    for method in METHODS:
        genes=example[(example.method==method)&example.selected].gene_index.tolist();enrichment=module_enrichment(genes);matrix.append(-np.log10(np.maximum(enrichment['q'],1e-15)))
    image=ax.imshow(matrix,cmap='YlGnBu',vmin=0,aspect='auto');ax.set(xticks=range(4),xticklabels=['SIM MODULE 01','SIM MODULE 02','SIM MODULE 03','SIM MODULE 04'],yticks=range(3),yticklabels=METHODS,title='Synthetic module enrichment in a prespecified fixture');ax.tick_params(axis='x',rotation=15);fig.colorbar(image,ax=ax,label='−log10 module q');fig.tight_layout();save(fig,'module_enrichment')
    print('Eight main scientific figures written')


if __name__=='__main__':main()
