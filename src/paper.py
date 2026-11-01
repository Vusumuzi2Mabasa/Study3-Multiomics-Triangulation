"""Research manuscript generated from every simulation cell, never selected outcomes."""
from pathlib import Path
import json
import pandas as pd
from analysis import wilson

ROOT=Path(__file__).resolve().parents[1]
SC=['null','regulation_only','abundance_only','convergent','interaction_only']
LABEL=['Null','Regulation only','Abundance only','Convergent','Interaction only']


def read(name):return pd.read_csv(ROOT/f'results/{name}.csv',keep_default_na=False,na_values=[''])


def table(head,rows):return '| '+' | '.join(head)+' |\n| '+' | '.join(['---']*len(head))+' |\n'+'\n'.join('| '+' | '.join(map(str,row))+' |' for row in rows)


def main():
    prediction=read('prediction_summary');tri=read('triangulation');cal=read('calibration_summary');diag=read('mofa_diagnostics');freq=read('selection_frequency');factors=read('factor_associations')
    methods=['MOFA factors','PCA scores','DIABLO','RNA ridge','Combined ridge','Clinical only']
    def score(s,m):return float(prediction[(prediction.scenario==s)&(prediction.n==480)&(prediction.method==m)&(prediction.metric=='balanced_accuracy')]['mean'].iloc[0])
    def probability(s,level='gene',reference='eQTL'):
        values=tri[(tri.scenario==s)&(tri.n==480)&(tri.level==level)&(tri.qtl_reference==reference)].any_triple.astype(int)
        k=int(values.sum());n=len(values);lo,hi=wilson(k,n)
        return f'{k}/{n} ({100*k/n:.0f}%; {100*lo:.0f}–{100*hi:.0f}%)'
    designrows=[[LABEL[i],'Yes' if s in ['regulation_only','convergent'] else 'No','Yes' if s in ['abundance_only','convergent'] else 'No','Yes' if s=='interaction_only' else 'No'] for i,s in enumerate(SC)]
    predrows=[[LABEL[i]]+[f'{score(s,m):.3f}' for m in methods] for i,s in enumerate(SC)]
    trirows=[[LABEL[i],probability(s),probability(s,'locus'),probability(s,'module'),probability(s,reference='Wild eQTL')] for i,s in enumerate(SC)]
    calrows=[]
    for s in ['null','abundance_only']:
        for n in [240,480]:
            cells=[]
            for method,family in [('Ordinary F','omnibus'),('Restricted wild F','omnibus'),('Ordinary F','interaction'),('Restricted wild F','interaction')]:
                row=cal[(cal.scenario==s)&(cal.n==n)&(cal.method==method)&(cal.family==family)].iloc[0]
                cells.append(f'{100*row.rate:.1f}% ({100*row.lower:.1f}–{100*row.upper:.1f})')
            calrows.append([LABEL[SC.index(s)],n]+cells)
    markerrows=[]
    for gene in range(6):
        markerrows.append([f'SIM_GENE_{gene+1:02}',f'SIM_LOCUS_{gene//2+1:02}']+[f'{100*float(freq[(freq.scenario=="convergent")&(freq.n==480)&(freq.method==method)&(freq.gene_index==gene)].selection_frequency.iloc[0]):.0f}%' for method in ['MOFA','DIABLO','eQTL']])
    nconv=int(diag.converged.sum());ndconv=int(diag.diablo_converged.sum())
    valid=factors[(factors.discovery_q<.05)&(factors.validation_q<.05)]
    text=r'''# Triangulation of Genotype and RNA Signals in Simulated Breast Cancer

**Vusumuzi Mabasa**

**Simulation study.** Every participant, genotype, LD pattern, gene, coordinate, expression measurement and module is synthetic. The molecular layers are genotype and RNA only. No patient data, proteins, actual breast cancer markers, real pathway enrichment or causal validation are included.

[Word report](report/Study3_Triangulation_Report.docx) · [Quarto manuscript](report/manuscript.qmd) · [Protocol](protocol/design.json) · [Complete marker tables](results/marker_selection.csv) · [Workflow](workflow.txt) · [Example dataset](data/example)

<!-- manuscript -->
## Abstract

We investigated when unsupervised molecular programmes, supervised subtype signatures and local genetic expression associations converge in synthetic breast cancer cases. MOFA+ factors were associated with TNBC using adjusted regression, DIABLO selected supervised genotype–RNA signatures, and regional cis-eQTL models tested common and subtype-dependent effects. Five scenarios separated regulatory associations, abundance shifts, convergent signals, an interaction mechanism and a molecular null. Twenty independent training/test pairs per scenario and training size yielded 200 pairs, each with 360 independent test participants. A further 2,000 fresh cohorts evaluated eQTL calibration. At training N = 480, mean independent-test balanced accuracy in the convergent scenario was '''+f'{score("convergent","MOFA factors"):.3f} for MOFA factors and {score("convergent","DIABLO"):.3f} for DIABLO.'+r''' Triple gene-level overlap occurred in '''+probability('convergent')+r''' convergent cohorts and '''+probability('regulation_only')+r''' regulation-only cohorts. Agreement is conditional on the mechanism and nomination rules; it is not independent biological validation or evidence of a universal method ranking.

## Introduction

Matched genotype and RNA measurements support distinct questions. Latent factors describe molecular variation; supervised signatures describe discrimination; local association models connect genotype to expression. Agreement may strengthen a candidate's interpretation, but a shared eQTL need not distinguish subtype, and subtype-associated RNA need not have a detectable local genetic association.

This simulation represents the genotype–RNA component of Study 3 rather than a complete genomic, transcriptomic and proteomic study. Three parallel analyses are followed by triangulation. MOFA+ does not supply inputs to DIABLO, and eQTL results do not screen DIABLO predictors. Synthetic measurements permit evaluation of agreement and informative disagreement without inventing discoveries about named human genes.

## Methods

### Design and synthetic participants

Training sizes were 240 and 480, with twenty independent replicates per scenario and size. Each training cohort had a separately seeded test cohort of 360 cases. The panels contained 48 dosages in twelve four-variant local regions and 24 RNA genes, two per region. Four synthetic modules each contained six genes. The same first six genes carried scenario-dependent planted effects.

Table 1. Prespecified mechanisms. An interaction changes conditional genotype–expression slopes and need not create a strong marginal mean difference or linearly predictive signature.

'''+table(['Scenario','Common eQTL','Abundance shift','Subtype interaction'],designrows)+r'''

![Parallel discovery signatures and local genetic associations](figures/study_design.png)

Figure 1. Three parallel analyses of the same molecular layers, followed by triangulation. Clinical covariates are adjustment or prediction inputs, not a third molecular assay. Arrows do not assert a causal chain.

Ancestry had three strata with probabilities 0.5, 0.3 and 0.2. Age, a cellular-composition covariate and three batches were observed. Within-region latent Gaussian correlation was 0.65 raised to variant-index distance. Binomial thresholds generated dosages with target allele frequencies 0.18, 0.25, 0.33 and 0.40, plus small ancestry shifts. This is synthetic LD, not an African or other population reference.

Subtype was generated before RNA:

$$
\Pr(T_i=1\mid C_i)=\operatorname{expit}\{-0.6+0.25A_i+0.20F_i+0.15I(S_i=2)\}.
$$

Here A denotes age, F the cellular-composition covariate and S ancestry. Genotype has no direct term in this subtype generator. The experiment concerns case-only subtype associations, not population breast cancer susceptibility or an established genotype-to-disease causal mechanism.

For gene g, one prespecified local variant and one synthetic module contributed to latent log2 abundance:

$$
L_{ig}=5+C_i^\top\eta+0.45H_{i,m(g)}+\beta_g\widetilde G_{iv(g)}
+\gamma_gT_i+\delta_g\widetilde G_{iv(g)}T_i+\epsilon_{ig}.
$$

Module programmes H were independent standard normals; residual standard deviation was 0.5. Genotypes were centred on twice their target allele frequency. Common effects were 0.5, abundance shifts 0.55, and interactions alternated between 0.55 and −0.55 across six genes. Other genes had none of these effects. Counts and analysed expression used the known whole-assay library s:

$$
K_{ig}\sim\operatorname{NB}(\mu_{ig},r=30),\quad
\mu_{ig}=s_i2^{L_{ig}}/10^6,\quad
Y_{ig}=\log_2(10^6K_{ig}/s_i+1).
$$

The library denominator was not the sum of the 24 selected genes. Checks covered valid dosages, finite matrices, integer nonnegative counts, positive libraries, dimensions and subtype representation. These idealised matrices do not represent FASTQ processing, variant calling, imputation or real RNA-seq QC.

### MOFA discovery and subtype association

Official mofapy2 fitted Gaussian genotype and RNA views:

$$
X^{(v)}=ZW^{(v)\top}+E^{(v)}.
$$

Training-only ancestry/batch adjustment and scaling excluded subtype labels. A 15% masked-entry validation selected four or eight factors using equal-weight mean squared reconstruction error across views. Each candidate had PCA and perturbed-PCA starts selected by ELBO. Masks were excluded from preprocessing estimation. The selected budget was refitted to complete training views; independent participants were projected using fixed training transformations, loadings and precisions [1,2].

The initial 150-iteration cap yielded convergence in 52/200 final fits. A disclosed amendment uniformly extended tuning, starts, final fits and fixed-fixture bootstrap fits to at most 1,000 iterations. Original-budget results remain available. No scenarios, effects, seeds, thresholds or factor grids changed. After extension, '''+str(nconv)+r'''/200 final fits met the medium convergence criterion. A numerical amendment zeroed factors with training-score SD below 10⁻⁸ before inference and prediction, preventing standardisation from amplifying ARD-pruned remnants. Test outcomes never determine factor activity. Inactive factor P values were one, retaining the nominal BH family size.

Each active factor was separately assessed with adjusted logistic regression:

$$
\operatorname{logit}\Pr(T_i=1)=\alpha+\theta_kZ_{ik}+C_i^\top\lambda.
$$

Likelihood-ratio tests compared each factor with the covariate-only model; BH adjustment used the within-cohort factor family. These are factor-wise effects, not mutually adjusted effects of all factors. Discovery-significant factors nominated their five largest absolute RNA loadings. The same fixed factors were tested in independent participants; a separate confirmed set required q below 0.05 in both cohorts. Effect-direction agreement was recorded. Loading ranks are not gene-level significance tests.

Total training reconstruction R² and single-factor reconstruction diagnostics describe the representation. Single-factor contributions can be negative and do not form additive explained-variance fractions. Modules are synthetic groups, not named biological pathways.

### DIABLO signatures and downstream prediction

Official R mixOmics DIABLO used genotype and RNA blocks plus a clinical-input block. Threefold stratified training-only validation selected one or two components and three or six molecular features per block per component. Transformations were refitted within folds. The clinical block retained all six covariates; block-design off-diagonal entries were 0.1. Nonzero final RNA loadings defined the DIABLO gene set; genotype selections were also released [3,4].

Independent participants supplied weighted-centroid class predictions. MOFA scores, PCA scores, RNA features, combined genotype/RNA features and clinical-only inputs were evaluated with standardised ridge logistic regression, fixed C = 1. PCA used the nominal selected MOFA budget; MOFA could retain fewer active factors. All references used identical clinical covariates and test participants.

$$
\operatorname{BA}=\tfrac12\{\operatorname{sensitivity}+\operatorname{specificity}\}.
$$

AUC and Brier scores were reported for probabilistic logistic references. DIABLO voting labels were not treated as calibrated probabilities. Prediction from MOFA factors is a downstream analysis, not its unsupervised fitting objective.

### Regional eQTL mapping and calibration

Four local variants were fitted jointly for each gene:

$$
Y_{ig}=\alpha_g+C_i^\top\eta_g+\tau_gT_i+
\sum_{v\in\mathcal C_g}b_{gv}G_{iv}+
\sum_{v\in\mathcal C_g}d_{gv}G_{iv}T_i+e_{ig}.
$$

The omnibus test jointly tested eight main and interaction coefficients against covariates plus subtype. The interaction test retained main genetic effects under its null and tested four interaction coefficients. Nested F statistics used observed design ranks. BH 0.05 was applied separately across 24 genes in each family. These are gene-level regional associations, not fine-mapped causal-variant findings.

Restricted wild-bootstrap sensitivity tests used reduced-model HC2 residuals and participant-level Rademacher multipliers:

$$
Y_i^*=\widehat Y_{0i}+\widehat e_{0i}w_i/\sqrt{1-h_{0ii}},\qquad
\Pr(w_i=1)=\Pr(w_i=-1)=\tfrac12,
$$

$$
\widehat p=\frac{1+\sum_{b=1}^{999}I(F^{*(b)}\ge F_{\mathrm{obs}})}{1000}.
$$

Multipliers were shared across genes/hypotheses within each cohort, and the same BH families were used [5]. Ordinary-F and bootstrap discovery sets were separate. Regional ordinary-F replication required omnibus q below 0.05 in both discovery and independent cohorts; that does not establish causal-allele replication.

A further 2,000 fresh cohorts assessed null and abundance-only scenarios, 500 per scenario and size. Under a wholly null family, the probability of any false gene equals expected FDP. Wilson 95% intervals describe these probabilities. No threshold was tuned from calibration outcomes, and no joint error guarantee is claimed across families or methods.

### Stability and triangulation

Selection frequency across twenty independent cohorts was distinguished from within-cohort bootstrap stability. For the fixed first replicate of each scenario/size, five resamples assessed MOFA loading agreement and marker overlap and DIABLO marker overlap. Transformations were refitted, but factor budgets and DIABLO tuning parameters were retained. MOFA bootstrap fits used one start. These conditional checks do not evaluate full repeated tuning or establish a stability winner.

$$
J(A,B)=|A\cap B|/|A\cup B|.
$$

Empty-union Jaccard was undefined. Locus overlap mapped each gene to its two-gene region. Module nominations required hypergeometric over-representation against the 24-gene background with BH across four six-gene modules. Three-way intersections were evaluated separately at gene, locus and module levels. Overlap from shared measurements is descriptive triangulation, not independent validation.

Truth was objective-specific: MOFA/DIABLO recovery used genes with planted abundance or interaction effects, eQTL recovery used genes with common or interaction effects, and shared-candidate recovery required both. Interaction-only effects can be poorly captured by linear discrimination. MOFA/DIABLO feature FDP is a synthetic recovery diagnostic, not a claim of controlled gene-selection FDR. DIABLO selects features even under a molecular null.

Means and t-based Monte Carlo intervals describe twenty main pairs per cell. Probability estimates use Wilson intervals. Zero overlap in twenty cohorts does not establish zero population probability. Every cell, unfavourable estimate and original-budget result is retained.

## Results

### Factor associations and interpretation

The prespecified convergent N = 480 first replicate illustrates discovery and held-out association. Loadings explain contributions to a factor; association concerns the score rather than automatic significance for every loading. Across the experiment, '''+str(len(valid))+r''' factor instances had q below 0.05 in both cohorts, and '''+str(int(valid.direction_consistent.sum()))+r''' had matching directions. These instances span different models and are not a common set of biological factors.

![MOFA factor associations reconstruction and RNA loadings](figures/mofa_factors.png)

Figure 2. Prespecified convergent fixture: discovery/independent factor q values, single-factor reconstruction and leading RNA weights for the strongest discovery-associated factor within this fixed fixture. Reconstruction contributions are not additive variance fractions. The heatmap colour scale starts at zero; numerical calculations retain any negative values. G identifiers are synthetic.

### Independent subtype prediction

Table 2 evaluates whether representations or signatures add information beyond simpler references. Predictive performance does not establish regulation, and successful eQTL mapping need not yield subtype discrimination. The regulation-only scenario deliberately separates those objectives.

Table 2. Mean independent-test balanced accuracy at training N = 480, twenty pairs per scenario. MOFA/PCA use downstream ridge classifiers; DIABLO uses supervised components and centroid voting.

'''+table(['Scenario','MOFA factors','PCA scores','DIABLO','RNA ridge','Combined ridge','Clinical only'],predrows)+r'''

![Independent subtype prediction across scenarios](figures/prediction.png)

Figure 3. Balanced accuracy and t-based 95% Monte Carlo intervals over twenty independent pairs per cell. The dotted line is 0.5 balanced accuracy, not an assumption that clinical covariates contain no information. Every test participant is synthetic.

### Marker nomination and stability

Table 3 retains the six prespecified signal genes for the convergent scenario. Complete tables include all 24 genes, selected null genes, q values and known truth. Frequencies across generated cohorts are not bootstrap confidence probabilities or clinical reproducibility estimates.

Table 3. Selection frequencies at convergent N = 480 over twenty independent training cohorts. All six genes have planted common regulatory effects and abundance shifts.

'''+table(['Synthetic gene','Synthetic locus','MOFA','DIABLO','eQTL'],markerrows)+r'''

![Selection frequency for all genes in signal scenarios](figures/marker_frequency.png)

Figure 4. N = 480 gene selection frequencies. G01–G06 carry scenario-dependent effects; the dotted line separates them from other genes. Regulation-only eQTL detection without subtype-signature agreement is expected. Null-scenario results remain in the released tables.

Bootstrap checks assess another aspect of reproducibility. MOFA loading agreement used absolute cosine similarity with optimal matching to account for factor permutation and sign. Empty marker-set comparisons were left undefined. Five conditional resamples provide limited precision and retain selected tuning parameters.

![Conditional bootstrap marker stability](figures/bootstrap_stability.png)

Figure 5. Jaccard across five resamples in each of five fixed scenario fixtures per size. Boxes pool those fixtures descriptively, not as twenty-cohort uncertainty estimates. Resampling strategies differ between MOFA and DIABLO.

### Convergence across genes loci and modules

Triple agreement was measured rather than required. Table 4 reports any shared nomination across MOFA, DIABLO and ordinary-F eQTL mapping and a wild-bootstrap sensitivity analysis. A shared module may arise through different genes, but a larger overlap is not automatically better discovery.

Table 4. Any triple overlap at N = 480: events/20, percentage and Wilson 95% interval. Module overlap requires module q below 0.05 separately for each method. Wild gene overlap replaces the ordinary-F eQTL set.

'''+table(['Scenario','Gene','Locus','Module','Wild gene'],trirows)+r'''

![Cross method agreement probabilities](figures/triangulation.png)

Figure 6. Any triple agreement across twenty pairs per cell, with Wilson intervals. Reg, Abund, Conv and Inter denote regulation-only, abundance-only, convergent and interaction-only. Colours distinguish eQTL references. Zero observed events retain nonzero upper bounds.

The fixed convergent fixture relates nominated genes to a known module. No real pathway database or human annotation was used. Synthetic module membership does not independently validate the generator that created it.

![Synthetic module enrichment](figures/module_enrichment.png)

Figure 7. BH-adjusted hypergeometric results in convergent N = 480 replicate zero. SIM MODULE 01 contains the planted signal genes. These are synthetic grouping checks, not breast cancer pathway discoveries.

### Calibration and computational checks

Table 5 reports every calibration cell for this generator and its 24-gene families. The earlier twelve-gene study's calibration is not transplanted as a guarantee. A confidence interval containing 5% does not prove exact control for real RNA.

Table 5. Any false gene after BH with Wilson 95% intervals, 500 fresh cohorts per scenario/size. Both test families are wholly null in these scenarios.

'''+table(['Scenario','N','F omnibus','Wild omnibus','F interaction','Wild interaction'],calrows)+r'''

![Independent eQTL error rate calibration](figures/calibration.png)

Figure 8. Calibration on 2,000 independent synthetic cohorts. BH is separate for each family; intervals are cell-specific rather than simultaneous guarantees.

All 200 pairs and 2,000 calibration cohorts passed data-validity checks. Independent R calculations replayed factor likelihood-ratio tests and seventeen explicit bootstrap statistics. Six unit tests covered generation, BH, regression refits, mapping, logistic fitting and separation of test labels from prediction fitting. Extended MOFA convergence was '''+str(nconv)+r'''/200; DIABLO's warning-based convergence check passed in '''+str(ndconv)+r'''/200 fits. Flags remain released; numerical completion, convergence and scientific validity are distinct.

## Discussion

Programme discovery, discrimination and local genetic association provide complementary evidence. A variant can regulate expression without distinguishing subtype, and subtype-associated RNA need not have a local genetic effect. Overlap depends on the mechanism, power, correlations and nomination rules. Requiring identical markers would discard valid findings about different questions.

For Study 3, MOFA factors can guide programme interpretation, DIABLO can identify signatures evaluated against simpler references, and eQTL mapping can supply regional genotype–RNA links and subtype interactions. Their agreement prioritises candidates, but does not identify causal variants, demonstrate RNA mediation or show that regulation causes TNBC.

Subtype precedes RNA in this generator and has no direct germline term. Supervised genotype selection may therefore reflect covariance with RNA or finite-sample selection rather than a true marginal genetic predictor. Classification accuracy cannot be attributed specifically to regulation without targeted ablations. Interaction-only changes may remain difficult for the linear classifiers evaluated here.

The panel contains only 24 genes, twelve synthetic LD regions, four modules and twenty main replicates per cell. Module agreement depends on the small universe and enrichment rules. Effect sizes are moderate relative to the chosen noise but not calibrated to a clinical assay. Cell composition is measured without error; missing assays, rare variants, relatedness and cross-country transportability were not modelled. A real study requires assay-specific, specimen, ancestry, multiplicity and power assessments.

Convergence and dormant-factor amendments address computational reliability and preserve original outputs. Shared measurements mean agreement is not three independent confirmations. Independent synthetic participants assess reproducibility within the generator. Patient replication, ancestry-appropriate fine-mapping and colocalisation remain separate tasks.

## Conclusion

The workflow implements MOFA discovery, supervised signatures and regional eQTL mapping as parallel components of a genotype–RNA Study 3 simulation. It quantifies convergence and disagreement while retaining objective-specific truth, independent evaluation and calibration uncertainty. It supports methodological study design rather than claiming validated breast cancer biomarkers.

## Reproducibility

The full execution order is in [workflow.txt](workflow.txt). Primary, calibration and test seeds are disjoint. Frozen scientific sources precede primary evaluation; two different-seed implementation pilots are excluded. Original-budget sensitivity results, raw matrices, models and every candidate gene are released. Temporary logs and environments are excluded. Existing records resume without silent replacement.

```bash
python -m pip install -r requirements.txt
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python -m unittest discover -s tests -v
python src/run.py --stage main --workers 2
Rscript --vanilla src/diablo.R .
python src/run.py --stage calibration --workers 2
python src/refine_mofa.py
python src/finalize.py
python src/analysis.py
python src/figures.py
python src/verify.py --fixture-only
Rscript --vanilla src/check_R.R .
python src/verify.py
python src/paper.py
python src/build_word.py
python src/audit_release.py
```

Python versions are pinned. R with mixOmics and jsonlite and Quarto are required. R package versions appear in individual outputs. An isolated full rerun can use `python src/reproduce.py --output /path/to/new-study --rscript /path/to/Rscript --quarto /path/to/quarto`. New scientific questions require documented amendments and distinct result banks.

## References

1. Argelaguet R, et al. MOFA+: a statistical framework for comprehensive integration of multi-modal single-cell data. Genome Biology. 2020;21:111. [Methods paper](https://doi.org/10.1186/s13059-020-02015-1).
2. MOFA developers. [Downstream analysis](https://github.com/bioFAM/MOFA2/blob/master/vignettes/downstream_analysis.Rmd).
3. Singh A, et al. DIABLO: an integrative approach for identifying key molecular drivers from multi-omics assays. Bioinformatics. 2019;35:3055–3062. [Methods paper](https://doi.org/10.1093/bioinformatics/bty1054).
4. mixOmics developers. [Official DIABLO documentation](https://mixomics.org/mixdiablo/).
5. Davidson R, Flachaire E. The wild bootstrap, tamed at last. Journal of Econometrics. 2008;146:162–169. [Methods paper](https://doi.org/10.1016/j.jeconom.2008.08.003).

## Data and code availability

The public repository contains synthetic measurements and annotations only, with no private proposal, clinical data, credentials or controlled-access reference. Released under MIT at [Study3 Multiomics Triangulation](https://github.com/Vusumuzi2Mabasa/Study3-Multiomics-Triangulation). The [earlier genotype–RNA study](https://github.com/Vusumuzi2Mabasa/Breast-Cancer-Genotype-RNA-eQTL) remains unchanged by this project.
'''
    interaction=read('interaction_recovery_summary')
    ipower=float(interaction[(interaction.scenario=='interaction_only')&(interaction.n==480)&(interaction.method=='Restricted wild F')&(interaction.metric=='sensitivity')]['mean'].iloc[0])
    text=text.replace('Table 2 evaluates whether representations or signatures add information beyond simpler references.',f'Abundance shifts supported subtype discrimination: at N = 480, DIABLO balanced accuracy was {score("abundance_only","DIABLO"):.3f} in the abundance-only scenario and {score("convergent","DIABLO"):.3f} in the convergent scenario. Corresponding MOFA-factor values were {score("abundance_only","MOFA factors"):.3f} and {score("convergent","MOFA factors"):.3f}. In the convergent scenario, RNA ridge reached {score("convergent","RNA ridge"):.3f}; these results do not show a general advantage from adding genotype. Table 2 evaluates whether representations or signatures add information beyond simpler references.')
    text=text.replace('### Marker nomination and stability',f'In contrast, the interaction-only scenario had mean wild-bootstrap interaction sensitivity of {100*ipower:.1f}% at N = 480, while DIABLO balanced accuracy was {score("interaction_only","DIABLO"):.3f}. Detectable subtype-dependent regulatory slopes therefore did not imply a useful linear classifier. These interaction results are retained in [the interaction summary](results/interaction_recovery_summary.csv).\n\n### Marker nomination and stability')
    text=text.replace('Triple agreement was measured rather than required.', 'At N = 480, all twenty convergent cohorts had at least one triple-supported gene, while only six had a triple-supported enriched module. The corresponding module proportion was 30%, not 100%. The stricter module over-representation rule and breadth of nominated gene sets make these different endpoints; pathway-level agreement is not automatically more frequent than gene overlap. Regulation-only and interaction-only each showed gene overlap in one of twenty cohorts, and null and abundance-only showed none. Triple agreement was measured rather than required.')
    text=text.replace('Table 5 reports every calibration cell for this generator and its 24-gene families.', 'Ordinary-F interaction rejection was elevated in the N = 240 molecular-null cell: 7.2% (95% interval 5.2–9.8%), versus 4.8% (3.2–7.0%) for restricted wild bootstrap. The bootstrap rate ranged from 4.0% to 5.6% across the eight calibration cells; this is encouraging within the generator, not universal calibration. Table 5 reports every calibration cell for this generator and its 24-gene families.')
    (ROOT/'README.md').write_text(text,encoding='utf-8',newline='\n')
    print('Research manuscript generated from all scenarios')


if __name__=='__main__':main()
