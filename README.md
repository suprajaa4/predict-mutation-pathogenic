
#  Mutation Effect Predictor 

**Predicting pathogenic vs benign missense variants using biological features**

Why missense mutation?

A missense mutation is a single nucleotide change that results in the substitution of one amino acid for another in a protein. Depending on where  and how drastic the change is, it may disrupt protein folding, affect catalytic or binding sites or have little to no functional effect. They are a major cause of inherited disease. However, many are benign. Accurately distinguishing **pathogenic** from **benign** variants is a central problem in medical genomics.

It's hard to predict as protein function depends on evolutionary conservation or physicochemical changes. No single feature is sufficient so models must combine multiple biological signals.

In this project, I built an end-to-end machine learning pipeline to:

* Construct a  missense variant dataset from ClinVar and extract protein-level features using Ensembl VEP
* Train and evaluate models to predict variant pathogenicity and interpret feature importance 

**Datasets used**

https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz

used the **ClinVar GRCh37 VCF** as the raw data source. ClinVar is a public database that aggregates variant interpretations submitted by clinical labs and expert panels.

Each variant includes: Genomic coordinates, Reference and alternate alleles and Clinical significance (e.g. Benign, Pathogenic)

After preprocessing:

* **~1.4 million labeled variants**
* **~12–15% pathogenic**, reflecting imbalance


**Features**

Features are extracted using **Ensembl Variant Effect Predictor (VEP)** on the GRCh37 assembly.

Protein-level features are extracted using **Ensembl Variant Effect Predictor (VEP)**:

* **Grantham distance** : how different the two amino acids are
* **Protein position** : amino acid index in the protein
* **SIFT score** : conservation-based deleteriousness score
* **PolyPhen score** : predicted structural/functional impact

Only variants annotated as "missense" are used for modeling.

**Models**

Two models are trained:

* **Logistic regression** using Grantham distance and protein position as a simple baseline
* **XGBoost** using all features to capture nonlinear interactions

Because pathogenic variants are rarer, evaluation focuses on both ROC-AUC and precision–recall.

**Results**

The baseline logistic regression performs reasonably well, showing that simple protein-level features already carry signal.

The XGBoost model performs much better, especially in precision–recall space, indicating that combining conservation scores with protein context improves prediction.

**Model performance**
               
 Logistic regression :
 ROC-AUC: ~0.77 ,  PR-AUC :~0.36     
 XGBoost:
 ROC-AUC:~0.94 ,  PR-AUC : ~0.79 

The nonlinear model provides a large performance gain, indicating strong feature interactions.


**Feature importance (XGBoost)**
           
**SIFT score** : ~0.63      
**PolyPhen score**:  ~0.17      
Protein position: ~0.12      
Grantham distance:  ~0.08  

Feature importance from XGBoost shows:

* SIFT is the most important feature-Conservation-based predictors dominate, as expected biologically
* PolyPhen adds additional signal (position importance)
* Protein position and Grantham distance contribute smaller but meaningful effects

Overall, the results align with biological expectations: conserved positions and functionally disruptive changes are more likely to be pathogenic.


**Limitations**

* Conservation features are derived from existing predictors, not raw MSAs
* Structural features (like, solvent accessibility) are not yet included
* Only missense variants are modeled
* ClinVar labels may contain noise


**References**

1. Landrum et al., *ClinVar: improving access to variant interpretations*, Nucleic Acids Research
2. McLaren et al., *The Ensembl Variant Effect Predictor*, Genome Biology
3. Grantham, *Amino acid difference formula*, Science (1974)
4. Ng & Henikoff, *SIFT*, Nucleic Acids Research
5. Adzhubei et al., *PolyPhen-2*, Nature Methods
6. Richards et al., *ACMG guidelines for variant interpretation*, Genetics in Medicine


Any feedback, please  reach out to Suprajaa V | suprajaav4@gmail.com


