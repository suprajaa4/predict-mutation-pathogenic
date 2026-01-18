
#  Mutation Effect Predictor 

**Predicting pathogenic vs benign missense variants using biologically interpretable features**

Single nucleotide variants (SNVs) that alter protein sequences (missense mutations) are a major cause of inherited disease. However, not all missense mutations are harmful — many are benign. Accurately distinguishing **pathogenic** from **benign** variants is a central problem in medical genomics.

In this project, I built an end-to-end machine learning pipeline to:

* Construct a  missense variant dataset from ClinVar and extract protein-level features using Ensembl VEP
* Train and evaluate models to predict variant pathogenicity and interpret feature importance in a biologically meaningful way

The emphasis is on **interpretability, correctness, and reproducibility**, rather than deep learning complexity.

##  Background (Beginner-Friendly)

### What is a missense mutation?

A missense mutation is a single nucleotide change that results in the substitution of one amino acid for another in a protein. Depending on where it occurs and how drastic the change is, it may disrupt protein folding, affect catalytic or binding sites or have little to no functional effect

### Why is this hard to predict?

Protein function depends on:

* **Evolutionary conservation** (important residues are conserved)
* **Physicochemical changes** (e.g. charge, size, hydrophobicity)
* **Protein context** (position within domains or motifs)

No single feature is sufficient — models must combine multiple biological signals.



## Datasets used

## ClinVar (VCF)

used the **ClinVar GRCh37 VCF** as the raw data source. ClinVar is a public database that aggregates variant interpretations submitted by clinical labs and expert panels.

Each variant includes: Genomic coordinates, Reference and alternate alleles and Clinical significance (e.g. Benign, Pathogenic)


### Label construction


* **Benign / Likely benign → label = 0**
* **Pathogenic / Likely pathogenic → label = 1**

Variants with conflicting or uncertain interpretations are excluded.

After preprocessing:

* **~1.4 million labeled SNVs**
* **~12–15% pathogenic**, reflecting imbalance


##  Feature Engineering

Features are extracted using **Ensembl Variant Effect Predictor (VEP)** on the GRCh37 assembly.

### Protein-level features used

| Feature               | Description                                         |
| --------------------- | --------------------------------------------------- |
| **Grantham distance** | Physicochemical severity of amino acid substitution |
| **Protein position**  | Amino acid index within the protein                 |
| **SIFT score**        | Conservation-based deleteriousness predictor        |
| **PolyPhen score**    | Structure/function impact predictor                 |



Only **missense variants with valid amino acid changes** are retained for modeling



##  Project Structure

```
mutation-effect-predictor/
│
├── data/
│   ├── raw/                # ClinVar VCF
│   ├── processed/          # Labeled SNVs, ML table
│   └── external/           # VEP-derived features
│
├── src/
│   ├── preprocess_vcf.py   # Build labeled SNV dataset
│   ├── build_features.py   # Call VEP, extract missense features
│   ├── merge_ml_table.py   # Join labels + features
│   └── train.py            # Train and evaluate models
│
├── results/
│   ├── figures/            # ROC, PR, calibration plots
│   └── metrics.json
│
└── README.md
```



##  Models



 1️. Logistic Regression 

* Features: Grantham distance + protein position
* Purpose: establish a simple, interpretable baseline

 2️. XGBoost (nonlinear)

* Features: all four protein-level features
* Captures nonlinear interactions (e.g. conserved + drastic change)

(Class imbalance is handled via weighting)


##  Evaluation Metrics

Because pathogenic variants are rare, reported:

* **ROC-AUC** (overall ranking ability)
* **Precision–Recall AUC** (more informative for imbalance)
* **Calibration curves** (probability reliability)

---

##  Results

### Model performance (missense-only dataset)

| Model               | ROC-AUC   | PR-AUC    |
| ------------------- | --------- | --------- |
| Logistic regression | ~0.77     | ~0.36     |
| **XGBoost**         | **~0.94** | **~0.79** |

The nonlinear model provides a large performance gain, indicating strong feature interactions.


### Feature importance (XGBoost)

| Feature            | Importance |
| ------------------ | ---------- |
| **SIFT score**     | ~0.63      |
| **PolyPhen score** | ~0.17      |
| Protein position   | ~0.12      |
| Grantham distance  | ~0.08      |

#### Interpretation

* Conservation-based predictors dominate, as expected biologically
* Protein position adds contextual signal
* Physicochemical severity alone is informative but insufficient

This ordering aligns with domain knowledge and increases confidence in the model.



##  Limitations

* Conservation features are derived from existing predictors, not raw MSAs
* Structural features (like, solvent accessibility) are not yet included
* Only missense variants are modeled
* ClinVar labels may contain noise


##  References

1. Landrum et al., *ClinVar: improving access to variant interpretations*, Nucleic Acids Research
2. McLaren et al., *The Ensembl Variant Effect Predictor*, Genome Biology
3. Grantham, *Amino acid difference formula*, Science (1974)
4. Ng & Henikoff, *SIFT*, Nucleic Acids Research
5. Adzhubei et al., *PolyPhen-2*, Nature Methods
6. Richards et al., *ACMG guidelines for variant interpretation*, Genetics in Medicine



Any feedback, please  reach out to Suprajaa V | suprajaav4@gmail.com


