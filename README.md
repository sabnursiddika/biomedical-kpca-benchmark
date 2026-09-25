# Biomedical KPCA Benchmark

This repository contains the code, final results, figures, and supplementary materials for the study:

**When Do Feature Representations Pay Off in Biomedical Classification? A Leakage-Controlled Repeated Nested Cross-Validation Study**

## Study overview

The study evaluates whether alternative feature representations improve biomedical binary classification relative to standardized raw predictors.

Six datasets are included:

- Wisconsin
- Pima
- Colon
- Gravier
- Prostate
- Leukemia

Four predictor representations are compared:

1. Raw standardized predictors
2. Principal component analysis (PCA)
3. Global RBF kernel PCA
4. Locally scaled kernel PCA

Each representation is evaluated with:

- Linear SVM
- RBF SVM

The primary performance metric is outer-test ROC-AUC.

## Validation design

The benchmark uses repeated nested cross-validation:

- 5 outer folds
- 3 outer repeats
- 3 inner folds
- 15 outer evaluations per method
- 8 representation-classifier combinations per dataset
- 6 datasets
- 720 total outer-test evaluations

All preprocessing, imputation, scaling, feature-representation fitting, hyperparameter tuning, and classification-threshold selection are performed using training data only.

## Repository structure

```text
biomedical-kpca-benchmark/
├── README.md
├── sessionInfo.txt
├── code/
│   ├── 01_main_benchmark.R
│   └── 02_final_figures.R
├── data/
│   └── README.md
├── results/
│   └── benchmark_results_final/
└── supplementary/
    ├── Supplementary_Table_S1.csv
    └── Supplementary_Table_S2.csv