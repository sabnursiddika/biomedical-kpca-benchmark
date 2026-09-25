# Data

This benchmark uses six publicly available biomedical binary-classification datasets.

Two datasets are supplied to the benchmark as local input files:

- `wdbc.data`
- `pima.csv`

The remaining four datasets are loaded directly from R packages.

## Wisconsin Breast Cancer Diagnostic

Expected filename:

`wdbc.data`

The benchmark expects:

- 569 observations
- 32 columns
- column 1: ID
- column 2: diagnosis
- columns 3–32: 30 numeric predictors

Diagnosis coding in the original data:

- `B` = benign
- `M` = malignant

For the benchmark:

- `0` = benign
- `1` = malignant

Source:

UCI Machine Learning Repository, Breast Cancer Wisconsin (Diagnostic).

DOI: `10.24432/C5DW2B`

## Pima Indians Diabetes

Expected filename:

`pima.csv`

The benchmark expects:

- 768 observations
- 9 columns
- no header row

Column order:

1. Pregnancies
2. Glucose
3. BloodPressure
4. SkinThickness
5. Insulin
6. BMI
7. DiabetesPedigree
8. Age
9. Outcome

Outcome coding:

- `0` = negative
- `1` = positive

Before modeling, zero values in the following variables are treated as missing:

- Glucose
- BloodPressure
- SkinThickness
- Insulin
- BMI

Missing-value imputation is performed using quantities estimated from the corresponding training data only.

## Package-provided datasets

The remaining datasets are loaded automatically by `code/01_main_benchmark.R`:

- Colon — `plsgenomics` package
- Gravier — `datamicroarray` package
- Prostate — `sda` package
- Leukemia — `SIS` package

For Gravier, the benchmark codes:

- `0` = good prognosis
- `1` = poor prognosis

## Raw data files

The local copies of `wdbc.data` and `pima.csv` are excluded from Git tracking by `.gitignore`.

Place both files in this directory before intentionally running the full benchmark:

```text
data/
├── README.md
├── wdbc.data
└── pima.csv