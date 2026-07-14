# Integrated Pre-GWAS Phenomic Pipeline & Core Collection Sampling (v7.2)

**Version:** 7.2
**Language:** R (≥ 4.3)  
**Purpose:** Phenomic preprocessing, core collection optimization, and statistical justification for whole-genome resequencing prior to GWAS.

---

## Overview

This pipeline performs an integrated pre-GWAS phenomic analysis and optimized core collection sampling for advanced backcross-derived populations.

The workflow combines:

- phenotypic preprocessing and robust missing-value imputation,
- phenotypic diversity assessment,
- PCA-based multivariate characterization,
- weighted Neyman allocation,
- stratified gradient sampling,
- exact sample-size enforcement through force-fill,
- statistical validation of sampling representativeness, and
- automatic generation of publication-quality figures and supplementary tables.

The pipeline was designed for projects requiring **exactly 200 individuals** for whole-genome resequencing while maximizing phenotypic diversity and minimizing family-structure bias.

---

# Workflow

```
Raw Phenotype Data
        │
        ▼
Stage 1
Data cleaning
Family-wise median imputation
Small-family filtering
Phenotypic clustering
        │
        ▼
Stage 2
Variance decomposition
Repeatability estimation
Principal component analysis
        │
        ▼
Weighted Neyman Allocation
        │
        ▼
Stratified Gradient Sampling
        │
        ▼
Phase 4 Force-Fill
(exact N = 200)
        │
        ▼
Sampling validation
(K-S test, cluster coverage,
allocation balance)
        │
        ▼
Publication-ready figures
+
Supplementary data tables
```

---

# Major Features

## 1. Robust Phenotypic Preprocessing

The pipeline performs hierarchical median imputation.

Missing observations are first imputed using the **within-family median**, followed by global median imputation if necessary.

This strategy preserves family-level phenotypic structure while avoiding unnecessary removal of individuals.

---

## 2. Phenotypic Diversity Assessment

Seven quantitative traits are jointly analyzed.

Current implementation includes:

- Plant height (Nov)
- Internode length (Nov)
- Plant height (Mar)
- Internode length (Mar)
- Branch number (Mar)
- Multifoliate score (Mar)
- Plant height (May)

Phenotypic diversity is evaluated using

- K-means clustering
- Silhouette analysis
- Principal component analysis (PCA)

---

## 3. Weighted Neyman Allocation

Family sampling sizes are determined using weighted Neyman allocation.

Sampling weight for each family was calculated as:

`Weight_h = N_h × S_h`

where `N_h` is the family size and `S_h` is the weighted phenotypic standard deviation estimated from the first two principal components.

The weighted variance was computed as:

`Var_h = Var(PC1) × PVE_PC1 + Var(PC2) × PVE_PC2`

This strategy preferentially allocates sequencing resources toward genetically informative families while maintaining representation across the population.

---

## 4. Stratified Gradient Sampling

Within each family, individuals are sampled using PC1 gradients.

The PC1 axis is divided into quantiles.

Random sampling is then performed independently within each gradient.

Compared with purely random sampling, this strategy preserves phenotypic diversity inside each family.

---

## 5. Exact Sample Size Enforcement

One practical limitation of constrained Neyman allocation is that integer rounding frequently produces fewer than the required number of samples.

Version 7.2 introduces a dedicated **Phase 4 Force-Fill** algorithm.

If

```
Selected < Target
```

the algorithm automatically

- identifies families with remaining individuals,
- prioritizes families according to Neyman weights,
- allows limited soft-cap relaxation,
- sequentially fills missing samples,

until

```
Selected == Target
```

is satisfied.

This guarantees the supervisor-required sequencing size of exactly **200 individuals**.

---

## 6. Statistical Validation

Sampling quality is evaluated using multiple independent metrics.

### Distribution preservation

Kolmogorov–Smirnov test

compares the PC1 distributions between

- complete population
- selected core collection

A non-significant result indicates representative sampling.

---

### Allocation balance

The Gini coefficient is calculated to quantify sampling inequality among families.

Lower values indicate more balanced sampling.

---

### Cluster coverage

Sampling ratios are calculated for every phenotypic cluster to verify that all major phenotypic groups remain represented.

---

# Output Structure

```
GWAS_Pipeline_Result_v7.2_YYYYMMDD_HHMM/

│
├──01_Scripts/
│
├──02_Data_Tables/
│     Imputed datasets
│     PCA results
│     Repeatability estimates
│     Core collection
│     Diagnostic reports
│
├──03_Figures/
│     Fig1
│     Fig2
│     Fig3
│     Supplementary figures
│
└──04_Plot_Data/
      Source data
      Plot-ready tables
      Force-fill logs
```

Each execution creates a timestamped directory, ensuring complete reproducibility without overwriting previous analyses.

---

# Publication Figures

The pipeline automatically generates publication-quality vector graphics.

Current outputs include

| Figure | Description |
|---------|-------------|
| Fig S1 | Phenotypic clustering |
| Fig 1 | Segregation analysis and repeatability |
| Fig 2 | PCA biplot |
| Fig 3 | Sampling justification and validation |

All figures

- use Arial fonts,
- adopt Nature-compatible formatting,
- use colorblind-friendly palettes,
- are exported as vector PDF via Cairo.

---

# Key Parameters

| Parameter | Default | Description |
|------------|---------|-------------|
| TARGET_TOTAL_N | 200 | Exact number of resequencing samples |
| MAX_PER_FAMILY | 7 | Soft cap during Neyman allocation |
| FILL_ALLOW_EXCESS | 2 | Maximum temporary cap increase during force-fill |
| MIN_FAMILY_SIZE | 3 | Minimum eligible family size |
| MIN_ALLOC | 1 | Minimum allocation for every eligible family |

---

# Required Packages

```
dplyr
tidyr
ggplot2
lme4
lmerTest
factoextra
readr
patchwork
cluster
viridis
ggsci
scales
RColorBrewer
```

All required packages are automatically installed if absent.

---

# Input Data

The pipeline expects one phenotype table

```
01_Raw_Data/
    data_202605_gwas.csv
```

Required variables include

- Family
- ID
- PH_Nov
- IN_Nov
- PH_Mar
- IN_Mar
- BN_Mar
- MF_Score_Mar
- PH_May

These variables are automatically renamed internally to standardized analysis names.

---

# Design Principles

Version 7.2 emphasizes three design objectives.

1. **Statistical efficiency**

Maximize GWAS power through variance-weighted sampling.

2. **Population representativeness**

Maintain phenotypic diversity and family representation.

3. **Reproducibility**

Generate all figures, tables, diagnostics, and logs automatically from raw phenotypic data.

---

# Version History

## v7.2

- Exact 200-individual guarantee
- Phase 4 Force-Fill algorithm
- Nature-style publication formatting
- Arial-only typography
- Colorblind-friendly palettes
- Cairo PDF export
- Timestamped output structure
- Complete diagnostic reporting

## v7.1

- Exact target sampling
- Balanced family allocation
- Minimum family contribution
- Force-fill framework

## v7.0

- Initial integrated pre-GWAS phenomic pipeline
- Neyman allocation
- PCA-based sampling
- Publication-ready visualization

---

# Citation

If this pipeline contributes to published work, please cite the associated manuscript describing the sampling strategy and phenomic optimization framework.

```
Integrated Pre-GWAS Phenomic Pipeline (v7.2)
Advanced Backcross Population Sampling for Whole-Genome Association Analysis.
```

---

# License

For academic research use.

Copyright © Your Laboratory.
All rights reserved.
