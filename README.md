# AlfalfaGeoSelect: Geometric Selection and Visualization for Phenotypic Screening

## Overview
**AlfalfaGeoSelect** is a lightweight phenotypic selection framework designed for rapid identification of elite individuals based on multi-trait geometric integration.

This module implements a **weighted geometric selection index** combining plant height and multifoliate score, accompanied by a **publication-ready contour visualization** to illustrate the selection boundary in trait space.

---

## Key Features

- **Geometric Selection Index**
  - Combines traits using a multiplicative model:
    ```
    Index = (0.5 × Height) × (0.5 × Multifoliate Score)
    ```
  - Emphasizes balanced performance across traits rather than single-trait dominance

- **Top-N Selection Strategy**
  - Automatically selects top-performing individuals (default: Top 50)
  - Defines an explicit **selection cutoff boundary**

- **Professional Visualization**
  - Contour-based selection landscape
  - Red dashed cutoff line for threshold interpretation
  - Jittered points to resolve overlap
  - Vector arrow indicating **selection gradient**

- **Flexible Data Input**
  - Supports both:
    - In-memory objects
    - External CSV files

---

## Input Data Format

### 1. Plant Height
- Matrix-like format (rows × columns)
- Example:
```

L1, L2, ...

```

### 2. Multifoliate Score
- Same structure as height
- Score range: typically 1–5

---

## Output

The script generates a timestamped directory:

```

Alfalfa_GeoViz_Pro_YYYYMMDD_HHMMSS/

```

Contents:

- `Plot_Geometric_Selection_Pro.png`
  - Publication-ready contour plot
- `01_Selection_List.csv`
  - Selected top individuals

---

## Methodological Notes

- The geometric index is equivalent to a **multiplicative selection surface**, which:
  - Penalizes imbalance between traits
  - Favors individuals performing well across all dimensions

- The contour plot represents:
  - Iso-selection lines (equal index values)
  - A continuous **fitness landscape**

---

## Use Case

This module is best suited for:

- Early-stage phenotypic screening
- Field trial visualization
- Rapid selection prior to genetic analysis

---

## Limitations

- Does not account for:
  - survival
  - environmental variation
  - multi-year dynamics

- Intended as a **first-pass selection tool**, not a final GWAS pipeline

---

## Related Module

For advanced multi-year selection and GWAS-ready candidate identification, see:
```
👉 `alfalfa-multiyear-gwas-pipeline.R`
```

---

# Alfalfa Multi-Year GWAS Pipeline: Dynamic Trait Integration and Elite Selection

## Overview

This repository contains a publication-ready R pipeline for:

* Genetic variation analysis in structured populations
* Variance component estimation using linear mixed models (LMM)
* Multivariate phenotypic structure analysis (PCA)
* Neyman optimal allocation for stratified sampling
* Quantile-based deterministic sampling within families
* Extraction of a representative core subset (n = 200) for GWAS

The pipeline is specifically designed for **two-year phenotypic datasets measured under the same environment**, with a focus on **family-based segregation populations**.

---

## Key Features

* **Repeatability estimation (R)** using variance components (family vs residual)
* **Multivariate sampling strategy** based on PC1–PC2 dispersion
* **Robust stratified sampling** combining:

  * Neyman allocation (between-family)
  * Quantile gradient sampling (within-family)
* **Deterministic + stochastic hybrid sampling** ensuring representativeness
* **Automatic export of all figure-level datasets** for reproducibility
* **Publication-quality visualization outputs (PDF, 600 dpi)**

---

## Input Data Requirements

### File 1: `data_2025.csv`

Required columns:

| Column Name | Description               |
| ----------- | ------------------------- |
| Parent      | Family identifier         |
| SampleID    | Individual plant ID       |
| Real_Count  | (Optional) count metadata |
| PlantHeight | Plant height (2025)       |
| NodeCount   | Internode number (2025)   |

### File 2: `data_2026.csv`

Required columns:

| Column Name        | Description              |
| ------------------ | ------------------------ |
| Group              | Family identifier        |
| Matrix             | Individual plant ID      |
| Plant_Height       | Plant height (2026)      |
| Internode          | Internode number         |
| Branch_Number      | Branch number            |
| Multifoliate_Score | Multifoliate trait score |
| Death_Code         | Survival status          |

### Notes

* `Family_ID` and `Plant_ID` must match across years for merging
* Missing or invalid family labels in 2026 are automatically filtered
* Only individuals labeled as **Alive** are retained for downstream analysis

---

## Output Structure

All results are written to:

```
GWAS_Sampling_Results_YYYYMMDD/
```

### Core Data Outputs

| File                                     | Description                           |
| ---------------------------------------- | ------------------------------------- |
| 01_Merged_Wide_Data_All.csv              | Merged two-year dataset               |
| 02_Variance_Components_Repeatability.csv | Variance components and repeatability |
| 03_Selected_200_GWAS.csv                 | Final selected GWAS subset            |

---

## Figure Outputs and Corresponding Data

### Figure 1 – Phenotypic Segregation

| File                          | Description       |
| ----------------------------- | ----------------- |
| Fig1_Segregation_Analysis.pdf | Density + boxplot |
| Fig1A_Data.csv                | Density plot data |
| Fig1B_Data.csv                | Boxplot data      |

---

### Figure 2 – PCA Structure

| File                   | Description                        |
| ---------------------- | ---------------------------------- |
| Fig2_PCA_Structure.pdf | Scree plot + PCA distribution      |
| Fig2A_Data.csv         | Eigenvalues and variance explained |
| Fig2B_Data.csv         | Individual PCA scores              |

---

### Figure 3 – Sampling Justification

| File                                 | Description                                   |
| ------------------------------------ | --------------------------------------------- |
| Fig3_Sampling_Justification.pdf      | Allocation + density + strip plot             |
| Fig3A_Data.csv                       | Neyman allocation results                     |
| Fig3B_Data.csv                       | Density comparison (before vs after sampling) |
| Fig3C_Data.csv                       | Within-family sampling distribution           |
| Fig3_Selected_200_Plants_Details.csv | Selected individuals with traits              |

---

## Methodological Workflow

### Step 1: Data Integration

* Merge 2025 and 2026 datasets using `Family_ID` and `Plant_ID`
* Define survival status
* Compute derived traits (e.g., ΔHeight)

### Step 2: Phenotypic Segregation Analysis

* Visualize within-family distributions (density + rug)
* Compare between-family variation (boxplot)

### Step 3: Variance Decomposition

Model:

```
y = μ + Family + ε
```

* Estimate:

  * Between-family variance
  * Within-family variance
* Compute repeatability:

```
R = σ²_family / (σ²_family + σ²_residual)
```

---

### Step 4: Multivariate Structure (PCA)

* Standardize traits
* Extract PC1 and PC2
* Use PCA space as sampling coordinate system

---

### Step 5: Neyman Optimal Allocation

For each family:

```
Weight_h = N_h × S_h
S_h = sqrt(Var(PC1) + Var(PC2))
```

* Allocate sample sizes proportionally
* Apply integer correction (`smart_round`) to ensure total = 200

---

### Step 6: Within-Family Sampling

* Divide individuals into **quantile-based PC1 gradients** (≤3 strata)
* Allocate samples proportionally within strata
* Perform random sampling within each gradient
* Fill shortfall if necessary

---

### Step 7: Representativeness Validation

* Compare PC1 density distributions (full vs sampled)
* Visualize coverage across families
* Confirm sampling spans full phenotypic gradients

---

## Reproducibility

* Random seed fixed: `set.seed(2026)`
* All intermediate datasets exported
* Figures generated directly from saved data

---

## Dependencies

Required R packages:

```
dplyr, tidyr, ggplot2, lme4, lmerTest,
factoextra, readr, patchwork, showtext
```

The script automatically installs missing packages.

---

## Recommended Citation (Methods Description)

This pipeline implements a **multivariate stratified sampling framework combining Neyman allocation and PCA-based quantile sampling**, designed to maximize genetic representativeness in GWAS subset selection from structured family populations.

---

## Contact / Notes

* Ensure consistent trait units across years
* Suitable for selfing-derived populations (e.g., EMS, inbred families)
* Easily extendable to additional traits or environments

---

## License

For academic and research use.
