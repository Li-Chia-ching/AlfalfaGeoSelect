
# 🧬 Integrated Pre-GWAS Phenomic Pipeline & Core Collection Sampling (v6.0)

An end-to-end R pipeline for **phenomic structure dissection and GWAS-oriented core subset selection** in *Medicago sativa* (alfalfa) or similar outcrossing populations.

This version integrates:

* **Multidimensional phenotypic structure analysis (PCA + clustering)**
* **Variance-aware Neyman allocation with stratified sampling**

while **removing spatial correction modules** to prevent statistical artifacts.

---

## 🚀 Key Features

### 1. Phenomic Structure Decomposition (Macro-scale)

* Robust **missing data imputation** (within-family median → global fallback)
* **K-means clustering** with silhouette-based K selection
* **PCA-based gradient analysis**
* Visualization:

  * Phenotypic clustering (Fig S1)
  * PCA biplot with trait loadings (Fig 2)

---

### 2. Genetic Signal Validation (Pre-GWAS Evidence)

* Within-family segregation visualization (density + boxplot)
* Variance component estimation using **LMM (lme4)**
* Repeatability calculation:  
  `R = σ²_family / (σ²_family + σ²_within)`

---

### 3. Core Collection Sampling (Micro-scale)

* **Weighted Neyman allocation**:

  `n_h ∝ N_h · S_h`

* Variance defined in **PCA space (PC1 + PC2 weighted)**

* Constraints:

  * Total sample size (`TARGET_TOTAL_N`)
  * Max per family (`MAX_PER_FAMILY`)

* Within-family sampling:

  * **PC1 quantile stratification**
  * Avoids extreme-value bias

---

### 4. Sampling Validation

* Distribution comparison (K-S test)
* Density overlap (population vs subset)
* Within-family sampling strip plot
* Cluster coverage summary

---

## ⚙️ Pipeline Structure

```
Stage 1:
  ├── Data loading & imputation
  ├── Phenotypic clustering (K-means)
  └── PCA structure visualization

Stage 2:
  ├── Segregation & variance validation
  ├── PCA gradient decomposition
  ├── Neyman allocation
  ├── Stratified sampling
  └── Sampling validation & coverage
```

---

## 📂 Input Data Requirements

### Required file:

```
data_202605.csv
```

### Required columns:

* `Family` → will be renamed to `Family_ID`
* `ID` → will be renamed to `Plant_ID`

### Trait variables:

```
Plant_Height_Nov
Internode_Nov
Plant_Height_Mar
Internode_Mar
Branch_Number_Mar
Multifoliate_Score_Mar
Plant_Height_May
```

---

## 📊 Output Structure

All outputs are written to:

```
GWAS_Integrated_Pipeline_YYYYMMDD/
```

### Core outputs:

| File                                       | Description                       |
| ------------------------------------------ | --------------------------------- |
| `00_Imputed_Dataset.csv`                   | Cleaned dataset                   |
| `01_Variance_Components_Repeatability.csv` | Variance decomposition            |
| `02_PCA_*`                                 | PCA eigenvalues, loadings, scores |
| `03_Selected_XXX_GWAS.csv`                 | Final selected individuals        |
| `03_Selected_XXX_WithTraits.csv`           | Selection with phenotypes         |
| `04_Cluster_Coverage_Summary.csv`          | Sampling representativeness       |

---

### Figures (publication-ready)

* `FigS1_Phenotypic_Clustering.pdf`
* `Fig1_Segregation_Analysis.pdf`
* `Fig2_PCA_Biplot_Enhanced.pdf`
* `Fig3_Sampling_Justification.pdf`

---

## 🔬 Methodological Highlights

### Why remove spatial correction?

Previous versions included spatial adjustment, but:

* Caused **overfitting**
* Led to **within-family variance collapse (SD → 0)**
* Invalidated downstream:

  * variance component estimation
  * Neyman allocation

👉 v6.0 uses **raw phenotypic structure + robust imputation**, ensuring:

* biologically meaningful variance
* stable sampling weights

---

### Why PCA-driven sampling?

* Captures **multivariate trait covariance**
* Reduces dimensionality
* Provides a **continuous sampling gradient (PC1)**

---

### Why Neyman allocation?

* Maximizes **sampling efficiency**
* Allocates more samples to:

  * larger families
  * more variable families

---

## 🧩 Dependencies

```r
dplyr, tidyr, ggplot2
lme4, lmerTest
factoextra
cluster
patchwork
viridis
readr
showtext
```

---

## ▶️ How to Run

```r
source("pipeline_v6.0.R")
```

Outputs will be automatically generated in a timestamped directory.

---

## 📌 Parameter Configuration

```r
TARGET_TOTAL_N <- 140   # total GWAS sample size
MAX_PER_FAMILY <- 20    # cap per family
```

---

## 🧰 Auxiliary Tool: 2D Trait Geometric Selection (V3.0)

A visualization-driven selection module for **bi-trait optimization**.

### Key Idea

Projects individuals into a **2D trait-performance surface**, enabling intuitive selection of top performers.

### Features

* Weighted geometric index
* Top-N extraction (default: 50)
* Contour-based selection boundary visualization
* Automatic field data reshaping

---

## 📌 Conceptual Workflow

```
Raw Field Data
      ↓
[Stage I]
Spatial Correction (2D Splines)
      ↓
Phenotypic Clustering (K-means)
      ↓
Macro Structure Definition
      ↓
[Stage II]
PCA Gradient Extraction
      ↓
Neyman Allocation
      ↓
Stratified Sampling
      ↓
GWAS Core Subset
```
