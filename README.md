---

# 🧬 Phenomic-to-GWAS Pipeline: Spatial Correction, Clustering & Core Subset Sampling (v5.4)

An end-to-end R pipeline for **field-scale phenomic preprocessing, macro-phenotypic structure discovery, and GWAS-oriented core subset sampling**.

Designed for **selfed family populations**, this framework integrates **2D spatial correction (P-splines), multivariate clustering, PCA-based gradient extraction, and variance-weighted Neyman allocation** to construct a statistically representative and biologically meaningful subset for Genome-Wide Association Studies (GWAS).

---

## 🚀 What's New in v5.4

### 🧭 1. Pre-GWAS Phenomic Module (NEW)

A dedicated upstream module has been introduced to **remove environmental noise and define phenotypic structure before sampling**:

* **2D Spatial Correction (P-splines via `sommer`)**
  Corrects field heterogeneity (row × column effects) and extracts **noise-reduced phenotypes (BLUP-like estimates)**.

* **Matrix-to-Coordinate Parsing**
  Converts field layout identifiers (e.g., A01, AA47) into **strict numerical spatial coordinates**, enabling spatial modeling.

* **Macro Phenotypic Clustering (K-means + Silhouette)**
  Identifies **population-level phenotypic strata**, with biologically constrained cluster number (K ≤ 3).

* **PCA Structural Visualization**
  Generates **low-dimensional representations of phenotypic architecture**, supporting downstream sampling justification.

* **Plot Data Export (Origin/GraphPad-ready)**
  Automatically exports **underlying plotting datasets**, ensuring full reproducibility and publication flexibility.

---

### 📊 2. Enhanced Sampling Validation

* **K-S Test (v5.3 retained)**
  Quantifies distributional similarity between total population and sampled subset.
* Now explicitly grounded in:

  > *“sampling from spatially corrected and cluster-informed phenotypic space”*

---

### ⚙️ 3. Unified Two-Stage Design

The pipeline is now explicitly structured as:

```
Stage I   : Spatial Correction + Phenotypic Clustering
Stage II  : PCA-driven Stratified Sampling (GWAS subset construction)
```

This separation aligns directly with **standard GWAS methodological logic**:

> noise removal → structure identification → representative sampling

---

## 🛠️ Key Features

### 🔬 Stage I: Phenomic Preprocessing

* [x] **Spatial Noise Removal:** 2D spline-based correction of field effects
* [x] **Robust LMM Fallback:** Automatic downgrade to standard LMM if spatial model fails
* [x] **Phenotypic Scaling & Cleaning:** Z-score normalization with NA handling
* [x] **Cluster-Constrained Stratification:** Avoids over-fragmentation of family-level variation

### 📈 Stage II: GWAS-Oriented Sampling

* [x] **Multi-Year Data Fusion:** Integrates 2025–2026 datasets and growth dynamics
* [x] **Repeatability Estimation (LMM):** Quantifies genetic signal strength
* [x] **Stepwise PCA Framework:** Extracts core phenotypic gradients (PC1 as growth axis)
* [x] **Capped Neyman Allocation:** Variance-weighted optimal sampling with family constraints
* [x] **Stratified Quantile Sampling:** Ensures coverage of extreme and median phenotypes

---

## 📦 Dependencies

Automatically installed if missing:

* **Data Wrangling:** `dplyr`, `tidyr`, `readr`
* **Modeling:** `sommer`, `lme4`, `lmerTest`
* **Multivariate Analysis:** `factoextra`, `cluster`
* **Visualization:** `ggplot2`, `patchwork`, `showtext`

> **⚠️ Font Note**
> The pipeline uses `SimSun` for publication-grade figures.
> On macOS/Linux, ensure availability of `STSong.ttf` or equivalent CJK font.

---

## 🚦 Usage

### 1. Input Data

Place the following files in your working directory:

* `data_2026.csv` → **Required for Stage I (spatial correction + clustering)**
* `data_2025.csv` → Optional but required for **two-year GWAS sampling**

---

### 2. Run Stage I (Pre-GWAS)

```r
source("PreGWAS_Phenomic_Pipeline.R")
```

Outputs:

* `data_2026_Corrected.csv`
* `data_2026_Corrected_with_Clusters.csv`
* Phenotypic diversity plots + plotting datasets

---

### 3. Configure Sampling Parameters

```r
TARGET_TOTAL_N <- 140
MAX_PER_FAMILY <- 20
```

---

### 4. Run Stage II (Sampling Pipeline)

```r
source("GWAS_Sampling_Pipeline_v5.4.R")
```

---

## 📂 Output Structure

### 📈 Phenomic Structure (Stage I)

* `01_Phenotypic_Diversity_Structure.pdf`

  * A: Optimal cluster number (Silhouette)
  * B: PCA-based phenotypic structure

* `PlotData_FigA_Silhouette_Scores.csv`

* `PlotData_FigB_PCA_Scatter.csv`

---

### 📈 GWAS Sampling Justification (Stage II)

* `Fig1_Segregation_Analysis.pdf`
* `Fig2_PCA_Stepwise_Biplots.pdf`
* `Fig3_Sampling_Justification.pdf`

  * A: Neyman allocation
  * B: K-S test distribution comparison
  * C: Within-family gradient coverage

---

### 📊 Core Data Outputs

* `data_2026_Corrected.csv` ⭐（关键输入）
* `03_Comprehensive_PCA_Scores.csv`
* **`04_Selected_140_GWAS.csv`** ⭐（最终GWAS样本）

---

## 🧠 Methodological Insight

### Why Spatial Correction First?

Field experiments introduce structured environmental variance:

* irrigation gradients
* soil heterogeneity
* row/column effects

Failing to remove these leads to:

> **phenotype ≠ genotype signal**

This pipeline explicitly models:

```
Phenotype = Genetic Effect + Spatial Noise
```

and extracts:

> **noise-reduced phenotypes for unbiased downstream analysis**

---

### Why Cluster Before Sampling?

Traditional GWAS sampling assumes random or variance-based selection.
This pipeline instead operates on:

> **cluster-informed phenotypic space**

Advantages:

* prevents over-sampling dominant families
* preserves **macro-level trait architecture**
* improves **allelic diversity capture**

---

## 💡 Interpreting the K-S Test (Fig 3B)

* **D-statistic** → maximum distributional deviation
* **p-value** → statistical difference

> A small D indicates strong representativeness.

⚠️ **Important:**
A significant p-value (*p < 0.05*) is expected.

Reason:

* Sampling intentionally enriches **phenotypic extremes**
* This increases GWAS detection power

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

---
