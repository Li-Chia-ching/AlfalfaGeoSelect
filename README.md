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

A reproducible R pipeline for **genetic variation analysis** and **GWAS core subset selection**, integrating:

- Multi-year phenotypic data
- Linear mixed models (LMM)
- Stepwise PCA decomposition
- Neyman optimal allocation
- Quantile-based stratified sampling

---

## 🚀 Key Features

- 📊 **Two-year phenotype integration (2025 + 2026)**
- 🧠 **Stepwise PCA framework**
  - Cross-year shared traits
  - Single-year expanded traits
- 📈 **Repeatability estimation (LMM)**
- 🎯 **Neyman optimal allocation**
- 🔬 **Deterministic within-family sampling**
- 📦 **Full export of publication-ready figures and data**

---

## 📂 Input Data Requirements

### `data_2025.csv`
| Column        | Description              |
|--------------|--------------------------|
| Parent       | Family ID                |
| SampleID     | Plant ID                 |
| PlantHeight  | Height                   |
| NodeCount    | Internode number         |

---

### `data_2026.csv`
| Column             | Description              |
|--------------------|--------------------------|
| Group              | Family ID                |
| Matrix             | Plant ID                 |
| Plant_Height       | Height                   |
| Internode          | Internode number         |
| Branch_Number      | Branch count             |
| Multifoliate_Score | Leaf complexity score    |
| Death_Code         | Survival status          |

---

## ⚙️ Installation

```r
install.packages(c(
  "dplyr", "tidyr", "ggplot2", "lme4", "lmerTest",
  "factoextra", "readr", "patchwork", "showtext"
))
````

---

## ▶️ Usage

```r
source("GWAS_sampling_v3.1.R")
```

Ensure input files are in the working directory:

```
data_2025.csv
data_2026.csv
```

---

## 🧬 Workflow Overview

### 1. Data Processing

* Merge 2025 and 2026 datasets
* Filter alive individuals
* Remove missing values

---

### 2. Phenotypic Analysis

* Density distribution (Fig1A)
* Boxplot comparison (Fig1B)
* LMM-based repeatability

---

### 3. PCA Analysis (Stepwise)

#### A. Cross-Year PCA

* Traits:

  * Height (2025 & 2026)
  * Internode (2025 & 2026)

#### B. 2026 PCA

* Traits:

  * Height
  * Internode
  * Branch
  * Multifoliate

#### C. Comprehensive PCA

* All 6 traits
* Used for:

  * PC1 extraction
  * Sampling

---

### 4. Sampling Strategy

#### Neyman Allocation

$$
n_h \propto N_h \cdot S_h
$$

Where:

* $S_h = \sqrt{\mathrm{Var}(PC1) + \mathrm{Var}(PC2)}$

#### Within-Family Sampling

* Stratified by **PC1 quantiles**
* Ensures gradient coverage:

  * High
  * Medium
  * Low

---

### 5. Output Generation

#### 📁 Main Outputs

```
01_Merged_Wide_Data_All.csv
02_Variance_Components_Repeatability.csv
03_Comprehensive_PCA_*.csv
04_Selected_200_GWAS.csv
```

#### 📊 Figures

* Fig1: Segregation analysis
* Fig2: PCA biplots
* Fig3: Sampling validation

All figures include:

* Source data (CSV)
* High-resolution PDF

---

## 📊 Output Directory Structure

```
GWAS_Sampling_Two-Year_YYYYMMDD/
├── Data tables (CSV)
├── PCA results
├── Sampling results
├── Figures (PDF)
```

---

## ⚠️ Notes

* Only **alive individuals** are analyzed
* Missing values in traits → removed
* PCA uses **standardization (scale = TRUE)**
* Sampling size fixed at **200 individuals**

---

## 🔧 Customization

You can modify:

```r
# Target sample size
smart_round(..., target = 200)

# Number of quantile groups
k_clusters <- 3
```

---

## 📌 Applications

* GWAS population design
* Core germplasm extraction
* Multi-environment trait integration
* Forage breeding research

---

## 📜 License

**MIT License**

---

## 👤 Author

Developed for advanced research in:

* Forage breeding
* Legume genetics
* Germplasm evaluation

---
