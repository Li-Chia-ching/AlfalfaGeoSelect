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

This pipeline implements a **multi-year, multi-trait selection framework** for identifying elite alfalfa individuals suitable for GWAS.

It integrates phenotypic data from:
- **Autumn 2025 (establishment stage)**
- **Spring 2026 (regrowth stage)**

The pipeline combines:
- **feature engineering (dynamic traits)**
- **Pareto optimality**
- **MGIDI-based ideotype selection**

to produce a robust set of candidate genotypes.

---

## Key Innovations

### 1. Dynamic Trait Engineering
Derived traits capturing temporal response:

- `Growth_Diff = Spring - Autumn`
- `Plasticity_Ratio = Spring / Autumn`

These traits quantify:
- compensatory growth
- seasonal plasticity

---

### 2. Reaction Norm Analysis
- Tracks **individual-level trajectories**
- Reveals:
  - stable genotypes
  - high plasticity genotypes
  - autumn-dormant vs spring-responsive types

---

### 3. Trait Correlation Network
- Constructs phenotype network using Pearson correlations
- Identifies **hub traits** influencing overall architecture
- Enables biological interpretation beyond pairwise correlations

---

### 4. Dual Selection Strategy

#### Pareto Front
- Non-dominated selection across multiple traits
- No weighting assumptions

#### MGIDI (Ideotype Distance)
- Selects individuals closest to the theoretical optimum

#### Final Selection:
```

Final Candidates = Pareto ∩ MGIDI

```

---

## Input Data

### Autumn Dataset (2025)
Required columns:
- `Parent`
- `SampleID`
- `PlantHeight`
- `NodeCount`
- `Group`

### Spring Dataset (2026)
Required columns:
- `Group`
- `Matrix`
- `Plant_Height`
- `Internode`
- `Branch_Number`
- `Multifoliate_Score`
- `Death_Code`

---

## Output Structure

```

gwas_selection_YYYY-MM-DD/
│
├── data/
│   ├── merged_data.csv
│   ├── correlation_matrix.csv
│   └── gwas_candidates.csv
│
├── plots/
│   ├── reaction_norms.png
│   └── trait_network.png

```

---

## Workflow

1. Data merging across years
2. Feature engineering (dynamic traits)
3. Survival filtering
4. Reaction norm visualization
5. Correlation network analysis
6. Pareto selection
7. MGIDI ranking
8. Final candidate extraction

---

## Use Case

This pipeline is designed for:

- GWAS population construction
- Multi-trait elite selection
- Phenotypic data mining across environments

---

## Advantages

- Integrates **temporal dynamics**
- Avoids arbitrary weighting (Pareto)
- Incorporates **ideotype theory (MGIDI)**
- Produces reproducible and structured outputs

---

## Limitations

- MGIDI implemented as Euclidean approximation (not full factor-based model)
- No mixed-model (BLUP) correction
- Environmental variance not explicitly modeled

---

## Recommended Extensions

- Replace MGIDI with `metan::mgidi()` implementation
- Incorporate mixed models (e.g., `lme4`, `sommer`)
- Use derived traits (plasticity) as GWAS phenotypes

---

## Relationship to GeoSelect

| Module | Purpose |
|------|--------|
| GeoSelect | Single-season geometric selection |
| Multi-Year Pipeline | Cross-season, GWAS-ready selection |

---

## Author Notes

This pipeline is designed to bridge **field phenotyping** and **genomic analysis**, providing a scalable and interpretable framework for modern alfalfa breeding.

---
