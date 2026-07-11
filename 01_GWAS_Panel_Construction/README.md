# GWAS Integrated Pipeline v7.1 – Exact-N Force-Fill Edition

**Version:** 7.1  
**Date:** 2026-07-11  
**Status:** Production‑ready for supervised GWAS core collection sampling  

---

## 1. Overview

This R pipeline performs integrated pre‑GWAS phenomic data processing and core collection selection. It implements a stratified Neyman allocation scheme followed by a **mandatory force‑fill phase** to guarantee a **strictly fixed total sample size** (N = 200), as required by the supervising investigator. The workflow includes data imputation, small‑family filtering, phenotypic clustering, PCA‑based gradient characterisation, variance‑component estimation (repeatability), and comprehensive diagnostic visualisations.

The pipeline is designed for self‑pollinated families (`Family_ID`) with multiple phenotypic records. It outputs a selected core set that balances family‑level genetic diversity, phenotypic coverage, and allocation equality, while respecting a per‑family soft cap to limit kinship inflation in downstream GWAS.

---

## 2. Key Changes from v6.0 → v7.1

| Feature | v6.0 | v7.1 | Rationale |
|---------|------|------|-----------|
| Target total N | 140 (flexible) | **200 (hard requirement)** | Supervisor directive; exact N is non‑negotiable |
| Per‑family soft cap | 20 | **7** | Tighter control over family‑wise kinship contribution |
| Minimum family size | none | **3** | Excludes micro‑families (<3 individuals) that cannot provide reliable within‑family variance estimates |
| Minimum allocation per family | 0 (some families may be omitted) | **1** | Guarantees representation of every eligible family, preserving rare alleles |
| Phase 4 – Force‑fill | absent | **implemented** | Detects shortfall after Neyman+stratified sampling and fills from large families up to a hard cap (soft cap + 2), ensuring exact 200 |
| Column mapping | implicit | explicit renaming of raw file fields (PH_Nov → Plant_Height_Nov, etc.) | Adapts to new data file format |
| Output directory | flat structure | hierarchical (`01_Scripts/`, `02_Data_Tables/`, `03_Figures/`, `04_Plot_Data/`) | Improved organisation for reproducibility |
| Diagnostic report | none | `05_GWASI_Diagnostic_Report.csv` with Gini, K‑S, allocation metrics | Enables post‑hoc quality review |
| Allocation algorithm | simple capped greedy | multi‑phase with floor guarantee, soft cap, and iterative fill | Reduces integer‑truncation gaps and ensures full utilisation of target |

---

## 3. System Requirements

- **R version:** ≥ 4.0.0
- **Packages:**  
  `dplyr`, `tidyr`, `ggplot2`, `lme4`, `lmerTest`, `factoextra`, `readr`, `patchwork`, `showtext`, `cluster`, `viridis`  
  *These are automatically installed if missing.*
- **Fonts:** The script attempts to load SimSun (Chinese) for PDF labels; if unavailable, falls back to default sans‑serif. No functional impact.

---

## 4. Input Data Specification

The pipeline expects a **single CSV file** located at:  
`01_Raw_Data/data_202605_gwas.csv`

**Required columns (exact names):**

| Original column | Mapped internal name | Description |
|-----------------|-----------------------|-------------|
| `Family`        | `Family_ID`           | Selfed family identifier (character) |
| `ID`            | `Plant_ID`            | Individual plant identifier (character) |
| `PH_Nov`        | `Plant_Height_Nov`    | November plant height (numeric) |
| `IN_Nov`        | `Internode_Nov`       | November internode length (numeric) |
| `PH_Mar`        | `Plant_Height_Mar`    | March plant height (numeric) |
| `IN_Mar`        | `Internode_Mar`       | March internode length (numeric) |
| `BN_Mar`        | `Branch_Number_Mar`   | March branch count (numeric) |
| `MF_Score_Mar`  | `Multifoliate_Score_Mar` | March multifoliate score (numeric) |
| `PH_May`        | `Plant_Height_May`    | May plant height (numeric) – **used as key trait** |

All missing values are imputed by **within‑family median**; any remaining NAs (families with all missing) are filled by global median.

---

## 5. Key Parameters (Configurable)

Set these variables at the top of the script (lines ~50–60):

| Parameter | Value | Description |
|-----------|-------|-------------|
| `TARGET_TOTAL_N` | 200 | Hard target – final core set size (assertion enforced) |
| `MAX_PER_FAMILY` | 7 | Soft cap during Neyman allocation phase |
| `MIN_FAMILY_SIZE` | 3 | Families with fewer individuals are excluded from allocation |
| `MIN_ALLOC` | 1 | Minimum number of plants taken from any eligible family |
| `FILL_ALLOW_EXCESS` | 2 | Maximum extra plants per family during Phase 4 force‑fill (hard cap = `MAX_PER_FAMILY + FILL_ALLOW_EXCESS` = 9) |
| `all_traits` | (vector of 7 trait names) | Traits used in PCA and clustering – **do not modify unless column mapping is updated** |

---

## 6. Execution Workflow

1. **Package installation & font initialisation**
2. **Stage 1 – Data imputation & clustering**
   - Read CSV, map columns, impute missing values
   - Exclude families with < `MIN_FAMILY_SIZE` individuals
   - Perform K‑means clustering (optimal k by silhouette) and PCA visualisation (Fig S1)
3. **Stage 2 – Sampling pipeline**
   - **2.1** Compute variance components and repeatability (Fig 1)
   - **2.2** Run PCA on all traits; extract PC1/PC2 scores (Fig 2)
   - **2.3** Weighted Neyman allocation (Phases 1–3) with floor, soft cap, and integer adjustment
   - **2.3b** Within‑family stratified sampling by PC1 quantiles
   - **★ 2.3c Phase 4 – Force‑fill**  
     - Detect gap: `gap = TARGET_TOTAL_N - nrow(selected_plants)`  
     - If gap > 0, iteratively select additional plants from **remaining unselected** individuals in large families, prioritised by Neyman weight  
     - Allow per‑family count up to `MAX_PER_FAMILY + FILL_ALLOW_EXCESS`  
     - Continue until exact target is reached; assertion fails if impossible
   - **2.4** Sampling validation: K‑S test, density comparison, family‑wise strip plot (Fig 3)
   - **2.5** Cluster coverage cross‑tabulation
   - **2.6** Generate diagnostic report (GWASI)
4. **Output** – all tables, figures, and logs saved under a timestamped root directory.

---

## 7. Output Structure

All outputs are written to:  
`GWAS_Pipeline_Result_YYYYMMDD_HHMM/`

Subdirectories:

| Directory | Contents |
|-----------|----------|
| `01_Scripts/` | Copy of the executed R script (`GWAS_Pipeline_v7.1_Final.R`) |
| `02_Data_Tables/` | All core data files (imputed dataset, variance components, PCA eigenvalues/loadings/scores, selected core lists with/without traits, cluster coverage, diagnostic report) |
| `03_Figures/` | Publication‑ready PDF figures: Fig S1 (clustering), Fig 1 (segregation), Fig 2 (PCA biplot), Fig 3 (sampling justification) |
| `04_Plot_Data/` | CSV files underlying each figure (silhouette data, PCA cluster data, density data, allocation table, force‑fill log, strip data) |

**Key output files (under `02_Data_Tables/`):**

- `03_Selected_200_GWAS.csv` – core set (Family_ID, Plant_ID only)  
- `03_Selected_200_WithTraits.csv` – core set with cluster, PC scores, and all traits  
- `05_GWASI_Diagnostic_Report.csv` – summary metrics (Gini, K‑S, min/max allocation, Phase 4 details, etc.)  
- `Phase4_ForceFill_Log.csv` – step‑by‑step record of added plants during force‑fill (if triggered)

---

## 8. Diagnostic Metrics & Interpretation

- **Gini coefficient** – measures equality of allocation across families. Values < 0.35 indicate balanced representation; v7.1 typically yields ≤ 0.10.
- **Kolmogorov‑Smirnov (K‑S) test** – compares PC1 distribution of the core set versus the total eligible population. A non‑significant p‑value (p > 0.05) indicates that the core set preserves the overall phenotypic gradient, reducing the risk of ascertainment bias.
- **Repeatability (R)** – intra‑class correlation coefficient for each trait; higher values suggest stronger genetic control and better trait suitability for GWAS.
- **Phase 4 fill** – documented to allow transparent review of any post‑allocation additions. The fill is minimal (typically 3–5 plants) and drawn from large families, so kinship structure remains essentially unchanged.

---

## 9. Important Notes

- The script **stops with an error** if the final selected count does not exactly equal `TARGET_TOTAL_N`. This is intentional to enforce the supervisor’s hard requirement.
- The `MAX_PER_FAMILY` soft cap is **not** a hard limit during Phase 4 – families may reach up to `MAX_PER_FAMILY + FILL_ALLOW_EXCESS` (9) to accommodate the fill. This trade‑off is considered acceptable because the fill volume is tiny and spread across multiple families.
- All random sampling uses a fixed seed (`set.seed(2026)`) to ensure reproducibility.
- The script assumes the input file path `01_Raw_Data/data_202605_gwas.csv` relative to the working directory. Adjust if necessary.

---

## 10. Execution

```bash
Rscript GWAS_Integrated_Pipeline_v7.1.R
```

The console will display progress messages and a final summary with key statistics. Total runtime is typically **2–4 minutes** for a dataset of ~800 individuals and 30 families.

---

## 11. Contact & Maintenance

This pipeline is maintained for internal use only. For questions or modifications, refer to the principal investigator. No external support is provided.
