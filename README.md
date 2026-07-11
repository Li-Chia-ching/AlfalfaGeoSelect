# GWAS‑Driven Alfalfa Breeding Pipeline — R Script Suite

**Repository:** Alfalfa-GWAS-Tools  
**Last Updated:** 2026-07-11  
**Status:** Production – validated for internal use

---

## 1. Repository Overview

This repository contains a **complete suite of R scripts** for phenotype‑driven core collection sampling, expert panel validation, crossing recommendation, and publication‑ready visualisation in an alfalfa (*Medicago sativa*) breeding programme. The pipeline is designed to be run sequentially, starting from raw field data and ending with recommended crossing plans and diagnostic figures.

The scripts are modular, self‑contained, and rely only on standard CRAN packages. They are intended for users with moderate R proficiency and are fully documented with inline comments.

---

## 2. Script Inventory

| Script | Purpose | Input | Output |
|--------|---------|-------|--------|
| `Alfalfa_Selection_V3.0.R` | Visualises selection index (height × multifoliate) via contour plot, selects top 50 individuals | CSV files or in‑memory data frames with row/column grid data | PNG figure + CSV of selected 50 |
| `GWAS_Integrated_Pipeline_v7.1.R` | Performs imputation, clustering, Neyman allocation, stratified sampling, and **force‑fill** to exact 200 individuals for GWAS | `data_202605_gwas.csv` (raw field data) | Hierarchical output directory with imputed data, PCA, variance components, selected core sets, figures, and diagnostic report |
| `Crossing_Recommendation.R` | Generates strong×weak and strong×strong crossing pairs for two traits from a core set of 200 | `GWAS_Pipeline_Robust200_*/02_Selected_200_GWAS_FullTraits.csv` | `03_RealData_Recommended_Crossing_Plan.csv` (top pairs) |
| `Expert_Panel_Diagnostics.R` | Compares original automatically selected 200 vs expert‑adjusted 212‑plant cohort; assesses balance and representativeness | Two Excel files + GWAS database CSV | PDF figure, family distribution table, added/removed ID lists, K‑S test result |

*Note: Script names in this table are placeholders – actual file names may vary. See each script’s header for exact naming.*

---

## 3. System Requirements

- **R version:** ≥ 4.0.0 (3.6+ may work but not tested)
- **Operating system:** Windows / Linux / macOS (path separators and font handling are cross‑platform aware)
- **Required R packages** (automatically installed if missing in most scripts):
  - Core: `dplyr`, `tidyr`, `ggplot2`, `readr`, `stringr`
  - Extended: `lme4`, `lmerTest`, `factoextra`, `patchwork`, `showtext`, `cluster`, `viridis`, `metR`, `readxl`, `extrafont`
- **Fonts:** The GWAS pipeline and diagnostics scripts attempt to use Arial (or SimSun for Chinese labels). If fonts are not available, the scripts fall back to system defaults without error.

---

## 4. Recommended Workflow

1. **Collect raw field data** in the required CSV format (see individual script documentation for exact column names).
2. **Run `GWAS_Integrated_Pipeline_v7.1.R`** to obtain a balanced core set of 200 individuals, along with all supporting data tables and figures.
3. **Perform expert review** of the 200‑set, adding/removing plants as needed, and save the final list as an Excel file.
4. **Run `Expert_Panel_Diagnostics.R`** to validate that the expert adjustments have not compromised representativeness.
5. **Run `Crossing_Recommendation.R`** on the final core set (or the original 200) to generate crossing plans for post‑GWAS validation.
6. **Optionally, use `Alfalfa_Selection_V3.0.R`** for a quick visual summary of the selection index if you are working with the original grid data.

All scripts are independent and can be executed separately; however, the intended order follows the steps above.

---

## 5. Detailed Script Descriptions

### 5.1. `Alfalfa_Selection_V3.0.R`
**Purpose:** Single‑file visualisation of a selection index derived from plant height and multifoliate score.  
**Key features:**
- Flexible data loading (prefers in‑memory objects over CSV)
- Automatic grid construction (L1–L50 × A–U)
- Index = (0.5 × Height) × (0.5 × Multi_Score)
- Top‑50 selection with contour plot, cutoff line, and gradient annotation
- Output: 300 DPI PNG + selection list CSV

### 5.2. `GWAS_Integrated_Pipeline_v7.1.R`
**Purpose:** Core collection sampling with exact‑N enforcement.  
**Key features:**
- Within‑family median imputation
- Exclusion of families with <3 individuals
- Phenotypic clustering (K‑means + silhouette) and PCA
- Variance component estimation (repeatability) for all traits
- Neyman allocation with per‑family soft cap (7) and floor (1)
- Stratified sampling by PC1 quantiles within each family
- **Phase 4 Force‑Fill:** automatically fills any shortfall to reach exactly 200 plants by drawing from large families (up to cap 9)
- Comprehensive output: imputed data, PCA tables, variance components, selected core sets (with/without traits), figures (S1, 1, 2, 3), and a detailed diagnostic report (GWASI)
- Hierarchical output structure: `01_Scripts/`, `02_Data_Tables/`, `03_Figures/`, `04_Plot_Data/`

### 5.3. `Crossing_Recommendation.R`
**Purpose:** Generate phenotype‑based crossing recommendations from a core set.  
**Key features:**
- Two strategies: Strong × Weak (extreme segregation) and Strong × Strong (allele pyramiding)
- For each strategy, selects top 5 pairs per trait (plant height and MF_Total)
- Strong×Weak: deduplicates parents to avoid overuse
- Strong×Strong: requires parents from different families
- Output: single CSV with all recommended pairs, including parental IDs, families, trait values, and phenotypic differences
- Hard‑coded input path – **must be edited** to match your output from the GWAS pipeline

### 5.4. `Expert_Panel_Diagnostics.R`
**Purpose:** Validate expert‑adjusted final cohort against the original automatically selected set.  
**Key features:**
- Reads original 200‑set and expert final list (expected 212 plants)
- Identifies added/removed plants
- Maps family IDs via manual dictionary, GWAS database, or original mapping
- Computes family distribution and Gini coefficient
- Performs Kolmogorov‑Smirnov test on plant height distributions
- Generates two‑panel figure: family bar chart + density comparison with K‑S annotation
- Outputs: PDF figure, family distribution table, added/removed lists, K‑S test result text file
- **Uses relative paths** – adapt to your folder structure

---

## 6. Data Formatting Guidelines

Each script expects specific column names and file formats. Please refer to the individual README files (or the script headers) for precise specifications. In summary:

| Script | Input file(s) | Required columns |
|--------|---------------|------------------|
| Alfalfa Selection | `Rawdata-PlantHeight_202504.csv`, `Rawdata-Multifoliate_202504.csv` | First column: row IDs (L1…L50); subsequent: column letters (A…U) |
| GWAS Pipeline | `data_202605_gwas.csv` | `Family`, `ID`, `PH_Nov`, `IN_Nov`, `PH_Mar`, `IN_Mar`, `BN_Mar`, `MF_Score_Mar`, `PH_May` |
| Crossing Recommendation | `GWAS_Pipeline_Robust200_*/02_Selected_200_GWAS_FullTraits.csv` | `ID`, `Family`, `Plant_Height`, `MF_Total` |
| Expert Diagnostics | `Original_Selected_200_WithTraits_Labeled.xlsx` (sheet `03_Selected_200_WithTraits_Labe`), `Latest_Selected_For_Height+Color.xlsx` (sheet `Latest_Selected`), `data_202605_gwas.csv` | See script details |

All scripts use `read_csv`/`read_excel` and assume files are in the working directory or specified relative paths. Edit the path variables if your structure differs.

---

## 7. Execution Instructions

1. **Clone the repository** to your local machine.
2. **Place raw data files** in the appropriate subdirectories (e.g., `./01_Raw_Data/` for the pipeline, `./data/` for the visualisation script).
3. **Open R or RStudio** and set the working directory to the repository root.
4. **Run each script** individually using `source("script_name.R")` or via command line `Rscript script_name.R`.
5. **Check console output** for progress messages and any warnings/errors.

All scripts create their own output directories with timestamps to avoid overwriting previous runs.

---

## 8. Customisation and Parameter Tuning

Key parameters are clearly marked at the top of each script. You can adjust:

- **GWAS pipeline:** `TARGET_TOTAL_N`, `MAX_PER_FAMILY`, `MIN_FAMILY_SIZE`, `MIN_ALLOC`, `FILL_ALLOW_EXCESS`
- **Crossing recommendation:** `pool_size`, `top_n_pairs` (inside function calls)
- **Alfalfa selection:** grid dimensions, contour resolution, selection threshold (top N)
- **Diagnostics:** output directory paths, colour palettes, expected final sample size (212)

Always review the parameters before running to match your experimental design.

---

## 9. Output Interpretation

- **Core selection sets** are saved as CSV files with full trait and PCA coordinates – use these for downstream GWAS.
- **Diagnostic reports** (GWASI) provide Gini coefficient, K‑S test results, allocation min/max, and Phase 4 fill details – use these to justify the sampling strategy in publications.
- **Crossing plans** give parent IDs and trait differences – use these to design field crosses for validation.
- **Figures** are publication‑ready (PDF/PNG at 300+ DPI) – include them directly in manuscripts.

---

## 10. Troubleshooting and Common Issues

| Issue | Suggested resolution |
|-------|----------------------|
| `Error: Input file not found` | Verify file paths; scripts use relative paths – adjust variables to match your directory layout. |
| `Package ‘xxx’ not available` | Run `install.packages("xxx")` manually; ensure you have a working internet connection. |
| Font warnings (SimSun/Arial) | Ignore – figures will use default system fonts; no functional impact. |
| `nrow(df_final) != 212` warning | The expert diagnostics script expects 212 plants; if yours differ, update the warning threshold or ignore if intentional. |
| `lmer` convergence failures | The pipeline handles these gracefully and records `NA` for repeatability; still proceeds. |
| Force‑fill unable to reach target | The script will stop with an error; this indicates insufficient remaining individuals – check family sizes or relax `MAX_PER_FAMILY`. |

---

## 11. Dependencies and Version Control

- All scripts are written in base R with tidyverse packages; they have been tested with R 4.2.3 on Windows 10 and Ubuntu 20.04.
- No external software or proprietary libraries are required.
- The repository is version‑controlled using Git. Tagged releases correspond to major pipeline versions (v6.0, v7.1, etc.).

---

## 12. License and Usage

This code is proprietary to the breeding programme and is provided for internal use only. Redistribution or commercial use is prohibited without explicit permission from the principal investigator.

---

## 13. Contact and Support

For questions, bug reports, or feature requests, contact the repository maintainer or the project lead. No external support is offered.
