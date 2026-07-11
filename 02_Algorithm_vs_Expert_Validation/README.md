# Expert‑Adjusted GWAS Core Panel Diagnostics

**Version:** 1.0  
**Date:** 2026-07-11  
**Status:** Validation script for post‑adjustment core collection assessment

---

## 1. Overview

This R script performs a **diagnostic evaluation** of a manually refined GWAS core panel. Starting from an automatically selected 200‑individual set, an expert (breeder) has made editorial changes—adding and removing specific plants—resulting in a final cohort of **212 individuals**. The script compares the final expert‑adjusted cohort against the original candidate pool to assess:

- **Compositional balance** – family‑wise distribution of the final set.
- **Representativeness** – whether the phenotypic distribution (plant height) has been significantly altered by the expert adjustments.
- **Changes in membership** – which plants were added and which were removed.

The output includes a publication‑ready two‑panel figure (family distribution bar chart + density comparison with Kolmogorov‑Smirnov test), summary statistics (Gini coefficient), and detailed data tables.

---

## 2. Input Data Requirements

The script expects three input files placed in the `./01_Raw_Data/` directory (relative to the working directory). All paths are **relative**—no absolute paths are used.

### File 1: `Original_Selected_200_WithTraits_Labeled.xlsx`
- **Sheet:** `03_Selected_200_WithTraits_Labe`
- **Required columns:** `Family_ID`, `Plant_ID`, `PC1`, `PC2`, `Plant_Height_May`
- **Description:** The original automatically selected 200 individuals with their family assignments, PCA coordinates, and May plant height.

### File 2: `Latest_Selected_For_Height+Color.xlsx`
- **Sheet:** `Latest_Selected`
- **Columns:** The script reads **columns 1 and 4**:
  - Column 1: `Plant_ID` (plant identifier)
  - Column 4: `Height_May` (plant height in May, numeric)
- **Description:** The expert‑adjusted final list of plants (expected 212 rows). Additional columns are ignored.

### File 3: `data_202605_gwas.csv`
- **Required columns:** First column = `Family_ID`, second column = `Plant_ID`
- **Description:** The master GWAS database used to look up family assignments for **newly added** plants that were not present in the original 200 set. The script assumes the first two columns are `Family_ID` and `Plant_ID`.

If the input files are not found or if the final list does not contain exactly 212 rows, a warning is issued (but the script continues).

---

## 3. Output Directory and Files

All outputs are saved to a timestamped subdirectory under `./05_Analysis_Outputs/`:

```
05_Analysis_Outputs/Panel_Diagnostics_212_YYYYMMDD_HHMM/
```

### Output files

| File | Content |
|------|---------|
| `Figure_Panel_Diagnostics_212.pdf` | Two‑panel figure (A: family distribution bar chart; B: density comparison of plant height with K‑S test annotation) |
| `Family_Distribution_212.csv` | Table of family counts in the final cohort, sorted descending |
| `01_Added_Plant_IDs.txt` | List of plant IDs added by the expert (present in final but not in original) |
| `01_Removed_Plant_IDs.txt` | List of plant IDs removed by the expert (present in original but not in final) |
| `KS_Test_Result.txt` | Full output of the Kolmogorov‑Smirnov test comparing the height distributions of the original and final panels |

---

## 4. Methodology

### 4.1. Data cleaning and family mapping
- Both plant ID lists are trimmed of whitespace.
- The final list is taken as‑is (no truncation) – expected to be 212 plants.
- For each plant in the final list, family assignment is resolved using a **priority cascade**:
  1. Manual dictionary (`map_manual`) – covers known special ID patterns (e.g., AL06, HY20).
  2. GWAS database look‑up (`map_gwas`) – for plants not in the original 200.
  3. Original mapping (`map_original`) – for plants that were already in the 200 set.
  4. Fallback: `"Unknown"` if none of the above apply.

### 4.2. Difference analysis
- **Added:** Plants in final but not in original.
- **Removed:** Plants in original but not in final.
- These lists are saved as plain text files for record‑keeping.

### 4.3. Family distribution and Gini coefficient
- The script tabulates the number of selected plants per family in the final cohort.
- The **Gini coefficient** is calculated as a measure of inequality in family representation. Values close to 0 indicate equal representation; values >0.35 suggest imbalance. This is computed but not plotted—use for internal diagnostics.

### 4.4. Kolmogorov‑Smirnov test
- The script compares the distribution of plant height (`Plant_Height_May` for original, `Height_May` for final) between the two cohorts using a two‑sample K‑S test.
- The test statistic (D) and p‑value are displayed directly on the density plot and saved to a text file.

### 4.5. Visualisation
- **Panel A:** Bar chart of family counts, with a dashed horizontal line at the mean count. The family with the highest count is annotated with its `n`.
- **Panel B:** Overlaid kernel density plots for the original (grey‑blue) and final (sand‑orange) cohorts, using colour‑blind‑friendly palettes. The K‑S test result is embedded in the plot.
- Both panels are combined using `patchwork` with panel labels (A, B).

---

## 5. Execution

### Prerequisites
- **R version** ≥ 4.0.0
- **Packages:** `tidyverse`, `readxl`, `patchwork`, `extrafont`
  - `extrafont` is used to ensure Arial font availability; if not installed, the script will still run but fall back to default fonts. On Windows, you may need to run `loadfonts(device = "win")` once.

### Running the script
Place the script in your project root (where `./01_Raw_Data/` and `./05_Analysis_Outputs/` are subdirectories). Then execute:

```bash
Rscript Expert_Adjusted_Panel_Diagnostics.R
```

The script will read all input files, perform the analyses, and create the output directory with all results. Console messages indicate progress and the final output location.

---

## 6. Customisation

### Adjusting input/output paths
The script uses relative paths defined at lines 45–47. To change the data or output root, modify:

```r
data_dir <- "./01_Raw_Data"
output_root <- "./05_Analysis_Outputs"
```

### Changing the expected final sample size
If your expert‑adjusted list contains a different number of plants (e.g., 215), update the warning threshold at line 75:

```r
if (nrow(df_final) != 212) {  # change 212 to your expected N
```

### Adding manual family mappings
To extend the manual dictionary, add entries to `map_manual` (lines 108–113) following the same `ID = "Family"` pattern.

### Modifying colours or theme
The colour palette `cb_palette` and theme `theme_gwas` are defined near the top; edit as desired.

---

## 7. Important Notes

- The script **does not** filter or truncate the final list; it uses all rows from the expert‑provided Excel sheet. Ensure that the sheet contains exactly the plants you wish to evaluate.
- The K‑S test compares **plant height distributions**, not PC1 or other traits. If you need to test other traits, modify the `height_orig` and `height_fin` extraction accordingly.
- No significance threshold is automatically applied; the p‑value is reported for the user to interpret.
- The script assumes that `Plant_Height_May` in the original file and `Height_May` in the final file are **comparable** (same measurement unit and protocol).

---

## 8. Interpretation of Results

- **Family distribution plot:** A wide spread indicates some families are over‑represented. The mean line gives a baseline for equal allocation. Large deviations (especially n > mean + 2×SD) may indicate bias introduced by expert selection.
- **Density comparison:** Overlap between the two density curves suggests the expert adjustments did not drastically shift the phenotypic distribution. A **non‑significant K‑S test** (p > 0.05) supports that the final cohort remains representative of the original pool. A significant result (p < 0.05) indicates that the expert intentionally shifted the distribution toward higher (or lower) heights.
- **Added/Removed lists:** Review these to confirm that the expert changes align with breeding objectives.

---

## 9. Limitations

- The script does not perform any multiple‑testing correction; the K‑S p‑value is reported as‑is.
- Family mapping relies on the manual dictionary; incomplete dictionaries may result in `"Unknown"` families.
- The Gini coefficient is calculated but not saved in the output table (only printed if you modify the script to show it). For record‑keeping, you may wish to add an explicit output for this value.
- The script uses `cairo_pdf` as the graphics device; on some Linux systems without Cairo support, you may need to change the `device` argument to `pdf()`.

---

## 10. Contact

This script is an internal diagnostic tool for the breeding programme. For modifications or questions, consult the project lead. No external support is provided.
