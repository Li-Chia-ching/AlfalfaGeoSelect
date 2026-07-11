# Automated Crossing Pair Recommendation Script (Post‑GWAS Validation)

**Version:** 1.0  
**Date:** 2026-07-11  
**Status:** Production – generates validated crossing plans from core GWAS collection

---

## 1. Overview

This R script reads a **core set of 200 genotyped and phenotyped individuals** (the output of a preceding GWAS sampling pipeline) and generates a **ranked list of recommended crossing pairs** for two target traits: plant height and multifoliate leaf score (MF_Total). The recommendations are based on two complementary genetic strategies:

- **Strong × Weak** – pairs extreme phenotypic performers to create segregating populations for QTL mapping and transgressive segregation studies.
- **Strong × Strong** – pairs high‑performing individuals from different families to pyramid favourable alleles and test allelism.

The script produces a single CSV file containing the top 5 pairs per strategy per trait, with full parental identifiers, family origins, trait values, and phenotypic differences. This plan is intended for **post‑GWAS validation crosses** and is designed for immediate use in field breeding programs.

---

## 2. Input Requirements

### Input file (hard‑coded path – **please adjust**)

The script expects a CSV file at:

```
GWAS_Pipeline_Robust200_20260512/02_Selected_200_GWAS_FullTraits.csv
```

**Columns required in this file:**

| Column name | Type | Description |
|-------------|------|-------------|
| `ID` | character | Unique plant identifier (must match field labels) |
| `Family` | character | Family (line) designation |
| `Plant_Height` | numeric | Phenotypic value for plant height (cm) |
| `MF_Total` | numeric | Phenotypic value for multifoliate leaf score (total) |

*Additional columns may be present but are ignored.*

If the file is not found at this location, the script stops with an error message. To adapt to your directory structure, **edit the `input_file` variable** at line 13 of the script.

---

## 3. Output

The script writes a single CSV file to the working directory:

```
03_RealData_Recommended_Crossing_Plan.csv
```

### Output columns

| Column | Description |
|--------|-------------|
| `Target_Trait` | Trait used for selection (`Plant_Height` or `MF_Total`) |
| `Cross_Type` | Strategy: `Strong x Weak` or `Strong x Strong` |
| `P1_ID` | Identifier of the first parent |
| `P1_Family` | Family of P1 |
| `P1_Value` | Phenotypic value of P1 for the target trait |
| `P2_ID` | Identifier of the second parent |
| `P2_Family` | Family of P2 |
| `P2_Value` | Phenotypic value of P2 for the target trait |
| `Pheno_Diff` | Absolute difference between P1_Value and P2_Value |

Each row represents one recommended crossing combination. The top 5 pairs are given for each `(Trait, Cross_Type)` combination, for a total of 20 recommendations (2 traits × 2 strategies × 5 pairs).

---

## 4. Algorithm Description

### 4.1. Pool construction
For each target trait, all 200 individuals are sorted by trait value in descending order. The **top 15** form the **strong pool**; the **bottom 15** form the **weak pool**.

### 4.2. Strategy A – Strong × Weak
- Generates all possible combinations between the strong and weak pools (15×15 = 225 combinations).
- Computes the absolute phenotypic difference for each pair.
- Sorts by difference descending.
- Deduplicates to ensure each **strong parent** appears only once and each **weak parent** appears only once in the final list (a critical step to avoid overusing a single parent).
- Keeps the top 5 pairs after deduplication.

### 4.3. Strategy B – Strong × Strong
- Generates all **unordered** combinations between the strong pool members (excluding self‑crosses, i.e., P1_ID < P2_ID).
- Filters out pairs coming from the **same family** to maintain genetic diversity.
- Computes the sum of the two parental values (`Combined_Strength`) and sorts descending.
- Keeps the top 5 pairs (no further deduplication, as each strong parent can appear in multiple strong×strong crosses).

### 4.4. Aggregation
Results from both strategies and both traits are merged into a single table. Column order is standardised for readability.

---

## 5. Execution

### Prerequisites
- **R version** ≥ 3.5.0
- **Packages:** `dplyr`, `tidyr`, `readr` (automatically loaded; will fail if not installed)

### Running
Place the script in your working directory with the input CSV file (or adjust the path). Then execute:

```bash
Rscript Crossing_Recommendation_Script.R
```

The script will:
1. Validate existence of the input file.
2. Read the data.
3. Generate recommendations for both traits.
4. Write the output CSV.
5. Print a preview to the console.

---

## 6. Customisation

You can adjust three parameters within the `generate_crosses()` function call (lines 59 and 62):

| Parameter | Default | Description |
|-----------|---------|-------------|
| `pool_size` | 15 | Number of top (and bottom) individuals to consider for each trait |
| `top_n_pairs` | 5 | Number of crossing pairs to output per strategy per trait |

For example, to increase the pool to 20 and output 8 pairs per strategy, change:

```r
crosses_height <- generate_crosses(df_200, trait_col = "Plant_Height", pool_size = 20, top_n_pairs = 8)
```

---

## 7. Important Notes

- The script **does not** check for missing trait values. Ensure your input data are complete for the two target traits.
- The strong×weak deduplication is **bidirectional** – it first retains distinct `P1_ID`s, then distinct `P2_ID`s. This may result in fewer than `top_n_pairs` if the pool is small; in practice, with 15 individuals each, 5 pairs are always attainable.
- For strong×strong, no deduplication is applied, so a single parent may appear in multiple recommended crosses – this is intentional to allow multiple uses of elite parents.
- The script assumes family and individual IDs are character strings with no leading/trailing spaces.
- All random elements are absent; the selection is deterministic based on trait ranking.

---

## 8. Output File Name

The output is saved as `03_RealData_Recommended_Crossing_Plan.csv` to distinguish it from any simulated or preliminary crossing lists. If this file already exists, it will be overwritten without warning.

---

## 9. Limitations

- Only two traits are currently hard‑coded (`Plant_Height`, `MF_Total`). To add more traits, extend the script by calling `generate_crosses()` with additional `trait_col` values.
- No consideration is given to genotype or marker data – the recommendations are purely phenotype‑based.
- The fixed input path may cause errors if the directory structure changes; manual editing of line 13 is required.

---

## 10. Contact

This script is an internal tool for the breeding programme. For modifications or questions, consult the project lead. No external support is provided.
