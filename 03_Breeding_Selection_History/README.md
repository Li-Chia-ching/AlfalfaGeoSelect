# Alfalfa Selection V3.0 – Publication‑Ready Visualization

**Version:** 3.0  
**Date:** 2026-07-11  
**Status:** Stable – single‑purpose visualization script for alfalfa breeding trials

---

## 1. Overview

This R script generates a **publication‑ready contour plot** that visualises the combined selection index for alfalfa (*Medicago sativa*) plants evaluated for two traits: plant height and multifoliate score. It reads raw field data from two CSV files, computes a weighted geometric index, selects the top 50 individuals, and produces a high‑resolution scatter‑contour graphic with clear annotations, a selection threshold line, and a gradient arrow.

The script is designed for **field grid experiments** (rows L1–L50, columns A–U) and performs all data cleaning, index calculation, and visualisation in a fully automated pipeline. No interactive input is required.

---

## 2. Key Features

- **Flexible data loading** – prefers in‑memory data frames (`initial_flowering_2025`, `Alfalfa_Multi_202504`) if present; otherwise reads from CSV files.
- **Robust data cleaning** – filters grid rows (`^L\d+$`), converts all columns to numeric, and reshapes from wide to long format.
- **Selection index** – defined as `(0.5 × Height) × (0.5 × Multi_Score)`. The top 50 individuals by this index are selected.
- **Contour surface** – generated on a 100×100 interpolation grid over the observed height range and a fixed multifoliate range (1–5).
- **Clear visual cues** – includes:
  - Colour‑coded contour lines (viridis palette) for the index surface.
  - Red dashed line marking the selection cutoff value.
  - Jittered data points coloured by selection status (selected vs. not selected), with size variation for emphasis.
  - Annotated arrow indicating the selection gradient direction.
  - Cutoff value label placed at the upper‑right corner.
- **Vector‑based output** – saves the plot as a 300 DPI PNG file (suitable for journal submission) and exports the selection list as a CSV.

---

## 3. Input Data Specification

The script expects **two CSV files** (or equivalent in‑memory data frames) with the following structure:

### File 1: Plant height data (`Rawdata-PlantHeight_202504.csv`)
- **First column**: Row identifiers (e.g., `L1`, `L2`, … `L50`).
- **Subsequent columns**: Column identifiers (e.g., `A`, `B`, … `U`) containing numeric height values (cm).
- Missing values are allowed and will be coerced to `NA`.

### File 2: Multifoliate data (`Rawdata-Multifoliate_202504.csv`)
- Identical layout: first column rows, subsequent columns for columns A–U.
- Values should be numeric multifoliate scores (expected range 1–5, but script does not enforce bounds).

**Important:** The script will first check if data frames named `initial_flowering_2025` and `Alfalfa_Multi_202504` exist in the global environment. If found, these are used preferentially, **bypassing file reading**. This facilitates integration with larger workflows. If not found, it falls back to the CSV files in the working directory.

---

## 4. Output

All outputs are saved in a **newly created timestamped directory** named:

```
Alfalfa_GeoViz_Pro_YYYYMMDD_HHMMSS/
```

### Files:

| File | Description |
|------|-------------|
| `geometric_selection.png` | 10 × 7 inch PNG figure at 300 DPI, ready for publication. |
| `selection_list.csv` | CSV containing the 50 selected plants with columns: `Row_ID`, `Col_ID`, `Plant_ID`, `Height`, `Multi_Score`, `Index`. |

---

## 5. Execution

### Prerequisites
- **R version** ≥ 3.6.0
- **Required packages** (automatically installed if missing):
  - `dplyr`, `tidyr`, `ggplot2`, `stringr`, `readr`, `metR`

### Running the script
Place the script in your working directory alongside the two input CSV files (or ensure the data frames are loaded in memory). Then execute:

```bash
Rscript Alfalfa_Selection_V3.0.R
```

The script will:
1. Check and install any missing packages.
2. Load the data (preferring in‑memory objects).
3. Clean and merge the data.
4. Compute indices and select the top 50.
5. Generate the contour plot.
6. Create the output directory and save the figure and selection list.

Console output will display a summary and the name of the output directory.

---

## 6. Configuration Parameters

The script contains a few hard‑coded values that can be adjusted by editing the source:

| Parameter | Value | Description |
|-----------|-------|-------------|
| `all_rows` | `paste0("L", 1:50)` | Row identifiers in the field grid. |
| `all_cols` | `LETTERS[1:21]` | Column identifiers (A–U). |
| `top_n` | `50` (embedded in `slice_head(n = 50)`) | Number of top‑performing individuals to select. |
| `grid_x` length | `100` | Resolution of the contour surface along the height axis. |
| `grid_y` range | `1` to `5` | Fixed multifoliate score range for the contour grid. |
| `contour bins` | `12` | Number of contour intervals. |

If you need to change these, edit the corresponding lines in the script.

---

## 7. Interpretation of the Figure

- **Contour lines** represent the selection index value across the trait space.
- **Red dashed line** indicates the minimum index required to be among the top 50.
- **Points** are individual plants; orange points are selected, grey points are not. Larger point size highlights selected individuals.
- **Arrow and label** indicate the direction of increasing selection index (towards higher height and higher multifoliate score).

The figure is self‑contained with a title, subtitle, axis labels, and a caption explaining the index formula.

---

## 8. Notes and Limitations

- The script **assumes** the row identifiers start with `L` followed by digits. Any rows not matching this pattern are discarded.
- Multifoliate scores are expected to be in the range 1–5; the contour grid uses this fixed range for consistency. If actual data exceed 5, the contour surface will still be drawn up to 5, which may clip some points. Adjust `grid_y` accordingly if needed.
- The plot uses `theme_minimal()` and is designed for clarity rather than decorative styling. The colour palette (viridis) is colour‑blind friendly.
- The script does **not** perform any spatial or statistical correction; it is purely a visualisation tool for selection outcomes.
- No log files are written; only the final output directory is created.

---

## 9. Dependencies & Versioning

This script is self‑contained and does not rely on external configuration files. It is intended for use with the **2025/2026 alfalfa field trial data** and has been tested with R 4.2+ on Windows and Linux.

For questions or modifications, refer to the original author. This script is not part of the integrated GWAS pipeline and serves a standalone reporting function.
