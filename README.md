# AlfalfaGeoSelect

**AlfalfaGeoSelect** is an R-based analytical and visualization pipeline for selecting superior alfalfa ( *Medicago sativa* ) hybrid candidates using a weighted geometric selection index derived from plant height and multifoliate score.

The workflow is designed for **breeding trials**, **phenotypic screening**, and **publication-ready visualization**.

---

## Features

- Automated phenotypic data cleaning and reshaping  
- Intelligent data loading (workspace object or CSV fallback)  
- Weighted geometric selection index construction  
- Top-N candidate identification with explicit cutoff  
- High-quality contour-based visualization suitable for journals  
- Fully reproducible and script-driven

---

## Selection Index Definition

The selection index is defined as:

\[
\text{Index} = (0.5 \times \text{Height}) \times (0.5 \times \text{Multifoliate Score})
\]

Plants ranked above the cutoff value (Top 50 by default) are selected as elite hybrid candidates.

---

## Workflow Overview

1. **Environment setup**  
   Required R packages are checked and installed automatically.

2. **Data loading & cleaning**  
   - Accepts either in-memory data frames or CSV files  
   - Converts wide field layouts into long-format plant-level data  

3. **Index calculation & selection**  
   - Computes weighted geometric index  
   - Classifies plants as *Selected*, *Not Selected*, or *Dead/Missing*

4. **Visualization**  
   - Contour map of the selection surface  
   - Explicit cutoff line  
   - Highlighted elite individuals  
   - Vector-based directional annotation

5. **Output**  
   - High-resolution PNG figure  
   - CSV list of selected plants  
   - Timestamped output directory

---

## Input Data Format

- **Rows**: L1–L50 (field rows)
- **Columns**: A–U (field columns)
- First column must contain row identifiers (e.g., `L1`, `L2`, …)

Supported inputs:
- `Rawdata-PlantHeight_YYYYMM.csv`
- `Rawdata-Multifoliate_YYYYMM.csv`

---

## Output Structure

```text
Alfalfa_GeoViz_Pro_YYYYMMDD_HHMMSS/
├── Plot_Geometric_Selection_Pro.png
└── 01_Selection_List.csv
```

------------

## Visualization Example

- Grey points: not selected
- Orange points: top-ranked candidates
- Red dashed contour: selection cutoff
- Background contours: index gradient

------------

## Dependencies

- dplyr
- tidyr
- ggplot2 (≥ 3.4.0 recommended)
- stringr
- readr

## Intended Use

- Alfalfa hybrid breeding programs
- Multi-trait phenotypic screening
- Trait integration method demonstration
- Reproducible research and publication figures

---

## License

MIT License
