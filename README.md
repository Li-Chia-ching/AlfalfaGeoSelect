-----

# 🧬 Genetic Variation Analysis & GWAS Core Subset Sampling Pipeline (v5.3)

An end-to-end R pipeline designed for multidimensional phenotypic gradient analysis and stratified sampling across two-year datasets. Specifically tailored for **selfed family populations**, this script utilizes Linear Mixed Models (LMM), Principal Component Analysis (PCA), and variance-weighted Neyman optimal allocation to extract a highly representative core subset for Genome-Wide Association Studies (GWAS).

## 🚀 What's New in v5.3

  * **📊 Quantitative Representativeness Validation (K-S Test):** Introduced a **Two-Sample Kolmogorov-Smirnov (K-S) test** to the final sampling validation step (Fig 3B). By calculating the *D*-statistic (maximum vertical distance) and *p*-value, the script now objectively quantifies the distributional overlap of the core phenotype (PC1) between the "Total Population" and the "Sampled Core Subset."
  * **⚙️ Dynamic Global Parameter Control:** Hardcoded limits have been replaced with a flexible configuration zone at the top of the script. You can now easily set `TARGET_TOTAL_N` (e.g., 140) and `MAX_PER_FAMILY` (e.g., 20). The script uses a D'Hondt-like smoothing algorithm to redistribute excess weights, ensuring the final sample size matches your target **exactly 100% of the time**.

-----

## 🛠️ Key Features

  - [x] **Multi-Year Data Fusion:** Automatically merges 2025 and 2026 phenotype data, calculates inter-annual growth dynamics, and filters out dead/missing records.
  - [x] **Phenotypic Segregation & Repeatability:** Leverages LMM (`lme4`) to extract variance components and calculate trait repeatability, generating publication-ready kernel density and boxplots (Fig 1).
  - [x] **Stepwise PCA Dimensionality Reduction:** Performs separate PCAs for cross-year shared traits and single-year specific traits (Fig 2 biplots). Extracts PC1 from a comprehensive 6-trait PCA as the quantitative baseline for overall growth performance.
  - [x] **Capped Neyman Optimal Allocation:** Distributes sample slots based on family population size and multivariate variance (PC1 & PC2), while strictly enforcing the user-defined per-family cap.
  - [x] **Stratified Quantile Sampling:** Divides individuals within each family into high/medium/low PC1 gradients, ensuring the deterministic capture of both extreme (tail variations) and median phenotypes.

-----

## 📦 Dependencies

The script automatically checks for and installs missing packages on its first run. Ensure you have an active internet connection.

  * **Data Wrangling:** `dplyr`, `tidyr`, `readr`
  * **Visualization:** `ggplot2`, `patchwork`, `showtext`, `factoextra`
  * **Modeling:** `lme4`, `lmerTest`

> **⚠️ Font Rendering Note:** The script uses `showtext` to render academic-compliant charts with the `SimSun` font. If you are on Linux/macOS and encounter font errors, ensure you have the appropriate TTF/TTC font files (like `STSong.ttf`) installed on your system.

-----

## 🚦 Usage

### 1\. Prepare Your Data

Place the following two CSV files in the same working directory as the R script:

  * `data_2025.csv` (Required columns: Parent, SampleID, PlantHeight, NodeCount, etc.)
  * `data_2026.csv` (Required columns: Group, Matrix, Plant\_Height, Internode, Branch\_Number, Multifoliate\_Score, Death\_Code, etc.)

### 2\. Configure Your Parameters

Open the R script and locate the `# 0. User Parameters` section at the very top. Adjust your sampling goals:

```r
TARGET_TOTAL_N <- 140  # Your target sequencing sample size for GWAS
MAX_PER_FAMILY <- 20   # Hard cap to prevent a single dominant family from monopolizing the subset
```

### 3\. Run the Pipeline

Simply source the script in your R console or RStudio:

```r
source("GWAS_Sampling_Pipeline_v5.3.R")
```

-----

## 📂 Output Structure

Upon completion, the script generates a timestamped directory (e.g., `GWAS_Sampling_Two-Year_YYYYMMDD`) containing three main categories of deliverables:

### 📈 Publication-Ready Plots (PDF, 600 DPI)

  * **`Fig1_Segregation_Analysis.pdf`**: Evidence of phenotypic segregation (density distributions and boxplots).
  * **`Fig2_PCA_Stepwise_Biplots.pdf`**: PCA biplots with vector arrows showing trait co-evolution.
  * **`Fig3_Sampling_Justification.pdf`**: **The definitive proof of sampling quality.**
      * A: Neyman allocation bar charts.
      * **B: K-S Test Probability Density Plot. Quantifies how well the 140-core subset mirrors the total population's distribution.**
      * C: Within-family strip plots showing the coverage across high/medium/low gradients.

### 📊 Result Matrices (CSV)

  * `01_Merged_Wide_Data_All.csv`: Cleaned and merged dataset.
  * `02_Variance_Components_Repeatability.csv`: LMM-derived variance and repeatability metrics.
  * `03_Comprehensive_PCA_Scores.csv`: PC scores for all living individuals.
  * **`04_Selected_140_GWAS.csv`**: The final roster of the 140 selected core individuals (ready to be sent to your sequencing provider).

### 🗄️ Raw Plotting Data

The folder also contains intermediate CSV files (`Fig1A_Data.csv` through `Fig3C_Data.csv`) so you can easily recreate or reformat the charts in external software like Origin or GraphPad Prism.

-----

> ### 💡 Pro-Tip: Interpreting the K-S Test Results in Fig 3B
>
> You will see a **$D$-statistic** and a **$p$-value** on the Fig 3B density plot. A smaller $D$ value indicates tighter alignment with the original population.
>
> **Don't panic if $p < 0.05$\!** In the context of GWAS sampling, our variance-weighted strategy *intentionally* enriches extreme phenotypes (the thick tails of the distribution). This slight statistical deviation from the original normal distribution is expected, highly beneficial for capturing maximum genetic variation, and ultimately improves the power to detect significant loci.

***

## 🧰 Auxiliary Tool: 2D Trait Geometric Selection (V3.0)

Included in this repository is an auxiliary script (`Alfalfa_Selection_V3.0.R`) designed for rapid, visualization-driven targeted selection. While the main GWAS pipeline focuses on multidimensional variance across the entire population, this supplementary tool isolates two specific agronomic traits (e.g., Plant Height and Multifoliate Score) to identify the absolute top performers.

### ✨ Overview & Features

This script calculates a combined performance index and generates a publication-ready contour plot to visualize the exact selection threshold.

* **Weighted Geometric Index:** Uses the formula `Index = (0.5 × Height) × (0.5 × Multifoliate)` to balance the contribution of both traits.
* **Top 50 Extraction:** Automatically isolates the top 50 individuals based on their calculated index scores.
* **Contour Visualization:** Generates a highly stylized 2D contour surface map. It overlays the total population scatter plot with the selection gradient, highlighting the selected top 50 plants and explicitly drawing the cutoff boundary line.
* **Smart Data Wrangling:** Automatically converts raw field-grid data (Rows L1-L50, Columns A-U) into a clean, analyzable long format.

### 🚀 How to Use the Auxiliary Tool

1. **Input Files:** Ensure your raw data files are in the working directory:
   * `Rawdata-PlantHeight_202504.csv`
   * `Rawdata-Multifoliate_202504.csv`
2. **Execution:** Source the script in R. It will automatically detect missing packages (like `metR` for contour generation) and install them.
3. **Outputs:** The script creates a timestamped folder (e.g., `Alfalfa_GeoViz_Pro_YYYYMMDD_HHMMSS`) containing:
   * `geometric_selection.png`: A high-resolution (300 DPI) visualization of the selection landscape.
   * `selection_list.csv`: The finalized list of the top 50 selected plants and their scores.
