# =============================================================================
# Pre-GWAS Phenomic Pipeline: Spatial Correction & Phenotypic Clustering
# Objective: Remove field spatial noise and explore macro phenotypic structure 
#            before downstream Neyman sampling.
# Input: data_2026.csv (Raw)
# Output: data_2026_Corrected.csv (Spatial-adjusted), Diversity Visualizations
# =============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

library(dplyr)
library(tidyr)
library(readr)
library(sommer)
library(factoextra)
library(cluster)
library(ggplot2)
library(showtext)

# ------------------- 1. Initialization & Data Loading -------------------
showtext_auto() 
if (.Platform$OS.type == "windows") {
  font_add("SimSun", "simsun.ttc") 
} else {
  font_add("SimSun", "STSong.ttf") 
}

out_dir <- paste0("Phenotypic_Diversity_", format(Sys.Date(), "%Y%m%d"))
dir.create(out_dir, showWarnings = FALSE)

academic_theme <- theme_bw(base_family = "SimSun") + 
  theme(plot.title = element_text(face = "bold", size = 14))

cat("Loading raw 2026 dataset...\n")
df_raw <- read_csv("data_2026.csv", show_col_types = FALSE) %>%
  filter(!is.na(Group) & trimws(Group) != "")

# ------------------- 2. Matrix Parsing & Spatial Coordinate Mapping -----------
# Convert alphanumeric Matrix identifiers (e.g., A01, AA47) into strict 2D numeric coordinates
col_letter_to_numeric <- function(letters_vec) {
  sapply(letters_vec, function(x) {
    if (is.na(x) || x == "") return(NA)
    chars <- strsplit(toupper(x), "")[[1]]
    sum(sapply(seq_along(chars), function(i) match(chars[i], LETTERS) * (26^(length(chars) - i))))
  })
}

df_spatial <- df_raw %>%
  mutate(
    Col_Num = col_letter_to_numeric(stringr::str_extract(Matrix, "^[A-Za-z]+")),
    Row_Num = as.numeric(stringr::str_extract(Matrix, "\\d+$"))
  ) %>%
  mutate(Row_Num = ifelse(is.na(Row_Num), 0, Row_Num), 
         Col_Num = ifelse(is.na(Col_Num), 0, Col_Num))

# ------------------- 3. Spatial BLUP Extraction (Environmental Filtering) -----
cat("Applying 2D spatial splines to extract noise-free individual phenotypes...\n")
trait_vars <- c("Plant_Height", "Internode", "Branch_Number", "Multifoliate_Score")

# Create a copy to store corrected phenotypic values
df_corrected <- df_spatial

for (trait in trait_vars) {
  cat(" - Processing architecture and yield trait:", trait, "\n")
  
  # Create temporary modeling dataset: retain only non-missing observations with valid coordinates
  df_model <- df_spatial %>%
    filter(!is.na(.data[[trait]])) %>%
    mutate(Row_Num = as.numeric(Row_Num), Col_Num = as.numeric(Col_Num))
  
  # If insufficient data, skip model fitting and retain original values
  if (nrow(df_model) < 30) {
    cat("   WARNING: Too few non-missing observations for", trait, ". Skipping spatial correction.\n")
    df_corrected[[trait]] <- df_spatial[[trait]]
    next
  }
  
  # Fit spatial model
  mix_mod <- tryCatch({
    mmer(fixed = as.formula(paste(trait, "~ 1")),
         random = ~ Group,
         rcov = ~ vsr(spl2D(Row_Num, Col_Num)),
         data = df_model,
         date.warning = FALSE)
  }, error = function(e) {
    cat("   WARNING: Spatial spline failed for", trait, "- falling back to standard LMM.\n")
    lme4::lmer(as.formula(paste(trait, "~ 1 + (1 | Group)")), data = df_model)
  })
  
  # Extract corrected values (only for modeled rows)
  if (inherits(mix_mod, "mmer")) {
    corrected_model <- fitted(mix_mod)$data[[trait]]
  } else {
    corrected_model <- predict(mix_mod)
  }
  
  # Backfill: initialize with original values, then overwrite modeled rows with corrected estimates
  df_corrected[[trait]] <- df_spatial[[trait]]
  df_corrected[[trait]][df_model %>% rownames() %>% as.integer()] <- corrected_model
}

# ------------------- 4. Phenotypic Diversity & Clustering ---------------------
cat("Exploring macro phenotypic diversity...\n")

# Extract trait matrix for clustering (base R subset to avoid select parsing issues)
pheno_mat <- df_corrected[, trait_vars, drop = FALSE]
pheno_scaled <- scale(pheno_mat)

# Remove rows with NaN values caused by missing data
na_rows <- apply(pheno_scaled, 1, function(x) any(is.na(x)))
if (any(na_rows)) {
  cat(" - Removing", sum(na_rows), "rows with missing values after scaling.\n")
  pheno_scaled <- pheno_scaled[!na_rows, , drop = FALSE]
  df_corrected_clean <- df_corrected[!na_rows, ]
} else {
  df_corrected_clean <- df_corrected
}

# 4.1 Determine optimal number of clusters K (with biological prior constraint)
# Limit k.max to 5 to avoid capturing trivial within-family variation
sil_res <- fviz_nbclust(pheno_scaled, kmeans, method = "silhouette", k.max = 5)
p_opt_k <- sil_res + academic_theme + labs(title = "A. Optimal number of phenotypic clusters")

best_k <- as.numeric(as.character(sil_res$data$clusters[which.max(sil_res$data$y)]))

# Enforced biological adjustment:
# Even within 1–5, if K is overly fragmented (e.g., 4 or 5),
# constrain to 3 macro-level gradients if biologically justified
if (length(best_k) == 0 || is.na(best_k) || best_k < 2) {
  best_k <- 3
} else if (best_k > 4) {
  cat(" - Silhouette suggested K =", best_k, "but capping at 3 for macro-gradient sampling.\n")
  best_k <- 3
}
cat(" - Final macro cluster number used:", best_k, "\n")

# 4.2 K-means clustering
set.seed(2026)
km_res <- kmeans(pheno_scaled, centers = best_k, nstart = 25)
df_corrected_clean$Phenotypic_Cluster <- paste0("Cluster_", km_res$cluster)

# Map cluster labels back to original dataset (index-aligned)
df_corrected$Phenotypic_Cluster <- NA
df_corrected$Phenotypic_Cluster[!na_rows] <- df_corrected_clean$Phenotypic_Cluster

# 4.3 PCA visualization of clustering structure
pca_res <- prcomp(pheno_scaled, center = FALSE, scale. = FALSE)
pca_df <- as.data.frame(pca_res$x[, 1:2])
pca_df$Cluster <- df_corrected_clean$Phenotypic_Cluster

p_cluster <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster, fill = Cluster)) +
  geom_point(shape = 21, size = 2.5, alpha = 0.7, stroke = 0.5) +
  stat_ellipse(geom = "polygon", alpha = 0.15, linetype = "dashed", na.rm = TRUE) +
  scale_color_viridis_d(option = "turbo") +
  scale_fill_viridis_d(option = "turbo") +
  labs(title = sprintf("B. Multidimensional Phenotypic Structure (K = %d)", best_k),
       subtitle = "Spatial-corrected plant architecture traits",
       x = "PC1", y = "PC2") +
  academic_theme

# Export diversity plots
library(patchwork)
ggsave(file.path(out_dir, "01_Phenotypic_Diversity_Structure.pdf"), 
       p_opt_k / p_cluster, width = 8, height = 10, dpi = 600)

# -------------------------------------------------------------------------
# [NEW] Export underlying plotting data for Fig A and Fig B to CSV
# -------------------------------------------------------------------------
cat("Exporting plotting data for external software (Origin/GraphPad)...\n")

# 1. Export underlying data for Fig A (silhouette scores)
# sil_res$data contains two columns: clusters (K) and y (silhouette score)
write_csv(sil_res$data, file.path(out_dir, "PlotData_FigA_Silhouette_Scores.csv"))

# 2. Export underlying data for Fig B (PCA scatter + cluster structure)
# pca_df already contains PC1, PC2, and Cluster labels
# Append plant IDs and corrected phenotypes for downstream validation
plot_data_pca <- pca_df %>%
  mutate(
    Plant_ID = df_corrected_clean$Matrix, 
    Family_ID = df_corrected_clean$Group
  ) %>%
  bind_cols(df_corrected_clean %>% select(all_of(trait_vars))) %>%
  # Reorder columns: place IDs and cluster labels first
  select(Family_ID, Plant_ID, Cluster, PC1, PC2, everything())

write_csv(plot_data_pca, file.path(out_dir, "PlotData_FigB_PCA_Scatter.csv"))

cat(" - Plotting data exported successfully to", out_dir, "\n")

# ------------------- 5. Export for Downstream Pipeline -------------------
cat("Exporting spatial-corrected dataset for downstream sampling...\n")

# Remove auxiliary coordinate columns, retain all other variables (including Phenotypic_Cluster)
keep_cols <- setdiff(names(df_corrected), c("Col_Num", "Row_Num"))
df_export <- df_corrected[, keep_cols, drop = FALSE]

# Save corrected dataset
write_csv(df_export, "data_2026_Corrected.csv")
write_csv(df_export, file.path(out_dir, "data_2026_Corrected_with_Clusters.csv"))

cat("Pre-processing complete. Field noise removed. You can now run the downstream script using 'data_2026_Corrected.csv'.\n")
