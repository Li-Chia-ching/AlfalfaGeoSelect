# ============================================================================
# Script: alfalfa-multiyear-gwas-pipeline.R
# Description: Multi-year (Autumn 2025 vs Spring 2026) phenotypic analysis
#              and GWAS candidate selection pipeline for alfalfa S1 families
#
# Key features:
#   1. Feature engineering: growth difference and plasticity
#   2. Individual-level reaction norms
#   3. Trait correlation network analysis
#   4. Integrated selection: Pareto front ∩ MGIDI
#   5. Structured outputs for publication and reproducibility
# ============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

# ------------------------------
# Load required packages
# ------------------------------
required_packages <- c("dplyr", "ggplot2", "tidyr", "rPref", "metan", 
                       "qgraph", "factoextra", "ggrepel", "corrplot")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# ------------------------------
# 1. Data loading and feature engineering
# ------------------------------
load_and_engineer_data <- function(file_25, file_26) {
  
  # ---- File validation ----
  if (!file.exists(file_25)) stop("File not found: ", file_25)
  if (!file.exists(file_26)) stop("File not found: ", file_26)
  
  # ---- Read data ----
  df25 <- read.csv(file_25)
  df26 <- read.csv(file_26)
  
  # ---- Column validation ----
  required_25 <- c("Parent", "SampleID", "PlantHeight", "NodeCount", "Group")
  required_26 <- c("Group", "Matrix", "Plant_Height", "Internode", 
                   "Branch_Number", "Multifoliate_Score", "Death_Code")
  
  if (!all(required_25 %in% colnames(df25))) {
    stop("Missing columns in 2025 data")
  }
  if (!all(required_26 %in% colnames(df26))) {
    stop("Missing columns in 2026 data")
  }
  
  # ---- Standardize column names ----
  df25 <- df25 %>%
    rename(Family = Parent,
           Height_Autumn = PlantHeight,
           Node_Autumn = NodeCount)
  
  df26 <- df26 %>%
    rename(Family = Group,
           SampleID = Matrix,
           Height_Spring = Plant_Height,
           Node_Spring = Internode)
  
  # ---- Merge datasets ----
  merged_df <- full_join(df25, df26, by = c("Family", "SampleID"))
  
  # ---- Survival encoding ----
  merged_df <- merged_df %>%
    mutate(
      Death_Code = ifelse(is.na(Death_Code) | Death_Code == "", 0, as.numeric(Death_Code)),
      Survival_Binary = ifelse(Death_Code == 1, 0, 1)
    )
  
  # ---- Feature engineering ----
  merged_df <- merged_df %>%
    mutate(
      Growth_Diff = Height_Spring - Height_Autumn,
      Plasticity_Ratio = Height_Spring / Height_Autumn,
      Plasticity_Ratio = ifelse(is.infinite(Plasticity_Ratio), NA, Plasticity_Ratio)
    )
  
  return(merged_df)
}

# ------------------------------
# 2. Reaction norm plot (individual level)
# ------------------------------
plot_reaction_norms <- function(data, output_dir) {
  
  rn_data <- data %>%
    filter(Survival_Binary == 1) %>%
    drop_na(Height_Autumn, Height_Spring) %>%
    pivot_longer(cols = c(Height_Autumn, Height_Spring),
                 names_to = "Season", values_to = "Height")
  
  p <- ggplot(rn_data, aes(x = Season, y = Height, group = SampleID, color = Family)) +
    geom_line(alpha = 0.3) +
    geom_point(size = 1) +
    theme_bw() +
    labs(title = "Reaction Norms (Autumn → Spring)",
         x = "Season", y = "Plant Height")
  
  ggsave(file.path(output_dir, "plots", "reaction_norms.png"), p, width = 8, height = 6)
}

# ------------------------------
# 3. Trait correlation network
# ------------------------------
plot_correlation_network <- function(data, output_dir) {
  
  traits <- data %>%
    select(Height_Autumn, Node_Autumn, Height_Spring, Node_Spring,
           Branch_Number, Multifoliate_Score, Growth_Diff, Plasticity_Ratio) %>%
    na.omit()
  
  cor_matrix <- cor(traits)
  
  write.csv(cor_matrix, file.path(output_dir, "data", "correlation_matrix.csv"))
  
  png(file.path(output_dir, "plots", "trait_network.png"), width = 2000, height = 2000)
  qgraph(cor_matrix, layout = "spring")
  dev.off()
}

# ------------------------------
# 4. GWAS candidate selection
# ------------------------------
select_candidates <- function(data, output_dir, n_select = 200) {
  
  traits <- c("Height_Spring", "Branch_Number", "Multifoliate_Score", "Plasticity_Ratio")
  
  df <- data %>%
    filter(Survival_Binary == 1) %>%
    drop_na(all_of(traits)) %>%
    mutate(ID = make.unique(paste(Family, SampleID)))
  
  # ---- Pareto selection ----
  pref <- high(Height_Spring) * high(Branch_Number) * high(Multifoliate_Score) * high(Plasticity_Ratio)
  pareto <- psel(df, pref, top = n_select * 1.5)
  
  # ---- MGIDI (distance to ideotype) ----
  mat <- scale(df[, traits])
  ideotype <- apply(mat, 2, max)
  dist <- sqrt(rowSums((mat - matrix(ideotype, nrow(mat), length(ideotype), byrow = TRUE))^2))
  
  df$MGIDI <- dist
  
  mgidi_top <- df %>% arrange(MGIDI) %>% head(n_select * 1.5)
  
  # ---- Intersection ----
  final <- pareto %>%
    inner_join(mgidi_top, by = "ID") %>%
    arrange(MGIDI) %>%
    head(n_select)
  
  write.csv(final, file.path(output_dir, "data", "gwas_candidates.csv"), row.names = FALSE)
  
  return(final)
}

# ------------------------------
# 5. Main pipeline
# ------------------------------
main <- function(file_25, file_26, n_select = 200) {
  
  out_dir <- paste0("gwas_selection_", Sys.Date())
  dir.create(out_dir, showWarnings = FALSE)
  dir.create(file.path(out_dir, "plots"), showWarnings = FALSE)
  dir.create(file.path(out_dir, "data"), showWarnings = FALSE)
  
  cat("Running pipeline...\n")
  
  df <- load_and_engineer_data(file_25, file_26)
  write.csv(df, file.path(out_dir, "data", "merged_data.csv"), row.names = FALSE)
  
  plot_reaction_norms(df, out_dir)
  plot_correlation_network(df, out_dir)
  
  candidates <- select_candidates(df, out_dir, n_select)
  
  cat("Done. Selected:", nrow(candidates), "samples\n")
}

# ------------------------------
# Run
# ------------------------------
main("file_25.csv", "file_26.csv", 200)
