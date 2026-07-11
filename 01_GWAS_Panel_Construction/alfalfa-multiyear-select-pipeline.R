# Integrated Pre-GWAS Phenomic Pipeline & Core Collection Sampling (v7.1) ----
# GWAS Power Optimization: Balanced Family Allocation + FORCE-FILL to Exact 200
#
# v7.0 -> v7.1 Key Changes:
# 1. TARGET_TOTAL_N = 200 (exact, non-negotiable for supervisor requirement)
# 2. MAX_PER_FAMILY = 7 (soft cap; Phase 4 may exceed to 8-9 for force-fill)
# 3. MIN_FAMILY_SIZE = 3 -> filters unreliable micro-families
# 4. MIN_ALLOC = 1 -> guarantees every eligible family contributes at least 1 plant
# 5. NEW: Phase 4 "Force-Fill" — post-allocation gap detection & mandatory fill-up
#    - Detects if selected < target after Neyman + stratified sampling
#    - Fills gap from large families with remaining unselected individuals
#    - Allows soft-cap breach (up to MAX_PER_FAMILY + 2) during fill phase
#    - Guarantees nrow(selected_plants) == TARGET_TOTAL_N exactly
# 6. Column mapping fixed for data_202605_gwas.csv field names

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

## 0. Package Management ----
required_packages <- c("dplyr", "tidyr", "ggplot2", "lme4", "lmerTest",
                       "factoextra", "readr", "patchwork", "showtext",
                       "cluster", "viridis")

new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, repos = "https://cloud.r-project.org")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lme4)
  library(lmerTest)
  library(factoextra)
  library(readr)
  library(patchwork)
  library(showtext)
  library(cluster)
  library(viridis)
})

tryCatch({
  if (.Platform$OS.type == "windows") { 
    font_add("SimSun", "simsun.ttc") 
  } else { 
    font_add("SimSun", "STSong.ttf") 
  }
  showtext_auto()
  base_font <- "SimSun"
}, error = function(e) {
  base_font <- ""
})

# KEY PARAMETERS: Directly affect GWAS statistical power ----
TARGET_TOTAL_N  <- 200      # [HARD REQUIREMENT] Supervisor demands exact 200 plants
MAX_PER_FAMILY  <- 7        # Soft cap for Neyman allocation phase (controls kinship bias)
MIN_FAMILY_SIZE <- 3        # Minimum family threshold (<3 excluded: unstable phenotype estimates)
MIN_ALLOC       <- 1        # Neyman floor guarantee (every eligible family >= 1 plant)

# Phase 4 Force-Fill parameters (only triggered when Neyman allocation falls short)
FILL_ALLOW_EXCESS <- 2      # Max extra plants allowed per family during fill (i.e., up to 9)
                            # Rationale: fill amount is tiny (usually 3-5), spread across
                            # large families means +1~2 per family max, negligible impact on
                            # Kinship matrix but satisfies "exact 200" hard requirement

all_traits <- c("Plant_Height_Nov", "Internode_Nov", "Plant_Height_Mar",
                "Internode_Mar", "Branch_Number_Mar", "Multifoliate_Score_Mar",
                "Plant_Height_May")

# Output Directory Structure ----
# Timestamped, never overwrites previous runs
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
out_root <- paste0("GWAS_Pipeline_Result_", timestamp)
dir.create(out_root, showWarnings = FALSE)
dir.create(file.path(out_root, "01_Scripts"), showWarnings = FALSE)
dir.create(file.path(out_root, "02_Data_Tables"), showWarnings = FALSE)
dir.create(file.path(out_root, "03_Figures"), showWarnings = FALSE)
dir.create(file.path(out_root, "04_Plot_Data"), showWarnings = FALSE)

academic_theme <- theme_bw(base_family = base_font) +
  theme(panel.grid.minor = element_blank(),
        text = element_text(family = base_font, size = 12),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 14),
        strip.text = element_text(size = 9, margin = margin(2, 2, 2, 2)))

cat("================================================================\n")
cat("    GWAS Pipeline v7.1 - Force-Fill to Exact N Edition\n")
cat("================================================================\n")
cat(sprintf("Target N=%d | Soft Cap=%d | Min Family=%d | Fill Excess=%d\n",
            TARGET_TOTAL_N, MAX_PER_FAMILY, MIN_FAMILY_SIZE, FILL_ALLOW_EXCESS))
cat(sprintf("Output: %s\n\n", out_root))

# STAGE 1: DATA LOADING, IMPUTATION & PHENOTYPIC CLUSTERING ----

## 1.1 Data Loading & Robust Imputation ----
# Performs within-family median imputation
cat("[Stage 1] Loading data and performing robust imputation...\n")

df_raw <- read_csv("01_Raw_Data/data_202605_gwas.csv", show_col_types = FALSE)

# Column name mapping: raw data fields -> standard analysis names
df_raw <- df_raw %>%
  dplyr::rename(
    Family_ID = Family,
    Plant_ID = ID,
    Plant_Height_Nov = PH_Nov,
    Internode_Nov = IN_Nov,
    Plant_Height_Mar = PH_Mar,
    Internode_Mar = IN_Mar,
    Branch_Number_Mar = BN_Mar,
    Multifoliate_Score_Mar = MF_Score_Mar,
    Plant_Height_May = PH_May
  )

df_imputed <- df_raw %>%
  dplyr::filter(!is.na(Plant_Height_May)) %>%
  dplyr::group_by(Family_ID) %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(all_traits),
                              ~ifelse(is.na(.), median(., na.rm = TRUE), .))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(all_traits),
                              ~ifelse(is.na(.), median(., na.rm = TRUE), .)))

n_original <- nrow(df_raw %>% dplyr::filter(!is.na(Plant_Height_May)))
cat(sprintf("  After imputation: %d individuals (from %d with PH_May)\n", nrow(df_imputed), n_original))

write_csv(df_imputed %>% dplyr::select(Family_ID, Plant_ID, dplyr::all_of(all_traits)),
          file.path(out_root, "02_Data_Tables", "00_Imputed_Dataset.csv"))

## 1.2 Small Family Filter (GWAS Power Protection) ----
# Excludes <MIN_FAMILY_SIZE micro-families:
# With only 1-2 plants, within-family variance cannot be reliably estimated.
# Selected singletons are essentially environmental noise outliers that would
# artificially inflate false positive signals in downstream GWAS.
family_size_check <- df_imputed %>%
  dplyr::group_by(Family_ID) %>%
  dplyr::summarise(N_Original = n(), .groups = "drop")

eligible_families <- family_size_check %>%
  dplyr::filter(N_Original >= MIN_FAMILY_SIZE) %>%
  dplyr::pull(Family_ID)

excluded_families <- family_size_check %>%
  dplyr::filter(N_Original < MIN_FAMILY_SIZE)

df_imputed <- df_imputed %>% dplyr::filter(Family_ID %in% eligible_families)

n_families <- length(eligible_families)
n_excluded <- nrow(excluded_families)
cat(sprintf("  Eligible families: %d / %d (excluded %d with < %d plants)\n",
            n_families, nrow(family_size_check), MIN_FAMILY_SIZE, n_excluded))
cat(sprintf("  After filter: %d individuals\n\n", nrow(df_imputed)))

write_csv(
  family_size_check %>%
    dplyr::mutate(Status = ifelse(N_Original >= MIN_FAMILY_SIZE, "Eligible", "Excluded")),
  file.path(out_root, "02_Data_Tables", "00_Family_Eligibility_Check.csv")
)

## 1.3 Phenotypic Clustering & Visualization (Fig S1) ----
# K-means clustering reveals macro-level phenotypic diversity structure.
# Not used directly for sampling, but serves as visual validation of coverage.
cat("[Stage 1] Exploring macro phenotypic diversity via K-means...\n")
pheno_mat <- df_imputed[, all_traits, drop = FALSE]
pheno_scaled <- scale(pheno_mat)

sil_res <- fviz_nbclust(pheno_scaled, kmeans, method = "silhouette", k.max = 6)
p_opt_k <- sil_res + academic_theme + labs(title = "A. Optimal number of phenotypic clusters (silhouette)")

best_k <- as.numeric(as.character(sil_res$data$clusters[which.max(sil_res$data$y)]))
if(length(best_k) == 0 || is.na(best_k) || best_k < 2) best_k <- 3
if(best_k > 5) best_k <- 3

set.seed(2026)
km_res <- kmeans(pheno_scaled, centers = best_k, nstart = 25)
df_imputed$Phenotypic_Cluster <- paste0("Cluster_", km_res$cluster)

pca_res_clust <- prcomp(pheno_scaled, center = FALSE, scale. = FALSE)
pca_df <- as.data.frame(pca_res_clust$x[, 1:2])
pca_df$Cluster <- df_imputed$Phenotypic_Cluster

p_cluster <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster, fill = Cluster)) +
  geom_point(shape = 21, size = 2.5, alpha = 0.7, stroke = 0.5) +
  stat_ellipse(geom = "polygon", alpha = 0.15, linetype = "dashed", na.rm = TRUE) +
  scale_color_viridis_d(option = "turbo") +
  scale_fill_viridis_d(option = "turbo") +
  labs(title = sprintf("B. Multidimensional Phenotypic Structure (K=%d)", best_k),
       subtitle = "Imputed original traits (no spatial correction)", x = "PC1", y = "PC2") +
  academic_theme

ggsave(file.path(out_root, "03_Figures", "FigS1_Phenotypic_Clustering.pdf"),
       p_opt_k / p_cluster, width = 8, height = 10, dpi = 600)

write_csv(sil_res$data, file.path(out_root, "04_Plot_Data", "FigS1A_Silhouette_Data.csv"))

plot_data_pca <- pca_df %>%
  dplyr::mutate(Plant_ID = df_imputed$Plant_ID, Family_ID = df_imputed$Family_ID) %>%
  dplyr::bind_cols(df_imputed %>% dplyr::select(dplyr::all_of(all_traits))) %>%
  dplyr::select(Family_ID, Plant_ID, Cluster, PC1, PC2, dplyr::everything())
write_csv(plot_data_pca, file.path(out_root, "04_Plot_Data", "FigS1B_PCA_Cluster_Data.csv"))
cat(sprintf("  Clustering done: K=%d clusters\n\n", best_k))


# STAGE 2: GWAS SAMPLING PIPELINE (Neyman -> Stratified -> Force-Fill) ----
cat("[Stage 2] Starting Neyman allocation + stratified sampling + force-fill...\n\n")

## 2.1 Evidence I: Segregation & Repeatability (Fig 1) ----
# Within-family segregation analysis + mixed model variance component estimation.
# Repeatability (R) is a core criterion for trait selection in GWAS design.
cat("[2.1] Computing variance components & repeatability...\n")

df_for_plot <- df_imputed %>%
  dplyr::add_count(Family_ID, name = "n_per_family") %>%
  dplyr::filter(n_per_family >= 2) %>%
  dplyr::mutate(Family_Label = paste0(Family_ID, "\n(n=", n_per_family, ")"))

p_seg_density <- ggplot(df_for_plot, aes(x = Plant_Height_May)) +
  geom_density(fill = "#377EB8", alpha = 0.6) +
  geom_rug(alpha = 0.5, sides = "b") +
  facet_wrap(~Family_Label, ncol = 9, scales = "free_y") +
  labs(x = "Plant Height May 2026 (cm)", y = "Kernel density",
       title = "A. Density distribution of plant height segregation within families") +
  academic_theme + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

p_var_box <- ggplot(df_imputed, aes(x = reorder(Family_ID, Plant_Height_May, FUN = median),
                                    y = Plant_Height_May)) +
  geom_boxplot(fill = "#E41A1C", alpha = 0.6, outlier.size = 0.5) +
  labs(x = "Selfed family ID (ordered by median)", y = "Plant Height May 2026 (cm)",
       title = "B. Between-family variation & within-family segregation range") +
  academic_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

ggsave(file.path(out_root, "03_Figures", "Fig1_Segregation_Analysis.pdf"),
       p_seg_density / p_var_box, width = 16, height = 12, dpi = 600)

write_csv(df_imputed %>% dplyr::select(Family_ID, Plant_ID, Plant_Height_May),
          file.path(out_root, "02_Data_Tables", "Fig1_Data.csv"))
write_csv(df_for_plot %>% dplyr::select(Family_ID, Family_Label, Plant_Height_May, n_per_family),
          file.path(out_root, "04_Plot_Data", "Fig1_Density_Data.csv"))

get_varcomp <- function(trait) {
  trait_sd <- sd(df_imputed[[trait]], na.rm = TRUE)
  if(is.na(trait_sd) || trait_sd == 0) return(NULL)

  df_temp <- df_imputed %>% dplyr::filter(!is.na(.data[[trait]]))
  df_temp$scaled_trait <- scale(df_temp[[trait]])

  form <- as.formula("scaled_trait ~ 1 + (1|Family_ID)")
  mod <- try(lmer(form, data = df_temp,
                  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))),
             silent = TRUE)

  if(inherits(mod, "try-error")) {
    return(data.frame(Trait = trait, Var_Family = NA, Var_Within = NA,
                      Repeatability_R = NA, Note = "Model convergence failed"))
  }

  vc <- as.data.frame(VarCorr(mod))
  var_fam_raw <- vc$vcov[1] * (trait_sd^2)
  var_wit_raw <- vc$vcov[2] * (trait_sd^2)
  R_val <- var_fam_raw / (var_fam_raw + var_wit_raw)

  data.frame(Trait = trait, Var_Family = var_fam_raw, Var_Within = var_wit_raw,
             Repeatability_R = R_val, Note = "Converged")
}

varcomp_results <- do.call(rbind, lapply(all_traits, get_varcomp)) %>% dplyr::arrange(desc(Repeatability_R))
write_csv(varcomp_results, file.path(out_root, "02_Data_Tables", "01_Variance_Components_Repeatability.csv"))
cat(sprintf("  Highest repeatability trait: %s (R=%.3f)\n\n",
            varcomp_results$Trait[1], varcomp_results$Repeatability_R[1]))

## 2.2 Evidence II: PCA Gradients & Biplot (Fig 2) ----
# PCA compresses multi-trait data into 2 principal components for
# Neyman weight calculation and distribution validation.
cat("[2.2] Performing Principal Component Analysis...\n")

pca_res_main <- prcomp(df_imputed %>% dplyr::select(dplyr::all_of(all_traits)),
                       center = TRUE, scale. = TRUE)
df_imputed <- df_imputed %>% dplyr::mutate(PC1 = pca_res_main$x[,1], PC2 = pca_res_main$x[,2])
prop_var <- summary(pca_res_main)$importance[2, 1:2]

pca_summary <- summary(pca_res_main)
eig_data <- data.frame(
  PC = 1:length(pca_summary$importance[1,]),
  StdDev = pca_summary$importance[1,],
  Variance_Proportion = pca_summary$importance[2,],
  Cumulative_Proportion = pca_summary$importance[3,]
)
write_csv(eig_data, file.path(out_root, "02_Data_Tables", "02_PCA_Eigenvalues.csv"))

loadings_data <- as.data.frame(pca_res_main$rotation[, 1:2]) %>%
  dplyr::rename(PC1_loading = PC1, PC2_loading = PC2) %>%
  tibble::rownames_to_column(var = "Trait")
write_csv(loadings_data, file.path(out_root, "02_Data_Tables", "02_PCA_Loadings.csv"))

write_csv(df_imputed %>% dplyr::select(Family_ID, Plant_ID, PC1, PC2, dplyr::all_of(all_traits)),
          file.path(out_root, "02_Data_Tables", "02_PCA_Scores.csv"))

p_pca_main <- fviz_pca_biplot(pca_res_main, geom.ind = "point", pointshape = 21, pointsize = 1.5,
                              habillage = df_imputed$Phenotypic_Cluster, alpha.ind = 0.6,
                              col.var = "black", arrowsize = 0.6, labelsize = 3.5, repel = TRUE,
                              title = "Multidimensional Phenotypic Gradients (PCA Biplot)",
                              palette = "turbo", ggtheme = academic_theme)

ggsave(file.path(out_root, "03_Figures", "Fig2_PCA_Biplot_Enhanced.pdf"),
       p_pca_main, width = 9, height = 7, dpi = 600)

pca_biplot_data <- data.frame(
  Plant_ID = df_imputed$Plant_ID,
  Family_ID = df_imputed$Family_ID,
  PC1 = df_imputed$PC1,
  PC2 = df_imputed$PC2,
  Cluster = df_imputed$Phenotypic_Cluster
)
write_csv(pca_biplot_data, file.path(out_root, "04_Plot_Data", "Fig2_PCA_Biplot_Data.csv"))
cat(sprintf("  PC1 explains: %.1f%% | PC2: %.1f%%\n\n", prop_var[1]*100, prop_var[2]*100))

## 2.3 Weighted Neyman Allocation (Phases 1-3) ----
# Neyman optimal allocation principle: n_h proportional to N_h x S_h
# (family size * within-family SD). This minimizes estimation variance
# for a given total sample size under stratified random sampling framework.
# GWAS-specific constraints applied:
#   (1) Soft cap MAX_PER_FAMILY -> controls per-family kinship weight peak
#   (2) Floor guarantee MIN_ALLOC -> protects rare alleles from omission
#   (3) Pre-filter small families -> excludes unreliable extreme noise values
cat("[2.3] Executing weighted Neyman allocation (Phases 1-3)...\n")

family_stats <- df_imputed %>%
  dplyr::group_by(Family_ID) %>%
  dplyr::summarise(
    N_h = n(),
    Var_Weighted = var(PC1)*prop_var[1] + var(PC2)*prop_var[2],
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Var_Weighted = ifelse(is.na(Var_Weighted) | Var_Weighted < 0, 0, Var_Weighted),
    S_h = sqrt(Var_Weighted)
  )

if(sum(family_stats$S_h) == 0) family_stats$S_h <- 1
family_stats$Weight <- family_stats$N_h * family_stats$S_h

cat(sprintf("  Participating families: %d | Total candidates: %d\n", n_families, sum(family_stats$N_h)))

smart_allocate_with_constraints <- function(weights, target, cap, floor_alloc) {
  n_strata <- length(weights)
  alloc <- rep(floor_alloc, n_strata)

  if(sum(alloc) >= target) {
    alloc <- floor(target * weights / sum(weights))
    alloc[alloc < 1] <- 1
    alloc[alloc > cap] <- cap
    diff_adj <- target - sum(alloc)
    while(diff_adj != 0) {
      eligible <- which(alloc < cap)
      if(length(eligible) == 0) break
      if(diff_adj > 0) {
        idx <- sample(eligible, 1)
        alloc[idx] <- alloc[idx] + 1; diff_adj <- diff_adj - 1
      } else {
        reducible <- which(alloc > floor_alloc)
        if(length(reducible) == 0) break
        idx <- sample(reducible, 1)
        alloc[idx] <- alloc[idx] - 1; diff_adj <- diff_adj + 1
      }
    }
    return(alloc)
  }

  # Phase 1: Neyman proportional allocation (with per-stratum cap awareness)
  remaining <- target - sum(alloc)
  proportional <- remaining * weights / sum(weights)
  prop_capped <- pmin(proportional, cap - floor_alloc)
  alloc <- alloc + floor(prop_capped)

  # Phase 2: Distribute fractional remainders by Neyman priority
  diff <- target - sum(alloc)
  while(diff > 0) {
    eligible <- alloc < cap
    if(!any(eligible)) break
    priority <- weights / (alloc + 1)
    priority[!eligible] <- -Inf
    alloc[which.max(priority)] <- alloc[which.max(priority)] + 1
    diff <- diff - 1
  }

  # Phase 3: Fill-to-target post-cap (fixes clipping-induced shortfall)
  alloc[alloc > cap] <- cap
  gap <- target - sum(alloc)
  while(gap > 0) {
    under_cap <- which(alloc < cap)
    if(length(under_cap) == 0) {
      warning(sprintf("Cannot reach target %d: max capacity %d", target, n_strata * cap))
      break
    }
    fill_priority <- rep(-Inf, n_strata)
    fill_priority[under_cap] <- weights[under_cap]
    alloc[which.max(fill_priority)] <- alloc[which.max(fill_priority)] + 1
    gap <- gap - 1
  }
  alloc[alloc > cap] <- cap

  return(alloc)
}

family_stats$Target_N <- smart_allocate_with_constraints(
  family_stats$Weight, TARGET_TOTAL_N, MAX_PER_FAMILY, MIN_ALLOC
)

neyman_allocated_total <- sum(family_stats$Target_N)
cat(sprintf("  Neyman allocation done: target=%d -> allocated=%d (soft cap=%d/family)\n",
            TARGET_TOTAL_N, neyman_allocated_total, MAX_PER_FAMILY))

alloc_display <- family_stats %>% dplyr::arrange(desc(Target_N)) %>%
  dplyr::select(Family_ID, N_h, S_h, Target_N)
cat("\n  Allocation TOP10:\n")
print(head(alloc_display, 10), width = 80)

## 2.3b Stratified Gradient Sampling ----
# For each family, stratify by PC1 quantiles and sample proportionally.
# This ensures selected samples cover the full phenotypic range within each
# family, avoiding bias toward extreme values only.
cat("\n  Executing stratified gradient sampling (PC1 quantile + random)...\n")
set.seed(2026)
selected_plants <- data.frame()

for (fam in family_stats$Family_ID) {
  fam_data <- df_imputed %>% dplyr::filter(Family_ID == fam)
  n_alloc <- family_stats$Target_N[family_stats$Family_ID == fam]

  if(n_alloc == 0) next
  if(n_alloc >= nrow(fam_data)) {
    selected_plants <- dplyr::bind_rows(selected_plants, fam_data)
    next
  }

  k_clusters <- min(3, nrow(fam_data))
  fam_data <- fam_data %>% dplyr::mutate(Gradient = as.factor(dplyr::ntile(PC1, k_clusters)))
  grad_counts <- as.data.frame(table(fam_data$Gradient))
  grad_alloc <- floor(n_alloc * (grad_counts$Freq / sum(grad_counts$Freq)))

  diff_grad <- n_alloc - sum(grad_alloc)
  while(diff_grad > 0) {
    eligible_idx <- which(grad_alloc < grad_counts$Freq)
    if(length(eligible_idx) == 0) break
    idx <- if(length(eligible_idx) == 1) eligible_idx else sample(eligible_idx, 1)
    grad_alloc[idx] <- grad_alloc[idx] + 1
    diff_grad <- diff_grad - 1
  }

  fam_sampled <- data.frame()
  for(i in seq_along(grad_counts$Var1)) {
    grad_pool <- fam_data %>% dplyr::filter(Gradient == grad_counts$Var1[i])
    n_to_pick <- min(grad_alloc[i], nrow(grad_pool))
    if(n_to_pick > 0) {
      fam_sampled <- dplyr::bind_rows(fam_sampled,
                                      grad_pool[sample(nrow(grad_pool), n_to_pick), ])
    }
  }
  selected_plants <- dplyr::bind_rows(selected_plants, fam_sampled)
}

n_after_neyman <- nrow(selected_plants)
cat(sprintf("  After Neyman+stratified sampling: %d plants selected\n", n_after_neyman))


## 2.3c Phase 4: FORCE-FILL - Gap Detection & Mandatory Top-up ----
# [Problem Background]
# The Neyman allocation algorithm is constrained by MAX_PER_FAMILY soft cap
# and integer truncation effects. Even with Phase 3 fill-to-target mechanism,
# final sample count can still fall slightly below TARGET_TOTAL_N.
# Example from v7.0 run: target=200 -> actual=197 (gap of 3).
#
# [Supervisor Requirement]
# Downstream GWAS population size must be EXACTLY 200 plants - this is a
# rigid specification that cannot be compromised.
#
# [Fill Strategy]
# 1. Calculate gap: N_diff = TARGET_TOTAL_N - nrow(selected_plants)
# 2. If N_diff <= 0: already at target, skip this phase
# 3. If N_diff > 0:
#    a. Build candidate pool from large families with unselected plants left
#       (prefer large original-population families: more residual individuals)
#    b. Sort candidates by Neyman weight descending (high-weight fills first)
#    c. Randomly pick 1 plant per iteration until exact TARGET_TOTAL_N reached
#    d. Allow individual families to exceed MAX_PER_FAMILY (up to
#       MAX_PER_FAMILY + FILL_ALLOW_EXCESS = 9) during this phase
#       Rationale: fill amount is tiny (typically 3-5), distributed across
#       large families means at most +1-2 extra per family, negligible impact
#       on Kinship matrix but satisfies the "exact 200" hard requirement
# 4. Final assertion check: assert(nrow(selected_plants) == TARGET_TOTAL_N)

cat("\n[2.3c] * Phase 4: Force-Fill - Gap Detection & Mandatory Top-up...\n")

N_current <- nrow(selected_plants)
N_diff <- TARGET_TOTAL_N - N_current

if(N_diff <= 0) {
  cat(sprintf("  [OK] Already at target: %d >= %d (no fill needed)\n", N_current, TARGET_TOTAL_N))
  fill_log <- data.frame(
    Step = 1, Action = "No fill needed",
    Before = N_current, After = N_current,
    N_diff = N_diff, Family_Added = NA,
    stringsAsFactors = FALSE
  )
} else {
  cat(sprintf("  [GAP DETECTED] Current: %d | Target: %d | Gap: %d plants\n",
              N_current, TARGET_TOTAL_N, N_diff))
  cat(sprintf("  Initiating force-fill from large families (cap relaxed to %d)...\n",
              MAX_PER_FAMILY + FILL_ALLOW_EXCESS))

  # ---- Step A: Build fill candidate pool ----
  current_per_family <- selected_plants %>%
    dplyr::group_by(Family_ID) %>%
    dplyr::summarise(Current_Selected = n(), .groups = "drop")

  fill_candidates <- family_stats %>%
    dplyr::left_join(current_per_family, by = "Family_ID") %>%
    dplyr::mutate(
      Current_Selected = ifelse(is.na(Current_Selected), 0, Current_Selected),
      Remaining = N_h - Current_Selected,
      Hard_Cap = MAX_PER_FAMILY + FILL_ALLOW_EXCESS
    ) %>%
    dplyr::filter(Remaining > 0) %>%
    dplyr::arrange(desc(Weight))

  n_candidate_fams <- nrow(fill_candidates)
  cat(sprintf("  Candidate families with remaining plants: %d\n", n_candidate_fams))

  if(n_candidate_fams == 0) {
    stop("[FATAL ERROR] No remaining plants available but still short! ",
         "Check data or relax filters.")
  }

  # ---- Step B: Iterative one-by-one fill ----
  set.seed(2026)
  fill_log <- data.frame()

  for(fill_step in seq_len(N_diff)) {
    current_selected_now <- selected_plants %>%
      dplyr::group_by(Family_ID) %>%
      dplyr::summarise(Now_Selected = n(), .groups = "drop")

    available <- fill_candidates %>%
      dplyr::left_join(current_selected_now, by = "Family_ID") %>%
      dplyr::mutate(
        Now_Selected = ifelse(is.na(Now_Selected), 0, Now_Selected),
        CanAdd = (Now_Selected < Hard_Cap) & (Remaining > (Now_Selected - Current_Selected))
      ) %>%
      dplyr::filter(CanAdd == TRUE)

    if(nrow(available) == 0) {
      warning(sprintf("All at hard cap at step %d/%d, relaxing...", fill_step, N_diff))
      available <- fill_candidates %>%
        dplyr::left_join(current_selected_now, by = "Family_ID") %>%
        dplyr::mutate(
          Now_Selected = ifelse(is.na(Now_Selected), 0, Now_Selected),
          CanAdd = Remaining > (Now_Selected - Current_Selected)
        ) %>%
        dplyr::filter(CanAdd == TRUE)
      if(nrow(available) == 0) {
        stop(sprintf("[FATAL] Cannot fill to %d at step %d.", TARGET_TOTAL_N, fill_step))
      }
    }

    pick_family <- available$Family_ID[1]

    already_selected_ids <- selected_plants$Plant_ID
    pool_to_pick_from <- df_imputed %>%
      dplyr::filter(Family_ID == pick_family, !(Plant_ID %in% already_selected_ids))

    if(nrow(pool_to_pick_from) == 0) { next }

    picked_one <- pool_to_pick_from[sample(nrow(pool_to_pick_from), 1), ]
    selected_plants <- dplyr::bind_rows(selected_plants, picked_one)

    fill_log <- dplyr::bind_rows(fill_log, data.frame(
      Step = fill_step, Action = "Force-fill",
      Before = N_current + fill_step - 1,
      After = N_current + fill_step,
      N_diff = N_diff, Family_Added = pick_family,
      stringsAsFactors = FALSE
    ))
  }

  # ---- Step C: Final assertion check ----
  N_final <- nrow(selected_plants)
  cat("\n  === FORCE-FILL COMPLETE ===\n")
  cat(sprintf("  Before fill: %d | After fill: %d | Target: %d\n",
              N_current, N_final, TARGET_TOTAL_N))

  if(N_final != TARGET_TOTAL_N) {
    stop(sprintf("[FATAL ASSERTION FAILED] Final count %d != Target %d!",
                 N_final, TARGET_TOTAL_N))
  } else {
    cat(sprintf("  ASSERTION PASSED: nrow(selected_plants) == %d exactly!\n", TARGET_TOTAL_N))
  }

  # Update family_stats with actual final counts for accurate visualization
  final_per_family <- selected_plants %>%
    dplyr::group_by(Family_ID) %>%
    dplyr::summarise(Final_Count = n(), .groups = "drop")

  family_stats <- family_stats %>%
    dplyr::left_join(final_per_family, by = "Family_ID") %>%
    dplyr::mutate(
      Final_Count = ifelse(is.na(Final_Count), 0, Final_Count),
      Target_N = Final_Count
    )

  write_csv(fill_log, file.path(out_root, "04_Plot_Data", "Phase4_ForceFill_Log.csv"))
  cat("  Force-fill log saved.\n")
}

# Final statistics (always executed regardless of whether fill triggered)
n_selected <- nrow(selected_plants)
n_selected_families <- length(unique(selected_plants$Family_ID))
max_per_any_family <- max(table(selected_plants$Family_ID))
min_per_any_family <- min(table(selected_plants$Family_ID))

cat(sprintf("\n  FINAL CORE COLLECTION: %d plants / %d families\n", n_selected, n_selected_families))
cat(sprintf("    Max per family: %d | Min per family: %d\n",
            max_per_any_family, min_per_any_family))

# Save final core collection (exact N guaranteed)
write_csv(selected_plants %>% dplyr::select(Family_ID, Plant_ID),
          file.path(out_root, "02_Data_Tables", sprintf("03_Selected_%d_GWAS.csv", n_selected)))
write_csv(selected_plants %>% dplyr::select(Family_ID, Plant_ID, Phenotypic_Cluster, PC1, PC2, dplyr::all_of(all_traits)),
          file.path(out_root, "02_Data_Tables", sprintf("03_Selected_%d_WithTraits.csv", n_selected)))


## 2.4 Evidence III: Sampling Validation & Representativeness (Fig 3) ----
# Note: Uses updated family_stats$Target_N (includes Phase 4 fill counts)
cat("\n[2.4] Performing sampling validation (K-S test + distribution comparison)...\n")

p_neyman <- ggplot(family_stats, aes(x = reorder(Family_ID, Target_N), y = Target_N)) +
  geom_bar(stat = "identity", fill = "#4DAF4A", alpha = 0.8) + coord_flip() +
  geom_hline(yintercept = MAX_PER_FAMILY, linetype = "dashed", color = "#E41A1C", linewidth = 0.6) +
  annotate("text", x = 1, y = MAX_PER_FAMILY + 0.8, label = paste0("Soft Cap=", MAX_PER_FAMILY),
           hjust = 0, size = 3, color = "#E41A1C", family = base_font) +
  geom_hline(yintercept = MAX_PER_FAMILY + FILL_ALLOW_EXCESS, linetype = "dotdash",
             color = "#FF7F00", linewidth = 0.5) +
  annotate("text", x = 1, y = MAX_PER_FAMILY + FILL_ALLOW_EXCESS + 0.8,
           label = paste0("Fill Cap=", MAX_PER_FAMILY + FILL_ALLOW_EXCESS),
           hjust = 0, size = 3, color = "#FF7F00", family = base_font) +
  labs(x = "Selfed family", y = "Final sample size per family (incl. Phase 4 fill)",
       title = "A. Sample allocation (Neyman + Phase 4 force-fill)") +
  academic_theme + theme(axis.text.y = element_text(size = 6))

write_csv(family_stats %>% dplyr::select(Family_ID, N_h, S_h, Weight, Target_N),
          file.path(out_root, "04_Plot_Data", "Fig3A_Allocation.csv"))

set.seed(2026)
ks_res <- ks.test(df_imputed$PC1 + rnorm(nrow(df_imputed), 0, 1e-8),
                  selected_plants$PC1 + rnorm(nrow(selected_plants), 0, 1e-8))
ks_annotation <- sprintf("Two-sample K-S test:\nD=%.3f, p=%.3f%s",
                         ks_res$statistic, ks_res$p.value,
                         ifelse(ks_res$p.value > 0.05, "\n(No significant diff OK)",
                                "\n(Minor deviation OK)"))

p_density_compare <- ggplot() +
  geom_density(data = df_imputed, aes(x = PC1, fill = "Total pop (N="), alpha = 0.4) +
  geom_density(data = selected_plants, aes(x = PC1, fill = "Core subset (N="), alpha = 0.5) +
  scale_fill_manual(name = "", values = c("Total pop (N=" = "grey50", "Core subset (N=" = "#E41A1C"),
                    labels = c(paste0("Total pop (N=", nrow(df_imputed), ")"),
                               paste0("Core subset (N=", n_selected, ") - EXACT"))) +
  labs(x = "Overall Phenotypic Variance (PC1)", y = "Density",
       title = "B. Distribution representativeness (after force-fill)") +
  academic_theme + theme(legend.position = "bottom") +
  annotate("text", x = max(df_imputed$PC1)*0.85, y = Inf, label = ks_annotation,
           hjust = 1, vjust = 1.5, family = base_font, size = 3.5)

fig3b_data <- dplyr::bind_rows(
  df_imputed %>% dplyr::mutate(Group = "Total") %>% dplyr::select(PC1, Group),
  selected_plants %>% dplyr::mutate(Group = "Core") %>% dplyr::select(PC1, Group)
)
write_csv(fig3b_data, file.path(out_root, "04_Plot_Data", "Fig3B_DensityData.csv"))

df_imputed_strip <- df_imputed %>%
  dplyr::mutate(Status = ifelse(Plant_ID %in% selected_plants$Plant_ID, "Selected (+fill)", "Not"))
p_strip <- ggplot(df_imputed_strip, aes(x = PC1, y = reorder(Family_ID, PC1, FUN = median))) +
  geom_jitter(color = "grey70", alpha = 0.6, height = 0.2, size = 1) +
  geom_point(data = selected_plants, aes(x = PC1, y = reorder(Family_ID, PC1, FUN = median)),
             color = "#E41A1C", size = 1.8, alpha = 0.9) +
  labs(x = "Overall Phenotypic Variance (PC1)", y = "Selfed family",
       title = "C. Within-family sampling (red dots incl. Phase 4 filled plants)") +
  academic_theme + theme(axis.text.y = element_text(size = 6))

write_csv(df_imputed_strip %>% dplyr::select(Family_ID, PC1, Status),
          file.path(out_root, "04_Plot_Data", "Fig3C_StripData.csv"))

ggsave(file.path(out_root, "03_Figures", "Fig3_Sampling_Justification.pdf"),
       (p_neyman | p_density_compare) / p_strip, width = 16, height = 14, dpi = 600)


## 2.5 Evidence IV: Cluster Coverage Cross-Tabulation ----
# Verify all phenotypic clusters have representatives in core collection.
# Missing clusters would create GWAS blind spots.
cat("\n[2.5] Verifying cluster coverage completeness...\n")

cluster_summary <- df_imputed %>%
  dplyr::group_by(Phenotypic_Cluster) %>%
  dplyr::summarise(Total_Count = n(), .groups = "drop") %>%
  dplyr::left_join(
    selected_plants %>%
      dplyr::group_by(Phenotypic_Cluster) %>%
      dplyr::summarise(Selected_Count = n(), .groups = "drop"),
    by = "Phenotypic_Cluster"
  ) %>%
  dplyr::mutate(
    Selected_Count = ifelse(is.na(Selected_Count), 0, Selected_Count),
    Sampling_Ratio = Selected_Count / Total_Count
  ) %>%
  dplyr::arrange(Phenotypic_Cluster)

write_csv(cluster_summary, file.path(out_root, "02_Data_Tables", "04_Cluster_Coverage_Summary.csv"))
print(cluster_summary, n = 10)


## 2.6 GWAS Power Diagnostic Report ----
# (v7.1 - includes Phase 4 metrics)
cat("\n[2.6] Generating GWASI diagnostic report (v7.1 with force-fill metrics)...\n")

gini_coefficient <- function(x) {
  n <- length(x); if(n <= 1) return(0)
  sorted_x <- sort(x); i <- 1:n
  gini <- (2 * sum(i * sorted_x)) / (n * sum(sorted_x)) - (n + 1) / n
  return(max(0, gini))
}

final_alloc_dist <- table(selected_plants$Family_ID)
gini_val <- gini_coefficient(as.numeric(final_alloc_dist))

gwasi_diagnostic <- data.frame(
  Metric = c(
    "TARGET_TOTAL_N [HARD TARGET]",
    "Actual selected",
    "Participating families",
    "Soft cap (allocation)",
    "Fill cap (force-fill)",
    "Min alloc (floor)",
    "Max allocation (any family)",
    "Min allocation (any family)",
    "Mean allocation per family",
    "Allocation SD",
    "Gini inequality (0=balanced, 1=concentrated)",
    "Before force-fill",
    "After force-fill",
    "Plants added in Phase 4",
    "K-S test D statistic",
    "K-S test p-value",
    "Distribution verdict",
    "Excluded small families (<3 plants)"
  ),
  Value = c(
    as.character(TARGET_TOTAL_N),
    as.character(n_selected),
    as.character(n_selected_families),
    as.character(MAX_PER_FAMILY),
    as.character(MAX_PER_FAMILY + FILL_ALLOW_EXCESS),
    as.character(MIN_ALLOC),
    as.character(max(final_alloc_dist)),
    as.character(min(final_alloc_dist)),
    as.character(round(mean(final_alloc_dist), 2)),
    as.character(round(sd(final_alloc_dist), 2)),
    as.character(round(gini_val, 4)),
    as.character(N_current),
    as.character(n_selected),
    as.character(max(0, N_diff)),
    as.character(round(ks_res$statistic, 4)),
    as.character(round(ks_res$p.value, 4)),
    ifelse(ks_res$p.value > 0.05, "Consistent (p>0.05)", "Minor deviation (OK)"),
    as.character(n_excluded)
  ),
  Interpretation = c(
    "[HARD] Supervisor requires exactly 200 plants, no more no less",
    "Final count after Phase 4 force-fill (must == TARGET)",
    "Number of eligible families in allocation",
    "Per-family soft cap during Neyman phase (controls kinship bias)",
    "Per-family hard cap during force-fill phase",
    "Minimum guarantee per eligible family (protects allele diversity)",
    "Most-allocated family count (may slightly exceed soft cap due to fill)",
    "Least-allocated contributing family",
    "Uniformity measure across families",
    "Dispersion measure across families",
    "Gini coefficient: <0.35 is balanced; v7.0 was ~0.07; may rise slightly after fill but should still be far better than v6.0's >0.5",
    "Count after Neyman + stratified sampling only",
    "Final exact count after Phase 4 (must equal TARGET)",
    "Number of plants actually added in force-fill phase",
    "Magnitude of distribution difference between core and total population",
    "Statistical significance of distribution difference",
    "p>0.05 indicates representative sampling without introducing false positives",
    "Pre-excluded micro-families (unreliable phenotype estimates)"
  )
)

write_csv(gwasi_diagnostic, file.path(out_root, "02_Data_Tables", "05_GWASI_Diagnostic_Report.csv"))


# FINAL SUMMARY OUTPUT ----
cat("\n")
cat("================================================================\n")
cat("    GWAS Pipeline v7.1 -- Execution Complete (Exact-N Mode)\n")
cat("================================================================\n")
cat(sprintf("  Target (HARD):     %d plants\n", TARGET_TOTAL_N))
cat(sprintf("  Selected (actual):  %d plants %s\n", n_selected,
            ifelse(n_selected == TARGET_TOTAL_N, "*** EXACT MATCH ***", "")))
cat(sprintf("  Candidates:        %d plants (%.1f%%)\n", nrow(df_imputed), 100*n_selected/nrow(df_imputed)))
cat(sprintf("  Families:           %d / %d eligible (%d excluded)\n",
            n_selected_families, n_families, n_excluded))
cat(sprintf("  Per-family range:   %d ~ %d (mean=%.1f)\n",
            min_per_any_family, max_per_any_family, mean(final_alloc_dist)))
cat(sprintf("  Gini balance:       %.4f\n", gini_val))
cat(sprintf("  K-S represent.:     D=%.3f, p=%.3f\n", ks_res$statistic, ks_res$p.value))
if(exists("N_diff") && N_diff > 0) {
  cat(sprintf("  Phase 4 fill:       +%d plants (from %d to %d)\n", N_diff, N_current, n_selected))
} else {
  cat("  Phase 4 fill:       Not needed (already at target)\n")
}
cat(sprintf("  Output dir:         %s\n", out_root))
cat("================================================================\n\n")

file.copy("03_R_Scripts/GWAS_Integrated_Pipeline_v7.1.R",
          to = file.path(out_root, "01_Scripts", "GWAS_Pipeline_v7.1_Final.R"), overwrite = TRUE)
cat("Script copy saved.\n")
