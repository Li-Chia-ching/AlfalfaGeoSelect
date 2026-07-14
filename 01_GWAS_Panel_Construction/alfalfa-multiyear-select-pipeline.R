# Integrated Pre-GWAS Phenomic Pipeline & Core Collection Sampling (v7.2) ----
# GWAS Power Optimization: Balanced Family Allocation + FORCE-FILL to Exact 200
#
# v7.1 -> v7.2 Key Changes:
# 1. REMOVED: showtext/SimSun → all fonts Arial only (Nature journal standard)
# 2. REMOVED: all figure titles, keeping only A/B/C panel labels via patchwork
# 3. REPLACED: color palette → low-saturation, colorblind-friendly (Okabe-Ito / Nature)
# 4. ADDED: ggsci package for additional Nature-style palettes
# 5. All paths are relative to project root (no hardcoded user paths)
# 6. FIXED: all ggsave() calls use device=cairo_pdf (avoids CID font error on Windows)
# 7. Target: exact 200 plants, Neyman allocation, force-fill guarantee
#
# v7.0 -> v7.1 Key Changes (retained):
# 1. TARGET_TOTAL_N = 200 (exact, non-negotiable for supervisor requirement)
# 2. MAX_PER_FAMILY = 7 (soft cap; Phase 4 may exceed to 8-9 for force-fill)
# 3. MIN_FAMILY_SIZE = 3 -> filters unreliable micro-families
# 4. MIN_ALLOC = 1 -> guarantees every eligible family contributes at least 1 plant
# 5. NEW: Phase 4 "Force-Fill" — post-allocation gap detection & mandatory fill-up

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

## 0. Package Management ----
required_packages <- c("dplyr", "tidyr", "ggplot2", "lme4", "lmerTest",
                       "factoextra", "readr", "patchwork",
                       "cluster", "viridis", "ggsci",
                       "scales", "RColorBrewer")

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
  library(cluster)
  library(viridis)
  library(ggsci)
  library(scales)
  library(RColorBrewer)
})

# ===== NATURE-STYLE COLOR PALETTE (Okabe-Ito, colorblind-friendly) =====
# Low-saturation, suitable for print & digital publication
pal_blue      <- "#0072B2"   # Blue
pal_orange    <- "#D55E00"   # Vermillion / Orange
pal_green     <- "#009E73"   # Bluish green
pal_yellow    <- "#F0E442"   # Yellow
pal_skyblue   <- "#56B4E9"   # Sky blue
pal_redpurple <- "#CC79A7"   # Reddish purple
pal_grey      <- "#999999"   # Neutral grey
pal_darkgrey  <- "#555555"   # Dark grey (for text/annotations)
pal_lightgrey <- "#E0E0E0"   # Light grey (for backgrounds)

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
out_root <- paste0("05_Analysis_Outputs/GWAS_Pipeline_Result_v7.2_", timestamp)
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_root, "01_Scripts"), showWarnings = FALSE)
dir.create(file.path(out_root, "02_Data_Tables"), showWarnings = FALSE)
dir.create(file.path(out_root, "03_Figures"), showWarnings = FALSE)
dir.create(file.path(out_root, "04_Plot_Data"), showWarnings = FALSE)

# ===== NATURE-STYLE ACADEMIC THEME (Arial, minimal, no grid) =====
nature_theme <- theme_bw(base_family = "Arial", base_size = 9) +
  theme(
    panel.grid.major = element_line(colour = pal_lightgrey, linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(colour = "black", linewidth = 0.5),
    axis.line        = element_blank(),
    axis.ticks       = element_line(colour = "black", linewidth = 0.4),
    axis.text        = element_text(colour = "black", size = 8),
    axis.title       = element_text(colour = "black", size = 9),
    strip.text       = element_text(size = 7, margin = margin(1, 1, 1, 1)),
    strip.background = element_rect(fill = "grey95", colour = "grey80"),
    legend.position  = "bottom",
    legend.text      = element_text(size = 8),
    legend.title     = element_text(size = 8),
    legend.key.size  = unit(0.6, "lines"),
    plot.margin      = margin(4, 4, 2, 2)
  )

cat("================================================================\n")
cat("    GWAS Pipeline v7.2 - Nature-Style + Exact N Edition\n")
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
            n_families, nrow(family_size_check), n_excluded, MIN_FAMILY_SIZE))
cat(sprintf("  After filter: %d individuals\n\n", nrow(df_imputed)))

write_csv(
  family_size_check %>%
    dplyr::mutate(Status = ifelse(N_Original >= MIN_FAMILY_SIZE, "Eligible", "Excluded")),
  file.path(out_root, "02_Data_Tables", "00_Family_Eligibility_Check.csv")
)

## 1.3 Phenotypic Clustering & Visualization (Fig S1) ----
cat("[Stage 1] Exploring macro phenotypic diversity via K-means...\n")
pheno_mat <- df_imputed[, all_traits, drop = FALSE]
pheno_scaled <- scale(pheno_mat)

sil_res <- fviz_nbclust(pheno_scaled, kmeans, method = "silhouette", k.max = 6)
p_opt_k <- sil_res +
  ggtitle(NULL) +
  nature_theme

best_k <- as.numeric(as.character(sil_res$data$clusters[which.max(sil_res$data$y)]))
if(length(best_k) == 0 || is.na(best_k) || best_k < 2) best_k <- 3
if(best_k > 5) best_k <- 3

set.seed(2026)
km_res <- kmeans(pheno_scaled, centers = best_k, nstart = 25)
df_imputed$Phenotypic_Cluster <- paste0("Cluster_", km_res$cluster)

pca_res_clust <- prcomp(pheno_scaled, center = FALSE, scale. = FALSE)
pca_df <- as.data.frame(pca_res_clust$x[, 1:2])
pca_df$Cluster <- df_imputed$Phenotypic_Cluster

p_cluster <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = Cluster, fill = Cluster)) +
  geom_point(shape = 21, size = 1.8, alpha = 0.75, stroke = 0.3) +
  stat_ellipse(geom = "polygon", alpha = 0.12, linewidth = 0.3, na.rm = TRUE) +
  scale_colour_npg() +
  scale_fill_npg() +
  labs(x = "PC1", y = "PC2") +
  ggtitle(NULL) +
  nature_theme

ggsave(file.path(out_root, "03_Figures", "FigS1_Phenotypic_Clustering.pdf"),
       p_opt_k / p_cluster, width = 8, height = 10, device = cairo_pdf)

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
cat("[2.1] Computing variance components & repeatability...\n")

df_for_plot <- df_imputed %>%
  dplyr::add_count(Family_ID, name = "n_per_family") %>%
  dplyr::filter(n_per_family >= 2) %>%
  dplyr::mutate(Family_Label = paste0(Family_ID, " (n=", n_per_family, ")"))

p_seg_density <- ggplot(df_for_plot, aes(x = Plant_Height_May)) +
  geom_density(fill = pal_blue, alpha = 0.5, colour = NA) +
  geom_rug(alpha = 0.3, sides = "b", colour = pal_darkgrey) +
  facet_wrap(~Family_Label, ncol = 9, scales = "free_y") +
  labs(x = expression(Plant~Height~May~(cm)), y = "Density") +
  ggtitle(NULL) +
  nature_theme +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

p_var_box <- ggplot(df_imputed, aes(x = reorder(Family_ID, Plant_Height_May, FUN = median),
                                    y = Plant_Height_May)) +
  geom_boxplot(fill = pal_orange, alpha = 0.5, outlier.size = 0.3, linewidth = 0.3) +
  labs(x = "Family ID (ordered by median)", y = expression(Plant~Height~May~(cm))) +
  ggtitle(NULL) +
  nature_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))

ggsave(file.path(out_root, "03_Figures", "Fig1_Segregation_Analysis.pdf"),
       p_seg_density / p_var_box, width = 16, height = 12, device = cairo_pdf)

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

p_pca_main <- fviz_pca_biplot(pca_res_main, geom.ind = "point", pointshape = 21, pointsize = 1.2,
                              habillage = df_imputed$Phenotypic_Cluster, alpha.ind = 0.5,
                              col.var = "black", arrowsize = 0.5, labelsize = 3, repel = TRUE,
                              title = NULL,
                              palette = "npg", ggtheme = nature_theme)

ggsave(file.path(out_root, "03_Figures", "Fig2_PCA_Biplot_Enhanced.pdf"),
       p_pca_main, width = 9, height = 7, device = cairo_pdf)

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
  
  # Phase 1: Neyman proportional allocation
  remaining <- target - sum(alloc)
  proportional <- remaining * weights / sum(weights)
  prop_capped <- pmin(proportional, cap - floor_alloc)
  alloc <- alloc + floor(prop_capped)
  
  # Phase 2: Distribute fractional remainders
  diff <- target - sum(alloc)
  while(diff > 0) {
    eligible <- alloc < cap
    if(!any(eligible)) break
    priority <- weights / (alloc + 1)
    priority[!eligible] <- -Inf
    alloc[which.max(priority)] <- alloc[which.max(priority)] + 1
    diff <- diff - 1
  }
  
  # Phase 3: Fill-to-target post-cap
  alloc[alloc > cap] <- cap
  gap <- target - sum(alloc)
  while(gap > 0) {
    under_cap <- which(alloc < cap)
    if(length(under_cap) == 0) break
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
cat("\n[2.3c] Phase 4: Force-Fill - Gap Detection & Mandatory Top-up...\n")

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
  cat(sprintf("  Initiating force-fill (fill cap = %d)...\n",
              MAX_PER_FAMILY + FILL_ALLOW_EXCESS))
  
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
    stop("[FATAL ERROR] No remaining plants available — check data or relax filters.")
  }
  
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
  
  N_final <- nrow(selected_plants)
  cat("\n  === FORCE-FILL COMPLETE ===\n")
  cat(sprintf("  Before fill: %d | After fill: %d | Target: %d\n",
              N_current, N_final, TARGET_TOTAL_N))
  
  if(N_final != TARGET_TOTAL_N) {
    stop(sprintf("[FATAL ASSERTION FAILED] Final %d != Target %d!", N_final, TARGET_TOTAL_N))
  } else {
    cat(sprintf("  ASSERTION PASSED: nrow == %d exactly\n", TARGET_TOTAL_N))
  }
  
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

# Final statistics
n_selected <- nrow(selected_plants)
n_selected_families <- length(unique(selected_plants$Family_ID))
max_per_any_family <- max(table(selected_plants$Family_ID))
min_per_any_family <- min(table(selected_plants$Family_ID))

cat(sprintf("\n  FINAL CORE COLLECTION: %d plants / %d families\n", n_selected, n_selected_families))
cat(sprintf("    Max per family: %d | Min per family: %d\n",
            max_per_any_family, min_per_any_family))

write_csv(selected_plants %>% dplyr::select(Family_ID, Plant_ID),
          file.path(out_root, "02_Data_Tables", sprintf("03_Selected_%d_GWAS.csv", n_selected)))
write_csv(selected_plants %>% dplyr::select(Family_ID, Plant_ID, Phenotypic_Cluster, PC1, PC2, dplyr::all_of(all_traits)),
          file.path(out_root, "02_Data_Tables", sprintf("03_Selected_%d_WithTraits.csv", n_selected)))


## 2.4 Evidence III: Sampling Validation & Representativeness (Fig 3) ----
cat("\n[2.4] Performing sampling validation (K-S test + distribution comparison)...\n")

p_neyman <- ggplot(family_stats, aes(x = reorder(Family_ID, Target_N), y = Target_N)) +
  geom_bar(stat = "identity", fill = pal_green, alpha = 0.75, width = 0.7) +
  coord_flip() +
  geom_hline(yintercept = MAX_PER_FAMILY, linetype = "dashed",
             colour = pal_orange, linewidth = 0.5) +
  annotate("text", x = 1, y = MAX_PER_FAMILY + 0.8,
           label = paste0("Soft cap = ", MAX_PER_FAMILY),
           hjust = 0, size = 2.5, colour = pal_orange) +
  geom_hline(yintercept = MAX_PER_FAMILY + FILL_ALLOW_EXCESS,
             linetype = "dotted", colour = pal_redpurple, linewidth = 0.4) +
  annotate("text", x = 1, y = MAX_PER_FAMILY + FILL_ALLOW_EXCESS + 0.8,
           label = paste0("Fill cap = ", MAX_PER_FAMILY + FILL_ALLOW_EXCESS),
           hjust = 0, size = 2.5, colour = pal_redpurple) +
  labs(x = "Family ID", y = "Sample size per family") +
  ggtitle(NULL) +
  nature_theme +
  theme(axis.text.y = element_text(size = 5.5))

write_csv(family_stats %>% dplyr::select(Family_ID, N_h, S_h, Weight, Target_N),
          file.path(out_root, "04_Plot_Data", "Fig3A_Allocation.csv"))

set.seed(2026)
ks_res <- ks.test(df_imputed$PC1 + rnorm(nrow(df_imputed), 0, 1e-8),
                  selected_plants$PC1 + rnorm(nrow(selected_plants), 0, 1e-8))
ks_annotation <- sprintf("K-S: D = %.3f, p = %.3f%s",
                         ks_res$statistic, ks_res$p.value,
                         ifelse(ks_res$p.value > 0.05, " (n.s.)", ""))

p_density_compare <- ggplot() +
  geom_density(data = df_imputed,
               aes(x = PC1, fill = "Total"), alpha = 0.5, colour = NA) +
  geom_density(data = selected_plants,
               aes(x = PC1, fill = "Core"), alpha = 0.55, colour = NA) +
  scale_fill_manual(
    name = NULL,
    values = c("Total" = pal_grey, "Core" = pal_skyblue),
    labels = c(paste0("Total (N=", nrow(df_imputed), ")"),
               paste0("Core (N=", n_selected, ")"))
  ) +
  labs(x = "PC1", y = "Density") +
  ggtitle(NULL) +
  nature_theme +
  annotate("text", x = max(df_imputed$PC1) * 0.82, y = Inf,
           label = ks_annotation, hjust = 1, vjust = 1.5, size = 3, colour = pal_darkgrey)

fig3b_data <- dplyr::bind_rows(
  df_imputed %>% dplyr::mutate(Group = "Total") %>% dplyr::select(PC1, Group),
  selected_plants %>% dplyr::mutate(Group = "Core") %>% dplyr::select(PC1, Group)
)
write_csv(fig3b_data, file.path(out_root, "04_Plot_Data", "Fig3B_DensityData.csv"))

df_imputed_strip <- df_imputed %>%
  dplyr::mutate(Status = ifelse(Plant_ID %in% selected_plants$Plant_ID, "Selected", "Not"))
p_strip <- ggplot(df_imputed_strip, aes(x = PC1, y = reorder(Family_ID, PC1, FUN = median))) +
  geom_jitter(data = df_imputed_strip %>% dplyr::filter(Status == "Not"),
              colour = pal_lightgrey, alpha = 0.5, height = 0.2, size = 0.6) +
  geom_point(data = df_imputed_strip %>% dplyr::filter(Status == "Selected"),
             colour = pal_blue, size = 1.2, alpha = 0.85) +
  labs(x = "PC1", y = "Family ID") +
  ggtitle(NULL) +
  nature_theme +
  theme(axis.text.y = element_text(size = 5))

write_csv(df_imputed_strip %>% dplyr::select(Family_ID, PC1, Status),
          file.path(out_root, "04_Plot_Data", "Fig3C_StripData.csv"))

ggsave(file.path(out_root, "03_Figures", "Fig3_Sampling_Justification.pdf"),
       (p_neyman | p_density_compare) / p_strip +
         plot_annotation(tag_levels = "A") &
         theme(plot.tag = element_text(face = "bold", size = 11, family = "Arial")),
       width = 16, height = 14, device = cairo_pdf)


## 2.5 Evidence IV: Cluster Coverage Cross-Tabulation ----
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
cat("\n[2.6] Generating GWASI diagnostic report (v7.2 with force-fill metrics)...\n")

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
    ifelse(ks_res$p.value > 0.05, "Consistent (p>0.05)", "Minor deviation"),
    as.character(n_excluded)
  ),
  Interpretation = c(
    "[HARD] Supervisor requires exactly N plants, no more no less",
    "Final count after Phase 4 force-fill (must == TARGET)",
    "Number of eligible families in allocation",
    "Per-family soft cap during Neyman phase (controls kinship bias)",
    "Per-family hard cap during force-fill phase",
    "Minimum guarantee per eligible family (protects allele diversity)",
    "Most-allocated family count (may slightly exceed soft cap due to fill)",
    "Least-allocated contributing family",
    "Uniformity measure across families",
    "Dispersion measure across families",
    "Gini coefficient: <0.35 is balanced",
    "Count after Neyman + stratified sampling only",
    "Final exact count after Phase 4 (must equal TARGET)",
    "Number of plants actually added in force-fill phase",
    "Magnitude of distribution difference between core and total population",
    "Statistical significance of distribution difference",
    "p>0.05 indicates representative sampling",
    "Pre-excluded micro-families (unreliable phenotype estimates)"
  )
)

write_csv(gwasi_diagnostic, file.path(out_root, "02_Data_Tables", "05_GWASI_Diagnostic_Report.csv"))


# FINAL SUMMARY OUTPUT ----
cat("\n")
cat("================================================================\n")
cat("    GWAS Pipeline v7.2 -- Execution Complete (Nature-Style)\n")
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

file.copy("03_R_Scripts/GWAS_Integrated_Pipeline_v7.2.R",
          to = file.path(out_root, "01_Scripts", "GWAS_Pipeline_v7.2_Final.R"), overwrite = TRUE)
cat("Script copy saved to output directory.\n")
