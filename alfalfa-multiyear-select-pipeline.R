# =============================================================================
# Integrated Pre-GWAS Phenomic Pipeline & Core Collection Sampling (v6.0)
# Stage 1: Data Imputation & Phenotypic Clustering (Macro-structure)
# Stage 2: Weighted Neyman Allocation & Stratified Sampling (Micro-structure)
# Note: Spatial correction removed due to over-fitting (all within-family SD became 0)
# =============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

# -------------------------------------------------------------------------
# 0. Initialization & Package Management
# -------------------------------------------------------------------------
required_packages <- c("dplyr", "tidyr", "ggplot2", "lme4", "lmerTest", 
                       "factoextra", "readr", "patchwork", "showtext", 
                       "cluster", "viridis")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(lme4)
  library(lmerTest); library(factoextra); library(readr)
  library(patchwork); library(showtext); library(cluster); library(viridis)
})

tryCatch({
  if (.Platform$OS.type == "windows") { font_add("SimSun", "simsun.ttc") } 
  else { font_add("SimSun", "STSong.ttf") }
  showtext_auto()
  base_font <- "SimSun"
}, error = function(e) {
  base_font <- ""
})

out_dir <- paste0("GWAS_Integrated_Pipeline_", format(Sys.Date(), "%Y%m%d"))
dir.create(out_dir, showWarnings = FALSE)

academic_theme <- theme_bw(base_family = base_font) + 
  theme(panel.grid.minor = element_blank(),
        text = element_text(family = base_font, size = 12),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 14),
        strip.text = element_text(size = 9, margin = margin(2, 2, 2, 2)))

TARGET_TOTAL_N <- 140  
MAX_PER_FAMILY <- 20   

all_traits <- c("Plant_Height_Nov", "Internode_Nov", "Plant_Height_Mar", 
                "Internode_Mar", "Branch_Number_Mar", "Multifoliate_Score_Mar", 
                "Plant_Height_May")

# =============================================================================
# STAGE 1: DATA LOADING, IMPUTATION & PHENOTYPIC CLUSTERING
# =============================================================================

# -------------------------------------------------------------------------
# 1.1 Data Loading & Robust Imputation
# -------------------------------------------------------------------------
cat("[Stage 1] Loading data and performing robust imputation...\n")
df_raw <- read_csv("data_202605.csv", show_col_types = FALSE) %>%
  dplyr::rename(Family_ID = Family, Plant_ID = ID) %>%
  dplyr::filter(!is.na(Plant_Height_May))

df_imputed <- df_raw %>%
  dplyr::group_by(Family_ID) %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(all_traits), ~ifelse(is.na(.), median(., na.rm = TRUE), .))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(all_traits), ~ifelse(is.na(.), median(., na.rm = TRUE), .)))

write_csv(df_imputed %>% dplyr::select(Family_ID, Plant_ID, dplyr::all_of(all_traits)), 
          file.path(out_dir, "00_Imputed_Dataset.csv"))

# -------------------------------------------------------------------------
# 1.2 Phenotypic Clustering & Visualization (Fig S1)
# -------------------------------------------------------------------------
cat("[Stage 1] Exploring macro phenotypic diversity via K-means...\n")
pheno_mat <- df_imputed[, all_traits, drop = FALSE]
pheno_scaled <- scale(pheno_mat)

sil_res <- fviz_nbclust(pheno_scaled, kmeans, method = "silhouette", k.max = 5)
p_opt_k <- sil_res + academic_theme + labs(title = "A. Optimal number of phenotypic clusters")

best_k <- as.numeric(as.character(sil_res$data$clusters[which.max(sil_res$data$y)]))
if (length(best_k) == 0 || is.na(best_k) || best_k < 2) {
  best_k <- 3
} else if (best_k > 4) {
  best_k <- 3
}

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
  labs(title = sprintf("B. Multidimensional Phenotypic Structure (K = %d)", best_k),
       subtitle = "Imputed original traits (no spatial correction)", x = "PC1", y = "PC2") +
  academic_theme

ggsave(file.path(out_dir, "FigS1_Phenotypic_Clustering.pdf"), 
       p_opt_k / p_cluster, width = 8, height = 10, dpi = 600)

write_csv(sil_res$data, file.path(out_dir, "FigS1A_Silhouette_Data.csv"))

plot_data_pca <- pca_df %>%
  dplyr::mutate(Plant_ID = df_imputed$Plant_ID, Family_ID = df_imputed$Family_ID) %>%
  dplyr::bind_cols(df_imputed %>% dplyr::select(dplyr::all_of(all_traits))) %>%
  dplyr::select(Family_ID, Plant_ID, Cluster, PC1, PC2, dplyr::everything())
write_csv(plot_data_pca, file.path(out_dir, "FigS1B_PCA_Cluster_Data.csv"))

# =============================================================================
# STAGE 2: GWAS SAMPLING PIPELINE (Using Imputed Original Data)
# =============================================================================
cat("\n[Stage 2] Commencing downstream validation and Neyman sampling...\n")

# -------------------------------------------------------------------------
# 2.1 Evidence I: Segregation & Repeatability (Fig 1)
# -------------------------------------------------------------------------
df_for_plot <- df_imputed %>%
  dplyr::add_count(Family_ID, name = "n_per_family") %>%
  dplyr::filter(n_per_family >= 2) %>%
  dplyr::mutate(Family_Label = paste0(Family_ID, "\n(n=", n_per_family, ")"))

p_seg_density <- ggplot(df_for_plot, aes(x = Plant_Height_May)) +
  geom_density(fill = "#377EB8", alpha = 0.6) +
  geom_rug(alpha = 0.5, sides = "b") +
  facet_wrap(~Family_Label, ncol = 9, scales = "free_y") +
  labs(x = "Plant Height in May 2026 (cm)", y = "Kernel density", 
       title = "A. Density distribution of plant height segregation") +
  academic_theme + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

p_var_box <- ggplot(df_imputed, aes(x = reorder(Family_ID, Plant_Height_May, FUN = median), 
                                    y = Plant_Height_May)) +
  geom_boxplot(fill = "#E41A1C", alpha = 0.6, outlier.size = 0.5) +
  labs(x = "Selfed family ID (ordered by median)", y = "Plant Height in May 2026 (cm)", 
       title = "B. Between-family variation and within-family segregation range") +
  academic_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

ggsave(file.path(out_dir, "Fig1_Segregation_Analysis.pdf"), 
       p_seg_density / p_var_box, width = 16, height = 12, dpi = 600)

write_csv(df_imputed %>% dplyr::select(Family_ID, Plant_ID, Plant_Height_May), 
          file.path(out_dir, "Fig1_Data.csv"))

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
  
  data.frame(Trait = trait, 
             Var_Family = var_fam_raw, 
             Var_Within = var_wit_raw, 
             Repeatability_R = R_val, 
             Note = "Converged")
}
varcomp_results <- do.call(rbind, lapply(all_traits, get_varcomp)) %>% dplyr::arrange(desc(Repeatability_R))
write_csv(varcomp_results, file.path(out_dir, "01_Variance_Components_Repeatability.csv"))

# -------------------------------------------------------------------------
# 2.2 Evidence II: PCA Gradients & PCA Visualization (Fig 2)
# -------------------------------------------------------------------------
pca_res_main <- prcomp(df_imputed %>% dplyr::select(dplyr::all_of(all_traits)), center = TRUE, scale. = TRUE)
df_imputed <- df_imputed %>% dplyr::mutate(PC1 = pca_res_main$x[,1], PC2 = pca_res_main$x[,2])
prop_var <- summary(pca_res_main)$importance[2, 1:2] 

pca_summary <- summary(pca_res_main)
eig_data <- data.frame(
  PC = 1:length(pca_summary$importance[1,]),
  StdDev = pca_summary$importance[1,],
  Variance_Proportion = pca_summary$importance[2,],
  Cumulative_Proportion = pca_summary$importance[3,]
)
write_csv(eig_data, file.path(out_dir, "02_PCA_Eigenvalues.csv"))

loadings_data <- as.data.frame(pca_res_main$rotation[, 1:2]) %>%
  dplyr::rename(PC1_loading = PC1, PC2_loading = PC2) %>%
  tibble::rownames_to_column(var = "Trait")
write_csv(loadings_data, file.path(out_dir, "02_PCA_Loadings.csv"))

write_csv(df_imputed %>% dplyr::select(Family_ID, Plant_ID, PC1, PC2, dplyr::all_of(all_traits)), 
          file.path(out_dir, "02_PCA_Scores.csv"))

p_pca_main <- fviz_pca_biplot(pca_res_main, geom.ind = "point", pointshape = 21, pointsize = 1.5, 
                              habillage = df_imputed$Phenotypic_Cluster, alpha.ind = 0.6,
                              col.var = "black", arrowsize = 0.6, labelsize = 3.5, repel = TRUE, 
                              title = "Multidimensional Phenotypic Gradients (PCA Biplot)",
                              palette = "turbo", ggtheme = academic_theme)

ggsave(file.path(out_dir, "Fig2_PCA_Biplot_Enhanced.pdf"), p_pca_main, width = 9, height = 7, dpi = 600)

# -------------------------------------------------------------------------
# 2.3 Weighted Neyman Allocation & Stratified Gradient Sampling
# -------------------------------------------------------------------------
family_stats <- df_imputed %>%
  dplyr::group_by(Family_ID) %>%
  dplyr::summarise(N_h = n(), 
                   Var_Weighted = var(PC1)*prop_var[1] + var(PC2)*prop_var[2], .groups = "drop") %>%
  dplyr::mutate(Var_Weighted = ifelse(is.na(Var_Weighted) | Var_Weighted < 0, 0, Var_Weighted),
                S_h = sqrt(Var_Weighted))

if(sum(family_stats$S_h) == 0) family_stats$S_h <- 1 
family_stats$Weight <- family_stats$N_h * family_stats$S_h

smart_allocate_with_cap <- function(weights, target, cap) {
  alloc <- floor(target * weights / sum(weights))
  alloc[alloc > cap] <- cap
  diff <- target - sum(alloc)
  while(diff > 0) {
    eligible <- alloc < cap
    if(!any(eligible)) break
    priority <- weights / (alloc + 1)
    priority[!eligible] <- -Inf
    alloc[which.max(priority)] <- alloc[which.max(priority)] + 1
    diff <- diff - 1
  }
  return(alloc)
}
family_stats$Target_N <- smart_allocate_with_cap(family_stats$Weight, TARGET_TOTAL_N, MAX_PER_FAMILY)

set.seed(2026) 
selected_plants <- data.frame()

for (fam in family_stats$Family_ID) {
  fam_data <- df_imputed %>% dplyr::filter(Family_ID == fam)
  n_alloc <- family_stats$Target_N[family_stats$Family_ID == fam]
  
  if (n_alloc == 0) next
  if (n_alloc >= nrow(fam_data)) { selected_plants <- dplyr::bind_rows(selected_plants, fam_data); next }
  
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
  for (i in seq_along(grad_counts$Var1)) {
    grad_pool <- fam_data %>% dplyr::filter(Gradient == grad_counts$Var1[i])
    n_to_pick <- min(grad_alloc[i], nrow(grad_pool))
    if(n_to_pick > 0) fam_sampled <- dplyr::bind_rows(fam_sampled, grad_pool[sample(nrow(grad_pool), n_to_pick), ])
  }
  selected_plants <- dplyr::bind_rows(selected_plants, fam_sampled)
}

write_csv(selected_plants %>% dplyr::select(Family_ID, Plant_ID), 
          file.path(out_dir, sprintf("03_Selected_%d_GWAS.csv", TARGET_TOTAL_N)))

write_csv(selected_plants %>% dplyr::select(Family_ID, Plant_ID, Phenotypic_Cluster, PC1, PC2, dplyr::all_of(all_traits)), 
          file.path(out_dir, sprintf("03_Selected_%d_WithTraits.csv", TARGET_TOTAL_N)))

# -------------------------------------------------------------------------
# 2.4 Evidence III: Sampling Validation & Cluster Coverage (Fig 3)
# -------------------------------------------------------------------------
p_neyman <- ggplot(family_stats, aes(x = reorder(Family_ID, Target_N), y = Target_N)) +
  geom_bar(stat = "identity", fill = "#4DAF4A", alpha = 0.8) + coord_flip() +
  labs(x = "Selfed family", y = "Allocated sample size (Neyman)", 
       title = "A. Sample allocation based on variance and size") +
  academic_theme + theme(axis.text.y = element_text(size = 6))

write_csv(family_stats %>% dplyr::select(Family_ID, N_h, S_h, Weight, Target_N), 
          file.path(out_dir, "Fig3A_Allocation.csv"))

set.seed(2026)
ks_res <- ks.test(df_imputed$PC1 + rnorm(nrow(df_imputed), 0, 1e-8), 
                  selected_plants$PC1 + rnorm(nrow(selected_plants), 0, 1e-8))
ks_annotation <- sprintf("Two-sample exact K-S test:\nD = %.3f, p = %.3f", ks_res$statistic, ks_res$p.value)

p_density_compare <- ggplot() +
  geom_density(data = df_imputed, aes(x = PC1, fill = "Total population"), alpha = 0.4) +
  geom_density(data = selected_plants, aes(x = PC1, fill = "Core subset"), alpha = 0.5) +
  scale_fill_manual(name = "", values = c("Total population" = "grey50", "Core subset" = "#E41A1C")) +
  labs(x = "Overall Phenotypic Variance (PC1)", y = "Density", title = "B. Distribution representativeness") +
  academic_theme + theme(legend.position = "bottom") +
  annotate("text", x = max(df_imputed$PC1)*0.9, y = Inf, label = ks_annotation, 
           hjust = 1, vjust = 1.5, family = base_font, size = 3.5)

fig3b_data <- dplyr::bind_rows(
  df_imputed %>% dplyr::mutate(Group = "Total") %>% dplyr::select(PC1, Group),
  selected_plants %>% dplyr::mutate(Group = "Core") %>% dplyr::select(PC1, Group)
)
write_csv(fig3b_data, file.path(out_dir, "Fig3B_DensityData.csv"))

df_imputed_strip <- df_imputed %>% dplyr::mutate(Status = ifelse(Plant_ID %in% selected_plants$Plant_ID, "Selected", "Not"))
p_strip <- ggplot(df_imputed_strip, aes(x = PC1, y = reorder(Family_ID, PC1, FUN = median))) +
  geom_jitter(color = "grey70", alpha = 0.6, height = 0.2, size = 1) +
  geom_point(data = selected_plants, aes(x = PC1, y = reorder(Family_ID, PC1, FUN = median)), 
             color = "#E41A1C", size = 1.8, alpha = 0.9) +
  labs(x = "Overall Phenotypic Variance (PC1)", y = "Selfed family", 
       title = "C. Within-family quantile sampling avoids extreme bias") +
  academic_theme + theme(axis.text.y = element_text(size = 6))

write_csv(df_imputed_strip %>% dplyr::select(Family_ID, PC1, Status), 
          file.path(out_dir, "Fig3C_StripData.csv"))

ggsave(file.path(out_dir, "Fig3_Sampling_Justification.pdf"), 
       (p_neyman | p_density_compare) / p_strip, width = 16, height = 14, dpi = 600)

# -------------------------------------------------------------------------
# 2.5 Evidence IV: Cluster Coverage Cross-Tabulation
# -------------------------------------------------------------------------
cluster_summary <- df_imputed %>%
  dplyr::group_by(Phenotypic_Cluster) %>%
  dplyr::summarise(Total_Count = n(), .groups = "drop") %>%
  dplyr::left_join(
    selected_plants %>% 
      dplyr::group_by(Phenotypic_Cluster) %>% 
      dplyr::summarise(Selected_Count = n(), .groups = "drop"),
    by = "Phenotypic_Cluster"
  ) %>%
  dplyr::mutate(Selected_Count = ifelse(is.na(Selected_Count), 0, Selected_Count),
                Sampling_Ratio = Selected_Count / Total_Count) %>%
  dplyr::arrange(Phenotypic_Cluster)

write_csv(cluster_summary, file.path(out_dir, "04_Cluster_Coverage_Summary.csv"))

cat("\n[Pipeline Complete] Integrated workflow finished successfully.\n")
cat("Output Directory: ", out_dir, "\n")
