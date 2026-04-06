# ==============================================================================
# Genetic Variation Analysis + Multivariate Gradient-Based Stratified Sampling for GWAS
# Features:
#   - Multivariate variance-based Neyman allocation
#   - Quantile-gradient deterministic sampling
#   - LMM-based repeatability estimation
#   - Publication-quality visualization and data export
# Application:
#   - Two-year common environment phenotyping dataset
#   - Extraction of a 200-individual core subset for GWAS
# ==============================================================================


# ==============================================================================
# 1. Dependency Check and Package Loading----
# ==============================================================================
required_packages <- c("dplyr", "tidyr", "ggplot2", "lme4", "lmerTest", 
                       "factoextra", "readr", "patchwork", "showtext")

new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(factoextra)
library(readr)
library(patchwork) 
library(showtext) 


# ==============================================================================
# 2. Global Font Configuration and Academic Plot Theme Initialization----
# ==============================================================================
showtext_auto() 

if (.Platform$OS.type == "windows") {
  font_add("SimSun", "simsun.ttc") 
} else {
  font_add("SimSun", "STSong.ttf") 
}

out_dir <- paste0("GWAS_Sampling_Results_", format(Sys.Date(), "%Y%m%d"))
dir.create(out_dir, showWarnings = FALSE)

academic_theme <- theme_bw(base_family = "SimSun") + 
  theme(
    panel.grid.minor = element_blank(),
    text = element_text(family = "SimSun", size = 12),
    axis.title = element_text(family = "SimSun", face = "bold"),
    plot.title = element_text(family = "SimSun", face = "bold", size = 14),
    strip.text = element_text(size = 9, margin = margin(2, 2, 2, 2)),
    panel.spacing = unit(1, "lines")
  )


# ==============================================================================
# 3. Data Processing (Data Lineage — Merge Integrity Safeguard)----
# ==============================================================================
cat("Reading and preprocessing two-year phenotypic data...\n")

df_2025 <- read_csv("data_2025.csv", show_col_types = FALSE) %>%
  rename(
    Family_ID = Parent,
    Plant_ID = SampleID,
    Height_2025 = PlantHeight,
    Internode_2025 = NodeCount
  ) %>%
  select(Family_ID, Plant_ID, Real_Count, Height_2025, Internode_2025)

df_2026 <- read_csv("data_2026.csv", show_col_types = FALSE) %>%
  filter(!is.na(Group) & trimws(Group) != "") %>%
  rename(
    Family_ID = Group,
    Plant_ID = Matrix,
    Height_2026 = Plant_Height,
    Internode_2026 = Internode,
    Branch_2026 = Branch_Number,
    Multi_2026 = Multifoliate_Score
  ) %>%
  select(Family_ID, Plant_ID, Height_2026, Internode_2026, Branch_2026, Multi_2026, Death_Code)

df_wide <- inner_join(df_2025, df_2026, by = c("Family_ID", "Plant_ID")) %>%
  mutate(
    Status = case_when(
      Death_Code == 1 ~ "Dead_in_2025",
      Death_Code == 2 ~ "Dead_in_2026",
      TRUE ~ "Alive"
    ),
    Delta_Height = Height_2026 - Height_2025
  )

write_csv(df_wide, file.path(out_dir, "01_Merged_Wide_Data_All.csv"))

df_alive <- df_wide %>%
  filter(Status == "Alive") %>%
  drop_na(Height_2026, Internode_2026, Branch_2026, Multi_2026)


# ==============================================================================
# 4. Phenotypic Segregation Evidence (LMM Reformulated as Repeatability)----
# ==============================================================================
cat("Generating phenotypic segregation patterns and variance structure...\n")

## ---- Fig.1A Data: Raw plant height values (no sample size annotation) ----
write_csv(
  df_alive %>% select(Family_ID, Height_2026),
  file.path(out_dir, "Fig1A_Data.csv")
)

### ---- Add per-family sample size to facet labels (fix density issue for n=1) ----
df_alive_for_plot <- df_alive %>%
  add_count(Family_ID, name = "n_per_family") %>%
  mutate(Family_Label = paste0(Family_ID, " (n = ", n_per_family, ")"))

### ---- Density plot + rug (explicit display of individual observations) ----
p_seg_density <- ggplot(df_alive_for_plot, aes(x = Height_2026)) +
  geom_density(fill = "#377EB8", alpha = 0.6) +
  geom_rug(alpha = 0.5, sides = "b") +
  facet_wrap(~Family_Label, ncol = 8, scales = "free_y") +
  labs(
    x = "Plant height at early flowering stage in 2026 (cm)",
    y = "Kernel density estimate",
    title = "A. Within-family phenotypic segregation density distribution"
  ) +
  academic_theme +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

## ---- Fig.1B Data ----
write_csv(
  df_alive %>% select(Family_ID, Height_2026),
  file.path(out_dir, "Fig1B_Data.csv")
)

p_var_box <- ggplot(df_alive, aes(
  x = reorder(Family_ID, Height_2026, FUN = median),
  y = Height_2026
)) +
  geom_boxplot(fill = "#E41A1C", alpha = 0.6, outlier.size = 0.3) +
  labs(
    x = "Inbred family (ordered by median)",
    y = "Plant height at early flowering stage in 2026 (cm)",
    title = "B. Between-family divergence vs within-family segregation amplitude"
  ) +
  academic_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))

ggsave(
  file.path(out_dir, "Fig1_Segregation_Analysis.pdf"),
  p_seg_density / p_var_box,
  width = 16, height = 12, dpi = 600
)

## ---- Variance components and repeatability estimation ----
traits <- c("Height_2026", "Internode_2026", "Branch_2026", "Multi_2026")

get_varcomp <- function(trait) {
  form <- as.formula(paste(trait, "~ 1 + (1|Family_ID)"))
  mod <- try(lmer(form, data = df_alive), silent = TRUE)
  
  if(inherits(mod, "try-error")) return(NULL)
  
  vc <- as.data.frame(VarCorr(mod))
  
  data.frame(
    Trait = trait,
    Var_Family = vc$vcov[1],
    Var_Within = vc$vcov[2],
    Repeatability_R = vc$vcov[1] / (vc$vcov[1] + vc$vcov[2])
  )
}

write_csv(
  do.call(rbind, lapply(traits, get_varcomp)),
  file.path(out_dir, "02_Variance_Components_Repeatability.csv")
)


# ==============================================================================
# 5. Multivariate Phenotypic Structure (PCA without Redundant Scaling)----
# ==============================================================================
pca_res <- prcomp(
  df_alive %>% select(all_of(traits)),
  center = TRUE,
  scale. = TRUE
)

df_alive <- df_alive %>%
  mutate(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2])

## ---- Fig.2A Data: Eigenvalues and variance explained ----
pca_summary <- summary(pca_res)

eig_data <- data.frame(
  PC = 1:length(pca_summary$importance[1,]),
  StdDev = pca_summary$importance[1,],
  Variance_Proportion = pca_summary$importance[2,],
  Cumulative_Proportion = pca_summary$importance[3,]
)

write_csv(eig_data, file.path(out_dir, "Fig2A_Data.csv"))

p_scree <- fviz_eig(
  pca_res,
  addlabels = TRUE,
  ylim = c(0, 100),
  main = "A. Variance explained by principal components (Scree plot)",
  ggtheme = academic_theme
)

## ---- Fig.2B Data: Individual scores ----
write_csv(
  df_alive %>% select(Family_ID, PC1, PC2, all_of(traits)),
  file.path(out_dir, "Fig2B_Data.csv")
)

p_pca_ind <- fviz_pca_ind(
  pca_res,
  geom = "point",
  col.ind = as.factor(df_alive$Family_ID),
  addEllipses = FALSE,
  legend.title = "Family",
  pointsize = 1.5,
  title = "B. PCA distribution of derived population",
  ggtheme = academic_theme
) + theme(legend.position = "none")

ggsave(
  file.path(out_dir, "Fig2_PCA_Structure.pdf"),
  p_scree + p_pca_ind,
  width = 12, height = 6, dpi = 600
)


# ==============================================================================
# 6. Neyman Optimal Allocation + Quantile-Based Deterministic Sampling----
# ==============================================================================
cat("Executing core algorithm: multivariate variance-weighted stratified sampling...\n")

family_stats <- df_alive %>%
  group_by(Family_ID) %>%
  summarise(
    N_h = n(),
    Var_PC1 = var(PC1, na.rm = TRUE),
    Var_PC2 = var(PC2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Var_PC1 = ifelse(is.na(Var_PC1), 0, Var_PC1),
    Var_PC2 = ifelse(is.na(Var_PC2), 0, Var_PC2),
    S_h_raw = sqrt(Var_PC1 + Var_PC2)
  )

base_weight <- median(family_stats$S_h_raw[family_stats$S_h_raw > 0], na.rm = TRUE) * 0.1

family_stats <- family_stats %>%
  mutate(
    S_h = ifelse(S_h_raw == 0, base_weight, S_h_raw),
    Weight = N_h * S_h
  )

smart_round <- function(x, target = 200) {
  res <- floor(x)
  diff <- target - sum(res)
  
  if(diff > 0) {
    idx <- order(x - res, decreasing = TRUE)
    idx <- idx[seq_len(min(diff, length(idx)))]
    res[idx] <- res[idx] + 1
  }
  
  return(res)
}

family_stats$Target_N <- smart_round(
  200 * family_stats$Weight / sum(family_stats$Weight)
)

set.seed(2026)

selected_plants <- data.frame()

for (fam in family_stats$Family_ID) {
  
  fam_data <- df_alive %>% filter(Family_ID == fam)
  n_alloc <- family_stats %>% filter(Family_ID == fam) %>% pull(Target_N)
  
  if (n_alloc == 0) next
  
  if (n_alloc >= nrow(fam_data)) {
    selected_plants <- bind_rows(selected_plants, fam_data)
    next
  }
  
  k_clusters <- min(3, nrow(fam_data))
  
  fam_data <- fam_data %>%
    mutate(Gradient = as.factor(ntile(PC1, k_clusters)))
  
  grad_counts <- as.data.frame(table(fam_data$Gradient))
  
  grad_alloc <- smart_round(
    n_alloc * (grad_counts$Freq / sum(grad_counts$Freq)),
    target = n_alloc
  )
  
  fam_sampled <- data.frame()
  
  for (i in 1:k_clusters) {
    grad_pool <- fam_data %>% filter(Gradient == grad_counts$Var1[i])
    n_to_pick <- min(grad_alloc[i], nrow(grad_pool))
    
    if(n_to_pick > 0) {
      fam_sampled <- bind_rows(
        fam_sampled,
        grad_pool[sample(nrow(grad_pool), n_to_pick), ]
      )
    }
  }
  
  shortfall <- n_alloc - nrow(fam_sampled)
  
  if(shortfall > 0) {
    rem_pool <- fam_data %>%
      filter(!Plant_ID %in% fam_sampled$Plant_ID)
    
    safe_shortfall <- min(shortfall, nrow(rem_pool))
    
    fam_sampled <- bind_rows(
      fam_sampled,
      rem_pool[sample(nrow(rem_pool), safe_shortfall), ]
    )
  }
  
  selected_plants <- bind_rows(selected_plants, fam_sampled)
}

write_csv(selected_plants, file.path(out_dir, "03_Selected_200_GWAS.csv"))

write_csv(
  selected_plants %>%
    select(Family_ID, Plant_ID, PC1, Height_2026, Internode_2026, Branch_2026, Multi_2026),
  file.path(out_dir, "Fig3_Selected_200_Plants_Details.csv")
)


# ==============================================================================
# 7. Sampling Representativeness Diagnostics (Publication Layout Optimized)----
# ==============================================================================
## ---- Fig.3A Data: Neyman allocation ----
write_csv(
  family_stats %>% select(Family_ID, N_h, S_h, Weight, Target_N),
  file.path(out_dir, "Fig3A_Data.csv")
)

p_neyman <- ggplot(family_stats, aes(
  x = reorder(Family_ID, Target_N),
  y = Target_N
)) +
  geom_bar(stat = "identity", fill = "#4DAF4A", alpha = 0.8) +
  coord_flip() +
  labs(
    x = "Inbred family",
    y = "Allocated sample size (Neyman)",
    title = "A. Allocation based on population size and multivariate dispersion"
  ) +
  academic_theme +
  theme(axis.text.y = element_text(size = 6))

## ---- Fig.3B Data: Density comparison ----
density_compare_data <- bind_rows(
  df_alive %>% mutate(Group = "Total population") %>% select(PC1, Group),
  selected_plants %>% mutate(Group = "Core subset (n = 200)") %>% select(PC1, Group)
)

write_csv(density_compare_data, file.path(out_dir, "Fig3B_Data.csv"))

p_density_compare <- ggplot() +
  geom_density(data = df_alive, aes(x = PC1, fill = "Total population"), alpha = 0.4) +
  geom_density(data = selected_plants, aes(x = PC1, fill = "Core subset (n = 200)"), alpha = 0.5) +
  scale_fill_manual(
    name = "Population",
    values = c("Total population" = "grey50", "Core subset (n = 200)" = "#E41A1C")
  ) +
  labs(
    x = "Principal Component 1 (PC1)",
    y = "Density",
    title = "B. Representativeness of sampled core subset"
  ) +
  academic_theme +
  theme(legend.position = "bottom")

## ---- Fig.3C Data: Within-family sampling distribution ----
strip_data <- df_alive %>%
  mutate(Selected = ifelse(Plant_ID %in% selected_plants$Plant_ID, "Selected", "Not selected")) %>%
  select(Family_ID, PC1, Selected)

write_csv(strip_data, file.path(out_dir, "Fig3C_Data.csv"))

p_strip <- ggplot(df_alive, aes(
  x = PC1,
  y = reorder(Family_ID, PC1, FUN=median)
)) +
  geom_jitter(color = "grey70", alpha = 0.6, height = 0.2, size = 1) +
  geom_point(
    data = selected_plants,
    aes(x = PC1, y = reorder(Family_ID, PC1, FUN=median)),
    color = "#E41A1C",
    size = 1.8,
    alpha = 0.9
  ) +
  labs(
    x = "Principal Component 1 (PC1)",
    y = "Inbred family",
    title = "C. Within-family quantile sampling across phenotypic gradients"
  ) +
  academic_theme +
  theme(axis.text.y = element_text(size = 6))

ggsave(
  file.path(out_dir, "Fig3_Sampling_Justification.pdf"),
  (p_neyman | p_density_compare) / p_strip,
  width = 16, height = 14, dpi = 600
)

cat("Pipeline completed. All outputs and figure-level datasets have been exported.\n")

# ==============================================================================
