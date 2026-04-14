# Genetic Variation Analysis + Stratified Within-Family GWAS Sampling R Script Based on Multidimensional Gradients (v5.3 - Modified)
# Features: Multivariate variance-based Neyman allocation, quantile-gradient deterministic sampling, 
# LMM-based repeatability estimation, and publication-ready visualization

# ------------------- 0. User Parameters (User parameter configuration section) -------------------
TARGET_TOTAL_N <- 140  # Requirement 1: Set the total number of individuals to be selected
MAX_PER_FAMILY <- 20   # Requirement 2: Set the maximum number of individuals selected per family (n <= 20)

# ------------------- 1. Check and load required packages -------------------
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

# ------------------- 2. Global font configuration and academic theme setup -------------------
showtext_auto() 
if (.Platform$OS.type == "windows") {
  font_add("SimSun", "simsun.ttc") 
} else {
  font_add("SimSun", "STSong.ttf") 
}

out_dir <- paste0("GWAS_Sampling_Two-Year_", format(Sys.Date(), "%Y%m%d"))
dir.create(out_dir, showWarnings = FALSE)

academic_theme <- theme_bw(base_family = "SimSun") + 
  theme(panel.grid.minor = element_blank(),
        text = element_text(family = "SimSun", size = 12),
        axis.title = element_text(family = "SimSun", face = "bold"),
        plot.title = element_text(family = "SimSun", face = "bold", size = 14),
        strip.text = element_text(size = 9, margin = margin(2, 2, 2, 2)),
        panel.spacing = unit(1, "lines"))

# ------------------- 3. Data processing (merge two-year phenotypes) -------------------
cat("Reading and preprocessing two-year datasets...\n")
df_2025 <- read_csv("data_2025.csv", show_col_types = FALSE) %>%
  rename(Family_ID = Parent, Plant_ID = SampleID, Height_2025 = PlantHeight, Internode_2025 = NodeCount) %>%
  select(Family_ID, Plant_ID, Real_Count, Height_2025, Internode_2025)

df_2026 <- read_csv("data_2026.csv", show_col_types = FALSE) %>%
  filter(!is.na(Group) & trimws(Group) != "") %>%
  rename(Family_ID = Group, Plant_ID = Matrix, Height_2026 = Plant_Height, 
         Internode_2026 = Internode, Branch_2026 = Branch_Number, Multi_2026 = Multifoliate_Score) %>%
  select(Family_ID, Plant_ID, Height_2026, Internode_2026, Branch_2026, Multi_2026, Death_Code)

df_wide <- inner_join(df_2025, df_2026, by = c("Family_ID", "Plant_ID")) %>%
  mutate(Status = case_when(Death_Code == 1 ~ "Dead_in_2025",
                            Death_Code == 2 ~ "Dead_in_2026",
                            TRUE ~ "Alive"),
         Delta_Height = Height_2026 - Height_2025)

write_csv(df_wide, file.path(out_dir, "01_Merged_Wide_Data_All.csv"))

# Define trait vectors by year
traits_2025 <- c("Height_2025", "Internode_2025")
traits_2026 <- c("Height_2026", "Internode_2026", "Branch_2026", "Multi_2026")
all_traits <- c(traits_2025, traits_2026)   # Total of 6 traits

# Filter alive individuals and ensure no missing values across all traits
df_alive <- df_wide %>% filter(Status == "Alive") %>%
  drop_na(all_of(all_traits))

# ------------------- 4. Phenotypic segregation evidence (LMM repeatability) -------------------
cat("Generating phenotypic segregation and variance structure...\n")

# Fig 1A data
write_csv(df_alive %>% select(Family_ID, Height_2026), 
          file.path(out_dir, "Fig1A_Data.csv"))

df_alive_for_plot <- df_alive %>%
  add_count(Family_ID, name = "n_per_family") %>%
  mutate(Family_Label = paste0(Family_ID, " (n = ", n_per_family, ")"))

p_seg_density <- ggplot(df_alive_for_plot, aes(x = Height_2026)) +
  geom_density(fill = "#377EB8", alpha = 0.6) +
  geom_rug(alpha = 0.5, sides = "b") +
  facet_wrap(~Family_Label, ncol = 8, scales = "free_y") +
  labs(x = "Plant height at early flowering in 2026 (cm)", y = "Kernel density", 
       title = "A. Density distribution of plant height segregation within selfed families") +
  academic_theme + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# Fig 1B data
write_csv(df_alive %>% select(Family_ID, Height_2026), 
          file.path(out_dir, "Fig1B_Data.csv"))

p_var_box <- ggplot(df_alive, aes(x = reorder(Family_ID, Height_2026, FUN = median), 
                                  y = Height_2026)) +
  geom_boxplot(fill = "#E41A1C", alpha = 0.6, outlier.size = 0.3) +
  labs(x = "Selfed family ID (ordered by median)", y = "Plant height at early flowering in 2026 (cm)", 
       title = "B. Between-family variation and within-family segregation range") +
  academic_theme + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))

ggsave(file.path(out_dir, "Fig1_Segregation_Analysis.pdf"), 
       p_seg_density / p_var_box, width = 16, height = 12, dpi = 600)

# Repeatability estimation
get_varcomp <- function(trait) {
  form <- as.formula(paste(trait, "~ 1 + (1|Family_ID)"))
  mod <- try(lmer(form, data = df_alive), silent = TRUE)
  if(inherits(mod, "try-error")) return(NULL)
  vc <- as.data.frame(VarCorr(mod))
  data.frame(Trait = trait, 
             Var_Family = vc$vcov[1], 
             Var_Within = vc$vcov[2], 
             Repeatability_R = vc$vcov[1] / (vc$vcov[1] + vc$vcov[2]))
}
varcomp_results <- do.call(rbind, lapply(all_traits, get_varcomp))
write_csv(varcomp_results, file.path(out_dir, "02_Variance_Components_Repeatability.csv"))

# ------------------- 5. Multidimensional phenotypic structure (stepwise PCA strategy) -------------------
cat("Performing stepwise PCA and rendering biplots...\n")

# 5.1 Fig A: PCA on cross-year shared traits
traits_common <- c("Height_2025", "Internode_2025", "Height_2026", "Internode_2026")
pca_common <- prcomp(df_alive %>% select(all_of(traits_common)), 
                     center = TRUE, scale. = TRUE)

pca_common_summary <- summary(pca_common)
eig_common <- data.frame(
  PC = 1:length(pca_common_summary$importance[1,]),
  StdDev = pca_common_summary$importance[1,],
  Variance_Proportion = pca_common_summary$importance[2,],
  Cumulative_Proportion = pca_common_summary$importance[3,]
)
write_csv(eig_common, file.path(out_dir, "Fig2A_Data_Common.csv"))

df_common_scores <- as.data.frame(pca_common$x[,1:2])
names(df_common_scores) <- c("PC1_common", "PC2_common")
df_common_export <- bind_cols(
  df_alive %>% select(Family_ID, Plant_ID, all_of(traits_common)), 
  df_common_scores
)
write_csv(df_common_export, file.path(out_dir, "Fig2B_Data_Common.csv"))

# 5.2 Fig B: PCA on all 2026 traits
traits_2026_all <- c("Height_2026", "Internode_2026", "Branch_2026", "Multi_2026")
pca_2026 <- prcomp(df_alive %>% select(all_of(traits_2026_all)), 
                   center = TRUE, scale. = TRUE)

pca_2026_summary <- summary(pca_2026)
eig_2026 <- data.frame(
  PC = 1:length(pca_2026_summary$importance[1,]),
  StdDev = pca_2026_summary$importance[1,],
  Variance_Proportion = pca_2026_summary$importance[2,],
  Cumulative_Proportion = pca_2026_summary$importance[3,]
)
write_csv(eig_2026, file.path(out_dir, "Fig2A_Data_2026.csv"))

df_2026_scores <- as.data.frame(pca_2026$x[,1:2])
names(df_2026_scores) <- c("PC1_2026", "PC2_2026")
df_2026_export <- bind_cols(
  df_alive %>% select(Family_ID, Plant_ID, all_of(traits_2026_all)), 
  df_2026_scores
)
write_csv(df_2026_export, file.path(out_dir, "Fig2B_Data_2026.csv"))

# 5.3 Comprehensive PCA (used for downstream sampling)
pca_all <- prcomp(df_alive %>% select(all_of(all_traits)), 
                  center = TRUE, scale. = TRUE)
df_alive <- df_alive %>% mutate(PC1 = pca_all$x[,1], PC2 = pca_all$x[,2])

pca_summary <- summary(pca_all)
eig_data <- data.frame(
  PC = 1:length(pca_summary$importance[1,]),
  StdDev = pca_summary$importance[1,],
  Variance_Proportion = pca_summary$importance[2,],
  Cumulative_Proportion = pca_summary$importance[3,]
)
write_csv(eig_data, file.path(out_dir, "03_Comprehensive_PCA_Eigenvalues.csv"))
write_csv(df_alive %>% select(Family_ID, Plant_ID, PC1, PC2, all_of(all_traits)), 
          file.path(out_dir, "03_Comprehensive_PCA_Scores.csv"))

# ------------------- 6. Neyman optimal allocation with limits & deterministic sampling -------------------
cat(sprintf("Executing allocation: Target %d samples, Max %d per family...\n", TARGET_TOTAL_N, MAX_PER_FAMILY))

family_stats <- df_alive %>%
  group_by(Family_ID) %>%
  summarise(N_h = n(), 
            Var_PC1 = var(PC1, na.rm = TRUE),
            Var_PC2 = var(PC2, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(Var_PC1 = ifelse(is.na(Var_PC1), 0, Var_PC1),
         Var_PC2 = ifelse(is.na(Var_PC2), 0, Var_PC2),
         S_h_raw = sqrt(Var_PC1 + Var_PC2))

base_weight <- median(family_stats$S_h_raw[family_stats$S_h_raw > 0], na.rm = TRUE) * 0.1
family_stats <- family_stats %>%
  mutate(S_h = ifelse(S_h_raw == 0, base_weight, S_h_raw),
         Weight = N_h * S_h)

# New allocation algorithm: integer allocation with upper bound (cap) handling
smart_allocate_with_cap <- function(weights, total_target, max_cap) {
  n <- length(weights)
  # Initialization: proportional floor allocation with cap constraint applied
  alloc <- floor(total_target * weights / sum(weights))
  alloc[alloc > max_cap] <- max_cap
  
  diff <- total_target - sum(alloc)
  
  # Distribute remaining slots using a D'Hondt-like priority adjustment
  while(diff > 0) {
    eligible <- alloc < max_cap
    if(!any(eligible)) {
      message("Warning: Cannot reach TARGET_TOTAL_N due to MAX_PER_FAMILY restrictions.")
      break
    }
    priority_score <- weights / (alloc + 1)
    priority_score[!eligible] <- -Inf
    best <- which.max(priority_score)
    alloc[best] <- alloc[best] + 1
    diff <- diff - 1
  }
  return(alloc)
}

# Apply allocation algorithm with cap
family_stats$Target_N <- smart_allocate_with_cap(family_stats$Weight, TARGET_TOTAL_N, MAX_PER_FAMILY)

set.seed(2026) 
selected_plants <- data.frame()
for (fam in family_stats$Family_ID) {
  fam_data <- df_alive %>% filter(Family_ID == fam)
  n_alloc <- family_stats %>% filter(Family_ID == fam) %>% pull(Target_N)
  
  if (n_alloc == 0) next
  if (n_alloc >= nrow(fam_data)) { 
    selected_plants <- bind_rows(selected_plants, fam_data); next 
  }
  
  k_clusters <- min(3, nrow(fam_data))
  fam_data <- fam_data %>% mutate(Gradient = as.factor(ntile(PC1, k_clusters)))
  
  grad_counts <- as.data.frame(table(fam_data$Gradient))
  # To ensure stratified sampling does not violate the cap, replace smart_round with simple correction so that subsamples sum to n_alloc
  grad_alloc <- floor(n_alloc * (grad_counts$Freq / sum(grad_counts$Freq)))
  grad_diff <- n_alloc - sum(grad_alloc)
  if(grad_diff > 0) {
    add_idx <- order(n_alloc * (grad_counts$Freq / sum(grad_counts$Freq)) - grad_alloc, decreasing = TRUE)[1:grad_diff]
    grad_alloc[add_idx] <- grad_alloc[add_idx] + 1
  }
  
  fam_sampled <- data.frame()
  for (i in 1:k_clusters) {
    grad_pool <- fam_data %>% filter(Gradient == grad_counts$Var1[i])
    n_to_pick <- min(grad_alloc[i], nrow(grad_pool))
    if(n_to_pick > 0) {
      fam_sampled <- bind_rows(fam_sampled, grad_pool[sample(nrow(grad_pool), n_to_pick), ])
    }
  }
  
  shortfall <- n_alloc - nrow(fam_sampled)
  if(shortfall > 0) {
    rem_pool <- fam_data %>% filter(!Plant_ID %in% fam_sampled$Plant_ID)
    safe_shortfall <- min(shortfall, nrow(rem_pool))
    fam_sampled <- bind_rows(fam_sampled, rem_pool[sample(nrow(rem_pool), safe_shortfall), ])
  }
  selected_plants <- bind_rows(selected_plants, fam_sampled)
}

write_csv(selected_plants, file.path(out_dir, sprintf("04_Selected_%d_GWAS.csv", TARGET_TOTAL_N)))
write_csv(selected_plants %>% select(Family_ID, Plant_ID, PC1, all_of(all_traits)), 
          file.path(out_dir, sprintf("Fig3_Selected_%d_Plants_Details.csv", TARGET_TOTAL_N)))

# ------------------- 7. Sampling representativeness validation -------------------
# Fig 3A data: Neyman allocation results
write_csv(family_stats %>% select(Family_ID, N_h, S_h, Weight, Target_N), 
          file.path(out_dir, "Fig3A_Data.csv"))

p_neyman <- ggplot(family_stats, aes(x = reorder(Family_ID, Target_N), y = Target_N)) +
  geom_bar(stat = "identity", fill = "#4DAF4A", alpha = 0.8) +
  coord_flip() +
  labs(x = "Selfed family", y = "Allocated sample size (Neyman)", 
       title = "A. Sample allocation based on population size and multivariate dispersion") +
  academic_theme + theme(axis.text.y = element_text(size = 6))

# ---------------------------------------------------------
# Fig 3B data: density comparison of PC1 + K-S Test
# ---------------------------------------------------------
subset_label <- paste(TARGET_TOTAL_N, "core subset")

# 1. Perform Kolmogorov-Smirnov test
ks_res <- ks.test(df_alive$PC1, selected_plants$PC1)
d_stat <- ks_res$statistic
p_val <- ks_res$p.value

# 2. Format test result text (automatically adjust p-value display)
p_text <- ifelse(p_val < 0.001, "p < 0.001", sprintf("p = %.3f", p_val))
ks_annotation <- sprintf("Two-sample K-S test:\nD = %.3f, %s", d_stat, p_text)

# 3. Prepare plotting data
density_compare_data <- bind_rows(
  df_alive %>% mutate(Group = "Total population") %>% select(PC1, Group),
  selected_plants %>% mutate(Group = subset_label) %>% select(PC1, Group)
)
write_csv(density_compare_data, file.path(out_dir, "Fig3B_Data.csv"))

fill_mapping <- c("Total population" = "grey50")
fill_mapping[subset_label] <- "#E41A1C"

# 4. Plot and add annotation
p_density_compare <- ggplot() +
  geom_density(data = df_alive, aes(x = PC1, fill = "Total population"), alpha = 0.4) +
  geom_density(data = selected_plants, aes(x = PC1, fill = subset_label), alpha = 0.5) +
  scale_fill_manual(name = "Population type", values = fill_mapping) +
  labs(x = "Cross-year composite growth (PC1)", y = "Density", 
       title = "B. Representativeness of sampled core subset") +
  academic_theme + 
  theme(legend.position = "bottom") +
  annotate("text", x = max(df_alive$PC1, na.rm = TRUE), 
           y = Inf, label = ks_annotation, 
           hjust = 1, vjust = 1.5, family = "SimSun", size = 3.5, color = "black")

# Fig 3C data: within-family sampling distribution
strip_data <- df_alive %>%
  mutate(Selected = ifelse(Plant_ID %in% selected_plants$Plant_ID, "Selected", "Not selected")) %>%
  select(Family_ID, PC1, Selected)
write_csv(strip_data, file.path(out_dir, "Fig3C_Data.csv"))

p_strip <- ggplot(df_alive, aes(x = PC1, y = reorder(Family_ID, PC1, FUN = median))) +
  geom_jitter(color = "grey70", alpha = 0.6, height = 0.2, size = 1) +
  geom_point(data = selected_plants, aes(x = PC1, y = reorder(Family_ID, PC1, FUN = median)), 
             color = "#E41A1C", size = 1.8, alpha = 0.9) +
  labs(x = "Cross-year composite growth (PC1)", y = "Selfed family", 
       title = "C. Within-family quantile sampling: selected individuals cover high/medium/low gradients") +
  academic_theme + theme(axis.text.y = element_text(size = 6))

ggsave(file.path(out_dir, "Fig3_Sampling_Justification.pdf"), 
       (p_neyman | p_density_compare) / p_strip, width = 16, height = 14, dpi = 600)

cat("Pipeline completed. Stepwise PCA biplots and full sampling results have been generated, with all figure data exported.\n")
