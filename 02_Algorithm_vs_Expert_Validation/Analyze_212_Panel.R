# Expert-Adjusted GWAS Core Panel Composition, Balance, and Representativeness Analysis (Final 212 Plants) ----

# Load required packages
library(tidyverse)
library(readxl)
library(patchwork)
library(extrafont)

# Ensure fonts are available (may require loadfonts(device = "win") on first run)
# loadfonts(device = "win")

# 0. Global Parameters and Theme Settings ----

# Low-saturation, colorblind-friendly palette
cb_palette <- c("#6E9CBD", "#D59C68", "#85B08A", "#D6C27A", "#A58FAD")
COLOR_200 <- cb_palette[2]   # Original candidate pool
COLOR_FINAL <- cb_palette[1] # Final confirmed cohort

theme_gwas <- theme_classic(base_size = 12, base_family = "Arial") +
  theme(
    plot.title = element_blank(),
    axis.text = element_text(color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    legend.position = "top",
    legend.title = element_blank()
  )

# 1. Directory and Path Configuration ----
# Use relative paths; ensure data files are placed under the input data directory
data_dir <- "./01_Raw_Data"                 # Input data directory
output_root <- "./05_Analysis_Outputs"      # Output root directory

# Input files (adjust according to actual file names)
original_file <- file.path(data_dir, "Original_Selected_200_WithTraits_Labeled.xlsx")
latest_file   <- file.path(data_dir, "Latest_Selected_For_Height+Color.xlsx")
gwas_db       <- file.path(data_dir, "data_202605_gwas.csv")

# Create a timestamped output subdirectory
timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
output_dir <- file.path(output_root, paste0("Panel_Diagnostics_212_", timestamp))
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 2. Data Loading and Cleaning ----

## 2.1 Load original 200-plant data (auto-selected results) ----
df_200 <- read_excel(original_file, sheet = "03_Selected_200_WithTraits_Labe") %>%
  select(Family_ID, Plant_ID, PC1, PC2, Plant_Height_May) %>%
  mutate(across(c(Family_ID, Plant_ID), str_trim))

## 2.2 Load latest expert-adjusted list ----
# Expected 212 plants; read all without truncation
df_latest_raw <- read_excel(latest_file, sheet = "Latest_Selected") %>%
  select(1, 4) %>%
  set_names(c("Plant_ID", "Height_May")) %>%
  filter(!is.na(Plant_ID) & str_trim(Plant_ID) != "") %>%
  mutate(
    Plant_ID = str_trim(Plant_ID),
    Height_May = as.numeric(Height_May)
  ) %>%
  distinct(Plant_ID, .keep_all = TRUE)

# Directly use all rows as the final cohort (assuming exactly 212 plants)
df_final <- df_latest_raw

# Validate sample size
if (nrow(df_final) != 212) {
  warning(sprintf("Warning: Final cohort sample size = %d, expected 212. Please verify the data!", nrow(df_final)))
}

# 3. Difference Comparison (Added vs. Removed) ----
common <- intersect(df_200$Plant_ID, df_final$Plant_ID)
added <- setdiff(df_final$Plant_ID, df_200$Plant_ID)
removed <- setdiff(df_200$Plant_ID, df_final$Plant_ID)

# 4. Family Mapping ----

# Extract mapping from the original 200 plants
map_original <- df_200 %>% select(Plant_ID, Family_ID) %>% deframe()

# Fill in families for newly added individual plants from the database
df_gwas <- read_csv(gwas_db, show_col_types = FALSE) %>%
  rename(Family_ID = 1, Plant_ID = 2) %>%
  mutate(across(c(Family_ID, Plant_ID), str_trim)) %>%
  filter(!is.na(Family_ID) & Family_ID != "")

map_gwas <- df_gwas %>% 
  filter(Plant_ID %in% added) %>%
  distinct(Plant_ID, .keep_all = TRUE) %>%
  select(Plant_ID, Family_ID) %>% 
  deframe()

# Manual mapping dictionary (expand as needed based on actual data)
map_manual <- c(
  AL06 = "Algonquin", AL07 = "Algonquin", AL08 = "Algonquin", AL10 = "Algonquin",
  AL18 = "Algonquin", AL2 = "Algonquin", AL3 = "Algonquin", AL9 = "Algonquin",
  EL8 = "Commercial_A", EM8 = "Commercial_A", EN8 = "Commercial_A",
  HY20 = "Huaiyang_Parent", HY201 = "Huaiyang_Parent", HY21 = "Huaiyang_Parent"
)

# Merge all dictionaries and apply mapping
df_final <- df_final %>%
  mutate(
    Family_ID = coalesce(
      map_manual[Plant_ID],
      map_gwas[Plant_ID],
      map_original[Plant_ID],
      "Unknown"
    )
  )

# 5. Statistical Calculations (Gini & K-S Test) ----

# Family distribution
family_dist <- df_final %>%
  count(Family_ID, name = "Count") %>%
  arrange(desc(Count))

# Gini coefficient calculation function
calc_gini <- function(x) {
  x <- x[x > 0]
  if (length(x) <= 1) return(0)
  x <- sort(x)
  n <- length(x)
  num <- sum((2 * 1:n - n - 1) * x)
  return(num / (n * sum(x)))
}
gini_val <- calc_gini(family_dist$Count)

# Kolmogorov-Smirnov test (compare plant height distribution between original 200 and final 212)
height_orig <- drop_na(df_200, Plant_Height_May)$Plant_Height_May
height_fin  <- drop_na(df_final, Height_May)$Height_May

ks_res <- ks.test(height_fin, height_orig)

# 6. Data Visualization (Patchwork) ----

# Panel A: Bar chart of family distribution
p_bar <- ggplot(family_dist, aes(x = reorder(Family_ID, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = COLOR_FINAL, alpha = 0.85,
           color = "white", linewidth = 0.3, width = 0.7) +
  geom_hline(yintercept = mean(family_dist$Count), color = COLOR_200,
             linetype = "dashed", linewidth = 1) +
  annotate("text", x = 1, y = family_dist$Count[1] + 1.5,
           label = paste("n =", family_dist$Count[1]),
           color = COLOR_FINAL, size = 3, fontface = "bold") +
  labs(x = "Family ID", y = "Number of Selected Individuals") +
  theme_gwas +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

# Panel B: Density distribution comparison (Original vs. Final)
df_density <- bind_rows(
  tibble(Height = height_orig, Panel = "Original Candidate Pool (n=200)"),
  tibble(Height = height_fin,  Panel = "Final Expert-Adjusted Cohort (n=212)")
) %>%
  mutate(Panel = factor(Panel, levels = c("Original Candidate Pool (n=200)",
                                          "Final Expert-Adjusted Cohort (n=212)")))

p_density <- ggplot(df_density, aes(x = Height, fill = Panel, color = Panel)) +
  geom_density(alpha = 0.35, linewidth = 1.2) +
  scale_fill_manual(values = c(COLOR_200, COLOR_FINAL)) +
  scale_color_manual(values = c(COLOR_200, COLOR_FINAL)) +
  labs(x = "Plant Height in May (cm)", y = "Density") +
  annotate("text", x = quantile(df_density$Height, 0.75, na.rm = TRUE),
           y = Inf, vjust = 2,
           label = sprintf("K-S test: D = %.4f, p = %.4f", ks_res$statistic, ks_res$p.value),
           size = 3.5, fontface = "italic", color = "grey30") +
  theme_gwas

# Combine plots (A | B)
combined_plot <- (p_bar | p_density) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 14, face = "bold", family = "Arial"))

# 7. Output and Archiving ----

# Save vector graphic
ggsave(file.path(output_dir, "Figure_Panel_Diagnostics_212.pdf"),
       plot = combined_plot, width = 12, height = 5.5, device = cairo_pdf)

# Save statistics table and difference lists
write_csv(family_dist, file.path(output_dir, "Family_Distribution_212.csv"))
write_lines(added,   file.path(output_dir, "01_Added_Plant_IDs.txt"))
write_lines(removed, file.path(output_dir, "01_Removed_Plant_IDs.txt"))

# Export K-S test results to a text file
sink(file.path(output_dir, "KS_Test_Result.txt"))
cat("Kolmogorov-Smirnov Test: Original (n=200) vs Final (n=212)\n")
print(ks_res)
sink()

message("[SUCCESS] Analysis complete! All outputs have been saved to: ", output_dir)
