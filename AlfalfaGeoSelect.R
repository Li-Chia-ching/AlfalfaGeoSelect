# ============================================================================
# Alfalfa Selection V3.0: Publication-Ready Visualization
# Clean implementation with standardized annotations and vector-based graphics
# ============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

# ------------------------------
# 1. Environment setup
# ------------------------------
required_packages <- c("dplyr", "tidyr", "ggplot2", "stringr", "readr", "metR")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages, repos = "https://cran.rstudio.com/")

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(readr)

# ------------------------------
# 2. Data loading and cleaning
# ------------------------------
clean_pheno_data <- function(df_raw, value_name) {
  colnames(df_raw)[1] <- "Row_ID"
  
  df_long <- df_raw %>%
    filter(str_detect(Row_ID, "^L\\d+$")) %>%
    mutate(across(-Row_ID, as.character)) %>%
    pivot_longer(cols = -Row_ID, names_to = "Col_ID", values_to = "Value") %>%
    mutate(
      Value = suppressWarnings(as.numeric(Value)),
      Plant_ID = paste0(Row_ID, "-", Col_ID)
    ) %>%
    rename(!!value_name := Value)
  
  return(df_long)
}

# ---- Load plant height data ----
if (exists("initial_flowering_2025")) {
  raw_height <- get("initial_flowering_2025")
} else if (file.exists("Rawdata-PlantHeight_202504.csv")) {
  raw_height <- read_csv("Rawdata-PlantHeight_202504.csv", show_col_types = FALSE)
} else {
  stop("Plant height data not found")
}

# ---- Load multifoliate data ----
if (exists("Alfalfa_Multi_202504")) {
  raw_multi <- get("Alfalfa_Multi_202504")
} else if (file.exists("Rawdata-Multifoliate_202504.csv")) {
  raw_multi <- read_csv("Rawdata-Multifoliate_202504.csv", show_col_types = FALSE)
} else {
  stop("Multifoliate data not found")
}

# ---- Construct full grid and merge ----
all_rows <- paste0("L", 1:50)
all_cols <- LETTERS[1:21]

full_grid <- expand.grid(Row_ID = all_rows, Col_ID = all_cols) %>%
  mutate(Plant_ID = paste0(Row_ID, "-", Col_ID))

df_analysis <- full_grid %>%
  left_join(clean_pheno_data(raw_height, "Height"), by = c("Row_ID", "Col_ID", "Plant_ID")) %>%
  left_join(clean_pheno_data(raw_multi, "Multi_Score"), by = c("Row_ID", "Col_ID", "Plant_ID"))

# ------------------------------
# 3. Index calculation and selection
# ------------------------------
df_analysis <- df_analysis %>%
  mutate(Index = (0.5 * Height) * (0.5 * Multi_Score))

selection_list <- df_analysis %>%
  filter(!is.na(Index)) %>%
  arrange(desc(Index)) %>%
  slice_head(n = 50)

cutoff_value <- min(selection_list$Index)

df_analysis <- df_analysis %>%
  mutate(
    Status = case_when(
      Plant_ID %in% selection_list$Plant_ID ~ "Selected (Top 50)",
      is.na(Index) ~ "Missing",
      TRUE ~ "Not Selected"
    )
  )

# ------------------------------
# 4. Contour visualization
# ------------------------------

# Background grid for contour surface
grid_x <- seq(min(df_analysis$Height, na.rm = TRUE),
              max(df_analysis$Height, na.rm = TRUE),
              length.out = 100)

grid_y <- seq(1, 5, length.out = 100)

grid_map <- expand.grid(Height = grid_x, Multi_Score = grid_y) %>%
  mutate(Index = (0.5 * Height) * (0.5 * Multi_Score))

p_contour <- ggplot(df_analysis %>% filter(Status != "Missing"),
                    aes(x = Height, y = Multi_Score)) +
  
  # Contour background
  geom_contour(data = grid_map,
               aes(z = Index, color = ..level..),
               bins = 12, linewidth = 0.3, alpha = 0.4) +
  scale_color_viridis_c(option = "viridis", name = "Selection Index") +
  
  # Selection cutoff line
  geom_contour(data = grid_map,
               aes(z = Index),
               breaks = cutoff_value,
               color = "#E41A1C",
               linewidth = 1,
               linetype = "longdash") +
  
  # Data points
  geom_jitter(aes(fill = Status, size = Status),
              shape = 21,
              stroke = 0.2,
              width = 0.4,
              height = 0.2,
              alpha = 0.85) +
  
  # Styling
  scale_fill_manual(values = c(
    "Not Selected" = "grey85",
    "Selected (Top 50)" = "#d35400"
  )) +
  scale_size_manual(values = c(
    "Not Selected" = 2,
    "Selected (Top 50)" = 4
  )) +
  
  # Selection gradient arrow
  annotate("segment",
           x = max(grid_x) * 0.85,
           xend = max(grid_x),
           y = 1.5,
           yend = 2.0,
           arrow = arrow(type = "closed", length = unit(0.3, "cm"))) +
  
  annotate("text",
           x = max(grid_x) * 0.92,
           y = 1.4,
           label = "Selection Gradient",
           size = 3.5,
           fontface = "bold") +
  
  # Cutoff label
  annotate("text",
           x = max(grid_x),
           y = 4.8,
           label = paste0("Cutoff: ", round(cutoff_value, 2)),
           color = "#E41A1C",
           hjust = 1,
           fontface = "italic") +
  
  labs(
    title = "Weighted Geometric Selection",
    subtitle = "Selection based on combined trait performance",
    x = "Plant Height (cm)",
    y = "Multifoliate Score",
    caption = "Index = (0.5 × Height) × (0.5 × Multifoliate)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.title = element_text(face = "bold")
  )

# ------------------------------
# 5. Output
# ------------------------------
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- paste0("Alfalfa_GeoViz_Pro_", timestamp)
dir.create(out_dir)

ggsave(file.path(out_dir, "geometric_selection.png"),
       p_contour, width = 10, height = 7, dpi = 300)

write_csv(selection_list,
          file.path(out_dir, "selection_list.csv"))

cat("\n========================================\n")
cat("Visualization completed\n")
cat("Output directory:", out_dir, "\n")
cat("========================================\n")
