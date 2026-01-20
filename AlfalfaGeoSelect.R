# ==============================================================================
# Alfalfa Selection V3.0: Publication-Ready Visualization
# 修正：移除特殊字符，使用专业矢量箭头和术语
# ==============================================================================

# 1. 环境准备 ------------------------------------------------------------------
required_packages <- c("dplyr", "tidyr", "ggplot2", "stringr", "readr", "metR") 
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, repos="https://cran.rstudio.com/")

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(readr)

# 2. 数据加载与清洗 ------------------------------------------------------------
clean_pheno_data <- function(df_raw, value_name) {
  colnames(df_raw)[1] <- "Row_ID"
  df_long <- df_raw %>%
    filter(str_detect(Row_ID, "^L\\d+$")) %>%
    mutate(across(-Row_ID, as.character)) %>%
    pivot_longer(cols = -Row_ID, names_to = "Col_ID", values_to = "Value") %>%
    mutate(Value = suppressWarnings(as.numeric(Value)), Plant_ID = paste0(Row_ID, "-", Col_ID)) %>%
    rename(!!value_name := Value)
  return(df_long)
}

# 智能加载
if(exists("initial_flowering_2025")) {
  raw_height <- get("initial_flowering_2025")
} else if(file.exists("Rawdata-PlantHeight_202504.csv")) {
  raw_height <- read_csv("Rawdata-PlantHeight_202504.csv", show_col_types = FALSE)
} else { stop("缺少株高数据") }

if(exists("Alfalfa_Multi_202504")) {
  raw_multi <- get("Alfalfa_Multi_202504")
} else if(file.exists("Rawdata-Multifoliate_202504.csv")) {
  raw_multi <- read_csv("Rawdata-Multifoliate_202504.csv", show_col_types = FALSE)
} else { stop("缺少多叶数据") }

# 合并
all_rows <- paste0("L", 1:50)
all_cols <- LETTERS[1:21]
full_grid <- expand.grid(Row_ID = all_rows, Col_ID = all_cols) %>% mutate(Plant_ID = paste0(Row_ID, "-", Col_ID))

df_analysis <- full_grid %>%
  left_join(clean_pheno_data(raw_height, "Height"), by = c("Row_ID", "Col_ID", "Plant_ID")) %>%
  left_join(clean_pheno_data(raw_multi, "Multi_Score"), by = c("Row_ID", "Col_ID", "Plant_ID"))

# 3. 计算指数与筛选 ------------------------------------------------------------
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
      is.na(Index) ~ "Dead/Missing",
      TRUE ~ "Not Selected"
    )
  )

# 4. 可视化：专业版等值线图 (Professional Contour Plot) ------------------------

# 生成背景网格
grid_x <- seq(min(df_analysis$Height, na.rm=T), max(df_analysis$Height, na.rm=T), length.out = 100)
grid_y <- seq(1, 5, length.out = 100)
grid_map <- expand.grid(Height = grid_x, Multi_Score = grid_y) %>%
  mutate(Index = (0.5 * Height) * (0.5 * Multi_Score))

p_contour <- ggplot(df_analysis %>% filter(Status != "Dead/Missing"), aes(x = Height, y = Multi_Score)) +
  
  # A. 背景等值线
  geom_contour(data = grid_map, aes(z = Index, color = ..level..), bins = 12, size = 0.3, alpha = 0.4) +
  scale_color_viridis_c(option = "viridis", name = "Selection Index") +
  
  # B. 红色阈值虚线
  geom_contour(data = grid_map, aes(z = Index), breaks = cutoff_value, color = "#E41A1C", size = 1, linetype = "longdash") +
  
  # C. 数据点
  geom_jitter(aes(fill = Status, size = Status), shape = 21, stroke = 0.2, width = 0.4, height = 0.2, alpha = 0.85) +
  
  # D. 样式调整
  scale_fill_manual(values = c("Not Selected" = "grey85", "Selected (Top 50)" = "#d35400")) +
  scale_size_manual(values = c("Not Selected" = 2, "Selected (Top 50)" = 4)) +
  
  # E. 专业注释 (替换原有的特殊字符)
  # 使用 geom_segment 绘制真实的矢量箭头
  annotate("segment", x = max(grid_x)*0.85, xend = max(grid_x), y = 1.5, yend = 2.0, 
           arrow = arrow(type = "closed", length = unit(0.3, "cm")), color = "black") +
  annotate("text", x = max(grid_x)*0.92, y = 1.4, label = "Selection Gradient", 
           color = "black", hjust = 0.5, size = 3.5, fontface = "bold") +
  
  # 标注阈值分
  annotate("text", x = max(grid_x), y = 4.8, label = paste0("Cutoff Score: ", round(cutoff_value, 2)), 
           color = "#E41A1C", hjust = 1, fontface = "italic") +
  
  labs(
    title = "Weighted Geometric Selection Analysis",
    subtitle = "Points above the red dashed line are selected based on Factor Analysis weights",
    x = "Plant Height (cm)",
    y = "Multifoliate Score (1-5)",
    caption = "Index = (0.5 * Height) * (0.5 * Multi)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.title = element_text(face = "bold", linewidth = 14)
  )

# 5. 保存结果 ------------------------------------------------------------------
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- paste0("Alfalfa_GeoViz_Pro_", timestamp)
dir.create(out_dir)

ggsave(file.path(out_dir, "Plot_Geometric_Selection_Pro.png"), p_contour, width = 10, height = 7, dpi = 300)
write_csv(selection_list, file.path(out_dir, "01_Selection_List.csv"))

cat("\n======================================================\n")
cat("专业版图表已生成！\n")
cat("修改点：移除特殊字符，增加 'Selection Gradient' 矢量箭头指示。\n")
cat("输出目录:", out_dir, "\n")
cat("======================================================\n")
