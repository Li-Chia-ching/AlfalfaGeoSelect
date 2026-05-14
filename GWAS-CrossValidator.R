# =============================================================================
# Automated Crossing Pair Recommendation Script (For Post-GWAS Validation)
# Input: The 200 samples selected from the previous pipeline
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

# 1. Read the 200 core samples that have already been selected
#    (Make sure the file name matches the one actually generated in your pipeline)
# df_200 <- read_csv("GWAS_Pipeline_Robust200_XXXX/02_Selected_200_GWAS_FullTraits.csv", show_col_types = FALSE)

# [Simulated data]
# If you want to test the code now, you can run the following block to generate mock data.
# If you are using real data, comment out this section.
set.seed(123)
df_200 <- data.frame(
  ID = paste0("P_", 1:200),
  Family = sample(paste0("Fam_", 1:20), 200, replace = TRUE),
  Plant_Height = rnorm(200, mean = 150, sd = 20),
  MF_Total = round(runif(200, 0, 50))
)

# 2. Define the crossing-pair generation function
generate_crosses <- function(df, trait_col, pool_size = 15, top_n_pairs = 5) {
  
  # Sort the data in descending order by the target trait
  df_sorted <- df %>% dplyr::arrange(desc(!!sym(trait_col)))
  
  # Extract the strong candidate pool and the weak candidate pool
  strong_pool <- head(df_sorted, pool_size)
  weak_pool <- tail(df_sorted, pool_size)
  
  # ==========================================================
  # Strategy A: Strong × Weak
  # Purpose: Build an extreme phenotype segregation population
  # ==========================================================
  sxw_pairs <- expand.grid(P1_ID = strong_pool$ID, P2_ID = weak_pool$ID, stringsAsFactors = FALSE) %>%
    dplyr::left_join(strong_pool %>% dplyr::select(ID, Family, !!sym(trait_col)), by = c("P1_ID" = "ID")) %>%
    dplyr::rename(P1_Family = Family, P1_Value = !!sym(trait_col)) %>%
    dplyr::left_join(weak_pool %>% dplyr::select(ID, Family, !!sym(trait_col)), by = c("P2_ID" = "ID")) %>%
    dplyr::rename(P2_Family = Family, P2_Value = !!sym(trait_col)) %>%
    # Calculate phenotypic difference
    dplyr::mutate(
      Pheno_Diff = abs(P1_Value - P2_Value),
      Cross_Type = "Strong x Weak",
      Target_Trait = trait_col
    ) %>%
    # Sort by phenotypic difference from largest to smallest
    dplyr::arrange(desc(Pheno_Diff)) %>%
    # [Fix: bidirectional deduplication]
    # Ensure that each strong parent (P1) and each weak parent (P2) appears only once in the final list
    dplyr::distinct(P1_ID, .keep_all = TRUE) %>%
    dplyr::distinct(P2_ID, .keep_all = TRUE) %>%  # Critical fix: prevent repeated use of weak parents
    head(top_n_pairs)
  
  # ==========================================================
  # Strategy B: Strong × Strong
  # Purpose: Allelic testing / aggregation of superior alleles
  # ==========================================================
  sxs_pairs <- expand.grid(P1_ID = strong_pool$ID, P2_ID = strong_pool$ID, stringsAsFactors = FALSE) %>%
    # Exclude self-crosses and remove duplicate combinations such as A×B and B×A
    dplyr::filter(P1_ID < P2_ID) %>% 
    dplyr::left_join(strong_pool %>% dplyr::select(ID, Family, !!sym(trait_col)), by = c("P1_ID" = "ID")) %>%
    dplyr::rename(P1_Family = Family, P1_Value = !!sym(trait_col)) %>%
    dplyr::left_join(strong_pool %>% dplyr::select(ID, Family, !!sym(trait_col)), by = c("P2_ID" = "ID")) %>%
    dplyr::rename(P2_Family = Family, P2_Value = !!sym(trait_col)) %>%
    # [Core biological constraint]: strictly require the two strong parents to come from different families
    dplyr::filter(P1_Family != P2_Family) %>%
    dplyr::mutate(
      Pheno_Diff = abs(P1_Value - P2_Value),
      Cross_Type = "Strong x Strong",
      Target_Trait = trait_col,
      # For strong × strong crosses, prioritize combinations with the largest summed trait values
      Combined_Strength = P1_Value + P2_Value 
    ) %>%
    dplyr::arrange(desc(Combined_Strength)) %>%
    dplyr::select(-Combined_Strength) %>%
    head(top_n_pairs)
  
  # Merge the results from both strategies
  result <- dplyr::bind_rows(sxw_pairs, sxs_pairs)
  return(result)
}

# 3. Generate recommended crossing lists for "Plant Height" and "Multi-leaf Number"
#    (Recommend 5 pairs per strategy for each trait)
cat("Generating recommended crossing combinations for Plant Height...\n")
crosses_height <- generate_crosses(df_200, trait_col = "Plant_Height", pool_size = 15, top_n_pairs = 5)

cat("Generating recommended crossing combinations for Multi-leaf Number...\n")
crosses_mf <- generate_crosses(df_200, trait_col = "MF_Total", pool_size = 15, top_n_pairs = 5)

# 4. Merge the final list and save it
final_crossing_plan <- dplyr::bind_rows(crosses_height, crosses_mf)

# Reorder the columns for easier reading
final_crossing_plan <- final_crossing_plan %>%
  dplyr::select(Target_Trait, Cross_Type, 
                P1_ID, P1_Family, P1_Value, 
                P2_ID, P2_Family, P2_Value, Pheno_Diff)

# Print a preview to the console
print(final_crossing_plan)

# Export as CSV for greenhouse / field staff
write_csv(final_crossing_plan, "03_Recommended_Crossing_Plan.csv")
cat("\nThe crossing plan has been saved to '03_Recommended_Crossing_Plan.csv'\n")
