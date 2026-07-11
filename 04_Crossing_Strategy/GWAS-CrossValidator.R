# Automated Crossing Pair Recommendation Script (For Post-GWAS Validation) ----
# Input: The REAL 200 samples selected from the previous pipeline

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

# 1. Read the Real 200 Core Samples ----
# Note: Ensure this folder and file actually exist in your working directory
input_file <- "GWAS_Pipeline_Robust200_20260512/02_Selected_200_GWAS_FullTraits.csv"

if (!file.exists(input_file)) {
  stop("Input file not found! Please verify your working directory and ensure the folder name 'GWAS_Pipeline_Robust200_20260512' matches the actual output from your pipeline.")
}

df_200 <- read_csv(input_file, show_col_types = FALSE)

# 2. Define the Crossing-Pair Generation Function ----
generate_crosses <- function(df, trait_col, pool_size = 15, top_n_pairs = 5) {
  
  # Sort samples in descending order based on the target trait
  df_sorted <- df %>% dplyr::arrange(desc(!!sym(trait_col)))
  
  # Extract the strong and weak candidate pools
  strong_pool <- head(df_sorted, pool_size)
  weak_pool <- tail(df_sorted, pool_size)
  
  ## Strategy A: Strong × Weak ----
  # Used to construct extreme phenotype segregation populations
  sxw_pairs <- expand.grid(P1_ID = strong_pool$ID, P2_ID = weak_pool$ID, stringsAsFactors = FALSE) %>%
    dplyr::left_join(strong_pool %>% dplyr::select(ID, Family, !!sym(trait_col)), by = c("P1_ID" = "ID")) %>%
    dplyr::rename(P1_Family = Family, P1_Value = !!sym(trait_col)) %>%
    dplyr::left_join(weak_pool %>% dplyr::select(ID, Family, !!sym(trait_col)), by = c("P2_ID" = "ID")) %>%
    dplyr::rename(P2_Family = Family, P2_Value = !!sym(trait_col)) %>%
    # Calculate phenotypic differences
    dplyr::mutate(
      Pheno_Diff = abs(P1_Value - P2_Value),
      Cross_Type = "Strong x Weak",
      Target_Trait = trait_col
    ) %>%
    # Sort by phenotypic difference in descending order
    dplyr::arrange(desc(Pheno_Diff)) %>%
    # [Bidirectional deduplication]
    # Ensure each strong parent (P1) and weak parent (P2) appears only once in the final list
    dplyr::distinct(P1_ID, .keep_all = TRUE) %>%
    dplyr::distinct(P2_ID, .keep_all = TRUE) %>%  
    head(top_n_pairs)
  
  ## Strategy B: Strong × Strong ----
  # Used for allelism testing and superior allele pyramiding
  sxs_pairs <- expand.grid(P1_ID = strong_pool$ID, P2_ID = strong_pool$ID, stringsAsFactors = FALSE) %>%
    # Exclude self-crosses and remove duplicate combinations (e.g., A×B and B×A)
    dplyr::filter(P1_ID < P2_ID) %>% 
    dplyr::left_join(strong_pool %>% dplyr::select(ID, Family, !!sym(trait_col)), by = c("P1_ID" = "ID")) %>%
    dplyr::rename(P1_Family = Family, P1_Value = !!sym(trait_col)) %>%
    dplyr::left_join(strong_pool %>% dplyr::select(ID, Family, !!sym(trait_col)), by = c("P2_ID" = "ID")) %>%
    dplyr::rename(P2_Family = Family, P2_Value = !!sym(trait_col)) %>%
    # [Core biological constraint]
    # Strictly require the two strong parents to originate from different families
    dplyr::filter(P1_Family != P2_Family) %>%
    dplyr::mutate(
      Pheno_Diff = abs(P1_Value - P2_Value),
      Cross_Type = "Strong x Strong",
      Target_Trait = trait_col,
      # For strong × strong combinations, prioritize pairs with the largest combined phenotypic value
      Combined_Strength = P1_Value + P2_Value 
    ) %>%
    dplyr::arrange(desc(Combined_Strength)) %>%
    dplyr::select(-Combined_Strength) %>%
    head(top_n_pairs)
  
  # Merge the results from both strategies
  result <- dplyr::bind_rows(sxw_pairs, sxs_pairs)
  return(result)
}

# 3. Generate Recommended Crossing Lists ----
cat("Generating recommended crossing combinations for Plant Height based on real data...\n")
crosses_height <- generate_crosses(df_200, trait_col = "Plant_Height", pool_size = 15, top_n_pairs = 5)

cat("Generating recommended crossing combinations for Multifoliate Score (MF_Total) based on real data...\n")
crosses_mf <- generate_crosses(df_200, trait_col = "MF_Total", pool_size = 15, top_n_pairs = 5)

# 4. Merge and Save the Final Crossing Plan ----
final_crossing_plan <- dplyr::bind_rows(crosses_height, crosses_mf)

# Reorder the columns for improved readability
final_crossing_plan <- final_crossing_plan %>%
  dplyr::select(Target_Trait, Cross_Type, 
                P1_ID, P1_Family, P1_Value, 
                P2_ID, P2_Family, P2_Value, Pheno_Diff)

# Print a preview to the console
print(final_crossing_plan)

# Save using a distinct file name to avoid overwriting previous simulated results
write_csv(final_crossing_plan, "03_RealData_Recommended_Crossing_Plan.csv")

cat("\n[SUCCESS] The crossing plan containing your real field IDs has been saved to '03_RealData_Recommended_Crossing_Plan.csv'\n")
