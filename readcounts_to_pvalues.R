# Functions to compute z-scores from read counts assuming a log-normal distribution.
# Also calculates p-values and FDR-adjusted values using the Benjamini–Hochberg method.

library(tidyverse)
library(readxl)
library(openxlsx)


## Calculating reads per million (normalization): Normalizes raw read counts to 
## reads per million (RPM) for each condition. This enables comparison across 
## samples with different sequencing depths.

calculate_rpm <- function(df, conditions){
  for(cond in names(conditions)){
    cols <- conditions[[cond]]
    df <- df |> mutate(across(all_of(cols), ~ .x / sum(.x) * 1e6, .names = "{.col}_rpm"))
  }
  return(df)
}

## Averaging replicates per condition: Averages RPM values across replicates 
## within each condition. Produces a single summary value per condition per guide.

average_replicates <- function(df, conditions){
  for(cond in names(conditions)){
    cols <- paste0(conditions[[cond]], "_rpm")
    df[[paste0(cond, "_avg")]] <- rowMeans(df[, cols])
  }
  return(df)
}

##  Log2 transformation of averaged read counts:Applies log2 transformation to 
## averaged values with a pseudocount of 1. Stabilizes variance and reduces 
## skewness in the data. 

log2_transform <- function(df){
  avg_cols <- grep("_avg$", colnames(df), value = TRUE)
  df |> mutate(across(all_of(avg_cols), ~ log2(.x + 1), .names = "{.col}_log"))
}

## Calculate LFC relative to reference: Computes log fold-change of each 
## condition relative to a reference condition. Captures relative 
## enrichment or depletion across conditions.

calculate_lfc <- function(df, reference_condition, conditions){
  ref_col <- paste0(reference_condition, "_avg_log")
  conds_to_calc <- setdiff(names(conditions), reference_condition)
  for(cond in conds_to_calc){
    df[[paste0(cond, "_lfc")]] <- df[[paste0(cond, "_avg_log")]] - df[[ref_col]]
  }
  df
}

## Filter outliers based on IQR:Removes extreme values based on the IQR of the 
## reference condition. Helps reduce the influence of outliers in downstream analysis.

filter_outliers <- function(df, reference_condition){
  ref_col <- paste0(reference_condition, "_avg_log")
  iqr_val <- IQR(df[[ref_col]], na.rm = TRUE)
  lower <- quantile(df[[ref_col]], 0.25, na.rm = TRUE) - 1.5 * iqr_val
  upper <- quantile(df[[ref_col]], 0.75, na.rm = TRUE) + 1.5 * iqr_val
  df |> filter(df[[ref_col]] >= lower & df[[ref_col]] <= upper)
}

## Calculate control statistics:Calculates mean and standard deviation of LFCs 
## for a control gene. These values define the null distribution for normalization. 

get_control_stats <- function(df, control_gene, lfc_cols){
  ctrl_df <- df |> filter(gene_names == control_gene)
  stats <- sapply(lfc_cols, function(col){
    c(mean = mean(ctrl_df[[col]], na.rm = TRUE), sd = sd(ctrl_df[[col]], na.rm = TRUE))
  })
  as.data.frame(t(stats))
}

# Calculate  z-scores: Converts LFC values into z-scores using control 
## statistics. Standardizes effects across conditions for comparison.
calculate_zscores <- function(df, control_stats, lfc_cols){
  for(col in lfc_cols){
    mean_val <- control_stats[col, "mean"]
    sd_val <- control_stats[col, "sd"]
    cond_name <- sub("_lfc$", "", col)
    df[[paste0(cond_name, "_zscore")]] <- (df[[col]] - mean_val)/sd_val
  }
  df
}

## add p_values and FDR_values:Computes one-sided upper-tail p-values from 
## z-scores and adjusts them using Benjamini–Hochberg correction. 
## Produces significance and FDR estimates for each condition. 

add_pvalues_and_fdr <- function(df){
  
  z_cols <- grep("_zscore$", colnames(df), value = TRUE)
  
  for(col in z_cols){
    
    condition <- sub("_zscore$", "", col)
    
    # one-sided upper tail p-value
    p_col <- paste0(condition, "_pvalue")
    df[[p_col]] <- 1 - pnorm(df[[col]])
    
    # FDR correction (Benjamini-Hochberg)
    fdr_col <- paste0(condition, "_FDR")
    df[[fdr_col]] <- p.adjust(df[[p_col]], method = "BH")
  }
  
  df
}

##  Wrapper function : Runs the full analysis pipeline from raw counts to 
## statistical outputs. Integrates normalization, transformation, scoring, 
## and significance testing into one workflow.

base_editing_pipeline <- function(raw_count_file, experiment_design){
  df <- read_excel(raw_count_file)
  df <- calculate_rpm(df, experiment_design$conditions)
  df <- average_replicates(df, experiment_design$conditions)
  df <- log2_transform(df)
  df <- calculate_lfc(df, experiment_design$reference_condition, experiment_design$conditions)
  df <- filter_outliers(df, experiment_design$reference_condition)
  
  lfc_cols <- grep("_lfc$", colnames(df), value = TRUE)
  control_stats <- get_control_stats(df, experiment_design$control_gene, lfc_cols)
  df <- calculate_zscores(df, control_stats, lfc_cols)
  df <- add_pvalues_and_fdr(df)
  
  return(df)
}




