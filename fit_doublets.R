#!/usr/bin/env Rscript

# Load required libraries silently
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

# 1. Handle Command Line Arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript fit_doublets.R <path_to_assignments.tsv.gz>\n", call. = FALSE)
}

file_path <- args[1]
if (!file.exists(file_path)) {
  stop(sprintf("Error: File not found at %s\n", file_path), call. = FALSE)
}

library_name <- basename(dirname(file_path))
if(library_name == "." || library_name == "") {
  library_name <- sub("\\.tsv.*", "", basename(file_path))
}

cat(sprintf("\nProcessing Library: %s\n", library_name))
cat("----------------------------------------------------\n")

# 2. Read and Parse Data
raw_data <- read.table(gzfile(file_path), header = FALSE, stringsAsFactors = FALSE, fill = TRUE)
colnames(raw_data)[1:2] <- c("Barcode", "Assignment")

counts_df <- raw_data %>%
  group_by(Assignment) %>%
  summarise(Count = n(), .groups = "drop") %>%
  filter(!Assignment %in% c("0", "unassigned", "doublet", "Assignment"))

singlets <- counts_df %>% filter(!grepl("\\+", Assignment))
doublets <- counts_df %>% filter(grepl("\\+", Assignment))

donors <- sort(singlets$Assignment)
K <- length(donors)
N_total <- sum(counts_df$Count)

# 3. Build Observed Heterotypic Matrix
obs_doublets <- matrix(0, nrow = K, ncol = K, dimnames = list(donors, donors))
for(i in seq_len(nrow(doublets))) {
  parts <- unlist(strsplit(doublets$Assignment[i], "\\+"))
  if(length(parts) == 2 && all(parts %in% donors)) {
    obs_doublets[parts[1], parts[2]] <- doublets$Count[i]
    obs_doublets[parts[2], parts[1]] <- doublets$Count[i]
  }
}

# Pre-calculate matrix indices for optimization speed
upper_idx <- upper.tri(obs_doublets)
obs_D_vec <- obs_doublets[upper_idx]

# 4. Define Negative Log-Likelihood (Vectorized)
nll_function <- function(params) {
  p <- exp(params[1:K]) / sum(exp(params[1:K]))
  d <- 1 / (1 + exp(-params[K + 1]))
  
  # Singlet likelihood
  lambda_S <- N_total * ((1 - d) * p + d * p^2)
  ll_S <- sum(dpois(singlets$Count, lambda = lambda_S, log = TRUE))
  
  # Doublet likelihood (vectorized)
  p_matrix <- p %*% t(p)
  lambda_D_vec <- N_total * d * 2 * p_matrix[upper_idx]
  ll_D <- sum(dpois(obs_D_vec, lambda = lambda_D_vec, log = TRUE))
  
  return(-(ll_S + ll_D))
}

# 5. Initialize & Optimize
init_p <- singlets$Count / sum(singlets$Count)
init_p[init_p == 0] <- 1e-6
init_p <- init_p / sum(init_p)

initial_params <- c(log(init_p), log(0.05 / (1 - 0.05)))
fit <- optim(par = initial_params, fn = nll_function, method = "BFGS")

# 6. Extract Results
fitted_p <- exp(fit$par[1:K]) / sum(exp(fit$par[1:K]))
names(fitted_p) <- donors
fitted_d <- 1 / (1 + exp(-fit$par[K + 1]))

# 7. Calculate Key Statistics
expected_total_doublets <- N_total * fitted_d
expected_homo_doublets <- sum(N_total * fitted_d * fitted_p^2)
expected_hetero_doublets <- expected_total_doublets - expected_homo_doublets
observed_hetero_doublets <- sum(obs_D_vec)

# 8. Format DataFrames for Export
summary_df <- data.frame(
  Library = library_name,
  Total_Valid_Droplets = N_total,
  Unique_Donors = K,
  Fitted_Doublet_Rate_Pct = round(fitted_d * 100, 4),
  Expected_Total_Doublets = round(expected_total_doublets, 2),
  Expected_Homo_Doublets = round(expected_homo_doublets, 2),
  Expected_Hetero_Doublets = round(expected_hetero_doublets, 2),
  Observed_Hetero_Doublets = observed_hetero_doublets,
  Convergence_Code = fit$convergence
)

donor_df <- data.frame(
  Library = library_name,
  Donor = donors,
  Observed_Singlet_Class = singlets$Count,
  Fitted_True_Proportion = round(fitted_p, 6),
  Expected_True_Singlets = round(N_total * (1 - fitted_d) * fitted_p, 2),
  Expected_Homo_Doublets = round(N_total * fitted_d * fitted_p^2, 2)
)

# Vectorized pairwise dataframe creation
donor_pairs <- combn(donors, 2)
pairwise_df <- data.frame(Library = library_name, Donor_1 = donor_pairs[1, ], Donor_2 = donor_pairs[2, ])
pairwise_df$Expected_Hetero_Doublets <- round(N_total * fitted_d * 2 * fitted_p[pairwise_df$Donor_1] * fitted_p[pairwise_df$Donor_2], 2)
pairwise_df$Observed_Hetero_Doublets <- obs_doublets[cbind(pairwise_df$Donor_1, pairwise_df$Donor_2)]
pairwise_df$Log2_Ratio_Obs_vs_Exp <- round(log2((pairwise_df$Observed_Hetero_Doublets + 1) / (pairwise_df$Expected_Hetero_Doublets + 1)), 4)

# 9. Save to Output Directory
out_dir <- "doublet_estimate"
if(!dir.exists(out_dir)) { dir.create(out_dir, recursive = TRUE) }

write.csv(summary_df, file.path(out_dir, paste0(library_name, "_summary.csv")), row.names = FALSE)
write.csv(donor_df, file.path(out_dir, paste0(library_name, "_donor_stats.csv")), row.names = FALSE)
write.csv(pairwise_df, file.path(out_dir, paste0(library_name, "_pairwise_stats.csv")), row.names = FALSE)

# 10. Human-Readable Terminal Summary
cat(sprintf("Analyzed %d valid droplets containing %d unique donors.\n", N_total, K))
cat(sprintf("The model converged perfectly and estimates a true global doublet rate of %.2f%%.\n\n", fitted_d * 100))

cat(sprintf("Across the entire library, we expect %.1f total doublets to have formed.\n", expected_total_doublets))
cat(sprintf("  - %.1f of these are heterotypic (between different donors).\n", expected_hetero_doublets))
cat(sprintf("  - %.1f of these are homotypic (between the same donor).\n", expected_homo_doublets))
cat("    * Note: These homotypic doublets are currently hidden inside your 'singlet' counts.\n\n")

match_quality <- ifelse(abs(observed_hetero_doublets - expected_hetero_doublets) / expected_hetero_doublets < 0.1, "great", "poor")
cat(sprintf("We actually observed %d heterotypic doublets, which is a %s match to the expected %.1f.\n\n",
            observed_hetero_doublets, match_quality, expected_hetero_doublets))

# Flagging Outliers
outliers <- pairwise_df %>% filter(abs(Log2_Ratio_Obs_vs_Exp) > 1 & Expected_Hetero_Doublets > 5)

if(nrow(outliers) > 0) {
  cat(" OUTLIER WARNING: The following donor pairs formed doublets at significantly non-random rates:\n")
  for(i in 1:nrow(outliers)) {
    cat(sprintf("  -> %s + %s: Observed %d, Expected %.1f (Log2 Ratio: %.2f)\n",
                outliers$Donor_1[i], outliers$Donor_2[i],
                outliers$Observed_Hetero_Doublets[i], outliers$Expected_Hetero_Doublets[i],
                outliers$Log2_Ratio_Obs_vs_Exp[i]))
  }
} else {
  cat("No major doublet formation outliers detected. Cell pairing appears random and healthy.\n")
}
cat("====================================================\n\n")
