# 1: Data Simulation and Generation of Summary Statistics

# Load required packages
library(dplyr)
library(mvtnorm)
library(data.table)

# Set seed for reproducibility
set.seed(2024)

# Define file paths
base_path <- "/Your/Path"
bhat_path <- file.path(base_path, "sim_bhat")
bedp_path <- file.path(base_path, "bfile")
X0_path <- file.path(base_path, "X0loc")
y_path <- file.path(base_path, "sim_y")
sumstat_path <- file.path(base_path, "sumstat")

# Create directories if they don't exist
dir.create(bhat_path, recursive = TRUE, showWarnings = FALSE)
dir.create(y_path, recursive = TRUE, showWarnings = FALSE)
dir.create(sumstat_path, recursive = TRUE, showWarnings = FALSE)
dir.create(X0_path, recursive = TRUE, showWarnings = FALSE)
dir.create(bedp_path, recursive = TRUE, showWarnings = FALSE)


# Load required data
load(file.path(base_path, "HDLL_LOC_snps.RData"))  # Loads 'NEWLOC' and 'nsnps.list.imputed'
load(file.path(base_path, "causal.snps.10.percent.list.sim.imputed.loc.RData"))  # Loads 'causal.snps.list'
NEWLOC = NEWLOC %>% filter(CHR == 22)

# Define simulation parameters
h1 <- 0.2
h2 <- 0.4
rg <- 0.5
num_simulations <- 100
N <- 335272
N1 <- N2 <- N / 2 #simulated two independent population 
index <- sample(c(rep(TRUE, ceiling(N1)), rep(FALSE, floor(N2))))
M <- length(causal.snps.list)  # Number of causal variants

# Function to generate covariance matrix for effect sizes
generate_cov_matrix <- function(h1, h2, rg, M) {
  h12 <- rg * sqrt(h1 * h2)
  genMat <- (1 / M) * matrix(c(h1, h12, h12, h2), ncol = 2)
  sigma_e1 <- 1 - h1
  sigma_e2 <- 1 - h2
  list(genMat = genMat, sigma_e1 = sigma_e1, sigma_e2 = sigma_e2)
}

# Generate covariance matrices
cov_matrices <- generate_cov_matrix(h1, h2, rg, M)
genMat <- cov_matrices$genMat
sigma_e1 <- cov_matrices$sigma_e1
sigma_e2 <- cov_matrices$sigma_e2

# Simulate true effect sizes (b)
simulate_effect_sizes <- function(M, genMat, num_simulations) {
  b_list <- vector("list", num_simulations)
  for (sim in 1:num_simulations) {
    b <- rmvnorm(M, sigma = genMat)
    rownames(b) <- causal.snps.list
    b_list[[sim]] <- b
  }
  return(b_list)
}

# Run simulations for effect sizes
b_list <- simulate_effect_sizes(M, genMat, num_simulations)

# Save the list of simulations
dir_path <- file.path(bhat_path, sprintf("h1_%s_h2_%s_rg_%s", h1, h2, rg))
dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
saveRDS(b_list, file = file.path(dir_path, "b_real_simulations.rds"))

# Initialize y1 and y2 matrices
y1_final <- matrix(0, nrow = N1, ncol = num_simulations)
y2_final <- matrix(0, nrow = N2, ncol = num_simulations)

# Function to process genotype data and accumulate genetic contributions
process_genotype_data <- function(chr, index, b_list, y1_final, y2_final) {
  num_pieces <- length(nsnps.list.imputed[[chr]])
  for (piece in 1:num_pieces) {
    X0_file <- file.path(X0_path, sprintf("ukb_chr%d.%d_X0_scaled.RData", chr, piece))
    if (!file.exists(X0_file)) next
    load(X0_file)  # Loads 'X0'
    causal_snps_in_X0 <- intersect(colnames(X0), causal.snps.list)
    if (length(causal_snps_in_X0) == 0) {
      rm(X0)
      next
    }
    X0_1_causal <- X0[index, causal_snps_in_X0, drop = FALSE]
    X0_2_causal <- X0[!index, causal_snps_in_X0, drop = FALSE]
    for (sim in 1:num_simulations) {
      b_causal <- b_list[[sim]][causal_snps_in_X0, ]
      y1_final[, sim] <- y1_final[, sim] + X0_1_causal %*% b_causal[, 1]
      y2_final[, sim] <- y2_final[, sim] + X0_2_causal %*% b_causal[, 2]
    }
    rm(X0, X0_1_causal, X0_2_causal)
    gc()
  }
  return(list(y1_final = y1_final, y2_final = y2_final))
}

# Accumulate genetic contributions across all chromosomes
#for (chr in 1:22) {
chr = 22 #adjust it when you need
  cat("Processing chromosome", chr, "\n")
  results <- process_genotype_data(chr, index, b_list, y1_final, y2_final)
  y1_final <- results$y1_final
  y2_final <- results$y2_final
#}

# Add residual errors to phenotypes
for (sim in 1:num_simulations) {
  y1_final[, sim] <- y1_final[, sim] + rnorm(N1, mean = 0, sd = sqrt(sigma_e1))
  y2_final[, sim] <- y2_final[, sim] + rnorm(N2, mean = 0, sd = sqrt(sigma_e2))
}

# Save final phenotypes
y_path2 <- file.path(y_path, sprintf("h1_%s_h2_%s_rg_%s", h1, h2, rg))
dir.create(y_path2, recursive = TRUE, showWarnings = FALSE)
save(y1_final, y2_final, file = file.path(y_path2, "y_final_y1_y2.rda"))

# Compute regression coefficients (bhat1 and bhat2)
compute_bhat <- function(chr, index, y1_final, y2_final, N1, N2) {
  num_pieces <- length(nsnps.list.imputed[[chr]])
  for (piece in 1:num_pieces) {
    X0_file <- file.path(X0_path, sprintf("ukb_chr%d.%d_X0_scaled.RData", chr, piece))
    if (!file.exists(X0_file)) next
    load(X0_file)  # Loads 'X0'
    X0_1 <- X0[index, , drop = FALSE]
    X0_2 <- X0[!index, , drop = FALSE]
    bhat1 <- crossprod(X0_1, y1_final) / N1
    bhat2 <- crossprod(X0_2, y2_final) / N2
    save(bhat1, bhat2, file = file.path(dir_path, sprintf("b.hat.chr%d.piece%d.rda", chr, piece)))
    rm(X0, X0_1, X0_2, bhat1, bhat2)
    gc()
  }
}

# Calculate bhat1 and bhat2 for all chromosomes
#for (chr in 1:22) {
  chr = 22 #adjust it when you need
  cat("Calculating bhat for chromosome", chr, "\n")
  compute_bhat(chr, index, y1_final, y2_final, N1, N2)
#}

# Combine bhat1 and bhat2 across all chromosomes
combine_bhat <- function() {
  bhat1_list <- list()
  bhat2_list <- list()
  snp_info_list <- list()
  #for (chr in 1:22) {
    chr = 22 #adjust it when you need
    num_pieces <- length(nsnps.list.imputed[[chr]])
    for (piece in 1:num_pieces) {
      bhat_file <- file.path(dir_path, sprintf("b.hat.chr%d.piece%d.rda", chr, piece))
      if (!file.exists(bhat_file)) next
      load(bhat_file)  # Loads 'bhat1' and 'bhat2'
      bhat1_list[[length(bhat1_list) + 1]] <- bhat1
      bhat2_list[[length(bhat2_list) + 1]] <- bhat2
      # Combine SNP info
      bim_file <- file.path(bedp_path, sprintf("ukb_chr%d.%d_n336000.imputed_clean.bim", chr, piece))
      if (!file.exists(bim_file)) next
      snps_info <- fread(bim_file, header = FALSE, data.table = FALSE)
      snps_selected <- snps_info[, c(2, 5, 6)]
      colnames(snps_selected) <- c("SNP", "A1", "A2")
      snp_info_list[[length(snp_info_list) + 1]] <- snps_selected
      rm(bhat1, bhat2, snps_info, snps_selected)
      gc()
    #}
  }
  bhat1_all <- do.call(rbind, bhat1_list)
  bhat2_all <- do.call(rbind, bhat2_list)
  munge_df1 <- rbindlist(snp_info_list)
  munge_df2 <- copy(munge_df1)
  # Validate dimensions
  num_snps <- nrow(munge_df1)
  if (nrow(bhat1_all) != num_snps || nrow(bhat2_all) != num_snps) {
    stop("Number of SNPs in SNP info does not match number of SNPs in bhat1_all or bhat2_all.")
  }
  cat("Number of SNPs:", num_snps, "\n")
  # Save combined data
  save(bhat1_all, bhat2_all, file = file.path(dir_path, "b.hat_all.rda"))
  save(munge_df1, munge_df2, file = file.path(sumstat_path, "mungedf12.rda"))
}

# Combine bhat values and SNP info
combine_bhat()

# Generate summary statistics
generate_sumstats <- function(bhat_all, N, sim_count, phenotype_name, munge_df) {
  sqrt_N <- sqrt(N)
  for (sim in 1:sim_count) {
    cat("Generating sumstats for", phenotype_name, "simulation", sim, "\n")
    Z <- bhat_all[, sim] * sqrt_N
    P <- pchisq(Z^2, df = 1, lower.tail = FALSE)
    sumstats <- data.frame(
      SNP = munge_df$SNP,
      A1 = munge_df$A1,
      A2 = munge_df$A2,
      N = N,
      Z = Z,
      P = P
    )
    sumstats_file <- file.path(sumstat_path, sprintf("%s.sim.%d.sum.stats.txt", phenotype_name, sim))
    fwrite(sumstats, file = sumstats_file, sep = "\t", row.names = FALSE, quote = FALSE)
  }
}

# Load combined bhat and SNP info data
load(file.path(dir_path, "b.hat_all.rda"))
load(file.path(sumstat_path, "mungedf12.rda"))

# Generate sumstats for y1 and y2
generate_sumstats(bhat1_all, N1, num_simulations, "y1", munge_df1)
generate_sumstats(bhat2_all, N2, num_simulations, "y2", munge_df2)
