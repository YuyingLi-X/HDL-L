# 3.Run LAVA
# Load required packages
library(dplyr)
library(LAVA)
library(parallel)

# Define file paths
base_path <- "/Your/Path"
sumstat_path <- file.path(base_path, "sumstat")
lava_path <- file.path(base_path, "LAVA")
ldsc_path <- file.path(base_path, "ldsc")

# Create LAVA result directory if it doesn't exist
dir.create(lava_path, recursive = TRUE, showWarnings = FALSE)

# Load required data
load(file.path(base_path, "HDLL_LOC_snps.RData"))  # Loads 'NEWLOC'
NEWLOC <- NEWLOC %>% filter(CHR == 22)

# Define LAVA analysis functions
run_lava_analysis <- function(sim) {
  # Prepare info file
  phenos <- sprintf("y%ssim%d", 1:2, sim)
  infofile_save <- data.frame(
    phenotype = phenos,
    cases = NA,
    controls = NA,
    filename = file.path(sumstat_path, sprintf("y%d.sim.%d.sum.stats.txt", 1:2, sim)),
    stringsAsFactors = FALSE
  )
  infofile_name <- sprintf("input.info.sim%s.txt", sim)
  infofile_path <- file.path(sumstat_path, infofile_name)
  write.table(infofile_save, file = infofile_path, sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Define paths for munge_sumstats.py and ldsc.py
  w_hm3_snplist <- file.path(ldsc_path, "eur_w_ld_chr", "w_hm3.snplist")
  munge_sumstats_py <- file.path(ldsc_path, "munge_sumstats.py")
  ldsc_py <- file.path(ldsc_path, "ldsc.py")
  
  # Munge sumstats for each phenotype
  for (phen in infofile_save$phenotype) {
    sumstats_file <- infofile_save$filename[infofile_save$phenotype == phen]
    output_prefix <- file.path(sumstat_path, phen)
    cmd <- sprintf(
      "python2 %s --sumstats %s --out %s --merge-alleles %s --ignore b,se --chunksize 500000",
      munge_sumstats_py,
      shQuote(sumstats_file),
      shQuote(output_prefix),
      shQuote(w_hm3_snplist)
    )
    system(cmd)
  }
  
  # Calculate genetic correlation using ldsc.py
  rg_output_prefix <- file.path(sumstat_path, sprintf("y1_y2_sim%s_rg", sim))
  ldsc_cmd <- sprintf(
    "python2 %s --rg %s,%s --ref-ld-chr %s --w-ld-chr %s --out %s",
    ldsc_py,
    shQuote(file.path(sumstat_path, sprintf("y1sim%s.sumstats.gz", sim))),
    shQuote(file.path(sumstat_path, sprintf("y2sim%s.sumstats.gz", sim))),
    shQuote(file.path(ldsc_path, "eur_w_ld_chr/")),
    shQuote(file.path(ldsc_path, "eur_w_ld_chr/")),
    shQuote(rg_output_prefix)
  )
  system(ldsc_cmd)
  
  # Process RG scores and create overlap matrix
  rg_log_file <- paste0(rg_output_prefix, ".log")
  if (file.exists(rg_log_file)) {
    rg_data <- readLines(rg_log_file)
    indices <- grep("gcov_int", rg_data)
    next_indices <- indices + 1
    relevant_lines <- c(rg_data[indices], rg_data[next_indices])
    relevant_order <- sort(c(indices, next_indices))
    relevant_lines_sorted <- rg_data[relevant_order]
    all_rg_file <- file.path(sumstat_path, sprintf("all.rg.sim%s.txt", sim))
    writeLines(relevant_lines_sorted, con = all_rg_file)
  } else {
    warning("rg log file not found. Skipping rg data aggregation.")
  }
  
  # Read and process RG scores
  rg_scores_file <- file.path(sumstat_path, sprintf("all.rg.sim%s.txt", sim))
  if (file.exists(rg_scores_file)) {
    scor <- read.table(rg_scores_file, header = TRUE, stringsAsFactors = FALSE)
    required_cols <- c("p1", "p2", "gcov_int")
    if (!all(required_cols %in% colnames(scor))) {
      stop("rg scores file does not contain the required columns.")
    }
    scor <- scor %>%
      mutate(
        p1 = basename(p1),
        p2 = basename(p2),
        p1 = gsub(".sumstats.gz", "", p1),
        p2 = gsub(".sumstats.gz", "", p2)
      )
    phen <- unique(c(scor$p1, scor$p2))
    n <- length(phen)
    mat <- diag(1, nrow = n, ncol = n)
    rownames(mat) <- colnames(mat) <- phen
    for (i in seq_len(nrow(scor))) {
      rna <- scor$p1[i]
      cna <- scor$p2[i]
      mat[rna, cna] <- scor$gcov_int[i]
      mat[cna, rna] <- scor$gcov_int[i]
    }
    mat <- round(cov2cor(mat), 5)
    overlap_file <- file.path(sumstat_path, sprintf("sample.overlap.sim%s.txt", sim))
    write.table(mat, file = overlap_file, sep = "\t", row.names = TRUE, quote = FALSE)
  } else {
    stop("rg scores file not found. Cannot proceed with overlap matrix.")
  }
  
  # Prepare LAVA input
  phenos <- sprintf("y%ssim%s", 1:2, sim)
  refdat <- file.path(ldsc_path, "g1000_eur", "g1000_eur")
  n_loc <- nrow(NEWLOC)
  input <- process.input(
    input.info.file = infofile_path,
    sample.overlap.file = overlap_file,
    ref.prefix = refdat,
    phenos = phenos
  )
  
  # Define output file names
  outname_univ <- file.path(lava_path, sprintf("results.univ.sim%s.txt", sim))
  outname_bivar <- file.path(lava_path, sprintf("results.bivar.sim%s.txt", sim))
  
  # Run LAVA analysis
  message(sprintf("Starting LAVA analysis for %s loci", n_loc))
  progress_points <- ceiling(seq(0, n_loc, length.out = 20))
  num_cores <- 1  # Adjust as needed
  out <- mclapply(
    1:n_loc,
    mc.cores = num_cores,
    function(i) {
      if (i %in% progress_points) {
        message(sprintf("Progress: %.1f%% (%d/%d)", (i / n_loc) * 100, i, n_loc))
      }
      locus <- process.locus(NEWLOC[i, ], input)
      if (is.null(locus)) return(NULL)
      loc_info <- data.frame(
        locus = locus$id,
        chr = locus$chr,
        start = locus$start,
        stop = locus$stop,
        n_snps = locus$n.snps,
        n_pcs = locus$K,
        stringsAsFactors = FALSE
      )
      univ <- tryCatch(run.univ(locus), error = function(e) NULL)
      bivar <- tryCatch(run.bivar(locus), error = function(e) NULL)
      list(
        univ = if (!is.null(univ)) cbind(loc_info, univ),
        bivar = if (!is.null(bivar)) cbind(loc_info, bivar)
      )
    }
  )
  univ_data <- do.call(rbind, lapply(out, `[[`, "univ"))
  bivar_data <- do.call(rbind, lapply(out, `[[`, "bivar"))
  if (!is.null(univ_data)) {
    write.table(univ_data, file = outname_univ, sep = "\t", row.names = FALSE, quote = FALSE)
  }
  if (!is.null(bivar_data)) {
    write.table(bivar_data, file = outname_bivar, sep = "\t", row.names = FALSE, quote = FALSE)
  }
}

# Run LAVA for each simulation
num_simulations <- 100  # Adjust if needed
for (sim in 1:num_simulations) {
  cat("Running LAVA for simulation", sim, "\n")
  run_lava_analysis(sim)
}