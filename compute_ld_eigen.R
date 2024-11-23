# An example of how to use this R script to split genomic bim file and calculate LD.
# snps = "the_vector_of_snps_file"
# 
# bim_file <- "path/to/your_file.bim"         # Path to your .bim file
# bed_file <- "path/to/your_file.bed"         # Path to your .bed file
# plink_path <- "plink"                       # Assuming PLINK is in your PATH
# output_prefix <- "output/ld_eigen"          # Prefix for output files
# folder_path <- "path/to/your/fortran_files" # Path to folder containing ldscore.so, bmult.so, Rbmult.r
# bandwidth <- 500                           # Bandwidth for LD calculation
# nval <- NULL                                # Number of eigenvalues/vectors to compute
# num_threads <- 4                            # Number of threads for PLINK
# 
# # Run the function
# result <- compute_ld_eigen(
#   snps = snps,
#   bim_file = bim_file,
#   bed_file = bed_file,
#   plink_path = plink_path,
#   output_prefix = output_prefix,
#   folder_path = folder_path,
#   bandwidth = bandwidth,
#   nval = nval,
#   num_threads = num_threads
# )



# Load required packages
library(RSpectra)
library(data.table)

# Define the function
compute_ld_eigen <- function(snps, bim_file, bed_file, plink_path, output_prefix, folder_path, bandwidth = 500, nval = NULL, num_threads = 1) {
  # Create a temporary directory for intermediate files
  temp_dir <- tempdir()
  
  # Write SNP list to a file
  snp_list_file <- file.path(temp_dir, "snps_to_extract.txt")
  write.table(snps, file = snp_list_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # Extract SNPs using PLINK
  extracted_bed_prefix <- file.path(temp_dir, "extracted_snps")
  plink_extract_cmd <- paste(
    plink_path,
    "--bfile", sub("\\.bim$", "", bim_file),
    "--extract", snp_list_file,
    "--make-bed",
    "--out", extracted_bed_prefix,
    "--threads", num_threads
  )
  system(plink_extract_cmd)
  
  # Compute banded LD matrix using PLINK
  ld_output_prefix <- file.path(temp_dir, "ld_matrix")
  plink_ld_cmd <- paste(
    plink_path,
    "--bfile", extracted_bed_prefix,
    "--r",
    "--ld-window", bandwidth + 1,
    "--ld-window-kb", 100000,
    "--out", ld_output_prefix,
    "--threads", num_threads
  )
  system(plink_ld_cmd)
  
  # Load required shared libraries and source files
  dyn.load(paste0(folder_path, "/fortran/ldscore.so"))
  dyn.load(paste0(folder_path, "/fortran/bmult.so"))
  source(paste0(folder_path, "/fortran/Rbmult.r"))
  
  # Read the banded LD matrix
  ld_matrix_file <- paste0(ld_output_prefix, ".ld")
  bandedR <- read.table(ld_matrix_file, header = TRUE)
  
  # Calculate nb and M
  nb <- sum(bandedR[, 2] == bandedR[1, 2])
  M <- length(unique(bandedR[, 2])) + 1
  pout <- bandedR[, "R"]
  
  # Read SNP info
  snps_info_file <- paste0(extracted_bed_prefix, ".bim")
  snps_info <- read.table(snps_info_file, header = FALSE)
  
  # Check if nval is specified, else set default
  if (is.null(nval)) {
    nval <- round(M * 0.5)
  }
  
  # Perform eigen decomposition
  eigR <- eigs_sym(Rxfun, k = nval, n = M, which = "LA", args = list(nb = nb, pout = pout))
  lam <- eigR$values
  V <- eigR$vectors
  rownames(V) <- snps_info[, 2]
  
  # Compute LD scores
  LDsc <- Rldscore(M, nb, pout)
  
  # Save the results
  save(nb, LDsc, lam, V, file = paste0(output_prefix, "_LDSVD.rda"))
  
  # Clean up temporary files
  unlink(temp_dir, recursive = TRUE)
  
  # Return the results
  return(list(
    nb = nb,
    LDsc = LDsc,
    eigenvalues = lam,
    eigenvectors = V,
    snps = snps_info[, 2]
  ))
}


