library(dplyr)
library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)
fn <- gsub(x = args[grep(x = args, pattern = "gwas.file=")], pattern = "gwas.file=", replacement = "")
LD.path <- gsub(x = args[grep(x = args, pattern = "LD.path=")], pattern = "LD.path=", replacement = "")
GWAS.type <- gsub(x = args[grep(x = args, pattern = "GWAS.type=")], pattern = "GWAS.type=", replacement = "")
output.file <- gsub(x = args[grep(x = args, pattern = "output.file=")], pattern = "output.file=", replacement = "")
log.file <- gsub(x = args[grep(x = args, pattern = "log.file=")], pattern = "log.file=", replacement = "")

SNP <- gsub(x = args[grep(x = args, pattern = "SNP=")], pattern = "SNP=", replacement = "")
A1 <- gsub(x = args[grep(x = args, pattern = "A1=")], pattern = "A1=", replacement = "")
A2 <- gsub(x = args[grep(x = args, pattern = "A2=")], pattern = "A2=", replacement = "")
N <- gsub(x = args[grep(x = args, pattern = "N=")], pattern = "N=", replacement = "")
b <- gsub(x = args[grep(x = args, pattern = "b=")], pattern = "b=", replacement = "")
se <- gsub(x = args[grep(x = args, pattern = "se=")], pattern = "se=", replacement = "")
Z <- gsub(x = args[grep(x = args, pattern = "Z=")], pattern = "Z=", replacement = "")

if(length(output.file) == 0)
  output.file <- fn

if(length(log.file) != 0){
  log.file <- paste(log.file, "txt", sep = ".")
  if(file.exists(log.file) == T){
    system(paste0("rm ",log.file))
  }
}

library(data.table)
data.table.version <- packageVersion("data.table")
if(data.table.version < "1.12.1"){
  message("Searching for the updated version of package 'data.table' in user's library:")
  detach("package:data.table", unload=TRUE)
  library(data.table, lib.loc = Sys.getenv("R_LIBS_USER"))
}


smart.reader <- function(path){
  path.split <- unlist(strsplit(path, split = "\\."))
  file.type <- path.split[length(path.split)]
  if(file.type == "rds"){
    return(readRDS(path))
  } else if(file.type == "gz" | file.type == "bgz"){
    options(datatable.fread.input.cmd.message=FALSE)
    return(fread(input = paste("zcat < ",path)))
  } else{
    try_error <- try(return(fread(path)))
    if(class(try_error) == "try-error"){
      error.message <- "This file type is not supported by fread function in data.table package. Please reformat it to .txt, .csv or .tsv."
      if(output.file != ""){
        cat(error.message, file = output.file, append = T)
      }
      stop(error.message)
    }
  }
}

time.start <- date()
cat("Program starts on",time.start,"\n")
cat("Loading GWAS summary statistics from",fn,"\n")

if(length(log.file) != 0){
  cat("Program starts on",time.start,"\n", file = log.file, append = T)
  cat("Loading GWAS summary statistics from",fn,"\n", file = log.file, append = T)
}
gwas.all <- smart.reader(fn)

cat("Data are loaded successfully. Data wrangling starts. \n")
if(length(log.file) != 0){
  cat("Data are loaded successfully. Data wrangling starts. \n", file = log.file, append = T)
}

# Load SNP list and counter files with error handling
    snp_list_file <- paste0(LD.path, "UKB_snp_list_imputed_vector.RData")
    snp_counter_file <- paste0(LD.path, "UKB_snp_counter_imputed.RData")
    # Load the snp list and counter files
    if (!is.null(snp_list_file)||!is.null(snp_counter_file)) {
        load(snp_list_file)
        load(snp_counter_file)
    } else {
        stop("No snp informations' files found in LD.path.")
}


# Check LD files exsitence
all_files <- list.files(LD.path, pattern = "\\.rda$", full.names = TRUE)

expected_patterns <- with(nsnps.df.imputed, paste0("ukb_", chr, "\\.", piece, "[\\._].*_LDSVD.rda$"))

# Function to check if any file matches the expected pattern
pattern_exists <- function(pattern) {
  sum(str_detect(all_files, regex(pattern, ignore_case = TRUE))) > 0
}

# Apply the function to each pattern
file_presence <- map_lgl(expected_patterns, pattern_exists)

# Combine results with the original DataFrame to see which files are missing
results <- nsnps.df.imputed %>%
  mutate(file_exists = file_presence, 
         expected_pattern = expected_patterns)


# To simply check if all files exist
all_files_exist <- all(file_presence)
if(!all_files_exist){
    error.message <- paste0("Lacking following LD file in LD.path:", "\n",
                     results[!results$file_exists, ], "\n",
                        "Please check your LD.path again.")
if (output.file != "") {
      cat(error.message, file = output.file, append = TRUE)
    }
    stop(error.message)
}


## the Neale's UKB GWAS format ##
if(length(GWAS.type) != 0){
  if(GWAS.type == "UKB.Neale"){
    load(file=paste0(LD.path, "snp.dictionary.imputed.rda"))
    gwas.hdl.df <- gwas.all %>%
      inner_join(snp.dictionary %>% filter(rsid %in% snps.list.imputed.vector), by = "variant")  %>%
      select(rsid, alt, ref, n_complete_samples, tstat) %>%
      rename(SNP = rsid, A1 = alt, A2 = ref, N = n_complete_samples, Z = tstat)
  }
}

## non built-in format
if(length(GWAS.type) == 0){
  if(length(Z) == 0 & (length(b) == 0 | length(se) == 0)){
    error.message <- "Z-score is not available and either b or se is missing. Please check."
    if(length(log.file) != 0){
      cat(error.message, file = log.file, append = T)
    }
    stop(error.message)
  }

  if(length(Z) != 0){
    gwas.hdl.df <- gwas.all %>%
      rename(SNP = SNP, A1 = A1, A2 = A2, N = N, Z = Z) %>%
      filter(SNP %in% snps.list.imputed.vector)
  } else{
    gwas.hdl.df <- gwas.all %>%
      rename(SNP = SNP, A1 = A1, A2 = A2, N = N, b = b, se = se) %>%
      filter(SNP %in% snps.list.imputed.vector)
    
    # b is likely to be OR in stead of log(OR)
    if(abs(median(gwas.hdl.df$b) - 1) < 0.1){
      cat("Taking log(b) because b is likely to be OR in stead of log(OR). \n")
      if(length(log.file) != 0){
        cat("Taking log(b) because b is likely to be OR in stead of log(OR). \n", file = log.file, append = T)
      }
      gwas.hdl.df <- gwas.hdl.df %>% mutate(b = log(b)) %>%
        mutate(Z = (b/se)) %>% select(SNP, A1, A2, N, Z)
    }
  }
}
gwas.hdl.df$Z[!is.finite(gwas.hdl.df$Z)] <- NA
cat("Data wrangling completed. \n")
if(length(log.file) != 0){
  cat("Data wrangling completed. \n", file = log.file, append = T)
}


gwas.hdl.df$A1 <- toupper(gwas.hdl.df$A1)
gwas.hdl.df$A2 <- toupper(gwas.hdl.df$A2)

fn.rds <- paste0(output.file, ".hdl.rds")
saveRDS(gwas.hdl.df, fn.rds)
cat("The output is saved to", fn.rds, "\n")
if(length(log.file) != 0){
  cat("The output is saved to", fn.rds, "\n", file = log.file, append = T)
}

if(length(log.file) != 0){
  cat("The log is saved to", log.file, "\n")
  cat("The log is saved to", log.file, "\n", file = log.file, append = T)
}
