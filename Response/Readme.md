# HDL-L Response Simulations

This repository contains R scripts used in our response to de Leeuw *et al.* regarding the modeling assumptions underlying LAVA [(Werme *et al.*, 2022)] and HDL-L [(Li *et al.*, 2025)]. Our work addresses key statistical and biological differences between the fixed-effects model used by LAVA and the random-effects model underlying HDL-L, demonstrating why the latter is both statistically more appropriate and biologically more interpretable for modeling genetic effects.

## Contents

- **`Figure1R.R`**  
  R code to reproduce **Figure 1R** in our response.  
  This script simulates the discrepancy between *forced* subject-level genetic correlations and the corresponding *realized* variant-level genetic correlations under varying proportions of causal variants.

- **`Figure2R.R`**  
  R code to reproduce **Figure 2R** in our response.  
  This script demonstrates how LAVA’s type-I error inflation under the random-effects model arises as a statistical artifact, particularly when the proportion of causal variants is small.

## Data Availability

The simulations use individual-level genotype data from the UK Biobank (first region of chromosome 22), available via application at:  
<https://www.ukbiobank.ac.uk>

## Usage

1. Ensure you have access to the required genotype dataset from UK Biobank.  
2. Place the genotype data files in your working directory (adjust file paths in scripts as necessary).  
3. Run the scripts in R:

```r
# For Figure 1R
source("Figure1R.R")

# For Figure 2R
source("Figure2R.R")
