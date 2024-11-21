# Simulation Pipeline for HDL-L and LAVA Analyses

This repository contains a simulation pipeline for performing HDL-L and LAVA analyses using simulated genetic data. The pipeline is divided into three scripts:

1. **Data Simulation and Summary Statistics Generation**
2. **HDL-L Analysis**
3. **LAVA Analysis**

The simulation currently focuses on chromosome 22 but can be adjusted in the code for whole-genome simulation.

## Table of Contents

- [Introduction](#introduction)
- [Requirements](#requirements)
- [Data Preparation](#data-preparation)
- [Scripts Overview](#scripts-overview)
  - [Script 1: Data Simulation and Summary Statistics Generation](#script-1-data-simulation-and-summary-statistics-generation)
  - [Script 2: HDL-L Analysis](#script-2-hdl-l-analysis)
  - [Script 3: LAVA Analysis](#script-3-lava-analysis)
- [Running the Pipeline](#running-the-pipeline)
- [Adjusting for Whole-Genome Simulation](#adjusting-for-whole-genome-simulation)
- [Outline of the Simulation Process](#outline-of-the-simulation-process)
- [References](#references)
- [License](#license)

## Introduction

This simulation pipeline is designed to generate synthetic genetic data and perform HDL-L and LAVA analyses. It allows researchers to study genetic correlations and shared genetic architectures between two traits using simulated datasets.

## Requirements

### Programming Languages and Tools

- **R** (version ≥ 4.0)
- **Python 2.7** (required for LDSC tools)

### R Packages

- `dplyr`
- `mvtnorm`
- `data.table`
- `LAVA` (install from GitHub)
- `parallel`

### External Tools

- **LDSC**: [Linkage Disequilibrium Score Regression](https://github.com/bulik/ldsc)
- **LAVA**: [Local Analysis of Variants and Annotations](https://github.com/josefin-werme/LAVA)
- **HDL-L**

### Data Files

- **Sample `X0` Data**: Two sample `X0` data files are provided in HDL-L/data directory for testing purposes.
- **LD Files**: Download required LD files following HDL-L instructions.
- **SNP Information**: Obtain SNP information files.
- **bim Files**: Ensure bim files corresponding to the genotype data are available.

## Data Preparation

Before running the scripts, make sure you have prepared the necessary data:

1. **Sample `X0` Data**:

   - Place the provided sample `X0` data files in the `X0loc` directory.

2. **LD Files**:

   - Download LD files required by HDL-L and LAVA and place them in appropriate directories as specified by each tool's instructions.

3. **SNP Information and bim Files**:

   - Download SNP information and bim files needed for the analyses and place them in the `bfile` directory.

## Scripts Overview

### Script 1: Data Simulation and Summary Statistics Generation

**Filename:** `script1_data_simulation.R`

**Description:**

This script performs the following tasks:

- **Effect Size Simulation**:
  - Generates simulated effect sizes (`b`) for two traits based on specified heritabilities (`h1`, `h2`) and genetic correlation (`rg`).

- **Phenotype Simulation**:
  - Simulates phenotypes (`y1` and `y2`) by combining genetic contributions and residual errors.

- **Regression Coefficients Calculation**:
  - Computes regression coefficients (`bhat1` and `bhat2`) from the simulated data.

- **Summary Statistics Generation**:
  - Calculates Z-scores and p-values to create summary statistics files for each trait.

- **Outputs**:
  - Saves the summary statistics to files for further analysis.

### Script 2: HDL-L Analysis

**Filename:** `script2_hdl_analysis.R`

**Description:**

This script performs HDL-L analysis using the simulated summary statistics:

- **Loads Summary Statistics**:
  - Reads the summary statistics generated in Script 1.

- **Runs HDL-L Analysis**:
  - Estimates genetic correlations and other parameters between the two traits using HDL-L.

- **Outputs**:
  - Saves the HDL-L analysis results to output files.

### Script 3: LAVA Analysis

**Filename:** `script3_lava_analysis.R`

**Description:**

This script performs LAVA analysis using the simulated summary statistics:

- **Prepares Input Files**:
  - Munges summary statistics and computes sample overlap using LDSC.

- **Runs LAVA Analysis**:
  - Performs univariate and bivariate LAVA analyses across genomic loci.

- **Outputs**:
  - Saves the LAVA analysis results to output files.

## Running the Pipeline

### Step 1: Ensure All Requirements Are Met

- **Install the necessary R packages**:

  ```R
  install.packages(c("dplyr", "mvtnorm", "data.table", "parallel"))
  install.packages("devtools")
  ```
  
- Install LDSC, LAVA and HDL-L by following their GitHub instructions:

LDSC: [Installation Guide](https://github.com/bulik/ldsc)
HDL-L: [Installation Guide](https://github.com/YuyingLi-X/HDL-L)
LAVA: [Installation Guide](https://github.com/josefin-werme/LAVA)

### Step 2: Prepare the Data
Place the sample X0 data files in the X0loc directory.
Download and place the required LD files, SNP information in LD_path, and BIM files in bedp_path directories.

### Step 3: Run Script 1
```bash
Rscript 1.Generation_of_Summary_Statistics.R
```

### Step 4: Run HDL-L and LAVA
```bash
Rscript HDL-L.R
Rscript LAVA.R
```
## Adjusting for Whole-Genome Simulation

The current simulation is set up for chromosome 22. Adjust the simulation for the whole genome when you need.

In the scripts, modify loops that iterate over chromosomes to include chromosomes 1 to 22:

```R
for (chr in 1:22) {
  # Existing code
}
```




## References
LDSC: https://github.com/bulik/ldsc
LAVA: https://github.com/josefin-werme/LAVA
