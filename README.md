# HDL-L - An enhanced framework for local genetic correlation analysis

## Introduction 
HDL-L is the local version of high-definition likelihood (HDL) (https://github.com/zhenin/HDL). It is specifically tailored for local heritability and genetic correlation analysis using GWAS summary statistics. As a specialized tool, HDL-L addresses the limitations of global correlation estimates by focusing on localized genetic signals which can differ in direction and magnitude across genomic loci. 

HDL-L is a full likelihood-based method, ensuring precision and computational efficiency. Notably, HDL-L achieves an approximate fifty-fold increase in computational speed compared to LAVA (Local Analysis of [co]Variant Association), facilitating more rapid analyses without sacrificing accuracy. This efficiency is critical for large-scale genomic studies where traditional methods may falter due to computational demands.

HDL-L not only enhances our ability to identify shared genetic pathways but also deepens our understanding of the biological mechanisms that underpin these relationships. By integrating both local and global correlation estimates, HDL-L provides a comprehensive view of pleiotropic genetic architectures, offering insights that are essential for both theoretical research and practical applications in human genomics.

Here is the preprint version of our paper: https://www.researchsquare.com/article/rs-4568593/v1

## System Requirements
HDL-L requires only a standard computer with enough RAM to support the in-memory operations.

### Software requirements

**OS Requirements**
This package is supported for macOS and Linux. 

**R Requirements**
R Dependencies: 
```R
dplyr
data.table
```

## Installation 
HDL-L can be easily installed from GitHub using the `remotes` package. If you don't already have `remotes` installed, the following commands will manage the installation for you:
```R
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")

remotes::install_github("YuyingLi-X/HDL-L")
library(HDLL)  #Notice: The name of the R package without `-`
```
This installation process only takes a few minutes.


## Quick vignette
You will need the reference panel and linkage disequilibrium (LD) data for each region pertaining to the European ancestry population from the UK Biobank. Access this data at [Zenodo](https://doi.org/10.5281/zenodo.11001214).

For a detailed, step-by-step tutorial on how to conduct your analysis using HDL-L, please refer to our comprehensive guide available below. This demo only takes a few minutes.

```R
# Load example GWAS summary statistics for basal metabolic rate
data(gwas1.example)

# Load example GWAS summary statistics for standing height
data(gwas2.example)

# Specify the paths to the LD and bim files
LD.path <- "/Path/to/LD/LD.path/"
bim.path <- "/Path/to/bim/bim.path/"

# Run HDL-L analysis
res.HDL <- HDL.L(
  gwas1.example, gwas2.example, 
  Trait1name = "Basal metabolic rate",  
  Trait2name = "Standing height", 
  LD.path = LD.path, 
  bim.path = bim.path, 
  chr = 1, 
  piece = 3, 
  N0 = 0
)

# View the results
print(res.HDL)

```

## Tutorial
Here, we provide a step-by-step tutorial for HDL-L and a real data example at the end. Before you begin your analysis, ensure that you have the necessary resources downloaded, so you need download reference panel and LD at first: 
### Step 1: Reference panel and local region definition
As with HDL, we already prepared the pre-computed reference panel and LD for each region of the European-ancestry population. You can download it from [Zenodo](https://doi.org/10.5281/zenodo.11001214).

In the "LD.path", it includes LD files, eigenvectors, and eigen matrixes for all local regions, end by "_LDSVD.rda"

In the "bim.path", it includes bim files for local regions, which helps to clean the summary statistics data and check if there are multiallelic or duplicated SNPs

In addition to the existing HDL-L reference panel, we have developed a new LD reference panel specifically tailored for cis-protein quantitative trait loci (cis-pQTL) analysis. This new panel was constructed using data from the Olink Proteomics proximity extension assay (PEA), focusing on approximately 3,000 proteins in the Olink® Explore panel. It can be downloaded from [Zenodo](https://zenodo.org/records/14209926).

If you want to build a reference panel for your predefined loci, you can use compute_ld_eigen.R

### Step 2: The format of summary statistics
To analyze your data using HDL-L, it is crucial to format your summary statistics correctly. Below are the required columns that your input data file must include:

- `SNP`: SNP ID  
- `A1`: Effect allele  
- `A2`: Reference allele  
- `N`: Sample size  
- `Z`: Z-score  

If `Z` is not available, alternatively, you may provide:  
- `b`: Estimate of marginal effect in GWAS  
- `se`: Standard error of the estimates of marginal effects in GWAS  

If the GWAS is based on logistic regression, `b` should be the logarithm of OR (odds ratio), and `se` is the standard error of log(OR). 

The summary statistics should look like (b and se can be absent in this example since Z is available):

```R
##          SNP A1 A2      N        b       se      Z
## 1  rs3131962  G  A 205475 0.001004 0.004590 0.2187
## 2 rs12562034  A  G 205475 0.005382 0.005011 1.0740
## 3 rs11240779  A  G 205475 0.002259 0.003691 0.6119
## 4 rs57181708  G  A 205475 0.005401 0.005114 1.0562
## 5  rs4422948  G  A 205475 0.005368 0.003604 1.4893
## 6  rs4970383  A  C 205475 0.004685 0.003582 1.3080
```

You can use `HDL.L.data.wrangling.R` to do data wrangling for data from [the Neale Lab round 2 GWAS of UK Biobank](https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit?ts=5b5f17db#gid=227859291) using commands.
```bash
Rscript /Path/to/HDL/HDL.L.data.wrangling.R \
gwas.file=/Path/to/gwas/data/datafile \
LD.path=/Path/to/LD.path/ \
GWAS.type=UKB.Neale \
output.file=/Path/to/gwas/gwas1 \
log.file=/Path/to/log/gwas1
```

If the GWAS is from other sources, you need to explicitly tell `HDL.L.data.wrangling.R` how to understand the variable names in the GWAS. Other than this, the syntax is the same as that in the previous section. For example, if your GWAS looks like this:
```R
##         rsid alt ref  tstat n_complete_samples     beta       se
## 1  rs3131962   G   A 0.2187             205475 0.001004 0.004590
## 2 rs12562034   A   G 1.0740             205475 0.005382 0.005011
## 3 rs11240779   A   G 0.6119             205475 0.002259 0.003691
## 4 rs57181708   G   A 1.0562             205475 0.005401 0.005114
## 5  rs4422948   G   A 1.4893             205475 0.005368 0.003604
## 6  rs4970383   A   C 1.3080             205475 0.004685 0.003582
```

You should use the below command to run HDL.L.data.wrangling.R
```bash
Rscript /Path/to/HDL/HDL.L.data.wrangling.R \
gwas.file=/Path/to/gwas/data/datafile \
LD.path=/Path/to/LD.path/ \
SNP=rsid A1=alt A2=ref N=n_complete_samples Z=tstat \
output.file=/Path/to/gwas/gwas1 \
log.file=/Path/to/log/gwas1
```

or

```bash
Rscript /Path/to/HDL/HDL.L.data.wrangling.R \
gwas.file=/Path/to/gwas/data/datafile \
LD.path=/Path/to/LD.path/ \
SNP=rsid A1=alt A2=ref N=n_complete_samples b=beta se=se \
output.file=/Path/to/gwas/gwas1 \
log.file=/Path/to/log/gwas1
```


### Step 3: Running HDL.L on each region

You can execute HDL.L on the entire genome or specific regions to obtain results. It is optimized for parallel processing, allowing users to configure the execution to utilize multiple cores. When utilizing a single core, the typical computational time for analyzing all local regions is estimated at approximately 1.5 hours. Users are encouraged to exploit the parallel processing capabilities of their systems by allocating additional cores, thereby reducing the overall computation time. Detailed procedural guidance is provided in the example below.

The HDL.L function takes several arguments as outlined below:

1. **gwas1.df**
The first formatted summary statistics data. The input data frame should include the following columns: 
- `SNP`: SNP ID
- `A1`: Effect allele
- `A2`: Reference allele
- `N`: Sample size
- `Z`: Z-score
Alternatively, if `Z` is not provided, you may include:
- `b`: Estimate of marginal effect in GWAS
- `se`: Standard error of the estimates of marginal effects in GWAS.

2. **gwas2.df**
For the second formatted summary statistics data, the columns should be the same as `gwas1.df`.

3. **Trait1name**
The trait name for **gwas1.df**.

4. **Trait2name**
The trait name for **gwas2.df**.

5. **LD.path**
Path to the directory where the decompressed LD.path.zip file is stored.

6. **bim.path**
Path to the directory where the decompressed bim.path.zip file is stored.

7. **Nref**
Sample size of the reference sample where LD is computed. If the default UK Biobank reference sample is used, Nref = 335272.

8. **N0**
Number of individuals included in both cohorts.

9. **output.file**
Location where the log and results should be written. If you do not specify a file, the log will be printed on the console.

10. **eigen.cut**
Specifies which eigenvalues and eigenvectors in each LD score matrix should be used for HDL. The default value is 0.95. Users are allowed to specify a numeric value between 0 and 1 for eigen.cut.

11. **intercept.output**
Logical, `FALSE` by default. Determines whether the intercept terms are included in the result `estimates.df` or not.

12. **fill.missing.N**
If `NULL` (default), SNPs with missing `N` are removed. You can specify "median", "min", or "max" so that the missing `N` will be filled accordingly. For example, "median" means the missing `N` are filled with the median `N` of the SNPs with available `N`.

13. **lim**
Tolerance limitation, default `lim = exp(-18)`.

14. **chr**
The chromosome to which the region belongs.

15. **piece**
The piece of the genome to which the region belongs. The whole genome is divided into 2,476 smaller, semi-independent blocks, each defined by LD calculated by Plink. The SNP information in each local region is included in these two data sets: "UKB_snp_counter_imputed.RData" and "UKB_snp_list_imputed_vector.RData".

## Real data example
We provide a real data example from UKBB:
```bash
cd /Path/to/Your/directory/
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/21001_irnt.gwas.imputed_v3.female.tsv.bgz -O 21001_irnt.gwas.imputed_v3.female.tsv.bgz

wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/I9_CHD.gwas.imputed_v3.male.tsv.bgz -O I9_CHD.gwas.imputed_v3.male.tsv.bgz


Rscript /Path/to/HDL/HDL.L.data.wrangling.R \
gwas.file=/Path/to/gwas/21001_irnt.gwas.imputed_v3.female.tsv.bgz \
LD.path=/Path/to/LD.path/ \
GWAS.type=UKB.Neale \
output.file=/Path/to/gwas/gwas1 \
log.file=/Path/to/log/gwas1

Rscript /Path/to/HDL/HDL.L.data.wrangling.R \
gwas.file=/Path/to/gwas/I9_CHD.gwas.imputed_v3.male.tsv.bgz \
LD.path=/Path/to/LD.path/ \
GWAS.type=UKB.Neale \
output.file=/Path/to/gwas/gwas2 \
log.file=/Path/to/log/gwas2

Rscript /Path/to/HDL/HDL.L.run.R \
gwas1.df=/Path/to/gwas/gwas1.hdl.rds \
gwas2.df=/Path/to/gwas/gwas2.hdl.rds \
Trait1name="gwas1" \
Trait2name="gwas2" \
LD.path=/Path/to/LD.path/ \
bim.path=/Path/to/bim.path/ \
N0=0 \
output.file=/Path/to/output/test.raw.gwas.Rout \
type="WG" \
cores=1 \
save.path=/Path/to/save/result/
```

If you only want to run on a specific region, you can run this command

```bash
Rscript /Path/to/HDL/HDL.L.run.R \
gwas1.df=/Path/to/gwas/gwas1.hdl.rds \
gwas2.df=/Path/to/gwas/gwas2.hdl.rds \
Trait1name="gwas1" \
Trait2name="gwas2" \
LD.path=/Path/to/LD.path/ \
bim.path=/Path/to/bim.path/ \
N0=0 \
output.file=/Path/to/output/test.raw.gwas.Rout \
chr=1 \
piece=147 \
cores=1 \
save.path=/Path/to/save/result/
```

The results of HDL.L are saved in the save.path as a rda file and here is the output file:

```R
Attaching package: ‘data.table’

The following objects are masked from ‘package:dplyr’:

    between, first, last

Loading GWAS1 ... 
Loading GWAS2 ... 
Processing chromosome 1 region 147 
Analysis starts on Mon Apr 22 00:58:23 2024 
0 SNPs were removed in GWAS 1 due to missing N or missing test statistic.  
0 SNPs were removed in GWAS 2 due to missing N or missing test statistic.  
329 out of 329 (100%) SNPs in reference panel are available in GWAS 1.  
329 out of 329 (100%) SNPs in reference panel are available in GWAS 2.  

Point estimates: 
Heritability of phenotype 1:  8.25e-05 
Heritability of phenotype 2:  2e-04 
Genetic Covariance:  7.40e-05 
Genetic Correlation:  0.6396 


Heritability of phenotype 1:  8.25e-05 ,P:0.0893 
Heritability of phenotype 2:  2e-04 , P:0.00e+00 
Genetic Correlation:  0.6396 (-0.4085,1) 
P:  0.2480275 

Analysis finished at Mon Apr 22 00:58:25 2024 
Saving results ... 
Finished!
```

## Citation
If you use the HDL-L software, please cite:

- Li, Y., Pawitan, Y., & Shen, X. *An enhanced framework for local genetic correlation analysis*. (2024)
- Ning, Z., Pawitan, Y. & Shen, X. *High-definition likelihood inference of genetic correlations across human complex traits*. Nat Genet (2020).

## License
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.

## Contact
- Report bugs by opening a new issue on this GitHub page.
- Send email to the authors: [Yuying Li](mailto:yuying.li@ki.se) or [Xia Shen](mailto:shenx@fudan.edu.cn).

## Future Plans
This repository will be integrated into the HDL software in the future.
