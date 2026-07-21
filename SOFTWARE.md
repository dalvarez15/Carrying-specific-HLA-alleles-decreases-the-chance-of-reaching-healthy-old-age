# Software requirements

This document covers all software needed to run the HLA imputation pipeline
(`hla_imputation/`) and the main association analysis (`main_script.R`).

---

## R environment

The imputation pipeline was run using the conda environment defined in
`environment.yml`. It includes R 4.3 and all required R packages.

### Recreating the environment (internet-connected machine)

```bash
conda env create -f environment.yml
conda activate r-hibag
```

### Offline transfer with conda-pack

If the target compute server has no internet access, the environment can be
packed on an internet-connected machine and transferred as a single archive.
This is how `r-hibag.tar.gz` (referred to as `r-hibag.gz`) was created:

```bash
# On the internet-connected machine:
conda install conda-pack          # install conda-pack if not already present
conda env create -f environment.yml
conda pack -n r-hibag -o r-hibag.tar.gz

# Transfer r-hibag.tar.gz to the offline server, then:
mkdir -p /path/to/conda/envs/r-hibag
tar -xzf r-hibag.tar.gz -C /path/to/conda/envs/r-hibag
source /path/to/conda/envs/r-hibag/bin/activate
conda-unpack
```

After unpacking, activate normally with `conda activate r-hibag`.

Note: the packed environment is OS-specific. An archive built on Linux x86-64
cannot be used on macOS or Windows.

---

## PLINK

Two versions of PLINK are used at different pipeline steps.

### PLINK 1.9 (`plink`)

Used in Steps 2 and 4 for `--extract`, `--update-map`, `--update-cm`,
and `--r2 square`.

PLINK 1.9 is required for `--r2 square` (not available in PLINK2) and for
`--update-cm`/`--update-map` (removed from PLINK2).

Download: [https://www.cog-genomics.org/plink/1.9/](https://www.cog-genomics.org/plink/1.9/)

The pipeline scripts call `plink` and expect it to be on `$PATH`. If it is
installed elsewhere, update the `system()` calls in the relevant scripts.

### PLINK 2 (`plink2`)

Used in Steps 1 and 4 for `--make-bed` (pgen→bed conversion) and
`--import-dosage` (for LD matrix computation).

Download: [https://www.cog-genomics.org/plink/2.0/](https://www.cog-genomics.org/plink/2.0/)

Similarly, `plink2` is expected on `$PATH`.

### Versions used in this study

| Tool | Version | Build |
|---|---|---|
| PLINK 1.9 | v1.90b7.2 | 64-bit Linux |
| PLINK 2 | v2.00a5.12 | 64-bit Linux (AVX2) |

---

## HIBAG

HIBAG is an R/Bioconductor package for HLA allele imputation. It is included
in the conda environment (`bioconductor-hibag`).

- Publication: Zheng et al. (2014) *Pharmacogenomics Journal*
  [doi:10.1038/tpj.2013.18](https://doi.org/10.1038/tpj.2013.18)
- Source: [https://github.com/zhengxwen/HIBAG](https://github.com/zhengxwen/HIBAG)
- Bioconductor page: [https://bioconductor.org/packages/HIBAG/](https://bioconductor.org/packages/HIBAG/)

HIBAG version used in this study: **2.0.1** (Bioconductor 3.18).

### Pre-fitted model

The imputation model is not included in this repository. Download the
ancestry- and platform-matched model from:

[https://hibag.s3.amazonaws.com/hlares_index.html](https://hibag.s3.amazonaws.com/hlares_index.html)

Model used in this study: `InfiniumGlobal-European-HLA4-hg19.RData`
(European ancestry, Illumina Infinium Global Screening Array v2.0, hg19).

---

## LiftOver chain file

The hg38→hg19 chain file used in Step 2 is not included in this repository.
Download from UCSC:

[https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz)

Decompress before use: `gunzip hg38ToHg19.over.chain.gz`

---

## Main analysis (`main_script.R`)

`main_script.R` uses only standard CRAN/Bioconductor R packages and does not
require PLINK or HIBAG. Install with:

```r
install.packages(c("readxl", "tidyr", "dplyr", "ggplot2", "corrplot", "data.table"))
```
