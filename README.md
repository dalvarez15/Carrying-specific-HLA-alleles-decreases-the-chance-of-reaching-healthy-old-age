# Carrying specific HLA alleles decreases the chance of reaching healthy old age

Code repository for:

> Álvarez Sirvent D, Tesi N, Hulsman M, Salazar AN, van Schoor NM, Huisman M, Pijnenburg Y, van der Flier WM, Strijbis EMM, Uitdehaag BMJ, Reinders MJT, van der Lee S, Holstege H. *Carrying specific HLA alleles decreases the chance of reaching healthy old age.* npj Aging (2025).

---

## Repository structure

```
.
├── main_script.R              # Association analysis and figure/table generation
├── README.md                  # This file
├── SOFTWARE.md                # Software versions, environment setup, PLINK/HIBAG install
├── environment.yml            # Conda environment for the imputation pipeline
├── inputs/                    # Input files for main_script.R (see below)
├── figures/                   # Output figures (Figures 2–6, Supplementary Figure 1)
├── tables/                    # Output tables (Table 1, Supplementary Data 1–2)
└── hla_imputation/            # Upstream pipeline: HLA imputation and regression
    ├── 1_make_input.R
    ├── 2_lift_over_to_hg19.R
    ├── 3_impute_hla_genotypes.R
    ├── 4_compute_regression.R
    └── hla_imputation_job.sh
```

---

## Part 1: HLA imputation pipeline (`hla_imputation/`)

This upstream pipeline imputes classical two-field HLA alleles (A, B, C, DRB1, DQB1, DQA1, DPB1) from SNP array genotype data using [HIBAG](https://github.com/zhengxwen/HIBAG), then computes the association statistics and LD matrix that feed into `main_script.R`. It must be run before the main analysis, and its outputs copied to `inputs/`.

### Steps

| Step | Script | Input | Output |
|---|---|---|---|
| 1 | `1_make_input.R` | PLINK2 chr6 fileset (hg38) | PLINK1 BED fileset |
| 2 | `2_lift_over_to_hg19.R` | PLINK1 BED fileset (hg38) | PLINK1 BED fileset (hg19) |
| 3 | `3_impute_hla_genotypes.R` | PLINK1 BED fileset (hg19) | `result_<LOCUS>.txt`, `summary_<LOCUS>.txt` per locus |
| 4 | `4_compute_regression.R` | HIBAG results + phenotypes + PCs | `regression_npj_aging_2025.csv`, `hla.ld`, `hla.bim` → copy to `inputs/` |

### Running Steps 1–3

Edit the paths in `hla_imputation_job.sh` and submit with:

```bash
sbatch hla_imputation_job.sh
```

Or run each step manually:

```bash
Rscript 1_make_input.R /path/to/chr6.dose /path/to/output/

Rscript 2_lift_over_to_hg19.R /path/to/output/ chr6_hla_imputation /path/to/hg38ToHg19.over.chain

Rscript 3_impute_hla_genotypes.R \
    /path/to/output/chr6_hla_imputation_updated \
    /path/to/InfiniumGlobal-European-HLA4-hg19.RData \
    /path/to/output/imputed_HLA_alleles/
```

### Running Step 4

Step 4 requires additional cohort-specific input files (phenotype file, PC file, allele dosage table for LD) that are not part of the imputation pipeline proper. Edit the `Configuration` block at the top of `4_compute_regression.R` and run:

```bash
Rscript 4_compute_regression.R
```

Copy the three output files to `inputs/` before running `main_script.R`:

```bash
cp /path/to/output/regression_npj_aging_2025.csv inputs/
cp /path/to/output/hla.ld inputs/
cp /path/to/output/hla.bim inputs/
```

### External files required

- **HIBAG pre-fitted model**: download the ancestry- and platform-matched `.RData` from [https://hibag.s3.amazonaws.com/hlares_index.html](https://hibag.s3.amazonaws.com/hlares_index.html). Model used in this study: `InfiniumGlobal-European-HLA4-hg19.RData` (European ancestry, Illumina Infinium Global Screening Array v2.0).
- **UCSC chain file**: `hg38ToHg19.over.chain`, from [https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/).

### Software and environment

See [`SOFTWARE.md`](SOFTWARE.md) for full details on required software versions,
R packages, PLINK installations, and the conda environment.

To recreate the R environment:

```bash
conda env create -f environment.yml
conda activate r-hibag
```

---

## Part 2: Association analysis (`main_script.R`)

Computes HLA allele associations with cognitively healthy centenarian (CHC) status, calculates the Centenarian Effect Ratio (CER) against published disease effect sizes, and produces all manuscript figures and tables.

### Running the script

Update the `setwd()` path at the top of the script to the repository root, then:

```bash
Rscript main_script.R
```

### Input files (`inputs/`)

| File | Description | Source |
|---|---|---|
| `regression_npj_aging_2025.csv` | Per-allele logistic regression results (CHC vs. controls) | This study (`4_compute_regression.R`) |
| `BUTLER-LAPORTE_2024.xlsx` | UK Biobank HLA–disease PheWAS summary statistics | Butler-Laporte et al. (2024) |
| `AD_GWAS_2022.xlsx` | Alzheimer's disease GWAS HLA summary statistics | Bellenguez et al. (2022) |
| `hla.bim` | HLA allele identifiers in PLINK BIM format | This study (`4_compute_regression.R`) |
| `hla.ld` | Pairwise LD (r²) matrix for HLA alleles | This study (`4_compute_regression.R`) |

### Output files

| Location | File | Description |
|---|---|---|
| `figures/` | `figure2_raw.pdf` | Allele frequencies (CHC vs. controls), by locus |
| | `figure3_raw.pdf` | LD matrix for the 6 FDR-significant alleles |
| | `figure4_raw.pdf` | CER corrplot for the 6 FDR-significant alleles |
| | `figure5_raw.pdf` | CER for all 194 allele-disease associations |
| | `figure6.pdf` | Literature / CHC effect sizes for AD- and FDR-significant alleles |
| | `suppfigure1_raw.pdf` | Literature ORs per allele-trait pair |
| `tables/` | `table1.txt` | Table 1: the 6 FDR-significant alleles |
| | `supptable1.csv` | Supplementary Data 1: full allele-trait association table |
| | `supptable2.csv` | Supplementary Data 2: sex-corrected CHC associations |

Figure 1 (study workflow schematic) is not generated by this script.

### Software

See [`SOFTWARE.md`](SOFTWARE.md) for R package requirements and installation instructions.

---

## Contact

For questions about the code or analysis, contact the corresponding author:
Henne Holstege — h.holstege@amsterdamumc.nl
