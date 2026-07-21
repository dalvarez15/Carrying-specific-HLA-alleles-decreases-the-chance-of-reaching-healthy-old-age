## =============================================================================
## 4_compute_regression.R — Step 4: logistic regression and LD matrix
##
## Produces the three inputs consumed by main_script.R:
##   inputs/regression_npj_aging_2025.csv  — per-allele logistic regression
##                                           results (CHC vs. controls)
##   inputs/hla.ld                         — pairwise r² matrix for HLA alleles
##   inputs/hla.bim                        — allele identifiers (PLINK BIM format)
##
## The script has two parts:
##
##   PART 1 — REGRESSION
##     Loads HIBAG imputation results (Step 3) for two genotyping batches,
##     merges with phenotype and PC files, builds per-sample allele dosages,
##     and runs six logistic regression models per allele:
##       1. cent_ctr_original             — PCs only
##       2. cent_ctr_sex_covariate        — PCs + sex
##       3. cent_ctr_sex_interaction_main — PCs + sex*allele (main effect)
##       4. cent_ctr_sex_interaction_term — PCs + sex*allele (interaction term)
##       5. cent_ctr_female_only          — females only, PCs
##       6. cent_ctr_male_only            — males only, PCs
##     Outcome: 1 = control, 0 = centenarian (CHC).
##
##   PART 2 — LD MATRIX
##     Converts a wide allele dosage table to transposed PLINK dosage format
##     and calls plink2 / plink to compute a square r² matrix.
##
## Usage:
##   Rscript 4_compute_regression.R
##   (paths are configured in the section below)
##
## Input files:
##   Batch 1 imputation results : result_<LOCUS>.txt per locus (from Step 3)
##                                sample.id = bare GWAS ID
##   Batch 2 imputation results : result_<LOCUS>.txt per locus (from Step 3)
##                                sample.id = compound string (PREFIX_GWASID)
##   Phenotype file             : tab-separated; columns must include
##                                ID_GWAS, ID_100plus, Study, diagnosis, chip, sex, age
##   Updated phenotype file     : used to recover missing sex values;
##                                columns must include ID_GWAS, ID_ORIGINAL, sex
##   PC file                    : PLINK2 eigenvec output (headerless);
##                                col 1 = FID, col 2 = IID, cols 3–7 = PC1–PC5
##   Allele dosage table        : semicolon-separated CSV for LD computation;
##                                col 1 = ID_GWAS, cols 2–14 = phenotype/PC fields,
##                                cols 15+ = allele dosage columns (adjust first_allele_col)
##
## Output (copy to inputs/ in the repository root after running):
##   regression_npj_aging_2025.csv, hla.ld, hla.bim
##
## Dependencies:
##   R packages : tidyr, dplyr, stringr, data.table
##   System     : plink (v1), plink2
## =============================================================================

library(tidyr)
library(dplyr)
library(stringr)
library(data.table)

## -----------------------------------------------------------------------------
## Configuration — set all paths here before running
## -----------------------------------------------------------------------------

# Batch 1: HIBAG result directory (result_A.txt, result_B.txt, etc.)
batch1_dir <- "/path/to/batch1/imputation_results"

# Batch 2: HIBAG result directory
# In this batch, sample.id is a compound string (PREFIX_GWASID);
# the GWAS ID is extracted as the second field after splitting on "_".
# Verify this matches your actual ID format before running.
batch2_dir <- "/path/to/batch2/imputation_results"

# Phenotype file (tab-separated)
phenotype_file <- "/path/to/phenotypes.txt"

# Updated phenotype file for recovering missing sex values
phenotype_file_new <- "/path/to/phenotypes_updated.txt"

# PLINK eigenvec file (headerless; col 2 = sample ID, cols 3–7 = PC1–PC5)
pc_file <- "/path/to/pca.eigenvec"

# Wide allele dosage CSV for LD computation (Part 2; semicolon-separated)
allele_dosage_file <- "/path/to/allele_dosage_table.csv"

# Chip label string used in the phenotype file to identify Batch 1 samples
# that appear in both batches (these are removed from Batch 2)
batch1_chip_label <- "GSA_with_custom_content"

# First column index containing allele dosages in allele_dosage_file
first_allele_col <- 15

# Output directory for all generated files
output_dir <- "/path/to/output"

## -----------------------------------------------------------------------------

loci <- c("A", "B", "C", "DRB1", "DQB1", "DQA1", "DPB1")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)


## =============================================================================
## PART 1 — REGRESSION
## =============================================================================

## -----------------------------------------------------------------------------
## Step 1: Load and standardise Batch 1 imputation results
##
##   result_<LOCUS>.txt files have columns: sample.id, allele1, allele2, prob,
##   matching. Allele names are prefixed with the locus name, and "HLA_" is
##   stripped to yield the canonical two-field form (e.g. "A*02:01", "DRB1*01:01")
##   used throughout the regression and in main_script.R.
## -----------------------------------------------------------------------------

cat("** Loading Batch 1 imputation results...\n")

batch1_list <- lapply(loci, function(locus) {
  df <- read.delim(file.path(batch1_dir, paste0("result_", locus, ".txt")))
  df$Locus <- ifelse(locus %in% c("A", "B", "C"), paste0("HLA_", locus), locus)
  df
})
imputed_hla <- do.call(rbind, batch1_list)
colnames(imputed_hla)[colnames(imputed_hla) == "sample.id"] <- "ID_GWAS"

imputed_hla <- imputed_hla %>%
  mutate(allele1 = gsub("^HLA_", "", paste0(Locus, "*", allele1)),
         allele2 = gsub("^HLA_", "", paste0(Locus, "*", allele2)))

## -----------------------------------------------------------------------------
## Step 2: Load and standardise Batch 2 imputation results
##
##   sample.id is a compound string (e.g. "PREFIX_18R1234"); the GWAS ID is
##   extracted as the second field after splitting on "_".
## -----------------------------------------------------------------------------

cat("** Loading Batch 2 imputation results...\n")

batch2_list <- lapply(loci, function(locus) {
  df <- read.table(file.path(batch2_dir, paste0("result_", locus, ".txt")), header = TRUE)
  df$ID_GWAS <- str_split_fixed(df$sample.id, "_", 2)[, 2]
  df <- df[, c("ID_GWAS", "allele1", "allele2", "prob", "matching")]
  df$Locus <- ifelse(locus %in% c("A", "B", "C"), paste0("HLA_", locus), locus)
  df
})
hla_batch2 <- do.call(rbind, batch2_list)

hla_batch2 <- hla_batch2 %>%
  mutate(allele1 = gsub("^HLA_", "", paste0(Locus, "*", allele1)),
         allele2 = gsub("^HLA_", "", paste0(Locus, "*", allele2)))

## -----------------------------------------------------------------------------
## Step 3: Merge with phenotype file and resolve cross-batch duplicates
##
##   Some samples are genotyped on both chips. Strategy: prefer Batch 2 calls
##   for overlapping samples.
##     a. Attach chip labels by merging with the phenotype file.
##     b. Remove from Batch 2 any sample carrying the Batch 1 chip label.
##     c. Remove from Batch 1 any sample that also appears in Batch 2.
##     d. Remove any sample with a duplicated ID_100plus (same participant
##        enrolled under two GWAS IDs); keep the first occurrence.
## -----------------------------------------------------------------------------

cat("** Merging with phenotype file...\n")
phenotypes <- read.delim(phenotype_file)

imputed_hla <- merge(imputed_hla, phenotypes, by = "ID_GWAS")
hla_batch2  <- merge(hla_batch2,  phenotypes, by = "ID_GWAS")

hla_batch2      <- hla_batch2[hla_batch2$chip != batch1_chip_label, ]
ids_in_batch2   <- unique(hla_batch2$ID_GWAS)
imputed_hla     <- imputed_hla[!imputed_hla$ID_GWAS %in% ids_in_batch2, ]

imputed_hla_combined <- rbind(imputed_hla, hla_batch2)

# Detect and remove duplicate ID_100plus values (one per participant × 7 loci = 7 rows;
# > 14 rows for the same ID_100plus indicates a true sample duplicate)
dup_100plus <- names(which(table(
  imputed_hla_combined$ID_100plus[!is.na(imputed_hla_combined$ID_100plus)]
) > 14))
if (length(dup_100plus) > 0) {
  cat("  Duplicated ID_100plus found:", paste(dup_100plus, collapse = ", "), "\n")
  for (dup_id in dup_100plus) {
    gwas_ids <- unique(imputed_hla_combined$ID_GWAS[
      imputed_hla_combined$ID_100plus == dup_id & !is.na(imputed_hla_combined$ID_100plus)
    ])
    if (length(gwas_ids) > 1)
      imputed_hla_combined <- imputed_hla_combined[imputed_hla_combined$ID_GWAS != gwas_ids[2], ]
  }
}
cat("  Unique samples after deduplication:", length(unique(imputed_hla_combined$ID_GWAS)), "\n")

## -----------------------------------------------------------------------------
## Step 4: Quality filter — retain calls with prob >= 0.5
##
##   Posterior probability >= 0.5 is the standard HIBAG threshold
##   (Zheng et al. 2014, Pharmacogenomics J).
## -----------------------------------------------------------------------------

imputed_hla_filtered <- imputed_hla_combined[
  !is.na(imputed_hla_combined$allele1) & imputed_hla_combined$prob >= 0.5, ]
cat("  Samples after prob >= 0.5 filter:", length(unique(imputed_hla_filtered$ID_GWAS)), "\n")

## -----------------------------------------------------------------------------
## Step 5: Load principal components and merge
## -----------------------------------------------------------------------------

cat("** Loading PCs...\n")
pcs <- read.table(pc_file, header = FALSE, stringsAsFactors = FALSE)
pcs <- pcs[, 2:7]
colnames(pcs) <- c("ID_GWAS", "PC1", "PC2", "PC3", "PC4", "PC5")

regression <- merge(imputed_hla_filtered, pcs, by = "ID_GWAS", all.x = TRUE)

## -----------------------------------------------------------------------------
## Step 6: Build per-sample allele dosage table
##
##   For each locus, pivot allele1/allele2 into a wide dosage matrix where each
##   column is one allele and the value is the number of copies carried (0/1/2).
##   This is the standard additive dosage encoding used in the regression.
##   NA values indicate the sample was not imputed at that locus.
## -----------------------------------------------------------------------------

cat("** Building allele dosage table...\n")

allele_dosages <- unique(regression[, c("ID_GWAS", "ID_100plus", "Study",
                                        "diagnosis", "chip", "sex", "age",
                                        "PC1", "PC2", "PC3", "PC4", "PC5")])

for (locus in unique(regression$Locus)) {
  pivot_wide <- subset(regression, Locus == locus) %>%
    pivot_longer(cols = c(allele1, allele2), values_to = "allele") %>%
    group_by(ID_GWAS, allele) %>%
    summarise(count = n(), .groups = "drop") %>%
    pivot_wider(names_from = allele, values_from = count, values_fill = 0)
  allele_dosages <- merge(allele_dosages, pivot_wide, by = "ID_GWAS", all.x = TRUE)
}
cat("  Samples in dosage table:", length(unique(allele_dosages$ID_GWAS)), "\n")

## -----------------------------------------------------------------------------
## Step 7: Define analysis subset and binary outcome
##
##   Centenarians (CHC) are coded 0; all control groups are coded 1.
##   Adjust case_label and control_labels to match your phenotype file.
## -----------------------------------------------------------------------------

case_label     <- "Centenarian"
control_labels <- c("Control_100plus", "Control_LASA", "Control_other_twin",
                    "Control_path", "SCD", "Control_MS")

allele_dosages_sub <- subset(allele_dosages, diagnosis %in% c(case_label, control_labels))
allele_dosages_sub$phenotype <- ifelse(allele_dosages_sub$diagnosis == case_label, 0, 1)
cat("  CHC (0):", sum(allele_dosages_sub$phenotype == 0),
    " Controls (1):", sum(allele_dosages_sub$phenotype == 1), "\n")

## -----------------------------------------------------------------------------
## Step 8: Fill in missing sex values from the updated phenotype file
##
##   First pass matches on ID_GWAS; second pass falls back to ID_ORIGINAL
##   for samples absent from the primary ID column.
## -----------------------------------------------------------------------------

cat("** Filling missing sex values...\n")
phenotypes_new <- read.delim(phenotype_file_new)

missing_idx <- which(!(allele_dosages_sub$sex %in% c("F", "M")) | is.na(allele_dosages_sub$sex))
allele_dosages_sub$sex[missing_idx] <- phenotypes_new$sex[
  match(allele_dosages_sub$ID_GWAS[missing_idx], phenotypes_new$ID_GWAS)]

missing_idx2 <- which(is.na(allele_dosages_sub$sex))
allele_dosages_sub$sex[missing_idx2] <- phenotypes_new$sex[
  match(allele_dosages_sub$ID_GWAS[missing_idx2], phenotypes_new$ID_ORIGINAL)]

cat("  Sex missing after fill:", sum(is.na(allele_dosages_sub$sex)), "\n")

## -----------------------------------------------------------------------------
## Step 9: Helper functions for regression
## -----------------------------------------------------------------------------

# Safely extract a named coefficient, SE, and p-value from a fitted glm.
# Returns NAs if the term is absent (e.g. perfect separation).
safe_coef <- function(model, coef_pattern) {
  coefs <- coef(summary(model))
  idx   <- grep(coef_pattern, rownames(coefs))
  if (length(idx) == 0) return(list(Estimate = NA, Std.Error = NA, P.value = NA))
  list(Estimate  = coefs[idx[1], "Estimate"],
       Std.Error = coefs[idx[1], "Std. Error"],
       P.value   = coefs[idx[1], "Pr(>|z|)"])
}

# Compute allele frequency statistics from a data frame with columns:
# dummy (allele dosage 0/1/2) and phenotype (0 = CHC, 1 = control).
allele_freq_stats <- function(df) {
  num_cases     <- sum(df$phenotype == 1)
  num_controls  <- sum(df$phenotype == 0)
  homo_cases    <- sum(df$phenotype == 1 & df$dummy == 2)
  homo_controls <- sum(df$phenotype == 0 & df$dummy == 2)
  het_cases     <- sum(df$phenotype == 1 & df$dummy == 1)
  het_controls  <- sum(df$phenotype == 0 & df$dummy == 1)
  count_cases   <- 2 * homo_cases    + het_cases
  count_controls<- 2 * homo_controls + het_controls
  list(num_cases      = num_cases,
       num_controls   = num_controls,
       c_cases        = sum(df$phenotype == 1 & df$dummy > 0),
       c_controls     = sum(df$phenotype == 0 & df$dummy > 0),
       homo_cases     = homo_cases,     homo_controls  = homo_controls,
       het_cases      = het_cases,      het_controls   = het_controls,
       count_cases    = count_cases,    count_controls = count_controls,
       freq_cases     = count_cases    / (2 * num_cases),
       freq_controls  = count_controls / (2 * num_controls),
       freq_total     = (count_cases + count_controls) / (2 * (num_cases + num_controls)))
}

# Build one results row from a glm coefficient result, frequency stats,
# allele metadata, and an association label.
make_allele_row <- function(locus, allele, coef_res, freq, label) {
  data.frame(
    locus          = locus,
    allele         = allele,
    beta           = coef_res$Estimate,
    standard_error = coef_res$Std.Error,
    odds_ratio     = exp(coef_res$Estimate),
    lower_ci       = exp(coef_res$Estimate - 1.96 * coef_res$Std.Error),
    upper_ci       = exp(coef_res$Estimate + 1.96 * coef_res$Std.Error),
    p_value        = coef_res$P.value,
    num_cases      = freq$num_cases,    num_controls   = freq$num_controls,
    c_cases        = freq$c_cases,      c_controls     = freq$c_controls,
    homo_cases     = freq$homo_cases,   homo_controls  = freq$homo_controls,
    het_cases      = freq$het_cases,    het_controls   = freq$het_controls,
    count_cases    = freq$count_cases,  count_controls = freq$count_controls,
    freq_cases     = freq$freq_cases,
    freq_controls  = freq$freq_controls,
    freq_total     = freq$freq_total,
    association    = label
  )
}

## -----------------------------------------------------------------------------
## Step 10: Logistic regression — six models per allele
##
##   Allele dosage is treated as a continuous additive predictor ("dummy").
##   PC1–5 correct for population stratification in all models.
## -----------------------------------------------------------------------------

cat("** Running logistic regression across all alleles...\n")

allele_columns  <- grep("\\*", colnames(allele_dosages_sub), value = TRUE)
coefficients_df <- data.frame()

for (allele in allele_columns) {
  locus <- strsplit(allele, "\\*")[[1]][1]

  # Model 1: PCs only
  df1      <- allele_dosages_sub[!is.na(allele_dosages_sub[[allele]]), ]
  df1$dummy <- df1[[allele]]
  m1       <- glm(phenotype ~ dummy + PC1 + PC2 + PC3 + PC4 + PC5, data = df1, family = "binomial")
  coefficients_df <- rbind(coefficients_df,
    make_allele_row(locus, allele, safe_coef(m1, "dummy"), allele_freq_stats(df1), "cent_ctr_original"))

  # Models 2–4: require non-missing sex
  df2      <- df1[df1$sex %in% c("F", "M"), ]

  # Model 2: sex as covariate
  m2       <- glm(phenotype ~ dummy + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = df2, family = "binomial")
  coefficients_df <- rbind(coefficients_df,
    make_allele_row(locus, allele, safe_coef(m2, "dummy"), allele_freq_stats(df2), "cent_ctr_sex_covariate"))

  # Models 3 & 4: sex interaction (main effect and interaction term saved separately)
  m3       <- glm(phenotype ~ dummy * sex + PC1 + PC2 + PC3 + PC4 + PC5, data = df2, family = "binomial")
  coefficients_df <- rbind(coefficients_df,
    make_allele_row(locus, allele, safe_coef(m3, "^dummy$"), allele_freq_stats(df2), "cent_ctr_sex_interaction_main"))

  interaction_term <- grep("^dummy:sex", rownames(coef(summary(m3))), value = TRUE)
  if (length(interaction_term) > 0)
    coefficients_df <- rbind(coefficients_df,
      make_allele_row(locus, allele, safe_coef(m3, interaction_term[1]),
                      allele_freq_stats(df2), "cent_ctr_sex_interaction_term"))

  # Model 5: females only
  df_f <- df2[df2$sex == "F", ]
  if (nrow(df_f) > 0) {
    m_f <- glm(phenotype ~ dummy + PC1 + PC2 + PC3 + PC4 + PC5, data = df_f, family = "binomial")
    coefficients_df <- rbind(coefficients_df,
      make_allele_row(locus, allele, safe_coef(m_f, "dummy"), allele_freq_stats(df_f), "cent_ctr_female_only"))
  }

  # Model 6: males only
  df_m <- df2[df2$sex == "M", ]
  if (nrow(df_m) > 0) {
    m_m <- glm(phenotype ~ dummy + PC1 + PC2 + PC3 + PC4 + PC5, data = df_m, family = "binomial")
    coefficients_df <- rbind(coefficients_df,
      make_allele_row(locus, allele, safe_coef(m_m, "dummy"), allele_freq_stats(df_m), "cent_ctr_male_only"))
  }
}

## -----------------------------------------------------------------------------
## Step 11: Write regression results
##
##   Semicolon-separated (write.csv2) to match the format read by main_script.R
##   (read.csv2). Copy the output file to inputs/ in the repository root.
## -----------------------------------------------------------------------------

regression_out <- file.path(output_dir, "regression_npj_aging_2025.csv")
cat("** Writing regression results to:", regression_out, "\n")
write.csv2(coefficients_df, file = regression_out, row.names = FALSE)
cat("** Part 1 complete.\n\n")


## =============================================================================
## PART 2 — LD MATRIX
##
##   Converts the wide allele dosage table to transposed PLINK dosage format
##   (.traw / .tfam) and computes a square pairwise r² matrix.
##   Outputs hla.ld (r² matrix) and hla.bim (allele identifiers) are used
##   by main_script.R for Figure 3 (LD heatmap) and Figure 6.
## =============================================================================

cat("** Starting Part 2: LD matrix computation...\n")

alleles_snps <- fread(allele_dosage_file, header = TRUE, sep = ";", stringsAsFactors = FALSE)

sb           <- alleles_snps[, c(1, first_allele_col:ncol(alleles_snps)), with = FALSE]
sample_names <- sb$ID_GWAS
sb$ID_GWAS   <- NULL

# Transpose: rows = alleles, columns = samples (required by plink --import-dosage)
sbt           <- data.frame(t(sb))
colnames(sbt) <- paste0("0_", sample_names)   # FID_IID format for id-delim=_

sbt$SNP <- rownames(sbt)
sbt$A1  <- "P"
sbt$A2  <- "A"
sbt     <- sbt[, c("SNP", "A1", "A2", setdiff(colnames(sbt), c("SNP", "A1", "A2")))]

traw_out    <- file.path(output_dir, "hla.traw")
tfam_out    <- file.path(output_dir, "hla.tfam")
pgen_prefix <- file.path(output_dir, "hla")

cat("** Writing transposed dosage file...\n")
write.table(sbt, traw_out, quote = FALSE, row.names = FALSE, sep = "\t")
write.table(
  data.frame(FID = 0, IID = sample_names, PAT = 0, MAT = 0, SEX = 0, PHENO = 0),
  tfam_out, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t"
)

cat("** Converting to PLINK2 pgen format...\n")
system(paste0("plink2 --import-dosage ", traw_out,
              " id-delim=_ format=1 --fam ", tfam_out,
              " --real-ref-alleles --make-pgen --out ", pgen_prefix))

cat("** Converting to PLINK1 BED format...\n")
system(paste0("plink2 --pfile ", pgen_prefix, " --make-bed --out ", pgen_prefix))

cat("** Computing pairwise r² matrix...\n")
system(paste0("plink --bfile ", pgen_prefix, " --r2 square --out ", pgen_prefix))

cat("** Part 2 complete.\n")
cat("   Copy", file.path(output_dir, "regression_npj_aging_2025.csv"), "to inputs/\n")
cat("   Copy", file.path(output_dir, "hla.ld"), "to inputs/\n")
cat("   Copy", file.path(output_dir, "hla.bim"), "to inputs/\n")
