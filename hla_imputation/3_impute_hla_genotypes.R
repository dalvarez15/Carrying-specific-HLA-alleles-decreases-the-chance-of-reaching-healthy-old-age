## =============================================================================
## 3_impute_hla_genotypes.R — Step 3: HIBAG HLA allele imputation
##
## Imputes two-field HLA alleles at seven classical loci (A, B, C, DRB1,
## DQB1, DQA1, DPB1) using HIBAG and a pre-fitted European-ancestry
## classifier. The PLINK1 BED fileset in hg19 coordinates produced by
## Step 2 is used as input.
##
## Usage:
##   Rscript 3_impute_hla_genotypes.R <plink_prefix> <model_file> <output_dir>
##
## Arguments:
##   plink_prefix : path prefix (no extension) for the hg19 BED/BIM/FAM
##                  fileset produced by Step 2
##                  e.g. /path/to/output/chr6_hla_imputation_updated
##   model_file   : HIBAG pre-fitted model (.RData); must be a named list
##                  with one hlaAttrBagObj per locus ($A, $B, $C, etc.)
##                  Download platform- and ancestry-matched models from:
##                  https://hibag.s3.amazonaws.com/hlares_index.html
##   output_dir   : directory for result and summary files (must exist)
##
## Output (in output_dir), one pair of files per locus:
##   result_<LOCUS>.txt  — per-sample allele calls and posterior probabilities
##                          columns: sample.id, allele1, allele2, prob
##   summary_<LOCUS>.txt — imputation quality summary (per-allele accuracy)
##
## Note: calls with posterior probability < 0.5 are excluded in Step 4
## (4_compute_regression.R), as recommended by HIBAG (Zheng et al. 2014).
##
## Dependencies:
##   R packages : HIBAG (Bioconductor)
## =============================================================================

library(HIBAG)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript 3_impute_hla_genotypes.R <plink_prefix> <model_file> <output_dir>")

plink_prefix <- args[1]
model_file   <- args[2]
output_path  <- args[3]

# Validate inputs
for (ext in c(".bed", ".bim", ".fam")) {
  if (!file.exists(paste0(plink_prefix, ext)))
    stop("Input file not found: ", plink_prefix, ext)
}
if (!file.exists(model_file))  stop("Model file not found: ", model_file)
if (!dir.exists(output_path))  stop("Output directory not found: ", output_path)

# Load genotype data and pre-fitted model
cat("** Loading genotype data...\n")
geno <- hlaBED2Geno(
  bed.fn = paste0(plink_prefix, ".bed"),
  fam.fn = paste0(plink_prefix, ".fam"),
  bim.fn = paste0(plink_prefix, ".bim")
)

cat("** Loading HIBAG model:", model_file, "\n")
mod_lst <- get(load(model_file))
cat("  Loci in model:", paste(names(mod_lst), collapse = ", "), "\n")

loci <- c("A", "B", "C", "DRB1", "DQB1", "DQA1", "DPB1")

for (locus in loci) {
  cat("\n** Imputing HLA-", locus, "...\n", sep = "")

  if (!(locus %in% names(mod_lst))) {
    warning("Locus ", locus, " not found in model; skipping.")
    next
  }

  # Load locus classifier and run parallel prediction (8 cores; match --cpus-per-task)
  classifier <- hlaModelFromObj(mod_lst[[locus]])
  result     <- hlaPredict(classifier, geno, cl = 8)

  summary(result)
  head(result$value)

  # Write per-sample allele calls (columns: sample.id, allele1, allele2, prob)
  write.table(result$value,
              file  = file.path(output_path, paste0("result_", locus, ".txt")),
              sep   = "\t", quote = FALSE, row.names = FALSE)

  # Write imputation quality summary
  write.table(capture.output(summary(result)),
              file      = file.path(output_path, paste0("summary_", locus, ".txt")),
              quote     = FALSE, row.names = FALSE, col.names = FALSE)
}

cat("\n** Step 3 complete. Results written to:", output_path, "\n")
