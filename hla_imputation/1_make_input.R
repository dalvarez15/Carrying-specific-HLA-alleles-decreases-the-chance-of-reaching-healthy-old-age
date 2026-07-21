## =============================================================================
## 1_make_input.R — Step 1: SNP extraction and PLINK2 → PLINK1 conversion
##
## Prepares chromosome 6 genotype data (PLINK2 pgen/pvar/psam) for HLA
## allele imputation. Fixes missing variant IDs (ID == "."), writes the
## selected variant list, and converts to PLINK1 BED format required by
## liftOver and HIBAG.
##
## Usage:
##   Rscript 1_make_input.R <chr6_plink2_prefix> <output_dir>
##
## Arguments:
##   chr6_plink2_prefix : PLINK2 fileset prefix for chromosome 6 (no extension)
##                        e.g. /path/to/chr6.dose
##   output_dir         : directory for all output files (created if absent)
##
## Output (in output_dir):
##   chr6_hla_imputation.bed/.bim/.fam  — PLINK1 fileset for Step 2
##   snps_hla_imputation.txt            — variant IDs passed to --extract
##
## Dependencies:
##   R packages : data.table
##   System     : plink2
## =============================================================================

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript 1_make_input.R <chr6_plink2_prefix> <output_dir>")

input_file    <- args[1]
output_folder <- args[2]

if (!dir.exists(output_folder)) {
  cat("** Creating output folder:", output_folder, "\n")
  dir.create(output_folder, recursive = TRUE)
}

# Copy PLINK2 files to the output folder so downstream steps can modify in-place
cat("** Copying input PLINK2 files...\n")
system(paste0("cp ", input_file, ".p* ", output_folder))

# Read the copied variant file
cat("** Reading .pvar...\n")
pvar_out <- file.path(output_folder, paste0(basename(input_file), ".pvar"))
data <- fread(pvar_out, showProgress = FALSE)

# Fix unnamed variants: assign CHROM:POS:REF:ALT to any variant with ID == "."
# These arise in imputed datasets and must be resolved before plink --extract
cat("** Fixing missing variant IDs...\n")
data$index <- seq_len(nrow(data))
missid <- data[data$ID == ".", ]
if (nrow(missid) > 0) {
  cat("  ", nrow(missid), "variants with ID == '.'; assigning CHROM:POS:REF:ALT\n")
  missid$ID <- paste(missid$"#CHROM", missid$POS, missid$REF, missid$ALT, sep = ":")
}
data <- rbind(data[data$ID != ".", ], missid)
data <- data[order(data$index), ]
data$index <- NULL

# Write updated .pvar, preserving the original header comment lines
cat("** Writing updated .pvar...\n")
pvar_tmp <- paste0(pvar_out, ".tmp")
system(paste0("grep '^#' ", pvar_out, " > ", pvar_tmp))
write.table(data, file = pvar_tmp, row.names = FALSE, col.names = FALSE,
            quote = FALSE, sep = "\t", append = TRUE)
system(paste0("mv ", pvar_tmp, " ", pvar_out))

# Re-read to confirm
data <- fread(pvar_out, showProgress = FALSE)
cat("  Total variants:", nrow(data), "\n")

# QC filters — disabled by default; uncomment to activate
# Filter 1: imputation quality (R2 > 0.8 recommended for HIBAG)
#   data$R2 <- as.numeric(sub(".*R2=([^;]+).*", "\\1", data$INFO))
#   data <- data[data$R2 > 0.8, ]
# Filter 2: restrict to HLA window on hg38 (chr6: 25–40 Mb)
#   data <- data[data$POS > 25000000 & data$POS < 40000000, ]

cat("** Variants selected for imputation:", nrow(data), "\n")

# Write variant list and convert to PLINK1 BED
snp_list <- file.path(output_folder, "snps_hla_imputation.txt")
write.table(data$ID, file = snp_list, row.names = FALSE, col.names = FALSE, quote = FALSE)

cat("** Running plink2 --make-bed...\n")
system(paste0(
  "plink2",
  " --pfile ", file.path(output_folder, basename(input_file)),
  " --extract ", snp_list,
  " --make-bed",
  " --out ", file.path(output_folder, "chr6_hla_imputation")
))

cat("** Step 1 complete.\n")
