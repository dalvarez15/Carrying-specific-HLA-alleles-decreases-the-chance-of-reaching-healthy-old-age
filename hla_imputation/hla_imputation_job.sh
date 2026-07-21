#!/bin/bash
#SBATCH --job-name=hla_imputation
#SBATCH --output=hla_imputation_%j.log
#SBATCH --error=hla_imputation_%j.err
#SBATCH --cpus-per-task=8

# =============================================================================
# hla_imputation_job.sh — SLURM submission script for the HLA imputation pipeline
#
# Runs Steps 1–3 sequentially. Step 4 (regression and LD matrix) requires
# additional input files (phenotype file, PC file, allele dosage table) and
# is configured separately inside 4_compute_regression.R; submit it as a
# second job once Steps 1–3 are complete.
#
# Edit the "Configuration" block below before submitting:
#   sbatch hla_imputation_job.sh
# =============================================================================

# --- Configuration -----------------------------------------------------------
SCRIPT_DIR="/path/to/hla_imputation"          # directory containing the R scripts
CHR6_INPUT="/path/to/chr6.dose"               # PLINK2 chr6 prefix (no extension)
OUTPUT_DIR="/path/to/imputation_results"      # intermediate and imputation outputs
BED_PREFIX="chr6_hla_imputation"              # prefix produced by Step 1
CHAIN_FILE="/path/to/hg38ToHg19.over.chain"  # UCSC hg38 → hg19 chain file
MODEL_FILE="/path/to/InfiniumGlobal-European-HLA4-hg19.RData"  # HIBAG model
# -----------------------------------------------------------------------------

conda activate r-hibag

mkdir -p "${OUTPUT_DIR}/imputed_HLA_alleles"

echo "=== Step 1: make input ==="
Rscript "${SCRIPT_DIR}/1_make_input.R" "${CHR6_INPUT}" "${OUTPUT_DIR}"
[ $? -ne 0 ] && { echo "Step 1 failed" >&2; exit 1; }

echo "=== Step 2: liftover ==="
Rscript "${SCRIPT_DIR}/2_lift_over_to_hg19.R" "${OUTPUT_DIR}" "${BED_PREFIX}" "${CHAIN_FILE}"
[ $? -ne 0 ] && { echo "Step 2 failed" >&2; exit 1; }

echo "=== Step 3: HIBAG imputation ==="
Rscript "${SCRIPT_DIR}/3_impute_hla_genotypes.R" \
    "${OUTPUT_DIR}/${BED_PREFIX}_updated" \
    "${MODEL_FILE}" \
    "${OUTPUT_DIR}/imputed_HLA_alleles"
[ $? -ne 0 ] && { echo "Step 3 failed" >&2; exit 1; }

echo "=== Steps 1–3 complete. Results in ${OUTPUT_DIR}/imputed_HLA_alleles ==="
echo "Next: configure and run 4_compute_regression.R to generate inputs/ files."
