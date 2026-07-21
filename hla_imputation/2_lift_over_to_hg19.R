## =============================================================================
## 2_lift_over_to_hg19.R — Step 2: LiftOver from hg38 to hg19
##
## Translates variant coordinates in the PLINK1 BED fileset produced by
## Step 1 from GRCh38 to GRCh37 (hg19), which is required by the HIBAG
## pre-fitted model used in Step 3. Variants that fail to map or produce
## duplicate IDs after liftover are excluded.
##
## Usage:
##   Rscript 2_lift_over_to_hg19.R <file_dir> <file_prefix> <chain_file>
##
## Arguments:
##   file_dir    : directory containing the BED/BIM/FAM input files;
##                 also used as the output directory
##                 e.g. /path/to/imputation_results/
##   file_prefix : BED/BIM/FAM basename without extension
##                 e.g. chr6_hla_imputation
##   chain_file  : path to the UCSC hg38ToHg19.over.chain file
##
## Output (in file_dir):
##   <prefix>_updated.bed/.bim/.fam          — hg19 BED fileset for Step 3
##   <prefix>_hg19.bed/.bim/.fam             — intermediate post-extract
##   <prefix>_failed_liftOver_snps.bim       — variants that could not be mapped
##   <prefix>lifted_over_snps.txt            — IDs of successfully mapped variants
##   <prefix>lifted_over_snps_set_for_update.txt — position table for plink
##                                               --update-map / --update-cm
##
## Dependencies:
##   R packages : data.table, liftOver (Bioconductor), GenomicRanges
##   System     : plink (v1)
##   Data       : hg38ToHg19.over.chain (UCSC)
## =============================================================================

library(data.table)
library(liftOver)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript 2_lift_over_to_hg19.R <file_dir> <file_prefix> <chain_file>")

file_folder     <- args[1]
file            <- args[2]
output_path     <- args[1]   # input and output share the same directory
chain_file_path <- args[3]

# Read BIM file (columns: chr, id, cm, pos, a1, a2 — no header)
cat("** Reading BIM file...\n")
positions_hg38 <- fread(file.path(file_folder, paste0(file, ".bim")))
colnames(positions_hg38) <- c("chr", "id", "cm", "pos", "a1", "a2")
cat("  Variants in BIM:", nrow(positions_hg38), "\n")

# Build GRanges from hg38 positions and apply the liftOver chain
cat("** Running liftOver...\n")
df_hg38 <- data.frame(
  chr   = paste0("chr", positions_hg38$chr),
  start = positions_hg38$pos,
  end   = positions_hg38$pos + 1
)
gr_hg38  <- makeGRangesFromDataFrame(df_hg38)
chain    <- import.chain(chain_file_path)
gr_hg19  <- liftOver(gr_hg38, chain)
df_hg19  <- as.data.frame(gr_hg19)
cat("  Variants successfully lifted:", nrow(df_hg19), "\n")

# Reconstruct mapping table joining hg19 positions back to original metadata.
# df_hg19$group is a 1-based index into the original GRanges.
sb <- data.frame(
  chr      = df_hg19$seqnames,
  pos_hg19 = df_hg19$start,
  pos_hg38 = df_hg38[df_hg19$group, "start"],
  rsid     = positions_hg38$id[df_hg19$group],
  a1       = positions_hg38$a1[df_hg19$group],
  a2       = positions_hg38$a2[df_hg19$group]
)

# Identify and remove variants that failed liftover
missing_snps <- positions_hg38[!(positions_hg38$id %in% sb$rsid), ]
cat("  Variants failed to map:", nrow(missing_snps), "\n")

# Remove post-liftover duplicates: one hg38 position mapping to multiple hg19
# positions (e.g. segmental duplications). All instances of a duplicated ID
# are excluded, as PLINK cannot handle duplicate variant identifiers.
dup_ids <- sb$rsid[duplicated(sb$rsid)]
cat("  Variants with duplicate IDs (excluded):", length(unique(dup_ids)), "\n")
sb <- sb[!(sb$rsid %in% dup_ids), ]
sb$cm_pos <- 0   # centimorgan values updated by plink --update-cm below
cat("  Variants retained:", nrow(sb), "\n")

# Write QC and mapping files
write.table(missing_snps,
            file.path(output_path, paste0(file, "_failed_liftOver_snps.bim")),
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

keep_out   <- file.path(output_path, paste0(file, "lifted_over_snps.txt"))
update_out <- file.path(output_path, paste0(file, "lifted_over_snps_set_for_update.txt"))
write.table(sb$rsid, keep_out,   quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(sb,      update_out, quote = FALSE, row.names = FALSE, sep = "\t")

# PLINK step a: extract mapped variants (still with hg38 positions)
cat("** Running plink --extract...\n")
system(paste0(
  "plink --bfile ", file.path(file_folder, file),
  " --extract ", keep_out,
  " --make-bed --out ", file.path(output_path, paste0(file, "_hg19"))
))

# PLINK step b: update bp positions and centimorgan values to hg19.
# Column indices in update_out: col 2 = pos_hg19, col 4 = rsid, col 7 = cm_pos
cat("** Running plink --update-map / --update-cm...\n")
system(paste0(
  "plink --bfile ", file.path(output_path, paste0(file, "_hg19")),
  " --update-cm ",  update_out, " 7 4",
  " --update-map ", update_out, " 2 4",
  " --make-bed --out ", file.path(output_path, paste0(file, "_updated"))
))

cat("** Step 2 complete. hg19 output:", file.path(output_path, paste0(file, "_updated")), "\n")
