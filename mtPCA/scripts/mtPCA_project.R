#!/usr/bin/env Rscript
# =============================================================================
# mtPCA_project.R
#
# Project external cohort mtDNA genotypes into a pre-computed 1000G PCA space.
#
# Use this script when variant overlap between external data and the reference
# is >= 95%. For lower overlap, use mtPCA_recompute.R instead.
#
# External cohort genotypes are read from PLINK BIM/BED/FAM.
#
# Variant harmonisation priority:
#   1. Exact match          (chr:pos:ref:alt identical)
#   2. REF/ALT swapped      (chr:pos:alt:ref) -> recode 0<->1
#   3. Strand complement    (chr:pos:comp(ref):comp(alt)) -> use as-is
#   4. Complement + swap    (chr:pos:comp(alt):comp(ref)) -> recode 0<->1
#   Ambiguous SNPs (A/T or C/G) are dropped with a warning.
#
# Usage:
#   Rscript mtPCA_project.R \
#     --ext_prefix <path to external cohort plink prefix> \
#     --out        <output prefix> \
#     [--ref_dir   <directory containing reference files, default '.'>]
#
# Required reference files (in --ref_dir):
#   mtPCA.variant_ids.txt       -- chr:pos:ref:alt, one per line, no header
#   mtPCA.mu.txt                -- reference mean AF, one per line, no header
#   mtPCA.loadings.PC1-20.txt  -- 317 x 20 matrix, no header or rownames
#   mtPCA.eigenvalues.txt       -- eigenvalues 1-20, no header
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
  library(BEDMatrix)
  library(dplyr)
})

# -----------------------------------------------------------------------------
# 0. Parse arguments
# -----------------------------------------------------------------------------
option_list <- list(
  make_option("--ext_prefix", type = "character",
              help = "Path to external cohort PLINK prefix (bim/bed/fam)"),
  make_option("--ref_dir",    type = "character", default = ".",
              help = "Directory containing reference files [default: .]"),
  make_option("--out",        type = "character", default = "mtPCA_projected",
              help = "Output file prefix [default: mtPCA_projected]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$ext_prefix)) stop("--ext_prefix is required. Run with --help for usage.")

cat("=== mtPCA Projection Workflow ===\n\n")

# -----------------------------------------------------------------------------
# Functions: strand complement and ambiguity check
# -----------------------------------------------------------------------------
complement <- function(allele) {
  chartr("ACGT", "TGCA", allele)
}

is_ambiguous <- function(ref, alt) {
  sorted <- paste(pmin(ref, alt), pmax(ref, alt), sep = "")
  sorted %in% c("AT", "CG")
}

# -----------------------------------------------------------------------------
# Helper: read PLINK dataset using BEDMatrix
# Returns list: geno (samples x variants, 0/1 numeric), bim, fam
# BEDMatrix returns 0/1/2; Recode heteroplasmic 1 > 0; recode 2->1 for haploid mt
# -----------------------------------------------------------------------------
read_plink <- function(prefix) {
  bed <- BEDMatrix(prefix)
  bim <- fread(paste0(prefix, ".bim"), header = FALSE,
               col.names = c("chrom", "snp_id", "cm", "pos", "allele1", "allele2"))
  fam <- fread(paste0(prefix, ".fam"), header = FALSE,
               col.names = c("fid", "iid", "pat", "mat", "sex", "phen"))

  # PLINK BIM: allele1 = ALT (counted allele), allele2 = REF
  bim$ref <- bim$allele2
  bim$alt <- bim$allele1

  # Load full matrix (small for mt SNPs)
  geno <- as.matrix(bed)
  rownames(geno) <- fam$iid

  # Recode heteroplasmic 1 > 0
  # Recode haploid: 2 -> 1
  geno[!is.na(geno) & geno == 1L] <- 0L
  geno[!is.na(geno) & geno == 2L] <- 1L


  # Name columns by chr:pos:ref:alt
  colnames(geno) <- paste(bim$chrom, bim$pos, bim$ref, bim$alt, sep = ":")

  list(geno = geno, bim = bim, fam = fam)
}

# -----------------------------------------------------------------------------
# Helper: harmonise external variants to reference variant universe
#
# For each reference variant (chr:pos:ref:alt), attempts matching in order:
#   1. Exact:           chr:pos:ref:alt
#   2. REF/ALT swap:    chr:pos:alt:ref           -> flip coding
#   3. Complement:      chr:pos:comp(ref):comp(alt)
#   4. Comp + swap:     chr:pos:comp(alt):comp(ref) -> flip coding
#
# Ambiguous SNPs (A/T or C/G) in the reference are dropped before matching.
# Returns data.frame with one row per (non-ambiguous) reference variant.
# -----------------------------------------------------------------------------
harmonise_variants <- function(target_ids, bim) {

  ## check if chr coded as "26" or "MT"
  if (!all(bim$chrom %in% c("26", "MT"))) {
    stop("Unexpected chromosome coding in BIM: expected '26' or 'MT'")
  }
  bim_chr <- case_when(
    all(bim$chrom == "26") ~ "26",
    all(bim$chrom == "MT") ~ "MT")
  


  # Parse target IDs
  target_parsed <- do.call(rbind, strsplit(target_ids, ":"))
  target_df <- data.frame(
    variant_id = paste(bim_chr, target_parsed[, 2], target_parsed[, 3], target_parsed[, 4], sep = ":"),
    chrom      = bim_chr,
    pos        = target_parsed[, 2],
    ref        = target_parsed[, 3],
    alt        = target_parsed[, 4],
    stringsAsFactors = FALSE
  )

  # Flag and drop ambiguous targets
  ambig <- is_ambiguous(target_df$ref, target_df$alt)
  if (any(ambig)) {
    cat(sprintf("  Dropping %d ambiguous SNP(s) (A/T or C/G) from reference list\n",
                sum(ambig)))
    cat(sprintf("  Dropped: %s\n", paste(target_ids[ambig], collapse = ", ")))
    target_df <- target_df[!ambig, ]
  }

  # Build lookups from external BIM
  ds_ref      <- bim$ref
  ds_alt      <- bim$alt
  ds_comp_ref <- complement(ds_ref)
  ds_comp_alt <- complement(ds_alt)
  ds_pos      <- as.character(bim$pos)
  ds_chr      <- bim_chr

  make_id <- function(chr, pos, r, a) paste(chr, pos, r, a, sep = ":")

  col_names        <- make_id(ds_chr, ds_pos, ds_ref, ds_alt)
  lookup_exact     <- setNames(col_names, col_names)
  lookup_swap      <- setNames(col_names, make_id(ds_chr, ds_pos, ds_alt,      ds_ref))
  lookup_comp      <- setNames(col_names, make_id(ds_chr, ds_pos, ds_comp_ref, ds_comp_alt))
  lookup_comp_swap <- setNames(col_names, make_id(ds_chr, ds_pos, ds_comp_alt, ds_comp_ref))

  # Results table
  result <- data.frame(
    target_id  = target_df$variant_id,
    matched_id = NA_character_,
    flip       = FALSE,
    match_type = NA_character_,
    found      = FALSE,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(target_df))) {
    tv <- target_df$variant_id[i]

    if (!is.na(lookup_exact[tv])) {
      result$matched_id[i] <- lookup_exact[tv]
      result$flip[i]       <- FALSE
      result$match_type[i] <- "exact"
      result$found[i]      <- TRUE
    } else if (!is.na(lookup_swap[tv])) {
      result$matched_id[i] <- lookup_swap[tv]
      result$flip[i]       <- TRUE
      result$match_type[i] <- "ref_alt_swap"
      result$found[i]      <- TRUE
    } else if (!is.na(lookup_comp[tv])) {
      result$matched_id[i] <- lookup_comp[tv]
      result$flip[i]       <- FALSE
      result$match_type[i] <- "strand_complement"
      result$found[i]      <- TRUE
    } else if (!is.na(lookup_comp_swap[tv])) {
      result$matched_id[i] <- lookup_comp_swap[tv]
      result$flip[i]       <- TRUE
      result$match_type[i] <- "complement_and_swap"
      result$found[i]      <- TRUE
    }
  }

  result
}

# -----------------------------------------------------------------------------
# 1. Load reference files
# -----------------------------------------------------------------------------
cat("Loading reference files...\n")

ref_variants <- fread(file.path(opt$ref_dir, "mtPCA.variant_ids.txt"),
                      header = FALSE, col.names = "variant_id")$variant_id
mu_ref       <- fread(file.path(opt$ref_dir, "mtPCA.mu.txt"),
                      header = FALSE)$V1
loadings     <- as.matrix(fread(file.path(opt$ref_dir, "mtPCA.loadings.PC1-20.txt"),
                                header = FALSE))
eigenvalues  <- fread(file.path(opt$ref_dir, "mtPCA.eigenvalues.txt"),
                      header = FALSE)$V1

stopifnot(length(ref_variants) == length(mu_ref))
stopifnot(nrow(loadings) == length(ref_variants))
stopifnot(ncol(loadings) == 20)

cat(sprintf("  Reference variants : %d\n", length(ref_variants)))
cat(sprintf("  Loadings matrix    : %d x %d\n", nrow(loadings), ncol(loadings)))

# -----------------------------------------------------------------------------
# 2. Load external cohort genotypes
# -----------------------------------------------------------------------------
cat("\nLoading external cohort genotypes...\n")
ext <- read_plink(opt$ext_prefix)
cat(sprintf("  Samples : %d\n", nrow(ext$geno)))
cat(sprintf("  Variants: %d\n", ncol(ext$geno)))

# -----------------------------------------------------------------------------
# 3. Harmonise external variants to reference universe
# -----------------------------------------------------------------------------
cat("\nHarmonising external variants to reference universe...\n")
ext_match <- harmonise_variants(ref_variants, ext$bim)

cat(sprintf("  Exact matches         : %d\n", sum(ext_match$match_type == "exact",               na.rm = TRUE)))
cat(sprintf("  REF/ALT swapped       : %d\n", sum(ext_match$match_type == "ref_alt_swap",        na.rm = TRUE)))
cat(sprintf("  Strand complement     : %d\n", sum(ext_match$match_type == "strand_complement",   na.rm = TRUE)))
cat(sprintf("  Complement + swap     : %d\n", sum(ext_match$match_type == "complement_and_swap", na.rm = TRUE)))
cat(sprintf("  Not found             : %d\n", sum(!ext_match$found)))

# Check overlap threshold
n_nonambig <- nrow(ext_match)
n_found    <- sum(ext_match$found)
pct_found  <- round(100 * n_found / n_nonambig, 1)

cat(sprintf("\n  Non-ambiguous reference variants: %d\n", n_nonambig))
cat(sprintf("  Matched in external             : %d (%.1f%%)\n", n_found, pct_found))

if (pct_found < 95) {
  cat(sprintf(
    "\n[WARNING] Variant overlap is %.1f%%, below the recommended 95%% threshold.\n",
    pct_found))
  cat("  Projection into the pre-computed reference space is not recommended.\n")
  cat("  Please use mtPCA_recompute.R with the provided 1000G reference data instead.\n")
  stop("Stopping: insufficient variant overlap for reliable projection.")
  # warning("Variant overlap below 95% threshold. Proceeding with caution.")

}

# -----------------------------------------------------------------------------
# 4. Build aligned genotype matrix (samples x ref_variants)
#    Unmatched variants imputed to reference allele (0) before centring,
#    so they contribute zero after centring with mu
# -----------------------------------------------------------------------------
cat("\nBuilding aligned genotype matrix...\n")

n_samples <- nrow(ext$geno)
G <- matrix(0.0, nrow = n_samples, ncol = n_nonambig)  # default 0 = ref allele
rownames(G) <- rownames(ext$geno)

for (i in seq_len(n_nonambig)) {
  if (!ext_match$found[i]) next
  col <- ext_match$matched_id[i]
  g   <- as.numeric(ext$geno[, col])
  if (ext_match$flip[i]) g <- 1 - g
  # Impute any within-sample missingness to reference allele (0)
  g[is.na(g)] <- 0
  G[, i] <- g
}

n_missing <- sum(is.na(ext$geno[, ext_match$matched_id[ext_match$found]]))
if (n_missing > 0)
  cat(sprintf("  Imputed %d missing genotypes to reference allele (0)\n", n_missing))

# -----------------------------------------------------------------------------
# 5. Align mu to non-ambiguous variants
#    mu_ref is ordered to match ref_variants; after ambiguous SNP removal
#    we need the corresponding subset of mu values
# -----------------------------------------------------------------------------
ambig       <- is_ambiguous(
  do.call(rbind, strsplit(ref_variants, ":"))[, 3],
  do.call(rbind, strsplit(ref_variants, ":"))[, 4]
)
mu_aligned      <- mu_ref[!ambig]
loadings_aligned <- loadings[!ambig, ]

stopifnot(length(mu_aligned) == n_nonambig)
stopifnot(nrow(loadings_aligned) == n_nonambig)

# -----------------------------------------------------------------------------
# 6. Centre using reference means and project
#    Do NOT use external allele frequencies for centring —
#    must use reference mu to preserve alignment with pre-computed loadings
# -----------------------------------------------------------------------------
cat("\nCentring and projecting...\n")
G_centred  <- sweep(G, 2, mu_aligned, "-")
PC_scores  <- G_centred %*% loadings_aligned
colnames(PC_scores) <- paste0("PC", 1:20)

cat(sprintf("  Projection complete: %d samples x 20 PCs\n", nrow(PC_scores)))

# -----------------------------------------------------------------------------
# 7. Output
# -----------------------------------------------------------------------------
cat("\nWriting output files...\n")

out_df <- as.data.frame(PC_scores)
write.table(out_df,
            paste0(opt$out, ".mtPCs.txt"),
            quote = FALSE, row.names = TRUE)
cat(sprintf("  Projected PCs   : %s.mtPCs.txt\n", opt$out))

write.table(
  data.frame(
    total_ref_variants = length(ref_variants),
    non_ambiguous      = n_nonambig,
    matched            = n_found,
    pct_overlap        = pct_found
  ),
  paste0(opt$out, ".overlap_summary.txt"),
  quote = FALSE, row.names = FALSE, sep = "\t"
)
cat(sprintf("  Overlap summary : %s.overlap_summary.txt\n", opt$out))

cat("\nSession info:\n")
print(sessionInfo())
cat("\n=== Done ===\n")


