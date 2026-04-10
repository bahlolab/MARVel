#!/usr/bin/env Rscript
# =============================================================================
# mtPCA_recompute.R
#
# Recompute mtDNA PCA from the intersection of high-quality variants shared
# between the 1000G reference panel and an external cohort.
#
# Use this script when variant overlap between external data and the reference
# is < 95%, making projection into the pre-computed space unreliable.
#
# Both 1000G and external cohort genotypes are read from PLINK BIM/BED/FAM.
#
# Variant harmonisation priority:
#   1. Exact match          (chr:pos:ref:alt identical)
#   2. REF/ALT swapped      (chr:pos:alt:ref) -> recode 0<->1
#   3. Strand complement    (chr:pos:comp(ref):comp(alt)) -> use as-is
#   4. Complement + swap    (chr:pos:comp(alt):comp(ref)) -> recode 0<->1
#   Ambiguous SNPs (A/T or C/G) are dropped with a warning.
#
# Usage:
#   Rscript mtPCA_recompute.R \
#     --ref_prefix <path to 1000G plink prefix> \
#     --ext_prefix <path to external cohort plink prefix> \
#     --out        <output prefix> \
#     [--n_pcs     <number of PCs, default 20>] \
#     [--ref_dir   <directory containing mtPCA.variant_ids.txt, default '.'>] \
#     [--save_loadings]
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
  make_option("--ref_prefix",    type = "character",
              help = "Path to 1000G reference PLINK prefix (bim/bed/fam)"),
  make_option("--ext_prefix",    type = "character",
              help = "Path to external cohort PLINK prefix (bim/bed/fam)"),
  make_option("--ref_dir",       type = "character", default = ".",
              help = "Directory containing mtPCA.variant_ids.txt [default: .]"),
  make_option("--n_pcs",         type = "integer",   default = 20L,
              help = "Number of PCs to compute [default: 20]"),
  make_option("--out",           type = "character", default = "mtPCA_recomputed",
              help = "Output file prefix [default: mtPCA_recomputed]"),
  make_option("--save_loadings", action = "store_true", default = FALSE,
              help = "Save recomputed loadings, mu, and variant IDs for future projection")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$ref_prefix)) stop("--ref_prefix is required. Run with --help for usage.")
if (is.null(opt$ext_prefix)) stop("--ext_prefix is required. Run with --help for usage.")

cat("=== mtPCA Recomputation Workflow ===\n\n")

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


  # Name columns by dataset's own chr:pos:ref:alt
  colnames(geno) <- paste(bim$chrom, bim$pos, bim$ref, bim$alt, sep = ":")

  list(geno = geno, bim = bim, fam = fam)
}

# -----------------------------------------------------------------------------
# Helper: harmonise a dataset's variants to a set of target variant IDs
#
# For each target variant (chr:pos:ref:alt), attempts matching in order:
#   1. Exact:           chr:pos:ref:alt
#   2. REF/ALT swap:    chr:pos:alt:ref           -> flip coding
#   3. Complement:      chr:pos:comp(ref):comp(alt)
#   4. Comp + swap:     chr:pos:comp(alt):comp(ref) -> flip coding
#
# Ambiguous SNPs (A/T or C/G) in the TARGET are dropped before matching.
# Returns data.frame with one row per (non-ambiguous) target variant.
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
# 1. Load reference variant universe
# -----------------------------------------------------------------------------
cat("Loading reference variant list...\n")
ref_variants <- fread(file.path(opt$ref_dir, "mtPCA.variant_ids.txt"),
                      header = FALSE, col.names = "variant_id")$variant_id
cat(sprintf("  High-quality reference variants: %d\n", length(ref_variants)))

# -----------------------------------------------------------------------------
# 2. Load genotype data
# -----------------------------------------------------------------------------
cat("\nLoading 1000G reference genotypes...\n")
ref <- read_plink(opt$ref_prefix)
cat(sprintf("  Samples : %d\n", nrow(ref$geno)))
cat(sprintf("  Variants: %d\n", ncol(ref$geno)))

cat("\nLoading external cohort genotypes...\n")
ext <- read_plink(opt$ext_prefix)
cat(sprintf("  Samples : %d\n", nrow(ext$geno)))
cat(sprintf("  Variants: %d\n", ncol(ext$geno)))

# -----------------------------------------------------------------------------
# 3. Harmonise both datasets to the reference variant universe
# -----------------------------------------------------------------------------
cat("\nHarmonising 1000G variants to reference universe...\n")
ref_match <- harmonise_variants(ref_variants, ref$bim)
cat(sprintf("  Exact matches         : %d\n", sum(ref_match$match_type == "exact",          na.rm = TRUE)))
cat(sprintf("  REF/ALT swapped       : %d\n", sum(ref_match$match_type == "ref_alt_swap",   na.rm = TRUE)))
cat(sprintf("  Strand complement     : %d\n", sum(ref_match$match_type == "strand_complement",  na.rm = TRUE)))
cat(sprintf("  Complement + swap     : %d\n", sum(ref_match$match_type == "complement_and_swap", na.rm = TRUE)))
cat(sprintf("  Not found             : %d\n", sum(!ref_match$found)))

cat("\nHarmonising external variants to reference universe...\n")
ext_match <- harmonise_variants(ref_variants, ext$bim)
cat(sprintf("  Exact matches         : %d\n", sum(ext_match$match_type == "exact",          na.rm = TRUE)))
cat(sprintf("  REF/ALT swapped       : %d\n", sum(ext_match$match_type == "ref_alt_swap",   na.rm = TRUE)))
cat(sprintf("  Strand complement     : %d\n", sum(ext_match$match_type == "strand_complement",  na.rm = TRUE)))
cat(sprintf("  Complement + swap     : %d\n", sum(ext_match$match_type == "complement_and_swap", na.rm = TRUE)))
cat(sprintf("  Not found             : %d\n", sum(!ext_match$found)))

# -----------------------------------------------------------------------------
# 4. Find shared intersection
# -----------------------------------------------------------------------------
cat("\nFinding shared variant intersection...\n")

shared_idx      <- which(ref_match$found & ext_match$found)
shared_variants <- ref_match$target_id[shared_idx]
n_shared        <- length(shared_variants)
# Denominator uses non-ambiguous variants (ambiguous already dropped in harmonise_variants)
n_nonambig      <- nrow(ref_match)
pct_shared      <- round(100 * n_shared / n_nonambig, 1)

cat(sprintf("  Non-ambiguous reference variants: %d\n", n_nonambig))
cat(sprintf("  Shared (intersection)           : %d / %d (%.1f%%)\n",
            n_shared, n_nonambig, pct_shared))

if (pct_shared < 95) {
  cat(sprintf(
    "\n[WARNING] Variant overlap is %.1f%%, below the recommended 95%% threshold.\n",
    pct_shared))
  cat("  PCA recomputation will proceed on the available intersection,\n")
  cat("  but results should be interpreted with caution.\n\n")
}
if (n_shared < 10)
  stop("Fewer than 10 variants in the intersection — cannot compute meaningful PCA.")

# -----------------------------------------------------------------------------
# 5. Build aligned genotype matrices (samples x shared_variants)
# -----------------------------------------------------------------------------
cat("\nBuilding aligned genotype matrices...\n")

build_aligned <- function(geno, match_df, shared_idx) {
  G <- matrix(NA_real_, nrow = nrow(geno), ncol = length(shared_idx))
  rownames(G) <- rownames(geno)
  for (i in seq_along(shared_idx)) {
    mi  <- shared_idx[i]
    col <- match_df$matched_id[mi]
    g   <- as.numeric(geno[, col])
    if (match_df$flip[mi]) g <- 1 - g
    G[, i] <- g
  }
  G
}

G_ref <- build_aligned(ref$geno, ref_match, shared_idx)
G_ext <- build_aligned(ext$geno, ext_match, shared_idx)

# Impute missing to reference allele (0)
impute_missing <- function(G, label) {
  n_miss <- sum(is.na(G))
  if (n_miss > 0) {
    cat(sprintf("  Imputing %d missing genotypes in %s to reference allele (0)\n",
                n_miss, label))
    G[is.na(G)] <- 0
  }
  G
}
G_ref <- impute_missing(G_ref, "1000G")
G_ext <- impute_missing(G_ext, "external")

# -----------------------------------------------------------------------------
# 6. Centre using 1000G means, run SVD, project external samples
#    Convention consistent with original PCA script:
#      mu      = colMeans(G_ref)
#      scores  = U %*% diag(d)
#      eigvals = d^2 / (n - 1)
#      project = G_ext_centred %*% V
# -----------------------------------------------------------------------------
cat("\nCentring using 1000G means...\n")
mu            <- colMeans(G_ref)
G_ref_centred <- sweep(G_ref, 2, mu, "-")

cat("Running SVD on 1000G reference...\n")
sv <- svd(G_ref_centred)
k  <- min(opt$n_pcs, length(sv$d))

ref_scores <- sv$u[, 1:k, drop = FALSE] %*% diag(sv$d[1:k], nrow = k, ncol = k)
loadings   <- sv$v[, 1:k, drop = FALSE]
eigvals    <- (sv$d^2) / (nrow(G_ref_centred) - 1)

colnames(ref_scores) <- paste0("PC", 1:k)
rownames(ref_scores) <- rownames(G_ref)
cat(sprintf("  SVD complete: %d reference samples x %d PCs\n", nrow(ref_scores), k))

cat("Projecting external samples...\n")
G_ext_centred        <- sweep(G_ext, 2, mu, "-")
ext_scores           <- G_ext_centred %*% loadings
colnames(ext_scores) <- paste0("PC", 1:k)
rownames(ext_scores) <- rownames(G_ext)
cat(sprintf("  Projection complete: %d external samples x %d PCs\n", nrow(ext_scores), k))

# -----------------------------------------------------------------------------
# 7. Write output files
# -----------------------------------------------------------------------------
cat("\nWriting output files...\n")

write.table(as.data.frame(ext_scores),
            sprintf("%s.scores.external.PC1-%d.txt", opt$out, k),
            quote = FALSE, row.names = TRUE)
cat(sprintf("  External PCs    : %s.scores.external.PC1-%d.txt\n", opt$out, k))

write.table(as.data.frame(ref_scores),
            sprintf("%s.scores.1000G.PC1-%d.txt", opt$out, k),
            quote = FALSE, row.names = TRUE)
cat(sprintf("  1000G PCs       : %s.scores.1000G.PC1-%d.txt\n", opt$out, k))

write.table(eigvals[1:k],
            paste0(opt$out, ".eigenvalues.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
cat(sprintf("  Eigenvalues     : %s.eigenvalues.txt\n", opt$out))

write.table(
  data.frame(
    total_ref_variants    = length(ref_variants),
    non_ambiguous         = n_nonambig,
    found_in_1000G        = sum(ref_match$found),
    found_in_external     = sum(ext_match$found),
    shared_variants       = n_shared,
    pct_overlap           = pct_shared
  ),
  paste0(opt$out, ".overlap_summary.txt"),
  quote = FALSE, row.names = FALSE, sep = "\t"
)
cat(sprintf("  Overlap summary : %s.overlap_summary.txt\n", opt$out))

if (opt$save_loadings) {
  write.table(loadings,
              sprintf("%s.loadings.PC1-%d.txt", opt$out, k),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(mu,
              paste0(opt$out, ".mu.txt"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(shared_variants,
              paste0(opt$out, ".variant_ids.txt"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  cat(sprintf("  Loadings        : %s.loadings.PC1-%d.txt\n", opt$out, k))
  cat(sprintf("  Means (mu)      : %s.mu.txt\n", opt$out))
  cat(sprintf("  Variant IDs     : %s.variant_ids.txt\n", opt$out))
}

cat("\nSession info:\n")
print(sessionInfo())
cat("\n=== Done ===\n")
