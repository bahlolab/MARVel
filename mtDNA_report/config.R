# ── paths ────────────────────────────────────────────────────────────────────
results_root <- "/vast/projects/bahlo_mtDNA/Results_ataxia"
report_dir   <- "/vast/projects/bahlo_mtDNA/Results_ataxia/mtDNA_report"
pca_dir      <- file.path(report_dir, "1000G")

# ID discovery: all subdirectories of results_root starting with "AA"
# Override by supplying a character vector, e.g.:
#   cohort_ids <- c("AA0001-01", "AA0002-01")
cohort_ids   <- NULL   # NULL = auto-discover

# ── thresholds ────────────────────────────────────────────────────────────────
min_median_depth  <- 100    # flag samples below this
contam_threshold  <- 0.02   # haplocheck contamination level
min_hom_vaf       <- 0.95   # VAF cut-off for homoplasmic calls
max_het_vaf       <- 0.95

# ── variant filtering ─────────────────────────────────────────────────────────
# Positions blacklisted as known artefact-prone sites
variant_blacklist <- c(301, 302, 310, 316, 3107, 5894, 10933, 16179,
                       16181:16183, 16188, 16189, 16192)
# FILTER tags that disqualify a variant (partial match)
filter_exclude    <- "strand_bias"

# ── labels ────────────────────────────────────────────────────────────────────
cohort_label <- "Ataxia cohort"