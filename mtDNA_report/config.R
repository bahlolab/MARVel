# ── paths ────────────────────────────────────────────────────────────────────
results_root <- "/path/to/mitoHPC/results"   # parent directory of per-sample folders
report_dir   <- "/path/to/output/directory"  # where the report and qc_summary.tsv are written
pca_dir      <- file.path(getwd(), "1000G")
# Assumes working directory is set to the mtDNA_report/ folder (see README).
# Override if needed: pca_dir <- "/path/to/mtDNA_report/1000G"

# ── sample discovery ──────────────────────────────────────────────────────────
# NULL = use all immediate subdirectories of results_root as sample IDs.
# Supply a character vector to specify samples explicitly:
#   cohort_ids <- c("SAMPLE001", "SAMPLE002")
cohort_ids <- NULL

# Optional regex to filter subdirectory names during auto-discovery.
# NULL = no filter (use all subdirs). Example: "^AA" matches IDs starting with "AA".
cohort_id_pattern <- NULL

# ── thresholds ────────────────────────────────────────────────────────────────
min_median_depth  <- 100    # flag samples below this median depth
contam_threshold  <- 0.02   # haplocheck contamination level threshold
min_hom_vaf       <- 0.95   # VAF cut-off for homoplasmic classification

# ── variant filtering ─────────────────────────────────────────────────────────
# Positions blacklisted as known artefact-prone sites
variant_blacklist <- c(301, 302, 310, 316, 3107, 5894, 16179,
                       16181:16183, 16188, 16189, 16192)
# FILTER tags that disqualify a variant (partial string match)
filter_exclude    <- "strand_bias"

# ── labels ────────────────────────────────────────────────────────────────────
cohort_label <- "My cohort"
