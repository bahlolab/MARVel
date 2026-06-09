# Module 01 — discover samples and load all mitoHPC outputs
#
# Expects in caller environment:
#   results_root  — path to the results directory
#   cohort_ids    — character vector of IDs, or NULL to auto-discover
#   min_median_depth, contam_threshold, min_hom_vaf (from config.R)
#
# Produces:
#   samples$meta  — data.frame, one row per sample, all scalar QC fields
#   samples$cvg   — named list of per-position depth data.frames (chrM, pos, depth)
#   samples$vcf   — named list of variant data.frames (one row per ALT allele)

# ── helpers ───────────────────────────────────────────────────────────────────

# Strip trailing suffixes to get the inner directory name used in mitoHPC output paths.
# Edit the regex if your naming convention differs. Examples:
#   "SAMPLE001-vcgs"    -> "SAMPLE001"
#   "SAMPLE001-batch2"  -> "SAMPLE001"
#   "SAMPLE001"         -> "SAMPLE001"  (no-op)
inner_id <- function(id) sub("-(vcgs|batch[0-9]*)$", "", id, ignore.case = TRUE)

# Build the path to a mitoHPC per-sample output file.
sample_path <- function(id, ext) {
  iid <- inner_id(id)
  file.path(results_root, id, "mitoHPC", "out",
            paste0(iid, ".merged"),
            paste0(iid, ".merged.", ext))
}

# Build the path to a mitoHPC cohort-level tab file.
tab_path <- function(id, name)
  file.path(results_root, id, "mitoHPC", "out", name)

# Read a TSV safely; return NULL with a warning on any failure.
read_tsv_safe <- function(path, quote = "\"", check.names = TRUE, ...) {
  if (!file.exists(path)) return(NULL)
  tryCatch(
    read.table(path, header = TRUE, sep = "\t",
               quote = quote, comment.char = "",
               stringsAsFactors = FALSE, check.names = check.names,
               na.strings = c("", "NA"), ...),
    error = function(e) { warning(path, ": ", e$message, call. = FALSE); NULL }
  )
}

# ── INFO parser ───────────────────────────────────────────────────────────────
parse_info <- function(info_str) {
  fields <- strsplit(info_str, ";")[[1]]
  out    <- list()
  for (f in fields) {
    if (grepl("=", f)) {
      kv           <- strsplit(f, "=")[[1]]
      out[[kv[1]]] <- paste(kv[-1], collapse = "=")
    } else {
      out[[f]] <- TRUE
    }
  }
  out
}

# ── VCF parser ────────────────────────────────────────────────────────────────
# Returns a data.frame with one row per ALT allele (multi-allelics expanded).
# Columns: CHROM POS REF ALT FILTER GT AF DP + all INFO fields
parse_vcf <- function(path) {
  if (!file.exists(path)) return(NULL)
  lines <- tryCatch(readLines(path), error = function(e) NULL)
  if (is.null(lines)) return(NULL)

  data_lines <- lines[!startsWith(lines, "##")]
  if (length(data_lines) < 2) return(NULL)

  con <- textConnection(data_lines)
  vcf <- tryCatch(
    read.table(con, header = TRUE, sep = "\t",
               quote = "", comment.char = "",
               stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL,
    finally = close(con)
  )
  if (is.null(vcf) || nrow(vcf) == 0) return(NULL)

  names(vcf)[1] <- "CHROM"
  sample_col    <- names(vcf)[10]

  info_list <- lapply(vcf$INFO, parse_info)
  all_keys  <- unique(unlist(lapply(info_list, names)))

  rows <- lapply(seq_len(nrow(vcf)), function(i) {
    alts <- strsplit(vcf$ALT[i], ",")[[1]]
    fmt  <- strsplit(vcf$FORMAT[i], ":")[[1]]
    vals <- strsplit(vcf[[sample_col]][i], ":")[[1]]
    fv   <- setNames(vals, fmt)

    afs <- if ("AF" %in% fmt)
             as.numeric(strsplit(fv["AF"], ",")[[1]])
           else rep(NA_real_, length(alts))
    dp  <- if ("DP" %in% fmt) as.integer(fv["DP"]) else NA_integer_
    gt  <- if ("GT" %in% fmt) fv["GT"]             else NA_character_

    info    <- info_list[[i]]
    info_df <- as.data.frame(
      lapply(setNames(all_keys, all_keys), function(k) {
        v <- info[[k]]; if (is.null(v)) NA else v
      }),
      stringsAsFactors = FALSE, check.names = FALSE
    )

    base_df  <- data.frame(
      CHROM  = vcf$CHROM[i],
      POS    = vcf$POS[i],
      REF    = vcf$REF[i],
      ALT    = alts,
      FILTER = vcf$FILTER[i],
      GT     = gt,
      AF     = afs,
      DP     = dp,
      stringsAsFactors = FALSE,
      row.names = NULL
    )
    info_row        <- info_df[rep(1, length(alts)), , drop = FALSE]
    info_row$SM     <- NULL
    info_row$NONSYN <- NULL
    info_row$HG     <- NULL
    info_row$NUMT   <- !is.na(info_row$NUMT)
    cbind(base_df, info_row, row.names = NULL)
  })
  do.call(rbind, rows)
}

# ── clade lookup (defined once) ───────────────────────────────────────────────
m_clade <- c("M", "C", "D", "G", "E", "Q", "Z")

get_macrohap <- function(hg) {
  if (is.na(hg)) return(NA_character_)
  if (grepl("^L", hg)) {
    m <- regmatches(hg, regexpr("^L[0-9]", hg))
    if (length(m)) m else "L"          # fallback for bare "L"
  } else {
    regmatches(hg, regexpr("^[A-Z]+", hg))[1]
  }
}

get_clade <- function(macrohap) {
  if (is.na(macrohap))              return(NA_character_)
  if (grepl("^L", macrohap))        return("L")
  if (macrohap %in% m_clade)        return("M")
  "N/R"
}

# ── ID discovery ──────────────────────────────────────────────────────────────
if (is.null(cohort_ids)) {
  dirs <- list.dirs(results_root, recursive = FALSE, full.names = FALSE)
  cohort_ids <- if (!is.null(cohort_id_pattern)) dirs[grepl(cohort_id_pattern, dirs)] else dirs
}
message(sprintf("[01_load_data] %d samples to load.", length(cohort_ids)))

# ── per-sample loading ────────────────────────────────────────────────────────
meta_list <- vector("list", length(cohort_ids))
cvg_list  <- vector("list", length(cohort_ids))
vcf_list  <- vector("list", length(cohort_ids))

for (k in seq_along(cohort_ids)) {
  id <- cohort_ids[k]
  if (k %% 50 == 0 || k == length(cohort_ids))
    message(sprintf("[01_load_data]   %d / %d", k, length(cohort_ids)))

  # -- coverage stats ----------------------------------------------------------
  cvg_stat <- read_tsv_safe(sample_path(id, "mutect2.cvg.stat"))
  if (is.null(cvg_stat))
    cvg_stat <- read_tsv_safe(sample_path(id, "cvg.stat"))

  depth_median <- if (!is.null(cvg_stat)) cvg_stat$median[1] else NA_real_
  depth_mean   <- if (!is.null(cvg_stat)) cvg_stat$mean[1]   else NA_real_
  depth_min    <- if (!is.null(cvg_stat)) cvg_stat$min[1]    else NA_real_
  depth_max    <- if (!is.null(cvg_stat)) cvg_stat$max[1]    else NA_real_

  # -- read counts + mtDNA-CN --------------------------------------------------
  cnt      <- read_tsv_safe(tab_path(id, "count.tab"), check.names = FALSE)
  all_reads <- if (!is.null(cnt)) cnt[["all_reads"]][1] else NA_integer_
  mt_reads  <- if (!is.null(cnt)) cnt[["MT_reads"]][1]  else NA_integer_
  mtdna_cn  <- if (!is.null(cnt)) cnt[["mtDNA-CN"]][1]  else NA_real_

  # -- haplogroup --------------------------------------------------------------
  hg         <- read_tsv_safe(sample_path(id, "mutect2.haplogroup"))
  haplogroup <- if (!is.null(hg)) hg$Haplogroup[1]          else NA_character_
  hg_quality <- if (!is.null(hg)) as.numeric(hg$Quality[1]) else NA_real_
  macrohap   <- get_macrohap(haplogroup)
  clade      <- get_clade(macrohap)

  # -- haplocheck (contamination) ----------------------------------------------
  hc            <- read_tsv_safe(sample_path(id, "mutect2.haplocheck"))
  contam_status <- if (!is.null(hc)) hc[["Contamination.Status"]][1] else NA_character_
  contam_level  <- if (!is.null(hc)) {
                     raw <- hc[["Contamination.Level"]][1]
                     if (!is.na(raw) && raw == "ND") 0 else as.numeric(raw)
                   } else NA_real_

  # -- VCF ---------------------------------------------------------------------
  v <- parse_vcf(sample_path(id, "mutect2.00.vcf"))
  vcf_list[k] <- list(v)

  n_variants <- if (!is.null(v)) nrow(v) else NA_integer_
  if (!is.null(v)) {
    pass   <- !grepl(filter_exclude, v$FILTER) & !(v$POS %in% variant_blacklist) &
              !is.na(v$AF) & v$AF >= 0.03
    n_pass <- sum(pass,                        na.rm = TRUE)
    n_hom  <- sum(pass & v$AF >= min_hom_vaf, na.rm = TRUE)
    n_het  <- sum(pass & v$AF <  min_hom_vaf, na.rm = TRUE)
  } else {
    n_pass <- n_hom <- n_het <- NA_integer_
  }

  # -- per-position coverage ---------------------------------------------------
  cvg_path <- sample_path(id, "mutect2.cvg")
  if (!file.exists(cvg_path)) cvg_path <- sample_path(id, "cvg")
  if (file.exists(cvg_path)) {
    cvg_df <- tryCatch(
      read.table(cvg_path, header = FALSE,
                 col.names = c("chrom", "pos", "depth")),
      error = function(e) NULL
    )
    cvg_list[k] <- list(cvg_df)
  }

  # -- QC flag -----------------------------------------------------------------
  qc_flag <- (!is.na(depth_median) && depth_median < min_median_depth) |
             (!is.na(contam_level) && contam_level  >= contam_threshold)

  meta_list[[k]] <- data.frame(
    ID            = id,
    inner_id      = inner_id(id),
    haplogroup    = haplogroup,
    hg_quality    = hg_quality,
    macrohap      = macrohap,
    clade         = clade,
    depth_median  = depth_median,
    depth_mean    = depth_mean,
    depth_min     = depth_min,
    depth_max     = depth_max,
    all_reads     = all_reads,
    mt_reads      = mt_reads,
    mtdna_cn      = mtdna_cn,
    contam_status = contam_status,
    contam_level  = contam_level,
    n_variants    = n_variants,
    n_pass        = n_pass,
    n_hom         = n_hom,
    n_het         = n_het,
    qc_flag       = qc_flag,
    stringsAsFactors = FALSE
  )
}

names(cvg_list) <- cohort_ids
names(vcf_list) <- cohort_ids

samples <- list(
  meta = do.call(rbind, meta_list),
  cvg  = cvg_list,
  vcf  = vcf_list
)

n_ok  <- sum(!samples$meta$qc_flag, na.rm = TRUE)
n_bad <- sum( samples$meta$qc_flag, na.rm = TRUE)
message(sprintf("[01_load_data] Loaded %d samples (%d pass QC, %d flagged).",
                nrow(samples$meta), n_ok, n_bad))

# ── incomplete job detection ───────────────────────────────────────────────────
# A sample is considered incomplete if all key output files are missing.
incomplete <- samples$meta[
  is.na(samples$meta$depth_median) &
  is.na(samples$meta$haplogroup)   &
  is.na(samples$meta$n_variants),
  "ID"
]

partial <- samples$meta[
  !is.na(samples$meta$depth_median) != !is.na(samples$meta$n_variants),
  "ID"
]

if (length(incomplete) > 0) {
  warning(sprintf(
    "\n[01_load_data] %d sample(s) appear to have NO mitoHPC output (job likely did not complete):\n%s\n\nRe-run the corresponding SLURM jobs, e.g.:\n%s",
    length(incomplete),
    paste0("  ", incomplete, collapse = "\n"),
    paste0("  sbatch ", file.path(sh_dir_hint <- Sys.getenv("SH_DIR", "<sh_scripts_dir>"),
                                   paste0(incomplete, ".sh")), collapse = "\n")
  ), call. = FALSE)
}

if (length(partial) > 0) {
  message(sprintf(
    "[01_load_data] WARNING: %d sample(s) have partial output (job may have been interrupted):\n%s\nCheck output logs and consider re-running these jobs.",
    length(partial),
    paste0("  ", partial, collapse = "\n")
  ))
}
