# Module 08 — QC summary table
#
# Expects in caller environment (from 01_load_data + config):
#   samples                — list with $meta
#   min_median_depth, contam_threshold, min_hom_vaf
#
# Optionally uses rel_data$pairs from 07_relatedness if available.
#
# Produces:
#   qc_table       — DT interactive table (HTML report)
#   qc_plots$flags — stacked bar of QC flag counts per category

library(ggplot2)
library(DT)

meta <- samples$meta

# ── build QC summary ──────────────────────────────────────────────────────────
qc <- data.frame(
  ID              = meta$ID,
  Haplogroup      = meta$haplogroup,
  HG_quality      = round(meta$hg_quality, 3),
  Clade           = meta$clade,
  Depth_median    = round(meta$depth_median, 0),
  Depth_mean      = round(meta$depth_mean,   0),
  mtDNA_CN        = round(meta$mtdna_cn,     1),
  Contam_level    = round(meta$contam_level, 4),
  Contam_status   = meta$contam_status,
  N_hom           = meta$n_hom,
  N_het           = meta$n_het,
  stringsAsFactors = FALSE
)

# ── QC flags (use raw meta values, not rounded display values) ────────────────
qc$flag_depth   <- !is.na(meta$depth_median) & meta$depth_median < min_median_depth
qc$flag_contam  <- !is.na(meta$contam_level) & meta$contam_level >= contam_threshold
qc$flag_hg      <- is.na(qc$Haplogroup)

# relatedness tag (Jaccard = 1 only) — informational, does not affect QC_pass
if (exists("rel_data") && !is.null(rel_data$pairs) && nrow(rel_data$pairs) > 0) {
  dup_pairs   <- rel_data$pairs[rel_data$pairs$jaccard >= 1.0, ]
  flagged_ids <- unique(c(dup_pairs$ID1, dup_pairs$ID2))
  qc$related  <- qc$ID %in% flagged_ids
} else {
  qc$related <- FALSE
}

# CV from module 02 (if available)
if (exists("cv_df") && !is.null(cv_df)) {
  qc <- merge(qc, cv_df[, c("ID", "CV")], by = "ID", all.x = TRUE)
  qc$CV <- round(qc$CV, 3)
} else {
  qc$CV <- NA_real_
}

# MSS from module 06 (if available)
if (exists("mss_df") && !is.null(mss_df)) {
  qc <- merge(qc, mss_df[, c("ID", "MSS")], by = "ID", all.x = TRUE)
  qc$MSS <- round(qc$MSS, 3)
} else {
  qc$MSS <- NA_real_
}

qc$N_flags <- rowSums(qc[, c("flag_depth", "flag_contam", "flag_hg")])
qc$QC_pass <- qc$N_flags == 0

# ── 1. QC flag summary bar chart ─────────────────────────────────────────────
flag_counts <- data.frame(
  Flag  = c("Low depth", "Contamination", "No haplogroup"),
  N     = c(sum(qc$flag_depth), sum(qc$flag_contam), sum(qc$flag_hg)),
  stringsAsFactors = FALSE
)
flag_counts$Flag <- factor(flag_counts$Flag,
                           levels = flag_counts$Flag[order(-flag_counts$N)])

qc_plots <- list()

qc_plots$flags <- ggplot(flag_counts, aes(x = Flag, y = N)) +
  geom_col(fill = "firebrick", width = 0.5) +
  geom_text(aes(label = N), vjust = -0.4, size = 3.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = sprintf("QC flags  |  %d / %d samples pass all QC",
                       sum(qc$QC_pass), nrow(qc)),
       x = NULL, y = "# samples flagged") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.x = element_blank())

# ── 2. interactive DT table ───────────────────────────────────────────────────
# colour-code flag columns; format numeric columns
display_cols <- c("ID", "Haplogroup", "HG_quality", "Clade",
                  "Depth_median", "Depth_mean", "mtDNA_CN",
                  "Contam_level", "Contam_status", "N_hom", "N_het", "CV", "MSS",
                  "QC_pass", "related",
                  "flag_depth", "flag_contam", "flag_hg")

qc_dt <- qc[, display_cols]

qc_table <- datatable(
  qc_dt,
  rownames  = FALSE,
  filter    = "top",
  extensions = "Buttons",
  options   = list(
    pageLength = 25,
    dom        = "Bfrtip",
    buttons    = c("csv", "excel"),
    scrollX    = TRUE,
    columnDefs = list(
      list(visible = FALSE,
           targets  = which(display_cols %in%
                              c("flag_depth","flag_contam","flag_hg")) - 1)
    )
  ),
  caption = htmltools::tags$caption(
    style = "caption-side: top; font-weight: bold;",
    sprintf("QC summary — %d samples  |  %d pass  |  %d flagged",
            nrow(qc), sum(qc$QC_pass), sum(!qc$QC_pass))
  )
) |>
  formatRound(c("HG_quality", "Contam_level", "CV", "MSS"), digits = 3) |>
  formatRound(c("Depth_median", "Depth_mean", "mtDNA_CN"), digits = 0) |>
  formatStyle("QC_pass",
              backgroundColor = styleEqual(c(TRUE, FALSE),
                                           c("#C8E6C9", "#FFCDD2"))) |>
  formatStyle("related",
              backgroundColor = styleEqual(c(TRUE, FALSE),
                                           c("#FFF9C4", "transparent"))) |>
  formatStyle("Depth_median",
              backgroundColor = styleInterval(min_median_depth,
                                              c("#FFCDD2", "transparent"))) |>
  formatStyle("Contam_level",
              backgroundColor = styleInterval(contam_threshold,
                                              c("transparent", "#FFCDD2")))

# ── export QC table ───────────────────────────────────────────────────────────
qc_out_path <- file.path(report_dir, "qc_summary.tsv")
write.table(qc, qc_out_path, sep = "\t", row.names = FALSE, quote = FALSE)
message(sprintf("[08_qc_summary] QC table written to %s", qc_out_path))

message(sprintf("[08_qc_summary] %d samples  |  %d pass  |  %d low depth  |  %d contaminated  |  %d no haplogroup  |  %d related (tag only).",
                nrow(qc), sum(qc$QC_pass),
                sum(qc$flag_depth), sum(qc$flag_contam),
                sum(qc$flag_hg),    sum(qc$related)))