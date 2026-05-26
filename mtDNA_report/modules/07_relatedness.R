# Module 07 — relatedness / duplicate check (mtDNA variant-based)
#
# Expects in caller environment (from 01_load_data + config):
#   samples, variant_blacklist, filter_exclude, min_hom_vaf
#
# Approach: pairwise Jaccard similarity on PASS variants (AF >= 0.03).
# Mother-offspring pairs and duplicates share near-identical variant profiles.
# Haplogroup-defining homoplasmic variants are included — same-haplogroup
# unrelated samples will be similar but duplicates/relatives will be near 1.0.
#
# Thresholds (adjust as needed):
#   flag_threshold >= 1.00  probable duplicate / identical profile
#   rel_threshold  >= 0.95  probable mother-offspring / close relative
#
# Produces:
#   rel_plots$ibs_distribution  — histogram of all pairwise Jaccard values
#   rel_plots$flagged_pairs     — dot plot of flagged pairs (if any)
#   rel_data$pairs              — data.frame of flagged pairs

library(ggplot2)
library(plotly)

rel_plots <- list()
rel_data  <- list()

flag_threshold  <- 1.0    # Jaccard = 1: duplicate or identical profile
rel_threshold   <- 0.95   # Relative?: Jaccard > 0.95 (~5% variants differ)

# ── helpers ───────────────────────────────────────────────────────────────────
pass_var <- function(v) {
  if (is.null(v) || nrow(v) == 0) return(NULL)
  v <- v[!grepl(filter_exclude, v$FILTER) & !(v$POS %in% variant_blacklist), ]
  if (nrow(v) == 0) return(NULL)
  v
}

# ── build variant key sets per sample ────────────────────────────────────────
vcf_ok <- Filter(Negate(is.null), samples$vcf)   # re-uses name from 06_variants

var_keys <- lapply(names(vcf_ok), function(id) {
  v <- pass_var(vcf_ok[[id]])
  if (is.null(v)) return(character(0))
  v <- v[order(v$POS, -v$AF), ]
  v <- v[!duplicated(v$POS), ]                   # one allele per position
  v <- v[!is.na(v$AF) & v$AF >= min_hom_vaf, ]
  paste(v$POS, v$REF, v$ALT, sep = ":")
})
names(var_keys) <- names(vcf_ok)
var_keys <- var_keys[sapply(var_keys, length) > 0]

n   <- length(var_keys)
ids <- names(var_keys)

if (n < 2) {
  message("[07_relatedness] Fewer than 2 samples with variants — skipping.")
} else {
  
  # ── binary presence/absence matrix (samples x variants) ──────────────────
  all_keys <- unique(unlist(var_keys))
  M <- matrix(0L, nrow = n, ncol = length(all_keys),
              dimnames = list(ids, all_keys))
  for (i in seq_len(n))
    M[i, var_keys[[ids[i]]]] <- 1L
  
  # ── pairwise metrics ──────────────────────────────────────────────────────
  inter     <- tcrossprod(M)                       # n x n intersection counts
  rs        <- rowSums(M)
  union_mat <- outer(rs, rs, "+") - inter
  jaccard   <- inter / union_mat
  
  diag(jaccard) <- NA
  
  # ── extract upper triangle ────────────────────────────────────────────────
  idx      <- which(upper.tri(jaccard), arr.ind = TRUE)
  pairs_df <- data.frame(
    ID1     = ids[idx[, 1]],
    ID2     = ids[idx[, 2]],
    jaccard = jaccard[idx],
    n_inter = inter[idx],
    n_mean  = (rs[idx[, 1]] + rs[idx[, 2]]) / 2,
    stringsAsFactors = FALSE
  )
  pairs_df <- pairs_df[!is.na(pairs_df$jaccard), ]
  
  pairs_df$label <- ifelse(pairs_df$jaccard >= flag_threshold, "Duplicate/Relative",
                           ifelse(pairs_df$jaccard >  rel_threshold,  "Relative?",
                                  "Unrelated"))
  
  flagged        <- pairs_df[pairs_df$label != "Unrelated", ]
  flagged        <- flagged[order(-flagged$jaccard), ]
  rel_data$pairs <- flagged
  
  n_dup_rel <- sum(pairs_df$label == "Duplicate/Relative")
  n_rel     <- sum(pairs_df$label == "Relative?")
  
  # ── 1. pairwise Jaccard distribution ─────────────────────────────────────
  rel_plots$ibs_distribution <- ggplot(pairs_df, aes(x = jaccard)) +
    geom_histogram(bins = 100, fill = "steelblue", colour = "white") +
    geom_vline(xintercept = flag_threshold,
               linetype = "dashed", colour = "firebrick", linewidth = 0.7) +
    annotate("text", x = flag_threshold, y = Inf,
             label = "duplicate/relative\n(Jaccard=1)",
             hjust = 1.1, vjust = 1.5, size = 3, colour = "firebrick") +
    scale_y_log10() +
    labs(title = "Pairwise mtDNA variant similarity (Jaccard)",
         x = "Jaccard similarity", y = "# pairs (log10)",
         caption = sprintf("%s pairs  |  %d duplicate/relative (Jaccard=1)  |  %d relative? (Jaccard>%.2f)",
                           format(nrow(pairs_df), big.mark = ","),
                           n_dup_rel, n_rel, rel_threshold)) +
    theme_minimal(base_size = 11)
  
  # ── 2. flagged pairs ──────────────────────────────────────────────────────
  if (nrow(flagged) > 0) {
    fp <- flagged[flagged$jaccard > 0.95, ]
    fp$pair_label <- factor(
      paste(fp$ID1, fp$ID2, sep = " — "),
      levels = rev(paste(fp$ID1, fp$ID2, sep = " — ")))
    
    p_fp <- ggplot(fp,
                   aes(x = jaccard, y = pair_label, colour = label,
                       text = paste0(ID1, " — ", ID2,
                                     "<br>Jaccard: ", round(jaccard, 4),
                                     "<br>", label))) +
      geom_point(size = 3) +
      geom_segment(aes(x = 0.95, xend = jaccard, yend = pair_label),
                   linewidth = 0.4, alpha = 0.5) +
      scale_colour_manual(values = c("Duplicate/Relative" = "firebrick",
                                     "Relative?"           = "orange"),
                          name = NULL) +
      scale_x_continuous(limits = c(0.95, 1.01)) +
      labs(title = "Flagged sample pairs",
           x = "Jaccard similarity", y = NULL) +
      theme_minimal(base_size = 11) +
      theme(legend.position  = "top",
            axis.text.y      = element_text(size = 8))
    
    plot_height <- max(400, nrow(fp) * 28)
    
    rel_plots$flagged_pairs <- ggplotly(p_fp, tooltip = "text",
                                        height = plot_height) |>
      layout(margin = list(l = 150))
  } else {
    message("[07_relatedness] No pairs flagged above thresholds.")
  }
  
  message(sprintf(
    "[07_relatedness] %d samples  |  %s pairs  |  %d duplicate/relative  |  %d relative.",
    n, format(nrow(pairs_df), big.mark = ","), n_dup_rel, n_rel))
}

