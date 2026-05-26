# Module 02 — coverage QC
#
# Expects in caller environment (from 01_load_data + config):
#   samples           — list with $meta and $cvg
#   min_median_depth  — depth flag threshold
#
# Produces:
#   cov_plots$depth_per_sample    — dot plot of median depth per sample
#   cov_plots$mtdna_cn            — histogram of mtDNA copy number
#   cov_plots$base_coverage       — per-position depth (all samples + cohort mean)
#   cov_plots$cumulative_coverage — % bases >= X depth per sample
#   cov_plots$cv_distribution     — histogram of per-sample coverage CV

library(ggplot2)
library(plotly)

cov_plots    <- list()
cv_threshold <- 0.2    # CV > this flagged as high in base_coverage plot

meta <- samples$meta

# ── 1. median depth per sample (interactive plotly) ───────────────────────────
meta_sorted        <- meta[order(meta$depth_median), ]
meta_sorted$ID_ord <- factor(meta_sorted$ID, levels = rev(meta_sorted$ID))
meta_sorted$colour <- ifelse(meta_sorted$qc_flag, "Flagged", "Pass")

n_flagged_depth <- sum(meta$qc_flag, na.rm = TRUE)

p <- ggplot(meta_sorted,
            aes(x = depth_median, y = ID_ord, colour = colour,
                text = paste0(ID, "<br>Depth: ", round(depth_median), "x"))) +
  geom_point(size = 2) +
  geom_vline(xintercept = min_median_depth, linetype = "dashed",
             colour = "firebrick", linewidth = 0.5) +
  scale_colour_manual(values = c("Pass" = "steelblue", "Flagged" = "firebrick"),
                      name = NULL) +
  scale_y_discrete(expand = expansion(add = 8)) +
  labs(title = sprintf("Median depth per sample  |  threshold = %d  |  %d / %d flagged",
                       min_median_depth, n_flagged_depth, nrow(meta)),
       x = "Median depth (post-Mutect2)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.y        = element_blank(),
        axis.ticks.y       = element_blank(),
        panel.grid.major.y = element_blank())

cov_plots$depth_per_sample <- ggplotly(p, tooltip = "text", height = 500) |>
  layout(margin = list(l = 20, t = 60, b = 60))

# ── 2. mtDNA copy number distribution ─────────────────────────────────────────
cn_df   <- meta[!is.na(meta$mtdna_cn), c("ID", "mtdna_cn")]
cn_vals <- cn_df$mtdna_cn

set.seed(42)
cn_df$jitter_y <- runif(nrow(cn_df), -0.4, 0.4)

p <- ggplot() +
  geom_histogram(data = cn_df, aes(x = mtdna_cn),
                 binwidth = 10, fill = "steelblue", colour = "white") +
  geom_point(data = cn_df,
             aes(x = mtdna_cn, y = jitter_y, text = paste0(ID, "<br>CN: ", round(mtdna_cn, 1))),
             colour = "steelblue", alpha = 0, size = 3) +
  geom_vline(xintercept = median(cn_vals), linetype = "dashed",
             colour = "firebrick", linewidth = 0.7) +
  geom_vline(xintercept = mean(cn_vals),   linetype = "dotted",
             colour = "firebrick", linewidth = 0.7) +
  labs(title = sprintf("mtDNA copy number  |  median = %.0f  |  mean = %.0f  |  range = %.0f–%.0f",
                       median(cn_vals), mean(cn_vals), min(cn_vals), max(cn_vals)),
       x = "mtDNA-CN", y = "# samples") +
  theme_minimal(base_size = 11)

cov_plots$mtdna_cn <- ggplotly(p, tooltip = "text")

# ── per-position coverage (plots 3–5) ─────────────────────────────────────────
cvg_ok <- Filter(Negate(is.null), samples$cvg)

if (length(cvg_ok) > 0) {
  
  # per-sample CV
  cv_vals <- sapply(names(cvg_ok), function(id) {
    d <- cvg_ok[[id]]$depth
    if (length(d) == 0 || mean(d) == 0) return(NA_real_)
    sd(d) / mean(d)
  })
  cv_df <- data.frame(ID = names(cv_vals), CV = cv_vals,
                      stringsAsFactors = FALSE)
  
  high_cv_ids <- cv_df$ID[!is.na(cv_df$CV) & cv_df$CV > cv_threshold]
  n_high_cv   <- length(high_cv_ids)
  
  # ── 3. per-position base coverage ──────────────────────────────────────────
  bin_size     <- 100
  pos_cov_rows <- lapply(names(cvg_ok), function(id) {
    d      <- cvg_ok[[id]]
    d$bin  <- floor((d$pos - 1) / bin_size) * bin_size + 1
    agg    <- aggregate(depth ~ bin, data = d, FUN = median)
    agg$ID <- id
    agg
  })
  pos_cov_df          <- do.call(rbind, pos_cov_rows)
  pos_cov_mean        <- aggregate(depth ~ bin, data = pos_cov_df, FUN = mean)
  pos_cov_df$high_cv  <- pos_cov_df$ID %in% high_cv_ids
  
  p_bc <- ggplot() +
    geom_line(data = pos_cov_df[!pos_cov_df$high_cv, ],
              aes(x = bin, y = depth, group = ID),
              colour = "grey70", linewidth = 0.2, alpha = 0.4) +
    geom_line(data = pos_cov_df[pos_cov_df$high_cv, ],
              aes(x = bin, y = depth, group = ID),
              colour = "#E07B3A", linewidth = 0.8, alpha = 0.7) +
    geom_line(data = pos_cov_mean,
              aes(x = bin, y = depth),
              colour = "firebrick", linewidth = 1) +
    scale_x_continuous(breaks = seq(0, 16569, by = 2000),
                       labels = function(x) paste0(x / 1000, "k")) +
    labs(title = "Base coverage across mtDNA",
         x = "mtDNA position", y = "Depth",
         caption = sprintf("Grey = normal CV  |  Orange = CV > %.1f (n=%d)  |  Red = cohort mean",
                           cv_threshold, n_high_cv)) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())
  
  cov_plots$base_coverage <- p_bc
  
  # ── 4. cumulative coverage curve ───────────────────────────────────────────
  depth_thresholds <- c(1, 10, 50, 100, 200, 500, 1000, 2000, 5000)
  
  cumcov_rows <- lapply(names(cvg_ok), function(id) {
    d   <- cvg_ok[[id]]$depth
    pct <- sapply(depth_thresholds, function(t) 100 * mean(d >= t))
    data.frame(ID        = id,
               threshold = depth_thresholds,
               pct_bases = pct,
               stringsAsFactors = FALSE)
  })
  cumcov_df  <- do.call(rbind, cumcov_rows)
  cumcov_med <- aggregate(pct_bases ~ threshold, data = cumcov_df, FUN = median)
  
  p_cum <- ggplot() +
    geom_line(data = cumcov_df,
              aes(x = threshold, y = pct_bases, group = ID),
              colour = "grey70", linewidth = 0.3, alpha = 0.4) +
    geom_line(data = cumcov_med,
              aes(x = threshold, y = pct_bases),
              colour = "firebrick", linewidth = 1.2) +
    scale_x_log10(breaks = depth_thresholds, labels = depth_thresholds) +
    labs(title = "Cumulative coverage (% mtDNA bases ≥ X depth)",
         x = "Depth threshold", y = "% bases covered") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())
  
  cov_plots$cumulative_coverage <- p_cum
  
  # ── 5. CV distribution ─────────────────────────────────────────────────────
  cv_ok <- cv_df[!is.na(cv_df$CV), ]
  
  set.seed(42)
  cv_ok$jitter_y <- runif(nrow(cv_ok), -0.4, 0.4)
  
  p_cv <- ggplot() +
    geom_histogram(data = cv_ok, aes(x = CV),
                   bins = 40, fill = "steelblue", colour = "white") +
    geom_point(data = cv_ok,
               aes(x = CV, y = jitter_y, text = paste0(ID, "<br>CV: ", round(CV, 3))),
               colour = "steelblue", alpha = 0, size = 3) +
    geom_vline(xintercept = median(cv_ok$CV), linetype = "dashed",
               colour = "firebrick", linewidth = 0.7) +
    geom_vline(xintercept = mean(cv_ok$CV),   linetype = "dotted",
               colour = "firebrick", linewidth = 0.7) +
    geom_vline(xintercept = cv_threshold,      linetype = "dashed",
               colour = "orange", linewidth = 0.7) +
    labs(title = sprintf("Coverage CV  |  median = %.3f  |  mean = %.3f  |  %d samples CV > %.1f",
                         median(cv_ok$CV), mean(cv_ok$CV), n_high_cv, cv_threshold),
         x = "CV (sd / mean depth)", y = "# samples") +
    theme_minimal(base_size = 11)
  
  cov_plots$cv_distribution <- ggplotly(p_cv, tooltip = "text")
  
}

message(sprintf("[02_coverage] %d samples  |  %d flagged (depth < %d)  |  %d high CV (> %.1f).",
                nrow(meta), n_flagged_depth, min_median_depth,
                if (length(cvg_ok) > 0) n_high_cv else 0, cv_threshold))