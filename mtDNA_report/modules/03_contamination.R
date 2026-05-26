# Module 03 — contamination QC (haplocheck)
#
# Expects in caller environment (from 01_load_data + config):
#   samples           — list with $meta
#   results_root, cohort_ids, contam_threshold
#
# Produces:
#   contam_plots$level_per_sample  — contamination level per sample dot plot
#   contam_plots$level_vs_depth    — contam level vs median depth scatter

library(ggplot2)
library(plotly)

contam_plots <- list()

# ── reload haplocheck with extra columns (Distance, SampleCoverage) ──────────
hc_list <- lapply(cohort_ids, function(id) {
  f <- tab_path(id, "mutect2.haplocheck.tab")
  if (!file.exists(f)) return(NULL)
  d <- tryCatch(read.table(f, header = TRUE, sep = "\t",
                           stringsAsFactors = FALSE),
                error = function(e) NULL)
  if (is.null(d)) return(NULL)
  d$ID <- id
  d$ContaminationLevel <- suppressWarnings(
    ifelse(d$ContaminationLevel == "ND", 0, as.numeric(d$ContaminationLevel))
  )
  d
})
hc_df <- do.call(rbind, Filter(Negate(is.null), hc_list))

if (is.null(hc_df) || nrow(hc_df) == 0) {
  message("[03_contamination] No haplocheck data found.")
} else {
  
  hc_df$flagged <- hc_df$ContaminationLevel >= contam_threshold
  
  # ── 1. contamination level per sample ──────────────────────────────────────
  hc_sorted        <- hc_df[hc_df$ContaminationLevel > 0, ]
  hc_sorted        <- hc_sorted[order(hc_sorted$ContaminationLevel), ]
  hc_sorted$ID_ord <- factor(hc_sorted$ID, levels = hc_sorted$ID)
  hc_sorted$colour <- ifelse(hc_sorted$flagged, "Flagged", "Pass")
  
  n_flagged   <- sum(hc_df$flagged, na.rm = TRUE)
  plot_height <- max(400, nrow(hc_sorted) * 14)
  
  p <- ggplot(hc_sorted,
              aes(x = ContaminationLevel, y = ID_ord, colour = colour,
                  text = paste0(ID, "<br>Contam: ",
                                round(ContaminationLevel * 100, 2), "%"))) +
    geom_point(size = 2) +
    geom_vline(xintercept = contam_threshold, linetype = "dashed",
               colour = "firebrick", linewidth = 0.5) +
    scale_colour_manual(values = c("Pass" = "steelblue", "Flagged" = "firebrick"),
                        name = NULL) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(title = sprintf("Contamination level per sample  |  threshold = %.0f%%  |  %d / %d flagged",
                         contam_threshold * 100, n_flagged, nrow(hc_df)),
         x = "Contamination level", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(axis.text.y        = element_text(size = 6),
          panel.grid.major.y = element_line(linewidth = 0.2))
  
  contam_plots$level_per_sample <- ggplotly(p, tooltip = "text",
                                            height = plot_height) |>
    layout(margin = list(l = 120))
  
  # ── 2. contamination level vs median depth ─────────────────────────────────
  meta <- samples$meta
  scatter_df <- merge(hc_df[, c("ID", "ContaminationLevel", "flagged")],
                      meta[, c("ID", "depth_median")],
                      by = "ID")
  
  scatter_df$colour <- ifelse(scatter_df$flagged, "Flagged", "Pass")
  
  p2 <- ggplot(scatter_df,
               aes(x = depth_median, y = ContaminationLevel, colour = colour,
                   text = paste0(ID, "<br>Depth: ", round(depth_median), "x",
                                 "<br>Contam: ", round(ContaminationLevel * 100, 2), "%"))) +
    geom_point(size = 2, alpha = 0.7) +
    geom_hline(yintercept = contam_threshold, linetype = "dashed",
               colour = "firebrick", linewidth = 0.5) +
    scale_colour_manual(values = c("Pass" = "steelblue", "Flagged" = "firebrick"),
                        name = NULL) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(title = "Contamination level vs median depth",
         x = "Median depth", y = "Contamination level") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "top")
  
  contam_plots$level_vs_depth <- ggplotly(p2, tooltip = "text")
  
}

message(sprintf("[03_contamination] %d samples  |  %d flagged (>= %.0f%% contamination).",
                nrow(hc_df), sum(hc_df$flagged, na.rm = TRUE),
                contam_threshold * 100))

