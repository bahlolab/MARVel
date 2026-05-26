# Module 04 — haplogroup distribution
#
# Expects in caller environment (from 01_load_data):
#   samples  — list with $meta (haplogroup, macrohap, clade, hg_quality)
#   m_clade  — character vector of M-clade macro-haplogroups
#
# Produces:
#   hg_plots$top_macro_haplogroups — all macro-haplogroups by frequency
#   hg_plots$quality               — haplogroup assignment quality histogram

library(ggplot2)
library(plotly)

hg_plots <- list()

meta <- samples$meta

clade_cols <- c("L"       = "#B39DDB",
                "M"       = "#A5C882",
                "N/R"     = "#80CBC4",
                "Unknown" = "grey70")

meta$clade[is.na(meta$clade)] <- "Unknown"

# ── 1. all macro-haplogroups ──────────────────────────────────────────────────
hg_df <- as.data.frame(table(Macrohap = meta$macrohap))
hg_df <- hg_df[!is.na(hg_df$Macrohap) & hg_df$Macrohap != "NA", ]
hg_df <- hg_df[order(-hg_df$Freq), ]
hg_df$pct   <- 100 * hg_df$Freq / nrow(meta)
hg_df$clade <- ifelse(grepl("^L", hg_df$Macrohap), "L",
                      ifelse(hg_df$Macrohap %in% m_clade,  "M", "N/R"))
hg_df$Macrohap <- factor(hg_df$Macrohap, levels = rev(hg_df$Macrohap))

hg_plots$top_macro_haplogroups <- ggplot(hg_df,
                                         aes(x = Freq, y = Macrohap, fill = clade)) +
  geom_col() +
  geom_text(aes(label = sprintf("%d (%.1f%%)", Freq, pct)),
            hjust = -0.1, size = 3) +
  scale_fill_manual(values = clade_cols, name = "Clade") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(title = "Macro-haplogroup distribution",
       x = "# samples", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank(),
        legend.position = "top")

# ── 2. haplogroup quality ─────────────────────────────────────────────────────
quality_df <- meta[!is.na(meta$hg_quality), ]

set.seed(42)
quality_df$jitter_y <- runif(nrow(quality_df), -0.4, 0.4)

p_q <- ggplot() +
  geom_histogram(data = quality_df, aes(x = hg_quality),
                 bins = 40, fill = "steelblue", colour = "white") +
  geom_point(data = quality_df,
             aes(x = hg_quality, y = jitter_y,
                 text = paste0(ID, "<br>Quality: ", round(hg_quality, 3))),
             colour = "steelblue", alpha = 0, size = 3) +
  geom_vline(xintercept = 0.5, linetype = "dashed",
             colour = "firebrick", linewidth = 0.7) +
  geom_vline(xintercept = 0.8, linetype = "dashed",
             colour = "orange", linewidth = 0.7) +
  annotate("text", x = 0.5, y = Inf, label = "50%",
           hjust = 1.2, vjust = 1.5, size = 3, colour = "firebrick") +
  annotate("text", x = 0.8, y = Inf, label = "80%",
           hjust = 1.2, vjust = 1.5, size = 3, colour = "orange") +
  labs(title = sprintf("Haplogroup assignment quality  |  median = %.3f  |  range = %.3f–%.3f",
                       median(quality_df$hg_quality),
                       min(quality_df$hg_quality),
                       max(quality_df$hg_quality)),
       x = "Quality score", y = "# samples") +
  theme_minimal(base_size = 11)

hg_plots$quality <- ggplotly(p_q, tooltip = "text")

message(sprintf("[04_haplogroup] %d samples  |  clades: %s",
                nrow(meta),
                paste(names(table(meta$clade)), table(meta$clade),
                      sep = "=", collapse = "  ")))