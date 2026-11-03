#### Supplementary Figure 5 ####
# Heteroplasmic variant VAF distribution and comparsion across tools
rm(list=ls())
library(ggstatsplot)
library(tidyr)
library(dplyr)
library(ggplot2)
library(rstatix)
library(ggpubr)
library(ggridges)
library(patchwork)
library(grid)


setwd("/vast/projects/bahlo_mtDNA/Analysis/heteroplasmy/allduos")

##### Figure 5a #####
data <- read.table("trio604_vaf03_het_newfilters.tsv", header = TRUE, sep = "\t")

overlaps <- data %>%
  group_by(SAMPLE, VAR) %>%
  mutate(count = n()) %>%
  filter(count == 4) %>%
  ungroup() %>%
  select(-count)

overlaps$pair <- paste0(overlaps$SAMPLE,"_",overlaps$VAR)
overlaps <- overlaps %>%
  mutate(TOOL = factor(TOOL, levels = c("mitoHPC", "mtSwirl", "mity", "mtDNA Server2")))

overlaps <- overlaps %>%
  mutate(pair = as.factor(pair))

overlaps_wide <- overlaps %>%
  select(SAMPLE, VAR, TOOL, VAF) %>%
  pivot_wider(
    names_from = TOOL,
    values_from = VAF
  )

overlaps_long <- overlaps_wide %>%
  pivot_longer(
    cols = c("mity", "mtDNA Server2", "mitoHPC", "mtSwirl"),
    names_to = "TOOL",
    values_to = "VAF"
  )

overlaps_long <- overlaps_long %>%
  mutate(pair_id = paste0(SAMPLE, "_", VAR))

tool_cols <- c(
  "mity" = "#F37520",
  "mitoHPC" = "#50CF4E", 
  "mtDNA Server2" = "#00BCE7", 
  "mtSwirl" = "#7D489C" 
)

p_vals <- overlaps_long %>%
  pairwise_wilcox_test(
    VAF ~ TOOL,
    paired = TRUE,
    id = "pair_id",
    p.adjust.method = "BH"
  )

eff_sizes <- overlaps_long %>%
  wilcox_effsize(
    VAF ~ TOOL, 
    paired = TRUE, 
    id = "pair_id"
  )

# 2. Join and format the label on one line
pairwise_res <- p_vals %>%
  left_join(eff_sizes, by = c("group1", "group2")) %>%
  add_significance("p.adj") %>%
  mutate(
    # Format: p=0.01, r=0.45
    label = paste0("p=", signif(p.adj, 2), ", r=", signif(effsize, 2)),
    y.position = max(overlaps_long$VAF, na.rm = TRUE) * 1.05 + 
      (dplyr::row_number() * 0.06 * max(overlaps_long$VAF, na.rm = TRUE))
  )

p_overlap2 <- ggplot(overlaps, aes(x = TOOL, y = VAF, fill = TOOL, color = TOOL)) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
  geom_boxplot(width = 0.6, colour = 'black', outlier.shape = NA, alpha = 0.6,staplewidth = 0.3) +
  stat_pvalue_manual(
    pairwise_res,
    label = "label",
    tip.length = 0.01,
    hide.ns = TRUE,
    size = 3.2    
  ) +
  labs(x = NULL, y = "VAF",title = "Comparison of Heteroplasmic VAFs across Tools") +
  theme_minimal(base_size = 12) + 
  theme(legend.position = "none") + 
  scale_fill_manual(values = tool_cols) + 
  scale_color_manual(values = tool_cols)

p_overlap2

##### Figure 5b #####
densityp <- ggplot(overlaps, aes(x = VAF, y = TOOL, fill = TOOL)) +
  geom_density_ridges(
    alpha = 0.7, 
    scale = 2.5,         # Values > 1.0 make the lanes "bleed" into each other/overlap
    rel_min_height = 0.01, # Cuts off the thin tails to keep it clean
    color = "black"
  ) +
  scale_fill_manual(values = tool_cols) +
  labs(
    x = "VAF",
    y = NULL,            # Removing "Tool" label as the y-axis text is clear
    title = "VAF Distributions Across Tools"
  ) +
  # theme_ridges() is optimized for this look
  theme_minimal(base_size = 12) + 
  theme(
    legend.position = "none",
    axis.title.x = element_text(hjust = 0.5)
  )

##### Combine into 1 plot #####
final_plot <- 
  (p_overlap2 | densityp) + 
  plot_annotation(tag_levels = "A") & theme(plot.tag=element_text(size=14,face="bold"))

ggsave(
  filename = "/vast/projects/bahlo_mtDNA/Analysis/final_figures/SupFig_4_final.tiff",
  plot = final_plot,
  device = "tiff",
  width = 12,          
  height = 8,         
  units = "in",
  dpi = 600,          
  compression = "lzw"
)

