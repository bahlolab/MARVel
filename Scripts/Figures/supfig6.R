#### Supplementary Figure 6 ####
# Heteroplasmic variant VAF comparsion across tools and visualisation in Haplogroups L and H
rm(list=ls())
library(tidyr)
library(dplyr)
library(ggplot2)
library(rstatix)
library(ggpubr)
library(ggridges)
library(patchwork)
library(grid)
library(circlize)
library(tidyverse)
library(scales)
library(png)

setwd("/vast/projects/bahlo_mtDNA/Analysis/heteroplasmy/allduos")

##### Figure 6a #####
df_orig = read.csv("/vast/projects/bahlo_mtDNA/Analysis/heteroplasmy/allduos/variants_3tools_newfilters_hapH.tsv",sep = "\t")

df<- df_orig %>%
  group_by(pos) %>%
  summarise(
    across(where(is.integer), ~ sum(.x, na.rm = TRUE)),
    across(where(~ !is.integer(.x)), 
           ~ paste(.x, collapse = ",")),
    .groups = "drop"
  )

df$max_calls <- apply(df[, c("mity", "mtDNA.Server2", "mitoHPC", "mtSwirl")], 1, max)

df$mity_missed          <- df$max_calls - df$mity
df$mtDNA_Server2_missed <- df$max_calls - df$mtDNA.Server2
df$mitoHPC_missed       <- df$max_calls - df$mitoHPC
df$mtSwirl_missed       <- df$max_calls - df$mtSwirl

df$max_calls <- NULL

# plot
png("/vast/projects/bahlo_mtDNA/Analysis/heteroplasmy/allduos/mitochondria_circos_dot_hapH.png",
    width = 2000, height = 2000,     # pixels
    res = 300)  

genome_len <- 16569
circos.clear()
circos.par(start.degree = 90, gap.degree = 2, track.margin = c(0.02, 0.02))

circos.initialize(factors = "chrM", xlim = c(1, genome_len))

circos.trackPlotRegion(
  factors = "chrM",
  ylim = c(0, 1),
  track.height = 0.05,
  bg.border = NA,
  panel.fun = function(x, y) {
    
    # MT-RNR2 region
    highlight_start <- 1671
    highlight_end   <- 3229
    circos.rect(
      xleft = highlight_start,
      ybottom = 0,
      xright = highlight_end,
      ytop = 1,
      col = "lightblue",
      border = NA
    )
    circos.text(
      x = (highlight_start + highlight_end)/2,
      y = 2,
      labels = "MT-RNR2",
      cex = 0.6,
      facing = "bending.inside",
      adj = c(0.5, 0.5)
    )
    
    # Control region
    start2 <- 16024
    end2 <- 576
    circos.rect(start2, 0, genome_len, 1, col = "orange", border = NA)
    circos.rect(1, 0, end2, 1, col = "orange", border = NA)
    
    L1 <- genome_len - start2 + 1
    L2 <- end2 - 1 + 1
    L  <- L1 + L2
    midpoint_from_start <- L / 2
    if(midpoint_from_start <= L1){
      midpoint_pos <- start2 + midpoint_from_start - 1
    } else {
      midpoint_pos <- (midpoint_from_start - L1)
    }
    
    circos.text(
      x = midpoint_pos,
      y = 2,
      labels = "Control-region",
      cex = 0.6,
      facing = "bending.inside",
      adj = c(0.5, 0.5)
    )
    
    # Mitochondrial genes
    mito_genes <- data.frame(
      gene  = c("MT-CO1","MT-TS1","MT-TD","MT-CO2","MT-TK","MT-ATP8","MT-ATP6"),
      start = c(5904, 7446, 7518, 7586, 8295, 8366, 8527),
      end   = c(7445, 7514, 7585, 8269, 8364, 8572, 9207),
      col   = c("indianred2","lightgreen","darkblue","mediumpurple","pink","yellow","cyan")
    )
    clockwise_genes <- c("MT-TS1","MT-TD","MT-TK","MT-ATP8")
    
    for(i in 1:nrow(mito_genes)){
      circos.rect(
        xleft   = mito_genes$start[i],
        ybottom = 0,
        xright  = mito_genes$end[i],
        ytop    = 1,
        col     = mito_genes$col[i],
        border  = NA
      )
      
      label_facing <- ifelse(mito_genes$gene[i] %in% clockwise_genes,
                             "clockwise",
                             "bending.outside")
      
      label_y <- ifelse(mito_genes$gene[i] %in% clockwise_genes, 3, 2)
      label_cex <- ifelse(mito_genes$gene[i] %in% clockwise_genes, 0.4, 0.6)
      
      circos.text(
        x      = (mito_genes$start[i] + mito_genes$end[i])/2,
        y      = label_y,
        labels = mito_genes$gene[i],
        cex    = label_cex,
        facing = label_facing,
        adj    = c(0.5, 0.5)
      )
    }
    
  } 
) 

# Axis track
circos.trackPlotRegion(
  factors = "chrM",
  ylim = c(0, 1),
  track.height = 0.06,
  bg.border = NA,
  panel.fun = function(x, y) {
    ticks <- seq(0, genome_len, by = 500)
    circos.axis(major.at = ticks,
                labels.cex = 0.5,
                labels.niceFacing = TRUE,
                direction = "inside",
                major.tick.length = mm_y(2))
  }
)

miss_mat <- as.matrix(df[, c("mity_missed",
                             "mtDNA_Server2_missed",
                             "mitoHPC_missed",
                             "mtSwirl_missed")])

tool_cols <- c(
  mity_missed         = "#F37520",
  mtDNA_Server2_missed = "#00BCE7",
  mitoHPC_missed      = "#50CF4E",
  mtSwirl_missed      = "#7D489C"
)

ring_y <- list(
  mity_missed          = 0.15,
  mtDNA_Server2_missed = 0.35,
  mitoHPC_missed       = 0.55,
  mtSwirl_missed       = 0.75
)

max_miss <- max(miss_mat, na.rm = TRUE)

circos.trackPlotRegion(
  factors = "chrM",
  ylim = c(0, 1),
  track.height = 0.35,
  bg.border = ,
  bg.col= "beige",
  track.margin = c(0.02,0.05),
  panel.fun = function(x, y) {
    for (tool in colnames(miss_mat)) {
      vals <- miss_mat[, tool]
      nonzero_idx <- which(vals > 0)
      
      # Scale sizes based on global max
      circos.points(
        x = df$pos[nonzero_idx],
        y = rep(ring_y[[tool]], length(nonzero_idx)),
        col = scales::alpha(tool_cols[tool], 0.7),
        pch = 16,
        cex = rescale(vals[nonzero_idx], to = c(0.5, 2), from = c(0, max_miss))
      )
    }
  }
)


# Legend
legend("bottomleft",
       legend = c("mity", "mtDNA Server2", "mitoHPC", "mtSwirl"),
       col = tool_cols,
       pch = 16,
       pt.cex = 1,
       cex = 0.75,
       bg = "white")

miss_cols <- c("mity_missed", "mtDNA_Server2_missed", "mitoHPC_missed", "mtSwirl_missed")
all_sizes <- unlist(df[miss_cols])
unique_sizes <- sort(unique(all_sizes))
size_values <- c(1,max_miss)#, unique_sizes[2],unique_sizes[6],unique_sizes[8]
size_cex <- rescale(size_values, to = c(0.5, 2), from = c(0, max_miss))
legend("bottomleft",  # change position as needed
       legend = size_values,
       pt.cex = size_cex,
       pch = 16,
       cex = 0.75,
       bg = "white",
       y.intersp = c(1,1,1,1,1.5),
       x.intersp = c(1.5,1.5,1.5,1.5,1.4),
       inset = c(0, 0.2))
# title("mtDNA Heteroplasmic Variants Detected By 3 out of 4 Tools", line = -0.25)

dev.off()

##### Figure 6b #####
data <- read.table("trio604_summary_het_newfilters_hapH.tsv", header = TRUE, sep = "\t")

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

pairwise_res <- p_vals %>%
  left_join(eff_sizes, by = c("group1", "group2")) %>%
  add_significance("p.adj") %>%
  mutate(
    # Format: p=0.01, r=0.45
    label = paste0("p=", signif(p.adj, 2), ", r=", signif(effsize, 2)),
    y.position = max(overlaps_long$VAF, na.rm = TRUE) * 1.05 + 
      (dplyr::row_number() * 0.06 * max(overlaps_long$VAF, na.rm = TRUE))
  )

p_overlap <- ggplot(overlaps, aes(x = TOOL, y = VAF, fill = TOOL, color = TOOL)) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
  geom_boxplot(width = 0.6, colour = 'black', outlier.shape = NA, alpha = 0.6,staplewidth = 0.3) +
  stat_pvalue_manual(
    pairwise_res,
    label = "label",
    tip.length = 0.01,
    hide.ns = TRUE,
    size = 3.2    
  ) +
  labs(x = NULL, y = "VAF") +
  theme_minimal(base_size = 12) + 
  theme(legend.position = "none") + 
  scale_fill_manual(values = tool_cols) + 
  scale_color_manual(values = tool_cols)

##### Figure 6c #####
df_origL = read.csv("/vast/projects/bahlo_mtDNA/Analysis/heteroplasmy/allduos/variants_3tools_newfilters_hapL.tsv",sep = "\t")

dfL<- df_origL %>%
  group_by(pos) %>%
  summarise(
    across(where(is.integer), ~ sum(.x, na.rm = TRUE)),
    across(where(~ !is.integer(.x)), 
           ~ paste(.x, collapse = ",")),
    .groups = "drop"
  )

dfL$max_calls <- apply(dfL[, c("mity", "mtDNA.Server2", "mitoHPC", "mtSwirl")], 1, max)

dfL$mity_missed          <- dfL$max_calls - dfL$mity
dfL$mtDNA_Server2_missed <- dfL$max_calls - dfL$mtDNA.Server2
dfL$mitoHPC_missed       <- dfL$max_calls - dfL$mitoHPC
dfL$mtSwirl_missed       <- dfL$max_calls - dfL$mtSwirl

dfL$max_calls <- NULL

png("/vast/projects/bahlo_mtDNA/Analysis/heteroplasmy/allduos/mitochondria_circos_dot_hapL.png",
    width = 2000, height = 2000,     # pixels
    res = 300)  

genome_len <- 16569
circos.clear()
circos.par(start.degree = 90, gap.degree = 2, track.margin = c(0.02, 0.02))

circos.initialize(factors = "chrM", xlim = c(1, genome_len))

circos.trackPlotRegion(
  factors = "chrM",
  ylim = c(0, 1),
  track.height = 0.05,
  bg.border = NA,
  panel.fun = function(x, y) {
    
    # MT-RNR2 region
    highlight_start <- 1671
    highlight_end   <- 3229
    circos.rect(
      xleft = highlight_start,
      ybottom = 0,
      xright = highlight_end,
      ytop = 1,
      col = "lightblue",
      border = NA
    )
    circos.text(
      x = (highlight_start + highlight_end)/2,
      y = 2,
      labels = "MT-RNR2",
      cex = 0.6,
      facing = "bending.inside",
      adj = c(0.5, 0.5)
    )
    
    # Control region
    start2 <- 16024
    end2 <- 576
    circos.rect(start2, 0, genome_len, 1, col = "orange", border = NA)
    circos.rect(1, 0, end2, 1, col = "orange", border = NA)
    
    L1 <- genome_len - start2 + 1
    L2 <- end2 - 1 + 1
    L  <- L1 + L2
    midpoint_from_start <- L / 2
    if(midpoint_from_start <= L1){
      midpoint_pos <- start2 + midpoint_from_start - 1
    } else {
      midpoint_pos <- (midpoint_from_start - L1)
    }
    
    circos.text(
      x = midpoint_pos,
      y = 2,
      labels = "Control-region",
      cex = 0.6,
      facing = "bending.inside",
      adj = c(0.5, 0.5)
    )
    
    # Mitochondrial genes
    mito_genes <- data.frame(
      gene  = c("MT-CO1","MT-TS1","MT-TD","MT-CO2","MT-TK","MT-ATP8","MT-ATP6"),
      start = c(5904, 7446, 7518, 7586, 8295, 8366, 8527),
      end   = c(7445, 7514, 7585, 8269, 8364, 8572, 9207),
      col   = c("indianred2","lightgreen","darkblue","mediumpurple","pink","yellow","cyan")
    )
    clockwise_genes <- c("MT-TS1","MT-TD","MT-TK","MT-ATP8")
    
    for(i in 1:nrow(mito_genes)){
      circos.rect(
        xleft   = mito_genes$start[i],
        ybottom = 0,
        xright  = mito_genes$end[i],
        ytop    = 1,
        col     = mito_genes$col[i],
        border  = NA
      )
      
      label_facing <- ifelse(mito_genes$gene[i] %in% clockwise_genes,
                             "clockwise",
                             "bending.outside")
      
      label_y <- ifelse(mito_genes$gene[i] %in% clockwise_genes, 3, 2)
      label_cex <- ifelse(mito_genes$gene[i] %in% clockwise_genes, 0.4, 0.6)
      
      circos.text(
        x      = (mito_genes$start[i] + mito_genes$end[i])/2,
        y      = label_y,
        labels = mito_genes$gene[i],
        cex    = label_cex,
        facing = label_facing,
        adj    = c(0.5, 0.5)
      )
    }
    
  } 
) 

# Axis track
circos.trackPlotRegion(
  factors = "chrM",
  ylim = c(0, 1),
  track.height = 0.06,
  bg.border = NA,
  panel.fun = function(x, y) {
    ticks <- seq(0, genome_len, by = 500)
    circos.axis(major.at = ticks,
                labels.cex = 0.5,
                labels.niceFacing = TRUE,
                direction = "inside",
                major.tick.length = mm_y(2))
  }
)

miss_mat <- as.matrix(dfL[, c("mity_missed",
                             "mtDNA_Server2_missed",
                             "mitoHPC_missed",
                             "mtSwirl_missed")])

tool_cols <- c(
  mity_missed         = "#F37520",
  mtDNA_Server2_missed = "#00BCE7",
  mitoHPC_missed      = "#50CF4E",
  mtSwirl_missed      = "#7D489C"
)

ring_y <- list(
  mity_missed          = 0.15,
  mtDNA_Server2_missed = 0.35,
  mitoHPC_missed       = 0.55,
  mtSwirl_missed       = 0.75
)

max_miss <- max(miss_mat, na.rm = TRUE)

circos.trackPlotRegion(
  factors = "chrM",
  ylim = c(0, 1),
  track.height = 0.35,
  bg.border = ,
  bg.col= "beige",
  track.margin = c(0.02,0.05),
  panel.fun = function(x, y) {
    for (tool in colnames(miss_mat)) {
      vals <- miss_mat[, tool]
      nonzero_idx <- which(vals > 0)
      
      # Scale sizes based on global max
      circos.points(
        x = dfL$pos[nonzero_idx],
        y = rep(ring_y[[tool]], length(nonzero_idx)),
        col = scales::alpha(tool_cols[tool], 0.7),
        pch = 16,
        cex = rescale(vals[nonzero_idx], to = c(0.5, 2), from = c(0, max_miss))
      )
    }
  }
)


# Legend
legend("bottomleft",
       legend = c("mity", "mtDNA Server2", "mitoHPC", "mtSwirl"),
       col = tool_cols,
       pch = 16,
       pt.cex = 1,
       cex = 0.75,
       bg = "white")

miss_cols <- c("mity_missed", "mtDNA_Server2_missed", "mitoHPC_missed", "mtSwirl_missed")
all_sizes <- unlist(dfL[miss_cols])
unique_sizes <- sort(unique(all_sizes))
size_values <- c(1,unique_sizes[3],max_miss)#, unique_sizes[2],unique_sizes[6],unique_sizes[8]
size_cex <- rescale(size_values, to = c(0.5, 2), from = c(0, max_miss))
legend("bottomleft",  # change position as needed
       legend = size_values,
       pt.cex = size_cex,
       pch = 16,
       cex = 0.75,
       bg = "white",
       y.intersp = c(1,1,1,1,1.5),
       x.intersp = c(1.5,1.5,1.5,1.5,1.4),
       inset = c(0, 0.2))
# title("mtDNA Heteroplasmic Variants Detected By 3 out of 4 Tools", line = -0.25)

dev.off()

##### Figure 6d #####
dataL <- read.table("trio604_summary_het_newfilters_hapL.tsv", header = TRUE, sep = "\t")

overlapsL <- dataL %>%
  group_by(SAMPLE, VAR) %>%
  mutate(count = n()) %>%
  filter(count == 4) %>%
  ungroup() %>%
  select(-count)

overlapsL$pair <- paste0(overlapsL$SAMPLE,"_",overlapsL$VAR)
overlapsL <- overlapsL %>%
  mutate(TOOL = factor(TOOL, levels = c("mitoHPC", "mtSwirl", "mity", "mtDNA Server2")))

tool_cols <- c(
  "mity" = "#F37520",
  "mitoHPC" = "#50CF4E",
  "mtDNA Server2" = "#00BCE7",
  "mtSwirl" = "#7D489C"
)

overlapsL <- overlapsL %>%
  mutate(pair = as.factor(pair))

overlaps_wideL <- overlapsL %>%
  select(SAMPLE, VAR, TOOL, VAF) %>%
  pivot_wider(
    names_from = TOOL,
    values_from = VAF
  )

overlaps_longL <- overlaps_wideL %>%
  pivot_longer(
    cols = c("mity", "mtDNA Server2", "mitoHPC", "mtSwirl"),
    names_to = "TOOL",
    values_to = "VAF"
  )

overlaps_longL <- overlaps_longL %>%
  mutate(pair_id = paste0(SAMPLE, "_", VAR))

p_valsL <- overlaps_longL %>%
  pairwise_wilcox_test(
    VAF ~ TOOL,
    paired = TRUE,
    id = "pair_id",
    p.adjust.method = "BH"
  )

eff_sizesL <- overlaps_longL %>%
  wilcox_effsize(
    VAF ~ TOOL, 
    paired = TRUE, 
    id = "pair_id"
  )

pairwise_resL <- p_valsL %>%
  left_join(eff_sizesL, by = c("group1", "group2")) %>%
  add_significance("p.adj") %>%
  mutate(
    # Format: p=0.01, r=0.45
    label = paste0("p=", signif(p.adj, 2), ", r=", signif(effsize, 2)),
    y.position = max(overlaps_longL$VAF, na.rm = TRUE) * 1.05 + 
      (dplyr::row_number() * 0.06 * max(overlaps_longL$VAF, na.rm = TRUE))
  )

p_overlapL <- ggplot(overlapsL, aes(x = TOOL, y = VAF, fill = TOOL, color = TOOL)) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
  geom_boxplot(width = 0.6, colour = 'black', outlier.shape = NA, alpha = 0.6,staplewidth = 0.3) +
  stat_pvalue_manual(
    pairwise_resL,
    label = "label",
    tip.length = 0.01,
    hide.ns = TRUE,
    size = 3.2    
  ) +
  labs(x = NULL, y = "VAF") +
  theme_minimal(base_size = 12) + 
  theme(legend.position = "none") + 
  scale_fill_manual(values = tool_cols) + 
  scale_color_manual(values = tool_cols)

##### Combine into 1 plot #####
hetL <- readPNG("/vast/projects/bahlo_mtDNA/Analysis/heteroplasmy/allduos/mitochondria_circos_dot_hapL.png")
hetH <- readPNG("/vast/projects/bahlo_mtDNA/Analysis/heteroplasmy/allduos/mitochondria_circos_dot_hapH.png")

hetL_grob <- rasterGrob(hetL, interpolate = TRUE)
hetH_grob <- rasterGrob(hetH, interpolate = TRUE)

p_hetL <- wrap_elements(full=hetL_grob)
p_hetH <- wrap_elements(full=hetH_grob)


final_plot <- 
  (p_hetH | p_overlap) /
  (p_hetL | p_overlapL) + 
  plot_annotation(tag_levels = "A") & theme(plot.tag=element_text(size=14,face="bold"))

ggsave(
  filename = "/vast/projects/bahlo_mtDNA/Analysis/final_figures/SupFig_5_final.tiff",
  plot = final_plot,
  device = "tiff",
  width = 10,          # inches (adjust as needed)
  height = 10,         # inches
  units = "in",
  dpi = 600,          # high resolution
  compression = "lzw"
)
