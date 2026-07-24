############ Figure 1 ##############
rm(list = ls())

library(ggstatsplot)
library(rstatix)
library(ggpubr)
library(patchwork)
library(circlize)
library(tidyverse)
library(scales)
library(png)
library(cowplot)
library(patchwork)
library(grid)

#### Figure 1a #####
setwd("/vast/projects/bahlo_mtDNA/Analysis/homoplasmy/allduos")
count <- read.table("trio604_count_hom.tsv", header = T, sep = "\t")
count$FID  <- factor(count$FID)
count$Tool <- factor(count$Tool, levels = c("mitoHPC", "mtSwirl", "mity", "mtDNA_Server2"))
count$HomN <- as.numeric(count$HomN)

stat_test <- count %>%
  pairwise_wilcox_test(
    HomN ~ Tool,
    paired = TRUE,
    p.adjust.method = "BH"
  ) %>%
  add_significance("p.adj") %>%
  mutate(
    y.position = max(count$HomN, na.rm = TRUE) * 1.05 + 
      dplyr::row_number() * 0.05 * max(count$HomN, na.rm = TRUE)
  )

tool_cols <- c(
  "mity" = "#F37520",#F8766D #F37520
  "mitoHPC" = "#50CF4E", #7CAE00 #49A942 #E01545
  "mtDNA_Server2" = "#00BCE7", #00BFC4 #2372B9
  "mtSwirl" = "#7D489C" #C77CFF #7D489C
)

p1 <- ggplot(count, aes(x = Tool, y = HomN, fill = Tool, color = Tool)) +
  
  geom_jitter(
    width = 0.15,
    size = 1.5,
    alpha = 0.6
  ) +
  
  geom_boxplot(
    width = 0.6,
    colour = 'black',
    outlier.shape = NA,
    staplewidth = 0.3,
    alpha = 0.6
  ) +
  
  stat_pvalue_manual(
    stat_test,
    label = "p.adj.signif",
    tip.length = 0.01,
    hide.ns = TRUE
  ) +
  
  labs( x = NULL,
        y = "Overlapping Homoplasmic variants (mother-offspring)"
  ) +
  
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  ) + scale_fill_manual(values = tool_cols) + scale_color_manual(values = tool_cols) + scale_x_discrete(labels = c(
    "mtDNA_Server2" = "mtDNA Server2",
    "mity" = "mity",
    "mitoHPC" = "mitoHPC",
    "mtSwirl" = "mtSwirl"
  ))
#### Figure 1b ####
setwd("/vast/projects/bahlo_mtDNA/")
df = read.csv("/vast/projects/bahlo_mtDNA/Analysis/homoplasmy/allduos/hom_diff_trio604_newdef_reduced.csv",
              sep = ",") # use a table with MNV issues removed
df$max_calls <- apply(df[, c("mity", "mtDNA_Server2", "mitoHPC", "mtSwirl")], 1, max)
df$mity_missed          <- df$max_calls - df$mity
df$mtDNA_Server2_missed <- df$max_calls - df$mtDNA_Server2
df$mitoHPC_missed       <- df$max_calls - df$mitoHPC
df$mtSwirl_missed       <- df$max_calls - df$mtSwirl
df$max_calls <- NULL
df$pos <- as.numeric(str_extract(df$VAR, "^[0-9]+"))

png("./Analysis/homoplasmy/allduos/mitochondria_circos_dot_final2.png",       
    width = 2000, height = 2000,     
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
    
    # Control region (wrap-around)
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
    # Other genes
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
      
      # Conditional facing
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
max_miss <- max(miss_mat, na.rm = TRUE) # Compute overall max across all tools

# Track for points, scaled globally
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
        cex = rescale(vals[nonzero_idx], to = c(0.75,3), from = c(0, max_miss))
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
size_values <- c(1, unique_sizes[15],unique_sizes[20],unique_sizes[22],max_miss)
size_cex <- rescale(size_values, to = c(0.75, 3), from = c(0, max_miss))
legend("bottomleft",  # change position as needed
       legend = size_values,
       pt.cex = size_cex,
       pch = 16,
       cex = 0.75,
       bg = "white",
       y.intersp = c(1,1,1,1,1.5),
       x.intersp = c(1.5,1.5,1.5,1.5,1.7),
       inset = c(0, 0.2))

dev.off()


#### Figure 1c ####
setwd("/vast/projects/bahlo_mtDNA/Analysis")
df_orig = read.csv("/vast/projects/bahlo_mtDNA/Analysis/heteroplasmy/allduos/variants_3tools_newfilters.tsv",sep = "\t")
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

png("./heteroplasmy/allduos/mitochondria_circos_dot_final2.png",       # file name
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
    
    # Control region (wrap-around)
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
size_values <- c(1, unique_sizes[3],unique_sizes[6],unique_sizes[8],max_miss)
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

#### Figure 1d ####
base_path <- "/vast/projects/bahlo_mtDNA/Simulation/1000G"
haplogroups <- read.table("/vast/projects/bahlo_mtDNA/Results_1000G/haplogroups_all.tsv",header = T)


stats_files <- list.files(
  path = base_path,
  pattern = "stats.tsv",
  recursive = TRUE,
  full.names = TRUE
)
remove_paths <- c(
  "/vast/projects/bahlo_mtDNA/Simulation/1000G/HG00428/stats.tsv",
  "/vast/projects/bahlo_mtDNA/Simulation/1000G/HG00739/stats.tsv"
)

stats_files <- setdiff(stats_files, remove_paths)

# Extract <id> from path

stats_df <- map_dfr(stats_files, function(f) {
  id <- basename(dirname(f))  # assumes .../<id>/stats.tsv
  
  read_tsv(f, show_col_types = FALSE) %>%
    mutate(id = id)
})

stats_df <- stats_df %>%
  left_join(
    haplogroups,
    by = c("id" = "SAMPLE")
  )

heatmap_df <- stats_df %>%
  select(
    id,
    Haplogroup,
    tool,
    VAF,
    F1,
    F1_thred
  ) %>%
  mutate(
    VAF_numeric = as.numeric(VAF) / 100,
    
    # Use F1 if VAF == 0.01, otherwise F1_thred
    F1_used = if_else(VAF_numeric == 0.01, F1, F1_thred),
    
    VAF = factor(VAF_numeric),
    tool = factor(tool, levels = c("mitoHPC", "mtSwirl", "mity", "server2")),
    id = factor(id),
    Haplogroup = factor(Haplogroup)
  )

tool_cols <- c(
  "mity" = "#F37520",
  "mitoHPC" = "#50CF4E",
  "server2" = "#00BCE7",
  "mtSwirl" = "#7D489C"
)


boxpl <- ggplot(heatmap_df, aes(x = tool, y = F1_used, fill = tool,color = tool)) +
  geom_jitter(
    #color = "black",
    width = 0.15,
    height = 0,
    size = 1.2,
    alpha = 1,
    stroke=0.85
  ) + 
  geom_boxplot(
    color = "black",
    linewidth = 0.6,
    alpha = 0.6,
    staplewidth = 0.3,
    outliers = FALSE
  ) +
  facet_wrap(~ VAF, nrow = 1) +
  labs(y = "F1",
       fill = "Tool",
       x = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) + scale_color_manual(values = tool_cols) + 
  scale_fill_manual(values = tool_cols) + scale_x_discrete(labels = c(
    "server2" = "mtDNA_Server2",
    "mity" = "mity",
    "mitoHPC" = "mitoHPC",
    "mtSwirl" = "mtSwirl"
  ))

#### Combine into one plot ####

hom <- readPNG("/vast/projects/bahlo_mtDNA/Analysis/homoplasmy/allduos/mitochondria_circos_dot_final2.png")
het <- readPNG("/vast/projects/bahlo_mtDNA/Analysis/heteroplasmy/allduos/mitochondria_circos_dot_final2.png")

hom_grob <- rasterGrob(hom, interpolate = TRUE)
het_grob <- rasterGrob(het, interpolate = TRUE)

p_hom <- wrap_elements(full=hom_grob)
p_het <- wrap_elements(full=het_grob)


final_plot <- 
  (p1 | p_hom) /
  (p_het | boxpl) +
  plot_annotation(tag_levels = "A") & theme(plot.tag=element_text(size=14,face="bold"))

ggsave(
  filename = "/vast/projects/bahlo_mtDNA/Analysis/final_figures/Figure_1_newcol3.tiff",
  plot = final_plot,
  device = "tiff",
  width = 10,          # inches (adjust as needed)
  height = 10,         # inches
  units = "in",
  dpi = 600,          # high resolution
  compression = "lzw"
)
