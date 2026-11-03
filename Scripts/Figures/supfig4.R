#### Supplementary Figure 4 ####
# Homoplasmic variants in Haplogroup L and H
rm(list = ls())
library(circlize)
library(tidyverse)
library(scales)
library(png)
library(cowplot)
library(patchwork)
library(grid)
setwd("/vast/projects/bahlo_mtDNA/")

#### Figure 4a ####
df = read.table("/vast/projects/bahlo_mtDNA/Analysis/homoplasmy/allduos/hom_diff_trio_H_reduced.tsv",
                sep = "\t",header = T) # use a table with MNV issues removed
df$max_calls <- apply(df[, c("mity", "mtDNA_Server2", "mitoHPC", "mtSwirl")], 1, max)

df$mity_missed          <- df$max_calls - df$mity
df$mtDNA_Server2_missed <- df$max_calls - df$mtDNA_Server2
df$mitoHPC_missed       <- df$max_calls - df$mitoHPC
df$mtSwirl_missed       <- df$max_calls - df$mtSwirl

df$max_calls <- NULL
df$pos <- as.numeric(str_extract(df$VAR, "^[0-9]+"))

# Plot
png("./Analysis/homoplasmy/allduos/mitochondria_circos_dot_hapH.png",       
    width = 2500, height = 2000,     
    res = 300)  

genome_len <- 16569
circos.clear()
circos.par(start.degree = 90, gap.degree = 2, track.margin = c(0.02, 0.02))

circos.initialize(factors = "chrM", xlim = c(1, genome_len))

# Highlight regions track (outside)
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
    
    # midpoint for text in circular space
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
    
    # Genes you want to rotate differently
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

# Colors with transparency
tool_cols <- c(
  mity_missed         = "#F37520",
  mtDNA_Server2_missed = "#00BCE7",
  mitoHPC_missed      = "#50CF4E",
  mtSwirl_missed      = "#7D489C"
)

# Radial positions for each tool
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
size_values <- c(1, unique_sizes[3],max_miss)
size_cex <- rescale(size_values, to = c(0.75, 3), from = c(0, max_miss))
legend("bottomleft",  # change position as needed
       legend = size_values,
       pt.cex = size_cex,
       pch = 16,
       cex = 0.75,
       bg = "white",
       y.intersp = c(1,1,1.2,1.2,1.7),
       x.intersp = c(1.5,1.5,1.4,1.5,1.5),
       inset = c(0, 0.2))

#title("Haplogroup H: Homoplasmic Variant Detection Gaps", line = -0.25)

dev.off()

#### Figure 4b ####
df = read.table("/vast/projects/bahlo_mtDNA/Analysis/homoplasmy/allduos/hom_diff_trio_L_reduced.tsv",
                sep = "\t",header = T) # use a table with MNV issues removed
df$max_calls <- apply(df[, c("mity", "mtDNA_Server2", "mitoHPC", "mtSwirl")], 1, max)

df$mity_missed          <- df$max_calls - df$mity
df$mtDNA_Server2_missed <- df$max_calls - df$mtDNA_Server2
df$mitoHPC_missed       <- df$max_calls - df$mitoHPC
df$mtSwirl_missed       <- df$max_calls - df$mtSwirl

df$max_calls <- NULL
df$pos <- as.numeric(str_extract(df$VAR, "^[0-9]+"))

# Plot
png("./Analysis/homoplasmy/allduos/mitochondria_circos_dot_hapL.png",       
    width = 2500, height = 2000,     
    res = 300)  

genome_len <- 16569
circos.clear()
circos.par(start.degree = 90, gap.degree = 2, track.margin = c(0.02, 0.02))

circos.initialize(factors = "chrM", xlim = c(1, genome_len))

# Highlight regions track (outside)
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
    
    # midpoint for text in circular space
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
    
    # Genes you want to rotate differently
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

# Colors
tool_cols <- c(
  mity_missed         = "#F37520", #FF000080
  mtDNA_Server2_missed = "#00BCE7", #0000FF80
  mitoHPC_missed      = "#50CF4E", #00FF0080
  mtSwirl_missed      = "#7D489C" #80008080
)


# Radial positions for each tool
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
        cex = rescale(vals[nonzero_idx], to = c(0.75,4), from = c(0, max_miss))
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
size_values <- c(1, unique_sizes[9],unique_sizes[14],unique_sizes[15],max_miss)
size_cex <- rescale(size_values, to = c(0.75, 3), from = c(0, max_miss))
legend("bottomleft",  # change position as needed
       legend = size_values,
       pt.cex = size_cex,
       pch = 16,
       cex = 0.75,
       bg = "white",
       y.intersp = c(1,1,1.2,1.2,1.7),
       x.intersp = c(1.5,1.5,1.4,1.5,1.5),
       inset = c(0, 0.2))

#title("Haplogroup L: Homoplasmic Variant Detection Gaps", line = -0.25)

dev.off()


#### Combine into 1 plot ####
homL <- readPNG("/vast/projects/bahlo_mtDNA/Analysis/homoplasmy/allduos/mitochondria_circos_dot_hapL.png")
homH <- readPNG("/vast/projects/bahlo_mtDNA/Analysis/homoplasmy/allduos/mitochondria_circos_dot_hapH.png")

homL_grob <- rasterGrob(homL, interpolate = TRUE)
homH_grob <- rasterGrob(homH, interpolate = TRUE)

p_homL <- wrap_elements(full=homL_grob)
p_homH <- wrap_elements(full=homH_grob)


final_plot <- 
  (p_homH | p_homL) +
  plot_annotation(tag_levels = "A") & theme(plot.tag=element_text(size=14,face="bold"))

ggsave(
  filename = "/vast/projects/bahlo_mtDNA/Analysis/final_figures/SupFig_3_final.tiff",
  plot = final_plot,
  device = "tiff",
  width = 10,          # inches (adjust as needed)
  height = 5,         # inches
  units = "in",
  dpi = 600,          # high resolution
  compression = "lzw"
)

