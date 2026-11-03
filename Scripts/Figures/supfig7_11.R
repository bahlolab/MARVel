rm(list = ls())
library(png)
library(patchwork)
library(grid)
library(ggplot2)
setwd("/vast/projects/bahlo_mtDNA/Scripts/final/mtdna-report-plots-updated-1102")

####### sup fig 7 ############
p1 <- readPNG("1_median_depth.png")
p2 <- readPNG("1_mtdna_cn.png")
p3 <- readPNG("1_coverage_cv.png")
p4 <- readPNG("1_base_coverage.png")

p1_grob <- rasterGrob(p1, interpolate = TRUE)
p2_grob <- rasterGrob(p2, interpolate = TRUE)
p3_grob <- rasterGrob(p3, interpolate = TRUE)
p4_grob <- rasterGrob(p4, interpolate = TRUE)

plot1 <- wrap_elements(full=p1_grob)
plot2 <- wrap_elements(full=p2_grob)
plot3 <- wrap_elements(full=p3_grob)
plot4 <- wrap_elements(full=p4_grob)

final_plot <- 
  (plot1 | plot2) /
  (plot3 | plot4) +
  plot_annotation(tag_levels = "A") & theme(plot.tag=element_text(size=14,face="bold"))

ggsave(
  filename = "/vast/projects/bahlo_mtDNA/Analysis/final_figures/SuppFig7_coverage_4panel.tiff",
  plot = final_plot,
  device = "tiff",
  width = 10,          # inches (adjust as needed)
  height = 7,         # inches
  units = "in",
  dpi = 500,          # high resolution
  compression = "lzw"
)

####### sup fig 8 ############
rm(list = ls())
library(png)
library(patchwork)
library(grid)
library(ggplot2)
setwd("/vast/projects/bahlo_mtDNA/Scripts/final/mtdna-report-plots-updated-1102/")


p1 <- readPNG("2_contam.png")
p2 <- readPNG("2_contam_vs_depth.png")

p1_grob <- rasterGrob(p1, interpolate = TRUE)
p2_grob <- rasterGrob(p2, interpolate = TRUE)


plot1 <- wrap_elements(full=p1_grob)
plot2 <- wrap_elements(full=p2_grob)

final_plot <- 
  (plot1 | plot2) +
  plot_annotation(tag_levels = "A") & theme(plot.tag=element_text(size=14,face="bold"))

ggsave(
  filename = "/vast/projects/bahlo_mtDNA/Analysis/final_figures/SuppFig8_contam_2panel.tiff",
  plot = final_plot,
  device = "tiff",
  width = 10,          # inches (adjust as needed)
  height = 4,         # inches
  units = "in",
  dpi = 600,          # high resolution
  compression = "lzw"
)

####### sup fig 9 ############
rm(list = ls())
library(png)
library(patchwork)
library(grid)
library(ggplot2)

p1 <- readPNG("/vast/projects/bahlo_mtDNA/Scripts/final/mtdna-report-plots-updated-1102/3_hap_distrib.png")
p2 <- readPNG("/vast/projects/bahlo_mtDNA/Scripts/final/mtdna-report-plots-updated-1102/3_hap_assign_qual.png")

p1_grob <- rasterGrob(p1, interpolate = TRUE)
p2_grob <- rasterGrob(p2, interpolate = TRUE)


plot1 <- wrap_elements(full=p1_grob)
plot2 <- wrap_elements(full=p2_grob)

final_plot <- 
  (plot1 | plot2) +
  plot_annotation(tag_levels = "A") & theme(plot.tag=element_text(size=14,face="bold"))

ggsave(
  filename = "/vast/projects/bahlo_mtDNA/Analysis/final_figures/SuppFig9_haplogroup_2panel.tiff",
  plot = final_plot,
  device = "tiff",
  width = 10,          # inches (adjust as needed)
  height = 4.5,         # inches
  units = "in",
  dpi = 600,          # high resolution
  compression = "lzw"
)

####### sup fig 11 ############
setwd("/vast/projects/bahlo_mtDNA/Scripts/final/mtdna-report-plots-updated-1102/")
p1 <- readPNG("4_variant_types.png")
p2 <- readPNG("4_var_conseq.png")
p3 <- readPNG("4_var_freq.png")
p4 <- readPNG("4_homoplasmy_pos.png")
p5 <- readPNG("4_heteroplasmy_pos.png")
p6 <- readPNG("4_mss.png")

p1_grob <- rasterGrob(p1, interpolate = TRUE)
p2_grob <- rasterGrob(p2, interpolate = TRUE)
p3_grob <- rasterGrob(p3, interpolate = TRUE)
p4_grob <- rasterGrob(p4, interpolate = TRUE)
p5_grob <- rasterGrob(p5, interpolate = TRUE)
p6_grob <- rasterGrob(p6, interpolate = TRUE)

plot1 <- wrap_elements(full=p1_grob)
plot2 <- wrap_elements(full=p2_grob)
plot3 <- wrap_elements(full=p3_grob)
plot4 <- wrap_elements(full=p4_grob)
plot5 <- wrap_elements(full=p5_grob)
plot6 <- wrap_elements(full=p6_grob)

final_plot <- 
  (plot1 | plot2) /
  (plot3 | plot4) /
  (plot5 | plot6) +
  plot_annotation(tag_levels = "A") & theme(plot.tag=element_text(size=14,face="bold"))

ggsave(
  filename = "/vast/projects/bahlo_mtDNA/Analysis/final_figures/SuppFig11_variants_6panel.tiff",
  plot = final_plot,
  device = "tiff",
  width = 10,          # inches (adjust as needed)
  height = 12,         # inches
  units = "in",
  dpi = 600,          # high resolution
  compression = "lzw"
)

