### Supplementary Figure 1 ####
# Heteroplasmy Upset Plot 
rm(list = ls())
library(UpSetR)
library(dplyr)
library(tidyr)
library(stringr)

setwd("/vast/projects/bahlo_mtDNA/Analysis/heteroplasmy/allduos")
data <- read.table("trio604_vaf03_het_newfilters.tsv", header = TRUE, sep = "\t")

variant_wide <- data %>%
  group_by(SAMPLE, VAR, TOOL) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = TOOL, values_from = count, values_fill = 0)

variant_wide$pos <- as.numeric(str_extract(variant_wide$VAR, "^[0-9]+"))

upsetdf <- variant_wide %>% 
  select(`mtDNA Server2`,mity,mitoHPC,mtSwirl)  %>%  # ensure 0/1
  as.data.frame()

png("/vast/projects/bahlo_mtDNA/Analysis/heteroplasmy/allduos/upset_plot.png", width = 4000, height = 4000, res = 600)

upset(upsetdf,sets = c("mtDNA Server2","mity","mitoHPC","mtSwirl"),order.by = "freq",
      text.scale=c(1.3, 1.3, 1, 1, 1.3, 1.2),point.size = 2.5,keep.order = TRUE,
      empty.intersections = "on")

dev.off()
