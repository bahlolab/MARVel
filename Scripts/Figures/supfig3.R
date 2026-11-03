### Supplementary Figure 3 ####
# Variants from non-subsampled Mutect2 and FreeBayes runs in MitoHPC.
rm(list=ls())
setwd("/vast/projects/bahlo_mtDNA/Test/MitoHPC_HG00405/")
library(vcfR)
library(ggplot2)
library(ggvenn)
library(patchwork)
library(grid)

#### Figure 3a ####

# check the different variants
freebayes.vcf <- read.vcfR("./out_nss_freebayes/freebayes.freebayes.03.vcf")
freebayes.fix <- data.frame(freebayes.vcf@fix)
# Extract AF value from INFO
freebayes.fix$AF <- sapply(freebayes.fix$INFO, function(info_string) {
  af_match <- regmatches(info_string, regexpr("AF=[^;]+", info_string))
  af_value <- sub("AF=", "", af_match)
  return(af_value)
})
freebayes.fix$AF <- as.numeric(freebayes.fix$AF)

mutect2.vcf <- read.vcfR("./out_nss_mutect2/mutect2.mutect2.03.vcf")
mutect2.fix <- data.frame(mutect2.vcf@fix)
mutect2.fix$AF <- sapply(mutect2.fix$INFO, function(info_string) {
  af_match <- regmatches(info_string, regexpr("AF=[^;]+", info_string))
  af_value <- sub("AF=", "", af_match)
  return(af_value)
})
mutect2.fix$AF <- as.numeric(mutect2.fix$AF)

AF.mutect2.freebayes <- merge(mutect2.fix[,-c(3,6:8)],freebayes.fix[,-c(3,6:8)], 
                              by = c("CHROM","POS","REF","ALT"), sort = F)
AF.mutect2.freebayes$mean_af <- (AF.mutect2.freebayes$AF.x + AF.mutect2.freebayes$AF.y) / 2
AF.mutect2.freebayes$diff_af <- AF.mutect2.freebayes$AF.x - AF.mutect2.freebayes$AF.y
AF.mutect2.freebayes[AF.mutect2.freebayes$diff_af != 0,] # 7

freebayes.fix$id <- paste(freebayes.fix$POS,freebayes.fix$REF, freebayes.fix$ALT,sep = "_")
mutect2.fix$id <- paste(mutect2.fix$POS,mutect2.fix$REF, mutect2.fix$ALT,sep = "_")

# Define sets
data <- list(
  mutect2 = mutect2.fix$id,
  freebayes = freebayes.fix$id
)
# Plot Venn diagram
p1 <- ggvenn(data, 
            auto_scale = TRUE,
            text_size = 0) + 
  annotate("text", x = -0.96, y = 0, 
           label = "3\n(1.9%)", size = 4) +
  
  annotate("text", x = -0.2, y = 0, 
           label = "49\n(31.2%)", size = 4) +
  
  annotate("text", x = 0.8, y = 0, 
           label = "105\n(66.9%)", size = 4) 
  # + theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))


#### Figure 3b ####
AF.mutect2.freebayes.all <- merge(mutect2.fix[,-c(3,6:8)],freebayes.fix[,-c(3,6:8)], 
                                  by = c("CHROM","POS","REF","ALT"), sort = F, all = T)
only_freebayes <- AF.mutect2.freebayes.all[(is.na(AF.mutect2.freebayes.all$AF.x) & !(is.na(AF.mutect2.freebayes.all$AF.y))),]
only_mutect2 <- AF.mutect2.freebayes.all[!is.na(AF.mutect2.freebayes.all$AF.x) & (is.na(AF.mutect2.freebayes.all$AF.y)),]

p2 <- ggplot(only_freebayes, aes(x = AF.y)) +
  geom_histogram(binwidth = 0.001, alpha = 0.6, position = "identity") +
  labs(title = "Histogram of AF.FreeBayes-only Variants",
       x = "Allele Frequency (AF)",
       y = "Count") +
  theme_minimal()

#### Combine into 1 plot ####
final_plot <- 
  (p1 | p2) + 
  plot_annotation(tag_levels = "A")

ggsave(
  filename = "/vast/projects/bahlo_mtDNA/Analysis/final_figures/SupFig_2.tiff",
  plot = final_plot,
  device = "tiff",
  width = 10,          # inches (adjust as needed)
  height = 5,         # inches
  units = "in",
  dpi = 600,          # high resolution
  compression = "lzw"
)