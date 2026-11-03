### Supplementary Figure 2 ####
# MitoHPC freebayes vs mutect2 comparison

rm(list = ls())
library(ggplot2)
library(vcfR)
library(cowplot)
library(patchwork)
library(grid)
setwd("/vast/projects/bahlo_mtDNA/Test/MitoHPC_HG00405/")



#### Figure 2a ####
SS.mutect2.count <- read.table("./out_ss_mutect2/count.tab", header = T)
NSS.mutect2.count <- read.table("./out_nss_mutect2/count.tab", header = T)
identical(SS.mutect2.count[1,], NSS.mutect2.count)
## same! use all reads
### coverage
SS.mutect2.cvg <- read.table("./out_ss_mutect2/cvg.tab", header = T)
NSS.mutect2.cvg <- read.table("./out_nss_mutect2/cvg.tab", header = T)
rbind(SS.mutect2.cvg[1,], NSS.mutect2.cvg)

### SNV counts
SS.mutect2.tab3 <- read.table("./out_ss_mutect2/mutect2.03.tab", header = T)
NSS.mutect2.tab3 <- read.table("./out_nss_mutect2/mutect2.03.tab", header = T)
rbind(SS.mutect2.tab3[1,],NSS.mutect2.tab3)
# different at Hp, hp, Sp, sp
SS.mutect2.mutect2.tab3 <- read.table("./out_ss_mutect2/mutect2.mutect2.03.tab", header = T)
NSS.mutect2.mutect2.tab3 <- read.table("./out_nss_mutect2/mutect2.mutect2.03.tab", header = T)
rbind(SS.mutect2.mutect2.tab3[1,],NSS.mutect2.mutect2.tab3)
# different at Hp, hp, Sp, sp

### check the different variants

vcf <- read.vcfR("./out_ss_mutect2/mutect2.mutect2.03.vcf")
SS.mutect2.03.vcf <- vcf[vcf@gt[,2] == "HG00405.final",]
SS.mutect2.03.fix <- data.frame(SS.mutect2.03.vcf@fix)
# Extract AF value from INFO
SS.mutect2.03.fix$AF <- sapply(SS.mutect2.03.fix$INFO, function(info_string) {
  af_match <- regmatches(info_string, regexpr("AF=[^;]+", info_string))
  af_value <- sub("AF=", "", af_match)
  return(af_value)
})
SS.mutect2.03.fix$AF <- as.numeric(SS.mutect2.03.fix$AF)

NSS.mutect2.03.vcf <- read.vcfR("./out_nss_mutect2/mutect2.mutect2.03.vcf")
NSS.mutect2.03.fix <- data.frame(NSS.mutect2.03.vcf@fix)
NSS.mutect2.03.fix$AF <- sapply(NSS.mutect2.03.fix$INFO, function(info_string) {
  af_match <- regmatches(info_string, regexpr("AF=[^;]+", info_string))
  af_value <- sub("AF=", "", af_match)
  return(af_value)
})
NSS.mutect2.03.fix$AF <- as.numeric(NSS.mutect2.03.fix$AF)

AF.SS.NSS <- merge(SS.mutect2.03.fix[,-c(3,6,8)],NSS.mutect2.03.fix[,-c(3,6,8)], 
                   by = c("CHROM","POS","REF","ALT"), sort = F)
AF.SS.NSS$mean_af <- (AF.SS.NSS$AF.x + AF.SS.NSS$AF.y) / 2
AF.SS.NSS$diff_af <- AF.SS.NSS$AF.x - AF.SS.NSS$AF.y
AF.SS.NSS[AF.SS.NSS$diff_af != 0,]

p1 <- ggplot(AF.SS.NSS, aes(x = AF.x, y = AF.y)) +
  geom_jitter(width = 0.01, height = 0.01, color = "blue", size = 2, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size = 1) +
  labs(title = "Comparison of Allele Frequencies (SS vs NSS)",
       x = "SS VAF", y = "NSS VAF") +
  theme_minimal()

#### Figure 2b ####
p2 <- ggplot(AF.SS.NSS, aes(x = mean_af, y = diff_af)) +
  geom_point(color = "darkgreen", size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 1) +
  labs(title = "Bland-Altman Plot: VAF Comparison",
       x = "Mean VAF", y = "Difference (SS - NSS)") +
  theme_minimal()

#### Figure 2c ####
SS.freebayes.freebayes.tab3 <- read.table("./out_ss_freebayes/freebayes.freebayes.03.tab", header = T)
NSS.freebayes.freebayes.tab3 <- read.table("./out_nss_freebayes/freebayes.freebayes.03.tab", header = T)
rbind(SS.freebayes.freebayes.tab3[1,],NSS.freebayes.freebayes.tab3)
# different at Hp, hp, Sp, sp

### check the different variants
vcf <- read.vcfR("./out_ss_freebayes/freebayes.freebayes.03.vcf")
SS.freebayes.03.vcf <- vcf[vcf@gt[,2] == "HG00405.final",]
SS.freebayes.03.fix <- data.frame(SS.freebayes.03.vcf@fix)
# Extract AF value from INFO
SS.freebayes.03.fix$AF <- sapply(SS.freebayes.03.fix$INFO, function(info_string) {
  af_match <- regmatches(info_string, regexpr("AF=[^;]+", info_string))
  af_value <- sub("AF=", "", af_match)
  return(af_value)
})
SS.freebayes.03.fix$AF <- as.numeric(SS.freebayes.03.fix$AF)

NSS.freebayes.03.vcf <- read.vcfR("./out_nss_freebayes/freebayes.freebayes.03.vcf")
NSS.freebayes.03.fix <- data.frame(NSS.freebayes.03.vcf@fix)
NSS.freebayes.03.fix$AF <- sapply(NSS.freebayes.03.fix$INFO, function(info_string) {
  af_match <- regmatches(info_string, regexpr("AF=[^;]+", info_string))
  af_value <- sub("AF=", "", af_match)
  return(af_value)
})
NSS.freebayes.03.fix$AF <- as.numeric(NSS.freebayes.03.fix$AF)

AF.SS.NSS <- merge(SS.freebayes.03.fix[,-c(3,6:8)],NSS.freebayes.03.fix[,-c(3,6:8)], 
                   by = c("CHROM","POS","REF","ALT"), sort = F)
AF.SS.NSS$mean_af <- (AF.SS.NSS$AF.x + AF.SS.NSS$AF.y) / 2
AF.SS.NSS$diff_af <- AF.SS.NSS$AF.x - AF.SS.NSS$AF.y
AF.SS.NSS.diff <- AF.SS.NSS[AF.SS.NSS$diff_af != 0 & AF.SS.NSS$mean_af>0.1,]

p3 <- ggplot(AF.SS.NSS, aes(x = AF.x, y = AF.y)) +
  geom_jitter(width = 0.01, height = 0.01, color = "blue", size = 2, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size = 1) +
  labs(title = "Comparison of Allele Frequencies (SS vs NSS)",
       x = "SS VAF", y = "NSS VAF") +
  theme_minimal()

#### Figure 2d ####
p4 <- ggplot(AF.SS.NSS, aes(x = mean_af, y = diff_af)) +
  geom_point(color = "darkgreen", size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 1) +
  labs(title = "Bland-Altman Plot: VAF Comparison",
       x = "Mean VAF", y = "Difference (SS - NSS)") +
  theme_minimal() 


#### Combine ####
final_plot <- 
  (p1 | p2) /
  (p3 | p4) +
  plot_annotation(tag_levels = "A")

ggsave(
  filename = "/vast/projects/bahlo_mtDNA/Analysis/final_figures/SupFig_1.tiff",
  plot = final_plot,
  device = "tiff",
  width = 10,          # inches (adjust as needed)
  height = 10,         # inches
  units = "in",
  dpi = 600,          # high resolution
  compression = "lzw"
)
