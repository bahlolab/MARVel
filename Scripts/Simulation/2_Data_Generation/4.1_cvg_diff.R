# rm(list = ls())
# ID <- "HG01881"
# VAF <- 10

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop('Usage: Rscript 4.1_cvg_diff.R "HG01881" 10')
}
ID <- args[1]
VAF <- as.integer(args[2])

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(vcfR)
})

work_dir <- paste0("/vast/projects/bahlo_mtDNA/Simulation/1000G/",ID,"/AF",VAF,"_DS2k/")

## --- Load inputs -------------------------------------------------------------

snv <- read.table(
  paste0(work_dir, "target.snv.", VAF, ".txt"), 
  header = FALSE, stringsAsFactors = FALSE)

indel <- read.table(
  paste0(work_dir, "target.indel.", VAF, ".txt"),
  sep = "\t", header = FALSE, fill = TRUE, stringsAsFactors = FALSE
)

incov  <- read.table(paste0(work_dir, "chrM.incover.txt"),  header = FALSE, stringsAsFactors = FALSE)
outcov <- read.table(paste0(work_dir, "chrM.outcover.txt"), header = FALSE, stringsAsFactors = FALSE)

## --- Coverage difference (incov - outcov) -----------------------------------

diff_df <- incov %>%
  transmute(pos = as.integer(V2), incov = as.numeric(V3)) %>%
  inner_join(
    outcov %>% transmute(pos = as.integer(V2), outcov = as.numeric(V3)),
    by = "pos"
  ) %>%
  mutate(diff = incov - outcov)

## --- Variant annotations -----------------------------------------------------

# SNVs: plot points on the actual diff curve
snv_df <- snv %>%
  transmute(pos = as.integer(V2)) %>%
  inner_join(diff_df %>% select(pos, diff), by = "pos") %>%
  rename(snv_diff = diff)

# INDELs: highlight regions (±10 bp padding)
indel_df <- indel %>%
  transmute(
    start = as.integer(V2) - 10L,
    end   = as.integer(V3) + 10L,
    type  = as.character(V5)
  ) %>%
  mutate(type = ifelse(is.na(type) | type == "", "INDEL", type))

## --- Dropped variant -----------------------------------------------------

# SNVs
snv_file <- paste0(work_dir, "/vcf/",ID,".chrM.40snv.af",VAF,".raw.addsnv.target.snv.",VAF,".vcf")
snv_vcf <- data.frame(read.vcfR(snv_file)@fix)
snv_drp <- snv_df[!(snv_df$pos %in% snv_vcf$POS),]

# INDELs
indel_file <- paste0(work_dir, "/vcf/",ID,".chrM.10indel.af",VAF,".raw.addindel.target.indel.",VAF,".vcf")
indel_vcf <- data.frame(read.vcfR(indel_file)@fix)
indel_drp <- indel_df %>%
  rowwise() %>%
  mutate(
    has_vcf = any(indel_vcf$CHROM == "chrM" &
                    indel_vcf$POS >= start &
                    indel_vcf$POS <= end)
  ) %>%
  ungroup() %>%
  filter(!has_vcf)

## --- Plot --------------------------------------------------------------------

cvg_diff <- ggplot(diff_df, aes(x = pos, y = diff)) +
  # INDEL regions (behind)
  geom_rect(
    data = indel_df,
    aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = type),
    inherit.aes = FALSE,
    alpha = 0.5
  ) +
  # dropped INDEL regions (purple outline)
  geom_rect(
    data = indel_drp,
    aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = NA,
    colour = "purple",
    linewidth = 0.5
  ) +
  # difference curve
  geom_line(colour = "grey30") +
  # zero line
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  # SNVs (on-curve)
  geom_point(
    data = snv_df,
    aes(x = pos, y = snv_diff),
    inherit.aes = FALSE,
    colour = "red",
    size = 1.2
  ) +
  geom_point(
    data = snv_drp,
    aes(x = pos, y = snv_diff),
    inherit.aes = FALSE,
    colour = "purple",
    size = 1.2
  ) +
  scale_fill_manual(
    values = c("INS" = "steelblue", "DEL" = "orange", "INDEL" = "grey70"),
    breaks = c("INS", "DEL")
  ) +
  labs(
    x = "chrM position",
    y = "Coverage difference (incov − outcov)",
    fill = "INDEL type",
    title = "Coverage difference with SNV and INDEL positions"
  ) +
  theme_bw()

ggsave(plot = cvg_diff, 
       filename = paste0(work_dir, "cvg.diff.png"), 
       width = 10, height = 3, 
       dpi = 300)


