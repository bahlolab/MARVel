#!/usr/bin/env Rscript

rm(list = ls())
setwd("/vast/scratch/users/wang.lo/1000G/PCA/")

suppressPackageStartupMessages({
  library(data.table)
})

# =====================
# User Inputs
# =====================
hg_file     <- "/vast/projects/bahlo_mtDNA/Results_1000G/haplogroups_all.tsv"
sample_file <- "mtPCA.sample_ids.txt"
scores_file <- "mtPCA.scores.1000G.PC1-20.txt"

thr <- 0.25
out_png <- "mtPCA_4panel_noL_splitPC1.png"

# =====================
# Load Haplogroup Table
# =====================
hg <- fread(hg_file)

# =====================
# Load Sample IDs
# =====================
sample_ids <- scan(sample_file, what="character", quiet=TRUE)
sample_ids <- sub("\\.final$", "", sample_ids)

# =====================
# Load PCA Scores
# =====================
scores_mat <- as.matrix(fread(scores_file, header=T))
colnames(scores_mat)[1:3] <- c("PC1","PC2","PC3")

scores_df <- as.data.frame(scores_mat)
scores_df$SAMPLE <- sample_ids

if (nrow(scores_df) != length(sample_ids)) {
  stop("Mismatch between scores rows and sample IDs")
}

# =====================
# Merge
# =====================
df <- merge(scores_df, hg, by="SAMPLE", all.x=TRUE)

cat("Total samples:", nrow(df), "\n")

# =====================
# Remove L*
# =====================
df2 <- df[!is.na(df$Haplogroup) & !grepl("^L", df$Haplogroup), ]
cat("After removing L*:", nrow(df2), "\n")
table(df2$Haplogroup)

# =====================
# Global Color Palette (scalable + colorblind-friendly)
# =====================
hgs_all <- sort(unique(df2$Haplogroup))

pal <- setNames(
  grDevices::hcl.colors(
    length(hgs_all),
    palette = "Dark 3"   # good contrast, avoids harsh red/green pairing
  ),
  hgs_all
)

df2$col <- adjustcolor(pal[as.character(df2$Haplogroup)], alpha.f = 0.7)
df2$col[is.na(df2$col)] <- adjustcolor("grey70", alpha.f = 0.7)

# =====================
# Split by PC1
# =====================
left  <- df2[df2$PC1 < thr, ]
right <- df2[df2$PC1 > thr, ]

cat("PC1 <", thr, ":", nrow(left), "\n")
cat("PC1 >", thr, ":", nrow(right), "\n")

set.seed(1)
if (nrow(left)  > 1) left  <- left[sample(nrow(left)), ]
if (nrow(right) > 1) right <- right[sample(nrow(right)), ]
table(left$Haplogroup)
table(right$Haplogroup)
# =====================
# Panel Plot Function
# =====================
plot_panel <- function(d, xvar, yvar, main_txt, panel_label) {
  plot(d[[xvar]], d[[yvar]],
       col = d$col, pch = 19, cex = 0.7,
       xlab = xvar, ylab = yvar,
       main = main_txt)
  
  mtext(panel_label, side = 3, adj = 0, line = 1, cex = 1.4, font = 2)
}

# =====================
# Legend prep (separate)
# =====================
get_legend <- function(d) {
  hgs <- sort(unique(d$Haplogroup))
  hg_n <- table(d$Haplogroup)
  list(
    hgs = hgs,
    labels = paste0(hgs, " (n=", as.integer(hg_n[hgs]), ")")
  )
}

leg_left  <- get_legend(left)   # M
leg_right <- get_legend(right)  # N

get_ncol <- function(n) {
  if (n <= 10) 1 else if (n <= 20) 2 else if (n <= 35) 3 else 4
}

# =====================
# Generate 4-panel Plot
# =====================
png(out_png, width = 2800, height = 2200, res = 220)

par(
  mfrow = c(2, 2),
  mar = c(5, 5, 4, 2),
  oma = c(0, 0, 0, 18),
  xpd = NA
)

# ---- A (M: PC1 vs PC2) ----
plot_panel(left, "PC1", "PC2", "Founder haplogroup M: PC1 vs PC2", "A")

legend(
  "topright",
  inset = c(-1.6, 0),
  legend = leg_left$labels,
  col = adjustcolor(pal[leg_left$hgs], alpha.f = 0.9),
  pch = 19,
  cex = 0.9,
  bty = "n",
  ncol = get_ncol(length(leg_left$hgs)),
  title = "M haplogroup"
)

# ---- B (M: PC2 vs PC3) ----
plot_panel(left, "PC2", "PC3", "Founder haplogroup M: PC2 vs PC3", "B")

# ---- C (N: PC1 vs PC2) ----
plot_panel(right, "PC1", "PC2", "Founder haplogroup N: PC1 vs PC2", "C")

legend(
  "topright",
  inset = c(-1.75, 0),
  legend = leg_right$labels,
  col = adjustcolor(pal[leg_right$hgs], alpha.f = 0.9),
  pch = 19,
  cex = 0.9,
  bty = "n",
  ncol = get_ncol(length(leg_right$hgs)),
  title = "N haplogroup"
)

# ---- D (N: PC2 vs PC3) ----
plot_panel(right, "PC2", "PC3", "Founder haplogroup N: PC2 vs PC3", "D")

dev.off()

cat("Saved:", out_png, "\n")

