############################################################
# PCA plot
############################################################
rm(list = ls())
setwd("/vast/scratch/users/wang.lo/1000G/PCA/")

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------
# Inputs
# -----------------------
hg_file     <- "/vast/projects/bahlo_mtDNA/Results_1000G/haplogroups_all.tsv"
sample_file <- "mtPCA.sample_ids.txt"
scores_file <- "mtPCA.scores.1000G.PC1-20.txt"
out_png     <- "mtPCA_phylogenetic_colors_tags.png"

# -----------------------
# Read data
# -----------------------
hg <- fread(hg_file)
scores_df <- fread(scores_file)

sample_ids <- scan(sample_file, what = "character", quiet = TRUE)
sample_ids <- sub("\\.final$", "", sample_ids)

# Ensure scores have PC1/PC2/PC3 column names
if (!all(c("PC1","PC2","PC3") %in% names(scores_df))) {
  setnames(scores_df, old = names(scores_df)[1:3], new = c("PC1","PC2","PC3"))
}

# Attach sample IDs (must match row order of scores file)
stopifnot(nrow(scores_df) == length(sample_ids))
scores_df[, SAMPLE := sample_ids]

# -----------------------
# Merge
# -----------------------
df <- merge(scores_df, hg, by = "SAMPLE", all.x = TRUE)

cat("Scores samples:", nrow(scores_df), "\n")
cat("Merged samples:", nrow(df), "\n")
cat("Missing Haplogroup rows:", sum(is.na(df$Haplogroup)), "\n")

# -----------------------
# Define haplogroup groups (single source of truth)
# -----------------------
grp_M <- c("M","C","D","G","Z","Q")
grp_N <- c("N","I","W","X","A","Y","S","R","H","HV","J","T","U","K","V","B","F","P")

assign_family <- function(hg) {
  if (is.na(hg)) return("Other")
  if (grepl("^L", hg)) return("L_family")
  if (hg %in% grp_M)   return("M_family")
  if (hg %in% grp_N)   return("N_family")
  "Other"
}

df$Family <- vapply(df$Haplogroup, assign_family, character(1))
print(table(df$Family, useNA = "ifany"))

# -----------------------
# Create label: L0–L6 + M/N/R tags
# -----------------------
df$label <- NA_character_

# L0–L6 exact labels
is_L0L6 <- !is.na(df$Haplogroup) & grepl("^L[0-6]$", df$Haplogroup)
df$label[is_L0L6] <- df$Haplogroup[is_L0L6]

# M/N/R family labels
df$label[df$Family == "M_family"] <- "M"
df$label[df$Family == "N_family"] <- "N"
df$label[df$Family == "R_family"] <- "R"

# -----------------------
# Colors: shades within each family (by Haplogroup)
# -----------------------
purple_base <- "#7B3294"  # L
green_base  <- "#1A9850"  # M family
blue_base <- "#2C7BB6"  # N family

make_shades <- function(hgs, light, dark) {
  hgs <- sort(unique(hgs))
  n <- length(hgs)
  if (n == 0) return(setNames(character(0), character(0)))
  cols <- colorRampPalette(c(light, dark))(n)
  setNames(cols, hgs)
}

col_map <- c(
  make_shades(df$Haplogroup[df$Family == "L_family"], "#E9D5F3", purple_base),
  make_shades(df$Haplogroup[df$Family == "M_family"], "#D5F5E3", green_base),
  make_shades(df$Haplogroup[df$Family == "N_family"], "#D6EAF8", blue_base)
)

# Anything else -> grey
all_hgs <- sort(unique(na.omit(df$Haplogroup)))
missing_in_map <- setdiff(all_hgs, names(col_map))
if (length(missing_in_map) > 0) {
  col_map <- c(col_map, setNames(rep("grey70", length(missing_in_map)), missing_in_map))
}

col_vec <- col_map[as.character(df$Haplogroup)]
col_vec[is.na(col_vec)] <- "grey80"

# Overplotting improvement
col_vec_alpha <- adjustcolor(col_vec, alpha.f = 0.55)
set.seed(1)
ord <- sample(seq_len(nrow(df)))

# -----------------------
# Plot + tags (centroid labels)
# -----------------------
png(out_png, width = 2400, height = 1200, res = 200)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))

# Single consistent font size
label_cex <- 1

# Helper to draw readable labels (white halo)
draw_label <- function(x, y, lab) {
  text(x, y, labels = lab, cex = label_cex, font = 2, col = "white")
  text(x, y, labels = lab, cex = label_cex, font = 2, col = "grey30")
}

# ---- PC1 vs PC2 ----
plot(df$PC1[ord], df$PC2[ord],
     col = col_vec_alpha[ord],
     pch = 19, cex = 0.6,
     xlab = "PC1", ylab = "PC2",
     main = "PC1 vs PC2")

mtext("A", side = 3, line = 1, adj = 0, cex = 1.5, font = 2)

if (any(!is.na(df$label))) {
  labs <- sort(unique(na.omit(df$label)))
  for (lab in labs) {
    idx <- which(df$label == lab)
    cx <- median(df$PC1[idx], na.rm = TRUE)
    cy <- median(df$PC2[idx], na.rm = TRUE)
    draw_label(cx, cy, lab)
  }
}

# ---- PC2 vs PC3 ----
plot(df$PC2[ord], df$PC3[ord],
     col = col_vec_alpha[ord],
     pch = 19, cex = 0.6,
     xlab = "PC2", ylab = "PC3",
     main = "PC2 vs PC3")

mtext("B", side = 3, line = 1, adj = 0, cex = 1.5, font = 2)

if (any(!is.na(df$label))) {
  labs <- sort(unique(na.omit(df$label)))
  for (lab in labs) {
    idx <- which(df$label == lab)
    cx <- median(df$PC2[idx], na.rm = TRUE)
    cy <- median(df$PC3[idx], na.rm = TRUE)
    draw_label(cx, cy, lab)
  }
}

dev.off()

cat("Saved:", out_png, "\n")

# -----------------------
# 3D PCA plot (interactive HTML)
# -----------------------
suppressPackageStartupMessages({
  library(plotly)
})

col_vec_alpha <- adjustcolor(col_vec, alpha.f = 0.85)
df$col3d <- col_vec_alpha  # reuse your colors with alpha

p3d <- plot_ly(
  data = df,
  x = ~PC1, y = ~PC2, z = ~PC3,
  type = "scatter3d", mode = "markers",
  marker = list(size = 4, color = ~col3d),
  text = ~paste0("SAMPLE: ", SAMPLE, "<br>HG: ", Haplogroup, "<br>Label: ", ifelse(is.na(label), "", label)),
  hoverinfo = "text"
) %>%
  layout(
    title = "3D PCA (PC1/PC2/PC3)",
    scene = list(
      xaxis = list(title = "PC1"),
      yaxis = list(title = "PC2"),
      zaxis = list(title = "PC3")
    )
  )

htmlwidgets::saveWidget(p3d, "mtPCA_3D_PC1_PC2_PC3.html", selfcontained = TRUE)
cat("Saved: mtPCA_3D_PC1_PC2_PC3.html\n")
