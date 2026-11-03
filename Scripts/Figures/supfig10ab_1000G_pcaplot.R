rm(list = ls())

library(tidyverse)
library(ggrepel)
library(data.table)
library(patchwork)
library(grid)


hg_file     <- "/vast/projects/bahlo_mtDNA/Results_1000G/haplogroups_all.tsv"
sample_file <- "/vast/projects/bahlo_mtDNA/Analysis/PCA/mtPCA.sample_ids.txt"
scores_file <- "/vast/projects/bahlo_mtDNA/Analysis/PCA/mtPCA.scores.1000G.PC1-20.txt"
out_png     <- "/vast/projects/bahlo_mtDNA/Analysis/final_figures/supfig_pca_1000G.tiff"
eigenvalues <- "/vast/projects/bahlo_mtDNA/Analysis/PCA/mtPCA.eigenvalues.txt"
# -----------------------
# Read data
# -----------------------
hg <- fread(hg_file)
scores_df <- fread(scores_file)
eigenval <- read.table(eigenvalues,header = FALSE)
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

centroids <- df %>%
  filter(!is.na(label)) %>%
  group_by(label) %>%
  summarise(
    PC1 = median(PC1),
    PC2 = median(PC2),
    PC3 = median(PC3),
    .groups = "drop"
  )

# -----------------------
# Theme
# -----------------------
pca_theme <- theme_minimal(base_size = 11) +
  theme(legend.position = "top",
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12.5),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12.5),
    plot.title = element_text(size = 14, face = "bold",hjust = 0.5)
  )

# -----------------------
# PC1 vs PC2 plot
# -----------------------
p1 <- ggplot(df, aes(x = PC1, y = PC2)) +
  
  # points (use your colour vector directly)
  geom_point(
    aes(colour = Haplogroup),
    size = 1.3,
    alpha = 0.6
  ) +
  
  scale_colour_manual(values = col_map, guide = "none") +
  
  # centroids
  geom_text_repel(
    data = centroids,
    aes(x = PC1, y = PC2, label = label),
    size = 4,
    fontface = "bold",
    box.padding = 0.3,
    point.padding = 0.25,
    max.overlaps = Inf
  ) +
  
  labs(
    title = "PC1 vs PC2",
    x = sprintf("PC1 (%.1f%%)", eigenval[1,1]/sum(eigenval)*100),
    y = sprintf("PC2 (%.1f%%)", eigenval[2,1]/sum(eigenval)*100)
  ) +
  
  pca_theme



# -----------------------
# PC2 vs PC3 plot
# -----------------------

p2 <- ggplot(df, aes(x = PC2, y = PC3)) +
  
  # points (use your colour vector directly)
  geom_point(
    aes(colour = Haplogroup),
    size = 1.3,
    alpha = 0.6
  ) +
  
  scale_colour_manual(values = col_map, guide = "none") +
  
  # centroids
  geom_text_repel(
    data = centroids,
    aes(x = PC2, y = PC3, label = label),
    size = 4,
    fontface = "bold",
    box.padding = 0.25,
    point.padding = 0.2,
    max.overlaps = Inf
  ) +
  
  labs(
    title = "PC2 vs PC3",
    x = sprintf("PC2 (%.1f%%)", eigenval[2,1]/sum(eigenval)*100),
    y = sprintf("PC3 (%.1f%%)", eigenval[3,1]/sum(eigenval)*100)
  ) +
  
  pca_theme


final_plot <- (p1 | p2) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 14, face = "bold"))

final_plot


ggsave(
  filename = "/vast/projects/bahlo_mtDNA/Analysis/final_figures/supfig_pca_1000G.tiff",
  plot = final_plot,
  device = "tiff",
  width = 12,          # inches (adjust as needed)
  height = 6,         # inches
  units = "in",
  dpi = 600,          # high resolution
  compression = "lzw"
)


####### RUN AFTER BOTH ATAXIA AND 1000G PLOTS CREATED ##########
library(tiff)
library(grid)
library(gridExtra)

img1 <- readTIFF("/vast/projects/bahlo_mtDNA/Analysis/final_figures/supfig_pca_1000G.tiff", native = FALSE)
img2 <- readTIFF("/vast/projects/bahlo_mtDNA/Analysis/final_figures/supfig_pca_ataxia.tiff", native = FALSE)

g1 <- rasterGrob(img1, interpolate = TRUE)
g2 <- rasterGrob(img2, interpolate = TRUE)

combined <- grid.arrange(g1, g2, ncol = 1)

# Save combined image
tiff("/vast/projects/bahlo_mtDNA/Analysis/final_figures/SuppFig10_pca_4panel.tiff", width = 12, height = 12, res = 600,units = "in",compression = "lzw")
grid.arrange(g1, g2, ncol = 1)
dev.off()