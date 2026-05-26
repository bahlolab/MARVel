# Module 05 — mtDNA PCA projection
#
# Expects in caller environment (from 01_load_data + config):
#   samples       — list with $meta and $vcf
#   pca_dir       — path to 1000G PCA reference files (set in config.R)
#   cohort_ids, filter_exclude, variant_blacklist
#   clade_cols    — named colour vector (defined in 04_haplogroup.R)
#
# 1000G PCA files:
#   mtPCA.variant_ids.txt        — chrM:POS:REF:ALT per line
#   mtPCA.mu.txt                 — per-variant mean AF (centering)
#   mtPCA.loadings.PC1-20.txt    — 317 variants x 20 PCs
#   mtPCA.eigenvalues.txt        — eigenvalues
#   1000G_mtPCApanel.bed/bim/fam — 1000G genotypes in PLINK format
#
# Produces:
#   pca_plots$pc1_pc2 — PC1 vs PC2, 1000G grey + cohort coloured by clade
#   pca_plots$pc2_pc3 — PC2 vs PC3
#   pca_plots$scree   — scree plot of variance explained

library(ggplot2)
library(plotly)
library(BEDMatrix)

pca_plots <- list()

# ── load PCA reference data ───────────────────────────────────────────────────
variant_ids <- readLines(file.path(pca_dir, "mtPCA.variant_ids.txt"))
mu          <- scan(file.path(pca_dir, "mtPCA.mu.txt"), quiet = TRUE)
loadings    <- as.matrix(read.table(file.path(pca_dir, "mtPCA.loadings.PC1-20.txt")))
eigenvalues <- scan(file.path(pca_dir, "mtPCA.eigenvalues.txt"), quiet = TRUE)

# parse variant IDs: chrM:POS:REF:ALT
var_split <- strsplit(variant_ids, ":")
var_df    <- data.frame(
  id  = variant_ids,
  pos = as.integer(sapply(var_split, `[`, 2)),
  ref = sapply(var_split, `[`, 3),
  alt = sapply(var_split, `[`, 4),
  stringsAsFactors = FALSE
)

# ── project 1000G samples ─────────────────────────────────────────────────────
bed      <- BEDMatrix(file.path(pca_dir, "1000G_mtPCApanel"), simple_names = TRUE)
geno_1kg <- as.matrix(bed)
geno_1kg[is.na(geno_1kg)] <- 0

centered_1kg <- sweep(geno_1kg, 2, mu, "-")
scores_1kg   <- centered_1kg %*% loadings
colnames(scores_1kg) <- paste0("PC", seq_len(ncol(scores_1kg)))

df_1kg       <- as.data.frame(scores_1kg)
df_1kg$group <- "1000G"

# ── project cohort samples ────────────────────────────────────────────────────
get_dosage <- function(vcf, var_df) {
  dosage <- rep(0L, nrow(var_df))
  if (is.null(vcf) || nrow(vcf) == 0) return(dosage)
  vcf <- vcf[!grepl(filter_exclude, vcf$FILTER) &
               !(vcf$POS %in% variant_blacklist), ]
  for (i in seq_len(nrow(var_df))) {
    hit <- vcf$POS == var_df$pos[i] &
      vcf$REF == var_df$ref[i] &
      vcf$ALT == var_df$alt[i]
    if (any(hit, na.rm = TRUE))
      dosage[i] <- 2L
  }
  dosage
}

meta      <- samples$meta
clade_map <- setNames(meta$clade,      meta$ID)
hg_map    <- setNames(meta$haplogroup, meta$ID)

cohort_scores <- lapply(cohort_ids, function(id) {
  dos      <- get_dosage(samples$vcf[[id]], var_df)
  centered <- dos - mu
  scores   <- as.numeric(centered %*% loadings)
  data.frame(as.list(setNames(scores, paste0("PC", seq_along(scores)))),
             ID         = id,
             clade      = ifelse(is.na(clade_map[id]), "Unknown", clade_map[id]),
             haplogroup = ifelse(is.na(hg_map[id]),    "Unknown", hg_map[id]),
             group      = "Cohort",
             stringsAsFactors = FALSE)
})
df_cohort <- do.call(rbind, cohort_scores)
df_cohort$tooltip <- paste0(df_cohort$ID,
                            "<br>Haplogroup: ", df_cohort$haplogroup,
                            "<br>Clade: ",      df_cohort$clade)

# ── scree plot ────────────────────────────────────────────────────────────────
var_exp  <- 100 * eigenvalues / sum(eigenvalues)
scree_df <- data.frame(PC = seq_along(var_exp), var_exp = var_exp)

pca_plots$scree <- ggplot(scree_df, aes(x = PC, y = var_exp)) +
  geom_col(fill = "steelblue", width = 0.7) +
  geom_line(colour = "firebrick", linewidth = 0.7) +
  geom_point(colour = "firebrick", size = 2) +
  scale_x_continuous(breaks = scree_df$PC) +
  labs(title = "Scree plot — mtDNA PCA",
       x = "PC", y = "% variance explained") +
  theme_minimal(base_size = 11)

# ── PC1 vs PC2 ────────────────────────────────────────────────────────────────
p12 <- ggplot() +
  geom_point(data = df_1kg,
             aes(x = PC1, y = PC2),
             colour = "grey80", size = 1, alpha = 0.5) +
  geom_point(data = df_cohort,
             aes(x = PC1, y = PC2, colour = clade, text = tooltip),
             size = 2, alpha = 0.8) +
  scale_colour_manual(values = clade_cols, name = "Clade") +
  labs(title = "mtDNA PCA — cohort projected onto 1000G space",
       x = sprintf("PC1 (%.1f%%)", var_exp[1]),
       y = sprintf("PC2 (%.1f%%)", var_exp[2]),
       caption = sprintf("1000G n=%d (grey)  |  Cohort n=%d (colour)",
                         nrow(df_1kg), nrow(df_cohort))) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")

pca_plots$pc1_pc2 <- ggplotly(p12, tooltip = "text", width = 650, height = 650) |>
  layout(margin = list(t = 80))

# ── PC2 vs PC3 ────────────────────────────────────────────────────────────────
p23 <- ggplot() +
  geom_point(data = df_1kg,
             aes(x = PC2, y = PC3),
             colour = "grey80", size = 1, alpha = 0.5) +
  geom_point(data = df_cohort,
             aes(x = PC2, y = PC3, colour = clade, text = tooltip),
             size = 2, alpha = 0.8) +
  scale_colour_manual(values = clade_cols, name = "Clade") +
  labs(title = "mtDNA PCA — PC2 vs PC3",
       x = sprintf("PC2 (%.1f%%)", var_exp[2]),
       y = sprintf("PC3 (%.1f%%)", var_exp[3])) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")

pca_plots$pc2_pc3 <- ggplotly(p23, tooltip = "text", width = 650, height = 650) |>
  layout(margin = list(t = 80))

message(sprintf("[05_pca] Projected %d cohort samples onto 1000G PCA space.",
                nrow(df_cohort)))