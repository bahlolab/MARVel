############################################################
# STEP 7 — Compute 1000G mtDNA PCA Reference from GT TSV
# Input: 1000G.mtPCA.GT.tsv (from bcftools query)
# Output:
#   mtPCA.variant_ids.txt
#   mtPCA.sample_ids.txt
#   mtPCA.mu.txt
#   mtPCA.loadings.PC1-20.txt
#   mtPCA.scores.1000G.PC1-20.txt
#   mtPCA.eigenvalues.txt
############################################################
rm(list = ls())
setwd("/vast/scratch/users/wang.lo/1000G/PCA/")
suppressPackageStartupMessages({
  library(data.table)
})

gt_tsv <- "1000G.mtPCA.GT.withHeader.tsv"
nPC <- 20

cat("Reading GT TSV...\n")
dt <- fread(gt_tsv, sep = "\t", header = T, data.table = FALSE)

variant_ids <- dt[[1]]
gt <- as.matrix(dt[, -1, drop = FALSE])  # variants x samples (strings)

cat("Variants:", length(variant_ids), "\n")
cat("Samples :", ncol(gt), "\n")

# Optional: save variant order (contract for projection)
write.table(variant_ids, "mtPCA.variant_ids.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

############################################################
# Convert GT strings to 0/1
# We treat:
#   "1", "1/1", "1|1" => 1
#   everything else (0, 0/0, ., ./., 0/1, etc.) => 0
############################################################
cat("Converting GT to 0/1...\n")

X01 <- matrix(0L, nrow = nrow(gt), ncol = ncol(gt))

X01[gt %in% c("1", "1/1", "1|1")] <- 1L
# Note: if any heteroplasmic "0/1" slipped in, we keep it 0 by default.
# If you prefer presence/absence (any ALT => 1), uncomment:
# X01[gt %in% c("0/1","1/0","0|1","1|0")] <- 1L

# Transpose to samples x variants
X <- t(X01)

############################################################
# Center using 1000G means (mu)
############################################################
cat("Centering...\n")
mu <- colMeans(X)
write.table(mu, "mtPCA.mu.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

Xc <- sweep(X, 2, mu, "-")

############################################################
# PCA via SVD
############################################################
cat("Running SVD...\n")
sv <- svd(Xc)

k <- min(nPC, ncol(sv$v), nrow(sv$u))
scores <- sv$u[, 1:k, drop = FALSE] %*% diag(sv$d[1:k], nrow = k, ncol = k)
loadings <- sv$v[, 1:k, drop = FALSE]

# Eigenvalues (variance of each PC)
eigvals <- (sv$d^2) / (nrow(Xc) - 1)

write.table(scores, sprintf("mtPCA.scores.1000G.PC1-%d.txt", k),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(loadings, sprintf("mtPCA.loadings.PC1-%d.txt", k),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(eigvals[1:k], "mtPCA.eigenvalues.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

scores_df <- as.data.frame(scores)
colnames(scores_df) <- paste0("PC",1:20)

write.table(scores_df, "mtPCA.scores.1000G.PC1-20.txt",
            quote = FALSE, row.names = FALSE)

cat("Done.\n")

