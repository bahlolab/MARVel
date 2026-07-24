# rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggstatsplot)   # ggwithinstats
  library(ggplot2)
})


# -----------------------------
# Helpers
# -----------------------------
### TP & sensitivity
# sensitivity = TP / (TP + FN)

TP_stat <- function(snv, indel, use_qc = TRUE, pass_label = "PASS"){
  
  tools <- c("server2","mity","mitoHPC","mtSwirl")
  
  n_snv <- nrow(snv)
  n_indel <- nrow(indel)
  n_total <- n_snv + n_indel
  
  get_fn <- function(df, suffix) {
    # suffix is "QC" or "VAF"
    sapply(tools, function(t) {
      col <- paste0(t, "_", suffix)
      if (!col %in% names(df)) return(nrow(df))  # column missing => all FN
      x <- df[[col]]
      
      if (suffix == "QC") {
        # FN = not PASS (including NA)
        sum(is.na(x) | x != pass_label)
      } else {
        # FN = missing VAF
        sum(is.na(x))
      }
    })
  }
  
  if (use_qc) {
    snv_FN   <- get_fn(snv, "QC")
    indel_FN <- get_fn(indel, "QC")
  } else {
    snv_FN   <- get_fn(snv, "VAF")
    indel_FN <- get_fn(indel, "VAF")
  }
  
  snv_TP   <- n_snv   - snv_FN
  indel_TP <- n_indel - indel_FN
  
  TP <- snv_TP + indel_TP
  
  snv_sensitivity <- round(snv_TP/n_snv, 3)
  indel_sensitivity <- round(indel_TP/n_indel, 3)
  
  sensitivity <- round(TP / n_total, 3)
  
  out <- data.frame(
    tool = tools,
    snv_TP = snv_TP,
    snv_sensitivity = if (n_snv > 0) round(snv_TP / n_snv, 3) else NA_real_,
    indel_TP = indel_TP,
    indel_sensitivity = if (n_indel > 0) round(indel_TP / n_indel, 3) else NA_real_,
    TP = TP,
    sensitivity = if (n_total > 0) round(TP / n_total, 3) else NA_real_
  )
  
  out
}

### FP

FP_stat <- function(FP_list){
  
  tools = c("server2","mity","mitoHPC","mtSwirl")
  
  # Expect FP_list order: server2, mity, mitoHPC, mtSwirl
  names(FP_list) <- tools
  
  safe_nrow <- function(x) {
    if (is.null(x)) return(0L)
    if (is.data.frame(x)) return(as.integer(nrow(x)))
    if (is.list(x)) {
      # sum rows of any data.frames inside
      return(as.integer(sum(vapply(x, function(z) if (is.data.frame(z)) nrow(z) else 0L, integer(1)))))
    }
    0L
  }
  
  FP <- vapply(FP_list, safe_nrow, integer(1))
  
  FP_in <- sapply(FP_list, function(x) {
    if (is.null(x) || !("FP_source" %in% names(x))) return(0L)
    sum(x$FP_source != "unknown", na.rm = TRUE)
  })
  
  # tool-specific thresholds
  FP_thred <- c(
    server2 = if (!is.null(FP_list$server2) && "VariantLevel" %in% names(FP_list$server2))
      sum(FP_list$server2$VariantLevel >= 0.02, na.rm = TRUE) else 0L,
    mity = if (!is.null(FP_list$mity) && "VAF" %in% names(FP_list$mity))
      sum(FP_list$mity$VAF >= 0.01, na.rm = TRUE) else 0L,
    mitoHPC = if (!is.null(FP_list$mitoHPC) && "VAF" %in% names(FP_list$mitoHPC))
      sum(FP_list$mitoHPC$VAF >= 0.03, na.rm = TRUE) else 0L,
    mtSwirl = if (!is.null(FP_list$mtSwirl) && "VAF" %in% names(FP_list$mtSwirl))
      sum(FP_list$mtSwirl$VAF >= 0.03, na.rm = TRUE) else 0L
  )
  
  out <- data.frame(
    tool = tools,
    FP = as.integer(FP),
    FP_in = as.integer(FP_in),
    FP_thred = as.integer(FP_thred[tools])
  )
  
  out
}

VAF_diff <- function(data, out=NULL){
  tools <- c("server2", "mity", "mitoHPC", "mtSwirl")
  
  # compute diffs in wide
  for (t in tools) {
    vaf_col  <- paste0(t, "_VAF")
    diff_col <- paste0(t, "_diff")
    
    v <- data[[vaf_col]]
    v[is.na(v)] <- 0
    data[[diff_col]] <- pmax(0, data$minVAF - v, v - data$maxVAF)
  }
  
  # summary
  diff_cols <- paste0(tools, "_diff")
  diff_mean <- round(colMeans(data[, diff_cols, drop = FALSE], na.rm = TRUE), 3)
  stat <- data.frame(tool = tools, diff_mean = as.numeric(diff_mean))

  if (!is.null(out)) {
    diff_long <- data %>%
      select(all_of(diff_cols)) %>%
      pivot_longer(everything(), names_to = "tool", values_to = "diff") %>%
      mutate(
        tool = sub("_diff$", "", tool),
        tool = factor(tool, levels = tools),
        diff = as.numeric(diff)
      )

    levels(diff_long$tool) <- c("mtDNA-server-2","mity", "MitoHPC", "mtSwirl")
    
    p <- ggwithinstats(
      data = diff_long,
      x = tool,
      y = diff,
      type = "np",
      plot.type = "violin",
      pairwise.comparisons = TRUE,
      p.adjust.method = "BH",
      centrality.type = "median",
      point.args = list(size = 1.2, alpha = 0.6),
      violin.args = list(alpha = 0.4),
      xlab = NULL,
      ylab = "Out-of-range distance",
      title = "VAF out-of-range distance by tool"
    ) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
    
    ggsave(filename = paste0(out, "_VAF_cmp.png"), plot = p,
           width = 9, height = 7, dpi = 300)
  }
  
  stat
}

################################################################################################
# -----------------------------
# Main
# -----------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop('Usage: Rscript ${scp_dir}/5_stats.R "HG01881"')
}
ID <- args[1]

base_dir <- paste0("/vast/projects/bahlo_mtDNA/Simulation/1000G/", ID)

results_all <- list()
VAFs <- c(1, 5, 10, 50)

for (VAF in VAFs) {
  vdir <- paste0(base_dir, "/AF", VAF, "_DS2k")
  
  snv   <- read.table(paste0(vdir, "/all_tools_snv.tsv"), header = TRUE, sep = "\t")
  indel <- read.table(paste0(vdir, "/all_tools_indel.tsv"), header = TRUE, sep = "\t")
  
  FP_list <- readRDS(paste0(vdir, "/FP_list.rds"))
  
  TP <- TP_stat(snv, indel, use_qc = TRUE)
  FP <- FP_stat(FP_list)
  
  res <- merge(TP, FP, by = "tool", sort = FALSE)
  
  # precision = TP / (TP + FP)
  res$precision <- round(res$TP / (res$TP + res$FP), 3)
  res$precision_thred <- round(res$TP / (res$TP + res$FP_thred), 3)
  
  # F1 = 2x(precision x sensitivity)/(precision + sensitivity)
  res$F1 <- round(2 * (res$precision * res$sensitivity) / (res$precision + res$sensitivity), 3)
  res$F1_thred <- round(2 * (res$precision_thred * res$sensitivity) /
                          (res$precision_thred + res$sensitivity), 3)
  
  VAF_snv   <- VAF_diff(snv,   out = paste0(vdir, "/SNV"))
  VAF_indel <- VAF_diff(indel, out = paste0(vdir, "/INDEL"))
  
  VAF_stat <- merge(VAF_snv, VAF_indel, by = "tool", suffixes = c("_SNV", "_INDEL"), sort = FALSE)
  
  res2 <- merge(res, VAF_stat, by = "tool", sort = FALSE)
  res2$VAF <- VAF
  
  results_all[[as.character(VAF)]] <- res2
  
}

stats_df <- bind_rows(results_all)
write.table(stats_df, file = paste0(base_dir, "/stats.tsv"),
            row.names = FALSE, quote = FALSE, sep = "\t")

