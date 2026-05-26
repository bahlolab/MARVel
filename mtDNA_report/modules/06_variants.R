# Module 06 — variant summary
#
# Expects in caller environment (from 01_load_data + config):
#   samples           — list with $meta and $vcf
#   variant_blacklist, filter_exclude, min_hom_vaf
#
# Produces:
#   var_plots$vaf_distribution      — VAF histogram, hom vs het
#   var_plots$consequence_breakdown — MLC_consq category bar chart
#   var_plots$mlc_score_dist        — MLC_score distribution histogram
#   var_plots$mss_distribution      — MSS (sum heteroplasmic MLC) per sample
#   var_plots$hom_per_sample_clade  — homoplasmy count histogram coloured by clade
#   var_plots$het_per_sample_clade  — heteroplasmy count histogram coloured by clade
#   var_plots$freq_spectrum         — unique variants by cohort frequency bin, coloured by hom/het class
#   var_plots$summary_pies          — 4 pie charts: bases, SNV type, variant class, variant calls

library(ggplot2)
library(plotly)
library(patchwork)

var_plots <- list()

# ── apply custom filter to all per-sample VCFs ───────────────────────────────
pass_variants <- function(v) {
  if (is.null(v) || nrow(v) == 0) return(v)
  v[!grepl(filter_exclude, v$FILTER) & !(v$POS %in% variant_blacklist), ]
}

# ── classify variants ─────────────────────────────────────────────────────────
classify_variant <- function(v) {
  v$zygosity <- ifelse(v$AF >= min_hom_vaf, "Homoplasmic", "Heteroplasmic")
  v
}

# ── build cohort-level variant table ─────────────────────────────────────────
# rbind_fill: like dplyr::bind_rows but base-R — pads missing columns with NA
rbind_fill <- function(lst) {
  lst       <- Filter(Negate(is.null), lst)
  if (length(lst) == 0) return(NULL)
  all_cols  <- unique(unlist(lapply(lst, names)))
  lst       <- lapply(lst, function(df) {
    missing <- setdiff(all_cols, names(df))
    df[missing] <- NA
    df[all_cols]
  })
  do.call(rbind, lst)
}

vcf_ok <- Filter(Negate(is.null), samples$vcf)

all_vars <- lapply(names(vcf_ok), function(id) {
  v <- pass_variants(vcf_ok[[id]])
  if (is.null(v) || nrow(v) == 0) return(NULL)
  v <- classify_variant(v)
  v$ID <- id
  v
})
all_vars <- rbind_fill(all_vars)

# For multi-allelic sites, keep one row per sample×position (highest AF allele)
# Then apply minimum AF threshold: heteroplasmic >= 0.03, homoplasmic >= min_hom_vaf
if (!is.null(all_vars) && nrow(all_vars) > 0) {
  all_vars <- all_vars[order(all_vars$ID, all_vars$POS, -all_vars$AF), ]
  all_vars <- all_vars[!duplicated(all_vars[, c("ID", "POS")]), ]
  all_vars <- all_vars[!is.na(all_vars$AF) & all_vars$AF >= 0.03, ]
}

meta    <- samples$meta
all_ids <- data.frame(ID = names(vcf_ok), stringsAsFactors = FALSE)

# ── 1. VAF distribution ───────────────────────────────────────────────────────
if (!is.null(all_vars) && nrow(all_vars) > 0) {
  
  var_plots$vaf_distribution <- ggplot(all_vars,
                                       aes(x = AF, fill = zygosity, colour = zygosity)) +
    geom_histogram(binwidth = 0.02, position = "stack", alpha = 0.8) +
    scale_fill_manual(values  = c("Homoplasmic"   = "steelblue",
                                  "Heteroplasmic" = "#E07B3A"), name = NULL) +
    scale_colour_manual(values = c("Homoplasmic"  = "steelblue",
                                   "Heteroplasmic"= "#E07B3A"), name = NULL) +
    geom_vline(xintercept = c(0.03, 0.95), linetype = "dashed",
               colour = "grey30", linewidth = 0.6) +
    annotate("text", x = 0.03, y = Inf, label = "0.03",
             hjust = -0.15, vjust = 1.5, size = 3, colour = "grey30") +
    annotate("text", x = 0.95, y = Inf, label = "0.95",
             hjust = -0.15, vjust = 1.5, size = 3, colour = "grey30") +
    labs(title = "VAF distribution (PASS variants)",
         x = "Allele frequency", y = "# variants") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "top")
  
  
  
  # ── per-unique-variant class (used in sections 5 and 9) ──────────────────
  has_hom <- aggregate(AF ~ POS + REF + ALT, data = all_vars,
                       FUN = function(x) any(x >= min_hom_vaf))
  names(has_hom)[4] <- "any_hom"
  
  has_het <- aggregate(AF ~ POS + REF + ALT, data = all_vars,
                       FUN = function(x) any(x < min_hom_vaf))
  names(has_het)[4] <- "any_het"
  
  class_levels <- c("heteroplasmic only", "homoplasmic only", "homoplasmic/heteroplasmic")
  class_cols   <- c("heteroplasmic only"       = "grey70",
                    "homoplasmic only"          = "#5B9BD5",
                    "homoplasmic/heteroplasmic" = "#1A3A6B")
  
  var_class_df <- merge(has_hom, has_het, by = c("POS", "REF", "ALT"))
  var_class_df$var_class <- factor(
    ifelse( var_class_df$any_hom &  var_class_df$any_het, "homoplasmic/heteroplasmic",
            ifelse( var_class_df$any_hom & !var_class_df$any_het, "homoplasmic only",
                    "heteroplasmic only")),
    levels = class_levels)
  
  # ── 5. consequence breakdown (MLC_consq) ─────────────────────────────────
  # severity ranking: pick most severe when multiple consequences are present
  consq_severity <- c("stop_gained", "frameshift_variant", "splice_site_variant",
                      "start_lost", "stop_lost", "missense_variant",
                      "stop_retained_variant", "synonymous_variant",
                      "tRNA", "rRNA", "intron_variant",
                      "upstream_gene_variant", "downstream_gene_variant",
                      "intergenic_variant")
  
  pick_worst <- function(x) {
    terms <- unique(trimws(strsplit(x, ",")[[1]]))
    hit   <- consq_severity[consq_severity %in% terms]
    if (length(hit)) hit[1] else terms[1]
  }
  
  if ("MLC_consq" %in% names(all_vars)) {
    consq_df <- all_vars[!is.na(all_vars$MLC_consq), ]
    consq_df <- consq_df[!duplicated(consq_df[, c("POS", "REF", "ALT")]), ]
    consq_df$MLC_consq <- sapply(consq_df$MLC_consq, pick_worst)
    
    # reclassify non_coding_transcript_exon_variant -> tRNA / rRNA using flag cols
    nc_idx <- consq_df$MLC_consq == "non_coding_transcript_exon_variant"
    if (any(nc_idx, na.rm = TRUE)) {
      consq_df$MLC_consq[nc_idx & !is.na(consq_df$TRN)] <- "tRNA"
      consq_df$MLC_consq[nc_idx & !is.na(consq_df$RNR)] <- "rRNA"
    }
    
    consq_df <- merge(consq_df, var_class_df[, c("POS", "REF", "ALT", "var_class")],
                      by = c("POS", "REF", "ALT"), all.x = TRUE)
    
    if (nrow(consq_df) > 0) {
      consq_tbl <- as.data.frame(table(Consequence = consq_df$MLC_consq,
                                       var_class   = consq_df$var_class))
      consq_tbl <- consq_tbl[consq_tbl$Freq > 0, ]
      
      # order consequences by total count
      consq_order <- aggregate(Freq ~ Consequence, data = consq_tbl, FUN = sum)
      consq_order <- consq_order[order(consq_order$Freq), ]
      consq_tbl$Consequence <- factor(consq_tbl$Consequence,
                                      levels = consq_order$Consequence)
      consq_tbl$var_class <- factor(consq_tbl$var_class, levels = class_levels)
      
      var_plots$consequence_breakdown <- ggplot(consq_tbl,
                                                aes(x = Freq, y = Consequence, fill = var_class)) +
        geom_col() +
        scale_fill_manual(values = class_cols, name = NULL) +
        scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
        labs(title = "Variant consequence breakdown (PASS)",
             x = "# unique variants", y = NULL) +
        theme_minimal(base_size = 11) +
        theme(panel.grid.major.y = element_blank(),
              legend.position = "top")
    }
  }
  
  # ── 6. MLC_score distribution ────────────────────────────────────────────
  if ("MLC_score" %in% names(all_vars)) {
    score_df <- all_vars[!is.na(all_vars$MLC_score), ]
    score_df$MLC_score <- as.numeric(score_df$MLC_score)
    score_df <- score_df[!is.na(score_df$MLC_score), ]
    
    if (nrow(score_df) > 0) {
      var_plots$mlc_score_dist <- ggplot(score_df, aes(x = MLC_score)) +
        geom_histogram(bins = 40, fill = "steelblue", colour = "white") +
        geom_vline(xintercept = median(score_df$MLC_score, na.rm = TRUE),
                   linetype = "dashed", colour = "firebrick", linewidth = 0.7) +
        labs(title = "MLC score distribution (PASS variants)",
             x = "MLC score", y = "# variants",
             caption = sprintf("median = %.3f", median(score_df$MLC_score, na.rm = TRUE))) +
        theme_minimal(base_size = 11)
    }
  }
  
  # ── 7. MSS (mitochondrial severity score) per sample ─────────────────────
  # MSS = sum of MLC_score across heteroplasmic PASS variants per sample
  if ("MLC_score" %in% names(all_vars)) {
    het_vars <- all_vars[all_vars$AF < min_hom_vaf, ]
    het_vars$MLC_score <- suppressWarnings(as.numeric(het_vars$MLC_score))
    het_vars <- het_vars[!is.na(het_vars$MLC_score), ]
    
    if (nrow(het_vars) > 0) {
      mss_df <- aggregate(MLC_score ~ ID, data = het_vars, FUN = sum)
      names(mss_df)[2] <- "MSS"
      mss_df <- merge(all_ids, mss_df, by = "ID", all.x = TRUE)
      mss_df$MSS[is.na(mss_df$MSS)] <- 0
      mss_df <- merge(mss_df, meta[, c("ID", "clade")], by = "ID", all.x = TRUE)
      mss_df$clade[is.na(mss_df$clade)] <- "Unknown"
      mss_df$clade <- factor(mss_df$clade,
                             levels = intersect(c("N/R", "M", "L", "Unknown"),
                                                unique(mss_df$clade)))
      
      set.seed(42)
      mss_df$jitter_y <- runif(nrow(mss_df), -0.4, 0.4)
      
      p_mss <- ggplot() +
        geom_histogram(data = mss_df, aes(x = MSS, fill = clade, colour = clade),
                       bins = 40, alpha = 0.5, position = "identity") +
        geom_point(data = mss_df,
                   aes(x = MSS, y = jitter_y,
                       text = paste0(ID, "<br>Clade: ", clade,
                                     "<br>MSS: ", round(MSS, 3))),
                   colour = "black", alpha = 0, size = 3) +
        geom_vline(xintercept = median(mss_df$MSS), linetype = "dashed",
                   colour = "grey30", linewidth = 0.6) +
        geom_vline(xintercept = mean(mss_df$MSS), linetype = "dotted",
                   colour = "grey30", linewidth = 0.6) +
        scale_fill_manual(values  = clade_cols, name = "Clade") +
        scale_colour_manual(values = clade_cols, name = "Clade") +
        labs(title = "Mitochondrial Local Constraint Score sum (MSS) per sample",
             x = "MSS (sum of heteroplasmic MLC scores)", y = "# samples") +
        theme_minimal(base_size = 11) +
        theme(legend.position = "top")
      
      var_plots$mss_distribution <- ggplotly(p_mss, tooltip = "text") |>
        layout(margin = list(t = 80, b = 80),
               annotations = list(
                 list(text = sprintf("mean = %.3f  |  median = %.3f  |  range = %.3f–%.3f",
                                     mean(mss_df$MSS), median(mss_df$MSS),
                                     min(mss_df$MSS),  max(mss_df$MSS)),
                      xref = "paper", yref = "paper",
                      x = 1, y = -0.11, xanchor = "right", yanchor = "top",
                      showarrow = FALSE,
                      font = list(size = 10, color = "grey40"))))
    }
  }
  
  # ── 8. homoplasmy count per sample coloured by clade ─────────────────────
  hom_df <- all_vars[all_vars$AF >= min_hom_vaf, ]
  hom_cnt <- as.data.frame(table(ID = hom_df$ID))
  names(hom_cnt)[2] <- "n_hom"
  
  hom_cnt <- merge(all_ids, hom_cnt, by = "ID", all.x = TRUE)
  hom_cnt$n_hom[is.na(hom_cnt$n_hom)] <- 0
  
  hom_cnt <- merge(hom_cnt, meta[, c("ID", "clade")], by = "ID", all.x = TRUE)
  hom_cnt$clade[is.na(hom_cnt$clade)] <- "Unknown"
  
  clade_levels <- intersect(c("N/R", "M", "L", "Unknown"), unique(hom_cnt$clade))
  hom_cnt$clade <- factor(hom_cnt$clade, levels = clade_levels)
  
  set.seed(42)
  hom_cnt$jitter_y <- runif(nrow(hom_cnt), -0.4, 0.4)
  
  p_hom <- ggplot() +
    geom_histogram(data = hom_cnt, aes(x = n_hom, fill = clade, colour = clade),
                   binwidth = 1, alpha = 0.5, position = "identity") +
    geom_point(data = hom_cnt,
               aes(x = n_hom, y = jitter_y,
                   text = paste0(ID, "<br>Clade: ", clade,
                                 "<br>Hom variants: ", n_hom)),
               colour = "black", alpha = 0, size = 3) +
    geom_vline(xintercept = median(hom_cnt$n_hom), linetype = "dashed",
               colour = "grey30", linewidth = 0.6) +
    geom_vline(xintercept = mean(hom_cnt$n_hom), linetype = "dotted",
               colour = "grey30", linewidth = 0.6) +
    scale_fill_manual(values  = clade_cols, name = "Clade") +
    scale_colour_manual(values = clade_cols, name = "Clade") +
    labs(title = "Homoplasmy positions / sample",
         x = "# homoplasmic variants", y = "# samples") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "top")
  
  var_plots$hom_per_sample_clade <- ggplotly(p_hom, tooltip = "text") |>
    layout(margin = list(t = 80, b = 80),
           annotations = list(
             list(text = sprintf("mean = %.1f  |  median = %.0f  |  range = %d–%d",
                                 mean(hom_cnt$n_hom), median(hom_cnt$n_hom),
                                 min(hom_cnt$n_hom),  max(hom_cnt$n_hom)),
                  xref = "paper", yref = "paper",
                  x = 1, y = -0.11, xanchor = "right", yanchor = "top",
                  showarrow = FALSE,
                  font = list(size = 10, color = "grey40"))))
  
  # ── 8. heteroplasmy count per sample coloured by clade ───────────────────
  het_df  <- all_vars[all_vars$AF < min_hom_vaf, ]
  het_cnt <- as.data.frame(table(ID = het_df$ID))
  names(het_cnt)[2] <- "n_het"
  
  het_cnt <- merge(all_ids, het_cnt, by = "ID", all.x = TRUE)
  het_cnt$n_het[is.na(het_cnt$n_het)] <- 0
  
  het_cnt <- merge(het_cnt, meta[, c("ID", "clade")], by = "ID", all.x = TRUE)
  het_cnt$clade[is.na(het_cnt$clade)] <- "Unknown"
  
  het_cnt$clade <- factor(het_cnt$clade,
                          levels = intersect(c("N/R", "M", "L", "Unknown"),
                                             unique(het_cnt$clade)))
  
  set.seed(42)
  het_cnt$jitter_y <- runif(nrow(het_cnt), -0.4, 0.4)
  
  p_het <- ggplot() +
    geom_histogram(data = het_cnt, aes(x = n_het, fill = clade, colour = clade),
                   binwidth = 1, alpha = 0.5, position = "identity") +
    geom_point(data = het_cnt,
               aes(x = n_het, y = jitter_y,
                   text = paste0(ID, "<br>Clade: ", clade,
                                 "<br>Het variants: ", n_het)),
               colour = "black", alpha = 0, size = 3) +
    geom_vline(xintercept = median(het_cnt$n_het), linetype = "dashed",
               colour = "grey30", linewidth = 0.6) +
    geom_vline(xintercept = mean(het_cnt$n_het), linetype = "dotted",
               colour = "grey30", linewidth = 0.6) +
    scale_fill_manual(values  = clade_cols, name = "Clade") +
    scale_colour_manual(values = clade_cols, name = "Clade") +
    labs(title = "Heteroplasmy positions / sample",
         x = "# heteroplasmic variants", y = "# samples") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "top")
  
  var_plots$het_per_sample_clade <- ggplotly(p_het, tooltip = "text") |>
    layout(margin = list(t = 80, b = 80),
           annotations = list(
             list(text = sprintf("mean = %.1f  |  median = %.0f  |  range = %d–%d",
                                 mean(het_cnt$n_het), median(het_cnt$n_het),
                                 min(het_cnt$n_het),  max(het_cnt$n_het)),
                  xref = "paper", yref = "paper",
                  x = 1, y = -0.11, xanchor = "right", yanchor = "top",
                  showarrow = FALSE,
                  font = list(size = 10, color = "grey40"))))
  
  # ── 9. variant frequency spectrum (site frequency) ───────────────────────
  n_samples <- length(vcf_ok)
  
  carrier_cnt <- aggregate(ID ~ POS + REF + ALT, data = all_vars, FUN = length)
  names(carrier_cnt)[4] <- "n_carriers"
  
  var_freq <- merge(carrier_cnt, var_class_df, by = c("POS", "REF", "ALT"))
  
  freq        <- var_freq$n_carriers / n_samples
  freq_levels <- c("singleton", "doubleton", "doubleton-1%", "1-5%", "5-10%", ">10%")
  var_freq$freq_cat <- factor(
    ifelse(var_freq$n_carriers == 1, "singleton",
           ifelse(var_freq$n_carriers == 2, "doubleton",
                  ifelse(freq <  0.01,             "doubleton-1%",
                         ifelse(freq <= 0.05,             "1-5%",
                                ifelse(freq <= 0.10,             "5-10%",  ">10%"))))),
    levels = freq_levels)
  
  # proportion of unique variants per frequency category, coloured by hom/het
  freq_tbl      <- as.data.frame(table(freq_cat  = var_freq$freq_cat,
                                       var_class = var_freq$var_class))
  total_vars    <- sum(freq_tbl$Freq)
  freq_tbl$prop <- freq_tbl$Freq / total_vars
  
  var_plots$freq_spectrum <- ggplot(freq_tbl, aes(x = freq_cat, y = prop, fill = var_class)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = class_cols, name = NULL) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(title = "Variant frequency spectrum",
         x = "Allele frequency", y = "Proportion of unique variants") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.x = element_blank(),
          legend.position = "top")
  
  # ── 10. summary pie charts ────────────────────────────────────────────────
  mt_genome_size <- 16569
  ti_pairs       <- c("AG", "GA", "CT", "TC")
  
  # pie helper
  make_pie <- function(df, fill_col, n_col, colours, title) {
    df$pct <- 100 * df[[n_col]] / sum(df[[n_col]])
    ggplot(df, aes(x = "", y = .data[[n_col]], fill = .data[[fill_col]])) +
      geom_col(width = 1, colour = "white", linewidth = 0.5) +
      coord_polar("y", start = 0) +
      geom_text(aes(label = sprintf("%.0f%%", pct)),
                position = position_stack(vjust = 0.5),
                size = 3.5, colour = "white", fontface = "bold") +
      scale_fill_manual(values = colours, name = NULL) +
      labs(title = title) +
      theme_void(base_size = 11) +
      theme(legend.position  = "right",
            plot.title = element_text(face = "bold", size = 10, hjust = 0.5))
  }
  
  # pie 1: mtDNA bases with/without variant
  n_var_bases <- length(unique(all_vars$POS))
  pie1_df <- data.frame(
    category = c("variant", "no variant"),
    n        = c(n_var_bases, mt_genome_size - n_var_bases)
  )
  pie1 <- make_pie(pie1_df, "category", "n",
                   c("variant" = "grey40", "no variant" = "grey80"),
                   sprintf("%d mtDNA bases", mt_genome_size))
  
  # pie 2: unique variants — SNV transition / transversion / indel
  uniq_vars <- all_vars[!duplicated(all_vars[, c("POS", "REF", "ALT")]),
                        c("POS", "REF", "ALT")]
  uniq_vars$snv_type <- ifelse(
    nchar(uniq_vars$REF) != 1 | nchar(uniq_vars$ALT) != 1, "indel",
    ifelse(paste0(uniq_vars$REF, uniq_vars$ALT) %in% ti_pairs,
           "SNV transition", "SNV transversion"))
  type_cnt <- as.data.frame(table(snv_type = uniq_vars$snv_type))
  names(type_cnt)[2] <- "n"
  type_cnt$snv_type <- factor(type_cnt$snv_type,
                              levels = c("SNV transition", "SNV transversion", "indel"))
  pie2 <- make_pie(type_cnt, "snv_type", "n",
                   c("SNV transition"   = "#2E7D32",
                     "SNV transversion" = "#81C784",
                     "indel"            = "grey70"),
                   sprintf("%s unique variants", format(nrow(uniq_vars), big.mark = ",")))
  
  # pie 3: unique variants — homoplasmic only / hom+het / heteroplasmic only
  class_cnt <- as.data.frame(table(var_class = var_class_df$var_class))
  names(class_cnt)[2] <- "n"
  pie3 <- make_pie(class_cnt, "var_class", "n",
                   c("heteroplasmic only"         = "grey70",
                     "homoplasmic only"            = "#5B9BD5",
                     "homoplasmic/heteroplasmic"   = "#1A3A6B"),
                   sprintf("%s unique variants", format(nrow(var_class_df), big.mark = ",")))
  
  # pie 4: variant calls — homoplasmic vs heteroplasmic
  calls_cnt <- as.data.frame(table(zygosity = all_vars$zygosity))
  names(calls_cnt)[2] <- "n"
  pie4 <- make_pie(calls_cnt, "zygosity", "n",
                   c("Heteroplasmic" = "grey70",
                     "Homoplasmic"   = "#5B9BD5"),
                   sprintf("%s variant calls", format(sum(calls_cnt$n), big.mark = ",")))
  
  var_plots$summary_pies <- (pie1 + pie2) / (pie3 + pie4)
  
}

message(sprintf("[06_variants] Variant plots built. %d PASS variants across %d samples.",
                if (!is.null(all_vars)) nrow(all_vars) else 0,
                length(vcf_ok)))