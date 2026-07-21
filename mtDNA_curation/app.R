library(shiny)
library(shinythemes)
library(DT)
library(data.table)
suppressPackageStartupMessages(library(tidyverse))

# ---------------------------------------------------------------------------
# Database paths — edit these when deploying to HPC
# ---------------------------------------------------------------------------
# Set MTDNA_DB_DIR environment variable to override (e.g. on HPC).
# Otherwise defaults to database/ sibling folder next to this app.
# DB_DIR <- "/stornext/Bioinf/data/lab_bahlo/ref_db/human/mtDNA/"
APP_DIR <- tryCatch(
  normalizePath(dirname(sys.frame(1)$ofile), mustWork = FALSE),
  error = function(e) getwd()
)
DB_DIR <- Sys.getenv("MTDNA_DB_DIR",
                     unset = file.path(APP_DIR, "database"))

GNOMAD_FILE  <- file.path(DB_DIR, "gnomad", "gnomad_chrM_hap_AF.tsv.gz")  # from prepare_databases.R
MITOTIP_FILE <- file.path(DB_DIR, "MitoTIP", "mitotip_scores_27_04_2020_VAR.txt")

# MITOMAP filenames include download dates — discover by pattern rather than hardcoding
find_mitomap <- function(db_dir, pattern) {
  hits <- list.files(file.path(db_dir, "MITOMAP"), pattern = pattern,
                     full.names = TRUE, ignore.case = TRUE)
  if (length(hits) == 0) return(NULL)
  if (length(hits) > 1) message("Multiple MITOMAP files matched '", pattern, "', using: ", hits[1])
  hits[1]
}
MITOMAP_CDS  <- find_mitomap(DB_DIR, "MutationsCodingControl.*_VAR\\.csv$")
MITOMAP_TRNA <- find_mitomap(DB_DIR, "MutationstRNA.*_VAR\\.csv$")

BLACKLIST_POS <- c(301, 302, 310, 316, 3107, 5894, 16179,
                   16181, 16182, 16183, 16188, 16189, 16192)

# Longer codes before shorter ones so prefix matching doesn't short-circuit
# (e.g. HV must be checked before H, L0-L5 before L)
HAP_ORDER <- c("HV","L0","L1","L2","L3","L4","L5",
               "A","B","C","D","E","F","G","H","I","J","K",
               "M","N","P","R","T","U","V","W","X","Y","Z")

# ---------------------------------------------------------------------------
# Load reference databases at startup
# ---------------------------------------------------------------------------

# gnomAD mtDNA: pre-processed by prepare_databases.R
# Columns: VAR, gnomAD_AF_hom, gnomAD_AF_het, gnomAD_AC_hom, gnomAD_AC_het,
#          gnomAD_AN, gnomAD_max_hl, gnomAD_filters,
#          hap_AF_hom_{A..Z}, hap_AF_het_{A..Z}
load_gnomad <- function(db_dir) {
  f <- file.path(db_dir, "gnomad", "gnomad_chrM_hap_AF.tsv.gz")
  if (!file.exists(f)) { message("gnomAD not found: ", f); return(NULL) }
  message("gnomAD: loading ", f)
  as.data.frame(fread(f))
}

# MITOMAP coding: Locus, Disease, Status, GB Freq FL (CR)*‡, GB Seqs FL (CR)*
# MITOMAP tRNA:   Position, Disease, Status, Homoplasmy, Heteroplasmy, MitoTIP†
# Merge both into one table keyed by VAR; coding and tRNA columns are the same
# conceptually so we harmonise to: MITOMAP_locus, MITOMAP_disease,
# MITOMAP_status, MITOMAP_gb_freq, MITOMAP_gb_seqs
load_mitomap <- function(db_dir) {
  mitomap_cds  <- find_mitomap(db_dir, "MutationsCodingControl.*_VAR\\.csv$")
  mitomap_trna <- find_mitomap(db_dir, "MutationstRNA.*_VAR\\.csv$")
  rows <- list()
  
  if (file.exists(mitomap_cds)) {
    message("MITOMAP CDS: loading ", mitomap_cds)
    cds <- as.data.frame(fread(mitomap_cds))
    cds <- cds[!is.na(cds$VAR), ]
    cds$MITOMAP_locus   <- cds$Locus
    cds$MITOMAP_disease <- cds$Disease
    cds$MITOMAP_status  <- cds$Status
    # "Plasmy Reports(Homo/Hetero)" already in Homo/Het format e.g. "+/-"
    plasmy_col <- grep("Plasmy", colnames(cds), value = TRUE)[1]
    cds$MITOMAP_plasmy  <- if (!is.na(plasmy_col)) cds[[plasmy_col]] else NA
    gb_freq_col <- grep("GB Freq", colnames(cds), value = TRUE)[1]
    gb_seq_col  <- grep("GB Seqs", colnames(cds), value = TRUE)[1]
    cds$MITOMAP_gb_freq <- if (!is.na(gb_freq_col)) cds[[gb_freq_col]] else NA
    cds$MITOMAP_gb_seqs <- if (!is.na(gb_seq_col))  cds[[gb_seq_col]]  else NA
    rows[["cds"]] <- cds[, c("VAR","MITOMAP_locus","MITOMAP_disease","MITOMAP_status",
                             "MITOMAP_plasmy","MITOMAP_gb_freq","MITOMAP_gb_seqs")]
  }
  
  if (file.exists(mitomap_trna)) {
    message("MITOMAP tRNA: loading ", mitomap_trna)
    trna <- as.data.frame(fread(mitomap_trna))
    trna <- trna[!is.na(trna$VAR), ]
    trna$MITOMAP_locus   <- trna$Locus
    trna$MITOMAP_disease <- trna$Disease
    trna$MITOMAP_status  <- trna$Status
    # Separate Homoplasmy / Heteroplasmy columns — unify to "Homo/Het" format
    hom_col <- grep("^Homoplasmy$", colnames(trna), value = TRUE)[1]
    het_col <- grep("^Heteroplasmy$", colnames(trna), value = TRUE)[1]
    trna$MITOMAP_plasmy <- if (!is.na(hom_col) && !is.na(het_col))
      paste0(trna[[hom_col]], "/", trna[[het_col]]) else NA
    gb_freq_col <- grep("GB Freq", colnames(trna), value = TRUE)[1]
    gb_seq_col  <- grep("GB Seqs", colnames(trna), value = TRUE)[1]
    trna$MITOMAP_gb_freq <- if (!is.na(gb_freq_col)) trna[[gb_freq_col]] else NA
    trna$MITOMAP_gb_seqs <- if (!is.na(gb_seq_col))  trna[[gb_seq_col]]  else NA
    rows[["trna"]] <- trna[, c("VAR","MITOMAP_locus","MITOMAP_disease","MITOMAP_status",
                               "MITOMAP_plasmy","MITOMAP_gb_freq","MITOMAP_gb_seqs")]
  }
  
  if (length(rows) == 0) return(NULL)
  # Multiple MITOMAP entries per VAR are possible; collapse with " | "
  combined <- bind_rows(rows)
  combined %>%
    group_by(VAR) %>%
    summarise(across(everything(), ~paste(unique(na.omit(.)), collapse = " | ")),
              .groups = "drop") %>%
    as.data.frame()
}

# MitoTIP: Position, rCRS, Alt, MitoTIP_Score, Quartile, Mitomap_Status, VAR
load_mitotip <- function(db_dir) {
  f <- file.path(db_dir, "MitoTIP", "mitotip_scores_27_04_2020_VAR.txt")
  if (!file.exists(f)) { message("MitoTIP not found: ", f); return(NULL) }
  message("MitoTIP: loading ", f)
  d <- as.data.frame(fread(f))
  d <- d[!is.na(d$VAR), c("VAR","MitoTIP_Score","Quartile")]
  colnames(d)[2:3] <- c("mitoTIP_score","mitoTIP_quartile")
  d
}

# PhyloP: produced by prepare_databases.R, columns: POS, PhyloP

# Databases are loaded reactively inside the server when the user clicks Load

# ---------------------------------------------------------------------------
# Parse annotated VCF (mutect2.{03|05|10}.vcf)
# Returns data.frame with one row per variant + a VAR join key
# ---------------------------------------------------------------------------
parse_vcf <- function(vcf_file, sample_id) {
  lines    <- readLines(vcf_file)
  hdr      <- which(startsWith(lines, "#CHROM"))[1]
  data_lines <- lines[(hdr + 1):length(lines)]
  if (length(data_lines) == 0) return(data.frame())
  
  vcf <- read.table(text = data_lines, sep = "\t", header = FALSE,
                    stringsAsFactors = FALSE, quote = "")
  colnames(vcf) <- c("CHROM","POS","rsID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")
  
  parse_info <- function(info_str) {
    fields <- strsplit(info_str, ";")[[1]]
    out <- list()
    for (f in fields) {
      if (grepl("=", f)) {
        kv <- strsplit(f, "=", fixed = TRUE)[[1]]
        out[[kv[1]]] <- paste(kv[-1], collapse = "=")
      } else {
        out[[f]] <- TRUE
      }
    }
    out
  }
  get_field <- function(lst, key) {
    v <- lst[[key]]
    if (is.null(v) || isFALSE(v)) NA_character_ else as.character(v)
  }
  get_flag <- function(lst, key) !is.null(lst[[key]]) && !isFALSE(lst[[key]])
  
  info_list <- lapply(vcf$INFO, parse_info)
  
  # GT:DP:AF are in FORMAT/SAMPLE columns, not INFO
  fmt_keys <- strsplit(vcf$FORMAT[1], ":", fixed = TRUE)[[1]]
  samp_vals <- strsplit(vcf$SAMPLE, ":", fixed = TRUE)
  get_fmt <- function(vals, key) {
    i <- match(key, fmt_keys)
    if (is.na(i)) return(NA_character_)
    sapply(vals, function(v) if (length(v) >= i) v[i] else NA_character_)
  }
  vcf$GT  <- get_fmt(samp_vals, "GT")
  vcf$DP  <- as.integer(get_fmt(samp_vals, "DP"))
  vcf$VAF <- as.numeric(get_fmt(samp_vals, "AF"))
  vcf$CDS     <- sapply(info_list, get_field, "CDS")
  vcf$COMPLEX <- sapply(info_list, get_field, "COMPLEX")
  vcf$RNR     <- sapply(info_list, get_field, "RNR")
  vcf$TRN     <- sapply(info_list, get_field, "TRN")
  vcf$AP      <- sapply(info_list, get_field, "AP")
  vcf$APS     <- as.numeric(sapply(info_list, get_field, "APS"))
  vcf$MLC_consq  <- sapply(info_list, get_field, "MLC_consq")
  vcf$MLC_score  <- as.numeric(sapply(info_list, get_field, "MLC_score"))
  vcf$MCC        <- as.numeric(sapply(info_list, get_field, "MCC"))
  vcf$Hypervariable <- sapply(info_list, get_flag, "Hypervariable")
  dloop_flag        <- sapply(info_list, get_flag, "DLOOP")
  vcf$Homopolymer <- sapply(info_list, get_flag, "Homopolymer")
  vcf$is_INDEL    <- sapply(info_list, get_flag, "INDEL")
  
  vcf$Region <- case_when(
    !is.na(vcf$CDS) ~ paste0("CDS:", vcf$CDS),
    !is.na(vcf$RNR) ~ paste0("rRNA:", vcf$RNR),
    !is.na(vcf$TRN) ~ paste0("tRNA:", vcf$TRN),
    dloop_flag      ~ "D-loop",
    TRUE            ~ "Other"
  )
  vcf$Gene     <- coalesce(vcf$CDS, vcf$RNR, vcf$TRN)
  vcf$VAR      <- paste0(vcf$POS, "_", vcf$REF, "_", vcf$ALT)
  vcf$SampleID <- sample_id
  
  vcf[, c("SampleID","VAR","POS","REF","ALT","rsID","FILTER","VAF","DP","GT",
          "Region","Gene","COMPLEX","MLC_consq","AP","APS",
          "MLC_score","MCC","Hypervariable",
          "Homopolymer","is_INDEL")]
}

# ---------------------------------------------------------------------------
# Load one sample directory
# ---------------------------------------------------------------------------
load_sample <- function(sample_folder, sample_id = NULL) {
  if (is.null(sample_id)) sample_id <- basename(sample_folder)
  hits <- list.files(file.path(sample_folder, "out"),
                     pattern = "\\.mutect2\\.mutect2\\.00\\.vcf$",
                     recursive = TRUE, full.names = TRUE)
  if (length(hits) == 0) stop(paste("No mutect2.00 VCF found under:", sample_folder))
  if (length(hits) > 1)  message(sprintf("Multiple VCFs found for %s, using first: %s", sample_id, hits[1]))
  vcf_file <- hits[1]
  
  vars <- parse_vcf(vcf_file, sample_id)
  
  hg_file <- file.path(sample_folder, "out", "mutect2.haplogroup.tab")
  haplogroup <- if (file.exists(hg_file)) fread(hg_file)$haplogroup[1] else NA_character_
  
  cvg_file <- file.path(sample_folder, "out", "count.tab")
  coverage  <- if (file.exists(cvg_file)) as.data.frame(fread(cvg_file)) else NULL
  
  hc_file   <- file.path(sample_folder, "out", "mutect2.haplocheck.tab")
  haplocheck <- if (file.exists(hc_file)) as.data.frame(fread(hc_file)) else NULL
  
  list(vars = vars, haplogroup = haplogroup,
       coverage = coverage, haplocheck = haplocheck,
       sample_id = sample_id)
}

# ---------------------------------------------------------------------------
# Load all samples from manifest TSV (SampleID, FamilyID, SampleFolder)
# ---------------------------------------------------------------------------
load_all_samples <- function(manifest_path) {
  mf <- as.data.frame(fread(manifest_path))
  if (!all(c("SampleID","SampleFolder") %in% colnames(mf)))
    stop("Manifest requires at least SampleID and SampleFolder columns.")
  samples <- list()
  for (i in seq_len(nrow(mf))) {
    sid <- mf$SampleID[i]
    s   <- tryCatch(load_sample(mf$SampleFolder[i], sid),
                    error = function(e) { message(e$message); NULL })
    if (!is.null(s)) samples[[sid]] <- s
  }
  samples
}

# ---------------------------------------------------------------------------
# Stack all samples into long format: one row per sample x variant
# ---------------------------------------------------------------------------
merge_variants <- function(samples_list) {
  if (length(samples_list) == 0) return(data.frame())
  
  all_vars <- lapply(names(samples_list), function(sid) {
    s   <- samples_list[[sid]]
    dat <- s$vars
    dat$haplogroup <- s$haplogroup
    dat$genotype   <- case_when(
      is.na(dat$VAF)   ~ NA_character_,
      dat$VAF >= 0.95  ~ "hom",
      dat$VAF >= 0.03  ~ "het",
      TRUE             ~ "low_het"
    )
    dat$blacklist <- dat$POS %in% BLACKLIST_POS
    dat
  })
  
  merged <- bind_rows(all_vars)
  
  # Cohort-level allele frequencies per VAR
  n_samples <- length(samples_list)
  cohort_ac <- merged |>
    group_by(VAR) |>
    summarise(
      cohort_AC_hom = sum(genotype == "hom", na.rm = TRUE),
      cohort_AC_het = sum(genotype == "het", na.rm = TRUE),
      .groups = "drop"
    )
  
  left_join(merged, cohort_ac, by = "VAR")
}

# ---------------------------------------------------------------------------
# Join pre-loaded reference databases
# ---------------------------------------------------------------------------
annotate_variants <- function(vars, db) {
  if (!is.null(db$gnomad)) {
    gnomad <- db$gnomad[, !(colnames(db$gnomad) %in% c("POS","REF","ALT")), drop = FALSE]
    vars <- left_join(vars, gnomad, by = "VAR")
  }
  if (!is.null(db$mitomap)) vars <- left_join(vars, db$mitomap, by = "VAR")
  if (!is.null(db$mitotip)) vars <- left_join(vars, db$mitotip, by = "VAR")
  
  if (!is.null(db$gnomad)) {
    # For each row look up its own haplogroup-specific AF
    match_hap <- function(hap) {
      if (is.na(hap)) return(NA_character_)
      for (h in HAP_ORDER) if (startsWith(toupper(hap), toupper(h))) return(h)
      NA_character_
    }
    hap_matched <- sapply(vars$haplogroup, match_hap)
    
    vars$gnomAD_hap_AF_hom <- mapply(function(hap, i) {
      col <- paste0("hap_AF_hom_", hap)
      if (!is.na(hap) && col %in% colnames(vars)) vars[[col]][i] else NA_real_
    }, hap_matched, seq_len(nrow(vars)))
    
    vars$gnomAD_hap_AF_het <- mapply(function(hap, i) {
      col <- paste0("hap_AF_het_", hap)
      if (!is.na(hap) && col %in% colnames(vars)) vars[[col]][i] else NA_real_
    }, hap_matched, seq_len(nrow(vars)))
    
    # Drop raw 58-column arrays
    vars <- vars[, !(colnames(vars) %in% c(paste0("hap_AF_hom_", HAP_ORDER),
                                           paste0("hap_AF_het_", HAP_ORDER)))]
  }
  
  # Clickable links
  vars$gnomAD_link <- sprintf(
    '<a href="https://gnomad.broadinstitute.org/variant/M-%d-%s-%s?dataset=gnomad_r3" target="_blank">gnomAD</a>',
    vars$POS, vars$REF, vars$ALT)
  
  vars
}

# ---------------------------------------------------------------------------
# Pedigree from PLINK FAM file
# ---------------------------------------------------------------------------
load_fam <- function(fam_path) {
  fam <- as.data.frame(fread(fam_path, header = FALSE))
  colnames(fam) <- c("FamilyID","IndividualID","FatherID","MotherID","Sex","Affected")
  fam$Sex    <- ifelse(fam$Sex == 1, "male", ifelse(fam$Sex == 2, "female", "unknown"))
  fam$Status <- case_when(fam$Affected == 2 ~ "affected",
                          fam$Affected == 1 ~ "unaffected",
                          TRUE              ~ "unknown")
  fam
}

# ---------------------------------------------------------------------------
# Family inheritance filter
# mode: "all" | "shared" | "de_novo"
#
# Comparators (used for shared / de novo):
#   - Mother if present in the dataset
#   - Siblings: same FamilyID + same MotherID as proband (if mother known),
#     otherwise all other loaded family members
# At least one comparator must be present in the loaded data.
# ---------------------------------------------------------------------------
# Family filter
# Determines the family member set (proband + mother + siblings on maternal line),
# restricts vars to those members, then applies mode:
#
#   all      — all variants for all family members
#   shared   — variants present in every family member (intersection)
#   de_novo  — variants not shared by all members (symmetric difference:
#               unique to proband OR unique to one or more other members)
#
# vars: variant table already filtered by other criteria (all samples present)
family_filter <- function(vars, fam, family_id, proband_id, mode) {
  fam_sub     <- fam[fam$FamilyID == family_id, ]
  proband_row <- fam_sub[fam_sub$IndividualID == proband_id, ]
  if (nrow(proband_row) == 0) stop("Proband not found in FAM file.")
  
  mother_id  <- proband_row$MotherID
  has_mother <- !is.null(mother_id) && !is.na(mother_id) &&
    mother_id != "0" && mother_id %in% vars$SampleID
  
  # Siblings: same maternal lineage as proband
  if (has_mother) {
    sib_ids <- fam_sub$IndividualID[
      !is.na(fam_sub$MotherID) & fam_sub$MotherID == mother_id &
        fam_sub$IndividualID != proband_id &
        fam_sub$IndividualID %in% vars$SampleID
    ]
  } else {
    sib_ids <- fam_sub$IndividualID[
      fam_sub$IndividualID != proband_id &
        fam_sub$IndividualID %in% vars$SampleID
    ]
  }
  
  member_ids <- unique(c(proband_id, if (has_mother) mother_id, sib_ids))
  if (length(member_ids) == 0) stop("No family members found in the loaded dataset.")
  message("Family members: ", paste(member_ids, collapse = ", "))
  
  # Restrict to family members
  family_vars <- vars[vars$SampleID %in% member_ids, ]
  if (nrow(family_vars) == 0) stop("No variants found for any family member.")
  if (mode == "all") return(family_vars)
  
  # Per-member VAR sets — basis for shared / de_novo
  vars_per_member <- tapply(family_vars$VAR, family_vars$SampleID,
                            function(v) unique(v[!is.na(v)]))
  
  if (mode == "shared") {
    # Present in proband AND at least one comparator (mother or sibling)
    comparator_ids  <- setdiff(member_ids, proband_id)
    proband_vars    <- vars_per_member[[proband_id]]
    comparator_vars <- unique(unlist(vars_per_member[comparator_ids]))
    shared_vars     <- intersect(proband_vars, comparator_vars)
    return(family_vars[family_vars$VAR %in% shared_vars, ])
  }
  
  if (mode == "de_novo") {
    # Present in proband AND absent from all comparators (mother or sibling)
    comparator_ids  <- setdiff(member_ids, proband_id)
    proband_vars    <- vars_per_member[[proband_id]]
    comparator_vars <- unique(unlist(vars_per_member[comparator_ids]))
    denovo_vars     <- setdiff(proband_vars, comparator_vars)
    return(family_vars[family_vars$VAR %in% denovo_vars, ])
  }
  
  family_vars
}

# ---------------------------------------------------------------------------

# ===========================================================================
# Shared DT options
# ===========================================================================
dt_options <- list(
  dom        = "Blfrtip",
  scrollX    = TRUE,
  autoWidth  = TRUE,
  buttons    = c("copy","csv","excel","pdf","print"),
  lengthMenu = list(c(25, 50, -1), c(25, 50, "All"))
)

# ===========================================================================
# UI
# ===========================================================================
ui <- navbarPage(
  title = "mtDNA Variant Curation",
  theme = shinytheme("cerulean"),
  
  # --- Home ---
  tabPanel(icon("home"),
           h1("mtDNA Variant Curation"),
           br(),
           h3("Overview"),
           p("Interactive viewer for mitochondrial DNA variant calls from short-read sequencing.
      Variants are called with Mutect2 using the GATK mitochondrial pipeline.
      Reference annotations are loaded from the local database at startup:
      gnomAD v3.1, MITOMAP, and mitoTIP (when available)."),
           
           h4("Input files"),
           p(strong("Manifest (TSV, required)"), "— columns:",
             code("SampleID"), ",", code("FamilyID"), "(optional),",
             code("SampleFolder"), "(absolute path). The app reads",
             code("*.mutect2.mutect2.00.vcf"), "found recursively under", code("out/"), "in each folder."),
           p(strong("FAM file (optional)"), "— PLINK format (6 columns: FamilyID, IndividualID,
      FatherID, MotherID, Sex, Affected). Required for maternal inheritance analysis."),
           
           h4("Variant filters"),
           tags$ul(
             tags$li(strong("AF threshold"), "— Mutect2 calling threshold: 03/05/10 (%)"),
             tags$li(strong("Min heteroplasmy"), "— further AF cutoff across all samples"),
             tags$li(strong("PASS only"), "— require FILTER = PASS in all samples"),
             tags$li(strong("Exclude D-loop / Hypervariable"), "— remove control region variants"),
             tags$li(strong("Region"), "— CDS, rRNA, tRNA, D-loop, or All"),
             tags$li(strong("Gene"), "— comma-separated (e.g.", code("ND1,COX1"), ") or All")
           ),
           
           h4("Family analysis"),
           p("mtDNA is strictly maternally inherited — no paternal contribution, no compound het.
      Comparators are mother (if loaded) and/or siblings (same family, same maternal line).
      Works with mother-only, sibling-only, or mixed families. Three modes:"),
           tags$ul(
             tags$li(strong("All"), "— all variants in proband (no family filter)"),
             tags$li(strong("Shared"), "— present in proband AND at least one comparator (mother or sibling, AF > 0)"),
             tags$li(strong("De novo"), "— present in proband, absent in all comparators (AF = 0 or not called)")
           ),
           br(), br()
  ),
  
  # --- Variants ---
  tabPanel("Variants",
           fluidPage(
             headerPanel("mtDNA Variant List"),
             sidebarPanel(width = 3,
                          h4("Input"),
                          textInput("db_dir", "Database directory", value = DB_DIR),
                          textInput("manifest_file", "Manifest file (path)",
                                    placeholder = "/path/to/manifest.tsv"),
                          textInput("fam_file", "FAM file (path, optional)",
                                    placeholder = "/path/to/cohort.fam"),
                          hr(),
                          h4("Variant filters"),
                          numericInput("min_af", "Min heteroplasmy level / VAF", value = 0.03, min = 0, max = 1, step = 0.01),
                          numericInput("max_cohort_ac", "Max cohort AC (hom + het)", value = 9999, min = 0, step = 1),
                          numericInput("max_gnomad_af", "Max gnomAD AF (leave 1 to skip)", value = 1, min = 0, max = 1, step = 0.001),
                          checkboxInput("pass_only",        "PASS only",             value = FALSE),
                          checkboxInput("excl_strand_bias", "Exclude strand bias",  value = TRUE),
                          checkboxInput("exclude_hv",       "Exclude hypervariable", value = FALSE),
                          checkboxInput("exclude_blacklist", "Exclude blacklist",        value = TRUE),
                          checkboxInput("exclude_synonymous", "Exclude synonymous",     value = FALSE),
                          checkboxInput("exclude_intergenic", "Exclude intergenic",     value = FALSE),
                          selectInput("region_filter", "Region",
                                      choices = c("All","CDS","rRNA","tRNA","D-loop"), selected = "All"),
                          textInput("gene_filter",   "Gene (e.g. ND1,COX1 or All)", value = "All"),
                          textInput("sample_filter", "Sample ID (leave blank for all)", value = ""),
                          hr(),
                          h4("Family analysis"),
                          textInput("family_id",  "Family ID",  value = ""),
                          textInput("proband_id", "Proband ID", value = ""),
                          selectInput("family_mode", "Family filter",
                                      choices = c("All"      = "all",
                                                  "Shared"   = "shared",
                                                  "De novo"  = "de_novo"),
                                      selected = "all"),
                          hr(),
                          actionButton("view_btn", "Load / Refresh", class = "btn-primary")
             ),
             mainPanel(width = 9,
                       DT::dataTableOutput("var_table")
             )
           )
  ),
  
  
)

# ===========================================================================
# Server
# ===========================================================================
server <- function(input, output, session) {
  
  db_data <- eventReactive(input$view_btn, {
    dir <- trimws(input$db_dir)
    db  <- list(
      gnomad  = load_gnomad(dir),
      mitomap = load_mitomap(dir),
      mitotip = load_mitotip(dir)
    )
    if (is.null(db$gnomad))  showNotification("gnomAD not found — gnomAD columns will be empty.",  type = "warning", duration = 10)
    if (is.null(db$mitomap)) showNotification("MITOMAP not found — MITOMAP columns will be empty.", type = "warning", duration = 10)
    if (is.null(db$mitotip)) showNotification("MitoTIP not found — MitoTIP columns will be empty.", type = "warning", duration = 10)
    db
  })
  
  samples_data <- eventReactive(input$view_btn, {
    req(nchar(trimws(input$manifest_file)) > 0)
    samps <- withProgress(message = "Loading samples...", {
      load_all_samples(trimws(input$manifest_file))
    })
    updateNumericInput(session, "max_cohort_ac", value = length(samps))
    samps
  })
  
  fam_data <- eventReactive(input$view_btn, {
    path <- trimws(input$fam_file)
    if (nchar(path) == 0 || !file.exists(path)) return(NULL)
    load_fam(path)
  })
  
  variants_merged <- eventReactive(input$view_btn, {
    samps <- samples_data()
    req(length(samps) > 0)
    vars <- merge_variants(samps)
    annotate_variants(vars, db_data())
  })
  
  variants_filtered <- reactive({
    vars <- variants_merged()
    fam  <- fam_data()
    
    # Min heteroplasmy filter
    if (input$min_af > 0)
      vars <- vars[!is.na(vars$VAF) & vars$VAF >= input$min_af, ]
    
    if (input$pass_only)
      vars <- vars[!is.na(vars$FILTER) & vars$FILTER == "PASS", ]
    
    if ("cohort_AC_hom" %in% colnames(vars) && "cohort_AC_het" %in% colnames(vars))
      vars <- vars[vars$cohort_AC_hom <= input$max_cohort_ac &
                     vars$cohort_AC_het <= input$max_cohort_ac, ]
    
    # gnomAD MAF filter — for hom: both gnomAD_AF_hom and gnomAD_hap_AF_hom must pass
    #                     for het: both gnomAD_AF_het and gnomAD_hap_AF_het must pass
    if (input$max_gnomad_af < 1) {
      thresh <- input$max_gnomad_af
      cn     <- colnames(vars)
      is_hom <- !is.na(vars$genotype) & vars$genotype == "hom"
      af_pass <- function(col_hom, col_het) {
        if (!(col_hom %in% cn) || !(col_het %in% cn)) return(rep(TRUE, nrow(vars)))
        af <- rep(NA_real_, nrow(vars))
        af[is_hom]  <- vars[[col_hom]][is_hom]
        af[!is_hom] <- vars[[col_het]][!is_hom]
        is.na(af) | af <= thresh
      }
      vars <- vars[af_pass("gnomAD_AF_hom",     "gnomAD_AF_het") &
                     af_pass("gnomAD_hap_AF_hom",  "gnomAD_hap_AF_het"), ]
    }
    
    if (input$excl_strand_bias)
      vars <- vars[!grepl("strand_bias", vars$FILTER, ignore.case = TRUE), ]
    
    if (input$exclude_hv)        vars <- vars[!vars$Hypervariable %in% TRUE, ]
    if (input$exclude_blacklist)  vars <- vars[!vars$blacklist %in% TRUE, ]
    if (input$exclude_synonymous) vars <- vars[is.na(vars$MLC_consq) | vars$MLC_consq != "synonymous_variant", ]
    if (input$exclude_intergenic) vars <- vars[is.na(vars$MLC_consq) | vars$MLC_consq != "intergenic_variant", ]
    
    if (input$region_filter != "All")
      vars <- vars[startsWith(vars$Region, input$region_filter), ]
    
    gene_input <- trimws(input$gene_filter)
    if (gene_input != "All" && nchar(gene_input) > 0) {
      genes <- trimws(strsplit(gene_input, ",")[[1]])
      vars  <- vars[!is.na(vars$Gene) & vars$Gene %in% genes, ]
    }
    
    sample_input <- trimws(input$sample_filter)
    if (nchar(sample_input) > 0)
      vars <- vars[vars$SampleID == sample_input, ]
    
    # Family filter
    if (!is.null(fam) &&
        nchar(trimws(input$family_id))  > 0 &&
        nchar(trimws(input$proband_id)) > 0) {
      vars <- tryCatch(
        family_filter(vars,
                      fam        = fam,
                      family_id  = trimws(input$family_id),
                      proband_id = trimws(input$proband_id),
                      mode       = input$family_mode),
        error = function(e) { showNotification(e$message, type = "error"); vars }
      )
    }
    vars
  })
  
  output$var_table <- DT::renderDataTable({
    vars <- variants_filtered()
    req(nrow(vars) > 0)
    
    # Column display order:
    #   SampleID | core anno | AF DP FILTER genotype haplogroup hap_AF
    #   | gnomAD | MITOMAP | mitoTIP | PhyloP | links
    core_cols <- c("SampleID","VAR","POS","REF","ALT","rsID","FILTER","VAF","DP",
                   "genotype","haplogroup","Region","Gene","COMPLEX",
                   "MLC_consq","AP","APS","MLC_score","MCC",
                   "Hypervariable","Homopolymer","is_INDEL","blacklist")
    db_cols   <- c("cohort_AC_hom","cohort_AC_het",
                   "gnomAD_AF_hom","gnomAD_AF_het","gnomAD_AC_hom","gnomAD_AC_het",
                   "gnomAD_max_hl",
                   "gnomAD_hap_AF_hom","gnomAD_hap_AF_het",
                   "MITOMAP_locus","MITOMAP_disease","MITOMAP_status","MITOMAP_plasmy",
                   "mitoTIP_score","mitoTIP_quartile")
    link_cols <- c("gnomAD_link")
    ord <- c(core_cols, db_cols, link_cols)
    vars <- vars[, ord[ord %in% colnames(vars)]]
    
    DT::datatable(vars,
                  caption    = "mtDNA variant list",
                  filter     = "top",
                  rownames   = FALSE,
                  escape     = FALSE,
                  extensions = "Buttons",
                  options    = dt_options)
  })
  
  
}

shinyApp(ui = ui, server = server)
