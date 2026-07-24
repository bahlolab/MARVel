#### target variants (generate candidate positions that avoid existing calls + homopolymers + breakpoints)
# rm(list = ls())
# ID <- "HG00407"
# seed <- 123

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop('Usage: Rscript 1_target_variant.R "HG01881" [seed]')
}
ID <- args[1]

seed <- if (length(args) >= 2) as.integer(args[2]) else 123
if (is.na(seed)) stop("Seed must be an integer.")
set.seed(seed)

suppressPackageStartupMessages({
  library(vcfR)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(readxl)
  library(Biostrings)
})

resDIR <- "/vast/projects/bahlo_mtDNA/Results_1000G/"
workDIR <- file.path("/vast/projects/bahlo_mtDNA/Simulation/1000G", ID)
if (!dir.exists(workDIR)) dir.create(workDIR, recursive = TRUE)
ref_fa <- readDNAStringSet("/vast/projects/bahlo_mtDNA/Simulation/1000G/rCRS.fasta")
ref_seq <- ref_fa[[1]]

# -------------------------------
# Helper: expand multi-allelic ALT and compute VCF-consistent end per allele
# end = POS + max(nchar(REF), nchar(ALT)) - 1
# (Insertion with REF length 1 => end == POS, deletion spans REF length, etc.)
# -------------------------------
vcf_expand_end <- function(df, pos_col = "POS", ref_col = "REF", alt_col = "ALT") {
  df[[pos_col]] <- as.integer(df[[pos_col]])
  
  df %>%
    tidyr::separate_rows(.data[[alt_col]], sep = ",") %>%
    dplyr::mutate(
      start = .data[[pos_col]],
      .ref_len = nchar(.data[[ref_col]]),
      .alt_len = nchar(.data[[alt_col]]),
      .span = pmax(.ref_len, .alt_len),
      end = start + .span - 1L
    ) %>%
    dplyr::select(-.ref_len, -.alt_len, -.span)
}

# -------------------------------
# Input files
# -------------------------------
mitoHPC_file  <- file.path(resDIR, ID, "mitoHPC", paste0(ID, ".final.mutect2.mutect2.00.vcf")) # edit vcf to vcf.gz
mity_file     <- file.path(resDIR, ID, "mity", paste0(ID, ".mity.report.xlsx"))
server2_file  <- file.path(resDIR, ID, "mtdna-server-2", "results", "variants.annotated.txt")
mtSwirl_file  <- file.path(resDIR, ID, "mtSwirl", paste0(ID, ".self.ref.split.selfToRef.final.vcf"))

# -------------------------------
# 1) Existing variants from tools (convert to start/end spans)
# -------------------------------

## mitoHPC (VCF)
mitoHPC_vcf <- suppressWarnings(read.vcfR(mitoHPC_file, verbose = FALSE))
mitoHPC_fix <- data.frame(mitoHPC_vcf@fix, stringsAsFactors = FALSE)
mitoHPC <- vcf_expand_end(mitoHPC_fix, pos_col = "POS", ref_col = "REF", alt_col = "ALT") %>%
  dplyr::select(start, end)

## mtSwirl (VCF)
mtSwirl_vcf <- suppressWarnings(read.vcfR(mtSwirl_file, verbose = FALSE))
mtSwirl_fix <- data.frame(mtSwirl_vcf@fix, stringsAsFactors = FALSE)
mtSwirl <- vcf_expand_end(mtSwirl_fix, pos_col = "POS", ref_col = "REF", alt_col = "ALT") %>%
  dplyr::select(start, end)

## mity (XLSX)
mity_raw <- readxl::read_excel(mity_file)
# Expect at least POS / REF / ALT columns; keep robust handling for multi-allelic ALT too
mity <- mity_raw %>%
  dplyr::mutate(POS = as.integer(POS)) %>%
  tidyr::separate_rows(ALT, sep = ",") %>%
  dplyr::mutate(
    start = POS,
    end = start + pmax(nchar(REF), nchar(ALT)) - 1L
  ) %>%
  dplyr::select(start, end)

## mtDNA-Server 2 (TSV)
server2 <- read.table(server2_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
server2 <- server2 %>%
  dplyr::mutate(
    start = as.integer(Pos),
    end = start + nchar(Variant) - 1L
  ) %>%
  dplyr::select(start, end)

# Combine and deduplicate existing variant spans
vars <- dplyr::bind_rows(server2, mity, mitoHPC, mtSwirl) %>%
  dplyr::mutate(id = paste0(start, ":", end)) %>%
  dplyr::distinct(id, .keep_all = TRUE) %>%
  dplyr::select(start, end)

# -------------------------------
# 2) Homopolymer tracts (+/- 2bp padding)
# -------------------------------
hpl <- read.table("/vast/projects/bahlo_mtDNA/1000G/reference/homopolymer_tracts.txt", header = TRUE)
hpl <- hpl %>%
  dplyr::mutate(
    start = start - 2L,
    end   = end + 2L
  )

# -------------------------------
# 3) Artificial breakpoints to avoid
# -------------------------------
brk <- c(1:150, 16569 + (-149:0))

# -------------------------------
# All possible sites (SNV-style single-base positions)
# -------------------------------
VAR.all <- data.frame(chr = "chrM", start = 1:16569, end = 1:16569, VAF = 0.1)

inside_hpl <- sapply(
  VAR.all$start,
  function(p) any(p >= hpl$start & p <= hpl$end)
)

inside_vars <- sapply(
  VAR.all$start,
  function(p) any(p >= vars$start & p <= vars$end)
)

# Final candidate set
VAR <- VAR.all[!inside_hpl &
                 !inside_vars &
                 VAR.all$start != 3107 &
                 !(VAR.all$start %in% brk) &
                 !(VAR.all$start %in% c(301:319, 16179:16194)), ]

# Write out (no header, tab-separated, no quotes)
out_file <- file.path(workDIR, "target.vars.txt")
write.table(VAR, out_file, row.names = FALSE, col.names = T, quote = FALSE, sep = "\t")

message("Wrote: ", out_file, " (n = ", nrow(VAR), " candidate sites)")

# -------------------------------
# Pick SNVs + Indels
# -------------------------------

genome_len <- 16569
min_dist   <- 200

# =========================
# 1) Pick SNVs (n_var = 40)
# =========================
n_var <- 40

VAR_binned <- VAR %>%
  mutate(
    bin = cut(
      start,
      breaks = seq(1, genome_len + 1, length.out = n_var + 1),
      include.lowest = TRUE,
      labels = FALSE
    )
  )

# one candidate per bin
candidates <- VAR_binned %>%
  group_by(bin) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  arrange(start)

# enforce spacing greedily
keep <- logical(nrow(candidates))
last_pos <- -Inf
for (i in seq_len(nrow(candidates))) {
  if (candidates$start[i] - last_pos >= min_dist) {
    keep[i] <- TRUE
    last_pos <- candidates$start[i]
  }
}
VAR_SNV <- candidates[keep, ]

# back-fill if < n_var (random order; still enforce min_dist)
if (nrow(VAR_SNV) < n_var) {
  extras <- VAR_binned %>%
    filter(!start %in% VAR_SNV$start) %>%
    distinct(start, .keep_all = TRUE) %>%  # avoid duplicates
    slice_sample(prop = 1)                 # shuffle
  
  for (p in extras$start) {
    if (all(abs(p - VAR_SNV$start) >= min_dist)) {
      VAR_SNV <- bind_rows(VAR_SNV, extras[extras$start == p, ])
      if (nrow(VAR_SNV) == n_var) break
    }
  }
}

VAR_SNV <- VAR_SNV %>% arrange(start)

if (nrow(VAR_SNV) < n_var) {
  stop(sprintf("Could only pick %d SNVs with min_dist=%d. Reduce n_var or min_dist.", nrow(VAR_SNV), min_dist))
}

write.table(
  VAR_SNV[, 1:4],
  file = file.path(workDIR, "target.40snv.txt"),
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

# =========================
# 2) Pick indels (5 DEL + 5 INS), all >= min_dist from SNVs and each other
# =========================

# positions allowed for variant placement (already filtered for hpl/artifacts etc.)
allowed_pos <- sort(unique(VAR$start))
allowed_set <- allowed_pos

# helper: enforce >= min_dist away from an existing set of positions
is_far_enough <- function(pos, others, min_dist) {
  if (length(others) == 0) return(TRUE)
  all(abs(pos - others) >= min_dist)
}

# for deletions: require that ALL bases in [start, end] are allowed
del_candidates <- function(L) {
  cand <- allowed_pos[allowed_pos + (L - 1) <= genome_len]
  cand[sapply(cand, function(p) all((p:(p + L - 1)) %in% allowed_set))]
}

# sample one position from candidates subject to distance constraints
pick_one <- function(candidates, forbidden_positions, min_dist = 200) {
  for (p in sample(candidates)) {
    if (is_far_enough(p, forbidden_positions, min_dist)) return(p)
  }
  stop("Failed to find a valid position under constraints. Try relaxing min_dist or increasing candidate space.")
}

# Start forbidden set with SNV positions
forbidden <- sort(unique(VAR_SNV$start))

# ---- Deletions: 1,3,5,7,9 bp ----
del_lengths <- c(1, 3, 5, 7, 9)
DEL <- vector("list", length(del_lengths))

for (i in seq_along(del_lengths)) {
  L <- del_lengths[i]
  cand <- del_candidates(L)
  
  # also avoid starts that are too close to already chosen indels
  p <- pick_one(cand, forbidden_positions = forbidden, min_dist = min_dist)
  # DEL[[i]] <- data.frame(
  #   chr = "chrM",
  #   start = p,
  #   end = p + L,
  #   VAF = 0.1,
  #   TYPE = "DEL",
  #   SEQ = "",
  #   stringsAsFactors = FALSE
  # )
  del_seq <- as.character(subseq(ref_seq, start = p+1, end = p + L))
  
  DEL[[i]] <- data.frame(
    chr = "chrM",
    start = p,
    end = p + L,
    VAF = 0.1,
    TYPE = "DEL",
    SEQ = del_seq,
    stringsAsFactors = FALSE
  )
  # update: add reference sequence for deletions
  forbidden <- c(forbidden, p)
}
DEL <- bind_rows(DEL)

# ---- Insertions: 1,3,5,7,9 bp sequences ----
ins_seqs <- c("C", "ATG", "ATGGC", "ATGGCAT", "ATGGCATGC")
INS <- vector("list", length(ins_seqs))

for (i in seq_along(ins_seqs)) {
  seqi <- ins_seqs[i]
  
  p <- pick_one(allowed_pos, forbidden_positions = forbidden, min_dist = min_dist)
  INS[[i]] <- data.frame(
    chr = "chrM",
    start = p,
    end = p,      # INS: start=end
    VAF = 0.1,
    TYPE = "INS",
    SEQ = seqi,
    stringsAsFactors = FALSE
  )
  forbidden <- c(forbidden, p)
}
INS <- bind_rows(INS)

# combine indels
VAR_indel <- bind_rows(DEL, INS) %>% arrange(start)

# sanity checks
stopifnot(nrow(VAR_indel) == 10)
stopifnot(min(diff(sort(VAR_indel$start))) >= min_dist)
# also ensure far from SNVs
stopifnot(all(sapply(VAR_indel$start, function(p) is_far_enough(p, VAR_SNV$start, min_dist))))
stopifnot(all(sapply(VAR_indel$end, function(p) is_far_enough(p, VAR_SNV$start, min_dist))))

# write BamSurgeon addindel varfile
# (chr start end vaf TYPE [SEQ if INS])
out_indel <- VAR_indel %>%
  mutate(VAF = sprintf("%.2f", VAF)) %>%
  select(chr, start, end, VAF, TYPE, SEQ)

write.table(
  out_indel,
  file = file.path(workDIR, "target.10indel.txt"),
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE,
  sep = "\t"
)


