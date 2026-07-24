##### Haplogroup-specific homoplasmy analysis (Haplogroup H and L) #####

rm(list = ls())
setwd("/vast/projects/bahlo_mtDNA/Results_1000G/")
library(readxl)
library(vcfR)
library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(ggstatsplot)
library(tidyr)
library(ggplot2)

artifact_prone_sites <- c(301,302,310,316,3107,5894,10933,16179,16181,16182,16183,16188,16189,16192)
allids = read.csv("/vast/projects/bahlo_mtDNA/Scripts_Kelly/sample.tsv",sep = "\t")
all_ids = allids$library_name
pedigree = read.csv("/vast/projects/bahlo_mtDNA/Scripts_Kelly/pedigree.tsv",sep = "\t")
ped <- pedigree %>% 
  filter(MotherID != 0) %>% 
  filter(entity.pedigree_id %in% all_ids) %>% 
  filter(MotherID %in% all_ids) %>% 
  select(entity.pedigree_id, MotherID) %>% 
  rename(mother = MotherID,
         proband = entity.pedigree_id)


base_dir <- "/vast/projects/bahlo_mtDNA/Results_1000G"
id_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
id_dirs <- id_dirs[grepl("^(HG|NA)", basename(id_dirs))]
ids = sort(c(ped$mother,ped$proband))
id_dirs <- id_dirs[basename(id_dirs) %in% ids] 

haplogroups <- read.table("/vast/projects/bahlo_mtDNA/Results_1000G/haplogroups_all.tsv",header = T)

## only run 1 depending on the haplogroup to be analysed
ids_with_L <- haplogroups %>%
  filter(grepl("^L", Haplogroup)) %>% 
  pull(SAMPLE)
ped_L <- ped %>%
  filter(mother %in% ids_with_L | proband %in% ids_with_L)
id_dirs_L <- id_dirs[basename(id_dirs) %in% ids_with_L]
id_dirs <- id_dirs_L
ped <- ped_L

ids_with_H <- haplogroups %>%
  filter(Haplogroup == "H") %>%
  pull(SAMPLE)  
ped_H <- ped %>%
  filter(mother %in% ids_with_H | proband %in% ids_with_H)
id_dirs_H <- id_dirs[basename(id_dirs) %in% ids_with_H]
id_dirs <- id_dirs_H
ped <- ped_H
## ---------- mtDNA-Server2 ----------
process_annotated <- function(id_dir) {
  file_path <- file.path(id_dir, "/mtdna-server-2/results/variants.annotated.txt")
  
  if (!file.exists(file_path)) {
    message("Missing file for ID: ", file_path)
    return(NULL)}
  
  df <- readr::read_tsv(file_path, col_types = readr::cols(.default = "c"))
  
  df %>%
    mutate(
      Pos = as.numeric(Pos),
      VAF = suppressWarnings(as.numeric(VariantLevel)),
      SAMPLE = sub("\\.chrM.*", "", ID),
      VAR = paste(Pos, Ref, Variant, sep = "_")
    ) %>%
    filter(!Pos %in% artifact_prone_sites) %>%
    filter(!grepl("strand_bias", Filter))
}

mtdna_server <- map_dfr(id_dirs, process_annotated)

mtdna_server_hom <- mtdna_server %>%
  filter(VAF >= 0.95) %>%
  select(SAMPLE, VAR,VAF)

## ---------- mity ----------
process_mity <- function(id_dir) {
  file_path <- file.path(id_dir, paste0("/mity/",basename(id_dir), ".mity.report.xlsx"))
  if (!file.exists(file_path)) {
    message("Missing file for ID: ", file_path)
    return(NULL)}
  df <- read_excel(file_path)
  
  df %>%
    filter(FILTER == "PASS" | (FILTER == "FAIL" & POS_FILTER == "0" & SBR_FILTER == "1" & SBA_FILTER == "1" & MQMR_FILTER == "1" & AQR_FILTER == "1")) %>%
    # filter(FILTER == "PASS") %>% 
    mutate(
      POS = as.numeric(POS),
      VAF = suppressWarnings(as.numeric(`VARIANT HETEROPLASMY`)),
      VAR = paste(POS, REF, ALT, sep = "_")
    ) %>%
    filter(!POS %in% artifact_prone_sites)
}

mity <- map_dfr(id_dirs, process_mity)
mity_concise <- mity %>% select(SAMPLE, HGVS, `GENE/LOCUS`, `VARIANT HETEROPLASMY`, CHR, POS, REF, ALT, FILTER, VAF, VAR,GT_FORMAT,FILTER,INFO,FORMAT)

mity_hom <- mity %>%
  filter(VAF >= 0.95) %>%  
  select(SAMPLE, VAR,VAF)

## ---------- mitoHPC (Mutect2) ----------
process_vcf <- function(id_dir) {
  file_path <- file.path(id_dir, "/mitoHPC/mutect2.mutect2.03.vcf")
  if (!file.exists(file_path)) {
    message("Missing file for ID: ", file_path)
    return(NULL)}
  
  vcf <- read.vcfR(file_path)
  fix_df <- as.data.frame(vcf@fix, stringsAsFactors = FALSE)
  gt_df  <- as.data.frame(vcf@gt,  stringsAsFactors = FALSE)
  
  fix_df$SAMPLE <- sub("\\.final", "", gt_df[,2])
  
  af_from_info <- function(info) {
    m <- str_match(info, "(?:^|;)AF=([^;]+)")
    suppressWarnings(as.numeric(m[,2]))
  }
  
  fix_df %>%
    transmute(
      CHROM = CHROM,
      POS   = as.numeric(POS),
      ID    = ID,
      REF   = REF,
      ALT   = ALT,
      QUAL  = QUAL,
      FILTER= FILTER,
      INFO  = INFO,
      SAMPLE= SAMPLE,
      VAF   = af_from_info(INFO),
      VAR   = paste(POS, REF, ALT, sep = "_")
    ) %>%
    filter(!POS %in% artifact_prone_sites) %>%
    filter(!grepl("strand_bias", FILTER)) %>%
    filter(!is.na(VAF))
}

mitoHPC <- map_dfr(id_dirs, process_vcf)

mitoHPC_hom <- mitoHPC %>%
  filter(VAF >= 0.95) %>%
  select(SAMPLE, VAR)

## ---------- mtSwirl ----------
read_mtSwirl_one <- function(path) {
  v <- read.vcfR(path)
  fix <- as.data.frame(v@fix, stringsAsFactors = FALSE)
  gt  <- as.data.frame(v@gt,  stringsAsFactors = FALSE)
  fmt <- gt$FORMAT[1]
  fields <- strsplit(fmt, ":")[[1]]
  af_idx <- match(TRUE, fields %in% c("AF","VAF","VF"))
  if (is.na(af_idx)) af_idx <- 2
  samp_col <- colnames(gt)[2]
  vafs <- suppressWarnings(as.numeric(sapply(gt[[samp_col]], function(x) strsplit(x, ":")[[1]][af_idx])))
  tibble(
    POS = as.numeric(fix$POS),
    REF = fix$REF,
    ALT = fix$ALT,
    FILTER = fix$FILTER,
    VAF = vafs,
    VAR = paste(fix$POS, fix$REF, fix$ALT, sep = "_")
  ) %>%
    filter(!POS %in% artifact_prone_sites) %>%
    filter(!grepl("strand_bias|blacklisted_site|base_qual|position", FILTER))
}

mtSwirl_hom_list <- lapply(seq_len(nrow(ped)), function(i) {
  mid <- ped$mother[i]; iid <- ped$proband[i]
  mid_df <- read_mtSwirl_one(paste0("/vast/projects/bahlo_mtDNA/Results_1000G/", mid, "/mtSwirl/",mid,".self.ref.split.selfToRef.final.vcf")) %>% filter(VAF >= 0.95) %>% mutate(SAMPLE = mid)
  iid_df <- read_mtSwirl_one(paste0("/vast/projects/bahlo_mtDNA/Results_1000G/", iid, "/mtSwirl/",iid,".self.ref.split.selfToRef.final.vcf")) %>% filter(VAF >= 0.95) %>% mutate(SAMPLE = iid)
  bind_rows(mid_df, iid_df) %>% select(SAMPLE, VAR,VAF)
})
mtSwirl_hom <- bind_rows(mtSwirl_hom_list)

## ---------- overlap counting ----------
tools <- c("mity","mtDNA_Server2","mitoHPC","mtSwirl")
overlap <- tibble(
  FID  = rep(1:54, each = length(tools)), # change number depending on mother-offspring pairs
  Tool = rep(tools, times = 54),
  HomN = NA_integer_,
  HomVariants = NA_character_
)

for (i in seq_len(nrow(ped))) {
  mid <- ped$mother[i]; iid <- ped$proband[i]
  
  # build per-tool sets
  s <- list(
    mtDNA_Server2 = intersect(
      mtdna_server_hom %>% filter(SAMPLE == mid) %>% pull(VAR),
      mtdna_server_hom %>% filter(SAMPLE == iid) %>% pull(VAR)
    ),
    mity = intersect(
      mity_hom %>% filter(SAMPLE == mid) %>% pull(VAR),
      mity_hom %>% filter(SAMPLE == iid) %>% pull(VAR)
    ),
    mitoHPC = intersect(
      mitoHPC_hom %>% filter(SAMPLE == mid) %>% pull(VAR),
      mitoHPC_hom %>% filter(SAMPLE == iid) %>% pull(VAR)
    ),
    mtSwirl = intersect(
      mtSwirl_hom %>% filter(SAMPLE == mid) %>% pull(VAR),
      mtSwirl_hom %>% filter(SAMPLE == iid) %>% pull(VAR)
    )
  )
  
  for (t in names(s)) {
    idx <- which(overlap$FID == i & overlap$Tool == t)
    vs  <- s[[t]]
    overlap$HomN[idx] <- length(vs)
    overlap$HomVariants[idx] <- if (length(vs) == 0) NA_character_ else paste(vs, collapse = ",")
  }
}

setwd("/vast/projects/bahlo_mtDNA/Analysis/homoplasmy/allduos")
count <- overlap
count <- count %>%
  mutate(
    FID  = factor(FID),
    Tool = factor(Tool, levels = c("mity","mtDNA_Server2","mitoHPC","mtSwirl")),
    HomN = as.numeric(HomN)
  )
long <- count %>%
  mutate(HomVariants = strsplit(HomVariants, ",")) %>%
  unnest(HomVariants)

variant_counts <- long %>%
  group_by(HomVariants, Tool) %>%
  summarise(count = n(), .groups = "drop")

variant_matrix <- variant_counts %>%
  pivot_wider(names_from = Tool, values_from = count, values_fill = 0)

variant_matrix <- variant_matrix %>%
  mutate(Pos = as.numeric(str_extract(HomVariants, "^[0-9]+"))) %>%
  arrange(Pos) %>%
  select(-Pos)

variant_matrix_df <- variant_matrix %>% as.data.frame()  
rownames(variant_matrix_df) <- variant_matrix_df$HomVariants
variant_matrix_df <- variant_matrix_df[, -which(names(variant_matrix_df) == "HomVariants")]

variant_diff <- variant_matrix_df[rowSums(variant_matrix_df != variant_matrix_df[, 1]) > 0, ]

write.table(variant_diff, "hom_diff_trio_hapH.tsv", quote = F, sep = "\t") # change file name depending on haplogroup
