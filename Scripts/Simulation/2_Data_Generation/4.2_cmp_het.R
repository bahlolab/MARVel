# rm(list = ls())
# ID <- "HG01881"
# VAF <- 50

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop('Usage: Rscript 4.2_cmp_het.R "HG01881" 10')
}
ID <- args[1]
VAF <- as.integer(args[2])

suppressPackageStartupMessages({
  library(dplyr)
  library(vcfR)
  library(stringr)
  library(tidyr)
})

work_dir <- paste0("/vast/projects/bahlo_mtDNA/Simulation/1000G/",ID,"/AF",VAF,"_DS2k/")

## --- Load inputs -------------------------------------------------------------

snv <- read.table(
  paste0(work_dir,"target.snv.",VAF,".txt"), 
  header = FALSE, stringsAsFactors = FALSE)

indel <- read.table(
  paste0(work_dir,"target.indel.",VAF,".txt"),
  sep = "\t", header = FALSE, fill = TRUE, stringsAsFactors = FALSE
)
indel$base <- indel$V3 - indel$V2
indel$base[indel$V5 == "INS"] <- nchar(indel$V6[indel$V5 == "INS"])

incov  <- read.table(paste0(work_dir,"chrM.incover.txt"),  header = FALSE, stringsAsFactors = FALSE)
outcov <- read.table(paste0(work_dir,"chrM.outcover.txt"), header = FALSE, stringsAsFactors = FALSE)

## --- Coverage difference (incov - outcov) -----------------------------------

diff_df <- incov %>%
  transmute(pos = as.integer(V2), incov = as.numeric(V3)) %>%
  inner_join(
    outcov %>% transmute(pos = as.integer(V2), outcov = as.numeric(V3)),
    by = "pos"
  ) %>%
  mutate(diff = incov - outcov)

## --- Variant annotations -----------------------------------------------------

snv1 <- snv %>%
  transmute(pos = as.integer(V2)) %>%
  inner_join(diff_df %>% select(pos, diff, incov, outcov), by = "pos")

snv1$minVAF <- round( (snv1$incov*VAF/100 - snv1$diff)/snv1$outcov, 3)
snv1$maxVAF <- round( (snv1$incov*VAF/100)/snv1$outcov, 3)

indel1 <- indel %>%
  transmute(
    start = as.integer(V2) - 1L,
    end   = as.integer(V3) + 1L,
    type  = as.character(V5)
  ) %>%
  mutate(type = ifelse(is.na(type) | type == "", "INDEL", type))

## --- Remove dropped variant -----------------------------------------------------

library(vcfR)
# SNVs
snv_file <- paste0(work_dir,"/vcf/",ID,".chrM.40snv.af",VAF,".raw.addsnv.target.snv.",VAF,".vcf")
snv_vcf <- data.frame(read.vcfR(snv_file)@fix)
snv_df <- snv1[(snv1$pos %in% snv_vcf$POS),]
  
# INDELs
indel_file <- paste0(work_dir,"/vcf/",ID,".chrM.10indel.af",VAF,".raw.addindel.target.indel.",VAF,".vcf")
indel_vcf <- data.frame(read.vcfR(indel_file)@fix)
indel_df <- indel1 %>%
  rowwise() %>%
  mutate(
    has_vcf = any(indel_vcf$CHROM == "chrM" &
                    indel_vcf$POS >= start &
                    indel_vcf$POS <= end)
  ) %>%
  ungroup() %>%
  filter(has_vcf)

## --- Files -------------------------------------------------------------------

sim_dir <- "/vast/projects/bahlo_mtDNA/Simulation/"
mity_file        <- paste0(sim_dir, "mity/",ID,"_AF",VAF,"/",ID,".mity.report.xlsx")
server2_file     <- paste0(sim_dir, "mtdna-server-2/",ID,"_AF",VAF,"/results/variants.annotated.txt")
mitoHPC_file     <- paste0(sim_dir, "mitoHPC/",ID,"_AF",VAF,"/out/",ID,".af",VAF,"/",ID,".af",VAF,".mutect2.mutect2.00.vcf")
mtSwirl_file <- Sys.glob(paste0(sim_dir,"mtSwirl/",ID,"_AF",VAF,"/*_MitochondriaPipeline/out/final_vcf/",ID,".self.ref.split.selfToRef.final.vcf"))

artifact_prone_sites <- c(301, 302, 310, 316, 3107, 16179:16183)

## Helper: indel-window selector (±1bp around each indel span)
in_indel_window <- function(pos_vec, indel_start, indel_end, pad = 1L) {
  Reduce(`|`, Map(function(s, e) pos_vec >= (s - pad) & pos_vec <= (e + pad), indel_start, indel_end))
}

### mtDNA-Server2 --------------------------------------------------------------

server2 <- read.table(server2_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

server2$server2_VAF <- as.numeric(server2$VariantLevel)
server2$server2_VAR <- paste(server2$Pos, server2$Ref, server2$Variant, sep = "_")
server2$server2_CVG <- server2$Coverage

server2$server2_QC <- NA_character_
keep1 <- !(server2$Pos %in% artifact_prone_sites) & !grepl("strand_bias", server2$Filter)
server2$server2_QC[keep1] <- "PASS"

server2_snv <- merge(
  snv_df,
  server2[, c("Pos", "server2_VAF", "server2_CVG", "server2_QC")],
  by.x = "pos", by.y = "Pos", all.x = TRUE
)

keep2 <- in_indel_window(server2$Pos, as.integer(indel$V2), as.integer(indel$V3), pad = 1L)
server2_indel <- server2[keep2, c("server2_VAR", "Pos", "server2_VAF", "server2_CVG", "server2_QC")]

server2_indel_no_snv <- server2_indel %>%
  separate(server2_VAR, into = c("pos2","REF","ALT"), sep = "_", remove = FALSE) %>%
  filter(nchar(REF) != nchar(ALT)) %>%   # removes substitutions like 12351_T_A
  select(-pos2, -REF, -ALT)

server2_indel <- server2_indel_no_snv

### mity ----------------------------------------------------------------------

mity <- readxl::read_excel(mity_file)

mity$mity_VAF <- round(as.numeric(mity$`VARIANT HETEROPLASMY`), 3)   # numeric for now
mity$mity_VAR <- paste(mity$POS, mity$REF, mity$ALT, sep = "_")
mity$mity_CVG <- mity$`TOTAL LOCUS DEPTH`

mity$mity_QC <- NA_character_
keep1 <- !(mity$POS %in% artifact_prone_sites) & (mity$FILTER == "PASS")
mity$mity_QC[keep1] <- "PASS"

snv_df <- server2_snv
mity_snv <- merge(
  snv_df,
  mity[, c("POS", "mity_VAF", "mity_CVG", "mity_QC")],
  by.x = "pos", by.y = "POS", all.x = TRUE
)

keep2 <- in_indel_window(mity$POS, as.integer(indel$V2), as.integer(indel$V3), pad = 1L)
mity_indel <- mity[keep2, c("mity_VAR", "POS", "mity_VAF", "mity_CVG", "mity_QC")]

mity_indel_no_snv <- mity_indel %>%
  separate(mity_VAR, into = c("pos2","REF","ALT"), sep = "_", remove = FALSE) %>%
  filter(
    (nchar(REF) == 1 & nchar(ALT) > 1) |
    (nchar(REF) > 1 & nchar(ALT) == 1)
  ) %>%
  select(-pos2, -REF, -ALT) # update: edit filtering step to include reads where ref or alt is 1 (not both either)

mity_indel <- mity_indel_no_snv

### mitoHPC -------------------------------------------------------------------

mitoHPC_vcf <- read.vcfR(mitoHPC_file)
mitoHPC <- data.frame(mitoHPC_vcf@fix, stringsAsFactors = FALSE)
mitoHPC_gt  <- data.frame(mitoHPC_vcf@gt,  stringsAsFactors = FALSE)

mitoHPC$POS <- as.integer(mitoHPC$POS)
mitoHPC$mitoHPC_VAF <- round(as.numeric(sapply(mitoHPC_gt[, 2], function(x) strsplit(x, ":")[[1]][3])),3)
mitoHPC$mitoHPC_VAR <- paste(mitoHPC$POS, mitoHPC$REF, mitoHPC$ALT, sep = "_")
mitoHPC$mitoHPC_CVG <- as.numeric(sapply(mitoHPC_gt[, 2], function(x) strsplit(x, ":")[[1]][2]))

mitoHPC$mitoHPC_QC <- NA_character_
keep1 <- !(mitoHPC$POS %in% artifact_prone_sites) & !grepl("strand_bias", mitoHPC$FILTER)
mitoHPC$mitoHPC_QC[keep1] <- "PASS"

snv_df <- mity_snv
mitoHPC_snv <- merge(
  snv_df,
  mitoHPC[, c("POS", "mitoHPC_VAF", "mitoHPC_CVG", "mitoHPC_QC")],
  by.x = "pos", by.y = "POS", all.x = TRUE
)

keep2 <- in_indel_window(as.integer(mitoHPC$POS), as.integer(indel$V2), as.integer(indel$V3), pad = 1L)
mitoHPC_indel <- mitoHPC[keep2, c("mitoHPC_VAR", "POS", "mitoHPC_VAF","mitoHPC_CVG", "mitoHPC_QC")]

### mtSwirl -------------------------------------------------------------------

mtSwirl_vcf <- read.vcfR(mtSwirl_file)
mtSwirl <- data.frame(mtSwirl_vcf@fix, stringsAsFactors = FALSE)[, c(1, 2, 4, 5, 7)]
mtSwirl.gt <- data.frame(mtSwirl_vcf@gt, stringsAsFactors = FALSE)

# Columns should be: CHROM POS REF ALT FILTER (depending on vcfR version)
mtSwirl$POS <- as.integer(mtSwirl$POS)
mtSwirl$mtSwirl_VAF <- round(as.numeric(sapply(mtSwirl.gt[, 2], function(x) strsplit(x, ":")[[1]][2])),3)
mtSwirl$mtSwirl_VAR <- paste(mtSwirl$POS, mtSwirl$REF, mtSwirl$ALT, sep = "_")
mtSwirl$mtSwirl_CVG <- as.numeric(sapply(mtSwirl.gt[, 2], function(x) strsplit(x, ":")[[1]][3]))

mtSwirl$mtSwirl_QC <- NA_character_
keep1 <- !(mtSwirl$POS %in% artifact_prone_sites) & !grepl("strand_bias|blacklisted_site|base_qual|position|FAIL|map_qual", mtSwirl$FILTER) # update: add FAIL and map_qual filter
mtSwirl$mtSwirl_QC[keep1] <- "PASS"

snv_df <- mitoHPC_snv
mtSwirl_snv <- merge(
  snv_df,
  mtSwirl[, c("POS", "mtSwirl_VAF", "mtSwirl_CVG", "mtSwirl_QC")],
  by.x = "pos", by.y = "POS", all.x = TRUE
)

keep2 <- in_indel_window(mtSwirl$POS, as.integer(indel$V2), as.integer(indel$V3), pad = 1L)
mtSwirl_indel <- mtSwirl[keep2, c("mtSwirl_VAR", "POS", "mtSwirl_VAF","mtSwirl_CVG", "mtSwirl_QC")]
mtSwirl_indel_no_snv <- mtSwirl_indel %>%
  separate(mtSwirl_VAR, into = c("pos2","REF","ALT"), sep = "_", remove = FALSE) %>%
  filter(nchar(REF) != nchar(ALT)) %>%   
  select(-pos2, -REF, -ALT) # update: added SNV filting step for mtSwirl

mtSwirl_indel <- mtSwirl_indel_no_snv

write.table(mtSwirl_snv,  paste0(work_dir,"all_tools_snv.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## --- Combine indels (NO coverage) into LONG table ----------------------------

all_indel_long <- bind_rows(
  server2_indel %>%
    transmute(tool = "server2",
              POS  = as.integer(Pos),
              VAR  = server2_VAR,
              VAF  = sprintf("%.3f", as.numeric(server2_VAF)),
              CVG  = as.numeric(server2_CVG),
              QC   = server2_QC),

  mity_indel %>%
    transmute(tool = "mity",
              POS  = as.integer(POS),
              VAR  = mity_VAR,
              VAF  = sprintf("%.3f", as.numeric(mity_VAF)),
              CVG  = as.numeric(mity_CVG),
              QC   = mity_QC),
  
  mitoHPC_indel %>%
    transmute(tool = "mitoHPC",
              POS  = as.integer(POS),
              VAR  = mitoHPC_VAR,
              VAF  = sprintf("%.3f", as.numeric(mitoHPC_VAF)),
              CVG  = as.numeric(mitoHPC_CVG),
              QC   = mitoHPC_QC),
  
  mtSwirl_indel %>%
    transmute(tool = "mtSwirl",
              POS  = as.integer(POS),
              VAR  = mtSwirl_VAR,
              VAF  = sprintf("%.3f", as.numeric(mtSwirl_VAF)),
              CVG  = as.numeric(mtSwirl_CVG),
              QC   = mtSwirl_QC),
) %>%
  arrange(POS, tool)

## --- Filter for real variants  ------------------
# update: added this filter
indeltab <- read.table(
  paste0(work_dir,"target.indel.",VAF,".txt"),
  sep = "\t", header = FALSE, fill = TRUE, stringsAsFactors = FALSE
)

target_indel_truth <- indeltab %>%
  transmute(
    POS = as.integer(V2),
    type = V5,  # "DEL" or "INS"
    
    # deletion
    true_del_len = ifelse(V5 == "DEL", as.integer(V3 - V2), NA_integer_),
    true_del_seq = ifelse(V5 == "DEL" & V6 != "", V6, NA_character_),
    
    # insertion
    true_ins_len = ifelse(V5 == "INS" & V6 != "", nchar(V6), NA_integer_),
    true_ins_seq = ifelse(V5 == "INS" & V6 != "", V6, NA_character_)
  )

all_indel_long <- all_indel_long %>%
  separate(
    VAR,
    into = c("pos_tmp", "REF_tmp", "ALT_tmp"),
    sep = "_",
    remove = FALSE
  ) %>%
  mutate(
    POS = as.integer(POS),
    
    ref_len = nchar(REF_tmp),
    alt_len = nchar(ALT_tmp),
    
    # deletion payload (REF minus anchor)
    del_seq = ifelse(ref_len > alt_len, substring(REF_tmp, 2), NA_character_),
    
    # insertion payload (ALT minus anchor)
    ins_seq = ifelse(alt_len > ref_len, substring(ALT_tmp, 2), NA_character_)
  ) %>%
  left_join(target_indel_truth, by = "POS") %>%
  
  group_by(tool, POS) %>%
  filter(
    # keep single calls unchanged
    n() == 1 |
      
      # multi-allelic case → validate against truth
      (n() > 1 & (
        
        ## ---- DELETIONS ----
        (type == "DEL" & (
          (!is.na(true_del_seq) & del_seq == true_del_seq) |
            (is.na(true_del_seq) & !is.na(true_del_len) & (ref_len-1) == true_del_len)
        )) |
          
          ## ---- INSERTIONS ----
        (type == "INS" & (
          (!is.na(true_ins_seq) & ins_seq == true_ins_seq) |
            (is.na(true_ins_seq) & !is.na(true_ins_len) & (alt_len - 1) == true_ins_len)
        ))
        
      ))
  ) %>%
  ungroup() %>%
  select(
    -pos_tmp, -REF_tmp, -ALT_tmp,
    -ref_len, -alt_len,
    -del_seq, -ins_seq,
    -true_del_len, -true_del_seq,
    -true_ins_len, -true_ins_seq,
    -type
  )


## --- Merge by POS into WIDE table --------------------------------------------

all_indel_wide <- all_indel_long %>%
  pivot_wider(
    id_cols = POS,
    names_from = tool,
    values_from = c(VAR, VAF, CVG, QC),
    names_glue = "{tool}_{.value}"
  ) %>%
  arrange(POS)

# all_indel <- all_indel_wide[,c("POS","server2_VAR","server2_VAF","server2_CVG","server2_QC",
#                                "mity_VAR","mity_VAF","mity_CVG","mity_QC",
#                                "mitoHPC_VAR","mitoHPC_VAF","mitoHPC_CVG","mitoHPC_QC",
#                                "mtSwirl_VAR","mtSwirl_VAF","mtSwirl_CVG","mtSwirl_QC")]
# update: fills missing columns with NA to deal with cases when tool doesn't pick up any variant calls (e.g. HG00428 VAF=0.01)
all_indel <- all_indel_wide %>%
  dplyr::select(
    POS,
    dplyr::any_of(c(
      "server2_VAR","server2_VAF","server2_CVG","server2_QC",
      "mity_VAR","mity_VAF","mity_CVG","mity_QC",
      "mitoHPC_VAR","mitoHPC_VAF","mitoHPC_CVG","mitoHPC_QC",
      "mtSwirl_VAR","mtSwirl_VAF","mtSwirl_CVG","mtSwirl_QC"
    ))
  )
required_cols <- c(
  "server2_VAR","server2_VAF","server2_CVG","server2_QC",
  "mity_VAR","mity_VAF","mity_CVG","mity_QC",
  "mitoHPC_VAR","mitoHPC_VAF","mitoHPC_CVG","mitoHPC_QC",
  "mtSwirl_VAR","mtSwirl_VAF","mtSwirl_CVG","mtSwirl_QC"
)

missing_cols <- setdiff(required_cols, names(all_indel))
all_indel[missing_cols] <- NA

close_pairs <- all_indel %>%
  select(POS, server2_VAR, mity_VAR, mitoHPC_VAR, mtSwirl_VAR) %>%
  arrange(POS) %>%
  mutate(next_POS = lead(POS),
         dist_to_next = next_POS - POS) %>%
  filter(!is.na(dist_to_next) & dist_to_next <= 10) 

if(nrow(close_pairs) >0){
  check_pair <- all_indel[all_indel$POS %in% c(close_pairs$POS, close_pairs$next_POS),
                          c("POS", "server2_VAR", "mity_VAR", "mitoHPC_VAR", "mtSwirl_VAR")]
  # rows where ONLY ONE tool has a non-NA (and non-empty) value
  x <- check_pair[, -1, drop = FALSE]
  nonempty <- !is.na(x) & x != ""
  n_nonNA <- rowSums(nonempty)
  
  check_pair_single <- check_pair[n_nonNA == 1, , drop = FALSE]
  
  if(nrow(check_pair_single) >0) all_indel <- all_indel[!(all_indel$POS %in% check_pair_single$POS),]
}

indel_df <- diff_df[diff_df$pos %in% all_indel$POS,]
indel_df$minVAF <- round( (indel_df$incov*VAF/100 - indel_df$diff)/indel_df$outcov, 3)
indel_df$maxVAF <- round( (indel_df$incov*VAF/100)/indel_df$outcov, 3)

indel_anno <- merge(indel_df, all_indel, by.x = "pos", by.y = "POS")

write.table(indel_anno,  paste0(work_dir,"all_tools_indel.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

