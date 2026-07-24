# rm(list = ls())
# ID <- "HG01881"
# VAF <- 10

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop('Usage: Rscript 4.3_FP.R "HG01881" 10')
}
ID <- args[1]
VAF <- as.integer(args[2])

suppressPackageStartupMessages({
  library(vcfR)
})

work_dir <- paste0("/vast/projects/bahlo_mtDNA/Simulation/1000G/",ID,"/AF",VAF,"_DS2k/")
sim_dir <- "/vast/projects/bahlo_mtDNA/Simulation/"
res_dir <- "/vast/projects/bahlo_mtDNA/Results_1000G/"

snv <- read.table(paste0(work_dir,"all_tools_snv.tsv"), header = T, sep = "\t")
indel <- read.table(paste0(work_dir,"all_tools_indel.tsv"), header = T, sep = "\t")

artifact_prone_sites <- c(301, 302, 310, 316, 3107, 16179:16183)

## --- Files -------------------------------------------------------------------

mity_file        <- paste0(sim_dir, "mity/",ID,"_AF",VAF,"/",ID,".mity.report.xlsx")
server2_file     <- paste0(sim_dir, "mtdna-server-2/",ID,"_AF",VAF,"/results/variants.annotated.txt")
mitoHPC_file     <- paste0(sim_dir, "mitoHPC/",ID,"_AF",VAF,"/out/",ID,".af",VAF,"/",ID,".af",VAF,".mutect2.mutect2.00.vcf")
mtSwirl_file <- Sys.glob(paste0(sim_dir,"mtSwirl/",ID,"_AF",VAF,"/*_MitochondriaPipeline/out/final_vcf/",ID,".self.ref.split.selfToRef.final.vcf"))

## --- Original Files -------------------------------------------------------------------

mity_file0 <- paste0(sim_dir, "1000G/",ID,"/mity/",ID,".mity.report.xlsx")
server2_file0 <- paste0(sim_dir, "1000G/",ID,"/mtdna-server-2/results/variants.annotated.txt")
mitoHPC_file0 <- paste0(res_dir,ID,"/mitoHPC/",ID,".final.mutect2.mutect2.00.vcf")
mtSwirl_file0 <- paste0(res_dir,ID,"/mtSwirl/",ID,".self.ref.split.selfToRef.final.vcf")

## Helper: indel-window selector (±1bp around each indel span)
in_indel_window <- function(pos_vec, indel_start, indel_end, pad = 10L) {
  Reduce(`|`, Map(function(s, e) pos_vec >= (s - pad) & pos_vec <= (e + pad), indel_start, indel_end))
}

indel_df <- read.table(
  paste0(work_dir,"target.indel.",VAF,".txt"),
  sep = "\t", header = FALSE, fill = TRUE, stringsAsFactors = FALSE
)

### mtDNA-Server2 --------------------------------------------------------------

server2_new <- read.table(server2_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
server2_org <- read.table(server2_file0, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

server2_add <- server2_new[!(server2_new$Mutation %in% server2_org$Mutation),]
server2_add$VAR <- paste(server2_add$Pos, server2_add$Ref, server2_add$Variant, sep = "_")
server2 <- server2_add[!(server2_add$Pos %in% snv$pos) & !(server2_add$VAR %in% indel$server2_VAR),]
in_indel <- in_indel_window(server2$Pos, as.integer(indel_df$V2), as.integer(indel_df$V3), pad = 10L)
in_snv <- in_indel_window(server2$Pos, as.integer(snv$pos), as.integer(snv$pos), pad = 10L)
server2$FP_source <- "unknown"; server2$FP_source[in_indel] <- "indel"; server2$FP_source[in_snv] <- "snv"

keep1 <- !(server2$Pos %in% artifact_prone_sites) & !grepl("strand_bias", server2$Filter)
server2_FP <- server2[keep1, c("VAR","VariantLevel","FP_source")]

### mity ----------------------------------------------------------------------

mity_new <- readxl::read_excel(mity_file)
mity_org <- readxl::read_excel(mity_file0)

mity_add <- mity_new[!(mity_new$HGVS %in% mity_org$HGVS),]
mity_add$VAR <- paste(mity_add$POS, mity_add$REF, mity_add$ALT, sep = "_")
mity <- mity_add[!(mity_add$POS %in% snv$pos) & !(mity_add$VAR %in% indel$mity_VAR),]
mity$VAF <- round(as.numeric(mity$`VARIANT HETEROPLASMY`),3)
in_indel <- in_indel_window(mity$POS, as.integer(indel_df$V2), as.integer(indel_df$V3), pad = 10L)
in_snv <- in_indel_window(mity$POS, as.integer(snv$pos), as.integer(snv$pos), pad = 10L)
mity$FP_source <- "unknown"; mity$FP_source[in_indel] <- "indel"; mity$FP_source[in_snv] <- "snv"

keep2 <- !(mity$POS %in% artifact_prone_sites) & (mity$FILTER == "PASS")
mity_FP <- mity[keep2, c("VAR","VAF","FP_source")]

### mitoHPC -------------------------------------------------------------------

mitoHPC_vcf <- read.vcfR(mitoHPC_file)
mitoHPC_new <- data.frame(mitoHPC_vcf@fix, stringsAsFactors = FALSE)
mitoHPC_org <- data.frame(read.vcfR(mitoHPC_file0)@fix, stringsAsFactors = FALSE)
mitoHPC_new$VAR <- paste(mitoHPC_new$POS, mitoHPC_new$REF, mitoHPC_new$ALT, sep = "_")
mitoHPC_org$VAR <- paste(mitoHPC_org$POS, mitoHPC_org$REF, mitoHPC_org$ALT, sep = "_")
mitoHPC_gt  <- data.frame(mitoHPC_vcf@gt,  stringsAsFactors = FALSE)
mitoHPC_new$VAF <- round(as.numeric(sapply(mitoHPC_gt[, 2], function(x) strsplit(x, ":")[[1]][3])),3)

mitoHPC_add <- mitoHPC_new[!(mitoHPC_new$VAR %in% mitoHPC_org$VAR),]
mitoHPC <- mitoHPC_add[!(mitoHPC_add$POS %in% snv$pos) & !(mitoHPC_add$VAR %in% indel$mitoHPC_VAR),]

if(nrow(mitoHPC)){
  in_indel <- in_indel_window(mitoHPC$POS, as.integer(indel_df$V2), as.integer(indel_df$V3), pad = 10L)
  in_snv <- in_indel_window(mitoHPC$POS, as.integer(snv$pos), as.integer(snv$pos), pad = 10L)
  mitoHPC$FP_source <- "unknown"; mitoHPC$FP_source[in_indel] <- "indel"; mitoHPC$FP_source[in_snv] <- "snv"
  
  keep3 <- !(mitoHPC$POS %in% artifact_prone_sites) & !grepl("strand_bias", mitoHPC$FILTER)
  mitoHPC_FP <- mitoHPC[keep3, c("VAR","VAF","FP_source")]
}else{
  mitoHPC_FP <- NA
}

### mtSwirl -------------------------------------------------------------------

mtSwirl_vcf <- read.vcfR(mtSwirl_file)
mtSwirl_new <- data.frame(mtSwirl_vcf@fix, stringsAsFactors = FALSE)[, c(1, 2, 4, 5, 7)]
mtSwirl_org <- data.frame(read.vcfR(mtSwirl_file0)@fix, stringsAsFactors = FALSE)[, c(1, 2, 4, 5, 7)]
mtSwirl_new$VAR <- paste(mtSwirl_new$POS, mtSwirl_new$REF, mtSwirl_new$ALT, sep = "_")
mtSwirl_org$VAR <- paste(mtSwirl_org$POS, mtSwirl_org$REF, mtSwirl_org$ALT, sep = "_")
mtSwirl.gt <- data.frame(mtSwirl_vcf@gt, stringsAsFactors = FALSE)
mtSwirl_new$VAF <- round(as.numeric(sapply(mtSwirl.gt[, 2], function(x) strsplit(x, ":")[[1]][2])),3)

mtSwirl_add <- mtSwirl_new[!(mtSwirl_new$VAR %in% mtSwirl_org$VAR),]
mtSwirl <- mtSwirl_add[!(mtSwirl_add$POS %in% snv$pos) & !(mtSwirl_add$VAR %in% indel$mtSwirl_VAR),]
in_indel <- in_indel_window(mtSwirl$POS, as.integer(indel_df$V2), as.integer(indel_df$V3), pad = 10L)
in_snv <- in_indel_window(mtSwirl$POS, as.integer(snv$pos), as.integer(snv$pos), pad = 10L)
mtSwirl$FP_source <- "unknown"; mtSwirl$FP_source[in_indel] <- "indel"; mtSwirl$FP_source[in_snv] <- "snv"

keep4 <- !(mtSwirl$POS %in% artifact_prone_sites) & !grepl("FAIL|strand_bias|blacklisted_site|base_qual|position|map_qual", mtSwirl$FILTER) # update: add map_qual filter
mtSwirl_FP <- mtSwirl[keep4, c("VAR","VAF","FP_source")]

### FP list --------------------------------------------------------------

FP_list <- list(server2_FP, mity_FP, mitoHPC_FP, mtSwirl_FP)
saveRDS(FP_list, file = paste0(work_dir,"FP_list.rds"))

