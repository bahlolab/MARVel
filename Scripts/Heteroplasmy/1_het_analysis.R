##### Part A - Extract heteroplasmic variants #####
rm(list = ls())
setwd("/vast/projects/bahlo_mtDNA/Results_1000G/")
library(readxl)
library(vcfR)
library(dplyr)
library(stringr)
library(readr)
library(purrr)

artifact_prone_sites <- c(301,302,310,316,3107,5894,10933,16179,16181,16182,16183,16188,16189,16192)

base_dir <- "/vast/projects/bahlo_mtDNA/Results_1000G"
id_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
id_dirs <- id_dirs[grepl("^(HG|NA)", basename(id_dirs))]
all_ids = basename(id_dirs)

allids = read.csv("/vast/projects/bahlo_mtDNA/Scripts_Kelly/sample.tsv",sep = "\t")
all_ids = allids$library_name
pedigree = read.csv("/vast/projects/bahlo_mtDNA/Scripts_Kelly/pedigree.tsv",sep = "\t")
ped <- pedigree %>% 
  filter(MotherID != 0) %>% 
  filter(entity.pedigree_id %in% all_ids) %>% 
  filter(MotherID %in% all_ids) %>% 
  select(entity.pedigree_id, MotherID) %>% 
  dplyr::rename(mother = MotherID,
                proband = entity.pedigree_id)
ids = sort(c(ped$mother,ped$proband))
idsu = unique(ids)
id_dirs <- id_dirs[basename(id_dirs) %in% ids] 

#---------- mtDNA-server2 ----------#
process_annotated <- function(id_dir) {
  file_path <- file.path(id_dir, "/mtdna-server-2/results/variants.annotated.txt")
  if (!file.exists(file_path)) {
    message("Missing file for ID: ", file_path)
    return(NULL)}
  
  df <- readr::read_tsv(file_path, col_types = readr::cols(.default = "c"))
}

mtdna_server <- map_dfr(id_dirs, process_annotated)
mtdna_server <- mtdna_server[!(mtdna_server$Pos %in% artifact_prone_sites),]
mtdna_server$VAF <- as.numeric(mtdna_server$VariantLevel)
mtdna_server$SAMPLE <- sub("\\.chrM.*", "", mtdna_server$ID)
mtdna_server$VAR <- paste(mtdna_server$Pos,mtdna_server$Ref,mtdna_server$Variant,sep = "_")
mtdna_server <- mtdna_server %>%
  filter(!grepl("strand_bias", Filter))

#------------ mity ------------#
process_mity <- function(id_dir) {
  file_path <- file.path(id_dir, paste0("/mity/",basename(id_dir), ".mity.report.xlsx"))
  if (!file.exists(file_path)) {
    message("Missing file for ID: ", file_path)
    return(NULL)}
  df <- read_excel(file_path)
}

mity <- map_dfr(id_dirs, process_mity)
mity <- mity %>% 
  filter(FILTER == "PASS" | (FILTER == "FAIL" & POS_FILTER == "0" & SBR_FILTER == "1" & SBA_FILTER == "1" & MQMR_FILTER == "1" & AQR_FILTER == "1"))
mity <- mity[!(mity$POS %in% artifact_prone_sites),]
mity$VAF <- as.numeric(mity$`VARIANT HETEROPLASMY`)
mity$VAR <- paste(mity$POS,mity$REF,mity$ALT,sep = "_")

mity_concise <- mity %>% select(SAMPLE, HGVS, `GENE/LOCUS`, `VARIANT HETEROPLASMY`, 
                                CHR, POS, REF, ALT, FILTER, VAF, VAR,GT_FORMAT,
                                FILTER,INFO,FORMAT)

#------------ mitoHPC ------------#
process_vcf <- function(id_dir) {
  file_path <- file.path(id_dir, "mitoHPC", "mutect2.mutect2.03.vcf")
  if (!file.exists(file_path)) return(NULL)
  
  vcf <- read.vcfR(file_path)
  vcf_fix <- data.frame(vcf@fix, stringsAsFactors = FALSE)
  vcf_gt  <- data.frame(vcf@gt, stringsAsFactors = FALSE)
  
  df <- cbind(vcf_fix, SAMPLE = sub("\\.final", "", vcf_gt[,2])) %>%
    mutate(
      POS  = as.numeric(POS),
      VAF  = as.numeric(str_extract(INFO, "(?<=AF=)[^;]+")),
      VAR  = paste(POS, REF, ALT, sep = "_")
    ) %>%
    filter(!POS %in% artifact_prone_sites) %>%
    filter(!grepl("strand_bias", FILTER)) %>% 
    select("POS","REF","ALT","FILTER","INFO","SAMPLE","VAF","VAR")
  
  return(df)
}

mitoHPC <- map_dfr(id_dirs, process_vcf)

#------------ tool-defined heteroplasmy ------------#
mtdna_server_het <- mtdna_server[ (mtdna_server$VAF< 0.95) 
                                  & mtdna_server$Ref != mtdna_server$Variant,
                                  c("SAMPLE","VAR","VAF","Filter")]
mity_het <- mity[mity$VAF < 0.95 & (mity$REF != mity$ALT) ,c("SAMPLE","VAR","VAF","FILTER")]

mitoHPC_het <- mitoHPC[mitoHPC$VAF < 0.95 & (mitoHPC$REF != mitoHPC$ALT), c("SAMPLE","VAR","VAF","FILTER")]

data <- data.frame(SAMPLE=NA, VAR=NA, VAF=NA, FILTER=NA)
nsample = 1200 
for(i in 1:nsample){
  
  id <- idsu[i]
  
  ### mtSwirl
  mtSwirl_id_vcf <- read.vcfR(paste0("/vast/projects/bahlo_mtDNA/Results_1000G/", id, "/mtSwirl/",id,".self.ref.split.selfToRef.final.vcf"))
  mtSwirl_id <- data.frame(mtSwirl_id_vcf@fix)[,c(1,2,4,5,7)]
  
  mtSwirl_id.gt <- data.frame(mtSwirl_id_vcf@gt)
  mtSwirl_id$VAF <- as.numeric(sapply(mtSwirl_id.gt[,2], function(x) { strsplit(x, ":")[[1]][2]}))
  mtSwirl_id$VAR <- paste(mtSwirl_id$POS,mtSwirl_id$REF,mtSwirl_id$ALT,sep = "_")
  mtSwirl_id <- mtSwirl_id[!(mtSwirl_id$POS %in% artifact_prone_sites),]
  
  mtSwirl_id_het <- mtSwirl_id[mtSwirl_id$VAF < 0.95 & mtSwirl_id$REF != mtSwirl_id$ALT,c("VAR","VAF","FILTER")]
  
  if(nrow(mtSwirl_id_het)>0) {
    mtSwirl_id_het$SAMPLE <- id
    data <- rbind(data, mtSwirl_id_het[,c("SAMPLE","VAR","VAF","FILTER")])
  }
}

mtSwirl_het <- data[-1,]
mtSwirl_het <- mtSwirl_het %>%
  filter(!grepl("strand_bias|blacklisted_site|base_qual|position|FAIL|map_qual", FILTER)) # update: add FAIL & map_qual filters
mtSwirl_het$TOOL <- "mtSwirl"
mtdna_server_het$TOOL <- "mtDNA Server2"
colnames(mtdna_server_het)[4] <- "FILTER"
mity_het$TOOL <- "mity"
mitoHPC_het$TOOL <- "mitoHPC"

### table of heteroplasmic variants across 4 tools
sum <- rbind(mity_het, mtdna_server_het, mitoHPC_het, mtSwirl_het)
rownames(sum) <- NULL
write.table(sum, "/vast/projects/bahlo_mtDNA/Analysis/heteroplasmy/allduos/trio604_summary_het_newfilters.tsv", 
            row.names = F, quote = F, sep = "\t") 

##### Part B - Heteroplasmy Analysis ######
rm(list = ls())
setwd("/vast/projects/bahlo_mtDNA/Analysis/heteroplasmy/allduos")
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggstatsplot)
library(stringr)

data <- read.table("trio604_summary_het_newfilters.tsv", header = TRUE, sep = "\t")
setwd("/vast/projects/bahlo_mtDNA/Analysis/heteroplasmy/allduos")

# extract heteroplasmic variants using definition as 3% <= VAF < 95%
dat03 <- data[data$VAF>=0.03 & data$VAF < 0.95,]
dat03 <- na.omit(dat03)
rownames(dat03) <- NULL
write.table(dat03, "./trio604_vaf03_het_newfilters.tsv", row.names = F, quote = F, sep = "\t") # table of VAF>0.03 het variants

# Check heteroplasmic MNVs in mity 
mity_het_mnv <- mity_het %>%
  separate(VAR, into = c("pos", "ref", "alt"), sep = "_", remove = FALSE) %>%
  filter(nchar(ref) >= 2 & nchar(alt) >= 2)

# Count number of tools detecting each variant per individual
tool_counts <- dat03 %>%
  group_by(SAMPLE, VAR) %>%
  summarise(
    tools_detected = n_distinct(TOOL),
    tools_list = paste(sort(unique(TOOL)), collapse = ", "),
    .groups = "drop"
  )

# Filter variant–sample combinations with exactly 3 tools
variants_3tools <- tool_counts %>%
  filter(tools_detected == 3) %>%
  select(SAMPLE, VAR)  


variant_wide <- dat03 %>%
  semi_join(variants_3tools, by = c("SAMPLE", "VAR")) %>%
  group_by(VAR, TOOL) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = TOOL, values_from = count, values_fill = 0)

variant_wide$pos <- as.numeric(str_extract(variant_wide$VAR, "^[0-9]+"))
variant_wide <- variant_wide %>% arrange(pos) %>% 
  select(VAR,mity,`mtDNA Server2`,mitoHPC,mtSwirl,pos)

sample_list <- variants_3tools %>%
  group_by(VAR) %>%
  summarise(
    samples = paste(sort(unique(SAMPLE)), collapse = ", "),
    .groups = "drop"
  )

# table including sample ID
variant_wide2 <- variant_wide %>%
  left_join(sample_list, by = "VAR") %>%
  select(VAR, samples, mity, `mtDNA Server2`, mitoHPC, mtSwirl, pos)

write.table(variant_wide, "./variants_3tools_newfilters.tsv", row.names = F, quote = F, sep = "\t")
write.table(variant_wide2, "./variants_3tools_withsample_newfilters.tsv", row.names = F, quote = F, sep = "\t")
write.table(mity_het_mnv, "./mity_het_mnv_newfilters.tsv",row.names = F, quote = F, sep = "\t")

#### Part C - Statistics ####
rm(list=ls())
library(dplyr)
library(rstatix)

setwd("/vast/projects/bahlo_mtDNA/Analysis/heteroplasmy/allduos")
data <- read.table("trio604_vaf03_het_newfilters.tsv", header = TRUE, sep = "\t")

overlaps <- data %>%
  group_by(SAMPLE, VAR) %>%
  mutate(count = n()) %>%
  filter(count == 4) %>%
  ungroup() %>%
  select(-count)

overlaps$pair <- paste0(overlaps$SAMPLE,"_",overlaps$VAR)
overlaps <- overlaps %>%
  mutate(TOOL = factor(TOOL, levels = c("mitoHPC", "mtSwirl", "mity", "mtDNA Server2")))

overlaps <- overlaps %>%
  mutate(pair = as.factor(pair))

friedman_res <- overlaps %>%
  friedman_test(VAF ~ TOOL | pair)
friedman_res