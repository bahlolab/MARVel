rm(list = ls())
library(dplyr)

setwd("/vast/scratch/users/chen.k")

bam_dirs <- c(
  "/stornext/Bioinf/data/lab_bahlo/projects/ataxia/GRREAT/bams",
  "/stornext/Bioinf/data/lab_bahlo/projects/ataxia/GRREAT/bams/macrogen/bam"
)

# get bam files
bam <- unlist(lapply(
  bam_dirs,
  list.files,
  pattern = "\\.bam$",
  full.names = TRUE
))

bam_names <- basename(bam)

sample <- sub("\\.merged\\.bam$", "", bam_names)

bai <- paste0(bam, ".bai") # expected bai paths
bai[!file.exists(bai)] <- NA # replace missing bai files with NA


input <- data.frame(
  sample = sample,
  bam = bam,
  bai = bai,
  stringsAsFactors = FALSE
)

no_bai_files <- input %>% filter(is.na(bai)) # check if any bai files are not present

dup_samples <- unique(input$sample[duplicated(input$sample) | duplicated(input$sample, fromLast = TRUE)])

input$sample <- ifelse(
  input$sample %in% dup_samples &
    !grepl("macrogen", input$bam, ignore.case = TRUE),
  paste0(input$sample, "-vcgs"),
  input$sample
)
filtered_input <- input %>% filter(!is.na(bai))

write.table(filtered_input, "input_list_complete.txt", row.names = F, quote = F, sep = "\t")
