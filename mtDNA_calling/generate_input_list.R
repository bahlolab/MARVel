rm(list = ls())
library(dplyr)

setwd(paste0("/vast/scratch/users/", Sys.getenv("USER")))

bam_dirs <- c("/path/to/bam/1", "/path/to/bam/2")

bam <- unlist(lapply(bam_dirs, list.files, pattern = "\\.bam$", full.names = TRUE))
bam_names <- basename(bam)
sample <- sub("\\.merged\\.bam$", "", bam_names)

bai <- paste0(bam, ".bai")
bai[!file.exists(bai)] <- NA

input <- data.frame(sample = sample, bam = bam, bai = bai, stringsAsFactors = FALSE)

# BAI check
no_bai_files <- input %>% filter(is.na(bai))
if (nrow(no_bai_files) > 0) {
  message("WARNING: ", nrow(no_bai_files), " BAM file(s) are missing index (.bai). Generate them with:\n",
          paste0("  samtools index ", no_bai_files$bam, collapse = "\n"),
          "\nThese samples will be excluded from the output.")
}

# Duplicate check
dup_samples <- unique(input$sample[duplicated(input$sample) | duplicated(input$sample, fromLast = TRUE)])
if (length(dup_samples) > 0) {
  dup_info <- input %>% filter(sample %in% dup_samples) %>% select(sample, bam)
  message("WARNING: ", length(dup_samples), " duplicate sample name(s) detected:\n",
          paste0("  ", dup_info$sample, " -> ", dup_info$bam, collapse = "\n"),
          "\nAction required: either exclude the duplicate samples or assign unique names before proceeding.")
}

filtered_input <- input %>% filter(!is.na(bai))
write.table(filtered_input, "input_list_complete.txt", row.names = FALSE, quote = FALSE, sep = "\t")
