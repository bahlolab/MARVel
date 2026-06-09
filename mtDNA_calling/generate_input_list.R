rm(list = ls())
library(dplyr)

# ================================
# INPUTS — edit these
# ================================
bam_dirs <- c("/path/to/bam/dir1", "/path/to/bam/dir2")
bam_pattern <- "\\.bam$"          # regex to match BAM files
sample_strip <- "\\.merged\\.bam$" # suffix to strip for sample name
outfile <- "input_list_complete.txt"

# ================================
# Build input table
# ================================
bam <- unlist(lapply(bam_dirs, list.files, pattern = bam_pattern, full.names = TRUE))

if (length(bam) == 0) stop("No BAM files found in the specified directories.")

sample <- sub(sample_strip, "", basename(bam))
bai <- paste0(bam, ".bai")
bai[!file.exists(bai)] <- NA

input <- data.frame(sample = sample, bam = bam, bai = bai, stringsAsFactors = FALSE)

# BAI check
no_bai <- input |> filter(is.na(bai))
if (nrow(no_bai) > 0) {
  message("WARNING: ", nrow(no_bai), " BAM file(s) missing index (.bai). Generate with:\n",
          paste0("  samtools index ", no_bai$bam, collapse = "\n"),
          "\nThese samples will be excluded.")
}

# Duplicate sample name check
dup_samples <- unique(input$sample[duplicated(input$sample) | duplicated(input$sample, fromLast = TRUE)])
if (length(dup_samples) > 0) {
  dup_info <- input |> filter(sample %in% dup_samples) |> select(sample, bam)
  message("WARNING: ", length(dup_samples), " duplicate sample name(s):\n",
          paste0("  ", dup_info$sample, " -> ", dup_info$bam, collapse = "\n"),
          "\nResolve before proceeding.")
}

filtered_input <- input |> filter(!is.na(bai))
message(nrow(filtered_input), " samples written to ", outfile)
write.table(filtered_input, outfile, row.names = FALSE, quote = FALSE, sep = "\t")
