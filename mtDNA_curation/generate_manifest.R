rm(list = ls())

# ================================
# INPUTS — edit these
# ================================
results_dir       <- "/path/to/mitoHPC/results"  # parent directory of per-sample folders
out_file          <- "manifest.tsv"               # output path
sample_id_pattern <- NULL                         # optional regex to filter sample dirs
                                                  # e.g. "^AA" — NULL = use all subdirs

# ================================
# Discover samples
# ================================
dirs <- list.dirs(results_dir, recursive = FALSE, full.names = FALSE)
if (!is.null(sample_id_pattern)) dirs <- dirs[grepl(sample_id_pattern, dirs)]
dirs <- sort(dirs)
message(sprintf("Found %d sample director%s.", length(dirs), ifelse(length(dirs) == 1, "y", "ies")))

rows <- lapply(dirs, function(id) {
  sample_folder <- file.path(results_dir, id, "mitoHPC")

  vcf_hits <- list.files(file.path(sample_folder, "out"),
                         pattern   = "\\.mutect2\\.00\\.vcf$",
                         recursive = TRUE,
                         full.names = FALSE)
  if (length(vcf_hits) == 0) {
    message(sprintf("  SKIP %s — VCF not found under %s/out/", id, sample_folder))
    return(NULL)
  }

  data.frame(SampleID     = id,
             FamilyID     = NA_character_,
             SampleFolder = sample_folder,
             stringsAsFactors = FALSE)
})

manifest <- do.call(rbind, Filter(Negate(is.null), rows))
message(sprintf("%d / %d samples included.", nrow(manifest), length(dirs)))

write.table(manifest, out_file, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
message("Written: ", out_file)
