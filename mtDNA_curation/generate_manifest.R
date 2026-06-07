# Generate manifest.tsv for mtDNA curation Shiny app
# Scans RESULTS_DIR for AA* folders and checks that the expected VCF exists
rm(list = ls())
RESULTS_DIR <- "/vast/projects/bahlo_mtDNA/Results_ataxia"
OUT_FILE    <- file.path(RESULTS_DIR, "manifest.tsv")

dirs <- list.dirs(RESULTS_DIR, recursive = FALSE, full.names = FALSE)
aa_dirs <- sort(dirs[startsWith(dirs, "AA")])
message(sprintf("Found %d AA* directories.", length(aa_dirs)))

rows <- lapply(aa_dirs, function(id) {
  sample_folder <- file.path(RESULTS_DIR, id, "mitoHPC")
  
  # Check VCF exists
  vcf_hits <- list.files(file.path(sample_folder, "out"),
                         pattern    = "\\.mutect2\\.mutect2\\.00\\.vcf$",
                         recursive  = TRUE,
                         full.names = FALSE)
  if (length(vcf_hits) == 0) {
    message(sprintf("  SKIP %s — VCF not found", id))
    return(NULL)
  }
  
  data.frame(SampleID     = id,
             FamilyID     = NA_character_,
             SampleFolder = sample_folder,
             stringsAsFactors = FALSE)
})

manifest <- do.call(rbind, rows)
message(sprintf("%d / %d samples included.", nrow(manifest), length(aa_dirs)))

write.table(manifest, OUT_FILE, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
message("Written: ", OUT_FILE)
