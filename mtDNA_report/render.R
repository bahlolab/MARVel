#!/usr/bin/env Rscript
# Run this script to render the mtDNA report.
# Always call from the mtDNA_report/ directory:
#
#   cd /path/to/mtDNA_report
#   Rscript render.R
#
# Or from R:
#   setwd("/path/to/mtDNA_report")
#   source("render.R")

# ── ensure working directory is the mtDNA_report folder ──────────────────────
script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(0)$ofile)),   # when source()'d
  error = function(e) getwd()                   # when run via Rscript
)
setwd(script_dir)

# ── load config to get report_dir and cohort_label ────────────────────────────
source("config.R")

if (!dir.exists(report_dir)) {
  dir.create(report_dir, recursive = TRUE)
  message("Created report_dir: ", report_dir)
}

# ── render ────────────────────────────────────────────────────────────────────
output_file <- file.path(
  report_dir,
  paste0("mtDNA_report_", cohort_label, "_", Sys.Date(), ".html")
)
output_file <- gsub("[[:space:]]", "_", output_file)  # no spaces in filename

message("Rendering report...")
message("  Config:  ", normalizePath("config.R"))
message("  Output:  ", output_file)

rmarkdown::render(
  input       = "report.Rmd",
  output_file = output_file,
  params      = list(
    config = normalizePath("config.R"),
    title  = paste0(cohort_label, " — mtDNA Variant Report")
  ),
  envir = new.env(parent = globalenv())
)

message("Done: ", output_file)
