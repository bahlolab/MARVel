rm(list = ls())

# ================================
# INPUTS — edit these
# ================================
input_file  <- "input_list_complete.txt"   # output of generate_input_list.R
sh_dir      <- "sh_scripts"                # directory containing per-sample .sh scripts
output_base <- "/path/to/output"           # base directory for mitoHPC output (same as run_multisample.R)

# ================================
# Check jobs
# ================================
files <- read.table(input_file, header = TRUE, sep = "\t")

# Derive bam_base from the bam column (same logic as run_multisample.R)
files$bam_base <- sub("\\.bam$", "", basename(files$bam))

# Expected mitoHPC output files per sample (relative to <output_base>/<sample>/mitoHPC/).
# Uses bam_base (BAM filename without .bam extension) — no hardcoded suffixes.
# A job is considered complete if ALL of these exist.
expected_files_for <- function(bam_base) {
  c(
    sprintf("out/%s/%s.mutect2.cvg.stat",  bam_base, bam_base),
    sprintf("out/%s/%s.mutect2.00.vcf",    bam_base, bam_base),
    sprintf("out/%s/%s.mutect2.haplogroup",bam_base, bam_base),
    sprintf("out/%s/%s.mutect2.haplocheck",bam_base, bam_base)
  )
}

check_sample <- function(sample, bam_base) {
  base_dir <- file.path(output_base, sample, "mitoHPC")
  paths    <- expected_files_for(bam_base)
  full     <- file.path(base_dir, paths)
  exists   <- file.exists(full)
  list(complete = all(exists), missing = paths[!exists])
}

results  <- mapply(check_sample, files$sample, files$bam_base, SIMPLIFY = FALSE)
complete <- sapply(results, `[[`, "complete")

n_total    <- nrow(files)
n_complete <- sum(complete)
n_failed   <- sum(!complete)

message(sprintf("Job check: %d / %d complete, %d incomplete.", n_complete, n_total, n_failed))

if (n_failed > 0) {
  failed_samples <- files$sample[!complete]

  message(sprintf("\n%d sample(s) with missing output:", n_failed))
  for (s in failed_samples) {
    miss <- results[[which(files$sample == s)]]$missing
    message(sprintf("  %s — missing: %s", s, paste(basename(miss), collapse = ", ")))
  }

  sh_paths <- file.path(sh_dir, paste0(failed_samples, ".sh"))
  existing <- file.exists(sh_paths)

  rerun_script <- "rerun_failed.sh"
  lines <- c("#!/bin/bash", "")
  for (i in seq_along(failed_samples)) {
    if (existing[i]) {
      lines <- c(lines, paste("sbatch", sh_paths[i]))
    } else {
      lines <- c(lines, paste0("# No script found for: ", failed_samples[i]))
    }
  }
  writeLines(lines, rerun_script)
  Sys.chmod(rerun_script, mode = "0755")

  message(sprintf(
    "\nRe-run script written to: %s\nSubmit with:\n  bash %s",
    rerun_script, rerun_script
  ))
} else {
  message("All jobs completed successfully.")
}
