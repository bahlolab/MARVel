rm(list = ls())
library(readr)
library(dplyr)
library(fs)

# ================================
# INPUTS — edit these
# ================================
input_file  <- "input_list_complete.txt"   # output of generate_input_list.R
sh_dir      <- "sh_scripts"                # directory for per-sample shell scripts
mitohpc_dir <- "/path/to/MitoHPC/scripts" # HP_SDIR: path to MitoHPC scripts directory
output_base <- "/path/to/output"           # base directory for MitoHPC output per sample

# SBATCH settings
sbatch_time  <- "24:00:00"
sbatch_mem   <- "8G"
sbatch_nodes <- 1
sbatch_cpus  <- 8
sbatch_email <- "user@example.com"

# ================================
# Setup
# ================================
files <- read.table(input_file, header = TRUE, sep = "\t")
dir_create(sh_dir)
dir_create(file.path(sh_dir, "logs"))

combined_sh <- file.path(sh_dir, "run_all.sh")

# ================================
# Per-sample shell scripts
# ================================
for (i in seq_len(nrow(files))) {
  row         <- files[i, ]
  sample_name <- row$sample
  sample_dir  <- file.path(output_base, sample_name, "mitoHPC")
  bam_base    <- sub("\\.bam$", "", basename(row$bam))
  sh_path     <- file.path(sh_dir, paste0(sample_name, ".sh"))

  sh_content <- c(
    "#!/bin/bash",
    paste0("#SBATCH --time=",          sbatch_time),
    paste0("#SBATCH --mem=",           sbatch_mem),
    paste0("#SBATCH --nodes=",         sbatch_nodes),
    paste0("#SBATCH --cpus-per-task=", sbatch_cpus),
    paste0("#SBATCH --job-name=",      sample_name),
    paste0("#SBATCH --error=logs/",    sample_name, ".err"),
    paste0("#SBATCH --output=logs/",   sample_name, ".out"),
    paste0("#SBATCH --mail-user=",     sbatch_email),
    "",
    paste0("mkdir -p ", sample_dir),
    paste0("cd ", sample_dir),
    "",
    paste0("HP_SDIR=", mitohpc_dir),
    "cp -n $HP_SDIR/init.sh .",
    paste0("sed -i 's|^export HP_ADIR=.*|export HP_ADIR=", dirname(row$bam), "|' init.sh"),
    "sed -i 's|^export HP_L=.*|export HP_L=|' init.sh",
    ". ./init.sh",
    paste0(
      "printf '", bam_base, "\\t", row$bam, "\\tout/", bam_base, "/", bam_base, "\\n'",
      " > in.txt"
    ),
    "$HP_SDIR/run.sh > run.all.sh",
    "bash ./run.all.sh > output.log 2>&1"
  )

  writeLines(sh_content, sh_path)
  file_chmod(sh_path, "u+x")
}

# ================================
# Combined sbatch script
# ================================
combined_content <- c(
  "#!/bin/bash",
  paste0("cd ", normalizePath(sh_dir, mustWork = FALSE)),
  ""
)
for (i in seq_len(nrow(files))) {
  sh_path <- file.path(sh_dir, paste0(files$sample[i], ".sh"))
  combined_content <- c(combined_content, paste("sbatch", sh_path))
}

writeLines(combined_content, combined_sh)
file_chmod(combined_sh, "u+x")

message("Done. Submit all jobs with:\n  bash ", combined_sh)
