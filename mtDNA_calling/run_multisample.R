rm(list = ls())
library(readr)
library(dplyr)
library(fs)

# ================================
# INPUTS
# ================================
files <- read.table(paste0("/vast/scratch/users/", Sys.getenv("USER"), "/input_list_complete.txt"), header = T, sep = "\t")
sh_dir <- paste0("/vast/scratch/users/", Sys.getenv("USER"), "/sh_scripts")  # folder to write per-row shell scripts
if(!(dir.exists(sh_dir))) dir.create(sh_dir)
combined_sh <- file.path(sh_dir, "run_all_full.sh")  # combined sbatch script
download_base <- "/path/to/download"  # base folder for downloads
if(!(dir.exists(download_base))) dir.create(download_base)

# SBATCH settings
sbatch_time <- "24:00:00"
sbatch_mem <- "8G"
sbatch_nodes <- 1
sbatch_cpus <- 8
sbatch_email <- "<email address>"

# ================================
# Create per-row shell scripts
# ================================
for (i in seq_len(nrow(files))) {
  row <- files[i, ]
  sample_name <- row$sample
  sample_dir <- file.path(download_base, sample_name)
  bam_base <- sub("\\.bam$", "", basename(row$bam))
  
  # path to per-row shell script
  sh_path <- file.path(sh_dir, paste0(sample_name, ".sh"))
  
  # start shell script content
  sh_content <- c(
    "#!/bin/bash",
    paste0("#SBATCH --time=", sbatch_time),
    paste0("#SBATCH --mem=", sbatch_mem),
    paste0("#SBATCH --nodes=", sbatch_nodes),
    paste0("#SBATCH --cpus-per-task=", sbatch_cpus),
    paste0("#SBATCH --job-name=", sample_name),
    paste0("#SBATCH --error=", sample_name, ".err"),
    paste0("#SBATCH --output=", sample_name, ".out"),
    paste0("#SBATCH --mail-user=", sbatch_email),
    "",
    paste0("mkdir -p ", sample_dir)
  )
  
  # run mitohpc
  sh_content <- c(sh_content,
                  "",
                  "####################### Run mitoHPC #######################",
                  paste0("cd ",sample_dir),
                  "mkdir -p mitoHPC",
                  "cd mitoHPC",
                  "HP_SDIR=/stornext/Bioinf/data/lab_bahlo/software/apps/MitoHPC/scripts",
                  "cp -i $HP_SDIR/init.sh .",
                  paste0("sed -i 's|^export HP_ADIR=.*|export HP_ADIR=",dirname(row$bam),"|' init.sh"),
                  "sed -i 's|^export HP_L=.*|export HP_L=|' init.sh",
                  #"sed -i 's|^export HP_M=.*|export HP_M=mutect2|' init.sh", # edit if want to use different caller
                  ". ./init.sh",
                  # write in.txt
                  paste0(
                    "echo -e '",
                    bam_base, "\\t",
                    row$bam, "\\t",
                    "out/", bam_base, "/", bam_base,
                    "' > in.txt"
                  ),
                  "$HP_SDIR/run.sh > run.all.sh",
                  "bash ./run.all.sh > output.log 2>&1"
  )
  # write the shell script
  writeLines(sh_content, sh_path)
  system(paste("chmod +x", sh_path))
}

# ================================
# Create combined sbatch run script
# ================================
combined_content <- c(
  "#!/bin/bash",
  paste0("cd ",sh_dir,"/logs"),
  ""
)

# add sbatch commands for each per-row script
for (i in seq_len(nrow(files))) {
  sample_name <- files$sample[i]
  sh_path <- file.path(sh_dir, paste0(sample_name, ".sh"))
  combined_content <- c(combined_content, paste("sbatch", sh_path))
}

# write combined script
writeLines(combined_content, combined_sh)
system(paste("chmod +x", combined_sh))

#cat("Done! Run the combined script with:\n")
#cat(paste0("sbatch ", combined_sh, "\n"))
