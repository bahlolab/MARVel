rm(list = ls())
setwd("/vast/projects/bahlo_mtDNA/Simulation/")

lines <- readLines("./1000G/rCRS.fasta")
seq <- paste(lines[!grepl("^>", lines)], collapse = "")

# Find homopolymer runs of length >= 5 (A/C/G/T)
#m <- gregexpr("(A{5,}|C{5,}|G{5,}|T{5,})", seq, perl = TRUE)[[1]]
m <- gregexpr("(A{4,}|C{4,}|G{4,}|T{4,})", seq, perl = TRUE)[[1]]

if (length(m) == 1 && m[1] == -1) {
  runs <- data.frame()
} else {
  run_seq <- regmatches(seq, list(m))[[1]]
  start <- as.integer(m)
  end <- start + nchar(run_seq) - 1
  
  runs <- data.frame(
    start = start,
    end = end,
    base = substr(run_seq, 1, 1),
    length = nchar(run_seq),
    seq = run_seq,
    stringsAsFactors = FALSE
  )
}

# 5 repeat: 114; 4 repeat: 325
write.table(runs, "/vast/projects/bahlo_mtDNA/1000G/reference/homopolymer_tracts.txt",
            row.names = F, sep = "\t", quote = F)
