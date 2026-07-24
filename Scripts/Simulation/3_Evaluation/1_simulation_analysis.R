######## Part A - Table of Results ##########
rm(list = ls())

library(dplyr)
library(readr)
library(purrr)
library(stringr)

base_path <- "/vast/projects/bahlo_mtDNA/Simulation/1000G"

stats_files <- list.files(
  path = base_path,
  pattern = "stats.tsv",
  recursive = TRUE,
  full.names = TRUE
)
remove_paths <- c(
  "/vast/projects/bahlo_mtDNA/Simulation/1000G/HG00428/stats.tsv",
  "/vast/projects/bahlo_mtDNA/Simulation/1000G/HG00739/stats.tsv"
)
files <- setdiff(stats_files, remove_paths)

read_stats <- function(file) {
  df <- read_tsv(file, show_col_types = FALSE)
  
  df %>%
    mutate(
      F1_used = if_else(VAF == 1, F1, F1_thred),
      precision = if_else(VAF == 1, precision, precision_thred),
      FP = if_else(VAF == 1, FP, FP_thred),
      Tool = case_when(
        tool == "server2" ~ "server2",
        tool == "mity" ~ "mity",
        tool == "mitoHPC" ~ "MitoHPC",
        tool == "mtSwirl" ~ "mtSwirl",
        TRUE ~ tool
      )
    ) %>%
    select(
      VAF,
      Tool,
      TPSNV = snv_TP,
      TPINDEL = indel_TP,
      TP,
      FP,
      Sens = sensitivity,
      Prec = precision,
      F1 = F1_used,
      snv_sensitivity,
      indel_sensitivity
    ) %>% 
    mutate(
      Tool = factor(Tool, levels = c("MitoHPC","mtSwirl","mity","server2"))
    )
}
all_stats <- map_dfr(files, read_stats)

# Summarise (mean + sd)
final_summary <- all_stats %>%
  group_by(VAF, Tool) %>%
  summarise(
    TPSNV_mean = mean(TPSNV, na.rm = TRUE),
    TPSNV_sd   = sd(TPSNV, na.rm = TRUE),
    
    TPINDEL_mean = mean(TPINDEL, na.rm = TRUE),
    TPINDEL_sd   = sd(TPINDEL, na.rm = TRUE),
    
    TP_mean = mean(TP, na.rm = TRUE),
    TP_sd   = sd(TP, na.rm = TRUE),
    
    FP_mean = mean(FP, na.rm = TRUE),
    FP_sd   = sd(FP, na.rm = TRUE),
    
    Sens_mean = mean(Sens, na.rm = TRUE),
    Sens_sd   = sd(Sens, na.rm = TRUE),
    
    Prec_mean = mean(Prec, na.rm = TRUE),
    Prec_sd   = sd(Prec, na.rm = TRUE),
    
    F1_mean = mean(F1, na.rm = TRUE),
    F1_sd   = sd(F1, na.rm = TRUE),
    
    snv_sensitivity_mean = mean(snv_sensitivity,na.rm = TRUE),
    snv_sensitivity_sd = sd(snv_sensitivity,na.rm = TRUE),
    indel_sensitivity_mean = mean(indel_sensitivity,na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  arrange(VAF, Tool)

final_summary <- final_summary %>%
  mutate(
    Sens_mean = Sens_mean * 100,
    Prec_mean = Prec_mean * 100,
    Sens_sd = Sens_sd* 100,
    Prec_sd = Prec_sd* 100
  )

final_summary_output <- final_summary %>% 
  mutate(
    across(where(is.numeric), ~ round(.x, 2))
  ) %>% 
  dplyr::select(VAF,Tool, TPSNV_mean,TPINDEL_mean, TP_mean, FP_mean, Sens_mean, Prec_mean,F1_mean,snv_sensitivity_mean,TPSNV_sd, TPINDEL_sd,TP_sd,FP_sd,Sens_sd,Prec_sd,F1_sd,snv_sensitivity_sd)

write.table(
  final_summary_output,
  "/vast/projects/bahlo_mtDNA/Analysis/simulation_final_results.tsv",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE)


######## Part B - Statistics ########
# Pairwise significance between tools
stat_test <- all_stats %>%
  group_by(VAF) %>%
  pairwise_wilcox_test(
    F1 ~ Tool,
    paired = TRUE,
    p.adjust.method = "BH"
  )

effectsize <- all_stats %>% group_by(VAF) %>% rstatix::wilcox_effsize(F1 ~ Tool, paired = TRUE)

combined_results <- stat_test %>%
  left_join(
    effectsize,
    by = c("VAF", "group1", "group2")
  ) %>% 
  select("VAF", "group1", "group2","statistic","p","p.adj","p.adj.signif","effsize","magnitude")

write.table(
  combined_results,
  "/vast/projects/bahlo_mtDNA/Analysis/simulation/pairwise_results.tsv",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE)

