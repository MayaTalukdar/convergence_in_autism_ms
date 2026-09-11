library(tidyverse)

file_path <- "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/DEG_Results/WithManualAdjustment/"
genes <- unname(unlist(read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/final_targets_to_carry_forward.csv", header = FALSE)))
sgRNAs <- unname(sapply(genes, function(x) gsub("_", "-", x)))
all_DESeq_res <- lapply(sgRNAs, function(x) read.table(paste0(file_path, x, "_vs_NTC_withDedup_manual_adjustment_Cutoff_Of_10.txt")))
names(all_DESeq_res) <- sgRNAs
gene_lists <- lapply(all_DESeq_res, function(x) row.names(x %>% filter(padj < 0.05)))
basemean_df <-  read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/baseMeanInNTCs.csv", row.names = 1)
universe <- row.names(basemean_df %>% filter(baseMean > 10))
gene_lists[["All"]] <- universe

#convert to df
max_length <- max(lengths(gene_lists))
padded_list <- lapply(gene_lists, function(x) {
  if (length(x) < max_length) {
    c(x, rep("", max_length - length(x)))
  } else {
    x
  }
})
df <- as.data.frame(t(as.data.frame(matrix(unlist(padded_list), nrow=length(gene_lists), byrow=TRUE))))
colnames(df) <- names(gene_lists)
row.names(df) <- NULL
df <- df %>% relocate(All)

write.table(df, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/magic_analysis/input_for_MAGIC.csv", sep = ",", row.names = F, col.names = T, quote = F)
