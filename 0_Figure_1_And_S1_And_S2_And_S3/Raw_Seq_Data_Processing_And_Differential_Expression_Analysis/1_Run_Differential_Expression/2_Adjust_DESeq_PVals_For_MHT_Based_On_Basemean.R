library(tidyverse)

#taken from the NTC 1 vs NTC2/3 comparisons 
baseMean_in_NTCs <- read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/baseMeanInNTCs.csv", row.names = 1)
genes_to_keep <- baseMean_in_NTCs %>% filter(baseMean >= 10) %>% rownames()

path <- "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/DEG_Results/OriginalResultsBeforeManualAdjustment/"
DESeq_res_list <- lapply(list.files(path), function(x) read.table(paste0(path, x), header = TRUE, row.names = 1))
DESeq_res_list_filtered <- lapply(DESeq_res_list, function(x) x %>% filter(!row.names(x) %in% genes_to_keep))
DESeq_res_list_filtered <- lapply(DESeq_res_list_filtered, function(x) x %>% mutate(padj_DEseq = padj) %>% mutate(padj = NA) %>% select(baseMean, log2FoldChange, lfcSE, stat, pvalue, padj_DEseq, padj))
DESeq_res_list <- lapply(DESeq_res_list, function(x) x %>% filter(row.names(x) %in% genes_to_keep) %>% arrange(padj))
DESeq_res_list <- lapply(DESeq_res_list, function(x) x %>% mutate(padj_DEseq = padj) %>% mutate(padj = p.adjust(pvalue, method = "BH")) %>% select(baseMean, log2FoldChange, lfcSE, stat, pvalue, padj_DEseq, padj))
DESeq_res_all_list <- lapply(seq_along(DESeq_res_list), function(x) rbind(DESeq_res_list[[x]], DESeq_res_list_filtered[[x]]))
DESeq_res_all_list <- lapply(DESeq_res_all_list, function(x) x %>% arrange(padj))
names(DESeq_res_all_list) <- sapply(list.files(path), function(x) gsub(".txt", "", x))

lapply(seq_along(DESeq_res_all_list), function(x) write.table(DESeq_res_all_list[[x]], paste0("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/DEG_Results/WithManualAdjustment/", names(DESeq_res_all_list)[x], "_manual_adjustment_Cutoff_Of_10.txt")))
