###############
#I/O
###############
.libPaths(c("/n/groups/walsh/indData/Maya/RLibs", .libPaths()))
library(tidyverse)
library(DESeq2)
library(data.table)
library(pheatmap)
library(RColorBrewer)
library(stringr)
library(openxlsx)

#read in the raw conts
processed_mat_merged <- fread("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/shRNA_perturbation_data_from_plasmidsaurus/K2Y3JM_results/K2Y3JM-expression-matrix.tsv") %>% as.data.frame()
processed_mat_merged <- processed_mat_merged[-which(duplicated(processed_mat_merged$gene_name)),]
row.names(processed_mat_merged) <- processed_mat_merged$gene_name
processed_mat_merged <- processed_mat_merged[,which(grepl("_count", colnames(processed_mat_merged)))]

#read in sample key 
sample_key <- read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/shRNA_perturbation_data_from_plasmidsaurus/Sample_Key_Plasmidsaurus.txt", header = TRUE)

#correct colnames
tube_numbers <- str_extract(colnames(processed_mat_merged), "(?<=K2Y3JM_)\\d+(?=_count)")
lookup <- setNames(sample_key$Sample, sample_key$Tube)
colnames(processed_mat_merged) <- lookup[tube_numbers]

#exclude ctrl2_lacZ1to3
processed_mat_merged <- processed_mat_merged[,-which(grepl("ctrl2_lacZ1to3", colnames(processed_mat_merged)))]

# create metadata
metadata_merged <- lapply(colnames(processed_mat_merged), function(x) {
  strsplit(x, "_")[[1]][1]}) %>%
  unlist() %>%
  data.frame(condition = .)

row.names(metadata_merged) <- colnames(processed_mat_merged)

metadata_merged$condition <- sapply(metadata_merged$condition, function(x) ifelse(grepl("ctrl", x), "ctrl", x))

# #############
# #RUN DESEQ2
# #############
# #run DESeq2
# dds <- DESeqDataSetFromMatrix(
#  countData = round(processed_mat_merged),
#   colData = metadata_merged,
#   design= ~ condition)
# dds <- DESeq(dds)
# saveRDS(dds, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Revisions/Analysis_of_shRNA_Data/0_Process_And_Run_Differential_Expression/differential_expression_analysis/dds.RDS") 
dds <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Revisions/Analysis_of_shRNA_Data/0_Process_And_Run_Differential_Expression/differential_expression_analysis/dds.RDS") 

##############
#GET DIFF EXPR RESULTS
##############
#get results 
control <- "ctrl"
path <- "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Revisions/Analysis_of_shRNA_Data/0_Process_And_Run_Differential_Expression/differential_expression_analysis/DEG_Results/OriginalResultsBeforeManualAdjustment"
test_conditions <- setdiff(unique(metadata_merged$condition), "ctrl")

for (condition in test_conditions)
{
	cleaned_condition <- gsub("_", "-", condition)
    write.table(results(dds, contrast = c("condition", condition, control)), paste0(path, "/", cleaned_condition, "_vs_", control, ".txt"))
    print(paste0(condition, " completed!"))
}

###############
#GET BASE MEAN OF GENE EXPR IN SHRNA NPCS
###############
baseMean_df <- read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Revisions/Analysis_of_shRNA_Data/0_Process_And_Run_Differential_Expression/differential_expression_analysis/DEG_Results/OriginalResultsBeforeManualAdjustment/BAZ2B_vs_ctrl.txt", header = TRUE) %>% dplyr::select(baseMean)
universe_shRNA <- row.names(baseMean_df %>% filter(baseMean > 10))

##############
#MANUALLY ADJUST PVALS
##############
path <- "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Revisions/Analysis_of_shRNA_Data/0_Process_And_Run_Differential_Expression/differential_expression_analysis/DEG_Results/OriginalResultsBeforeManualAdjustment/"
DESeq_res_list <- lapply(list.files(path), function(x) read.table(paste0(path, x), header = TRUE, row.names = 1))
DESeq_res_list_filtered <- lapply(DESeq_res_list, function(x) x %>% filter(!row.names(x) %in% universe_shRNA))
DESeq_res_list_filtered <- lapply(DESeq_res_list_filtered, function(x) x %>% mutate(padj_DEseq = padj) %>% mutate(padj = NA) %>% select(baseMean, log2FoldChange, lfcSE, stat, pvalue, padj_DEseq, padj))
DESeq_res_list <- lapply(DESeq_res_list, function(x) x %>% filter(row.names(x) %in% universe_shRNA) %>% arrange(padj))
DESeq_res_list <- lapply(DESeq_res_list, function(x) x %>% mutate(padj_DEseq = padj) %>% mutate(padj = p.adjust(pvalue, method = "BH")) %>% select(baseMean, log2FoldChange, lfcSE, stat, pvalue, padj_DEseq, padj))
DESeq_res_all_list <- lapply(seq_along(DESeq_res_list), function(x) rbind(DESeq_res_list[[x]], DESeq_res_list_filtered[[x]]))
DESeq_res_all_list <- lapply(DESeq_res_all_list, function(x) x %>% arrange(padj))
names(DESeq_res_all_list) <- sapply(list.files(path), function(x) gsub(".txt", "", x))
lapply(seq_along(DESeq_res_all_list), function(x) write.table(DESeq_res_all_list[[x]], paste0("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Revisions/Analysis_of_shRNA_Data/0_Process_And_Run_Differential_Expression/differential_expression_analysis/DEG_Results/WithManualAdjustment/", names(DESeq_res_all_list)[x], "_manual_adjustment_Cutoff_Of_10.txt")))

#write out as excel file
wb <- createWorkbook()

for (i in seq_along(DESeq_res_all_list)) {
  sheet_name <- names(DESeq_res_all_list)[i]
  
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet = sheet_name, DESeq_res_all_list[[i]], rowNames = TRUE)
}
saveWorkbook(
  wb,
  file = "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Revisions/Analysis_of_shRNA_Data/0_Process_And_Run_Differential_Expression/differential_expression_analysis/DEG_Results/WithManualAdjustment/Supp_Table_shRNA_Diff_Expr_Results.xlsx",
  overwrite = TRUE
)

