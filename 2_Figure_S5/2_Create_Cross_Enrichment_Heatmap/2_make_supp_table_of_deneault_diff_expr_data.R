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

wb <- createWorkbook()

###############
#WRITE OUT IPSC DATA
###############
wd <- "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_denault_perturbation_data/differential_expression_analysis/deg_results/ipscs/WithManualAdjustment/"
all_res_list <- lapply(list.files(wd), function(x) read.table(paste0(wd, x), header = TRUE))
genes <- unname(sapply(list.files(wd), function(x) gsub("_vs_control_DEG_res_manual_adjustment_Cutoff_Of_10.txt", "", x)))
names(all_res_list) <- genes

for (gene in genes) 
{
  sheet_name <- paste0("IPSC_", gene)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet = sheet_name, all_res_list[[gene]], rowNames = TRUE)
}

###############
#WRITE OUT NEURON DATA
###############
wd <- "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_denault_perturbation_data/differential_expression_analysis/deg_results/neurons/WithManualAdjustment/"
all_res_list <- lapply(list.files(wd), function(x) read.table(paste0(wd, x), header = TRUE))
genes <- unname(sapply(list.files(wd), function(x) gsub("_vs_control_DEG_res_manual_adjustment_Cutoff_Of_10.txt", "", x)))
names(all_res_list) <- genes

for (gene in genes) 
{
  sheet_name <- paste0("Neuron_", gene)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet = sheet_name, all_res_list[[gene]], rowNames = TRUE)
}

saveWorkbook(
  wb,
  file = "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Revisions/Analysis_of_Deneault_Data/Output/Supp_Table_Deneault_Diff_Expr_Results.xlsx",
  overwrite = TRUE
)


