###################
#I/O
###################
.libPaths(c("/n/groups/walsh/indData/Maya/RLibs", .libPaths()))
library(pheatmap)
library(Seurat)
library(tidyverse)
library(RColorBrewer)

#read in comparisons
comparisons <- unname(unlist(read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/analysis_of_mo_dias_15q_data/comparisons.txt", header = F))) #In /n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/analysis_of_mo_dias_15q_data/comparisons.txt
genes_to_analyze <- c(
  "YTHDC2", "ZNF711", "LRP2", "DDX3X", "SMC1A", "USP9X",
  "TOP2B", "HNRNPU", "SRSF11", "OPHN1", "PHF3", "ZFX"
) #order defined by current 6a

###################
#CREATE SINGLE HEATMAP ROW SHOWING THE EXPRESSION OF GENES TO ANALYZE IN EXN ASD CASES VS CONTROLS
###################
#first, read in plot df used to create fig 6a so we can appropriately set our color legend
plot_df_curr_fig_6a <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Revisions/Redo_Figure_6_And_S6_To_Add_Decasien_Data/Figure_6A_S6A_Heatmap_Expression_Candidate_Contributors_FPE/fig_6a_plot_df_for_heatmap.RDS")

#now, read in exn asd cases vs control results 
##read in deg results of interest 
log2fc_mat_list <- list()
sig_mat_list <- list()
for (comp in comparisons)
{
  print('***************')
  print(comp)
  print('***************')
  deg_res_df <- readRDS(paste0("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/analysis_of_mo_dias_15q_data/Output/DEG_Results/MAST/Case_vs_Control_Male/", comp, "/Neurons_DEG_res.RDS"))
  print(deg_res_df[which(row.names(deg_res_df) == "ZFX"),])
  print(paste(
  "These genes are not expressed:",
  paste(setdiff(genes_to_analyze, row.names(deg_res_df)), collapse = ", ")
))
  deg_res_df_targeted <- deg_res_df[row.names(deg_res_df) %in% genes_to_analyze,,drop = FALSE]
  
  ##get log2fc
  log2fc_mat <- as.data.frame(t(as.data.frame(deg_res_df_targeted[,"avg_log2FC",drop = FALSE])))
  log2fc_mat_list[[comp]] <- log2fc_mat

  ##get significance 
  padj_mat <- as.data.frame(t(as.data.frame(deg_res_df_targeted[,"p_val_adj",drop = FALSE])))
  sig_mat <- as.data.frame(t(as.data.frame(apply(padj_mat, 1, function(x) ifelse(x < 0.05, '*', '')))))
  sig_mat_list[[comp]] <- sig_mat
}
full_log2fc_mat <- do.call(rbind, log2fc_mat_list)
full_sig_mat <- do.call(rbind, sig_mat_list)
full_sig_mat[abs(full_log2fc_mat) < 0.1] <- ""

#plot heatmap - note: manually confirmed that 6a has the overall global max and min
paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
global_min <- min(plot_df_curr_fig_6a, na.rm = TRUE)
#[1] -1.001275
global_max <- max(plot_df_curr_fig_6a, na.rm = TRUE)
#0.1938674
plot_limit <- max(abs(global_min), abs(global_max))
myBreaks <- c(seq(-plot_limit, 0, length.out = ceiling(paletteLength / 2) + 1),
              seq(plot_limit / paletteLength, plot_limit, length.out = floor(paletteLength / 2)))

pheatmap(
  full_log2fc_mat,
  border_color = "black",
  legend = TRUE,
  na_col = "lightgray",
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  color = myColor,
  breaks = myBreaks,
  fontsize_col = 5,
  fontsize_row = 5,
  filename = "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Revisions/Analysis_of_15q_and_Idiopathic_Autism_Data/Output/Plots/exn_all_case_vs_control_heatmap_for_fig_6a.pdf",
  display_numbers = full_sig_mat,
  width = 4,
  height = 3,
) 
