##############################
#I/O
##############################
library(tidyverse)
library(RColorBrewer)
library(pheatmap)

num_perturb_cutoff <- as.numeric(commandArgs(trailingOnly=TRUE)[1])

#read in differentially expressed genes across all perturbations
#get list of genes to use
genes <- unname(unlist(read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/final_targets_to_carry_forward.csv", header = FALSE))) #In Input_Files_Not_Generated_By_Scripts
genes <- unname(sapply(genes, function(x) gsub("_", "-", x)))
path <- "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/DEG_Results/WithManualAdjustment/" 
all_DESeq_res_list <- lapply(genes, function(x) read.table(paste0(path, x, "_vs_NTC_withDedup_manual_adjustment_Cutoff_Of_10.txt"), header = TRUE))
names(all_DESeq_res_list) <- genes
DEG_res_list <- setNames(lapply(all_DESeq_res_list, function(x) row.names(x %>% filter(padj < 0.05))), genes)
names(DEG_res_list) <- genes

#read in high confidence SFARI genes 
sfari <- read.csv("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/sfari.csv", header = TRUE) #In Input_Files_Not_Generated_By_Scripts
all_hc_SFARI_genes <- sfari$gene.symbol[which(sfari$gene.score == 1)]

##############################
#PLOT SETUP 
##############################
#determine the set of SFARI genes that are DEGs in at least one perturbation 
hc_sfari_count <- sapply(all_hc_SFARI_genes, function(x) length(which(sapply(DEG_res_list, function(y) x %in% y))))
all_hc_SFARI_genes_DEGs <- names(which(hc_sfari_count >= num_perturb_cutoff))

#get all relevant information for these genes that needs to be plotted (ex: log fold change, number of perturbations in which DEG, etc.)
SFARI_DEG_res_list <- setNames(lapply(names(all_DESeq_res_list), function(x) all_DESeq_res_list[[x]] %>% filter(row.names(all_DESeq_res_list[[x]]) %in% all_hc_SFARI_genes_DEGs) %>% mutate(Target = x) %>% rownames_to_column(var = "Gene") %>% arrange(Gene) %>% mutate(log2FC_to_plot = log2FoldChange)), genes)

#extract p-values for each of the DEG-perturbation combinations, as we will use these to label a heatmap box with '*' if < 0.05
label_df <- lapply(SFARI_DEG_res_list, function(x) x$padj)
label_df <- as.data.frame(do.call(cbind, label_df))
colnames(label_df) <- genes
row.names(label_df) <- all_hc_SFARI_genes_DEGs
label_df <- as.data.frame(t(as.data.frame(label_df)))
label_df <- apply(label_df, 2, function(x) ifelse(x < 0.05, "*", ""))

#create plot df
SFARI_DEG_res_list <- lapply(SFARI_DEG_res_list, function(x) x$log2FC_to_plot)
SFARI_DEG_res_df <- as.data.frame(do.call(cbind, SFARI_DEG_res_list))
colnames(SFARI_DEG_res_df) <- genes
row.names(SFARI_DEG_res_df) <- all_hc_SFARI_genes_DEGs
SFARI_DEG_res_df <- as.data.frame(t(as.data.frame(SFARI_DEG_res_df)))
plot_df <- SFARI_DEG_res_df

##############################
#CREATE PLOT 
##############################
#plot  
paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
plot_min <- min(plot_df, na.rm = TRUE)
plot_max <- max(plot_df, na.rm = TRUE)
plot_limit <- max(abs(plot_min), abs(plot_max))
myBreaks <- c(seq(-plot_limit, 0, length.out = ceiling(paletteLength / 2) + 1),
              seq(plot_limit / paletteLength, plot_limit, length.out = floor(paletteLength / 2)))

output_path <- paste0("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Figure_1/minNumberPerturbs=", num_perturb_cutoff, "/")
dir.create(output_path)
pheatmap(SFARI_DEG_res_df, border_color = "black", legend = TRUE, na_col = "lightgray", cluster_cols = TRUE, cluster_rows = TRUE, color=myColor, breaks = myBreaks, fontsize_col = 5, fontsize_row = 5, filename = paste0(output_path, "num_perturb_heatmap_including_perturbed_genes_with_dendrogram.pdf"), display_numbers = label_df, width = 8, height = 3)
outs <- pheatmap(SFARI_DEG_res_df, border_color = "black", legend = TRUE, na_col = "lightgray", cluster_cols = TRUE, cluster_rows = TRUE, color=myColor, breaks = myBreaks, fontsize_col = 5, fontsize_row = 5, display_numbers = label_df, width = 10, height = 3)





