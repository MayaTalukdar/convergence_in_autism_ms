####################
#I/O
###################
.libPaths(c("/n/groups/walsh/indData/Maya/RLibs", .libPaths()))
library(tidyverse)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(ggplot2)

#read in SFARI genes
sfari <- read.csv("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/sfari.csv", header = TRUE)
SFARI_genes <- unique(sfari$gene.symbol)
all_SFARI_genes <- SFARI_genes

#read in differentially expressed genes across all perturbations
genes <- unname(unlist(read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/final_targets_to_carry_forward.csv", header = FALSE)))
genes <- unname(sapply(genes, function(x) gsub("_", "-", x)))
path <- "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/DEG_Results/WithManualAdjustment/"
all_DESeq_res_list <- lapply(genes, function(x) read.table(paste0(path, x, "_vs_NTC_withDedup_manual_adjustment_Cutoff_Of_10.txt"), header = TRUE))
names(all_DESeq_res_list) <- genes
all_DESeq_res_list <- lapply(all_DESeq_res_list, function(x) x %>% filter(!is.na(padj)))
crispr_universe <- row.names(all_DESeq_res_list[[1]])
all_DEG_res_list <- lapply(all_DESeq_res_list, function(x) row.names((x %>% filter(padj < 0.05) %>% filter(log2FoldChange < 0))))

path <- "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Revisions/Analysis_of_shRNA_Data/0_Process_And_Run_Differential_Expression/differential_expression_analysis/DEG_Results/WithManualAdjustment/"
all_DESeq_res_list_shRNA <- lapply(list.files(path, pattern = ".txt"), function(x) read.table(paste0(path, x), header = TRUE))
names(all_DESeq_res_list_shRNA) <- sapply(list.files(path, pattern = ".txt"), function(x) gsub("_vs_ctrl_manual_adjustment_Cutoff_Of_10.txt", "", x))
all_DESeq_res_list_shRNA <- lapply(all_DESeq_res_list_shRNA, function(x) x %>% filter(!is.na(padj)))
shrna_universe <- row.names(all_DESeq_res_list_shRNA[[1]])
all_DEG_res_list_shrna <- lapply(all_DESeq_res_list_shRNA, function(x) row.names((x %>% filter(padj < 0.05) %>% filter(log2FoldChange < 0) )))

#get common perturbations
common_perturbations <- intersect(names(all_DESeq_res_list), names(all_DESeq_res_list_shRNA))
universe <- intersect(shrna_universe, crispr_universe)

generate_overlap_p_val_and_or <- function(list1, list2, universe)
{
      list1 <- unique(toupper(list1))
      list2 <- unique(toupper(list2))
      universe <- unique(toupper(universe))
	list1 <- intersect(list1, universe)
	list2 <- intersect(list2, universe)


	inList1AndList2 <- length(intersect(list1, list2))
	inList1AndNotList2 <- length(setdiff(list1, list2))
	inList2AndNotList1 <- length(setdiff(list2, list1))
	inNeither <- length(setdiff(universe, c(list1, list2)))


	mat <- matrix(c(inList1AndList2, inList1AndNotList2, inList2AndNotList1, inNeither), nrow = 2)
	pval <- fisher.test(mat, alternative = "greater")$p.value
	estimate <- unname(fisher.test(mat, alternative = "greater")$estimate)


	return(c(pval, estimate, fisher.test(mat, alternative = "greater")$conf.int[1]))
}

################
#RUN OVERLAP 
################
plot_df <- as.data.frame(setNames(lapply(common_perturbations, function(x) generate_overlap_p_val_and_or(all_DEG_res_list[[x]], all_DEG_res_list_shrna[[x]], universe)), common_perturbations)) %>% t() %>% as.data.frame()
colnames(plot_df) <- c("pval", "OR", "lower.ci")
plot_df$gene <- row.names(plot_df)
plot_df$p_adj <- p.adjust(plot_df$pval, method = "BH")
plot_df <- plot_df %>% arrange(desc(OR))
plot_df$gene <- factor(plot_df$gene, levels = unique(plot_df$gene))
plot_df$sig <- sapply(plot_df$p_adj, function(x) ifelse(x <0.0001, "< 0.0001", ifelse(x < 0.001, "< 0.001", ifelse(x < 0.01, "< 0.01", ifelse(x < 0.05, "< 0.05", "Not Significant")))))
plot_df$sig <- factor(as.character(plot_df$sig), levels = c("< 0.0001", "< 0.001", "< 0.01", "< 0.05", "Not Significant"))

#plot
cols <- c(brewer.pal(8, "Greens")[c(8, 6)], "lightgray")
pdf("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Actually_Submitted/For_Rebuttal/Output/Figures/shRNA_Cas13_KD_Downreg_DEGs_Overlap/downreg_deg_overlap_barplot.pdf", height = 3, width = 5)
sfari_plot <- (ggplot(data = plot_df, aes(x = OR, y = forcats::fct_rev(gene), fill = sig, width = 0.8, color = "black")) + geom_bar(stat = "identity", color = "black") + 
theme_minimal() +
ylab("Target") + 
xlab("Odds Ratio") +     
geom_errorbar(aes(xmin = lower.ci, xmax = OR),  # only go down to lower CI
                  width = 0,                      # horizontal cap width
                  color = "black") +
theme(axis.text=element_text(size=7), axis.title=element_text(size=8,face="bold")) + 
scale_fill_manual(values = cols) + 
guides(fill=guide_legend(title="Adjusted p-value")) + 
theme(legend.text=element_text(size=7), legend.title=element_text(size=8,face="bold")) + 
guides(color = "none"))
print(sfari_plot)
dev.off()
