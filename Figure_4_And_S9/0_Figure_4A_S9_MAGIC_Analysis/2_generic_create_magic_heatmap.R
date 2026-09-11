#import libraries
library(tidyverse)
library(readxl)
library(pheatmap)
library(RColorBrewer)
library(data.table)
library(viridis) 
library(data.table)
library(fgsea)

generate_overlap_p_val <- function(list1, list2, universe)
{
    list1 <- intersect(list1, universe)
    list2 <- intersect(list2, universe)

    inList1AndList2 <- length(intersect(list1, list2))
    inList1AndNotList2 <- length(setdiff(list1, list2))
    inList2AndNotList1 <- length(setdiff(list2, list1))
    inNeither <- length(setdiff(universe, c(list1, list2)))

    mat <- matrix(c(inList1AndList2, inList1AndNotList2, inList2AndNotList1, inNeither), nrow = 2)
    pval <- fisher.test(mat, alternative = "greater")$p.value

    return(pval)
}

######################
#JUST TOP REGULATORS
#######################
#read in MAGIC results
targets <- unname(unlist(read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/final_targets_to_carry_forward.csv")))
MAGIC_res_list <- lapply(targets, function(x) read_excel(paste0("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/magic_analysis/Output/", x, "/", x, "_Summary.xls")))
names(MAGIC_res_list) <- targets

#subset to significant regulators
MAGIC_res_list_sig <- lapply(MAGIC_res_list, function(x) x[which(x[,"Corrected p"] < 0.05),])
MAGIC_res_list_sig_regulators <- lapply(MAGIC_res_list_sig, function(x) x$Factor)
MAGIC_res_list_sig_regulators_summary <- as.data.frame(table(do.call(c, MAGIC_res_list_sig_regulators))) %>% arrange(desc(Freq))

#regulators for figure
top_regs <- table(do.call(c, MAGIC_res_list_sig_regulators))
top_regs <- rev(sort(top_regs))
top_regs <- top_regs[grepl("\\w$", names(top_regs))]
cutoff <- 16
targets <- names(top_regs[top_regs >= cutoff])
regs_for_figures <- targets

#create df for heatmap
MAGIC_heatmap_df <- lapply(regs_for_figures, function(x) {
  sapply(MAGIC_res_list, function(y) {
    p_val <- y %>% filter(Factor == x) %>% pull(`Corrected p`)
    if (length(p_val) == 0) 1 else p_val
  })
}) %>% as.data.frame()

#add annotation row for whether a regulator is itself implicated in ASD 
sfari <- read.csv("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/sfari.csv", header = TRUE)
genes <- unique(sfari$gene.symbol)
annotation_row_df <- as.data.frame(sapply(regs_for_figures, function(x) ifelse(x %in% genes, "Yes", "No")))

#sfari analysis 
##read in sfari genes 
sfari <- read.csv("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/sfari.csv")
sfari <- sfari %>% filter(gene.score == 1)
sfari_og <- sfari$gene.symbol

##perturbation analysis
baseMean_in_NTCs <- read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/baseMeanInNTCs.csv", row.names = 1)
genes_to_keep <- row.names(baseMean_in_NTCs %>% filter(baseMean >= 10))
processed_mat_merged <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/processed_mat_merged.RDS")
processed_mat_merged <- processed_mat_merged[genes_to_keep,]
universe <- genes_to_keep
sfari <- sfari_og[sfari_og %in% universe]

##rank of regulator cor to sfari genes relative to all other genes 
corrmat <- cor(t(processed_mat_merged))
corrmat_sfari <- corrmat[,sfari]
sfari_corr <- apply(corrmat_sfari, 1, median)
sfari_corr <- rev(sfari_corr[order(sfari_corr)])
sfari_corr_percentile <- setNames(ecdf(sfari_corr)(sfari_corr), names(sfari_corr))
sfari_corr_annot_perturbation <- sapply(regs_for_figures, function(x) ifelse(x %in% names(sfari_corr_percentile), sfari_corr_percentile[x], NA))

##rank of regulator cor to sfari genes vs all other gene sets 
p_val_vec <- c()
for (reg in regs_for_figures)
{
    print(reg)
    if (reg %in% row.names(processed_mat_merged))
    {
    	true_med <- median(corrmat[reg, colnames(corrmat) %in% sfari])
    	res <- unlist(lapply(seq_len(1000), function(x) median(corrmat[reg, colnames(corrmat) %in% sample(genes_to_keep, length(sfari))])))
    	p_val <- 1 - ecdf(res)(true_med)
   		p_val_vec <- c(p_val_vec, p_val)
    } else
    {
    	p_val_vec <- c(p_val_vec, NA)
    }
}
names(p_val_vec) <- regs_for_figures
reg_corr_annot_perturbation <- p_val_vec
reg_corr_annot_perturbation_sig <- sapply(reg_corr_annot_perturbation, function(x) ifelse(x < 0.05, "Yes", "No"))

annotation_row_df <- cbind(annotation_row_df, reg_corr_annot_perturbation_sig)
colnames(annotation_row_df) <- c("Mutations Associated with ASD?", "Sig. Co-expression with SFARI Genes? (Perturbations)")
annotation_colors_diff_expr <- c("Yes" = "black" , "No" = "white")
annot_colors <- list("Mutations Associated with ASD?" = annotation_colors_diff_expr, "Sig. Co-expression with SFARI Genes? (Perturbations)" = annotation_colors_diff_expr)

colnames(MAGIC_heatmap_df) <- regs_for_figures
saveRDS(MAGIC_heatmap_df, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/magic_analysis/MAGIC_heatmap_df_all_sig_regs.RDS")
saveRDS(annotation_row_df, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/magic_analysis/annotation_row_df_all_sig_regs.RDS")
breaks <- c(0, 0.0001, 0.001, 0.01, 0.05, 1)
colors <- c(brewer.pal(8, "Greens")[c(8, 6, 4, 2)], "lightgray")
pheatmap(MAGIC_heatmap_df, color = colors, breaks = breaks, border_color = "black", legend = FALSE, annotation_legend = TRUE, annotation_col = annotation_row_df,  annotation_colors = annot_colors, filename =  "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Figure_4/heatmap_of_MAGIC_top_regulators.pdf", width = 14, height = 6)

