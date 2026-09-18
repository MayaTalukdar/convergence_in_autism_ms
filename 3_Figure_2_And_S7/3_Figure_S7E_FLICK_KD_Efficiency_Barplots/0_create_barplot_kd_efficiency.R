##############################
#I/O
##############################
.libPaths(c("/n/groups/walsh/indData/Maya/RLibs", .libPaths()))
library(tidyverse)
library(data.table)
library(Seurat)
library(RColorBrewer)

##############################
#PLOT EXPRESSION OF TARGET GENES IN D30 RGS
##############################
sgRNAs <- c("BAZ2B", "CLASP1", "EHMT1", "NR2F1-AS1", "PPP3CA",  "ST7", "WDFY3")
input_path <- "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas/Output/DEG_Results/day30_L1/"
RG_DEG_res_list <- setNames(lapply(sgRNAs, function(x) readRDS(paste0(input_path, x, "_KD/RG_DEG_res.RDS"))), sgRNAs)
sgRNA_kd_eff_df <- do.call(rbind, lapply(sgRNAs, function(x) RG_DEG_res_list[[x]][which(row.names(RG_DEG_res_list[[x]]) == x),]))

sgRNA_kd_eff_df$sgRNA <- row.names(sgRNA_kd_eff_df)

sgRNA_kd_eff_df$p_adj_just_7_targets <- p.adjust(sgRNA_kd_eff_df$p_val, method = "BH")

sgRNA_kd_eff_df$sig <- sapply(sgRNA_kd_eff_df$p_adj_just_7_targets, function(x) ifelse(is.na(x), "Not Applicable", ifelse(x <0.0001, "< 0.0001", ifelse(x < 0.001, "< 0.001", ifelse(x < 0.01, "< 0.01", ifelse(x < 0.05, "< 0.05", "Not Significant"))))))
sgRNA_kd_eff_df$sig <- factor(as.character(sgRNA_kd_eff_df$sig), levels = c("< 0.0001", "< 0.001", "< 0.01", "< 0.05", "Not Significant", "Not Applicable"))

sgRNA_kd_eff_df <- sgRNA_kd_eff_df %>% mutate(exp_avg_log2FC = 2^avg_log2FC) %>% arrange(desc(exp_avg_log2FC))

sgRNA_kd_eff_df$sgRNA <- factor(sgRNA_kd_eff_df$sgRNA, levels = rev(unique(sgRNA_kd_eff_df$sgRNA)))

#create plot (for supplement)
cols <- c(brewer.pal(8, "Greens")[c(8, 6, 4, 2)], "lightgray", "white")
pdf("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Revisions/Analysis_of_KD_Efficiency_In_LT_Organoids/knockdown_efficiency_barplot_using_adj_pvalue.pdf", height = 3, width = 10)
ggplot(data = sgRNA_kd_eff_df, aes(x = sgRNA, y = exp_avg_log2FC, fill = sig, width = 0.8)) + geom_bar(stat = "identity", color = "black") + 
theme_minimal() +
xlab("Gene") + 
ylab("Relative Expression of Target") + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
geom_hline(yintercept = 1, linetype = "dashed") + 
theme(axis.text=element_text(size=7), axis.title=element_text(size=8,face="bold")) + 
scale_fill_manual(values = cols) + 
guides(fill=guide_legend(title="Adj. P-value")) + 
theme(legend.text=element_text(size=7), legend.title=element_text(size=8,face="bold")) + 
guides(color = "none")
dev.off()