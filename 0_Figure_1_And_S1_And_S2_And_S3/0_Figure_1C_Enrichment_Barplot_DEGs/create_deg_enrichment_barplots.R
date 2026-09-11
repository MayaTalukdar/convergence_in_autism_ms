###################
#I/O
###################
library(tidyverse)
library(RColorBrewer)
library(patchwork)library(dplyr)
library(patchwork)


#read in differentially expressed genes across all perturbations
genes <- unname(unlist(read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/final_targets_to_carry_forward.csv", header = FALSE)))
genes <- unname(sapply(genes, function(x) gsub("_", "-", x)))
path <- "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/DEG_Results/WithManualAdjustment/"
all_DESeq_res_list <- lapply(genes, function(x) read.table(paste0(path, x, "_vs_NTC_withDedup_manual_adjustment_Cutoff_Of_10.txt"), header = TRUE))
names(all_DESeq_res_list) <- genes

##filter out targeted gene from all_DESeq_res_list
for (gene in genes)
{
	res <- all_DESeq_res_list[[gene]]
	if (gene %in% row.names(res))
	{
		res <- res %>% filter(!row.names(res) == gene)
	}
	all_DESeq_res_list[[gene]] <- res
}
DEG_res_list <- lapply(all_DESeq_res_list, function(x) x %>% filter(padj < 0.05))
names(DEG_res_list) <- genes


#read in SFARI genes 
sfari <- read.csv("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/sfari.csv", header = TRUE)
SFARI_genes <- unique(sfari$gene.symbol)
all_SFARI_genes <- SFARI_genes


#read in ID genes 
id_genes <- read.csv("/n/groups/walsh/indData/becky/SysNDD_figures_250823/Input/SysNDD_ID_genes_definitive.csv", header = TRUE)
all_id_genes <- setdiff(id_genes$symbol, all_SFARI_genes)


###################
#PLOT DATAFRAME SETUP
###################
plot_df <- as.data.frame(t(as.data.frame(lapply(DEG_res_list, function(x) length(intersect(row.names(x), all_SFARI_genes)))))) %>% arrange(desc(V1))
colnames(plot_df) <- "Number_Of_SFARI_DEGs"
plot_df$Total_Number_Of_DEGs <- sapply(row.names(plot_df), function(x) as.numeric(nrow(DEG_res_list[[gsub(".", "-", x, fixed = TRUE)]])))
plot_df$Proportion_DEGs_SFARI_Genes <- plot_df$Number_Of_SFARI_DEGs/plot_df$Total_Number_Of_DEGs

#calculate enrichment p-value of DEGs and SFARI high confidence genes
calc_enrichment_test_p_value <- function(all_SFARI_genes, all_DESeq_res, DEG_res) 
{
	universe <- row.names(all_DESeq_res %>% filter(!is.na(padj)))
	DEGs <- row.names(DEG_res)
	all_SFARI_genes <- all_SFARI_genes[which(all_SFARI_genes %in% universe)]
	not_DEGs <- setdiff(universe, DEGs)
	not_all_SFARI_genes <- setdiff(universe, all_SFARI_genes)


	is_SFARI_gene_and_is_DEG <- length(intersect(DEGs, all_SFARI_genes))
	is_SFARI_gene_and_is_not_DEG <- length(intersect(not_DEGs, all_SFARI_genes))
	is_not_SFARI_gene_is_DEG <- length(intersect(DEGs, not_all_SFARI_genes))
	is_not_SFARI_gene_is_not_DEG <- length(intersect(not_DEGs, not_all_SFARI_genes))


	tab <- matrix(data = c(is_SFARI_gene_and_is_DEG, is_SFARI_gene_and_is_not_DEG, is_not_SFARI_gene_is_DEG, is_not_SFARI_gene_is_not_DEG), nrow = 2)


	return (c(fisher.test(tab, alternative = "greater")$p.value, fisher.test(tab, alternative = "greater")$estimate, fisher.test(tab, alternative = "greater")$conf.int[1]))
}
row.names(plot_df) <- sapply(row.names(plot_df), function(x) gsub(".", "-", x, fixed = TRUE))
plot_df$pval <- sapply(row.names(plot_df), function(x) calc_enrichment_test_p_value(all_SFARI_genes, all_DESeq_res_list[[x]], DEG_res_list[[x]])[1])
plot_df <- plot_df %>% arrange(pval)


#adjust for MHT for enrichment p-values 
plot_df$p_adj <- p.adjust(plot_df$pval, method = "BH")
plot_df <- plot_df %>% mutate(Pconvert = -log10(p_adj))
plot_df$gene <- row.names(plot_df)
plot_df$gene <- factor(as.character(plot_df$gene), levels = as.character(plot_df$gene))
plot_df$sig <- sapply(plot_df$p_adj, function(x) ifelse(x <0.0001, "< 0.0001", ifelse(x < 0.001, "< 0.001", ifelse(x < 0.01, "< 0.01", ifelse(x < 0.05, "< 0.05", "Not Significant")))))
plot_df$sig <- factor(as.character(plot_df$sig), levels = c("< 0.0001", "< 0.001", "< 0.01", "< 0.05", "Not Significant"))


#add OR
plot_df$OR <- sapply(row.names(plot_df), function(x) calc_enrichment_test_p_value(all_SFARI_genes, all_DESeq_res_list[[x]], DEG_res_list[[x]])[2])
plot_df <- plot_df %>% arrange(desc(OR))
plot_df$lower.ci <- sapply(row.names(plot_df), function(x) calc_enrichment_test_p_value(all_SFARI_genes, all_DESeq_res_list[[x]], DEG_res_list[[x]])[3])
plot_df$gene <- factor(as.character(plot_df$gene), levels = as.character(plot_df$gene))
sfari_plot_df <- plot_df

xlimit<-max(sfari_plot_df$OR)

###################
#CREATE PLOT
###################
#create plot
cols <- c(brewer.pal(8, "Greens")[c(8, 6)], "lightgray")
pdf("/n/groups/walsh/indData/becky/SysNDD_figures_250823/Fig1c/SFARI_gene_enrichment_barplot.pdf", height = 3, width = 5)
sfari_plot <- (ggplot(data = plot_df, aes(x = OR, y = forcats::fct_rev(gene), fill = sig, width = 0.8, color = "black")) + geom_bar(stat = "identity", color = "black") + 
theme_minimal() +
xlim(0,xlimit) +
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


###################
#PLOT DATAFRAME SETUP
###################
plot_df <- as.data.frame(t(as.data.frame(lapply(DEG_res_list, function(x) length(intersect(row.names(x), all_id_genes)))))) %>% arrange(desc(V1))
colnames(plot_df) <- "Number_Of_ID_DEGs"
plot_df$Total_Number_Of_DEGs <- sapply(row.names(plot_df), function(x) as.numeric(nrow(DEG_res_list[[gsub(".", "-", x, fixed = TRUE)]])))
plot_df$Proportion_DEGs_ID_Genes <- plot_df$Number_Of_ID_DEGs/plot_df$Total_Number_Of_DEGs


#calculate enrichment p-value of DEGs and id high confidence genes
calc_enrichment_test_p_value <- function(all_id_genes, all_DESeq_res, DEG_res) 
{
	universe <- row.names(all_DESeq_res %>% filter(!is.na(padj)))
	DEGs <- row.names(DEG_res)
	all_id_genes <- all_id_genes[which(all_id_genes %in% universe)]
	not_DEGs <- setdiff(universe, DEGs)
	not_all_id_genes <- setdiff(universe, all_id_genes)


	is_id_gene_and_is_DEG <- length(intersect(DEGs, all_id_genes))
	is_id_gene_and_is_not_DEG <- length(intersect(not_DEGs, all_id_genes))
	is_not_id_gene_is_DEG <- length(intersect(DEGs, not_all_id_genes))
	is_not_id_gene_is_not_DEG <- length(intersect(not_DEGs, not_all_id_genes))


	tab <- matrix(data = c(is_id_gene_and_is_DEG, is_id_gene_and_is_not_DEG, is_not_id_gene_is_DEG, is_not_id_gene_is_not_DEG), nrow = 2)


	return (c(fisher.test(tab, alternative = "greater")$p.value, fisher.test(tab, alternative = "greater")$estimate, fisher.test(tab, alternative = "greater")$conf.int[1]))
}
row.names(plot_df) <- sapply(row.names(plot_df), function(x) gsub(".", "-", x, fixed = TRUE))
plot_df$pval <- sapply(row.names(plot_df), function(x) calc_enrichment_test_p_value(all_id_genes, all_DESeq_res_list[[x]], DEG_res_list[[x]])[1])
plot_df <- plot_df %>% arrange(pval)


#adjust for MHT for enrichment p-values 
plot_df$p_adj <- p.adjust(plot_df$pval, method = "BH")
plot_df <- plot_df %>% mutate(Pconvert = -log10(p_adj))
plot_df$gene <- row.names(plot_df)
plot_df$gene <- factor(as.character(plot_df$gene), levels = as.character(plot_df$gene))
plot_df$sig <- sapply(plot_df$p_adj, function(x) ifelse(x <0.0001, "< 0.0001", ifelse(x < 0.001, "< 0.001", ifelse(x < 0.01, "< 0.01", ifelse(x < 0.05, "< 0.05", "Not Significant")))))
plot_df$sig <- factor(as.character(plot_df$sig), levels = c("< 0.0001", "< 0.001", "< 0.01", "< 0.05", "Not Significant"))


#add OR
plot_df$OR <- sapply(row.names(plot_df), function(x) calc_enrichment_test_p_value(all_id_genes, all_DESeq_res_list[[x]], DEG_res_list[[x]])[2])
plot_df <- plot_df %>% arrange(desc(OR))
plot_df$lower.ci <- sapply(row.names(plot_df), function(x) calc_enrichment_test_p_value(all_id_genes, all_DESeq_res_list[[x]], DEG_res_list[[x]])[3])
plot_df$gene <- factor(as.character(plot_df$gene), levels = as.character(sfari_plot_df$gene))
id_plot_df <- plot_df 


###################
#CREATE PLOT
###################
#create plot
cols <- "lightgray"
pdf("/n/groups/walsh/indData/becky/SysNDD_figures_250823/Fig1c/id_gene_enrichment_barplot.pdf", height = 3, width = 5)
id_plot <- (ggplot(data = plot_df, aes(x = OR, y = forcats::fct_rev(gene), fill = sig, width = 0.8, color = "black")) + geom_bar(stat = "identity", color = "black") + 
theme_minimal() +
xlim(0,xlimit) +
ylab("Target") + 
xlab("Odds Ratio") + geom_errorbar(aes(xmin = lower.ci, xmax = OR),  # only go down to lower CI
                  width = 0,                      # horizontal cap width
                  color = "black") +
theme(axis.text=element_text(size=7), axis.title=element_text(size=8,face="bold")) + 
scale_fill_manual(values = cols) + 
guides(fill=guide_legend(title="Adjusted p-value")) + 
theme(legend.text=element_text(size=7), legend.title=element_text(size=8,face="bold")) + 
guides(color = "none"))
print(id_plot)
dev.off()

###################
#CREATE COMBINED PLOT
###################
pdf("/n/groups/walsh/indData/becky/SysNDD_figures_250823/Fig1c/sfari_id_gene_enrichment_barplot.pdf", height = 3, width = 10)
wrap_plots(sfari_plot, id_plot, nrow = 1)
dev.off()

###################
#CREATE COMBINED SUPP TABLE
###################
# Add "_ID" to id_plot_df columns (if not already there)
colnames(id_plot_df) <- ifelse(grepl("ID", colnames(id_plot_df)),
                               colnames(id_plot_df),
                               paste0(colnames(id_plot_df), "_ID"))

# Add "_SFARI" to sfari_plot_df columns (if not already there)
colnames(sfari_plot_df) <- ifelse(grepl("SFARI", colnames(sfari_plot_df)),
                                  colnames(sfari_plot_df),
                                  paste0(colnames(sfari_plot_df), "_SFARI"))

# Merge on gene (make sure both have the same column name for merging)
merged_df <- merge(sfari_plot_df, id_plot_df, 
                   by.x = "gene_SFARI", 
                   by.y = "gene_ID", 
                   all = TRUE)

# Move gene column to the front
merged_df <- merged_df %>% dplyr::relocate(gene_SFARI, .before = 1)

merged_df <- merged_df %>%
  dplyr::mutate(Total_Number_Of_DEGs = Total_Number_Of_DEGs_SFARI) %>%
  dplyr::select(gene_SFARI, Total_Number_Of_DEGs, dplyr::everything(), 
                -Total_Number_Of_DEGs_SFARI, -Total_Number_Of_DEGs_ID)

# Rename gene column back to just "gene"
merged_df <- merged_df %>%
  dplyr::rename(gene = gene_SFARI) %>% arrange(desc(OR_SFARI))

write.table(merged_df, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Figure_1_And_S1/Figure_1C_Enrichment_Barplot_DEGs/supp_table_fig_1c.txt", sep = "\t", row.names = F, col.names = T, quote = F)