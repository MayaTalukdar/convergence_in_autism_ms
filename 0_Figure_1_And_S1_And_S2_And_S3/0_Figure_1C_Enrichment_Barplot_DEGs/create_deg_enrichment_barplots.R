###################
#I/O
###################
.libPaths(c("/n/groups/walsh/indData/Maya/RLibs", .libPaths()))
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(dplyr)

#read in differentially expressed genes across all perturbations
genes <- unname(unlist(read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/final_targets_to_carry_forward.csv", header = FALSE))) #In Input_Files_Not_Generated_By_Scripts 
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
sfari <- read.csv("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/sfari.csv", header = TRUE) #In Input_Files_Not_Generated_By_Scripts 
SFARI_genes <- unique(sfari$gene.symbol)
all_SFARI_genes <- SFARI_genes

#read in ID genes 
id_genes <- read.csv("/n/groups/walsh/indData/becky/SysNDD_figures_250823/Input/SysNDD_ID_genes_definitive.csv", header = TRUE) #In Input_Files_Not_Generated_By_Scripts 
all_id_genes <- setdiff(id_genes$symbol, all_SFARI_genes)

###################
#PLOT DATAFRAME SETUP -- SFARI
###################
plot_df <- as.data.frame(t(as.data.frame(lapply(DEG_res_list, function(x) length(intersect(row.names(x), all_SFARI_genes)))))) %>% arrange(desc(V1))
colnames(plot_df) <- "Number_Of_SFARI_DEGs"
plot_df$Total_Number_Of_DEGs <- sapply(row.names(plot_df), function(x) as.numeric(nrow(DEG_res_list[[gsub(".", "-", x, fixed = TRUE)]])))
plot_df$Proportion_DEGs_SFARI_Genes <- plot_df$Number_Of_SFARI_DEGs/plot_df$Total_Number_Of_DEGs

calc_enrichment_test_p_value <- function(gene_set, all_DESeq_res, DEG_res) 
{
	universe <- row.names(all_DESeq_res %>% filter(!is.na(padj)))
	DEGs <- row.names(DEG_res)
	gene_set_in_universe <- gene_set[which(gene_set %in% universe)]
	not_DEGs <- setdiff(universe, DEGs)
	not_gene_set <- setdiff(universe, gene_set_in_universe)

	#list1 = gene_set (SFARI or ID genes), list2 = DEGs

	InList1AndList2       <- length(intersect(DEGs, gene_set_in_universe))       # top-left
	InList1AndNotList2    <- length(intersect(not_DEGs, gene_set_in_universe))   # bottom-left
	InList2AndNotList1    <- length(intersect(DEGs, not_gene_set))               # top-right
	InNeither <- length(intersect(not_DEGs, not_gene_set))           # bottom-right

	tab <- matrix(data = c(InList1AndList2, InList1AndNotList2, InList2AndNotList1, InNeither),
	              nrow = 2,
	              dimnames = list(c("DEG", "Not_DEG"), c("GeneSet", "Not_GeneSet")))

	ft <- fisher.test(tab, alternative = "greater")

	return(c(pval = ft$p.value,
	         OR = unname(ft$estimate),
	         lower.ci = ft$conf.int[1],
	         InList1AndList2 = InList1AndList2,
	         InList1AndNotList2 = InList1AndNotList2,
	         InList2AndNotList1 = InList2AndNotList1,
	         InNeither = InNeither))
}

row.names(plot_df) <- sapply(row.names(plot_df), function(x) gsub(".", "-", x, fixed = TRUE))

enrichment_stats <- as.data.frame(t(sapply(row.names(plot_df), function(x)
	calc_enrichment_test_p_value(all_SFARI_genes, all_DESeq_res_list[[x]], DEG_res_list[[x]]))))
enrichment_stats <- enrichment_stats[row.names(plot_df), ]

plot_df <- cbind(plot_df, enrichment_stats)
plot_df <- plot_df %>% arrange(pval)

#adjust for MHT for enrichment p-values 
plot_df$p_adj <- p.adjust(plot_df$pval, method = "BH")
plot_df <- plot_df %>% mutate(Pconvert = -log10(p_adj))
plot_df$gene <- row.names(plot_df)
plot_df$gene <- factor(as.character(plot_df$gene), levels = as.character(plot_df$gene))
plot_df$sig <- sapply(plot_df$p_adj, function(x) ifelse(x <0.0001, "< 0.0001", ifelse(x < 0.001, "< 0.001", ifelse(x < 0.01, "< 0.01", ifelse(x < 0.05, "< 0.05", "Not Significant")))))
plot_df$sig <- factor(as.character(plot_df$sig), levels = c("< 0.0001", "< 0.001", "< 0.01", "< 0.05", "Not Significant"))

plot_df <- plot_df %>% arrange(desc(OR))
plot_df$gene <- factor(as.character(plot_df$gene), levels = as.character(plot_df$gene))
sfari_plot_df <- plot_df

xlimit<-max(sfari_plot_df$OR)

###################
#CREATE PLOT -- SFARI
###################
cols <- c(brewer.pal(8, "Greens")[c(8, 6)], "lightgray")
pdf("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Actually_Submitted/For_Third_Submission_And_What_We_Published/0_CLEANED_0_Figure_1_And_S1_And_S2_And_S3/0_Figure_1C_Enrichment_Barplot_DEGs/SFARI_gene_enrichment_barplot.pdf", height = 3, width = 5)
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
#PLOT DATAFRAME SETUP -- ID
###################
plot_df <- as.data.frame(t(as.data.frame(lapply(DEG_res_list, function(x) length(intersect(row.names(x), all_id_genes)))))) %>% arrange(desc(V1))
colnames(plot_df) <- "Number_Of_ID_DEGs"
plot_df$Total_Number_Of_DEGs <- sapply(row.names(plot_df), function(x) as.numeric(nrow(DEG_res_list[[gsub(".", "-", x, fixed = TRUE)]])))
plot_df$Proportion_DEGs_ID_Genes <- plot_df$Number_Of_ID_DEGs/plot_df$Total_Number_Of_DEGs

#calculate enrichment p-value of DEGs and id high confidence genes (same logic as above,
#returns the relevant 2x2 contingency table counts alongside pval/OR/CI)
calc_enrichment_test_p_value <- function(gene_set, all_DESeq_res, DEG_res) 
{
	universe <- row.names(all_DESeq_res %>% filter(!is.na(padj)))
	DEGs <- row.names(DEG_res)
	gene_set_in_universe <- gene_set[which(gene_set %in% universe)]
	not_DEGs <- setdiff(universe, DEGs)
	not_gene_set <- setdiff(universe, gene_set_in_universe)

	#list1 = gene_set (SFARI or ID genes), list2 = DEGs

	InList1AndList2       <- length(intersect(DEGs, gene_set_in_universe))       # top-left
	InList1AndNotList2    <- length(intersect(not_DEGs, gene_set_in_universe))   # bottom-left
	InList2AndNotList1    <- length(intersect(DEGs, not_gene_set))               # top-right
	InNeither <- length(intersect(not_DEGs, not_gene_set))           # bottom-right

	tab <- matrix(data = c(InList1AndList2, InList1AndNotList2, InList2AndNotList1, InNeither),
	              nrow = 2,
	              dimnames = list(c("DEG", "Not_DEG"), c("GeneSet", "Not_GeneSet")))

	ft <- fisher.test(tab, alternative = "greater")

	return(c(pval = ft$p.value,
	         OR = unname(ft$estimate),
	         lower.ci = ft$conf.int[1],
	         InList1AndList2 = InList1AndList2,
	         InList1AndNotList2 = InList1AndNotList2,
	         InList2AndNotList1 = InList2AndNotList1,
	         InNeither = InNeither))
}

row.names(plot_df) <- sapply(row.names(plot_df), function(x) gsub(".", "-", x, fixed = TRUE))

enrichment_stats <- as.data.frame(t(sapply(row.names(plot_df), function(x)
	calc_enrichment_test_p_value(all_id_genes, all_DESeq_res_list[[x]], DEG_res_list[[x]]))))
enrichment_stats <- enrichment_stats[row.names(plot_df), ]

plot_df <- cbind(plot_df, enrichment_stats)
plot_df <- plot_df %>% arrange(pval)

#adjust for MHT for enrichment p-values 
plot_df$p_adj <- p.adjust(plot_df$pval, method = "BH")
plot_df <- plot_df %>% mutate(Pconvert = -log10(p_adj))
plot_df$gene <- row.names(plot_df)
plot_df$gene <- factor(as.character(plot_df$gene), levels = as.character(plot_df$gene))
plot_df$sig <- sapply(plot_df$p_adj, function(x) ifelse(x <0.0001, "< 0.0001", ifelse(x < 0.001, "< 0.001", ifelse(x < 0.01, "< 0.01", ifelse(x < 0.05, "< 0.05", "Not Significant")))))
plot_df$sig <- factor(as.character(plot_df$sig), levels = c("< 0.0001", "< 0.001", "< 0.01", "< 0.05", "Not Significant"))

plot_df <- plot_df %>% arrange(desc(OR))
plot_df$gene <- factor(as.character(plot_df$gene), levels = as.character(sfari_plot_df$gene))
id_plot_df <- plot_df 

###################
#CREATE PLOT -- ID
###################
cols <- "lightgray"
pdf("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Actually_Submitted/For_Third_Submission_And_What_We_Published/0_CLEANED_0_Figure_1_And_S1_And_S2_And_S3/0_Figure_1C_Enrichment_Barplot_DEGs/id_gene_enrichment_barplot.pdf", height = 3, width = 5)
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
pdf("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Actually_Submitted/For_Third_Submission_And_What_We_Published/0_CLEANED_0_Figure_1_And_S1_And_S2_And_S3/0_Figure_1C_Enrichment_Barplot_DEGs/sfari_id_gene_enrichment_barplot.pdf", height = 3, width = 10)
wrap_plots(sfari_plot, id_plot, nrow = 1)
dev.off()

###################
#CREATE COMBINED SUPP TABLE
###################
colnames(id_plot_df) <- ifelse(grepl("ID", colnames(id_plot_df)),
                               colnames(id_plot_df),
                               paste0(colnames(id_plot_df), "_ID"))

# Add "_SFARI" to sfari_plot_df columns (if not already there) -- same automatic pickup
# of the new matrix columns
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

# Rename gene column back to just "gene", and sort by SFARI OR (unchanged from before)
merged_df <- merged_df %>%
  dplyr::rename(gene = gene_SFARI) %>% arrange(desc(OR_SFARI))

###################
#RENAME COLUMNS TO PLAIN-ENGLISH / SUPPLEMENTARY-TABLE-READY NAMES
###################
# note: matrix columns (InList1AndList2 / InList1AndNotList2 / InList2AndNotList1 / InNeither)
# are left out of this map on purpose so they keep their raw name + gene set suffix
# (e.g. InList1AndList2_SFARI) instead of a made-up English phrase. InList1AndList2 is the
# same number as SFARI DEGs / ID DEGs above, just reported again in the matrix nomenclature
friendly_names <- c(
  "gene"                            = "Gene",
  "Total_Number_Of_DEGs"            = "Total DEGs",

  #--- SFARI columns ---
  "Number_Of_SFARI_DEGs"            = "SFARI DEGs",
  "Proportion_DEGs_SFARI_Genes"     = "Proportion SFARI",
  "OR_SFARI"                        = "SFARI Odds Ratio",
  "lower.ci_SFARI"                  = "SFARI Odds Ratio Lower CI",
  "pval_SFARI"                      = "SFARI P-value",
  "p_adj_SFARI"                     = "SFARI Adjusted P-value",
  "Pconvert_SFARI"                  = "SFARI -log10(Adjusted P-value)",
  "sig_SFARI"                       = "SFARI Significance",

  #--- ID columns ---
  "Number_Of_ID_DEGs"               = "ID DEGs",
  "Proportion_DEGs_ID_Genes"        = "Proportion ID",
  "OR_ID"                           = "ID Odds Ratio",
  "lower.ci_ID"                     = "ID Odds Ratio Lower CI",
  "pval_ID"                         = "ID P-value",
  "p_adj_ID"                        = "ID Adjusted P-value",
  "Pconvert_ID"                     = "ID -log10(Adjusted P-value)",
  "sig_ID"                          = "ID Significance"
)

matched <- match(names(friendly_names), colnames(merged_df))
stopifnot(!any(is.na(matched)))
colnames(merged_df)[matched] <- unname(friendly_names)

write.table(merged_df, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Actually_Submitted/For_Third_Submission_And_What_We_Published/0_CLEANED_0_Figure_1_And_S1_And_S2_And_S3/0_Figure_1C_Enrichment_Barplot_DEGs/supp_table_fig_1c.txt", sep = "\t", row.names = F, col.names = T, quote = F)