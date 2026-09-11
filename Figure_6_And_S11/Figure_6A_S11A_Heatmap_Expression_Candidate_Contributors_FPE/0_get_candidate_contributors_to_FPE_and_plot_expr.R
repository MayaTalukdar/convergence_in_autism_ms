#################
#I/O
#################
.libPaths(c("/n/groups/walsh/indData/Maya/RLibs", .libPaths()))
library(tidyverse)
library(readxl)
library(data.table)
library(org.Hs.eg.db)
library(pheatmap)
library(RColorBrewer)
library(biomaRt)

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

  print(mat)

	return(pval)
}

generate_overlap_or <- function(list1, list2, universe)
{
	list1 <- intersect(list1, universe)
	list2 <- intersect(list2, universe)

	inList1AndList2 <- length(intersect(list1, list2))
	inList1AndNotList2 <- length(setdiff(list1, list2))
	inList2AndNotList1 <- length(setdiff(list2, list1))
	inNeither <- length(setdiff(universe, c(list1, list2)))

	mat <- matrix(c(inList1AndList2, inList1AndNotList2, inList2AndNotList1, inNeither), nrow = 2)
	stat <- fisher.test(mat, alternative = "greater")$estimate

	return(stat)
}

#read in SFARI genes 
sfari <- read.csv("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/sfari.csv", header = TRUE)
SFARI_genes <- unique(sfari$gene.symbol)

#define the universe as all npc genes expressed with baseMean >= 10 
basemean_df <-  read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/baseMeanInNTCs.csv", row.names = 1)
universe <- row.names(basemean_df %>% filter(baseMean > 10))

#read in perturbation data
file_path <- "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/DEG_Results/WithManualAdjustment/"
sgRNAs <- unname(unlist(read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/final_targets_to_carry_forward.csv")))
print(length(sgRNAs) == 18)
all_DESeq_res <- lapply(sgRNAs, function(x) read.table(paste0(file_path, x, "_vs_NTC_withDedup_manual_adjustment_Cutoff_Of_10.txt")))
names(all_DESeq_res) <- sgRNAs

#extract only the psDEGs
psDEG_list <- list()
all_psDEGs <- lapply(all_DESeq_res, function(x) row.names(x %>% filter(padj < 0.05)))

#read in constraint scores
constraint_scores <- fread("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/supplementary_dataset_11_full_constraint_metrics.tsv", header = TRUE)
constraint_scores <- as.data.frame(constraint_scores)
constraint_scores <- constraint_scores %>% as.data.frame() %>% dplyr::select(oe_lof_upper, gene) %>% filter(!is.na(oe_lof_upper))
constraint_scores <- constraint_scores %>% group_by(gene) %>% summarize(median_oe_lof_upper_across_transcripts = median(oe_lof_upper))
ploeuf_leq_0p35_genes <- constraint_scores %>% filter(median_oe_lof_upper_across_transcripts < 0.35) %>% pull(gene) %>% as.character()

#################
#IDENTIFY FEMALE-BIASED GENES
#################
#read in female biased genes 
kissel_2024_data <- read.csv("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/figures_for_manuscript/get_sex_assoc_regs/X_Linked_Gene_Classes/kissel_2024_supp_table_4.csv", header = TRUE)

#read in naqvi conserved sex bias data
naqvi_s4_data <- read.csv("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/figures_for_manuscript/get_sex_assoc_regs/X_Linked_Gene_Classes/naqvi_table_s4.csv", header = TRUE)

#read in oliva gtex sex biased data
oliva_data <- fread("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/figures_for_manuscript/get_sex_assoc_regs/X_Linked_Gene_Classes/GTEx_Analysis_v8_sbgenes 2/signif.sbgenes.txt", header = TRUE) %>% as.data.frame()

#read in decasien data
decasien_data <- read.csv("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Revisions/Analysis_of_Decasien_Data/Raw_Data/decasien_table_s13_ubiquitously_sex_biased_genes.csv", header = TRUE)

#set up list 
potential_contributors_to_the_fpe_list <- list()

#KISSEL
##get UCLA sex-biased genes
ucla_fem_biased_genes <- kissel_2024_data %>% filter(adj.P.Val_UCLA < 0.1) %>% filter(logFC_UCLA < 0) %>% filter(Biotype != "pseudogene")
ucla_male_biased_genes <- kissel_2024_data %>% filter(adj.P.Val_UCLA < 0.1) %>% filter(logFC_UCLA > 0)  %>% filter(Biotype != "pseudogene")

##get BV sex-biased genes
BV_fem_biased_genes <- kissel_2024_data %>% filter(adj.P.Val_BV < 0.1) %>% filter(logFC_BV < 0)   %>% filter(Biotype != "pseudogene")
BV_male_biased_genes <- kissel_2024_data %>% filter(adj.P.Val_BV < 0.1) %>% filter(logFC_BV > 0)   %>% filter(Biotype != "pseudogene")

##get inv sex biased genes
kissel_2024_data %>%  filter(Biotype != "pseudogene") %>% filter(invNorm_adjPval < 0.1) %>% filter(logFC_UCLA < 0 & logFC_BV > 0 | logFC_UCLA > 0 & logFC_BV < 0) %>% dplyr::select(gene_name, logFC_BV, adj.P.Val_BV, logFC_UCLA, adj.P.Val_UCLA, invNorm_adjPval)
discordant_genes <- c("UCK2", "PCSK1")
meta_fem_biased_genes <- kissel_2024_data %>%  filter(Biotype != "pseudogene") %>% filter(invNorm_adjPval < 0.1) %>% filter(logFC_UCLA < 0 | logFC_BV < 0) %>% filter(!gene_name %in% discordant_genes)

##get fem-biased genes
fem_biased_genes <- unique(c(ucla_fem_biased_genes %>% pull(gene_name) %>% as.character(), 
BV_fem_biased_genes %>% pull(gene_name) %>% as.character(), 
meta_fem_biased_genes %>% pull(gene_name) %>% as.character()))

#NAQVI 
conserved_female_biased_genes <- naqvi_s4_data %>% filter(Brain == -1) %>% pull(Row.names) %>% as.character()
human_gained_female_biased_genes <-  naqvi_s4_data %>% filter(Brain == -6) %>% pull(Row.names) %>% as.character()

#OLIVA GTEX
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gtex_female_biased_genes_ensembl <- unname(sapply(oliva_data %>% filter(grepl("Cortex", tissue)) %>% filter(grepl("Brain", tissue)) %>% filter(effsize > 0) %>% pull(gene) %>% as.character(), function(x) strsplit(x, ".", fixed = TRUE)[[1]][1]))
gtex_female_biased_genes <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                    filters = "ensembl_gene_id", 
                    values = gtex_female_biased_genes_ensembl, 
                    mart = ensembl) %>% filter(external_gene_name != "") %>% pull(external_gene_name) %>% as.character()

#DECASIEN
decasien_female_biased_genes <- decasien_data %>% filter(Bias == "female") %>% pull(Unique.gene.list)

potential_contributors_to_the_fpe_list[["secondary_contributors"]] <- sort(unique(c(fem_biased_genes, conserved_female_biased_genes, gtex_female_biased_genes, decasien_female_biased_genes)))
saveRDS(potential_contributors_to_the_fpe_list, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Revisions/Redo_Figure_6_And_S6_To_Add_Decasien_Data/Figure_6A_S6A_Heatmap_Expression_Candidate_Contributors_FPE/female_biased_genes.RDS")

#################
#GENERATE LIST OF CANDIDATE ASD RISK GENES
#################
#filter to geq 3 perturbations
freq_table_of_genes_as_psDEGs <- table(do.call(c,all_psDEGs)) %>% as.data.frame() %>% arrange(desc(Freq)) %>% filter(Freq >= 3) %>% pull(Var1) %>% as.character()

#filter to genes that show at least 80% concordance in direction of sig DEGs across all ASD KDs
plot_df <- sapply(all_DESeq_res, function(x) (x[freq_table_of_genes_as_psDEGs,])$log2FoldChange)
row.names(plot_df) <- freq_table_of_genes_as_psDEGs
plot_df <- as.data.frame(t(plot_df))

sig_df <- sapply(all_DESeq_res, function(x) (x[freq_table_of_genes_as_psDEGs,])$padj)
row.names(sig_df) <- freq_table_of_genes_as_psDEGs
sig_df <- apply(sig_df,2, function(x) ifelse(x < 0.05, "*", ""))
sig_df <- as.data.frame(t(sig_df))

consistent_downreg_DEG_in_asd_kd_genes <- c()
for (gene in colnames(plot_df)) 
{
  logfc_vals <- plot_df[[gene]]
  sig_vals <- sig_df[[gene]]

  sig_mask <- sig_vals == "*"
  sig_logfc <- logfc_vals[sig_mask]

  if (length(sig_logfc) == 0) next

  prop_up <- mean(sig_logfc > 0)
  prop_down <- mean(sig_logfc < 0)

  if (prop_down >= 0.8) 
  {
    consistent_downreg_DEG_in_asd_kd_genes <- c(consistent_downreg_DEG_in_asd_kd_genes, gene)
  }
}

#################
#GENERATE FINAL LIST OF POTENTIAL CONTRIBUTORS TO THE FPE
#################
consistent_downreg_DEG_in_asd_kd_genes_and_potential_contributors_to_the_fpe <- intersect(universe, unique(intersect(consistent_downreg_DEG_in_asd_kd_genes, unique(do.call(c, potential_contributors_to_the_fpe_list)))))
saveRDS(consistent_downreg_DEG_in_asd_kd_genes_and_potential_contributors_to_the_fpe, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Revisions/Redo_Figure_6_And_S6_To_Add_Decasien_Data/Figure_6A_S6A_Heatmap_Expression_Candidate_Contributors_FPE/prioritized_regs_consistent_downreg_DEG_in_asd_kd_genes_and_potential_contributors_to_the_fpe.RDS")
consistent_downreg_DEG_in_asd_kd_genes_and_potential_contributors_to_the_fpe_constrained <- intersect(consistent_downreg_DEG_in_asd_kd_genes_and_potential_contributors_to_the_fpe, ploeuf_leq_0p35_genes)
consistent_downreg_DEG_in_asd_kd_genes_and_potential_contributors_to_the_fpe_constrained_sfari <- intersect(consistent_downreg_DEG_in_asd_kd_genes_and_potential_contributors_to_the_fpe_constrained, SFARI_genes)
saveRDS(consistent_downreg_DEG_in_asd_kd_genes_and_potential_contributors_to_the_fpe_constrained , "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Revisions/Redo_Figure_6_And_S6_To_Add_Decasien_Data/Figure_6A_S6A_Heatmap_Expression_Candidate_Contributors_FPE/prioritized_regs_consistent_downreg_DEG_in_asd_kd_genes_and_potential_contributors_to_the_fpe_constrained.RDS")

#run enrichment test of candidate asd risk genes and sfari genes
print(length(consistent_downreg_DEG_in_asd_kd_genes_and_potential_contributors_to_the_fpe_constrained_sfari))
print(length(consistent_downreg_DEG_in_asd_kd_genes_and_potential_contributors_to_the_fpe_constrained))
universe_constraint <- intersect(constraint_scores %>% filter(median_oe_lof_upper_across_transcripts < 0.35) %>% pull(gene) %>% as.character(), universe) 
print(paste0('p-value of candidate asd genes & sfari genes, using a universe of NPC-expressed genes with pLOEUF <0.35: ', generate_overlap_p_val(SFARI_genes, consistent_downreg_DEG_in_asd_kd_genes_and_potential_contributors_to_the_fpe_constrained, universe_constraint)))
print(paste0('or of candidate asd genes & sfari genes, using a universe of NPC-expressed genes with pLOEUF <0.35: ', generate_overlap_or(SFARI_genes, consistent_downreg_DEG_in_asd_kd_genes_and_potential_contributors_to_the_fpe_constrained, universe_constraint)))

#################
#PLOT HEATMAP - ALL 
#################
regs_for_figures <- consistent_downreg_DEG_in_asd_kd_genes_and_potential_contributors_to_the_fpe_constrained

#get matrices to plot 
plot_df <- sapply(all_DESeq_res, function(x) (x[regs_for_figures,])$log2FoldChange)
row.names(plot_df) <- regs_for_figures
plot_df <- as.data.frame(t(plot_df))

sig_df <- sapply(all_DESeq_res, function(x) (x[regs_for_figures,])$padj)
row.names(sig_df) <- regs_for_figures
sig_df <- apply(sig_df,2, function(x) ifelse(x < 0.05, "*", ""))
sig_df <- as.data.frame(t(sig_df))

#get information to nominate autism risk genes - coexpression 
sfari <- read.csv("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/sfari.csv")
sfari <- sfari %>% filter(gene.score == 1)
sfari_og <- sfari$gene.symbol

baseMean_in_NTCs <- read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/baseMeanInNTCs.csv", row.names = 1)
genes_to_keep <- row.names(baseMean_in_NTCs %>% filter(baseMean >= 10))

##perturbation analysis
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

#create annotation dataframe
annot_df <- sapply(regs_for_figures, function(x) ifelse(x %in% SFARI_genes, "TRUE", "FALSE")) %>% as.data.frame()
row.names(annot_df) <- colnames(plot_df)
annot_df <- cbind(annot_df, reg_corr_annot_perturbation_sig)
colnames(annot_df) <- c("Is SFARI?","Sig. Co-expression with SFARI Genes? (Perturbations)")

annotation_colors_diff_expr <- c("Yes" = "black" , "No" = "white")
ann_colors <- list(
  "Is SFARI?" = c("TRUE" = "black", "FALSE" = "white"),
  "Sig. Co-expression with SFARI Genes? (Perturbations)" = annotation_colors_diff_expr
)

#plot  
paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
plot_min <- min(plot_df, na.rm = TRUE)
plot_max <- max(plot_df, na.rm = TRUE)
plot_limit <- max(abs(plot_min), abs(plot_max))
myBreaks <- c(seq(-plot_limit, 0, length.out = ceiling(paletteLength / 2) + 1),
              seq(plot_limit / paletteLength, plot_limit, length.out = floor(paletteLength / 2)))

pheatmap(
  mat = plot_df,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "none",
  color = myColor,
  annotation_col = annot_df,
  annotation_colors = ann_colors,
  breaks = myBreaks,
  border_color = NA, 
  display_numbers = sig_df,
  fontsize_row = 8, fontsize_col = 8,
  filename = "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Revisions/Redo_Figure_6_And_S6_To_Add_Decasien_Data/Figure_6A_S6A_Heatmap_Expression_Candidate_Contributors_FPE/consistent_downreg_DEG_in_asd_kd_genes_and_potential_contributors_to_the_fpe_heatmap.pdf",
  width = 15, height = 6
)

saveRDS(plot_df, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Revisions/Redo_Figure_6_And_S6_To_Add_Decasien_Data/Figure_6A_S6A_Heatmap_Expression_Candidate_Contributors_FPE/fig_s6a_plot_df_for_heatmap.RDS")


#################
#PLOT HEATMAP - JUST SFARI GENES
#################
regs_for_figures <- intersect(consistent_downreg_DEG_in_asd_kd_genes_and_potential_contributors_to_the_fpe_constrained, SFARI_genes)

#get matrices to plot 
plot_df <- sapply(all_DESeq_res, function(x) (x[regs_for_figures,])$log2FoldChange)
row.names(plot_df) <- regs_for_figures
plot_df <- as.data.frame(t(plot_df))

sig_df <- sapply(all_DESeq_res, function(x) (x[regs_for_figures,])$padj)
row.names(sig_df) <- regs_for_figures
sig_df <- apply(sig_df,2, function(x) ifelse(x < 0.05, "*", ""))
sig_df <- as.data.frame(t(sig_df))

#get information to nominate autism risk genes - coexpression 
sfari <- read.csv("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/sfari.csv")
sfari <- sfari %>% filter(gene.score == 1)
sfari_og <- sfari$gene.symbol

baseMean_in_NTCs <- read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/baseMeanInNTCs.csv", row.names = 1)
genes_to_keep <- row.names(baseMean_in_NTCs %>% filter(baseMean >= 10))

##perturbation analysis
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

#create annotation dataframe
annot_df <- sapply(regs_for_figures, function(x) ifelse(x %in% SFARI_genes, "TRUE", "FALSE")) %>% as.data.frame()
row.names(annot_df) <- colnames(plot_df)
annot_df <- cbind(annot_df, reg_corr_annot_perturbation_sig)
colnames(annot_df) <- c("Is SFARI?","Sig. Co-expression with SFARI Genes? (Perturbations)")

annotation_colors_diff_expr <- c("Yes" = "black" , "No" = "white")
ann_colors <- list(
  "Is SFARI?" = c("TRUE" = "black", "FALSE" = "white"),
  "Sig. Co-expression with SFARI Genes? (Perturbations)" = annotation_colors_diff_expr
)

#plot  
paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
plot_min <- min(plot_df, na.rm = TRUE)
plot_max <- max(plot_df, na.rm = TRUE)
plot_limit <- max(abs(plot_min), abs(plot_max))
myBreaks <- c(seq(-plot_limit, 0, length.out = ceiling(paletteLength / 2) + 1),
              seq(plot_limit / paletteLength, plot_limit, length.out = floor(paletteLength / 2)))

pheatmap(
  mat = plot_df,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "none",
  color = myColor,
  annotation_col = annot_df,
  annotation_colors = ann_colors,
  breaks = myBreaks,
  border_color = NA, 
  display_numbers = sig_df,
  fontsize_row = 8, fontsize_col = 8,
  filename = "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Revisions/Redo_Figure_6_And_S6_To_Add_Decasien_Data/Figure_6A_S6A_Heatmap_Expression_Candidate_Contributors_FPE/consistent_downreg_DEG_in_asd_kd_genes_and_potential_contributors_to_the_fpe_heatmap_just_SFARI.pdf",
  width = 15, height = 6
)

saveRDS(plot_df, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Revisions/Redo_Figure_6_And_S6_To_Add_Decasien_Data/Figure_6A_S6A_Heatmap_Expression_Candidate_Contributors_FPE/fig_6a_plot_df_for_heatmap.RDS")

#################
#GET JUST DOWNREG GENES P-VALUES
#################
