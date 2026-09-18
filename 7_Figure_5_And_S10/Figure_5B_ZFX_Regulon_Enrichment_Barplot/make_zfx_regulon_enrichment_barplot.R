##################
#I/O
##################
.libPaths(c("/n/groups/walsh/indData/Maya/RLibs", .libPaths()))
library(tidyverse)
library(readxl)    
library(ggplot2)
library(data.table)
library(ggrepel)
library(RColorBrewer)
library(igraph)

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


	#also return the matrix values so they can be reported as source data for reviewers
	return(c(pval, estimate, fisher.test(mat, alternative = "greater")$conf.int[1], inList1AndList2, inList1AndNotList2, inList2AndNotList1, inNeither))
}


#read in network 
network <- read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/full_aracne_network_analysis/Output/All_TFs_And_Epigenetic_Regulators/aracne_output/network_formatted.txt", header = FALSE)
colnames(network) <- c('Regulator', 'Regulatee', 'weight') #note: weight is MI
all_network_genes <- unique(c(network$Regulatee, network$Regulator))

#read in candidate list
final_candidates <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/full_aracne_network_analysis/Output/Analysis_Results/GenePrioritization/final_candidates.rds")

#set up pathway list
#read in SFARI genes 
sfari <- read.csv("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/sfari.csv", header = TRUE) #In Input_Files_Not_Generated_By_Scripts
SFARI_genes <- unique(sfari$gene.symbol)
high_conf_SFARI_genes <- unique((sfari %>% filter(gene.score == 1))$gene.symbol)

#read in ID genes 
id_genes <- read.csv("/n/groups/walsh/indData/becky/SysNDD_figures_250823/Input/SysNDD_ID_genes_definitive.csv", header = TRUE) #In Input_Files_Not_Generated_By_Scripts
just_id_genes <- setdiff(id_genes$symbol, SFARI_genes)

##constraint scores
constraint_scores <- fread("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/supplementary_dataset_11_full_constraint_metrics.tsv", header = TRUE) #In Input_Files_Not_Generated_By_Scripts
constraint_scores <- as.data.frame(constraint_scores)
constraint_scores <- constraint_scores %>% as.data.frame() %>% dplyr::select(oe_lof_upper, gene) %>% filter(!is.na(oe_lof_upper))
constraint_scores <- constraint_scores %>% group_by(gene) %>% summarize(median_oe_lof_upper_across_transcripts = median(oe_lof_upper))
ploeuf_leq_0p35_genes <- constraint_scores %>% filter(median_oe_lof_upper_across_transcripts < 0.35) %>% pull(gene) %>% as.character()

#consistently downregulated genes
consistent_downreg_DEG_in_asd_kd_genes <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/figures_for_manuscript/get_sex_assoc_regs/0_get_consistent_ASD_DEGs/consistent_downreg_DEG_in_asd_kd_genes.RDS")
consistent_downreg_DEG_in_asd_kd_genes_constrained <- intersect(consistent_downreg_DEG_in_asd_kd_genes, ploeuf_leq_0p35_genes)

#female-biased downDEGs
female_biased_downDEG <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Revisions/Redo_Figure_6_And_S6_To_Add_Decasien_Data/Figure_6A_S6A_Heatmap_Expression_Candidate_Contributors_FPE/prioritized_regs_consistent_downreg_DEG_in_asd_kd_genes_and_potential_contributors_to_the_fpe.RDS")

#put ontologies in a single list
pathway_list <- list()
pathway_list[["SFARI"]] <- SFARI_genes
pathway_list[["HCSFARI"]] <- high_conf_SFARI_genes
pathway_list[["Only ID"]] <- just_id_genes
pathway_list[["Candidate ASD Risk"]] <- setdiff(consistent_downreg_DEG_in_asd_kd_genes_constrained, SFARI_genes)
pathway_list[["Downregulated DEGs"]] <- consistent_downreg_DEG_in_asd_kd_genes
pathway_list[["Top Drivers"]] <- final_candidates
pathway_list[["Female-Biased Down. DEGs"]] <- female_biased_downDEG

#get regulon list
regulon_list <- setNames(
  lapply(c("ZFX"), function(x) {
    first_degree_reg <- unique(network %>% filter(Regulator == x) %>% pull(Regulatee))
    second_degree_reg <- unique(network %>% filter(Regulator %in% first_degree_reg) %>% pull(Regulatee))
    unique(c(first_degree_reg, second_degree_reg))
  }),
  c("ZFX")
)

#get dataframe (list1 = ZFX regulon, list2 = each gene set/ontology in pathway_list)
enr_df <- sapply(pathway_list, function(x) generate_overlap_p_val_and_or(regulon_list[["ZFX"]], x, all_network_genes)) %>% as.data.frame() %>% t() %>% as.data.frame()
enr_df$regulator <- "ZFX"
colnames(enr_df) <- c("pval", "or", "lower.ci", "InList1AndList2", "InList1AndNotList2", "InList2AndNotList1", "InNeither", "regulator")
enr_df$padj <- p.adjust(enr_df$pval, method = "BH")
enr_df <- enr_df %>% arrange(desc(or))
enr_df$ontology <- sapply(row.names(enr_df), function(x) gsub("1", "", gsub("2", "", x, fixed = TRUE), fixed = TRUE))
row.names(enr_df) <- NULL

#write out the matrix values as source data for the reviewer comment
write.table(enr_df, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Actually_Submitted/For_Third_Submission_And_What_We_Published/7_CLEANED_7_Figure_5_And_S10/Figure_5B_ZFX_Regulon_Enrichment_Barplot/ZFX_regulon_enr_source_data.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


#set up plot
enr_df$ontology <- factor(enr_df$ontology, levels = unique(as.character(enr_df$ontology)))
enr_df$sig <- sapply(enr_df$padj, function(x) ifelse(x <0.0001, "< 0.0001", ifelse(x < 0.001, "< 0.001", ifelse(x < 0.01, "< 0.01", ifelse(x < 0.05, "< 0.05", "Not Significant")))))
enr_df$sig <- factor(as.character(enr_df$sig), 
                                levels = c("< 0.0001", "< 0.001", "< 0.01", "< 0.05", "Not Significant"))
cols <- c(brewer.pal(8, "Greens")[c(8, 6, 4, 2)], "lightgray")
names(cols) <- c("< 0.0001", "< 0.001", "< 0.01", "< 0.05", "Not Significant")
enr_df$regulator <- forcats::fct_rev(factor(enr_df$regulator, levels = c("ZFX")))
print(enr_df)


pdf("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Actually_Submitted/For_Third_Submission_And_What_We_Published/7_CLEANED_7_Figure_5_And_S10/Figure_5B_ZFX_Regulon_Enrichment_Barplot/ZFX_regulon_enr.pdf", width = 6, height = 3)
ggplot(enr_df, aes(y = forcats::fct_rev(ontology), x = or, fill = sig, group = regulator)) + 
  geom_bar(stat = "identity", 
           position = position_dodge(width = 0.8), 
           color = "black", 
           width = 0.7) + 
  geom_errorbar(aes(xmin = lower.ci, xmax = or),
                position = position_dodge(width = 0.8),  # <-- same dodge here
                width = 0,
                color = "black") +
  theme_minimal() +
  ylab("Ontology") + 
  xlab("Odds Ratio") + 
  theme(axis.text = element_text(size = 7), 
        axis.title = element_text(size = 8, face = "bold")) + 
  scale_fill_manual(values = cols) + 
  guides(fill = guide_legend(title = "Adjusted p-value")) + 
  theme(legend.text = element_text(size = 7), 
        legend.title = element_text(size = 8, face = "bold"))
dev.off()