##########################
#I/O
##########################
.libPaths(c("/n/groups/walsh/indData/Maya/RLibs", .libPaths()))
library(viper)
library(data.table)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(grid)
library(igraph)
library(fgsea)
library(readxl)
library(RColorBrewer)
library(readxl)
library(Orthology.eg.db)
library(org.Mm.eg.db)
library(org.Hs.eg.db)


#functions
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


#read in target genes & degs
targets <- unname(unlist(read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/final_targets_to_carry_forward.csv"))) #Input_Files_Not_Generated_By_Scripts
file_path <- "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/DEG_Results/WithManualAdjustment/"
all_DESeq_res <- lapply(targets, function(x) read.table(paste0(file_path, x, "_vs_NTC_withDedup_manual_adjustment_Cutoff_Of_10.txt")))
names(all_DESeq_res) <- targets
all_psDEGs <- lapply(all_DESeq_res, function(x) row.names(x %>% filter(padj < 0.05)))

#read in chd8 data
chd8_data <- read.csv("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/sugathan_data_chd8.csv", header = TRUE) #In Input_Files_Not_Generated_By_Scripts
colnames(chd8_data)[12] <- "has_CHIP_seq_binding_site"
colnames(chd8_data)[6] <- "log2FC"
colnames(chd8_data)[10] <- "padj"
chd8_data_sugathan <- chd8_data 
chd8_targets <- (chd8_data %>% filter(has_CHIP_seq_binding_site == 1))$gene.symbol


#read in network 
network <- read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/full_aracne_network_analysis/Output/All_TFs_And_Epigenetic_Regulators/aracne_output/network_formatted.txt", header = FALSE)
colnames(network) <- c('Regulator', 'Regulatee', 'weight') #note: weight is MI
all_network_genes <- unique(c(network$Regulatee, network$Regulator))


#read in universe
basemean_df <-  read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/baseMeanInNTCs.csv", row.names = 1) #In Input_Files_Not_Generated_By_Scripts
universe <- row.names(basemean_df %>% filter(baseMean > 10))


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


#set up function to read in gmt ontology file
pathway_list <- list()
pathway_list[["SysID"]] <- just_id_genes
pathway_list[["SFARI"]] <- SFARI_genes
pathway_list[["HCSFARI"]] <- high_conf_SFARI_genes
pathway_list[["KMT2E_DEGs"]] <-  all_psDEGs[["KMT2E"]]
pathway_list[["SETD2_DEGs"]] <-  all_psDEGs[["SETD2"]]
pathway_list[["SETD5_DEGs"]] <-  all_psDEGs[["SETD5"]]


##########################
#CREATE GRAPH FROM NETWORK
##########################
graph <- graph.data.frame(network, directed = TRUE)
pathway_list[["CHD8_First_Deg_Regulon"]] <-  network %>% filter(Regulator == "CHD8") %>% pull(Regulatee)
pathway_list[["CHD8_Second_Deg_Regulon"]] <-  unique(c(pathway_list[["CHD8_First_Deg_Regulon"]], unique(do.call(c, lapply(pathway_list[["CHD8_First_Deg_Regulon"]], function(x) network %>% filter(Regulator == x) %>% pull(Regulatee))))))

##########################
#PLOT ENRICHMENT
##########################
#run enrichment
enr_df <- sapply(pathway_list, function(x) generate_overlap_p_val_and_or(x, chd8_targets, universe)) %>% as.data.frame() %>% t() %>% as.data.frame()
colnames(enr_df) <- c("pval", "OR", "lower.ci", "InList1AndList2", "InList1AndNotList2", "InList2AndNotList1", "InNeither")
enr_df$pathway <- row.names(enr_df)
enr_df$padj <- p.adjust(enr_df$pval, method = "BH")
enr_df <- enr_df %>% arrange(desc(OR))
enr_df$pathway <- factor(enr_df$pathway, levels = unique(enr_df$pathway))
enr_df$sig <- sapply(enr_df$padj, function(x) ifelse(x <0.0001, "< 0.0001", ifelse(x < 0.001, "< 0.001", ifelse(x < 0.01, "< 0.01", ifelse(x < 0.05, "< 0.05", "Not Significant")))))
enr_df$sig <- factor(as.character(enr_df$sig), 
                              levels = c("< 0.0001", "< 0.001", "< 0.01", "< 0.05", "Not Significant"))

#write out the matrix values as source data for the reviewer comment (list1 = pathway/gene set, list2 = chd8 targets)
write.table(enr_df, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Actually_Submitted/For_Third_Submission_And_What_We_Published/6_CLEANED_6_Figure_4_And_S9/5_Figure_4F_CHD8_Analysis/CHD8_enr_barplot_source_data.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


#plot
cols <- c(brewer.pal(8, "Greens")[c(8, 6, 4, 2)], "lightgray")
names(cols) <- c("< 0.0001", "< 0.001", "< 0.01", "< 0.05", "Not Significant")

pdf("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Actually_Submitted/For_Third_Submission_And_What_We_Published/6_CLEANED_6_Figure_4_And_S9/5_Figure_4F_CHD8_Analysis/CHD8_enr_barplot.pdf", width = 6, height =3)
ggplot(data = enr_df, aes(x = OR, y = forcats::fct_rev(pathway), fill = sig, width = 0.8, color = "black")) + 
        geom_bar(stat = "identity", color = "black") + 
        theme_minimal() +
        ylab("Gene Set") + 
        xlab("Odds Ratio") +   
geom_errorbar(aes(xmin = lower.ci, xmax = OR),  # only go down to lower CI
                  width = 0.2,                      # horizontal cap width
                  color = "black") +
        theme(axis.text = element_text(size = 7), 
              axis.title = element_text(size = 8, face = "bold")) + 
        scale_fill_manual(values = cols) + 
        guides(fill = guide_legend(title = "Adjusted p-value")) + 
        theme(legend.text = element_text(size = 7), 
              legend.title = element_text(size = 8, face = "bold")) + 
        guides(color = "none") 
dev.off()