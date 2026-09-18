##############################
#I/O
##############################
.libPaths(c("/n/groups/walsh/indData/Maya/RLibs", .libPaths()))
library(tidyverse)
library(data.table)
library(pheatmap)
library(RColorBrewer)

#read in crispr kd data 
file_path <- "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/DEG_Results/WithManualAdjustment/"
sgRNAs <- unname(unlist(read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/final_targets_to_carry_forward.csv"))) #In Input_Files_Not_Generated_By_Scripts
all_res_list <- lapply(sgRNAs, function(x) read.table(paste0(file_path, x, "_vs_NTC_withDedup_manual_adjustment_Cutoff_Of_10.txt")))
names(all_res_list) <- sgRNAs
all_DEGs_list_npc <- lapply(all_res_list, function(x) row.names((x %>% filter(padj < 0.05) %>% filter(log2FoldChange < 0))))

##define the universe as all genes expressed with baseMean >= 10 
basemean_df <-  read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/baseMeanInNTCs.csv", row.names = 1)
universe_npc <- row.names(basemean_df %>% filter(baseMean > 10))
npc_genes <- sgRNAs

#read in shRNA kd data
path <- "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Revisions/Analysis_of_shRNA_Data/0_Process_And_Run_Differential_Expression/differential_expression_analysis/DEG_Results/WithManualAdjustment/"
all_DESeq_res_list <- lapply(list.files(path, pattern = ".txt"), function(x) read.table(paste0(path, x), header = TRUE))
names(all_DESeq_res_list) <- sapply(list.files(path, pattern = ".txt"), function(x) gsub("_vs_NTC_manual_adjustment_Cutoff_Of_10.txt", "", x))
all_DEGs_list_shRNA <- lapply(all_DESeq_res_list, function(x) row.names(x %>% filter(padj < 0.05) %>% filter(log2FoldChange < 0)))
shRNA_genes <- paste0(sapply(names(all_DEGs_list_shRNA), function(x) gsub("_vs_ctrl_manual_adjustment_Cutoff_Of_10.txt", "", x)),"_shRNA")

##define the universe as all genes expressed with baseMean >= 10 
baseMean_df <- read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Revisions/Analysis_of_shRNA_Data/0_Process_And_Run_Differential_Expression/differential_expression_analysis/DEG_Results/OriginalResultsBeforeManualAdjustment/BAZ2B_vs_ctrl.txt", header = TRUE) %>% dplyr::select(baseMean)
universe_shRNA <- row.names(baseMean_df %>% filter(baseMean > 10))

length(intersect(universe_shRNA, universe_npc))

##########################
#HELPER FUNCTIONS
##########################
#run enrichment test 
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

##########################
#iPSCs
##########################
#get DEGs per perturbation
wd <- "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_denault_perturbation_data/differential_expression_analysis/deg_results/ipscs/WithManualAdjustment/"
all_res_list <- lapply(list.files(wd), function(x) read.table(paste0(wd, x), header = TRUE))
genes <- sapply(list.files(wd), function(x) gsub("_vs_control_withDedup_manual_adjustment_Cutoff_Of_10.txt", "", x))
all_DEGs_list_ipscs <- lapply(all_res_list, function(x) row.names((x %>% filter(padj < 0.05) %>% filter(log2FoldChange < 0))))
basemean <- read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_denault_perturbation_data/differential_expression_analysis/deg_results/ipscs/baseMean.csv", sep = ",", row.names = 1)
universe_ipscs <- row.names(basemean %>% filter(V2 >= 10))
ipsc_genes <- unname(sapply(genes, function(x) gsub("_vs_control_DEG_res_manual_adjustment_Cutoff_Of_10.txt", "", x)))

##########################
#Neurons
##########################
#get DEGs per perturbation
wd <- "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_denault_perturbation_data/differential_expression_analysis/deg_results/neurons/WithManualAdjustment/"
all_res_list <- lapply(list.files(wd), function(x) read.table(paste0(wd, x), header = TRUE))
genes <- sapply(list.files(wd), function(x) gsub("_vs_control_withDedup_manual_adjustment_Cutoff_Of_10.txt", "", x))
all_DEGs_list_neurons <- lapply(all_res_list, function(x) row.names((x %>% filter(padj < 0.05) %>% filter(log2FoldChange < 0))))
basemean <- read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_denault_perturbation_data/differential_expression_analysis/deg_results/neurons/baseMean.csv", sep = ",", row.names = 1)
universe_neurons <- row.names(basemean %>% filter(V2 >= 10))
neuron_genes <- unname(sapply(genes, function(x) gsub("_vs_control_DEG_res_manual_adjustment_Cutoff_Of_10.txt", "", x)))

##########################
#Cross-enrichment
##########################
all_DEGs_list <- c(all_DEGs_list_ipscs, all_DEGs_list_neurons, all_DEGs_list_npc, all_DEGs_list_shRNA)
names(all_DEGs_list) <- c(sapply(ipsc_genes, function(x) paste0(x, "_ipsc")), sapply(neuron_genes, function(x) paste0(x, "_neuron")),sapply(npc_genes, function(x) paste0(x, "_npc")), shRNA_genes)
universe_labels <- c(rep("ipsc", length(all_DEGs_list_ipscs)), rep("neuron", length(all_DEGs_list_neurons)), rep("npc", length(all_DEGs_list_npc)), rep("shRNA", length(all_DEGs_list_shRNA)))
universe_map <- list("ipsc" = universe_ipscs, "neuron" = universe_neurons, "npc" = universe_npc, "shRNA" = universe_shRNA)
cross_enrichment_res <- sapply(seq_along(all_DEGs_list), function(i)
{
	sapply(seq_along(all_DEGs_list), function(j)
	{
		x <- all_DEGs_list[[i]]
		y <- all_DEGs_list[[j]]

		#get relevant universes
		Ui <- universe_map[[universe_labels[i]]]
		Uj <- universe_map[[universe_labels[j]]]

		#pairwise intersected universe
		joint_universe_ij <- intersect(Ui, Uj)

		length(intersect(x,y))

		generate_overlap_p_val(x, y, joint_universe_ij)

	})
})
row.names(cross_enrichment_res) <- c(sapply(ipsc_genes, function(x) paste0(x, "_ipsc")), sapply(neuron_genes, function(x) paste0(x, "_neuron")),sapply(npc_genes, function(x) paste0(x, "_npc")), shRNA_genes)
colnames(cross_enrichment_res) <- c(sapply(ipsc_genes, function(x) paste0(x, "_ipsc")), sapply(neuron_genes, function(x) paste0(x, "_neuron")), sapply(npc_genes, function(x) paste0(x, "_npc")), shRNA_genes)
write.table(cross_enrichment_res, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Revisions/Analysis_of_Deneault_Data/Output/Just_Downregulated_DEGs/raw_cross_enr_res.txt", row.names = TRUE, sep = "\t", col.names = TRUE,quote = F)

#correct p-values 
cross_enrichment_res <- apply(cross_enrichment_res, 2, function(x) p.adjust(x, method = "BH"))
row.names(cross_enrichment_res) <- c(sapply(ipsc_genes, function(x) paste0(x, "_ipsc")), sapply(neuron_genes, function(x) paste0(x, "_neuron")),sapply(npc_genes, function(x) paste0(x, "_npc")),shRNA_genes)
colnames(cross_enrichment_res) <- c(sapply(ipsc_genes, function(x) paste0(x, "_ipsc")), sapply(neuron_genes, function(x) paste0(x, "_neuron")), sapply(npc_genes, function(x) paste0(x, "_npc")),shRNA_genes)
cross_enrichment_res[lower.tri(cross_enrichment_res)] = t(cross_enrichment_res)[lower.tri(cross_enrichment_res)]
write.table(cross_enrichment_res, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Revisions/Analysis_of_Deneault_Data/Output/Just_Downregulated_DEGs/cross_enrichment_res.txt", sep = "\t", row.names= TRUE,col.names= TRUE,quote= F)

##########################
#Create plot
##########################
cross_enrichment_res <- read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Revisions/Analysis_of_Deneault_Data/Output/Just_Downregulated_DEGs/cross_enrichment_res.txt")

breaks <- c(0, 0.0001, 0.001, 0.01, 0.05, 1)
colors <- c(brewer.pal(8, "Greens")[c(8, 6, 4, 2)], "lightgray")
row_annotations <- setNames(sapply(row.names(cross_enrichment_res), function(x) strsplit(x, "_")[[1]][2]), row.names(cross_enrichment_res))
row_annotations_df <- as.data.frame(row_annotations)
colnames(row_annotations_df) <- "Cell Type"
row.names(row_annotations_df) <- names(row_annotations)
col_annotations <- setNames(sapply(colnames(cross_enrichment_res), function(x) strsplit(x, "_")[[1]][2]), colnames(cross_enrichment_res))
col_annotations_df <- as.data.frame(col_annotations)
colnames(col_annotations_df) <- "Cell Type"
row.names(col_annotations_df) <- names(col_annotations)

pheatmap(cross_enrichment_res, na_col = "white", annotation_row = row_annotations_df, annotation_col = col_annotations_df, cluster_rows = TRUE, cluster_cols = TRUE, color = colors, breaks = breaks, number_color = "black", border_color = "black", legend = FALSE, filename =  "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Revisions/Analysis_of_Deneault_Data/Output/Just_Downregulated_DEGs/cross_enrichment_plot_with_clustering.pdf", width = 10, height = 8)
pheatmap(cross_enrichment_res, na_col = "white", cluster_cols = FALSE,cluster_rows = FALSE,annotation_row = row_annotations_df, annotation_col = col_annotations_df, color = colors, breaks = breaks, number_color = "black", border_color = "black", legend = FALSE, filename =  "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Revisions/Analysis_of_Deneault_Data/Output/Just_Downregulated_DEGs/cross_enrichment_plot_without_clustering.pdf", width = 10, height = 8)
ph <- pheatmap(cross_enrichment_res, na_col = "white", annotation_row = row_annotations_df, annotation_col = col_annotations_df, cluster_rows = TRUE, cluster_cols = TRUE, color = colors, breaks = breaks, number_color = "black", border_color = "black", legend = FALSE, filename =  "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Revisions/Analysis_of_Deneault_Data/Output/Just_Downregulated_DEGs/cross_enrichment_plot_with_clustering.pdf", width = 10, height = 8)
row_order <- ph$tree_row$order
col_order <- ph$tree_col$order
cross_enrichment_res_clustered <- cross_enrichment_res[row_order, col_order]
write.table(cross_enrichment_res_clustered, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Revisions/Analysis_of_Deneault_Data/Output/Just_Downregulated_DEGs/cross_enrichment_res_clustered.txt", sep = "\t", row.names= TRUE,col.names= TRUE,quote= F)

##########################
#Examine cross enrichment with SFARI and ID genes
##########################
#read in SFARI genes 
sfari <- read.csv("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/sfari.csv", header = TRUE) #In Input_Files_Not_Generated_By_Scripts
SFARI_genes <- unique(sfari$gene.symbol)
high_conf_SFARI_genes <- unique((sfari %>% filter(gene.score == 1))$gene.symbol)

#read in ID genes 
id_genes <- read.csv("/n/groups/walsh/indData/becky/SysNDD_figures_250823/Input/SysNDD_ID_genes_definitive.csv", header = TRUE) #In Input_Files_Not_Generated_By_Scripts
just_id_genes <- setdiff(id_genes$symbol, SFARI_genes)

#get enrichment pvalues 
cross_enrichment_res_sfari <- setNames(lapply(seq_along(all_DEGs_list), function(i) generate_overlap_p_val(all_DEGs_list[[i]], SFARI_genes,  universe_map[[universe_labels[i]]])), names(all_DEGs_list)) %>% as.data.frame() %>% t() %>% as.data.frame()
cross_enrichment_res_sysID <- setNames(lapply(seq_along(all_DEGs_list), function(i) generate_overlap_p_val(all_DEGs_list[[i]], just_id_genes,  universe_map[[universe_labels[i]]])), names(all_DEGs_list)) %>% as.data.frame() %>% t() %>% as.data.frame()
cross_enrichment_res_sfari_sysID <- cbind(cross_enrichment_res_sfari , cross_enrichment_res_sysID)
colnames(cross_enrichment_res_sfari_sysID) <- c("SFARI_p", "SysID_p")
cross_enrichment_res_sfari_sysID %>% as.data.frame() %>% arrange(SysID_p)
cross_enrichment_res_sfari_sysID_adjusted <- apply(cross_enrichment_res_sfari_sysID, 2, function(x) p.adjust(x, method = "BH"))
colnames(cross_enrichment_res_sfari_sysID_adjusted) <- c("SFARI_padj", "SysID_padj")
cross_enrichment_res_sfari_sysID_adjusted <- cross_enrichment_res_sfari_sysID_adjusted %>% as.data.frame() %>% mutate(is_SFARI_sig = SFARI_padj < 0.05) %>% mutate(is_SysID_sig = SysID_padj < 0.05) 
cross_enrichment_res_sfari_sysID_adjusted$class <- sapply(row.names(cross_enrichment_res_sfari_sysID_adjusted), function(x) strsplit(x, "_")[[1]][2])
table(cross_enrichment_res_sfari_sysID_adjusted$class, cross_enrichment_res_sfari_sysID_adjusted$is_SFARI_sig)
table(cross_enrichment_res_sfari_sysID_adjusted$class, cross_enrichment_res_sfari_sysID_adjusted$is_SysID_sig)