##########################
#I/O
##########################
library(viper)
library(data.table)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(UpSetR)
library(grid)
library(igraph)
library(fgsea)
library(lsa)
library(readxl)
library(RColorBrewer)
library(readxl)
library(Orthology.eg.db)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(ggrepel)
library(biomaRt)

#functions
#write function to continually cluster until meeting max.size criteria
perform_iterative_walktrap_clustering <- function(graph, current_name = "initialCluster", max_module_size = 300)
{
  print(paste0("STARTING!"))
  finalized_clusters <- list()
  hasConverged <- FALSE 
  current_graph <- graph


  # Perform walktrap clustering 
  clusterWalktrap_clust <- igraph::cluster_walktrap(current_graph)


  # Determine size of clusters - we have converged for a cluster if it has < 300 genes
  membership_tab <- table(membership(clusterWalktrap_clust)) 
  names(membership_tab) <- paste0(current_name, "_", names(membership_tab))


  # Extract community membership
  communities_list <- communities(clusterWalktrap_clust)
  names(communities_list) <- paste0(current_name, "_", names(communities_list))


  # If a community has already converged in size, add it to the finalized_clusters list
  converged_clusters <- names(membership_tab)[which(membership_tab <= max_module_size)]
  nonconverged_clusters <- names(membership_tab)[which(membership_tab > max_module_size)]
  finalized_clusters <- c(finalized_clusters, communities_list[converged_clusters])


  # If there are still nonconverged clusters, then rerun the algorithm on each of them 
  for (nc in nonconverged_clusters)
  {
    # Create graph of these nodes 
    new_graph <- subgraph(current_graph, unique(do.call(c, communities_list[nc])))
    finalized_clusters <- c(finalized_clusters, perform_iterative_walktrap_clustering(new_graph, current_name))
  }


  return (finalized_clusters)
}


##generate hypergeometric p-value
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


#read in target genes & degs
targets <- unname(unlist(read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/final_targets_to_carry_forward.csv")))
file_path <- "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/DEG_Results/WithManualAdjustment/"
all_DESeq_res <- lapply(targets, function(x) read.table(paste0(file_path, x, "_vs_NTC_withDedup_manual_adjustment_Cutoff_Of_10.txt")))
names(all_DESeq_res) <- targets
all_psDEGs <- lapply(all_DESeq_res, function(x) row.names(x %>% filter(padj < 0.05)))


#read in network 
network <- read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/full_aracne_network_analysis/Output/All_TFs_And_Epigenetic_Regulators/aracne_output/network_formatted.txt", header = FALSE)
colnames(network) <- c('Regulator', 'Regulatee', 'weight') #note: weight is MI
all_network_genes <- unique(c(network$Regulator, network$Regulatee))
all_regulators <- unique(c(network$Regulator))


#get sig diff active regulons when considering all perturbations
mrs <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/full_aracne_network_analysis/Output/All_TFs_And_Epigenetic_Regulators/viper_results/all_perturbations/mrs.RDS")
mrs_sig_regulons <- sapply(mrs$regulon[ names(mrs$es$p.value[mrs$es$p.value< 0.05])], function(x) names(x$tfmode))


#read in SFARI genes 
sfari <- read.csv("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/sfari.csv", header = TRUE)
SFARI_genes <- unique(sfari$gene.symbol)
high_conf_SFARI_genes <- unique((sfari %>% filter(gene.score == 1))$gene.symbol)


#read in ID genes 
id_genes <- read.csv("/n/groups/walsh/indData/becky/SysNDD_figures_250823/Input/SysNDD_ID_genes_definitive.csv", header = TRUE)
just_id_genes <- setdiff(id_genes$symbol, SFARI_genes)
just_SFARI_genes <- setdiff(SFARI_genes, id_genes)
just_high_conf_SFARI_genes <- setdiff(high_conf_SFARI_genes, id_genes)


#set up function to read in gmt ontology file
#set up lists that will be used for fgsea 
pathway_list <- list()
pathway_list[["sfari"]] <- just_SFARI_genes
pathway_list[["sysID"]] <- just_id_genes
pathway_list[["full_sfari"]] <- SFARI_genes
pathway_list[["full_high_conf_sfari"]] <- high_conf_SFARI_genes
pathway_list[["high_conf_sfari"]] <- just_high_conf_SFARI_genes
pathway_list <- c(pathway_list, gmtPathways("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/full_aracne_network_analysis/Input/Ontologies/h.all.v2023.2.Hs.symbols.gmt"))


#read in expression matrix
dset <- as.data.frame(fread("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/full_aracne_network_analysis/Input/matrix.txt", header = TRUE))
dset <- dset[-which(is.na(dset$gene)),,drop = FALSE]
row.names(dset) <- dset$gene
dset$gene <- NULL
dset <- as.matrix(dset)


#read in network 
network <- read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/full_aracne_network_analysis/Output/All_TFs_And_Epigenetic_Regulators/aracne_output/network_formatted.txt", header = FALSE)
colnames(network) <- c('Regulator', 'Regulatee', 'weight') #note: weight is MI
all_network_genes <- unique(c(network$Regulatee, network$Regulator))


#read in universe
basemean_df <-  read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/baseMeanInNTCs.csv", row.names = 1)
universe <- row.names(basemean_df %>% filter(baseMean > 10))


#read in all NPC expressed TFs and epigenetic regulators
reg_universe <- unique(network$Regulator)


##########################
#CREATE GRAPH FROM NETWORK
##########################
graph <- graph.data.frame(network, directed = TRUE)


# ##########################
# #NOMINATE REGULATORS 
# ##########################
#criteria #1: sig diff active between controls and just samples in which genes associated with ASD had been perturbed
mrs <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/full_aracne_network_analysis/Output/All_TFs_And_Epigenetic_Regulators/viper_results/just_genes_associated_with_asd/mrs.RDS")
dars <- names(mrs$es$p.value)[which(mrs$es$p.value < 0.05)]


#criteria #2: regulon enriched for high pagerank scores
##run pagerank 
pagerank_res <- page_rank(graph)
pagerank_df <- as.data.frame(pagerank_res$vector)
colnames(pagerank_df) <- c("score")
pagerank_df$gene <- row.names(pagerank_df)
pagerank_df <- pagerank_df %>% arrange(desc(score)) 
pagerank_df$percentile <- rank(pagerank_df$score)/length(pagerank_df$score)
pagerank_df_pct_vec <- setNames(pagerank_df$percentile, pagerank_df$gene)


##run enrichment of regulons for high page rank scores 
regulon_list <- setNames(lapply(dars, function(x) (network %>% filter(Regulator == x))$Regulatee), dars)
ranks <- setNames(pagerank_df$score, pagerank_df$gene)
ranks <- ranks[-which(names(ranks) %in% targets)]
fgseaRes_pagerank <- fgsea(pathways = regulon_list, 
                  stats    = ranks,
                  minSize  = 50,
                  maxSize  = 5000,
                  scoreType = "pos")


fgseaRes_pagerank <- fgseaRes_pagerank %>% arrange(padj)
pagerank_enriched <- (fgseaRes_pagerank %>% filter(padj < 0.05))$pathway


#criteria #3: module enriched for genes associated with asd
reg_enrichment <- sapply(pagerank_enriched, function(x) generate_overlap_p_val((network %>% filter(Regulator == x))$Regulatee, pathway_list[["full_sfari"]], all_network_genes))

#final candidates
final_candidates <- names(reg_enrichment[which(reg_enrichment < 0.05)])

#run enrichment for sfari genes 
print("*****SFARI*****")
sfari_p_val <- generate_overlap_p_val(pathway_list[["full_sfari"]], final_candidates, reg_universe)
print(generate_overlap_p_val(pathway_list[["full_sfari"]], final_candidates, reg_universe))
print(length(intersect(pathway_list[["full_sfari"]], final_candidates)))
print("*****")


print("*****HCSFARI*****")
hcsfari_pval <- generate_overlap_p_val(pathway_list[["full_high_conf_sfari"]], final_candidates, reg_universe)
print(generate_overlap_p_val(pathway_list[["full_high_conf_sfari"]], final_candidates, reg_universe))
print(length(intersect(pathway_list[["full_high_conf_sfari"]], final_candidates)))
print("*****")

print("*****SYSID*****")
sysid_pval <- generate_overlap_p_val(pathway_list[["sysID"]], final_candidates, reg_universe)
print(generate_overlap_p_val(pathway_list[["sysID"]], final_candidates, reg_universe))
print(length(intersect(pathway_list[["sysID"]], final_candidates)))
print("*****")


print("*****ADJUSTED PVALS FOR SFARI, SYSID, AND HCSFARI*****")
p.adjust(c(sfari_p_val, hcsfari_pval, sysid_pval), method = "BH")

#also run enrichment of final candidate regulons against sysid
final_candidates <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/full_aracne_network_analysis/Output/Analysis_Results/GenePrioritization/final_candidates.rds")
final_cand_regulon_enrichment_sfari <- sapply(final_candidates, function(x) generate_overlap_p_val((network %>% filter(Regulator == x))$Regulatee, pathway_list[["full_sfari"]], all_network_genes))
final_cand_regulon_enrichment_sysID <- sapply(final_candidates, function(x) generate_overlap_p_val((network %>% filter(Regulator == x))$Regulatee, pathway_list[["sysID"]], all_network_genes))
final_cand_regulon_enrichment_df <- cbind(final_cand_regulon_enrichment_sfari, final_cand_regulon_enrichment_sysID)
#final_cand_regulon_enrichment_df <- apply(final_cand_regulon_enrichment_df, 2, function(x) p.adjust(x, method = "BH"))

##########################
#CREATE VISUALIZATION OF THESE REGULATORS (JUST BETWEEN CANDIDATES)
##########################
final_candidates <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/full_aracne_network_analysis/Output/Analysis_Results/GenePrioritization/final_candidates.rds")

#question #1: what is the average degree of various graphs?
med_degree <- list()


##just between our regulators
graph <- graph.data.frame(network %>% filter(Regulator %in% final_candidates))
med_degree[["just_btwn_candidate_regs"]] <- mean(degree(graph, normalized = TRUE))


##just between high conf sfari genes
graph <- graph.data.frame(network %>% filter(Regulator %in% high_conf_SFARI_genes))
med_degree[["just_btwn_high_conf_sfari_regs"]] <- mean(degree(graph, normalized = TRUE))


##just between sfari genes 
graph <- graph.data.frame(network %>% filter(Regulator %in% SFARI_genes))
med_degree[["just_btwn_sfari_regs"]] <- mean(degree(graph, normalized = TRUE))


##just between sys-id
graph <- graph.data.frame(network %>% filter(Regulator %in% just_id_genes))
med_degree[["just_btwn_sysid"]] <- mean(degree(graph, normalized = TRUE))


##just between all regulators 
graph <- graph.data.frame(network)
med_degree[["just_btwn_all_regs"]] <- mean(degree(graph, normalized = TRUE))


#question 2: is the high connectivity between our 70 candidate genes significant?
all_regulators <- unique(network$Regulator)
run_bootstrap_network <- function(network, num_regs) 
{
  #randomly select regulators
  selected_regs <- sample(unique(network$Regulator), num_regs, replace = FALSE)


  #subset network
  graph <- graph.data.frame(network %>% filter(Regulator %in% selected_regs))


  return(mean(degree(graph, normalized = TRUE)))
}


bootstrap_TopRegs <- sapply(seq(1:1000000), function(x) run_bootstrap_network(network, length(intersect(final_candidates,all_regulators))))
table(med_degree[["just_btwn_candidate_regs"]] > bootstrap_TopRegs)
saveRDS(bootstrap_TopRegs, "/n/groups/walsh/indData/becky/SysNDD_figures_250823/Fig5d/bootstrap_TopRegs.RDS")


bootstrap_HCSFARI <- sapply(seq(1:1000000), function(x) run_bootstrap_network(network, length(intersect(high_conf_SFARI_genes,all_regulators))))
table(med_degree[["just_btwn_high_conf_sfari_regs"]] > bootstrap_HCSFARI)
saveRDS(bootstrap_HCSFARI, "/n/groups/walsh/indData/becky/SysNDD_figures_250823/Fig5d/bootstrap_HCSFARI.RDS")


bootstrap_SFARI <- sapply(seq(1:1000000), function(x) run_bootstrap_network(network, length(intersect(SFARI_genes,all_regulators))))
table(med_degree[["just_btwn_sfari_regs"]] > bootstrap_SFARI)
saveRDS(bootstrap_SFARI, "/n/groups/walsh/indData/becky/SysNDD_figures_250823/Fig5d/bootstrap_SFARI.RDS")


bootstrap_sysID <- sapply(seq(1:1000000), function(x) run_bootstrap_network(network, length(intersect(just_id_genes,all_regulators))))
table(med_degree[["just_btwn_sysid"]] > bootstrap_sysID)
saveRDS(bootstrap_sysID, "/n/groups/walsh/indData/becky/SysNDD_figures_250823/Fig5d/bootstrap_sysID.RDS")


#calculate relative med degree using mean bootstrap values
mean_degree_results<-as.data.frame(c("Top Drivers", "HC-SFARI Regulators", "SFARI Regulators", "Only ID Regulators"))
colnames(mean_degree_results)<-"gene_set"
mean_degree_results$rel_mean_norm_degree<-c(med_degree[["just_btwn_candidate_regs"]]/mean(bootstrap_TopRegs),
                                            med_degree[["just_btwn_high_conf_sfari_regs"]]/mean(bootstrap_HCSFARI),
                                            med_degree[["just_btwn_sfari_regs"]]/mean(bootstrap_SFARI),
                                            med_degree[["just_btwn_sysid"]]/mean(bootstrap_sysID))
mean_degree_results$pval<-c(length(which(bootstrap_TopRegs >= med_degree[["just_btwn_candidate_regs"]]))/length(bootstrap_TopRegs),
                        length(which(bootstrap_HCSFARI >= med_degree[["just_btwn_high_conf_sfari_regs"]]))/length(bootstrap_HCSFARI),
                        length(which(bootstrap_SFARI >= med_degree[["just_btwn_sfari_regs"]]))/length(bootstrap_SFARI),
                        length(which(bootstrap_sysID >= med_degree[["just_btwn_sysid"]]))/length(bootstrap_sysID))

##plot
mean_degree_results$Category <- row.names(mean_degree_results)
mean_degree_results <- mean_degree_results %>% arrange(mean_degree_results)
mean_degree_results$Category <- factor(mean_degree_results$Category, levels = as.character(mean_degree_results$Category))
pdf("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Figure_5_And_S5/Figure_5DE_S5B/RelMeanNormDegree.pdf", width = 4, height = 4)
ggplot(mean_degree_results, aes(x = rel_mean_norm_degree, y = Category)) +
  geom_bar(stat = "identity") +
  labs(y = "Category",
       x = "Rel. Mean Normalized Degree") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))
dev.off()


#visualize with community detection 
graph <- graph.data.frame(network %>% filter(Regulator %in% final_candidates) %>% filter(Regulatee %in% final_candidates))
V(graph)$color <- ifelse(V(graph)$name %in% pathway_list[["full_sfari"]], "goldenrod1", "cadetblue1")
finalized_walktrap_clusters_full_res <- igraph::cluster_walktrap(graph, steps = 8)
finalized_walktrap_clusters <- communities(finalized_walktrap_clusters_full_res)
names(finalized_walktrap_clusters) <- paste0("module_", seq_along(finalized_walktrap_clusters))

pdf("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Figure_5_And_S5/Figure_5DE_S5B/just_candidate_genes_with_comms.pdf", width = 12, height = 12)
plot(finalized_walktrap_clusters_full_res, graph, 
     layout = layout_with_fr(graph), 
     vertex.size = 5, 
     vertex.label.cex = 0.7, 
     vertex.label.color = "black", 
     vertex.label.font = 2, # bold text
     edge.arrow.size = 0.3, # small arrows
     edge.width = 0.5, # thinner edges
     edge.color = adjustcolor("gray", alpha.f = 0.3)) # translucent pale gray arrows
dev.off()


lapply(names(finalized_walktrap_clusters), function(x) write.table(finalized_walktrap_clusters[[x]], paste0("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/full_aracne_network_analysis/Output/Analysis_Results/Clustering/finalized_walktrap_clusters_", x, ".txt"), sep = "\t", row.names = F, col.names = F, quote = F))


##########################
#CREATE VISUALIZATION OF THESE REGULATORS IN PERTURBATION DATA (HEATMAP)
##########################
final_candidates <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/full_aracne_network_analysis/Output/Analysis_Results/GenePrioritization/final_candidates.rds")
regs_for_figures <- final_candidates

#read in perturbation data
file_path <- "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/DEG_Results/WithManualAdjustment/"
sgRNAs <- unname(unlist(read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/final_targets_to_carry_forward.csv")))
print(length(sgRNAs) == 18)
all_DESeq_res <- lapply(sgRNAs, function(x) read.table(paste0(file_path, x, "_vs_NTC_withDedup_manual_adjustment_Cutoff_Of_10.txt")))
names(all_DESeq_res) <- sgRNAs

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

##perturbation analysis
processed_mat_merged <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/processed_mat_merged.RDS")
baseMean_in_NTCs <- read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/baseMeanInNTCs.csv", row.names = 1)
genes_to_keep <- row.names(baseMean_in_NTCs %>% filter(baseMean >= 10))
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
  filename = "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Figure_5_And_S5/Figure_5DE_S5B/heatmapVis.pdf",
  width = 15, height = 6
)


##########################
#MAKE SUPP TABLE WITH FINAL CANDIDATE INFORMATION
##########################
final_candidate_df <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/full_aracne_network_analysis/Output/Analysis_Results/GenePrioritization/final_candidates.rds") %>% as.data.frame()
colnames(final_candidate_df) <- "Gene"
candidate_autism_risk_genes <- setdiff(readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/figures_for_manuscript/get_sex_assoc_regs/0_get_consistent_ASD_DEGs/candidate_asd_risk_genes.RDS"), SFARI_genes)
final_candidate_df$Category <- sapply(final_candidate_df$Gene, function(x) ifelse(x %in% high_conf_SFARI_genes, "HC-SFARI", ifelse(x %in% SFARI_genes, "SFARI", ifelse(x %in% candidate_autism_risk_genes, "Candidate Autism Risk Gene", "Other"))))
final_candidate_df <- final_candidate_df %>% arrange(Gene)
write.table(final_candidate_df, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Figure_5_And_S5/Figure_5DE_S5B_GRN_Properties_Visualization_Top_Candidates_Analysis/top_candidate_info_df.txt", sep = "\t", row.names = F, col.names = T, quote = F)