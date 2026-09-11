##########################
#I/O
##########################
.libPaths(c("/n/groups/walsh/indData/Maya/RLibs", .libPaths()))
library(data.table)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(fgsea)
library(igraph)


targets <- unname(unlist(read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/final_targets_to_carry_forward.csv")))


#read in network 
network <- read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/full_aracne_network_analysis/Output/All_TFs_And_Epigenetic_Regulators/aracne_output/network_formatted.txt", header = FALSE)
colnames(network) <- c('Regulator', 'Regulatee', 'weight') #note: weight is MI
all_network_genes <- unique(c(network$Regulator, network$Regulatee))
all_regulators <- unique(network$Regulator)


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


#set up function to read in gmt ontology file
#set up lists that will be used for fgsea 
pathway_list <- list()
pathway_list[["SysID"]] <- just_id_genes
pathway_list[["SFARI"]] <- SFARI_genes
pathway_list[["HCSFARI"]] <- high_conf_SFARI_genes
pathway_list <- c(pathway_list, gmtPathways("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/full_aracne_network_analysis/Input/Ontologies/h.all.v2023.2.Hs.symbols.gmt"))


##########################
#CREATE GRAPH FROM NETWORK
##########################
graph <- graph.data.frame(network, directed = TRUE)


##########################
#RUN PAGE-RANK ALGORITHM
##########################
#run pagerank 
pagerank_res <- page_rank(graph)
pagerank_df <- as.data.frame(pagerank_res$vector)
colnames(pagerank_df) <- c("score")
pagerank_df$gene <- row.names(pagerank_df)
pagerank_df <- pagerank_df %>% arrange(desc(score))
pagerank_df$percentile <- rank(pagerank_df$score)/length(pagerank_df$score)
pagerank_df_pct_vec <- setNames(pagerank_df$percentile, pagerank_df$gene)


#check for enrichment of genes ranked by pagerank score against ontologies of interest
ranks <- setNames(pagerank_df$score, pagerank_df$gene)
ranks <- ranks[-which(names(ranks) %in% targets)]
fgseaRes_pagerank <- fgsea(pathways = pathway_list, 
                  stats    = ranks,
                  minSize  = 20,
                  maxSize  = 5000,
                  scoreType = "pos")
fgseaRes_pagerank %>% filter(padj < 0.05)
fgseaRes_pagerank <- fgseaRes_pagerank %>% arrange(padj)
fgseaRes_pagerank$is_sig <- ifelse(fgseaRes_pagerank$padj < 0.05, "yes", "no")
fgseaRes_pagerank$neg_log_10_padj <- -1 * log10(fgseaRes_pagerank$padj)
fgseaRes_pagerank$cleaned_pathway <- sapply(fgseaRes_pagerank$pathway, function(x) gsub("_", " ", stringr::str_to_title(gsub("HALLMARK_", "", x))))
fgseaRes_pagerank$cleaned_pathway <- factor(fgseaRes_pagerank$cleaned_pathway, levels = fgseaRes_pagerank$cleaned_pathway)


#plot
dir.create("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Initial_Submission/Figure_5_And_S5/Figure_5C_PageRank_GSEA/")
pdf("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Initial_Submission/Figure_5_And_S5/Figure_5C_PageRank_GSEA/pagerank_gsea.pdf", width = 5, height = 10)
ggplot(data = fgseaRes_pagerank, aes(x = neg_log_10_padj, y = cleaned_pathway)) +
  geom_point() + 
  theme_minimal() + 
  xlab("-log10 Adj. P-Value") + 
  ylab("Pathway") + geom_vline(xintercept =  1.30103, color = "red")
dev.off()


pdf("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Initial_Submission/Figure_5_And_S5/Figure_5C_PageRank_GSEA/pagerank_gsea_sig.pdf", width = 3, height = 2)
ggplot(data = fgseaRes_pagerank %>% filter(padj < 0.05), aes(x = neg_log_10_padj, y = cleaned_pathway)) +
  geom_point() + 
  theme_minimal() + 
  xlab("-log10 Adj. P-Value") + 
  ylab("Pathway") + geom_vline(xintercept =  1.30103, color = "red") + xlim(0, 9)
dev.off()
write.table(pagerank_df, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Initial_Submission/Figure_5_And_S5/Figure_5C_PageRank_GSEA/pagerank_res.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
saveRDS(fgseaRes_pagerank, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Initial_Submission/Figure_5_And_S5/Figure_5C_PageRank_GSEA/pagerank_GSEA_res.RDS")


write.table(fgseaRes_pagerank %>% dplyr::select(pathway, NES, size, pval, padj), "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Initial_Submission/Figure_5_And_S5/Figure_5C_PageRank_GSEA/full_gsea_results_for_3d.txt", sep= "\t", quote = F, row.names = F, col.names = T)

