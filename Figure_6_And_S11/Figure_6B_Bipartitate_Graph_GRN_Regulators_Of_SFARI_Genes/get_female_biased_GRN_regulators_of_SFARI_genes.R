#################
#I/O
#################
.libPaths(c("/n/groups/walsh/indData/Maya/RLibs", .libPaths()))
library(tidyverse)
library(data.table)
library(pheatmap)
library(RColorBrewer)
library(igraph)
library(readxl)

#generate custom functions
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

	return(c(pval, estimate, fisher.test(mat, alternative = "greater")$conf.int[1]))
}

get_direction_of_regulation <- function(regulator, regulatee, processed_mat_merged)
{
    cor_val <- cor(unlist(processed_mat_merged[regulatee,]), unlist(processed_mat_merged[regulator,]), method = "spearman")
    if (cor_val > 0)
    {
        return ("POS")
    } else 
    {
        return ("NEG")
    }
}

#read in overall data matrix
baseMean_in_NTCs <- read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/baseMeanInNTCs.csv", row.names = 1)
universe <- row.names(baseMean_in_NTCs %>% filter(baseMean >= 10))
processed_mat_merged <- fread("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/full_aracne_network_analysis/Input/matrix.txt", header = TRUE) %>% as.data.frame() %>% filter(gene %in% universe)
row.names(processed_mat_merged) <- processed_mat_merged$gene
processed_mat_merged$gene <- NULL
rownames <- row.names(processed_mat_merged)
processed_mat_merged <- apply(processed_mat_merged, 2, function(x) scale(x))
row.names(processed_mat_merged) <- rownames

#read in network
network <- read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/full_aracne_network_analysis/Output/All_TFs_And_Epigenetic_Regulators/aracne_output/network_formatted.txt", header = FALSE)
colnames(network) <- c('Regulator', 'Regulatee', 'weight') #note: weight is MI
all_network_genes <- unique(c(network$Regulatee, network$Regulator))
universe <- all_network_genes
graph <- graph.data.frame(network, directed = TRUE)
full_network <- network 
full_network$cor <- apply(full_network, 1, function(x) get_direction_of_regulation(x["Regulator"], x["Regulatee"], processed_mat_merged))

#read in constraint scores
constraint_scores <- fread("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/supplementary_dataset_11_full_constraint_metrics.tsv", header = TRUE)
constraint_scores <- as.data.frame(constraint_scores)
constraint_scores <- constraint_scores %>% as.data.frame() %>% dplyr::select(oe_lof_upper, gene) %>% filter(!is.na(oe_lof_upper))
constraint_scores <- constraint_scores %>% group_by(gene) %>% summarize(median_oe_lof_upper_across_transcripts = median(oe_lof_upper))
ploeuf_leq_0p35_genes <- constraint_scores %>% filter(median_oe_lof_upper_across_transcripts < 0.35) %>% pull(gene) %>% as.character()

#read in sfari genes
sfari <- read.csv("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/sfari.csv", header = TRUE)
sfari_genes <- unique(sfari$gene.symbol)
high_conf_SFARI_genes <- unique((sfari %>% filter(gene.score == 1))$gene.symbol)

#read in regulators
all_prioritized_regs <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Revisions/Redo_Figure_6_And_S6_To_Add_Decasien_Data/Figure_6A_S6A_Heatmap_Expression_Candidate_Contributors_FPE/prioritized_regs_consistent_downreg_DEG_in_asd_kd_genes_and_potential_contributors_to_the_fpe_constrained.RDS")
all_prioritized_regs <- intersect(all_prioritized_regs, sfari_genes)
all_prioritized_regs <- intersect(all_prioritized_regs, network$Regulator)

#get regulons
regulon_list <- setNames(
  lapply(all_prioritized_regs, function(x) {
    first_degree_reg <- unique(network %>% filter(Regulator == x) %>% pull(Regulatee))
    second_degree_reg <- unique(network %>% filter(Regulator %in% first_degree_reg) %>% pull(Regulatee))
    unique(c(first_degree_reg, second_degree_reg))
  }),
  all_prioritized_regs
)

#################
#CREATE BARPLOT OF ENRICHMENT
#################
regulon_list_prioritized <- regulon_list[which(names(regulon_list) %in% all_prioritized_regs)]

#run enrichment for SFARI risk genes
regulon_enr_for_sfari_genes <- lapply(regulon_list_prioritized, function(x) generate_overlap_p_val_and_or(x, sfari_genes, universe)) %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame() 
colnames(regulon_enr_for_sfari_genes) <- c("pval", "OR", "lower.ci")
regulon_enr_for_sfari_genes$type <- "SFARI"
regulon_enr_for_sfari_genes$Gene <- row.names(regulon_enr_for_sfari_genes)

#run enrichment for HCSFARI risk genes
regulon_enr_for_hcsfari_genes <- lapply(regulon_list_prioritized, function(x) generate_overlap_p_val_and_or(x, high_conf_SFARI_genes, universe)) %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame()
colnames(regulon_enr_for_hcsfari_genes) <- c("pval", "OR", "lower.ci")
regulon_enr_for_hcsfari_genes$type <- "HCSFARI"
regulon_enr_for_hcsfari_genes$Gene <- row.names(regulon_enr_for_hcsfari_genes)

#combine dataframe 
regulon_enr_for_sfari_genes <- rbind(regulon_enr_for_sfari_genes, regulon_enr_for_hcsfari_genes)
regulon_enr_for_sfari_genes$padj <- p.adjust(regulon_enr_for_sfari_genes$pval, method = "BH")

#add constraint scores
regulon_enr_for_sfari_genes$pLOEUF <- as.numeric(sapply(row.names(regulon_enr_for_sfari_genes), function(x) constraint_scores$median_oe_lof_upper_across_transcripts[which(constraint_scores$gene == x)]))

#################
#CREATE BI PARTITE GRAPH
#################
#outer layer = all_prioritized regs 
outer_nodes <- names(regulon_list)

#inner layer = sfari genes in the regulons 
inner_nodes <- lapply(regulon_list, function(x) intersect(x, sfari_genes))

#create edges 
edges <- stack(inner_nodes)
colnames(edges) <- c("Regulatee", "Regulator")
edges <- edges %>% dplyr::select(Regulator, Regulatee) %>% filter(!Regulator == Regulatee)
inner_nodes <- setdiff(edges$Regulatee, outer_nodes)

#make graph
g <- graph_from_data_frame(edges, directed = TRUE)

#set layer attribute
V(g)$layer <- NA
V(g)$layer[V(g)$name %in% outer_nodes] <- 1
V(g)$layer[V(g)$name %in% inner_nodes] <- 2

#make layout
n1 <- sum(V(g)$layer == 1)
n2 <- sum(V(g)$layer == 2)

theta1 <- sample(seq(0, 2*pi, length.out = n1+1)[-1])
theta2 <- sample(seq(0, 2*pi, length.out = n2+1)[-1])
names(theta1) <- outer_nodes

r1 <- 2
r2 <- 1

layout <- matrix(NA, nrow = vcount(g), ncol = 2)
layout[V(g)$layer == 1, ] <- cbind(r1 * cos(theta1), r1 * sin(theta1))
layout[V(g)$layer == 2, ] <- cbind(r2 * cos(theta2), r2 * sin(theta2))

#set node labels
V(g)$label <- ifelse(V(g)$layer == 1, V(g)$name, NA)

#set edge styles
edge_ends <- ends(g, E(g), names = TRUE)

E(g)$color <- "black"

E(g)$width <- 0.5

E(g)$arrow.mode <- 0

E(g)$lty <- "solid"

#set node sizes
V(g)$size <- ifelse(
  V(g)$layer == 1, 20, 6
)

#set node colors
V(g)$color <- NA
V(g)$color[is.na(V(g)$color) & !V(g)$name %in% sfari_genes] <- "#0071C5"
V(g)$color[is.na(V(g)$color) & V(g)$name %in% high_conf_SFARI_genes] <- "#fddb27ff"
V(g)$color[is.na(V(g)$color) & V(g)$name %in% setdiff(sfari_genes, high_conf_SFARI_genes)] <- "#ee4266"
V(g)$frame.color <-  NA
V(g)$frame.width <- 0

#save plot
pdf("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Revisions/Redo_Figure_6_And_S6_To_Add_Decasien_Data/Figure_6B_Bipartitate_Graph_GRN_Regulators_Of_SFARI_Genes/bipartite_graph.pdf")

plot(
  g,
  layout = layout,
  vertex.size = V(g)$size,
  vertex.label.cex = 0.8,
  vertex.color = V(g)$color,
  edge.arrow.size = 0.3,
  edge.lty = E(g)$lty,
  edge.curved = 0, 
)

#close plot
dev.off()




