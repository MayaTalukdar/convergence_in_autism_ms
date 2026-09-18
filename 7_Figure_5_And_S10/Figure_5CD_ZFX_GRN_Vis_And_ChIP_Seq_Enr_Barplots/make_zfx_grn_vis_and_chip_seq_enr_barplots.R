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

#add helper functions 
generate_overlap_p_val <- function(list1, list2, universe)
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


	return(pval)
}


read_excel_allsheets <- function(filename, tibble = FALSE) {
    # I prefer straight data.frames
    # but if you like tidyverse tibbles (the default with read_excel)
    # then just pass tibble = TRUE
    sheets <- readxl::excel_sheets(filename)
    x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
    if(!tibble) x <- lapply(x, as.data.frame)
    names(x) <- sheets
    x
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


	#also return the matrix values so they can be reported as source data for reviewers
	return(c(pval, estimate, fisher.test(mat, alternative = "greater")$conf.int[1], inList1AndList2, inList1AndNotList2, inList2AndNotList1, inNeither))
}


#read in network
network <- read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/full_aracne_network_analysis/Output/All_TFs_And_Epigenetic_Regulators/aracne_output/network_formatted.txt", header = FALSE)
colnames(network) <- c('Regulator', 'Regulatee', 'weight') #note: weight is MI


#read in prioritized regulators
final_candidates <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/full_aracne_network_analysis/Output/Analysis_Results/GenePrioritization/final_candidates.rds")


#get regulon list
regulon_list <- setNames(
  lapply(c(final_candidates, "ZFX"), function(x) {
    first_degree_reg <- unique(network %>% filter(Regulator == x) %>% pull(Regulatee))
    second_degree_reg <- unique(network %>% filter(Regulator %in% first_degree_reg) %>% pull(Regulatee))
    unique(c(first_degree_reg, second_degree_reg))
  }),
  c(final_candidates, "ZFX")
)


#read in data from rhie et al., 2018 
##chip-seq - data obtained from san roman et al., 2024 
rhie_path <- "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/zfx_zfy_data/Rhie_Genome_Res_2018/" #In Input_Files_Not_Generated_By_Scripts
rhie_chip_seq_data <- lapply(list.files(rhie_path), function(x) read.table(paste0(rhie_path, x), header = TRUE))
names(rhie_chip_seq_data) <- sapply(list.files(rhie_path), function(x) gsub(".txt", "", gsub("ZFX_direct_targets_", "", x)))
rhie_targets <- lapply(rhie_chip_seq_data, function(x) (x %>% filter(padj < 0.05))$Gene)

##chip-seq - common targets across all cell types as defined by rhie 2018
common_zfx_targets <- unname(unlist(read.csv("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/figures_for_manuscript/get_sex_assoc_regs/4_perform_zfx_analysis/raw_data/targets_rhie_table_s3.csv", header = FALSE))) #In Input_Files_Not_Generated_By_Scripts

##whole rna-seq results (for defining universe)
rhie_targets_universe <- list()
rhie_targets_universe[["C42B"]] <- read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/figures_for_manuscript/get_sex_assoc_regs/4_perform_zfx_analysis/raw_data/rhie_2018_C42B_siZFX_experiment_gene_counts.txt", header = TRUE) %>%
  column_to_rownames("Gene") %>%
  filter_all(all_vars(!is.na(.))) %>%
  rownames() #In Input_Files_Not_Generated_By_Scripts
rhie_targets_universe[["MCF7"]] <- read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/figures_for_manuscript/get_sex_assoc_regs/4_perform_zfx_analysis/raw_data/rhie_2018_MCF7_siZFX_experiment_gene_counts.txt", header = TRUE) %>%
  column_to_rownames("Gene") %>%
  filter_all(all_vars(!is.na(.))) %>%
  rownames() #In Input_Files_Not_Generated_By_Scripts

##deg results from rna-seq on siRNA kd (includes all DEGs and just DEGs that also have ChIP evidence)
rhie_table_s4 <- read_excel_allsheets("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/figures_for_manuscript/get_sex_assoc_regs/4_perform_zfx_analysis/raw_data/rhie_table_s4_cleaned.xlsx")
rhie_degs <- list()
rhie_degs[["C42B"]] <- rhie_table_s4[["c42b_sirna"]]$Gene
rhie_degs[["MCF7"]] <- rhie_table_s4[["mcf7_sirna"]]$Gene

rhie_dual_support_genes <- list()
rhie_dual_support_genes[["C42B"]] <- rhie_table_s4[["c42b_dual_support"]]$Gene
rhie_dual_support_genes[["MCF7"]] <- rhie_table_s4[["mcf7_dual_support"]]$Gene

#read in ontologies
##constraint scores
constraint_scores <- fread("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/supplementary_dataset_11_full_constraint_metrics.tsv", header = TRUE) #In Input_Files_Not_Generated_By_Scripts
constraint_scores <- as.data.frame(constraint_scores)
constraint_scores <- constraint_scores %>% as.data.frame() %>% dplyr::select(oe_lof_upper, gene) %>% filter(!is.na(oe_lof_upper))
constraint_scores <- constraint_scores %>% group_by(gene) %>% summarize(median_oe_lof_upper_across_transcripts = median(oe_lof_upper))
ploeuf_leq_0p35_genes <- constraint_scores %>% filter(median_oe_lof_upper_across_transcripts < 0.35) %>% pull(gene) %>% as.character()

##sfari genes 
sfari <- read.csv("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/sfari.csv", header = TRUE) #In Input_Files_Not_Generated_By_Scripts
SFARI_genes <- unique(sfari$gene.symbol)
high_conf_SFARI_genes <- unique((sfari %>% filter(gene.score == 1))$gene.symbol)

##sysID
id_genes <- read.csv("/n/groups/walsh/indData/becky/SysNDD_figures_250823/Input/SysNDD_ID_genes_definitive.csv", header = TRUE) #In Input_Files_Not_Generated_By_Scripts
just_id_genes <- setdiff(id_genes$symbol, SFARI_genes)

#consistently downregulated genes
consistent_downreg_DEG_in_asd_kd_genes <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/figures_for_manuscript/get_sex_assoc_regs/0_get_consistent_ASD_DEGs/consistent_downreg_DEG_in_asd_kd_genes.RDS")
consistent_downreg_DEG_in_asd_kd_genes_constrained <- intersect(consistent_downreg_DEG_in_asd_kd_genes, ploeuf_leq_0p35_genes)

##put ontologies in a single list
pathway_list <- list()
pathway_list[["SFARI"]] <- SFARI_genes
pathway_list[["HCSFARI"]] <- high_conf_SFARI_genes
pathway_list[["SysID"]] <- just_id_genes
pathway_list[["Candidate_Autism_Risk_Genes"]] <- setdiff(consistent_downreg_DEG_in_asd_kd_genes_constrained, SFARI_genes)
pathway_list[["Downregulated DEGs"]] <- consistent_downreg_DEG_in_asd_kd_genes
pathway_list[["Top Drivers"]] <- final_candidates
pathway_list[["ZFX Regulon"]] <- regulon_list[["ZFX"]]
pathway_list[["Female Biased downDEG"]] <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Revisions/Redo_Figure_6_And_S6_To_Add_Decasien_Data/Figure_6A_S6A_Heatmap_Expression_Candidate_Contributors_FPE/prioritized_regs_consistent_downreg_DEG_in_asd_kd_genes_and_potential_contributors_to_the_fpe.RDS")

##################
#CREATE PLOT BETWEEN ZFX AND 70 CANDIDATE REGULATORS
##################
#ZFX is central to the graph
outer_nodes <- c("ZFX")

#show edges between nodes if they are in ZFX's regulon or vice versa 
nodes_in_ZFX_regulon <- list("ZFX" = intersect(regulon_list[["ZFX"]], final_candidates))
nodes_in_which_ZFX_is_in_regulon <- names(regulon_list[final_candidates])[which(sapply(regulon_list[final_candidates], function(x) "ZFX" %in% x))]

edges <- rbind(setNames(stack(nodes_in_ZFX_regulon), c("target", "regulator")), data.frame(target = "ZFX", regulator = nodes_in_which_ZFX_is_in_regulon, stringsAsFactors = FALSE))
colnames(edges) <- c("Regulatee", "Regulator")
edges <- edges %>% dplyr::select(Regulator, Regulatee) %>% filter(!Regulator == Regulatee)
edges$Regulator <- as.character(edges$Regulator)
edges$Regulatee <- as.character(edges$Regulatee)
inner_nodes <- unique(setdiff(c(edges$Regulator, edges$Regulatee), outer_nodes))


#make graph
g <- graph_from_data_frame(edges, directed = TRUE)


#set layer attribute
V(g)$layer <- NA
V(g)$layer[V(g)$name %in% outer_nodes] <- 1
V(g)$layer[V(g)$name %in% inner_nodes] <- 2


#make layout
n1 <- sum(V(g)$layer == 1)
n2 <- sum(V(g)$layer == 2)


theta1 <- seq(0, 2*pi, length.out = n1+1)[-1]
incoming <- edges %>% filter(Regulatee == "ZFX") %>% pull(Regulator)
outgoing <- edges %>% filter(Regulator == "ZFX") %>% pull(Regulatee)
bidirectional <- intersect(incoming, outgoing)
incoming_only <- setdiff(incoming, outgoing)
outgoing_only <- setdiff(outgoing, incoming)
sorted_inner_nodes <- c(bidirectional, outgoing_only, incoming_only)
theta2_vals <- seq(0, 2*pi, length.out = length(sorted_inner_nodes) + 1)[-1]
names(theta2_vals) <- sorted_inner_nodes
names(theta1) <- outer_nodes


r1 <- 1
r2 <- 1.06


layout <- matrix(NA, nrow = vcount(g), ncol = 2)
layout[V(g)$layer == 1, ] <- matrix(c(0, 0), nrow = 1)  # place ZFX at center
layout[V(g)$layer == 2, ] <- cbind(r2 * cos(theta2_vals[V(g)$name[V(g)$layer == 2]]),
                                   r2 * sin(theta2_vals[V(g)$name[V(g)$layer == 2]]))


#set node labels
V(g)$label <- V(g)$name


#set edge styles
E(g)$color <- "black"


E(g)$width <- 0.5


E(g)$arrow.mode <- 1


E(g)$lty <- "solid"


#set node sizes
V(g)$size <- 20


#set node colors
V(g)$color <- NA
V(g)$color[is.na(V(g)$color) & !V(g)$name %in% SFARI_genes] <- "#0071C5"
V(g)$color[is.na(V(g)$color) & V(g)$name %in% high_conf_SFARI_genes] <- "#fddb27ff"
V(g)$color[is.na(V(g)$color) & V(g)$name %in% setdiff(SFARI_genes, high_conf_SFARI_genes)] <- "#ee4266"
V(g)$frame.color <-  NA
V(g)$frame.width <- 0
V(g)$label.cex <- 0.8  # default size for all nodes
V(g)$label.cex[V(g)$name %in% inner_nodes] <- 0.5 
bidirectional_edges <- apply(as_data_frame(g, what = "edges")[, 1:2], 1, function(x) {
  any(edges$Regulator == x[2] & edges$Regulatee == x[1])
})
E(g)$curved <- ifelse(bidirectional_edges, 0.2, 0)


#save plot
pdf("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Actually_Submitted/For_Third_Submission_And_What_We_Published/7_CLEANED_7_Figure_5_And_S10/Figure_5CD_ZFX_GRN_Vis_And_ChIP_Seq_Enr_Barplots/zfx_top_candidates_graph.pdf")


plot(
  g,
  layout = layout,
  vertex.size = V(g)$size,
  vertex.color = V(g)$color,
  edge.arrow.size = 0.3,
  edge.arrow.mode = 2,
  edge.lty = E(g)$lty,
  edge.curved = E(g)$curved
)

#close plot
dev.off()

##################
#RUN ENRICHMENT ANALYSIS (ZFX TARGETS RE-ANALYZED BY SAN ROMAN 2024)
##################
#get dataframe (list1 = the ZFX ChIP target set for that cell type, list2 = each gene set/ontology in pathway_list)
enr_df <- sapply(pathway_list, function(x) generate_overlap_p_val_and_or(rhie_targets[["MCF7"]], x, universe = rhie_targets_universe[["MCF7"]])) %>% as.data.frame() %>% t() %>% as.data.frame()
enr_df$cellType <- "MCF7"
enr_df <- rbind(enr_df, sapply(pathway_list, function(x) generate_overlap_p_val_and_or(rhie_targets[["C42B"]], x, universe = rhie_targets_universe[["C42B"]])) %>% as.data.frame() %>% t() %>% as.data.frame() %>% mutate(cellType = "C42B"))
enr_df <- rbind(enr_df, sapply(pathway_list, function(x) generate_overlap_p_val_and_or(common_zfx_targets, x, universe = intersect(rhie_targets_universe[["MCF7"]], rhie_targets_universe[["C42B"]]))) %>% as.data.frame() %>% t() %>% as.data.frame() %>% mutate(cellType = "Common_Targets"))
colnames(enr_df) <- c("pval", "or", "lower.ci", "InList1AndList2", "InList1AndNotList2", "InList2AndNotList1", "InNeither", "cellType")
enr_df$padj <- p.adjust(enr_df$pval, method = "BH")
enr_df <- enr_df %>% arrange(desc(or))
enr_df$ontology <- sapply(row.names(enr_df), function(x) gsub("1", "", gsub("2", "", x, fixed = TRUE), fixed = TRUE))
row.names(enr_df) <- NULL

#write out the matrix values as source data for the reviewer comment
write.table(enr_df, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Actually_Submitted/For_Third_Submission_And_What_We_Published/7_CLEANED_7_Figure_5_And_S10/Figure_5CD_ZFX_GRN_Vis_And_ChIP_Seq_Enr_Barplots/zfx_ChIP_enrichment_source_data.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#set up plot
enr_df$ontology <- factor(enr_df$ontology, levels = unique(as.character(enr_df$ontology)))
enr_df$sig <- sapply(enr_df$padj, function(x) ifelse(x <0.0001, "< 0.0001", ifelse(x < 0.001, "< 0.001", ifelse(x < 0.01, "< 0.01", ifelse(x < 0.05, "< 0.05", "Not Significant")))))
enr_df$sig <- factor(as.character(enr_df$sig), 
                                levels = c("< 0.0001", "< 0.001", "< 0.01", "< 0.05", "Not Significant"))
cols <- c(brewer.pal(8, "Greens")[c(8, 6, 4, 2)], "lightgray")
names(cols) <- c("< 0.0001", "< 0.001", "< 0.01", "< 0.05", "Not Significant")
enr_df$cellType <- forcats::fct_rev(factor(enr_df$cellType, levels = c("C42B", "Common_Targets", "MCF7")))
print(enr_df)


pdf("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Actually_Submitted/For_Third_Submission_And_What_We_Published/7_CLEANED_7_Figure_5_And_S10/Figure_5CD_ZFX_GRN_Vis_And_ChIP_Seq_Enr_Barplots/zfx_ChIP_enrichment.pdf", width = 6, height = 3)
ggplot(enr_df, aes(y = forcats::fct_rev(ontology), x = or, fill = sig, group = cellType)) + 
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