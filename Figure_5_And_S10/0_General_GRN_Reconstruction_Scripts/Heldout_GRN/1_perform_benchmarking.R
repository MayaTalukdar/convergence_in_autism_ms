#!/usr/bin/env Rscript
reg_list <- "All_TFs_And_Epigenetic_Regulators"

##########################
#I/O
##########################
#set up libraries
library(viper)
library(data.table)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(UpSetR)
library(grid)
library(tidyverse)
library(readxl)
library(pheatmap)
library(RColorBrewer)
library(data.table)
library(viridis) 
library(clusterProfiler)
library(org.Hs.eg.db)
library(fgsea)
library(ggrepel)
library(patchwork)
library(cowplot)
library(igraph)

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

get_custom_degree_neighbors <- function(graph, node, degree)
{

    neighbors <-unique(do.call(c, lapply(seq_len(degree), function(x) names(unlist(ego(graph, order = x, nodes = node, mode = "out"))))))
}

get_second_degree_regs_that_are_pos_regs_of_first_degree_regs <- function(network, regulator)
{
    pos_regulatees <- network %>% filter(Regulator == regulator) %>% filter(cor == "POS") %>% pull(Regulatee)
    second_degree_regulon_pos <- unique(c(pos_regulatees, do.call(c, sapply(pos_regulatees, function(x) network %>% filter(Regulator == x) %>% filter(cor == "POS") %>% pull(Regulatee)))))

    neg_regulatees <- network %>% filter(Regulator == regulator) %>% filter(cor == "NEG") %>% pull(Regulatee)
    second_degree_regulon_neg <- unique(c(neg_regulatees, do.call(c, sapply(neg_regulatees, function(x) network %>% filter(Regulator == x) %>% filter(cor == "POS") %>% pull(Regulatee)))))

    return (list(pos_second_degree_reg = second_degree_regulon_pos, neg_second_degree_reg = second_degree_regulon_neg))
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

all_regs <- unname(unlist(read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/half_heldout_aracne_network_analysis/Output/regs_to_test.txt", header = FALSE)))
baseMean_in_NTCs <- read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/baseMeanInNTCs.csv", row.names = 1)
universe <- row.names(baseMean_in_NTCs %>% filter(baseMean >= 10))

#read in overall data matrix
processed_mat_merged <- fread("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/full_aracne_network_analysis/Input/matrix.txt", header = TRUE) %>% as.data.frame() %>% filter(gene %in% universe)
row.names(processed_mat_merged) <- processed_mat_merged$gene
processed_mat_merged$gene <- NULL
rownames <- row.names(processed_mat_merged)
processed_mat_merged <- apply(processed_mat_merged, 2, function(x) scale(x))
row.names(processed_mat_merged) <- rownames

#check what target genes were considered as regulators for aracne 
sgRNAs <- unname(unlist(read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/final_targets_to_carry_forward.csv")))
potential_regs <- unname(unlist(read.csv("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/half_heldout_aracne_network_analysis/Input/RegulatorsList/ProcessedLists/all_NPC_expressed_TFs_and_epigenetic_regulators.csv", header = FALSE)))
intersect(potential_regs, sgRNAs)

#get size of first and second degree regulon in the full network - becky approach
for (current_gene in all_regs)
{
    print("*********")
    print(current_gene)
    
    #read in network 
    full_network <- read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/full_aracne_network_analysis/Output/All_TFs_And_Epigenetic_Regulators/aracne_output/network_formatted.txt", header = FALSE)
    colnames(full_network) <- c('Regulator', 'Regulatee', 'weight') #note: weight is MI
    full_network <- full_network %>% filter(Regulator %in% universe) %>% filter(Regulatee %in% universe)
    first_degree_reg <- unique(full_network %>% filter(Regulator == current_gene) %>% pull(Regulatee))
    second_degree_reg <- unique(full_network %>% filter(Regulator %in% first_degree_reg) %>% pull(Regulatee))
    
    print(length(unique(c(first_degree_reg, second_degree_reg))))
    print("*********")
}

#####################
#DETERMINE IF DEGS ARE ENRICHED FOR REGULONS
#####################
enr_pval_list <- c()
plot_list <- list()

for (current_gene in all_regs)
{
    print(current_gene)
    
    #read in network results
    viper_res_dir <- paste0("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/half_heldout_aracne_network_analysis/","/Output/", current_gene, "/", reg_list, "/viper_results/")
    network <- read.table(paste0("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/half_heldout_aracne_network_analysis/Output/", current_gene, "/All_TFs_And_Epigenetic_Regulators/aracne_output/network_formatted.txt")) 
    colnames(network) <- c('Regulator', 'Regulatee', 'weight') #note: weight is MI
    network <- network %>% filter(Regulator %in% universe) %>% filter(Regulatee %in% universe)
    network$cor <- apply(network, 1, function(x) get_direction_of_regulation(x["Regulator"], x["Regulatee"], processed_mat_merged))
    graph <- graph.data.frame(network, directed = TRUE)
    held_out_genes <- current_gene
    identified_regulators <- unique(network$Regulator)
    identified_regulators_that_we_heldout <- intersect(identified_regulators, held_out_genes) 

    #define the universe as all npc genes expressed with baseMean >= 10 
    basemean_df <-  read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/baseMeanInNTCs.csv", row.names = 1)
    universe <- row.names(basemean_df %>% filter(baseMean > 10))

    #read in viper results
    mrs <- readRDS(paste0("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/half_heldout_aracne_network_analysis/","/Output/", current_gene, "/", reg_list, "/viper_results/all_perturbations/mrs.RDS"))

    #read in differential expression results
    ##read in complete DESeq results
    path <- "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/DEG_Results/WithManualAdjustment/"
    all_DESeq_res_list <- lapply(identified_regulators_that_we_heldout, function(x) read.table(paste0(path, x, "_vs_NTC_withDedup_manual_adjustment_Cutoff_Of_10.txt"), header = TRUE) %>% filter(!is.na(padj)) %>% arrange(desc(log2FoldChange)))
    names(all_DESeq_res_list) <- identified_regulators_that_we_heldout

    ##get only significant DEGs
    DEG_res_list <- lapply(all_DESeq_res_list, function(x) row.names(x %>% filter(padj < 0.05)))
    names(DEG_res_list) <- identified_regulators_that_we_heldout
    neg_DEG_res_list <- setNames(lapply(all_DESeq_res_list, function(x) row.names(x %>% filter(padj < 0.05) %>% filter(log2FoldChange < 0))), identified_regulators_that_we_heldout)

    #set up regulon 
    regulon_list_indirect_pos <- setNames(lapply(identified_regulators_that_we_heldout, function(x) get_second_degree_regs_that_are_pos_regs_of_first_degree_regs(network, x)[[1]]), identified_regulators_that_we_heldout)
    regulon_list_direct <- setNames(lapply(identified_regulators_that_we_heldout, function(x) (network %>% filter(Regulator == x))$Regulatee), identified_regulators_that_we_heldout)

    #are genes in these regulons actually ppDEGs?
    enr_pval_list[[current_gene]] <- sapply(identified_regulators_that_we_heldout, function(x) generate_overlap_p_val(regulon_list_indirect_pos[[x]], neg_DEG_res_list[[x]], universe))

    #create plot
    regulon_indirect_pos <- regulon_list_indirect_pos[[current_gene]]
    regulon_direct_pos <- regulon_list_direct[[current_gene]]
    DEG_res <- all_DESeq_res_list[[current_gene]]
    DEG_res$gene <- rownames(DEG_res)
    DEG_res <- DEG_res[!is.na(DEG_res$padj), ]

    DEG_res$regulon_type <- "Other"
    DEG_res$regulon_type[DEG_res$gene %in% regulon_indirect_pos] <- "Indirect"
    DEG_res$regulon_type[DEG_res$gene %in% regulon_direct_pos] <- "Direct"

    DEG_res$fill_color <- "lightgray"  # default
    DEG_res$fill_color[DEG_res$regulon_type != "Other" & DEG_res$padj >= 0.05] <- "gray"
    DEG_res$fill_color[DEG_res$regulon_type != "Other" & DEG_res$padj < 0.05 & DEG_res$log2FoldChange > 0] <- "red"
    DEG_res$fill_color[DEG_res$regulon_type != "Other" & DEG_res$padj < 0.05 & DEG_res$log2FoldChange < 0] <- "blue"

    indirect_layer <- DEG_res[DEG_res$regulon_type == "Indirect", ]
    direct_layer   <- DEG_res[DEG_res$regulon_type == "Direct", ]

    highlight_layer <- rbind(indirect_layer, direct_layer)
    x_max <- max(abs(highlight_layer$log2FoldChange), na.rm = TRUE)
    xlim_range <- c(-x_max, x_max)

    p <- ggplot() +
        # Indirect: plot first, no border
        geom_point(data = indirect_layer,
                aes(x = log2FoldChange, y = -log10(padj), fill = fill_color),
                shape = 21, stroke = NA, size = 1.5, alpha = 0.8) +
        
        # Direct: plot second, black border
        geom_point(data = direct_layer,
                aes(x = log2FoldChange, y = -log10(padj), fill = fill_color),
                shape = 21, color = "black", stroke = 0.6, size = 1.5, alpha = 0.8) +
        
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
        scale_fill_identity() +
        labs(title = paste("Volcano Plot:", current_gene),
            x = "Log2 Fold Change",
            y = "-log10 Adjusted P-value") +
        xlim(xlim_range) +
        theme_minimal()

    plot_list[[current_gene]] <- p
}

p.adjust(unlist(enr_pval_list), method = "BH")

pdf("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/half_heldout_aracne_network_analysis/Output/pos_regulon_volcano_plot.pdf", height = 9, width = 12)
wrap_plots(plot_list, nrow = 3)
dev.off()







