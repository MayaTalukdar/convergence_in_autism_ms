#!/usr/bin/env Rscript
.libPaths(c("/n/groups/walsh/indData/Maya/RLibs", .libPaths()))
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
library(igraph)
library(patchwork)

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

#perform fisher test with logs to prevent underflow
fisher_combined_p <- function(pvals) 
{
  # Ensure p-values are numeric
  pvals <- as.numeric(pvals)
  
  # Number of p-values
  k <- length(pvals)
  
  # -------------------------------
  # 1. Compute Fisher chi-squared statistic
  # Formula: X2 = -2 * sum(log(p_i))
  # -------------------------------
  X2 <- -2 * sum(log(pvals))
  
  # Degrees of freedom = 2 * number of p-values
  df <- 2 * k
  
  # -------------------------------
  # 2. Compute log of combined p-value
  # Use log.p = TRUE to avoid underflow for extremely small p-values
  # -------------------------------
  log_p_combined <- pchisq(X2, df, lower.tail = FALSE, log.p = TRUE)
  
  # -------------------------------
  # 3. Convert to 1e-notation
  # log10_p = log base 10 of p-value
  # exponent = floor(log10_p)
  # p_1e = leading factor (between 1 and 10)
  # -------------------------------
  log10_p <- log_p_combined / log(10)
  exponent <- floor(log10_p)
  p_1e <- 10^(log10_p - exponent)
  
  # -------------------------------
  # 4. Print nicely
  # Example: Combined p-value ≈ 2.34e-5700
  # -------------------------------
  cat("Combined p-value ≈", signif(p_1e, 3), "e", exponent, "\n")
  
  # -------------------------------
  # 5. Return detailed results invisibly
  # -------------------------------
  invisible(list(
    X2 = X2,
    df = df,
    log_p = log_p_combined,
    p_1e = p_1e,
    exponent = exponent
  ))
}

all_regs <- unname(unlist(read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/half_heldout_aracne_network_analysis/Output/regs_to_test.txt", header = FALSE)))
all_regs <- setdiff(all_regs, c("SETD1B", "EHMT1")) #these regulators removed bc their regulon size is < 250
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

    #set up regulon 
    regulon_list_indirect_pos <- setNames(lapply(identified_regulators_that_we_heldout, function(x) get_second_degree_regs_that_are_pos_regs_of_first_degree_regs(network, x)[[1]]), identified_regulators_that_we_heldout)
    regulon_list_direct <- setNames(lapply(identified_regulators_that_we_heldout, function(x) (network %>% filter(Regulator == x))$Regulatee), identified_regulators_that_we_heldout)

   #create plot
    regulon_indirect_pos <- regulon_list_indirect_pos[[current_gene]]
    regulon_direct_pos   <- regulon_list_direct[[current_gene]]

    DEG_res <- all_DESeq_res_list[[current_gene]]
    DEG_res$gene <- rownames(DEG_res)
    DEG_res <- DEG_res[!is.na(DEG_res$padj), ]

    DEG_res$regulon_type <- "Other"
    DEG_res$regulon_type[DEG_res$gene %in% regulon_indirect_pos] <- "Indirect"
    DEG_res$regulon_type[DEG_res$gene %in% regulon_direct_pos]   <- "Direct"

    DEG_res$fill_color <- "lightgray"  # default
    DEG_res$fill_color[DEG_res$regulon_type != "Other" &
                        DEG_res$padj >= 0.05] <- "gray"
    DEG_res$fill_color[DEG_res$regulon_type != "Other" &
                        DEG_res$padj < 0.05 &
                        DEG_res$log2FoldChange > 0] <- "red"
    DEG_res$fill_color[DEG_res$regulon_type != "Other" &
                        DEG_res$padj < 0.05 &
                        DEG_res$log2FoldChange < 0] <- "blue"

    # Layers
    background_layer <- DEG_res[DEG_res$regulon_type == "Other", ]
    indirect_layer   <- DEG_res[DEG_res$regulon_type == "Indirect", ]
    direct_layer     <- DEG_res[DEG_res$regulon_type == "Direct", ]

    x_max <- max(abs(DEG_res$log2FoldChange), na.rm = TRUE)
    xlim_range <- c(-x_max, x_max)

    p <- ggplot() +

        # Background: all non-regulon genes
        geom_point(data = background_layer,
                aes(x = log2FoldChange,
                    y = -log10(padj)),
                color = "lightgray",
                size = 1.0,
                alpha = 0.5) +

        # Indirect regulon: no border
        geom_point(data = indirect_layer,
                aes(x = log2FoldChange,
                    y = -log10(padj),
                    fill = fill_color),
                shape = 21,
                stroke = NA,
                size = 1.8,
                alpha = 0.9) +

        # Direct regulon: black border
        geom_point(data = direct_layer,
                aes(x = log2FoldChange,
                    y = -log10(padj),
                    fill = fill_color),
                shape = 21,
                color = "black",
                stroke = 0.6,
                size = 1.8,
                alpha = 0.9) +

        geom_hline(yintercept = -log10(0.05),
                linetype = "dashed",
                color = "black") +
        geom_vline(xintercept = 0,
                linetype = "dashed",
                color = "black") +
        scale_fill_identity() +
        labs(title = paste("Volcano Plot:", current_gene),
            x = "Log2 Fold Change",
            y = "-log10 Adjusted P-value") +
        xlim(xlim_range) +
        theme_minimal()

    plot_list[[current_gene]] <- p
}

pdf("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Actually_Submitted/For_Rebuttal/Output/Figures/Volcano_Plots/grn_validation_volcano_plots.pdf", width = 20, height = 4)
wrap_plots(plot_list, nrow = 1)
dev.off()