library(tidyverse)
library(ggraph)
library(tidygraph)
library(igraph)

######################
#I/O
######################
#read in sfari final_targets 
sfari <- read.csv("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/sfari.csv", header = TRUE)
all_sfari_final_targets <- sfari$gene.symbol
high_confidence_sfari_final_targets <- (sfari %>% filter(gene.score == 1))$gene.symbol

#read in deg results
final_targets <- unname(unlist(read.csv("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/final_targets_to_carry_forward.csv", header = FALSE)))
path <- "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/DEG_Results/WithManualAdjustment/"
all_DESeq_res_list <- lapply(final_targets, function(x) read.table(paste0(path, x, "_vs_NTC_withDedup_manual_adjustment_Cutoff_Of_10.txt"), header = TRUE))
names(all_DESeq_res_list) <- final_targets
all_DEGs_list <- lapply(all_DESeq_res_list, function(x) row.names(x %>% filter(padj < 0.05)))

########################################################
#CREATE ADJACENCY DATAFRAMES & NODES FOR HOLOVIEWS
########################################################
#all sfari final_targets - green same direction, red opposite direction
##pull out only sfari gene DEGs 
deg_list_to_use <- lapply(all_DEGs_list, function(x) intersect(x, all_sfari_final_targets))
##remove target itself from this list
deg_list_to_use <- setNames(lapply(names(deg_list_to_use), function(x) setdiff(deg_list_to_use[[x]], c(x))), final_targets)
##create adjacency matrix
are_concordant_dir <- function(fc1, fc2)
{
    if(fc1 > 0 && fc2 > 0)
    {
        return (TRUE)
    } else if (fc1 < 0 && fc2 < 0)
    {
        return (TRUE)
    }

    return (FALSE)
}
create_adj_matrix <- function(target_1, target_2, deg_list_to_use)
{
    degs_to_consider <- intersect(deg_list_to_use[[target_1]], deg_list_to_use[[target_2]])
    if (length(degs_to_consider) > 0)
    {
        concordant_dir_vec <- sapply(degs_to_consider, function(x) are_concordant_dir(all_DESeq_res_list[[target_1]][x,"log2FoldChange"], all_DESeq_res_list[[target_2]][x,"log2FoldChange"]))
        num_same_dir <- length(which(concordant_dir_vec))
        num_opp_dir <- length(which(!concordant_dir_vec)) 
        dir_df_list <- list() 

        if (num_same_dir > 0)
        {
            same_dir_df <- as.data.frame(t(as.data.frame(c(target_1, target_2, num_same_dir, "same"))))
            row.names(same_dir_df) <- NULL
            colnames(same_dir_df) <- c("source", "target", "value", "type")
            dir_df_list[["same"]] <- same_dir_df
        } 
        if (num_opp_dir > 0)
        {
            opp_dir_df <- as.data.frame(t(as.data.frame(c(target_1, target_2, num_opp_dir, "opp"))))
            row.names(opp_dir_df) <- NULL
            colnames(opp_dir_df) <- c("source", "target", "value", "type")
            dir_df_list[["opp"]] <- opp_dir_df
        } 
        dir_df <- do.call(rbind, dir_df_list)
        return(dir_df)                                                   
    }
}
adj_matrix_df_sfari <- do.call(rbind, lapply(seq_along(final_targets), function(x) do.call(rbind, lapply(seq(x + 1, length(final_targets)), function(y) create_adj_matrix(final_targets[x], final_targets[y], deg_list_to_use)))))
adj_matrix_df_sfari <- adj_matrix_df_sfari %>% filter(!source == target)
write.table(adj_matrix_df_sfari, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/circos_plot_cross_enrichment/Input_Files/adj_matrix_sfari.csv", col.names = TRUE, row.names = FALSE, sep = ",", quote = F)

#high confidence sfari final_targets - green same direction, red opposite direction 
##pull out only sfari gene DEGs 
deg_list_to_use <- lapply(all_DEGs_list, function(x) intersect(x, high_confidence_sfari_final_targets))
##remove target itself from this list
deg_list_to_use <- setNames(lapply(names(deg_list_to_use), function(x) setdiff(deg_list_to_use[[x]], c(x))), final_targets)
##create adjacency matrix
adj_matrix_df_high_conf_sfari <- do.call(rbind, lapply(seq_along(final_targets), function(x) do.call(rbind, lapply(seq(x + 1, length(final_targets)), function(y) create_adj_matrix(final_targets[x], final_targets[y], deg_list_to_use)))))
adj_matrix_df_high_conf_sfari <- adj_matrix_df_high_conf_sfari %>% filter(!source == target)
write.table(adj_matrix_df_high_conf_sfari, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/circos_plot_cross_enrichment/Input_Files/adj_matrix_high_conf_sfari.csv", col.names = TRUE, row.names = FALSE, sep = ",", quote = F)

#all final_targets - green same direction, red opposite direction 
deg_list_to_use <- all_DEGs_list
##remove target itself from this list
deg_list_to_use <- setNames(lapply(names(deg_list_to_use), function(x) setdiff(deg_list_to_use[[x]], c(x))), final_targets)
##create adjacency matrix
adj_matrix_df_all <- do.call(rbind, lapply(seq_along(final_targets), function(x) do.call(rbind, lapply(seq(x + 1, length(final_targets)), function(y) create_adj_matrix(final_targets[x], final_targets[y], deg_list_to_use)))))
adj_matrix_df_all <- adj_matrix_df_all %>% filter(!source == target)
write.table(adj_matrix_df_all, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/circos_plot_cross_enrichment/Input_Files/adj_matrix_all.csv", col.names = TRUE, row.names = FALSE, sep = ",", quote = F)
