#!/usr/bin/env Rscript
heldout_gene <- commandArgs(trailingOnly=TRUE)[1]
print(heldout_gene)
library(viper)
library(data.table)
library(tidyverse)

all_targets <-  unname(unlist(read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/final_targets_to_carry_forward.csv", header = FALSE)))
all_targets <- setdiff(all_targets, heldout_gene)
dir.create(paste0("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/half_heldout_aracne_network_analysis/Input/", heldout_gene))
write.table(all_targets, paste0("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/half_heldout_aracne_network_analysis/Input/", heldout_gene, "/genes_to_keep.txt"), sep = "\t", row.names = F, col.names = F, quote = F)
genes_of_interest <- c(all_targets, "NTC1", "NTC2", "NTC3") 
basemean_df <-  read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/baseMeanInNTCs.csv", row.names = 1)
genes_to_keep <- row.names(basemean_df %>% filter(baseMean > 10))

#read in ASD pool key 
sample_key <- read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/run_alignment/ASD_Screen_Sample_Key_AllReps_updated.txt", header = TRUE)
sample_key <- sample_key %>% filter(target_gene %in% genes_of_interest)
 
#read in paths to Alevin output
pathsToAlevinOutput_withDedup <- setNames(unlist(read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/run_alignment/pathsToAlevinOutput_withDedup.txt", header = FALSE)), NULL)
#read in Alevin output
processed_output <- list()
counter <- 1 
for (path in pathsToAlevinOutput_withDedup)
{
	print(counter)
	#get pool ID
	Pool <- strsplit(path, "/", fixed = TRUE)[[1]][11]
	
	#get metadata specific to pool 
	metadata <- sample_key %>% filter(pool == Pool)

	#remove any samples that were not infected
	if (length(which(grepl("uninf", metadata$target_gene))) !=0)
	{
		metadata <- metadata[-which(grepl("uninf", metadata$target_gene)),]
	}

	#process mat to add perturbation names 
	mat <- read.table(paste0(path, "/quants_mat_merged_by_gene.txt"))
	processed_mat <- mat[,-which(!colnames(mat) %in% metadata$index_seq),drop = FALSE]
	colnames(processed_mat) <- sapply(colnames(processed_mat), function(x) metadata$target_gene[which(metadata$index_seq == x)])
	colnames(processed_mat) <- sapply(colnames(processed_mat), function(x) paste0(Pool, ".", x))

	processed_output[[counter]] <- list(processed_mat, metadata)

	counter <- counter + 1
}

#subset to only include intersected genes
library(tidyverse)
intersected_genes <- Reduce(intersect, lapply(lapply(processed_output, function(x) x[[1]]), row.names))
processed_mat_intersected <- list()
for (index in seq_along(processed_output))
{
        mat <- processed_output[[index]][[1]]
        mat <- mat[row.names(mat) %in% intersected_genes,,drop = FALSE]
        processed_mat_intersected[[index]] <- mat
}

#merge processed output for each screen sample
processed_mat_merged <- do.call(cbind, processed_mat_intersected)
metadata_merged <- do.call(rbind, lapply(processed_output, function(x) x[[2]]))

#subset to genes to keep 
processed_mat_merged <- processed_mat_merged[genes_to_keep,]
processed_mat_merged$gene <- row.names(processed_mat_merged)
processed_mat_merged <- processed_mat_merged %>% relocate(gene)

#ensure we have retained all sgRNAs that we want to 
unique_targets <- unique(sapply(colnames(processed_mat_merged), function(x) strsplit(x, ".", fixed = TRUE)[[1]][2]))
setdiff(all_targets, unique_targets)

write.table(processed_mat_merged, paste0("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/half_heldout_aracne_network_analysis/Input/", heldout_gene, "/matrix.txt"), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)






