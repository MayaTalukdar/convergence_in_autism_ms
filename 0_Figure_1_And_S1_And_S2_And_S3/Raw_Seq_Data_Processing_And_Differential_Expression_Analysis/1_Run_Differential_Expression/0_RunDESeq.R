library(tidyverse)
library(DESeq2)

#read in ASD pool key 
sample_key <- read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/run_alignment/ASD_Screen_Sample_Key_AllReps_updated.txt", header = TRUE)
all_perturbations <- unname(unlist(read.csv("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/all_perturbations.csv", header = FALSE))) 

target_genes <- c(all_perturbations, "NTC1", "NTC2", "NTC3")
sample_key <- sample_key %>% filter(target_gene %in% target_genes)

#read in paths to Alevin output
pathsToAlevinOutput_withDedup <- setNames(unlist(read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/run_alignment/pathsToAlevinOutput_withDedup.txt", header = FALSE)), NULL)
#read in Alevin output
processed_output <- list()
counter <- 1 
for (path in pathsToAlevinOutput_withDedup)
{
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
	processed_mat <- mat[,-which(!colnames(mat) %in% metadata$index_seq)]
	colnames(processed_mat) <- sapply(colnames(processed_mat), function(x) metadata$target_gene[which(metadata$index_seq == x)])
	colnames(processed_mat) <- sapply(colnames(processed_mat), function(x) paste0(Pool, ".", x))

	processed_output[[counter]] <- list(processed_mat, metadata)

	counter <- counter + 1
}

#subset to only include intersected genes
intersected_genes <- Reduce(intersect, lapply(lapply(processed_output, function(x) x[[1]]), row.names))
processed_mat_intersected <- list()
for (index in seq_along(processed_output))
{
        mat <- processed_output[[index]][[1]]
        mat <- mat[row.names(mat) %in% intersected_genes,]
        processed_mat_intersected[[index]] <- mat
}

#merge processed output for each screen sample
processed_mat_merged <- do.call(cbind, processed_mat_intersected)
saveRDS(processed_mat_merged, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/processed_mat_merged.RDS")
metadata_merged <- do.call(rbind, lapply(processed_output, function(x) x[[2]]))

#format metadata correctly
row.names(metadata_merged) <- sapply(seq_len(nrow(metadata_merged)), function(x) paste0(metadata_merged$pool[x], ".", metadata_merged$target_gene[x]))
metadata_merged <- metadata_merged[order(row.names(metadata_merged)),]
processed_mat_merged <- processed_mat_merged[,order(colnames(processed_mat_merged))]
colnames(processed_mat_merged) == row.names(metadata_merged)
metadata_merged$Batch <- sapply(metadata_merged$pool, function(x) as.numeric(substr(x, 1, 1)))
metadata_merged$Batch <- sapply(metadata_merged$Batch, function(x) paste0("batch_", x))
metadata_merged$Batch <- as.factor(metadata_merged$Batch)
metadata_merged$Batch <- relevel(metadata_merged$Batch, ref = "batch_1")
metadata_merged$target_gene <- sapply(metadata_merged$target_gene, function(x) gsub("-", "_", x))
metadata_merged$sgRNA_not_NTC_merged <- metadata_merged$target_gene
metadata_merged$sgRNA_NTC_1_2_merged <- sapply(metadata_merged$sgRNA_not_NTC_merged, function(x) ifelse(x %in% c("NTC1", "NTC2"), "NTC_1_2", x))
metadata_merged$sgRNA_NTC_1_3_merged <- sapply(metadata_merged$sgRNA_not_NTC_merged, function(x) ifelse(x %in% c("NTC1", "NTC3"), "NTC_1_3", x))
metadata_merged$sgRNA_NTC_2_3_merged <- sapply(metadata_merged$sgRNA_not_NTC_merged, function(x) ifelse(x %in% c("NTC2", "NTC3"), "NTC_2_3", x))
metadata_merged$target_gene <- sapply(metadata_merged$target_gene, function(x) ifelse(x %in% c("NTC1", "NTC2", "NTC3"), "NTC", x))
metadata_merged$target_gene <- as.factor(metadata_merged$target_gene)
metadata_merged$target_gene <- relevel(metadata_merged$target_gene, ref = "NTC")
metadata_merged$sgRNA_NTC_1_2_merged <- relevel(as.factor(metadata_merged$sgRNA_NTC_1_2_merged), ref = "NTC_1_2")
metadata_merged$sgRNA_NTC_1_3_merged <- relevel(as.factor(metadata_merged$sgRNA_NTC_1_3_merged), ref = "NTC_1_3")
metadata_merged$sgRNA_NTC_2_3_merged <- relevel(as.factor(metadata_merged$sgRNA_NTC_2_3_merged), ref = "NTC_2_3")
colnames(metadata_merged)[5] <- "sgRNA"

#run DESeq2
dds <- DESeqDataSetFromMatrix(
 countData = round(processed_mat_merged),
  colData = metadata_merged,
  design= ~ Batch + sgRNA)
dds <- DESeq(dds)
saveRDS(dds, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/dds.RDS") 
dds <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/dds.RDS") 

#get results 
control <- "NTC"

test_sgRNAs <- as.character(unique(metadata_merged$sgRNA)[-which(unique(metadata_merged$sgRNA) == control)])
path <- "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/DEG_Results/OriginalResultsBeforeManualAdjustment"

for (sgRNA in test_sgRNAs)
{
	cleaned_sgRNA <- gsub("_", "-", sgRNA)
       write.table(results(dds, contrast = c("sgRNA", sgRNA, control)), paste0(path, "/", cleaned_sgRNA, "_vs_", control, "_withDedup.txt"))
        print(paste0(sgRNA, " completed!"))
}

#run DESeq2 to get results for NTC 3
path <- "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/DEG_Results/OriginalResultsBeforeManualAdjustment"
library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = round(processed_mat_merged),
  colData = metadata_merged,
  design= ~Batch + sgRNA_NTC_1_2_merged)
dds <- DESeq(dds)
saveRDS(dds, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/dds_NTC_1_2_merged.RDS")                         

write.table(results(dds, contrast = c("sgRNA_NTC_1_2_merged", "NTC_1_2", "NTC3")), paste0(path, "/NTC-3_vs_NTC_1_2_withDedup.txt"))

#run DESeq2 to get results for NTC 2
dds <- DESeqDataSetFromMatrix(
  countData = round(processed_mat_merged),
  colData = metadata_merged,
  design= ~Batch + sgRNA_NTC_1_3_merged)
dds <- DESeq(dds)
saveRDS(dds, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/dds_NTC_1_3_merged.RDS")

write.table(results(dds, contrast = c("sgRNA_NTC_1_3_merged", "NTC_1_3","NTC2")), paste0(path, "/NTC-2_vs_NTC_1_3_withDedup.txt"))

#run DESeq2 to get results for NTC 1
dds <- DESeqDataSetFromMatrix(
  countData = round(processed_mat_merged),
  colData = metadata_merged,
  design= ~Batch + sgRNA_NTC_2_3_merged)
dds <- DESeq(dds)
saveRDS(dds, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/dds_NTC_2_3_merged.RDS")                         

write.table(results(dds, contrast = c("sgRNA_NTC_2_3_merged", "NTC1", "NTC_2_3")), paste0(path, "/NTC-1_vs_NTC_2_3_withDedup.txt"))




