library(tidyverse)
library(DESeq2)

#read in metadata matching sample name to perturbation
metadata <- read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_denault_perturbation_data/geo_files/sample_metadata.txt", header = FALSE) #In Input_Files_Not_Generated_By_Scripts
colnames(metadata) <- c("SampleName", "Description")
metadata$celltype <- sapply(metadata$Description, function(x) strsplit(x, "-")[[1]][2])
metadata$perturbation <- sapply(metadata$Description, function(x) strsplit(x, "-")[[1]][3])
row.names(metadata) <- metadata$SampleName

#######################################
#iPSCs
#######################################
wd <- '/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_denault_perturbation_data/geo_files/ipsc_files/'

#read in all perturbations 
perturbation_data <- lapply(list.files(wd), function(x) read.table(paste0(wd, x), header = TRUE))
names(perturbation_data) <- sapply(list.files(wd), function(x) strsplit(x, "_")[[1]][1])

#create single matrix
table(sapply(perturbation_data, nrow))
gene_names <- perturbation_data[[1]]$Gene
perturbation_data <- as.data.frame(do.call(cbind, lapply(perturbation_data, function(x) x$Counts)))
row.names(perturbation_data) <- gene_names
metadata_subset <- metadata %>% filter(celltype == "iPSC")

#run DESeq2
dds <- DESeqDataSetFromMatrix(
 countData = perturbation_data,
  colData = metadata_subset,
  design= ~perturbation)
dds <- DESeq(dds)
saveRDS(dds, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_denault_perturbation_data/differential_expression_analysis/deg_results/ipscs/dds.RDS")

#get results 
control <- "control"
test_sgRNAs <- as.character(unique(metadata_subset$perturbation)[-which(unique(metadata_subset$perturbation) == control)])
path <- "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_denault_perturbation_data/differential_expression_analysis/deg_results/ipscs/Original/"

for (sgRNA in test_sgRNAs)
{
	cleaned_sgRNA <- gsub("_", "-", sgRNA)
  write.table(results(dds, contrast = c("perturbation", sgRNA, control)), paste0(path, "/", cleaned_sgRNA, "_vs_", control, "_DEG_res.txt"))
  print(paste0(sgRNA, " completed!"))
}

#######################################
#Neurons
#######################################
wd <- '/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_denault_perturbation_data/geo_files/neuron_files/'

#read in all perturbations 
perturbation_data <- lapply(list.files(wd), function(x) read.table(paste0(wd, x), header = TRUE))
names(perturbation_data) <- sapply(list.files(wd), function(x) strsplit(x, "_")[[1]][1])

#create single matrix
table(sapply(perturbation_data, nrow))
gene_names <- perturbation_data[[1]]$Gene
perturbation_data <- as.data.frame(do.call(cbind, lapply(perturbation_data, function(x) x$Counts)))
row.names(perturbation_data) <- gene_names
metadata_subset <- metadata %>% filter(celltype == "Neuron")

#run DESeq2
dds <- DESeqDataSetFromMatrix(
 countData = perturbation_data,
  colData = metadata_subset,
  design= ~perturbation)
dds <- DESeq(dds)
saveRDS(dds, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_denault_perturbation_data/differential_expression_analysis/deg_results/neurons/dds.RDS")

#get results 
control <- "control"
test_sgRNAs <- as.character(unique(metadata_subset$perturbation)[-which(unique(metadata_subset$perturbation) == control)])
path <- "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_denault_perturbation_data/differential_expression_analysis/deg_results/neurons/Original/"

for (sgRNA in test_sgRNAs)
{
	cleaned_sgRNA <- gsub("_", "-", sgRNA)
  write.table(results(dds, contrast = c("perturbation", sgRNA, control)), paste0(path, "/", cleaned_sgRNA, "_vs_", control, "_DEG_res.txt"))
  print(paste0(sgRNA, " completed!"))
}


