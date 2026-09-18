library(tidyverse)

#read in metadata matching sample name to perturbation
metadata <- read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_denault_perturbation_data/geo_files/sample_metadata.txt", header = FALSE) #In Input_Files_Not_Generated_By_Scripts
colnames(metadata) <- c("SampleName", "Description")
metadata$celltype <- sapply(metadata$Description, function(x) strsplit(x, "-")[[1]][2])
metadata$perturbation <- sapply(metadata$Description, function(x) strsplit(x, "-")[[1]][3])
row.names(metadata) <- metadata$SampleName

#######################################
#iPSCs
#######################################
#get genes that are expressed >10 in control samples
wd <- '/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_denault_perturbation_data/geo_files/ipsc_files/'
perturbation_data <- lapply(list.files(wd), function(x) read.table(paste0(wd, x), header = TRUE))
names(perturbation_data) <- sapply(list.files(wd), function(x) strsplit(x, "_")[[1]][1])
table(sapply(perturbation_data, nrow))
gene_names <- perturbation_data[[1]]$Gene
perturbation_data <- as.data.frame(do.call(cbind, lapply(perturbation_data, function(x) x$Counts)))
row.names(perturbation_data) <- gene_names
metadata_subset <- metadata %>% filter(celltype == "iPSC")
control_samples <- (metadata_subset %>% filter(perturbation == "control"))$Sample
perturbation_data <- perturbation_data %>% select(control_samples)
baseMean <- apply(perturbation_data, 1, mean)
write.table(baseMean,"/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_denault_perturbation_data/differential_expression_analysis/deg_results/ipscs/baseMean.csv", sep = ",", quote = FALSE, col.names = FALSE)
genes_to_keep <- names(baseMean)[which(baseMean > 10)]

#manually adjust p-values 
path <- "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_denault_perturbation_data/differential_expression_analysis/deg_results/ipscs/Original/"
DESeq_res_list <- lapply(list.files(path), function(x) read.table(paste0(path, x), header = TRUE, row.names = 1))
DESeq_res_list_filtered <- lapply(DESeq_res_list, function(x) x %>% filter(!row.names(x) %in% genes_to_keep))
DESeq_res_list_filtered <- lapply(DESeq_res_list_filtered, function(x) x %>% mutate(padj_DEseq = padj) %>% mutate(padj = NA) %>% select(baseMean, log2FoldChange, lfcSE, stat, pvalue, padj_DEseq, padj))
DESeq_res_list <- lapply(DESeq_res_list, function(x) x %>% filter(row.names(x) %in% genes_to_keep) %>% arrange(padj))
DESeq_res_list <- lapply(DESeq_res_list, function(x) x %>% mutate(padj_DEseq = padj) %>% mutate(padj = p.adjust(pvalue, method = "BH")) %>% select(baseMean, log2FoldChange, lfcSE, stat, pvalue, padj_DEseq, padj))
DESeq_res_all_list <- lapply(seq_along(DESeq_res_list), function(x) rbind(DESeq_res_list[[x]], DESeq_res_list_filtered[[x]]))
DESeq_res_all_list <- lapply(DESeq_res_all_list, function(x) x %>% arrange(padj))
names(DESeq_res_all_list) <- sapply(list.files(path), function(x) gsub(".txt", "", x))

lapply(seq_along(DESeq_res_all_list), function(x) write.table(DESeq_res_all_list[[x]], paste0("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_denault_perturbation_data/differential_expression_analysis/deg_results/ipscs/WithManualAdjustment/", names(DESeq_res_all_list)[x], "_manual_adjustment_Cutoff_Of_10.txt")))

#######################################
#Neurons
#######################################
#get genes that are expressed >10 in control samples
wd <- '/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_denault_perturbation_data/geo_files/neuron_files/'
perturbation_data <- lapply(list.files(wd), function(x) read.table(paste0(wd, x), header = TRUE))
names(perturbation_data) <- sapply(list.files(wd), function(x) strsplit(x, "_")[[1]][1])
table(sapply(perturbation_data, nrow))
gene_names <- perturbation_data[[1]]$Gene
perturbation_data <- as.data.frame(do.call(cbind, lapply(perturbation_data, function(x) x$Counts)))
row.names(perturbation_data) <- gene_names
metadata_subset <- metadata %>% filter(celltype == "Neuron")
control_samples <- (metadata_subset %>% filter(perturbation == "control"))$Sample
perturbation_data <- perturbation_data %>% select(control_samples)
baseMean <- apply(perturbation_data, 1, mean)
write.table(baseMean,"/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_denault_perturbation_data/differential_expression_analysis/deg_results/neurons/baseMean.csv", sep = ",", quote = FALSE, col.names = FALSE)
genes_to_keep <- names(baseMean)[which(baseMean > 10)]

#manually adjust p-values 
path <- "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_denault_perturbation_data/differential_expression_analysis/deg_results/neurons/Original/"
DESeq_res_list <- lapply(list.files(path), function(x) read.table(paste0(path, x), header = TRUE, row.names = 1))
DESeq_res_list_filtered <- lapply(DESeq_res_list, function(x) x %>% filter(!row.names(x) %in% genes_to_keep))
DESeq_res_list_filtered <- lapply(DESeq_res_list_filtered, function(x) x %>% mutate(padj_DEseq = padj) %>% mutate(padj = NA) %>% select(baseMean, log2FoldChange, lfcSE, stat, pvalue, padj_DEseq, padj))
DESeq_res_list <- lapply(DESeq_res_list, function(x) x %>% filter(row.names(x) %in% genes_to_keep) %>% arrange(padj))
DESeq_res_list <- lapply(DESeq_res_list, function(x) x %>% mutate(padj_DEseq = padj) %>% mutate(padj = p.adjust(pvalue, method = "BH")) %>% select(baseMean, log2FoldChange, lfcSE, stat, pvalue, padj_DEseq, padj))
DESeq_res_all_list <- lapply(seq_along(DESeq_res_list), function(x) rbind(DESeq_res_list[[x]], DESeq_res_list_filtered[[x]]))
DESeq_res_all_list <- lapply(DESeq_res_all_list, function(x) x %>% arrange(padj))
names(DESeq_res_all_list) <- sapply(list.files(path), function(x) gsub(".txt", "", x))

lapply(seq_along(DESeq_res_all_list), function(x) write.table(DESeq_res_all_list[[x]], paste0("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_denault_perturbation_data/differential_expression_analysis/deg_results/neurons/WithManualAdjustment/", names(DESeq_res_all_list)[x], "_manual_adjustment_Cutoff_Of_10.txt")))
