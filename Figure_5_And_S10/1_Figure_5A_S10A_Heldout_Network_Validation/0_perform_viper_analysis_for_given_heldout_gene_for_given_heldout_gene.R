#!/usr/bin/env Rscript
current_gene <- commandArgs(trailingOnly=TRUE)[1]
reg_list <- "All_TFs_And_Epigenetic_Regulators"
library(viper)
library(data.table)
library(tidyverse)

print("*************************************")
print(current_gene)
print("*************************************")

###################################
#GENERATE REGULON OBJECT
###################################
#read in adjacency matrix
##modified based on https://support.bioconductor.org/p/123998/
adjfile <- as.data.frame(fread(paste0("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/half_heldout_aracne_network_analysis/","/Output/", current_gene, "/", reg_list, "/aracne_output/network.txt"), header = TRUE))
adjfile$pvalue <- NULL
colnames(adjfile) <- NULL
write.table(adjfile, paste0("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/half_heldout_aracne_network_analysis/","/Output/", current_gene, "/", reg_list, "/aracne_output/network_formatted.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#read in expression matrix
dset <- as.data.frame(fread(paste0("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/half_heldout_aracne_network_analysis/Input/", current_gene, "/matrix.txt"), header = TRUE))
row.names(dset) <- dset$gene
dset$gene <- NULL
dset <- as.matrix(dset)

#create regulon object
regul <- aracne2regulon(paste0("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/half_heldout_aracne_network_analysis/Output/", current_gene, "/", reg_list, "/aracne_output/network_formatted.txt"), dset, verbose = TRUE, format = "3col")

###################################
#GENERATE GENE EXPRESSION SIGNATURES & FIND MASTER REGS
###################################
viper_res_dir <- paste0("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/half_heldout_aracne_network_analysis/","/Output/", current_gene, "/", reg_list, "/viper_results/")
dir.create(viper_res_dir)

##across all perturbations vs all controls
all_perturb_path <- paste0(viper_res_dir, "all_perturbations/")
dir.create(all_perturb_path)

###subset expression matrix based on perturbations vs controls
control_cols <- colnames(dset)[which(grepl("NTC", colnames(dset)))]
dset_perturbations <- dset[,-which(colnames(dset) %in% control_cols)]

###extract controls
dset_controls <- dset[,control_cols]

###calculate signature
signature <- bootstrapTtest(dset_perturbations, dset_controls)

####calculate null distribution 
nullmodel <- ttestNull(dset_perturbations, dset_controls)

###find master regulators
mrs <- msviper(signature, regul, nullmodel, verbose = TRUE)
saveRDS(mrs, paste0(all_perturb_path, "mrs.RDS"))

