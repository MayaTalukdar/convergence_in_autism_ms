#!/usr/bin/env Rscript
reg_list <- "All_TFs_And_Epigenetic_Regulators"
library(viper)
library(data.table)
library(tidyverse)

print("*************************************")
print(reg_list)
print("*************************************")

###################################
#GENERATE REGULON OBJECT
###################################
#read in adjacency matrix
##modified based on https://support.bioconductor.org/p/123998/
adjfile <- as.data.frame(fread(paste0("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/full_aracne_network_analysis/Output/", reg_list, "/aracne_output/network.txt"), header = TRUE))
adjfile$pvalue <- NULL
colnames(adjfile) <- NULL
write.table(adjfile, paste0("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/full_aracne_network_analysis/Output/", reg_list, "/aracne_output/network_formatted.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#read in expression matrix
dset <- as.data.frame(fread("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/full_aracne_network_analysis/Input/matrix.txt", header = TRUE))
row.names(dset) <- dset$gene
dset$gene <- NULL
dset <- as.matrix(dset)

#create regulon object
regul <- aracne2regulon(paste0("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/full_aracne_network_analysis/Output/", reg_list, "/aracne_output/network_formatted.txt"), dset, verbose = TRUE, format = "3col")

###################################
#GENERATE GENE EXPRESSION SIGNATURES & FIND MASTER REGS
###################################
viper_res_dir <- paste0("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/full_aracne_network_analysis/Output/", reg_list, "/viper_results/")
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

##per perturbation
targets <-  unname(unlist(read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/final_targets_to_carry_forward.csv", header = FALSE)))
for (target in targets)
{
    current_perturb_path <- paste0(viper_res_dir, target, "/")
    dir.create(current_perturb_path)

    ###subset expression matrix based on perturbations vs controls
    control_cols <- colnames(dset)[which(grepl("NTC", colnames(dset)))]
    dset_perturbations <- dset[,-which(colnames(dset) %in% control_cols)]

    ###further filter to only take columns from current perturbation
    dset_perturbations <- dset_perturbations[,which(grepl(paste0(target, "$"), colnames(dset_perturbations)))]
    ncol(dset_perturbations) == 3

    ###extract controls
    dset_controls <- dset[,control_cols]

    ###calculate signature
    signature <- bootstrapTtest(dset_perturbations, dset_controls)

    ####calculate null distribution 
    nullmodel <- ttestNull(dset_perturbations, dset_controls)

    ###find master regulators
    mrs <- msviper(signature, regul, nullmodel, verbose = TRUE)
    saveRDS(mrs, paste0(current_perturb_path, "mrs.RDS"))

    print(paste0("Finished target ", target, "!"))
}

##just genes associated with ASD --> note: this is now the same as all perturbations (have maintained both code block sections so paths dont break)
current_perturb_path <- paste0(viper_res_dir, "just_genes_associated_with_asd/")
dir.create(current_perturb_path)
genes_to_use <- unname(unlist(read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/final_targets_to_carry_forward.csv", header = FALSE)))

###subset expression matrix based on perturbations vs controls
control_cols <- colnames(dset)[which(grepl("NTC", colnames(dset)))]
dset_perturbations <- dset[,-which(colnames(dset) %in% control_cols)]

###further filter to only take columns from current perturbation
dset_perturbations <- dset_perturbations[,which(grepl(paste0(paste(genes_to_use, collapse = "$|"), "$"), colnames(dset_perturbations)))]

###extract controls
dset_controls <- dset[,control_cols]

###calculate signature
signature <- bootstrapTtest(dset_perturbations, dset_controls)

####calculate null distribution 
nullmodel <- ttestNull(dset_perturbations, dset_controls)

###find master regulators
mrs <- msviper(signature, regul, nullmodel, verbose = TRUE)
saveRDS(mrs, paste0(current_perturb_path, "mrs.RDS"))
print("Finished just genes associated with ASD!")
