comparison <- commandArgs(trailingOnly = TRUE)[1] #can be ASD, CNV, or case 
path <- paste0("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/analysis_of_mo_dias_15q_data/Output/DEG_Results/MAST/Case_vs_Control_Male/", comparison, "/")
dir.create(path)
print("**********")
print(comparison)
print("**********")

###################
#I/O
###################
library(Seurat)
library(tidyverse)

seurat_obj <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/mo_dias_15q_data/seurat_obj_alisa_caroline_15q_data_from_chunhui_UMB4643_sexID_fixed.RDS")

###################
#SET UP OBJECT FOR DIFFERENTIAL EXPRESSION
###################
#subset by sex 
Idents(seurat_obj) <- "sex"
seurat_obj <- subset(seurat_obj, idents = "Male")
Idents(seurat_obj) <- "SampleID"

#subset by status
seurat_obj$case_control <- sapply(seurat_obj$group, function(x) ifelse(x == "CON", "control", "case"))

if (comparison == "case") #lumping together 15q and idiopathic asd as cases 
{
    Idents(seurat_obj) <- "case_control"
} else if (comparison == "CNV") #only comparing 15q to controls
{
    Idents(seurat_obj) <- "group"
    seurat_obj <- subset(seurat_obj, idents = "ASD", invert = TRUE)
} else if (comparison == "ASD") #only comparing idiopathic autism to controls
{
    Idents(seurat_obj) <- "group"
    seurat_obj <- subset(seurat_obj, idents = "CNV", invert = TRUE)
}

###################
#RUN DIFFERENTIAL EXPRESSION
###################
Idents(seurat_obj) <- "our_collapsed_broad_celltypes"
types <- unique(Idents(seurat_obj))

for (type in types)
{
    print("*********************************")
    print(paste0("Starting ", type, "!"))
        
    tryCatch({
            DEG_res <- FindMarkers(seurat_obj, ident.1 = "case", ident.2 = "control", group.by = "case_control",  test.use = "MAST", latent.vars = c("SampleID"), subset.ident = type, min.cells.group = 100, min.pct = 0.10, logfc.threshold = 0)
            saveRDS(DEG_res, paste0(path, type, "_DEG_res.RDS"))
    }, error = function(e) {
            cat("Error occurred while processing", type, ":", conditionMessage(e), "\n")
    })
}
