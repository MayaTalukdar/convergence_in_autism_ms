#!/usr/bin/env Rscript
curr_guide = commandArgs(trailingOnly=TRUE)[1]
day = commandArgs(trailingOnly=TRUE)[2]
level = commandArgs(trailingOnly=TRUE)[3]
day_level = paste0(day, "_", level)
col_name = paste0("CellType_", level)

print("=======================================")
print(paste0("Working on guide ", curr_guide, " for ", day, " and level ", level))
print("=======================================")

library(Seurat)
library(tidyverse)
library(rlang)
options(Seurat.object.assay.version = "v3")

#output path
path <- paste0("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas/Output/DEG_Results/", day_level, "/", curr_guide, "/")
dir.create(path)

#subset to only relevant day 
merged.Seurat.obj <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas/Output/Seurat_Objects/9_processed_seurat_obj_fullObj_integrated_harmony_ctrlOnlyAnnotations_qcFiltered_finalAnnots.RDS") 
Idents(merged.Seurat.obj) <- "day"
merged.Seurat.obj <- subset(merged.Seurat.obj, idents = day)

#subset to only relevant guides 
guides <- c(curr_guide, "NTC")
Idents(merged.Seurat.obj) <- "guide_assignment_NTC_merged"
merged.Seurat.obj <- subset(merged.Seurat.obj, idents = guides)

#set up differential expression 
types <- merged.Seurat.obj@meta.data %>%
  filter(guide_assignment == curr_guide) %>%
  pull(!!sym(col_name)) %>%
  unique()
Idents(merged.Seurat.obj) <- col_name
saveRDS(merged.Seurat.obj@meta.data, paste0(path, "metadata.RDS"))
for (type in types)
{
    print("*********************************")
    print(paste0("Starting ", type, "!"))
    
    tryCatch({
        DEG_res <- FindMarkers(merged.Seurat.obj, ident.1 = curr_guide, ident.2 = "NTC", group.by = "guide_assignment_NTC_merged", subset.ident = type, min.pct = 0.10, logfc.threshold = 0)
        saveRDS(DEG_res, paste0(path, type, "_DEG_res.RDS"))
    }, error = function(e) {
        cat("Error occurred while processing", type, ":", conditionMessage(e), "\n")
    })
} 
