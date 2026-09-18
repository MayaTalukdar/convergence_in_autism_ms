#!/usr/bin/env Rscript
curr_guide = commandArgs(trailingOnly=TRUE)[1]

print("=======================================")
print(paste0("Working on guide ", curr_guide, "!"))
print("=======================================")

library(Seurat)
library(tidyverse)
library(rlang)
options(Seurat.object.assay.version = "v3")

#output path
path <- paste0("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas_acuteKD/Output/DEG_Results/", curr_guide, "/")
dir.create(path)

#subset to only relevant guide
merged.Seurat.obj <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas_acuteKD/Output/Seurat_Objects/5_processed_seurat_obj_fullObj_integrated_harmony_organoidAtlasAnnotations_qcFiltered_finalAnnots.RDS")
guides <- c(curr_guide, "NTC")
Idents(merged.Seurat.obj) <- "guide_assignment_NTC_merged"
merged.Seurat.obj <- subset(merged.Seurat.obj, idents = guides)

#set up differential expression 
types <- unname(merged.Seurat.obj@meta.data %>%
  filter(guide_assignment == curr_guide) %>%
  dplyr::select(CellType_L1) %>%
  unique() %>% unlist())
Idents(merged.Seurat.obj) <- "CellType_L1"
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
