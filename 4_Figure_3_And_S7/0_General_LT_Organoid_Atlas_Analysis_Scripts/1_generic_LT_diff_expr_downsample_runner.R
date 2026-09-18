#!/usr/bin/env Rscript
curr_guide = commandArgs(trailingOnly=TRUE)[1]
day = commandArgs(trailingOnly=TRUE)[2]
level = commandArgs(trailingOnly=TRUE)[3]
numDownsampleCells = as.numeric(commandArgs(trailingOnly=TRUE)[4]) #note: this is the number of cells per guide, not total (total number is numDownsampleCells * 2); pass in -1 if you want to use all cells
day_level = paste0(day, "_", level)
col_name = paste0("CellType_", level)

print("=======================================")
print(paste0("Working on guide ", curr_guide, " for ", day, " and level ", level, " (Downsampling to ", numDownsampleCells, " cells!)"))
print("=======================================")

library(Seurat)
library(tidyverse)
library(rlang)
options(Seurat.object.assay.version = "v3")

#creat output path
dir.create(paste0("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas/Output/DEG_Results_Downsampled/", day_level, "/", curr_guide, "/"))
dir.create(paste0("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas/Output/DEG_Results_Downsampled/", day_level, "/", curr_guide, "/DownsamplingTo", numDownsampleCells, "Cells/"))
fullPath <- paste0("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas/Output/DEG_Results_Downsampled/", day_level, "/", curr_guide, "/DownsamplingTo", numDownsampleCells, "Cells/")

#subset to only relevant day 
merged.Seurat.obj <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas/Output/Seurat_Objects/9_processed_seurat_obj_fullObj_integrated_harmony_ctrlOnlyAnnotations_qcFiltered_finalAnnots.RDS")
#subset to only relevant guides 
guides <- c(curr_guide, "NTC")
Idents(merged.Seurat.obj) <- "guide_assignment_NTC_merged"
merged.Seurat.obj <- subset(merged.Seurat.obj, idents = guides)
merged.Seurat.obj.og <- merged.Seurat.obj

metadata <- merged.Seurat.obj.og@meta.data 
if (numDownsampleCells != -1)
{
    ntc_cells <- sample(row.names(metadata %>% filter(guide_assignment_NTC_merged == "NTC")), replace = FALSE, numDownsampleCells)
}

for (i in seq(1:5)) #run 5 replicates 
{
    print("*********************")
    paste0("Replicate: ", i)
    print("*********************")

    if (numDownsampleCells != -1)
    {
        #downsample to number of cells indicated, maintaining equivalent proportions of cell types per guide 
        metadata_list <- list()
        props <- table(metadata$guide_assignment_NTC_merged, metadata[[col_name]])
        props <- props[which(row.names(props) != "NTC"),,drop = FALSE]
        props <- ceiling(props/sum(props) * numDownsampleCells)

        for (row_ind in seq_len(nrow(props)))
        {
            temp_guide <- row.names(props)[row_ind]

            for (col_ind in seq_len(ncol(props)))
            {
                curr_metadata <- merged.Seurat.obj.og@meta.data

                curr_type <- colnames(props)[col_ind]
                num_cells <- props[row_ind, col_ind]

                curr_metadata <- curr_metadata %>% filter(guide_assignment_NTC_merged == temp_guide) %>% filter(!!sym(col_name) == curr_type) %>% sample_n(size = min(num_cells, n()))
                metadata_list[[length(metadata_list) + 1]] <- curr_metadata
            }
        }

        new_metadata <- do.call(rbind, metadata_list)
        new_metadata <- rbind(new_metadata, metadata[ntc_cells,])
        merged.Seurat.obj <- merged.Seurat.obj.og[,which(colnames(merged.Seurat.obj.og) %in% row.names(new_metadata))]
    } else 
    {
        merged.Seurat.obj <- merged.Seurat.obj.og
    }
    
    Idents(merged.Seurat.obj) <- col_name
    path <- paste0(fullPath, "Replicate=", i, "/")
    dir.create(path)
    saveRDS(merged.Seurat.obj@meta.data, paste0(path, "metadata.RDS"))

    types <- merged.Seurat.obj@meta.data %>%
    filter(guide_assignment_NTC_merged == curr_guide) %>%
    pull(!!sym(col_name)) %>%
    unique()
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
}





