###################
#I/O
###################
library(Seurat)
library(tidyverse)

seurat_obj <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/mo_dias_15q_data/seurat_obj_alisa_caroline_15q_data_from_chunhui_UMB4643_sexID_fixed.RDS") #from PMID: 39079538

###################
#CREATE AGGREGATED CELL TYPES
##################
#aggregate cell types
map_celltype <- function(ct) {
  case_when(
    ct %in% c("Astrocyte I", "Astrocyte II") ~ "Astrocytes",
    ct == "Endothelial" ~ "Endothelial",
    str_detect(ct, "^Inh-") ~ "IN",
    ct == "Microglia" ~ "Microglia",
    str_detect(ct, "^Neu") ~ "Neurons",
    ct %in% c("OL", "Oligodendrocyte") ~ "Oligos",
    ct == "OPC" ~ "OPCs",
    ct == "Pericyte" ~ "Pericytes",
    ct == "T Cells" ~ "TCells",
    TRUE ~ "Other"
  )
}
seurat_obj$our_collapsed_broad_celltypes <- map_celltype(seurat_obj$celltype)
table(seurat_obj$celltype, seurat_obj$our_collapsed_broad_celltypes)
types <- unique(seurat_obj$our_collapsed_broad_celltypes)

###################
#CONFIRM CORRECT SEX LABELING
##################
pdf("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Revisions/Analysis_of_15q_and_Idiopathic_Autism_Data/Output/Plots/SuppFig_PanelA_SexLinkedGeneExpr.pdf", width = 8, height = 8)
DotPlot(seurat_obj, features = c("XIST", "DDX3Y", "ZFY", "USP9Y"), cluster.idents = TRUE)
dev.off()