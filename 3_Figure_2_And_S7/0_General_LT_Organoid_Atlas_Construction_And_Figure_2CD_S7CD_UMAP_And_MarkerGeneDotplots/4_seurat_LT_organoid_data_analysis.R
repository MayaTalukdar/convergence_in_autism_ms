library(clusterProfiler)
library(cowplot)
library(DOSE)
library(dplyr)
library(enrichplot)
library(fgsea)
library(gplots)
library(ggplot2)
library(harmony)
library(limma)
library(lsa)
library(marray)
library(msigdbr)
library(org.Hs.eg.db)
library(patchwork)
library(pheatmap)
library(data.table)
library(reshape2)
library(RColorBrewer)
library(rjson)
library(R.utils)
library(scCustomize)
library(Seurat)
library(stringr)
library(tidyverse)
options(Seurat.object.assay.version = "v3")

#######################
#CREATE UNPROCESSED SEURAT OBJECT
#######################
list.of.Seurat.objs <- list()

#read in all 10X files as matrices 
cellbender_path <- "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/raw_data/day_60_organoids_scrna/Align_Data/3_RunCellBender/Output"
data.dirs <- list.files(cellbender_path, pattern = "LT")
for (dir in data.dirs)
{
  print(dir)
  #read in cellbender object
  raw_mat <- Read10X_h5(paste0(cellbender_path, "/", dir,  "/cellbender_filtered_seurat.h5"))
  expr <- raw_mat[['Gene Expression']]
  sgrnas <- raw_mat[['CRISPR Guide Capture']]

  #add sample id 
  colnames(expr) <- paste0(dir, "-", colnames(expr))
  colnames(sgrnas) <- paste0(dir, "-", colnames(sgrnas))

  #create the Seurat object
  seurat_obj <- CreateSeuratObject(
    counts = expr,
    min.cells = 0,
    min.features = 0,
    project = dir)

  #add the CRISPR Guide Capture assay
  seurat_obj[['Guides']] <- CreateAssayObject(counts = sgrnas)

  list.of.Seurat.objs[[dir]] <- seurat_obj
}

cellbender_path <- "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/raw_data/day_30_organoids_scrna/Align_Data/3_RunCellBender/Output"
data.dirs <- list.files(cellbender_path, pattern = "b4_d30")
for (dir in data.dirs)
{
  print(dir)
  #read in cellbender object
  raw_mat <- Read10X_h5(paste0(cellbender_path, "/", dir,  "/cellbender_filtered_seurat.h5"))
  expr <- raw_mat[['Gene Expression']]
  sgrnas <- raw_mat[['CRISPR Guide Capture']]

  #add sample id 
  colnames(expr) <- paste0(dir, "-", colnames(expr))
  colnames(sgrnas) <- paste0(dir, "-", colnames(sgrnas))

  #create the Seurat object
  seurat_obj <- CreateSeuratObject(
    counts = expr,
    min.cells = 0,
    min.features = 0,
    project = dir)

  #add the CRISPR Guide Capture assay
  seurat_obj[['Guides']] <- CreateAssayObject(counts = sgrnas)

  list.of.Seurat.objs[[dir]] <- seurat_obj
}

print("CREATED ALL SEURAT OBJECTS!")

#create merged seurat object 
merged.Seurat.obj <- merge(list.of.Seurat.objs[[1]], y = list.of.Seurat.objs[2:length(list.of.Seurat.objs)])
print("MERGED SEURAT OBJECT!")

#save merged seurat object 
saveRDS(merged.Seurat.obj, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas/Output/Seurat_Objects/0_unprocessed_Seurat_obj.RDS")
print("SAVE OBJECT!")

#####################
# DETERMINE DOUBLETS WITH SCRUBLET
######################
merged.Seurat.obj <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas/Output/Seurat_Objects/0_unprocessed_Seurat_obj.RDS")
merged.Seurat.obj$orig.ident <- sapply(row.names(merged.Seurat.obj@meta.data), function(x) ifelse(grepl("LT", x), paste0(strsplit(x, "-")[[1]][1:2], collapse = "-"), strsplit(x, "-")[[1]][1]))
lanes <- unique(merged.Seurat.obj@meta.data$orig.ident)
prop_doublets <- list()
num_doublets <- list()
num_singlets <- list()
cells_to_keep <- c()
for (lane in lanes)
{
    print(lane)
    if (grepl("LT", lane))
    {
      path <- "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/raw_data/day_60_organoids_scrna/Align_Data/2_PerformAlignment"
    } else 
    {
      path <- "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/raw_data/day_30_organoids_scrna/Align_Data/2_PerformAlignment"
    }

    #read in scrublet scores 
    scrublet_scores <- unname(unlist(read.csv(paste0("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas/Output/Scrublet/", lane, "_scrublet_predictions.csv"), header = TRUE)))
    wd <- paste0(path, "/", lane)
    names(scrublet_scores) <- unname(unlist(read.table(paste0(wd, "/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"))))
    prop_doublets[[lane]] <- length(which(scrublet_scores == "True"))/length(scrublet_scores)
    num_doublets[[lane]] <- length(which(scrublet_scores == "True"))
    num_singlets[[lane]] <- length(which(scrublet_scores == "False"))

    #get cells to keep 
    current_cells_to_keep <- sapply(names(scrublet_scores)[which(scrublet_scores == "False")], function(x) paste0(lane, "-", x))
    cells_to_keep <- c(cells_to_keep, current_cells_to_keep)
}

#make sure that we have cells to retain from all samples 
length(intersect(cells_to_keep, colnames(merged.Seurat.obj)))/ncol(merged.Seurat.obj) #~95% cells are retained

#get all putative doublets
cells_to_keep <- intersect(unname(cells_to_keep), colnames(merged.Seurat.obj))
saveRDS(cells_to_keep, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas/Output/Scrublet/cells_not_identified_as_doublets_by_scrublet.RDS")

#add this metadata
merged.Seurat.obj$is_putative_doublet <- (!row.names(merged.Seurat.obj@meta.data) %in% cells_to_keep)
saveRDS(merged.Seurat.obj, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas/Output/Seurat_Objects/1_processed_seurat_obj_with_scrublet_info.RDS")

###################
#GUIDE ASSIGNMENT 
###################
merged.Seurat.obj <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas/Output/Seurat_Objects/1_processed_seurat_obj_with_scrublet_info.RDS")
putative_doublet_vec <- setNames(merged.Seurat.obj@meta.data$is_putative_doublet, row.names(merged.Seurat.obj@meta.data))
dim(merged.Seurat.obj)

#read in 10X assignments 
lanes <- unique(merged.Seurat.obj@meta.data$orig.ident)
guide_df_10x <- list()
for (lane in lanes)
{
  print(lane)
  if (grepl("LT", lane))
  {
    guide_df <- read.csv(paste0("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/raw_data/day_60_organoids_scrna/Align_Data/2_PerformAlignment/", lane, "/outs/crispr_analysis/protospacer_calls_per_cell.csv"))
    guide_df$cell_barcode <- paste0(lane, "-", guide_df$cell_barcode)
    guide_df_10x[[lane]] <- guide_df
  } else if (grepl("b4", lane))
  {
    guide_df <- read.csv(paste0("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/raw_data/day_30_organoids_scrna/Align_Data/2_PerformAlignment/", lane, "/outs/crispr_analysis/protospacer_calls_per_cell.csv"))
    guide_df$cell_barcode <- paste0(lane, "-", guide_df$cell_barcode)
    guide_df_10x[[lane]] <- guide_df
  }
}
guide_df_10x <- do.call(rbind, guide_df_10x)
row.names(guide_df_10x) <- NULL

##get final assignments for multiple guides 
multiple_guide_assignments <- setNames(rep("No_Assigned_Guide", nrow(summary_table_multiple_guides)), summary_table_multiple_guides$cell_barcode)
multiple_guide_assignments[which(summary_table_multiple_guides$rescued_status == "Multiple: Rescued")] <- summary_table_multiple_guides$top_guide[which(summary_table_multiple_guides$rescued_status == "Multiple: Rescued")]
guide_assignment_vec <- c(guide_assignment_vec, multiple_guide_assignments)
num_guide_vec <- c(num_guide_vec, setNames(summary_table_multiple_guides$rescued_status, summary_table_multiple_guides$cell_barcode))

#add info to seurat object
merged.Seurat.obj@meta.data$guide_assignment <- guide_assignment_vec[row.names(merged.Seurat.obj@meta.data)]
merged.Seurat.obj@meta.data$num_guide_class <- num_guide_vec[row.names(merged.Seurat.obj@meta.data)]
saveRDS(merged.Seurat.obj@meta.data, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas/Output/Guide_Assignment/guide_assignment_metadata.RDS")
saveRDS(merged.Seurat.obj, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas/Output/Seurat_Objects/2_processed_seurat_obj_with_full_guide_info.RDS")

###################
#QUALITY CONTROL
###################
merged.Seurat.obj <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas/Output/Seurat_Objects/2_processed_seurat_obj_with_full_guide_info.RDS")
merged.Seurat.obj[["percent.mt"]] <- PercentageFeatureSet(merged.Seurat.obj, pattern = "^MT-")

#look at summary of key qc metrics based on predicted doublet status 
metadata <- merged.Seurat.obj@meta.data
metadata %>%
  group_by(is_putative_doublet) %>%
  summarize(
    Q1_nFeature = quantile(nFeature_RNA, 0.25),
    med_nFeature = median(nFeature_RNA),
    Q3_nFeature = quantile(nFeature_RNA, 0.75),
    Q1_nUMI = quantile(nCount_RNA, 0.25),
    med_nUMI = median(nCount_RNA),
    Q3_nUMI = quantile(nCount_RNA, 0.75)
  ) %>% as.data.frame() %>% t()

#filter out cells we did not assign a guide to 
Idents(merged.Seurat.obj) <- "num_guide_class"
merged.Seurat.obj <- subset(merged.Seurat.obj, idents = "Single_Guide")

#filter out putative doublets
Idents(merged.Seurat.obj) <- "is_putative_doublet"
merged.Seurat.obj <- subset(merged.Seurat.obj, idents = FALSE)

#filter cells 
merged.Seurat.obj<- subset(merged.Seurat.obj, subset = percent.mt < 20)
merged.Seurat.obj <- subset(merged.Seurat.obj, subset = nFeature_RNA > 500)
merged.Seurat.obj <- subset(merged.Seurat.obj, subset = nFeature_RNA < 6000)
saveRDS(merged.Seurat.obj, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas/Output/Seurat_Objects/3_processed_seurat_obj_filtered.RDS")
write.table(merged.Seurat.obj@meta.data, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas/Output/Seurat_Objects/3_processed_seurat_obj_filtered_metadata.csv", col.names = TRUE, row.names = TRUE, sep = ",", quote = FALSE)

#####################################
#NORMALIZATION & DIMENSIONALITY REDUCTION
##normalize and residualize cells 
#####################################
#convert to seurat assay v5
merged.Seurat.obj[["RNA5"]] <- as(object = merged.Seurat.obj[["RNA"]], Class = "Assay5")

#normalize data
DefaultAssay(merged.Seurat.obj) <- "RNA5"
Idents(merged.Seurat.obj) <- "orig.ident"
merged.Seurat.obj[["RNA5"]] <- split(merged.Seurat.obj[["RNA5"]], f = merged.Seurat.obj$orig.ident)
merged.Seurat.obj<- NormalizeData(merged.Seurat.obj)
merged.Seurat.obj <- FindVariableFeatures(merged.Seurat.obj, selection.method = "vst", nfeatures = 3000)
merged.Seurat.obj <- ScaleData(merged.Seurat.obj, vars.to.regress = c("percent_mito", "percent_ribo", "nCount_RNA", "nFeature_RNA"))

#perform dimensionality reduction
merged.Seurat.obj <- RunPCA(merged.Seurat.obj)
merged.Seurat.obj <- ProjectDim(object = merged.Seurat.obj)
merged.Seurat.obj <- FindNeighbors(merged.Seurat.obj, dims = 1:30, reduction = "pca")
merged.Seurat.obj <- FindClusters(merged.Seurat.obj, resolution = 1.2, cluster.name = "unintegrated_clusters.1.2")
merged.Seurat.obj <- FindClusters(merged.Seurat.obj, resolution = 0.6, cluster.name = "unintegrated_clusters.0.6")
merged.Seurat.obj <- RunUMAP(merged.Seurat.obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

#look at unintegrated cell types
DefaultAssay(merged.Seurat.obj) <- "RNA5"
pdf("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas/Output/UMAPs/All/unintegrated_umap.pdf", width = 6)
DimPlot(merged.Seurat.obj, reduction = "umap.unintegrated", group.by = "orig.ident", label = TRUE) + NoLegend()
dev.off()

#####################################
#INTEGRATE
#####################################
#integrate
merged.Seurat.obj <- IntegrateLayers(
 object = merged.Seurat.obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = TRUE
 )
print("Completed harmony!")

#cluster
merged.Seurat.obj <- JoinLayers(merged.Seurat.obj)
merged.Seurat.obj <- FindNeighbors(merged.Seurat.obj, reduction = "harmony", dims = 1:30)
merged.Seurat.obj <- FindClusters(merged.Seurat.obj, resolution = 1.2, cluster.name = "harmony_clusters_1.2")
merged.Seurat.obj <- FindClusters(merged.Seurat.obj, resolution = 0.6, cluster.name = "harmony_clusters_0.6")
merged.Seurat.obj <- RunUMAP(merged.Seurat.obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")

saveRDS(merged.Seurat.obj, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas/Output/Seurat_Objects/6_processed_seurat_obj_fullObj_integrated_harmony.RDS")

#####################################
#CELL TYPE LABEL TRANSFER FROM CONTROL ONLY OBJECT
#####################################
#run transfer
merged.Seurat.obj <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas/Output/Seurat_Objects/6_processed_seurat_obj_fullObj_integrated_harmony.RDS")

ref <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas/Output/Seurat_Objects/5_processed_seurat_obj_controlsOnly_integrated_harmony_xuyuAnnotations_finalAnnotations.RDS")
anchors <- FindTransferAnchors(reference = ref, query = merged.Seurat.obj, dims = 1:30,
    reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = ref$l2, dims = 1:30)
merged.Seurat.obj <- AddMetaData(merged.Seurat.obj, metadata = predictions)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
merged.Seurat.obj <- CellCycleScoring(merged.Seurat.obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

saveRDS(merged.Seurat.obj, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas/Output/Seurat_Objects/7_processed_seurat_obj_fullObj_integrated_harmony_ctrlOnlyAnnotations.RDS")

#####################################
#EXPLORATORY ANALYSIS 
#####################################
merged.Seurat.obj <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas/Output/Seurat_Objects/7_processed_seurat_obj_fullObj_integrated_harmony_ctrlOnlyAnnotations.RDS")
Idents(merged.Seurat.obj) <- "harmony_clusters_1.2"

#add in annotations from control only object 
ref <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas/Output/Seurat_Objects/5_processed_seurat_obj_controlsOnly_integrated_harmony_xuyuAnnotations_finalAnnotations.RDS")
control_annots <- setNames(ref$l2, row.names(ref@meta.data))
merged.Seurat.obj$control_annots <- setNames(control_annots[row.names(merged.Seurat.obj@meta.data)], row.names(merged.Seurat.obj@meta.data))
merged.Seurat.obj$control_annots[which(is.na(merged.Seurat.obj$control_annots))] <- "NotInControlObject"
merged.Seurat.obj$day <-sapply(merged.Seurat.obj$orig.ident, function(x) ifelse(grepl("LT", x), "day60", "day30"))

#look at quality control per cluster and remove low qual clusters 
merged.Seurat.obj@meta.data %>% 
  group_by(harmony_clusters_1.2) %>% 
  summarize(med = median(nFeature_RNA)) %>% arrange(desc(med)) %>% as.data.frame()

merged.Seurat.obj@meta.data %>% 
  group_by(harmony_clusters_1.2) %>% 
  summarize(med = median(percent.mt)) %>% arrange(desc(med)) %>% as.data.frame()

low_qual_clusters <- c(21, 28, 29, 31, 32)
Idents(merged.Seurat.obj) <- "harmony_clusters_1.2"
merged.Seurat.obj <- subset(merged.Seurat.obj, idents = low_qual_clusters, invert = TRUE)

saveRDS(merged.Seurat.obj, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas/Output/Seurat_Objects/8_processed_seurat_obj_fullObj_integrated_harmony_ctrlOnlyAnnotations_qcFiltered.RDS")

#####################################
#ADD IN FINAL ANNOTATIONS
#####################################
merged.Seurat.obj <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas/Output/Seurat_Objects/8_processed_seurat_obj_fullObj_integrated_harmony_ctrlOnlyAnnotations_qcFiltered.RDS")
annots_from_becky <- fread("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas/full_metadata_updated_CellType_final.txt", header = TRUE) %>% as.data.frame() #In Input_Files_Not_Generated_By_Scripts
row.names(annots_from_becky) <- annots_from_becky$rownames
annots_from_becky$rownames <- NULL
all.equal(row.names(annots_from_becky), row.names(merged.Seurat.obj@meta.data))
my_levels <- c("RG_Div","RG_nonDiv","IPC_Div","IPC_nonDiv","ExN_Immature","ExN_Mature","ExN_SP","IN")
Idents(merged.Seurat.obj) <- "CellType_L2"
Idents(merged.Seurat.obj) <- factor(Idents(merged.Seurat.obj), levels= my_levels)
merged.Seurat.obj@meta.data <- annots_from_becky
merged.Seurat.obj$guide_assignment_NTC_merged <- merged.Seurat.obj$guide_assignment
merged.Seurat.obj$guide_assignment_NTC_merged[which(grepl("NTC", merged.Seurat.obj$guide_assignment_NTC_merged))] <- "NTC"
merged.Seurat.obj$cleaned_guide_assignment <- gsub("_KD", "", merged.Seurat.obj$guide_assignment_NTC_merged)
Idents(merged.Seurat.obj) <- "cleaned_guide_assignment"
subset_merged.Seurat.obj <- subset(merged.Seurat.obj, idents = "NTC")
Idents(subset_merged.Seurat.obj) <- "CellType_L2"
Idents(subset_merged.Seurat.obj) <- factor(Idents(subset_merged.Seurat.obj), levels= my_levels)

#plot feature plot
markers <- c("MKI67", "ASPM", "HES1", "SOX2", "EOMES", "NHLH1", "NEUROD1", "SLC17A6", "NEUROD2", "SLC17A7", "TBR1", "BCL11B", "RORB", "SATB2", "RELN", "CALB2", "DLX2", "GAD2")
pdf("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Figure_2/MarkerGeneDotplot.pdf", width = 12, height = 5)
merged.Seurat.obj_dotplot <- DotPlot(subset_merged.Seurat.obj, features = markers) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(merged.Seurat.obj_dotplot + ggtitle("Marker Gene Expression"))
dev.off()

saveRDS(merged.Seurat.obj, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas/Output/Seurat_Objects/9_processed_seurat_obj_fullObj_integrated_harmony_ctrlOnlyAnnotations_qcFiltered_finalAnnots.RDS")

#####################################
#PLOT UMAPS
#####################################
merged.Seurat.obj <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas/Output/Seurat_Objects/9_processed_seurat_obj_fullObj_integrated_harmony_ctrlOnlyAnnotations_qcFiltered_finalAnnots.RDS")
merged.Seurat.obj$New_CellType_L2 <- gsub("_", "-", merged.Seurat.obj$CellType_L2)
new_cellType_L2_mapper <- c(
  "ExN-Immature" = "ExN-Imm", 
  "ExN-SP" = "ExN-Early", 
  "ExN-Mature" = "ExN-Mat"
)
merged.Seurat.obj$New_CellType_L2 <- sapply(merged.Seurat.obj$New_CellType_L2, function(x) ifelse(x %in% names(new_cellType_L2_mapper), new_cellType_L2_mapper[x], x))

type_color_vector <- c("IN" = "#6f3bbb", 
"ExN-Early" = "#1f63b4", 
"ExN-Mat" = "#12a2c4", 
"ExN-Imm" = "#78a641", 
"IPC-nonDiv" = "#ffd750", 
"IPC-Div" = "#ff7f0e", 
"RG-nonDiv" = "#d63a3a", 
"RG-Div" = "#ff8196")

guide_color_vector <- c(
  "NTC" = "#b4b4b4",
  "BAZ2B" = "#e75d62",
  "CLASP1" = "#e5a506",
  "EHMT1" = "#dedf61",
  "NR2F1-AS1" = "#217885",
  "PPP3CA" = "#4676b7",
  "ST7" = "#273e6b",
  "WDFY3" = "#93418d"
)

#subset to only new guides of interest (modified 7/2/25)
sgRNAs <- c("NTC", "BAZ2B", "CLASP1", "EHMT1", "NR2F1-AS1", "PPP3CA",  "ST7", "WDFY3")
merged.Seurat.obj$cleaned_guide_assignment <- gsub("_KD", "", merged.Seurat.obj$guide_assignment_NTC_merged)
Idents(merged.Seurat.obj) <- "cleaned_guide_assignment"
merged.Seurat.obj <- subset(merged.Seurat.obj, idents = sgRNAs)

#l2 
Idents(merged.Seurat.obj) <- "New_CellType_L2"

pdf("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Figure_2/5C_L2_UMAP.pdf", width = 4, height = 2.8)
DimPlot_scCustom(
  merged.Seurat.obj, 
  reduction = "umap.harmony", 
  colors_use = type_color_vector, 
  pt.size = 0.1, 
  shuffle = FALSE, 
  label = TRUE, 
  repel = FALSE) + ggtitle("Cell Type") + 
  NoLegend() +
  labs(x = "UMAP-1", y = "UMAP-2") +
  theme(
    axis.text = element_blank(),  # Remove tick labels
    axis.ticks = element_blank(), # Remove ticks
    axis.title = element_text(face = "bold")
  )
dev.off()

#day
Idents(merged.Seurat.obj) <- "day"

pdf("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Figure_2/S13C_day_UMAP_NoLegend.pdf", width = 4, height = 2.8)
DimPlot_scCustom(
  merged.Seurat.obj, 
  reduction = "umap.harmony", 
  pt.size = 0.1, ggplot_default_colors = TRUE,
  shuffle = TRUE, label = FALSE) + 
  NoLegend() + ggtitle("Timepoint") + 
  labs(x = "UMAP-1", y = "UMAP-2") +
  theme(
    axis.text = element_blank(),  # Remove tick labels
    axis.ticks = element_blank(), # Remove ticks
    axis.title = element_text(face = "bold")
  )
dev.off()

pdf("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Figure_2/S13C_day_UMAP_WithLegend.pdf", width = 4, height = 2.8)
DimPlot_scCustom(
  merged.Seurat.obj, 
  reduction = "umap.harmony", 
  pt.size = 0.1, ggplot_default_colors = TRUE,
  shuffle = TRUE, label = FALSE) + 
  labs(x = "UMAP-1", y = "UMAP-2") + ggtitle("Timepoint") +
  theme(
    axis.text = element_blank(),  # Remove tick labels
    axis.ticks = element_blank(), # Remove ticks
    axis.title = element_text(face = "bold")
  )
dev.off()

#guides
merged.Seurat.obj$cleaned_guide_assignment <- gsub("_KD", "", merged.Seurat.obj$guide_assignment_NTC_merged)
Idents(merged.Seurat.obj) <- "cleaned_guide_assignment"
levels(Idents(merged.Seurat.obj)) <- c("NTC", setdiff(sort(as.character(unique(Idents(merged.Seurat.obj)))), "NTC"))

pdf("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Figure_2/S13D_guide_UMAP_NoLegend.pdf", width = 4, height = 2.8)
DimPlot_scCustom(
  merged.Seurat.obj, 
  reduction = "umap.harmony", 
  pt.size = 0.1, colors_use = guide_color_vector,
  shuffle = TRUE, label = FALSE) + 
  NoLegend() + ggtitle("Guide") + 
  labs(x = "UMAP-1", y = "UMAP-2") +
  theme(
    axis.text = element_blank(),  # Remove tick labels
    axis.ticks = element_blank(), # Remove ticks
    axis.title = element_text(face = "bold")
  )
dev.off()

pdf("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Figure_2/S13D_guide_UMAP_WithLegend.pdf", width = 4, height = 2.8)
DimPlot_scCustom(
  merged.Seurat.obj, 
  reduction = "umap.harmony", 
  pt.size = 0.1, colors_use = guide_color_vector,
  shuffle = TRUE, label = FALSE) + 
  labs(x = "UMAP-1", y = "UMAP-2") + ggtitle("Guide") +
  theme(
    axis.text = element_blank(),  # Remove tick labels
    axis.ticks = element_blank(), # Remove ticks
    axis.title = element_text(face = "bold")
  )
dev.off()

