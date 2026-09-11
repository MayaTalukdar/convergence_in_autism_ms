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
data.dirs <- list.files(cellbender_path, pattern = "AD2")
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
saveRDS(merged.Seurat.obj, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas_acuteKD/Output/Seurat_Objects/0_unprocessed_Seurat_obj.RDS")
print("SAVE OBJECT!")

#####################
# DETERMINE DOUBLETS WITH SCRUBLET
######################
merged.Seurat.obj <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas_acuteKD/Output/Seurat_Objects/0_unprocessed_Seurat_obj.RDS")
merged.Seurat.obj$orig.ident <- sapply(row.names(merged.Seurat.obj@meta.data), function(x) ifelse(grepl("AD2", x), paste0(strsplit(x, "-")[[1]][1:2], collapse = "-"), strsplit(x, "-")[[1]][1]))
lanes <- unique(merged.Seurat.obj@meta.data$orig.ident)
prop_doublets <- list()
num_doublets <- list()
num_singlets <- list()
cells_to_keep <- c()
for (lane in lanes)
{
    print(lane)

    path <- "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/raw_data/day_60_organoids_scrna/Align_Data/2_PerformAlignment"

    #read in scrublet scores 
    scrublet_scores <- unname(unlist(read.csv(paste0("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas_acuteKD/Output/Scrublet/", lane, "_scrublet_predictions.csv"), header = TRUE)))
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
# [1] 0.9402987

#get all putative doublets
cells_to_keep <- intersect(unname(cells_to_keep), colnames(merged.Seurat.obj))
saveRDS(cells_to_keep, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas_acuteKD/Output/Scrublet/cells_not_identified_as_doublets_by_scrublet.RDS")

#add this metadata
merged.Seurat.obj$is_putative_doublet <- (!row.names(merged.Seurat.obj@meta.data) %in% cells_to_keep)
saveRDS(merged.Seurat.obj, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas_acuteKD/Output/Seurat_Objects/1_processed_seurat_obj_with_scrublet_info.RDS")

###################
#GUIDE ASSIGNMENT 
###################
merged.Seurat.obj <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas_acuteKD/Output/Seurat_Objects/1_processed_seurat_obj_with_scrublet_info.RDS")
putative_doublet_vec <- setNames(merged.Seurat.obj@meta.data$is_putative_doublet, row.names(merged.Seurat.obj@meta.data))
dim(merged.Seurat.obj)
# [1]  39107 250430

#read in 10X assignments 
lanes <- unique(merged.Seurat.obj@meta.data$orig.ident)
guide_df_10x <- list()
for (lane in lanes)
{
  guide_df <- read.csv(paste0("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/raw_data/day_60_organoids_scrna/Align_Data/2_PerformAlignment/", lane, "/outs/crispr_analysis/protospacer_calls_per_cell.csv"))
  guide_df$cell_barcode <- paste0(lane, "-", guide_df$cell_barcode)
  guide_df_10x[[lane]] <- guide_df
}
guide_df_10x <- do.call(rbind, guide_df_10x)
row.names(guide_df_10x) <- NULL
table(colnames(merged.Seurat.obj) %in% guide_df_10x$cell_barcode)
#  FALSE   TRUE 
# 144920 105510 
table(guide_df_10x$cell_barcode %in% colnames(merged.Seurat.obj))
# FALSE   TRUE 
#   1604 105510 

#identify cells with no assigned guides 
num_guide_vec <- setNames(sapply(setdiff(colnames(merged.Seurat.obj), guide_df_10x$cell_barcode), function(x) "No_Assigned_Guide"), setdiff(colnames(merged.Seurat.obj), guide_df_10x$cell_barcode))
guide_assignment_vec <- num_guide_vec

#identify cells with 1 assigned guide 
guide_df_single_guide <- guide_df_10x %>% filter(num_features == 1)
num_guide_vec <- c(num_guide_vec, setNames(sapply(guide_df_single_guide$feature_call, function(x) "Single_Guide"), guide_df_single_guide$cell_barcode))
guide_assignment_vec <- c(guide_assignment_vec, setNames(guide_df_single_guide$feature_call, guide_df_single_guide$cell_barcode))

#analyze cells that received multiple guides
guide_df_multiple_guides <- guide_df_10x %>% filter(num_features > 1)

##create summary table
create_summary_table_row <- function(row) 
{
  umi_values <- as.integer(unlist(str_split(row["num_umis"], "\\|")))
  guide_list <- unlist(str_split(row["feature_call"], "\\|"))
  
  max_umi <- max(umi_values)
  second_max_umi <- sort(umi_values, decreasing = TRUE)[2]
  
  ratio_of_top_umis <- max_umi / second_max_umi
  
  top_guide <- guide_list[which.max(umi_values)]
  
  return(c(row, max_umi = max_umi, second_max_umi = second_max_umi, 
           ratio_of_top_umis = ratio_of_top_umis, top_guide = top_guide))
}

summary_table_multiple_guides <- apply(guide_df_multiple_guides, 1, function(x) create_summary_table_row(x)) %>% as.data.frame() %>% t() %>% as.data.frame()
row.names(summary_table_multiple_guides) <- NULL
summary_table_multiple_guides$rescued_status <- 'Multiple: Rescued'
summary_table_multiple_guides$is_putative_doublet <- putative_doublet_vec[summary_table_multiple_guides$cell]

##any guides with a ratio of less than 4 should be not assigned 
summary_table_multiple_guides$ratio_of_top_umis <- as.numeric(summary_table_multiple_guides$ratio_of_top_umis)
summary(summary_table_multiple_guides$ratio_of_top_umis) 
summary_table_multiple_guides$rescued_status[which(summary_table_multiple_guides$ratio_of_top_umis < 4)] <- 'Muliple: Not Rescued - Failed Ratio Test'

##any guides in which the second max umi is greater than should be not assigned 
summary_table_multiple_guides$second_max_umi <- as.numeric(summary_table_multiple_guides$second_max_umi)
summary((summary_table_multiple_guides %>% filter(rescued_status != 'Not Rescued - Failed Ratio Test'))$second_max_umi) 
summary_table_multiple_guides$rescued_status[which(summary_table_multiple_guides$ratio_of_top_umis >= 5)] <- 'Multiple: Not Rescued - High Alt Guide UMI Count'

##get final assignments for multiple guides 
multiple_guide_assignments <- setNames(rep("No_Assigned_Guide", nrow(summary_table_multiple_guides)), summary_table_multiple_guides$cell_barcode)
multiple_guide_assignments[which(summary_table_multiple_guides$rescued_status == "Multiple: Rescued")] <- summary_table_multiple_guides$top_guide[which(summary_table_multiple_guides$rescued_status == "Multiple: Rescued")]
guide_assignment_vec <- c(guide_assignment_vec, multiple_guide_assignments)
num_guide_vec <- c(num_guide_vec, setNames(summary_table_multiple_guides$rescued_status, summary_table_multiple_guides$cell_barcode))

#add info to seurat object
merged.Seurat.obj@meta.data$guide_assignment <- guide_assignment_vec[row.names(merged.Seurat.obj@meta.data)]
merged.Seurat.obj@meta.data$num_guide_class <- num_guide_vec[row.names(merged.Seurat.obj@meta.data)]
saveRDS(merged.Seurat.obj@meta.data, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas_acuteKD/Output/Guide_Assignment/guide_assignment_metadata.RDS")

#examine qc metrics by these different assignments 
Idents(merged.Seurat.obj) <- "num_guide_class"

##doublets
meta_data <- merged.Seurat.obj@meta.data

proportion_by_guide <- meta_data %>%
  group_by(num_guide_class, is_putative_doublet) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(num_guide_class) %>%
  mutate(proportion = count / sum(count))

overall_proportion <- meta_data %>%
  group_by(is_putative_doublet) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(proportion = count / sum(count),
         num_guide_class = "Overall") # Label for overall bar

plot_data <- bind_rows(proportion_by_guide, overall_proportion)
full_names <- c(
  "Muliple: Not Rescued - Failed Ratio Test",
  "Multiple: Not Rescued - High Alt Guide UMI Count",
  "Multiple: Rescued",
  "No_Assigned_Guide",
  "Single_Guide",
  "Overall"
)
plot_data$num_guide_class <- factor(plot_data$num_guide_class, levels = full_names)

p <- ggplot(plot_data, aes(x = num_guide_class, y = proportion, fill = is_putative_doublet)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = "Guide Assignment",
    y = "Proportion",
    fill = "Putative Doublet",
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.minor = element_blank()
  )

##plot
pdf("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas_acuteKD/Output/Guide_Assignment/guide_assignment_doublet.pdf", width = 10, height = 8)
print(p)
dev.off()

saveRDS(merged.Seurat.obj, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas_acuteKD/Output/Seurat_Objects/2_processed_seurat_obj_with_full_guide_info.RDS")

###################
#QUALITY CONTROL
###################
merged.Seurat.obj <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas_acuteKD/Output/Seurat_Objects/2_processed_seurat_obj_with_full_guide_info.RDS")
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

table(merged.Seurat.obj$is_putative_doublet)
table(merged.Seurat.obj$num_guide_class)

#determine cells that will be lost by each qc filtering step 
cells_filtered_pct_mito <- row.names(merged.Seurat.obj@meta.data %>% filter(percent.mt >= 20))
cells_filtered_nFeatures <- row.names(merged.Seurat.obj@meta.data %>% filter(nFeature_RNA <= 500|nFeature_RNA >= 6000)) #3rd quartile of doublets 
filtering_list <- list("pct_mito" = cells_filtered_pct_mito, "nFeatures" = cells_filtered_nFeatures)
print("NUMBER OF CELLS FILTERED PER CRITERIA: ")
print(lapply(filtering_list, length))
print(length(unique(do.call(c, filtering_list))))
print(length(unique(do.call(c, filtering_list)))/ncol(merged.Seurat.obj))

#filter cells 
merged.Seurat.obj<- subset(merged.Seurat.obj, subset = percent.mt < 20)
merged.Seurat.obj <- subset(merged.Seurat.obj, subset = nFeature_RNA > 500)
merged.Seurat.obj <- subset(merged.Seurat.obj, subset = nFeature_RNA < 6000)
saveRDS(merged.Seurat.obj, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas_acuteKD/Output/Seurat_Objects/3_processed_seurat_obj_filtered.RDS")
write.table(merged.Seurat.obj@meta.data, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas_acuteKD/Output/Seurat_Objects/3_processed_seurat_obj_filtered_metadata.csv", col.names = TRUE, row.names = TRUE, sep = ",", quote = FALSE)

source("/n/groups/walsh/indData/Maya/microCHIP_AD_project/GOT/FCD_Analysis/Deprecated_For_Second_Submission/Final_Analysis_For_First_Submission/util.R")
ref <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas/Output/Seurat_Objects/9_processed_seurat_obj_fullObj_integrated_harmony_ctrlOnlyAnnotations_qcFiltered_finalAnnots.RDS")
Idents(ref) <- "guide_assignment_NTC_merged"
ref <- subset(ref, idents = "NTC")
Idents(ref) <- "day"
ref <- subset(ref, idents = "day60")

merged.Seurat.obj <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas_acuteKD/Output/Seurat_Objects/3_processed_seurat_obj_filtered.RDS")

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
merged.Seurat.obj <- ScaleData(merged.Seurat.obj, vars.to.regress = c("percent.mt", "percent_ribo", "nCount_RNA", "nFeature_RNA"))

#perform dimensionality reduction
merged.Seurat.obj <- RunPCA(merged.Seurat.obj)
merged.Seurat.obj <- ProjectDim(object = merged.Seurat.obj)
merged.Seurat.obj <- FindNeighbors(merged.Seurat.obj, dims = 1:30, reduction = "pca")
merged.Seurat.obj <- FindClusters(merged.Seurat.obj, resolution = 1.2, cluster.name = "unintegrated_clusters.1.2")
merged.Seurat.obj <- FindClusters(merged.Seurat.obj, resolution = 0.6, cluster.name = "unintegrated_clusters.0.6")
merged.Seurat.obj <- RunUMAP(merged.Seurat.obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

#look at unintegrated cell types
DefaultAssay(merged.Seurat.obj) <- "RNA5"
pdf("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas_acuteKD/Output/UMAPs/unintegrated_umap.pdf", width = 6)
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

#####################################
#CELL TYPE LABEL TRANSFER FROM OUR ORGANOID ATLAS
#####################################
anchors <- FindTransferAnchors(reference = ref, query = merged.Seurat.obj, dims = 1:30,
    reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = ref$CellType_L2, dims = 1:30)
merged.Seurat.obj <- AddMetaData(merged.Seurat.obj, metadata = predictions)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
merged.Seurat.obj <- CellCycleScoring(merged.Seurat.obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

saveRDS(merged.Seurat.obj, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas_acuteKD/Output/Seurat_Objects/4_processed_seurat_obj_fullObj_integrated_harmony_organoidAtlasAnnotations.RDS")

#####################################
#ADD FINAL ANNOTATIONS
#####################################
becky_metadata <- readRDS("/n/groups/walsh/indData/becky/organoids_acuteKD/AD2_metadata_updated_CellType_final.rds")
dim(becky_metadata)
dim(merged.Seurat.obj)
becky_metadata <- becky_metadata[row.names(merged.Seurat.obj@meta.data),]
all.equal(row.names(becky_metadata), row.names(merged.Seurat.obj@meta.data))
merged.Seurat.obj@meta.data <- becky_metadata

#add in a column to merge ntcs
merged.Seurat.obj$guide_assignment_NTC_merged <- merged.Seurat.obj$guide_assignment
merged.Seurat.obj$guide_assignment_NTC_merged[which(grepl("NTC", merged.Seurat.obj$guide_assignment_NTC_merged))] <- "NTC"
saveRDS(merged.Seurat.obj, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas_acuteKD/Output/Seurat_Objects/5_processed_seurat_obj_fullObj_integrated_harmony_organoidAtlasAnnotations_qcFiltered_finalAnnots.RDS")
merged.Seurat.obj <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas_acuteKD/Output/Seurat_Objects/5_processed_seurat_obj_fullObj_integrated_harmony_organoidAtlasAnnotations_qcFiltered_finalAnnots.RDS")

merged.Seurat.obj$cleaned_guide_assignment <- gsub("_KD", "", merged.Seurat.obj$guide_assignment_NTC_merged)
Idents(merged.Seurat.obj) <- "cleaned_guide_assignment"
subset_merged.Seurat.obj <- subset(merged.Seurat.obj, idents = "NTC")
new_cellType_L2_mapper <- c(
  "ExN-Immature" = "ExN-Imm", 
  "ExN-SP" = "ExN-Early", 
  "ExN-Mature" = "ExN-Mat"
)
subset_merged.Seurat.obj$New_CellType_L2 <- sapply(subset_merged.Seurat.obj$CellType_L2, function(x) ifelse(x %in% names(new_cellType_L2_mapper), new_cellType_L2_mapper[x], x))
my_levels <- c("RG_Div","RG_nonDiv","IPC_Div","IPC_nonDiv","ExN_Immature","ExN_Mature","ExN_SP","IN")
Idents(subset_merged.Seurat.obj) <- "New_CellType_L2"
Idents(subset_merged.Seurat.obj) <- factor(Idents(subset_merged.Seurat.obj), levels= my_levels)

#plot feature plot
markers <- c("MKI67", "ASPM", "HES1", "SOX2", "EOMES", "NHLH1", "NEUROD1", "SLC17A6", "NEUROD2", "SLC17A7", "TBR1", "BCL11B", "RORB", "SATB2", "RELN", "CALB2", "DLX2", "GAD2")
pdf("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Figure_3/MarkerGeneDotplot.pdf", width = 12, height = 5)
merged.Seurat.obj_dotplot <- DotPlot(subset_merged.Seurat.obj, features = markers) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(merged.Seurat.obj_dotplot + ggtitle("Marker Gene Expression"))
dev.off()

####################
#CREATE UMAPS
####################
merged.Seurat.obj$cleaned_guide_assignment <- gsub("_KD", "", merged.Seurat.obj$guide_assignment_NTC_merged)
Idents(merged.Seurat.obj) <- "cleaned_guide_assignment"
levels(Idents(merged.Seurat.obj)) <- c("NTC", setdiff(sort(as.character(unique(Idents(merged.Seurat.obj)))), "NTC"))

pdf("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/figures_for_manuscript/generate_fig_7b_s15a_s15b_umaps_and_grna_barplots/S15A_guide_UMAP_NoLegend.pdf", width = 4, height = 2.8)
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

pdf("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/figures_for_manuscript/generate_fig_7b_s15a_s15b_umaps_and_grna_barplots/S15A_guide_UMAP_WithLegend.pdf", width = 4, height = 2.8)
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

###########################
#CREATE CELL TYPE COMPOSITION PLOT
###########################
#read in metadata
metadata <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/figures_for_manuscript/generate_fig_7b_s15a_s15b_umaps_and_grna_barplots/metadata.RDS")
metadata$New_CellType_L2 <- gsub("_", "-", metadata$CellType_L2)
new_cellType_L2_mapper <- c(
  "ExN-Immature" = "ExN-Imm", 
  "ExN-SP" = "ExN-Early", 
  "ExN-Mature" = "ExN-Mat"
)
metadata$New_CellType_L2 <- sapply(metadata$New_CellType_L2, function(x) ifelse(x %in% names(new_cellType_L2_mapper), new_cellType_L2_mapper[x], x))

#set up color vector
color_vector <- c("IN" = "#6f3bbb", 
"ExN-Early" = "#1f63b4", 
"ExN-Mat" = "#12a2c4", 
"ExN-Imm" = "#78a641", 
"IPC-nonDiv" = "#ffd750", 
"IPC-Div" = "#ff7f0e", 
"RG-nonDiv" = "#d63a3a", 
"RG-Div" = "#ff8196")

curr_metadata <- metadata 
plot_df <- table(curr_metadata$guide_assignment_NTC_merged, curr_metadata$New_CellType_L2) %>% as.data.frame()
colnames(plot_df) <- c("guide", "type", "count")

plot_df <- plot_df %>%
  group_by(guide) %>%
  mutate(proportion = count / sum(count))

plot_df$guide <- factor(plot_df$guide, levels = c("NTC", setdiff(unique(plot_df$guide), "NTC")))
cell_type_order <- c("RG-Div", "RG-nonDiv",  "IPC-Div", "IPC-nonDiv",
                     "ExN-Imm", "ExN-Mat", "ExN-Early", "IN")
plot_df$type <- factor(plot_df$type, levels = rev(cell_type_order))

pdf("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/figures_for_manuscript/generate_fig_s15c/s15c_AD_cellTypeBarplot.pdf", width = 8, height = 6)
ggplot(plot_df, aes(x = guide, y = proportion, fill = type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_vector) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.title = element_text(face = "bold")  # Make title bold
  ) +
  labs(x = "Guide", y = "Proportion", fill = "Cell Type")
dev.off()



