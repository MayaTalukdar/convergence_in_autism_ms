###########################
#I/O
###########################
library(tidyverse)
library(data.table)
library(RColorBrewer)
library(ggplot2)
library(ggsci)
library(Seurat)
library(fgsea)
library(pheatmap)
library(cowplot)
library(patchwork)
library(scCustomize)
library(scProportionTest)
library(lme4)
library(sccomp)
library(forcats)
library(tidyr)
library(speckle)

#read in metadata
lt_metadata <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas/Output/Seurat_Objects/9_processed_seurat_obj_fullObj_integrated_harmony_ctrlOnlyAnnotations_qcFiltered_finalAnnots_metadata.RDS")
lt_metadata$New_CellType_L2 <- gsub("_", "-", lt_metadata$CellType_L2)
new_cellType_L2_mapper <- c(
  "ExN-Immature" = "ExN-Imm", 
  "ExN-SP" = "ExN-Early", 
  "ExN-Mature" = "ExN-Mat"
)
lt_metadata$New_CellType_L2 <- sapply(lt_metadata$New_CellType_L2, function(x) ifelse(x %in% names(new_cellType_L2_mapper), new_cellType_L2_mapper[x], x))

#set up color vector
color_vector <- c("IN" = "#6f3bbb", 
"ExN-Early" = "#1f63b4", 
"ExN-Mat" = "#12a2c4", 
"ExN-Imm" = "#78a641", 
"IPC-nonDiv" = "#ffd750", 
"IPC-Div" = "#ff7f0e", 
"RG-nonDiv" = "#d63a3a", 
"RG-Div" = "#ff8196")

sgRNAs <- c("NTC", sapply(c("BAZ2B", "CLASP1", "EHMT1", "NR2F1-AS1", "PPP3CA",  "ST7", "WDFY3"), function(x) paste0(x, '_KD')))
lt_metadata <- lt_metadata %>% filter(guide_assignment_NTC_merged %in% sgRNAs)

###########################
#CREATE LT DAY30 PLOT
###########################
curr_metadata <- lt_metadata %>% filter(day == "day30")
plot_df <- table(curr_metadata$guide_assignment_NTC_merged, curr_metadata$New_CellType_L2) %>% as.data.frame()
colnames(plot_df) <- c("guide", "type", "count")

plot_df <- plot_df %>%
  group_by(guide) %>%
  mutate(proportion = count / sum(count))

plot_df$guide <- factor(plot_df$guide, levels = c("NTC", setdiff(unique(plot_df$guide), "NTC")))
cell_type_order <- c("RG-Div", "RG-nonDiv",  "IPC-Div", "IPC-nonDiv",
                     "ExN-Imm", "ExN-Mat", "ExN-Early", "IN")
plot_df$type <- factor(plot_df$type, levels = rev(cell_type_order))

pdf("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Figure_2/5F_LT_day30_cellTypeBarplot.pdf", width = 8, height = 6)
ggplot(plot_df, aes(x = guide, y = proportion, fill = type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_vector) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.title = element_text(face = "bold")  # Make title bold
  ) +
  labs(x = "Guide", y = "Proportion", fill = "Cell Type", title = "LT Day 30")
dev.off()


###########################
#CREATE LT DAY60 PLOT
###########################
curr_metadata <- lt_metadata %>% filter(day == "day60")
plot_df <- table(curr_metadata$guide_assignment_NTC_merged, curr_metadata$New_CellType_L2) %>% as.data.frame()
colnames(plot_df) <- c("guide", "type", "count")

plot_df <- plot_df %>%
  group_by(guide) %>%
  mutate(proportion = count / sum(count))

plot_df$guide <- factor(plot_df$guide, levels = c("NTC", setdiff(unique(plot_df$guide), "NTC")))
cell_type_order <- c("RG-Div", "RG-nonDiv",  "IPC-Div", "IPC-nonDiv",
                     "ExN-Imm", "ExN-Mat", "ExN-Early", "IN")
plot_df$type <- factor(plot_df$type, levels = rev(cell_type_order))


pdf("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Figure_2/5F_LT_day60_cellTypeBarplot.pdf", width = 8, height = 6)
ggplot(plot_df, aes(x = guide, y = proportion, fill = type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_vector) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.title = element_text(face = "bold")  # Make title bold
  ) +
  labs(x = "Guide", y = "Proportion", fill = "Cell Type", title = "LT Day 60")
dev.off()


