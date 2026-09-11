####################
#I/O
####################
library(Seurat)
library(tidyverse)
library(rlang)
library(scCustomize)
library(patchwork)
library(cowplot)

options(Seurat.object.assay.version = "v3")

merged.Seurat.obj <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas_acuteKD/Output/Seurat_Objects/5_processed_seurat_obj_fullObj_integrated_harmony_organoidAtlasAnnotations_qcFiltered_finalAnnots.RDS")

sgRNAs <- c("NTC", "BAZ2B", "CLASP1", "EHMT1", "NR2F1-AS1", "PPP3CA",  "ST7", "WDFY3")
merged.Seurat.obj$cleaned_guide_assignment <- gsub("_KD", "", merged.Seurat.obj$guide_assignment_NTC_merged)
Idents(merged.Seurat.obj) <- "cleaned_guide_assignment"
merged.Seurat.obj <- subset(merged.Seurat.obj, idents = sgRNAs)

################################
#CREATE BARPLOT OF PROPORTIONS - S13A
################################
metadata <- merged.Seurat.obj@meta.data 
plot_df <- table(metadata$guide_assignment_NTC_merged) %>% as.data.frame()
plot_df$Prop <- plot_df$Freq/sum(plot_df$Freq)
plot_df$Var1 <- sapply(plot_df$Var1, function(x) gsub("_KD", "", x))
genes <- unique(plot_df$Var1)
plot_df <- plot_df %>% mutate(Var1 = factor(Var1, levels = c("NTC", sort(setdiff(unique(genes), "NTC")))))

#set up guide colors
color_vector <- c(
  "NTC" = "#b4b4b4",
  "BAZ2B" = "#e75d62",
  "CLASP1" = "#e5a506",
  "EHMT1" = "#dedf61",
  "NR2F1-AS1" = "#217885",
  "PPP3CA" = "#4676b7",
  "ST7" = "#273e6b",
  "WDFY3" = "#93418d",
)

pdf("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Figure_3/S15B_guide_abundance_AD.pdf", width = 3, height = 4)
ggplot(plot_df, aes(x = 1, y = Prop, fill = forcats::fct_rev(Var1))) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  scale_fill_manual(values = color_vector) + 
  labs(x = "Guide", y = "Proportion", fill = "Guide") +
  theme(
    axis.text.x = element_blank(),  # Remove x-axis labels
    axis.ticks.x = element_blank()  # Remove x-axis ticks
  ) +
  scale_y_continuous(labels = scales::percent_format())
dev.off()