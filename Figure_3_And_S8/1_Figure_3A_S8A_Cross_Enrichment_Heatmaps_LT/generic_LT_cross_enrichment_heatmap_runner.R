#!/usr/bin/env Rscript
cellType = commandArgs(trailingOnly=TRUE)[1]

##############################
#I/O
###############################
library(tidyverse)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(data.table)
library(pheatmap)

#create function to perform hypergeometric test 
generate_overlap_p_val <- function(deg_1, deg_2, upreg = FALSE, downreg = FALSE)
{
    universe <- intersect(row.names(deg_1), row.names(deg_2))
    deg_1 <- deg_1[universe,]
    deg_2 <- deg_2[universe,]
    if (upreg)
    {
      list1 <- row.names(deg_1 %>% filter(p_val_adj < 0.05) %>% filter(avg_log2FC > 0))
      list2 <- row.names(deg_2 %>% filter(p_val_adj < 0.05) %>% filter(avg_log2FC > 0))
    } else if (downreg)
    {
      list1 <- row.names(deg_1 %>% filter(p_val_adj < 0.05) %>% filter(avg_log2FC < 0))
      list2 <- row.names(deg_2 %>% filter(p_val_adj < 0.05) %>% filter(avg_log2FC < 0))
    } else
    {
      list1 <- row.names(deg_1 %>% filter(p_val_adj < 0.05))
      list2 <- row.names(deg_2 %>% filter(p_val_adj < 0.05))
    }
    

	inList1AndList2 <- length(intersect(list1, list2))
	inList1AndNotList2 <- length(setdiff(list1, list2))
	inList2AndNotList1 <- length(setdiff(list2, list1))
	inNeither <- length(setdiff(universe, c(list1, list2)))

	mat <- matrix(c(inList1AndList2, inList1AndNotList2, inList2AndNotList1, inNeither), nrow = 2)
	pval <- fisher.test(mat, alternative = "greater")$p.value

	return(pval)
}

#read in deg res 
# full_DEG_res_list <- list()
# for (level in c("L1"))
# {
#     for (day in c("day30", "day60"))
#     {
#         #read in deg results
#         path <- paste0("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas/Output/DEG_Results/", day, "_", level)

#         #read in degs 
#         DEG_res_list <- setNames(
#         lapply(
#         list.dirs(path, recursive = FALSE), 
#         function(x) {
#             setNames(
#             lapply(
#                 list.files(x, pattern = "DEG_res"), 
#                 function(y) readRDS(paste0(x, "/", y))), 
#             gsub("_DEG_res.RDS", "", list.files(x, pattern = "DEG_res")))}), 
#         gsub(
#         paste0(path, "/"), 
#         "", 
#         list.dirs(path, recursive = FALSE)))

#         #flatten nested structure 
#         flattened_DEG_res_list <- unlist(lapply(names(DEG_res_list), function(outer_name) {
#         inner_list <- DEG_res_list[[outer_name]]
#         setNames(inner_list, paste0(outer_name, "_", names(inner_list)))
#         }), recursive = FALSE)
#         names(flattened_DEG_res_list) <- sapply(names(flattened_DEG_res_list), function(x) paste0(day, "_", level, "_", x))
#         nm <- names(flattened_DEG_res_list)
#         flattened_DEG_res_list <- lapply(seq_along(flattened_DEG_res_list), function(x) {
#         gene_name <- strsplit(names(flattened_DEG_res_list)[x], "_")[[1]][1]
#         if (gene_name %in% rownames(flattened_DEG_res_list[[x]])) {
#             flattened_DEG_res_list[[x]] <- flattened_DEG_res_list[[x]][-which(rownames(flattened_DEG_res_list[[x]]) == gene_name), ]
#         }
#         return(flattened_DEG_res_list[[x]])
#         })
#         names(flattened_DEG_res_list) <- nm
#         full_DEG_res_list <- c(full_DEG_res_list, flattened_DEG_res_list)
#     }
# }

# saveRDS(full_DEG_res_list, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/figures_for_manuscript/generate_fig_6a_6c_s13b_s13c/full_DEG_res_list.RDS")
full_DEG_res_list <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/figures_for_manuscript/generate_fig_6a_6c_s13b_s13c/full_DEG_res_list.RDS")
full_DEG_res_list <- full_DEG_res_list[grepl(cellType, names(full_DEG_res_list))]
sgRNAs <- paste0(c("BAZ2B", "CLASP1", "EHMT1", "NR2F1-AS1", "PPP3CA",  "ST7", "WDFY3"), "_")
full_DEG_res_list <- full_DEG_res_list[grepl(paste0(sgRNAs, collapse = "|"), names(full_DEG_res_list))] #added 7/2/25
print(names(full_DEG_res_list))

###############################
#PERFORM CROSS ENRICHMENT 
###############################
cross_enrichment_res_adj_list <- list()

#all
cross_enrichment_res <- sapply(full_DEG_res_list, function(x) sapply(full_DEG_res_list, function(y) generate_overlap_p_val(x, y)))
cross_enrichment_res_adj <- apply(cross_enrichment_res, 2, function(x) p.adjust(x, method = "BH"))
cross_enrichment_res_adj[lower.tri(cross_enrichment_res_adj)] = t(cross_enrichment_res_adj)[lower.tri(cross_enrichment_res_adj)]
cross_enrichment_res_adj_list[["All"]] <- cross_enrichment_res_adj

###############################
#CREATE PLOTS 
###############################
guide_colors <- c(
  "BAZ2B" = "#e75d62",
  "CLASP1" = "#e5a506",
  "EHMT1" = "#dedf61",
  "NR2F1AS1" = "#217885",
  "PPP3CA" = "#4676b7",
  "ST7" = "#273e6b",
  "WDFY3" = "#93418d",
)
names(guide_colors) <- sapply(names(guide_colors), function(x) gsub("AS2", "-AS2", gsub("AS1", "-AS1", x)))

day_colors <- c("day60" = "#F8766D", 
"day30" = "#00BFC4")

#run plotting script
output_path <- "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Figure_3/"

for (type in names(cross_enrichment_res_adj_list))
{
  cross_enrichment_res_adj <- cross_enrichment_res_adj_list[[type]]

  #generate plot
  split_names <- strsplit(rownames(cross_enrichment_res_adj), "_")
  day_names <- sapply(split_names, `[`, 1)  # First part: Day
  guide_names <- sapply(split_names, `[`, 3)  # Third part: Guide
  cell_type_names <- sapply(split_names, `[`, 5)  # Fifth part: Cell Type
  print(cell_type_names)
  annotation_df <- data.frame(
    Day = day_names,
    Guide = guide_names,
    row.names = rownames(cross_enrichment_res_adj)
  )
  annotation_colors <- list(
    Day = day_colors,
    Guide = guide_colors
  )
  breaks <- c(0, 0.0001, 0.001, 0.01, 0.05, 1)
  colors <- c(brewer.pal(8, "Greens")[c(8, 6, 4, 2)], "lightgray")

  pheatmap(
    cross_enrichment_res_adj,
    na_col = "white",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = colors,
    breaks = breaks,
    number_color = "gray",
    border_color = "black",
    legend = FALSE,
    fontsize_row = 0.1,
    fontsize_col = 0.1,
    fontsize = 8,
    annotation_row = annotation_df,  # Apply row annotations
    annotation_col = annotation_df,  # Apply column annotations
    annotation_colors = annotation_colors,  # Apply color scales
    filename = paste0(output_path, "cross_perturbation_enrichment_heatmap_", cellType, ".pdf"),
    width = 10, height = 6
  )
}
