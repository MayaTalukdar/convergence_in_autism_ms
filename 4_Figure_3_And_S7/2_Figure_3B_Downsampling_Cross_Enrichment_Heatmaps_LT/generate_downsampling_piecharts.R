###############################
#I/O
###############################
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(Seurat)
library(patchwork)
library(cowplot)

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

#read in deg results
final_DEG_res_list <- list()
level <- "L1"

for (day in c("day60"))
{
  path <- paste0("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas/Output/DEG_Results_Downsampled/", day, "_", level, "/")
  guides <- unname(unlist(read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas/guides_downsample.txt", header = FALSE))) #In Input_Files_Not_Generated_By_Scripts
  numcells <- unname(unlist(read.table(paste0("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas/numcells_", day, "_downsample.txt"), header = FALSE)))
  numcells <- c(numcells, -1)
  all_deg_paths <- lapply(guides, function(x) sapply(numcells, function(y) paste0(path, x, "/DownsamplingTo", y, "Cells/Replicate=1")))
  all_deg_paths <- do.call(c, all_deg_paths)
  DEG_res_list <- lapply(all_deg_paths, function(x) setNames(lapply(list.files(x, pattern= "DEG_res"), function(y) readRDS(paste0(x,"/", y))), gsub("_DEG_res.RDS", "", list.files(x, pattern = "DEG_res"))))
  names(DEG_res_list) <- unname(sapply(all_deg_paths, function(x) gsub("Cells", "", gsub("DownsamplingTo", "", gsub("/", "-", gsub(path, "", gsub("/Replicate=1", "", x)), fixed = TRUE)))))

  #flatten nested structure 
  flattened_DEG_res_list <- unlist(lapply(names(DEG_res_list), function(outer_name) {
  inner_list <- DEG_res_list[[outer_name]]
  setNames(inner_list, paste0(outer_name, "_", names(inner_list)))
  }), recursive = FALSE)
  names(flattened_DEG_res_list) <- gsub("-1_", "All_", names(flattened_DEG_res_list), fixed = TRUE)
  names(flattened_DEG_res_list) <- paste(day, names(flattened_DEG_res_list), sep = "_")
  final_DEG_res_list <- c(final_DEG_res_list, flattened_DEG_res_list)
}

###############################
#PERFORM CROSS ENRICHMENT & PLOTTING - ALL 
###############################
#get relevant files 
n <- "All"
flattened_DEG_res_list <- final_DEG_res_list[which(grepl(paste0("-", n, "_"), names(final_DEG_res_list)))]

#run cross-enrichment
cross_enrichment_res_adj_list <- list()

##all
cross_enrichment_res <- sapply(flattened_DEG_res_list, function(x) sapply(flattened_DEG_res_list, function(y) generate_overlap_p_val(x, y)))
cross_enrichment_res_adj <- apply(cross_enrichment_res, 2, function(x) p.adjust(x, method = "BH"))
cross_enrichment_res_adj[lower.tri(cross_enrichment_res_adj)] = t(cross_enrichment_res_adj)[lower.tri(cross_enrichment_res_adj)]
cross_enrichment_res_adj_list[["All"]] <- cross_enrichment_res_adj

#run plotting 
pie_chart_list <- list()
for (type in names(cross_enrichment_res_adj_list))
{
  #set up matrix
  cross_enrichment_res_adj <- cross_enrichment_res_adj_list[[type]]
  breaks <- c(0, 0.0001, 0.001, 0.01, 0.05)
  colors <- c(brewer.pal(8, "Greens")[c(8, 6, 4, 2)])
  flat_df <- as.data.frame(as.table(as.matrix(cross_enrichment_res_adj)))
  colnames(flat_df) <- c("row", "col", "value")
  flat_df$combined <- paste(flat_df$row, flat_df$col, sep = "_")
  flat_df <- flat_df %>% filter(value < 0.05) 
  pathways_to_keep <- flat_df$combined

  #plot pie chart 
  flat_df$bin <- cut(flat_df$value, breaks = breaks, include.lowest = TRUE)
  bin_summary <- as.data.frame(table(flat_df$bin))
    colnames(bin_summary) <- c("bin", "count")

  # Plot pie chart
  p <- ggplot(bin_summary, aes(x = "", y = count, fill = bin)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar("y") +
  scale_fill_manual(values = colors, drop = FALSE) +
  theme_void() + ggtitle("All") + theme(plot.title = element_text(hjust = 0.5))
  pie_chart_list[["All"]] <-p
}

###############################
#PERFORM CROSS ENRICHMENT FOR DOWNSAMPLES
###############################
numcells <- setdiff(numcells, "-1")

for (n in numcells)
{
  print("=========")
  print(n)
  print("=========")

  flattened_DEG_res_list <- final_DEG_res_list[which(grepl(paste0("-", n, "_"), names(final_DEG_res_list)))]

  cross_enrichment_res_adj_list <- list()

  #all
  cross_enrichment_res <- sapply(flattened_DEG_res_list, function(x) sapply(flattened_DEG_res_list, function(y) generate_overlap_p_val(x, y)))
  cross_enrichment_res_adj <- apply(cross_enrichment_res, 2, function(x) p.adjust(x, method = "BH"))
  cross_enrichment_res_adj[lower.tri(cross_enrichment_res_adj)] = t(cross_enrichment_res_adj)[lower.tri(cross_enrichment_res_adj)]
  cross_enrichment_res_adj_list[["All"]] <- cross_enrichment_res_adj

  for (type in names(cross_enrichment_res_adj_list))
  {
    cross_enrichment_res_adj <- cross_enrichment_res_adj_list[[type]]

    #generate plot
    breaks <- c(0, 0.0001, 0.001, 0.01, 0.05, 1)
    colors <- c(brewer.pal(8, "Greens")[c(8, 6, 4, 2)], "lightgray")
    flat_df <- as.data.frame(as.table(as.matrix(cross_enrichment_res_adj)))
    colnames(flat_df) <- c("row", "col", "value")
    flat_df$combined <- paste(flat_df$row, flat_df$col, sep = "_")
    pathways_to_keep_custom <- gsub("-All", paste0('-', n), pathways_to_keep)
    flat_df <- flat_df %>% filter(combined %in% pathways_to_keep_custom)
    print(dim(flat_df))

    #plot pie chart 
    flat_df$bin <- cut(flat_df$value, breaks = breaks, include.lowest = TRUE)
    bin_summary <- as.data.frame(table(flat_df$bin))
        colnames(bin_summary) <- c("bin", "count")

    # Plot pie chart
    p <- ggplot(bin_summary, aes(x = "", y = count, fill = bin)) +
    geom_bar(stat = "identity", width = 1, color = "black") +
    coord_polar("y") +
    scale_fill_manual(values = colors, drop = FALSE) +
    theme_void() + ggtitle(paste('n=', n)) + theme(plot.title = element_text(hjust = 0.5))
    pie_chart_list[[paste0("n=", n)]] <-p
  }
}

###############################
#CREATE WRAPPED PLOT
###############################
pie_chart_list <- c(pie_chart_list[2:7], pie_chart_list[1])

pdf("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Figure_3/S14_DownsamplePieCharts.pdf", width = 12, height = 5)
wrap_plots(lapply(pie_chart_list, function(x) x + NoLegend()), nrow = 1)
dev.off()

