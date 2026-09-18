#!/usr/bin/env Rscript
use_universe = commandArgs(trailingOnly=TRUE)[1] #yes or no 
genes_to_use = commandArgs(trailingOnly=TRUE)[2] #all, SFARI, HCSFARI, CASD, SysID

print("**********************************")
print(paste0("USE UNIVERSE? ", use_universe))
print(paste0("GENES TO USE: ", genes_to_use))
print("**********************************")

###############################
#I/O
###############################
library(tidyverse)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(data.table)
library(pheatmap)
library(ggraph)
library(igraph)
library(scales)
library(grDevices)

cell_type_colors <- c(
    "ExN" = "#388E93",
    "IN" = "#6F3BBB", 
    "IPC" = "#FFAB2F",
    "RG" = "#F83C64"
)

#read in SFARI genes
sfari <- read.csv("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/sfari.csv", header = TRUE) #In Input_Files_Not_Generated_By_Scripts
SFARI_genes <- unique(sfari$gene.symbol)
high_conf_SFARI_genes <- unique((sfari %>% filter(gene.score == 1))$gene.symbol)

#read in deg res 
full_DEG_res_list <- list()
for (level in c("L1"))
{
    for (day in c("day30", "day60"))
    {
        #read in deg results
        path <- paste0("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas/Output/DEG_Results/", day, "_", level)

        #read in degs 
        DEG_res_list <- setNames(
        lapply(
        list.dirs(path, recursive = FALSE), 
        function(x) {
            setNames(
            lapply(
                list.files(x, pattern = "DEG_res"), 
                function(y) readRDS(paste0(x, "/", y))), 
            gsub("_DEG_res.RDS", "", list.files(x, pattern = "DEG_res")))}), 
        gsub(
        paste0(path, "/"), 
        "", 
        list.dirs(path, recursive = FALSE)))

        #flatten nested structure 
        flattened_DEG_res_list <- unlist(lapply(names(DEG_res_list), function(outer_name) {
        inner_list <- DEG_res_list[[outer_name]]
        setNames(inner_list, paste0(outer_name, "_", names(inner_list)))
        }), recursive = FALSE)
        names(flattened_DEG_res_list) <- sapply(names(flattened_DEG_res_list), function(x) paste0(day, "_", level, "_", x))
        nm <- names(flattened_DEG_res_list)
        flattened_DEG_res_list <- lapply(seq_along(flattened_DEG_res_list), function(x) {
        gene_name <- strsplit(names(flattened_DEG_res_list)[x], "_")[[1]][1]
        if (gene_name %in% rownames(flattened_DEG_res_list[[x]])) {
            flattened_DEG_res_list[[x]] <- flattened_DEG_res_list[[x]][-which(rownames(flattened_DEG_res_list[[x]]) == gene_name), ]
        }
        return(flattened_DEG_res_list[[x]])
        })
        names(flattened_DEG_res_list) <- nm
        full_DEG_res_list <- c(full_DEG_res_list, flattened_DEG_res_list)
    }
}

sgRNAs <- paste0(c("BAZ2B", "CLASP1", "EHMT1", "NR2F1-AS1", "PPP3CA",  "ST7", "WDFY3"), "_")
full_DEG_res_list <- full_DEG_res_list[grepl(paste0(sgRNAs, collapse = "|"), names(full_DEG_res_list))] #added 7/2/25
just_DEGs <- lapply(full_DEG_res_list, function(x) row.names(x %>% filter(p_val_adj < 0.05)))
just_DEGs <- just_DEGs[-which(sapply(just_DEGs, length) == 0)]
names(just_DEGs) <- sapply(names(just_DEGs), function(x) gsub("L1_", "", gsub("_KD", "", x)))

#create universe
universe_list <- lapply(full_DEG_res_list, function(x) row.names(x))
names(universe_list) <- names(just_DEGs)
universe <- Reduce(intersect, universe_list)
if (genes_to_use == "SFARI")
{
  universe <- intersect(universe, SFARI_genes)
} else if (genes_to_use == "HCSFARI")
{
  universe <- intersect(universe, high_conf_SFARI_genes)
} 

###############################
#ADJ. MATRIX SETUP
#main idea here is to create an edge between two nodes if a gene is diff expr in multiple perturbations
###############################
#create base adj matrix
##idea here is that we have a row for each comparison (besides self comparisons) and edge weight is the number of shared degs 
adj_matrix <- bind_rows(lapply(names(just_DEGs), function(x) {
  bind_rows(lapply(names(just_DEGs), function(y) {
    if (x != y) {  #no self-loops
      shared_genes <- intersect(just_DEGs[[x]], just_DEGs[[y]])

      if (use_universe == "yes")
      {
        shared_genes <- intersect(shared_genes, universe)
      } 

      if (genes_to_use == "SFARI")
      {
        shared_genes <- intersect(shared_genes, SFARI_genes)
      } else if (genes_to_use == "HCSFARI")
      {
        shared_genes <- intersect(shared_genes, high_conf_SFARI_genes)
      } else if (genes_to_use == "CASD")
      {
        shared_genes <- intersect(shared_genes, casd)
      } else if (genes_to_use == "SysID")
      {
        shared_genes <- intersect(shared_genes, id_genes)
      } 

      if (length(shared_genes) > 0) {
        tibble(x_name = x, 
               y_name = y, 
               edge_weight_all = length(shared_genes))
      }
    }
  }))
}))

###############################
#CREATE GRAPH (without edge width scaling)
#main idea here is to create an edge between two nodes if a gene is diff expr in multiple perturbations
###############################
#add cell type names
adj_matrix <- adj_matrix %>%
  mutate(
    x_type = sapply(x_name, function(n) strsplit(n, "_")[[1]][3]),
    y_type = sapply(y_name, function(n) strsplit(n, "_")[[1]][3]),
    typeOfEdge = ifelse(x_type == y_type, "Intra", "Inter")
  )

#rescale edge width
adj_matrix <- adj_matrix %>%
  mutate(edge_weight_capped = edge_weight_all) %>%
  mutate(edge_weight_scaled = rescale(edge_weight_capped, to = c(2, 6))) %>% as.data.frame()

#sort edges in descending order by edge weight
adj_matrix <- adj_matrix %>%
  arrange(desc(edge_weight_all))

#shuffle rows within each edge_weight_all group, keeping all columns intact
adj_matrix <- adj_matrix %>%
  group_by(edge_weight_all) %>%
  mutate(row_id = row_number()) %>%  # Create a row identifier
  ungroup() %>%
  arrange(edge_weight_all, sample(row_id)) %>%  # Shuffle rows within each group
  select(-row_id) %>% as.data.frame()

#examine whether there are more shared degs between cell types than within cell type
edge_stats_df <- adj_matrix %>% group_by(typeOfEdge) %>% summarize(medianNumSharedDEGs = median(edge_weight_all)) %>% mutate(universe = use_universe) %>% mutate(genes = genes_to_use)
print(edge_stats_df)
ratio <- edge_stats_df %>%
    summarize(ratio = medianNumSharedDEGs[typeOfEdge == "Intra"] / medianNumSharedDEGs[typeOfEdge == "Inter"]) %>% dplyr::select(ratio) %>% as.numeric()

##bootstrap 
n_boot <- 1000
boot_ratios <- c()
for (i in seq_len(n_boot))
{
  sample_adj_matrix <- adj_matrix %>% sample_n(n(), replace = TRUE)
  
  boot_ratio <- sample_adj_matrix %>%
    group_by(typeOfEdge) %>%
    summarize(medianNumSharedDEGs = median(edge_weight_all), .groups = "drop") %>%
    mutate(universe = use_universe, genes = genes_to_use) %>%
    summarize(ratio = medianNumSharedDEGs[typeOfEdge == "Intra"] / medianNumSharedDEGs[typeOfEdge == "Inter"]) %>% dplyr::select(ratio) %>% as.numeric()
    boot_ratios <- c(boot_ratios, boot_ratio)
}
ci_95 <- quantile(boot_ratios, probs = c(0.025, 0.975))
edge_stats_df <- edge_stats_df %>% mutate(intraToInterRatio = ratio, lower.ci = ci_95[1], upper.ci = ci_95[2])

#create graph
g <- graph_from_data_frame(adj_matrix, directed = FALSE)

##set up edge attributes
E(g)$width <- adj_matrix$edge_weight_scaled  
E(g)$color <- ifelse(adj_matrix$typeOfEdge == "Intra", "#0072B5", "#C57633")

##set up vertex attributes
V(g)$cell_type <- sapply(V(g)$name, \(n) strsplit(n, "_")[[1]][3])
V(g)$day <- sapply(V(g)$name, \(n) strsplit(n, "_")[[1]][1])
V(g)$shape <- ifelse(V(g)$day == "day30", "circle", "square")

##plot
set.seed(123)  
layout <- layout_with_fr(g, 
                         niter = 500,  
                         dim = 2,  
                         weights = E(g)$weight)  

output_path <- "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Figure_3/"
pdf(paste0(output_path, "sharedDEGPlot_useUniverse=", use_universe, "_genesToUse=", genes_to_use, ".pdf"), width = 20, height = 20)

plot(g, layout = layout,
     edge.color = E(g)$color,  
     vertex.color = cell_type_colors[V(g)$cell_type],
     vertex.shape = V(g)$shape,
     vertex.size = 5,  # Increased vertex size
     vertex.label = NA,
     vertex.frame.color = "black",  # Dark black outline for vertices
     vertex.frame.width = 7,  # Increased outline width (up to 7)
     main = paste0("use_universe=", use_universe, "; genes_to_use=", genes_to_use))

legend("topright", legend = c("Intra-Cell Type", "Inter-Cell Type"), 
       col = c("#0072B5", "#C57633"), lty = 1, lwd = 2)
legend("bottomright", legend = c("Day 30", "Day 60"), 
       pch = c(21, 22), pt.bg = "white", col = "black", pt.cex = 1.5)
legend("bottomleft", legend = unique(V(g)$cell_type), 
       fill = cell_type_colors[unique(V(g)$cell_type)],  
       title = "Cell Type", cex = 0.8)  

##add edge stats in top left corner
edge_text <- paste(edge_stats_df$typeOfEdge, ":", round(edge_stats_df$medianNumSharedDEGs, 2))
text(x = par("usr")[1] + 0.05 * diff(par("usr")[1:2]),  
     y = par("usr")[4] - 0.05 * diff(par("usr")[3:4]),  
     labels = paste(edge_text, collapse = "\n"),  
     adj = c(0, 1),  
     cex = 1,  
     font = 2,  
     col = "black")  

##add edge width legend based on edge_weight_all
legend_labels <- c(
  paste("Min number of shared DEGs =", min(adj_matrix$edge_weight_all)),
  paste("Max number of shared DEGs =",  max(adj_matrix$edge_weight_all))
)

legend("topleft", legend = legend_labels, 
       lty = 1, lwd = c(0.1, 5), 
       col = "black", title = "Edge Width")

dev.off()




