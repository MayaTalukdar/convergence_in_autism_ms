##########################
#I/O
##########################
library(viper)
library(data.table)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(UpSetR)
library(grid)
library(igraph)
library(fgsea)
library(lsa)
library(readxl)
library(RColorBrewer)
library(readxl)
library(Orthology.eg.db)
library(org.Mm.eg.db)
library(org.Hs.eg.db)

#read in SFARI genes 
sfari <- read.csv("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/sfari.csv", header = TRUE)
SFARI_genes <- unique(sfari$gene.symbol)
high_conf_SFARI_genes <- unique((sfari %>% filter(gene.score == 1))$gene.symbol)

##########################
#DOWNSAMPLE NETWORK
##########################
#create network
network <- read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/full_aracne_network_analysis/Output/All_TFs_And_Epigenetic_Regulators/aracne_output/network_formatted.txt", header = FALSE)
colnames(network) <- c('Regulator', 'Regulatee', 'weight') #note: weight is MI
network <- network %>% filter(Regulatee %in% unique(c(unique(network$Regulator), SFARI_genes, high_conf_SFARI_genes)))
network_og <- network
all_regulators <- unique(network$Regulator)

#downsample it
set.seed(42)
total_edges <- nrow(network)
pdf("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Figure_5_And_S5/Figure_5B_GRN_Network_Vis/network_vis.pdf")
for (prop_downsample in c(1))
{
    print(prop_downsample)

    network <- network_og
    num_to_keep <- round(prop_downsample * total_edges)
    network$probability <- network$weight / sum(network$weight)
    sampled_edges <- network[sample(1:nrow(network), size = num_to_keep, prob = network$probability, replace = FALSE), ]
    network <- sampled_edges[, c("Regulator", "Regulatee", "weight")]

    ##########################
    #CREATE GRAPH FROM NETWORK
    ##########################
    graph <- graph.data.frame(network, directed = TRUE)

    ##########################
    #VISUALIZE NETWORK
    ##########################
    # Set up a force-directed layout
    layout <- layout_with_fr(graph)

    # Add information on node coloring 
    ##in high conf sfari - #540d6e
    ##in all_casd_genes - #ee4266
    ##not in either - #ffb400
    node_colors <- ifelse(
    V(graph)$name %in% high_conf_SFARI_genes, "#fddb27ff", # High confidence SFARI genes
    ifelse(
        V(graph)$name %in% SFARI_genes, "#ee4266",      # All consensus ASD genes
        "lightgray"                                         # Not in either
    )
    )

    # Add information on size 
    vertex_sizes <- ifelse(
    V(graph)$name %in% all_regulators, 1.5,  # Larger size for regulators
    0.75                                     # Smaller size for regulatees
    )

    # Plot the graph using the force-directed layout
    plot(graph, 
    layout = layout,
    vertex.frame.color="black",
    vertex.color = node_colors,          
    vertex.frame.width = 0.2,
    vertex.size = vertex_sizes,             
    vertex.label = NA,  
    vertex.label.cex = NA,      
    edge.arrow.size = 0,
    edge.width = 0.1,  # This controls edge thickness, set to a small value
    edge.color = "black", 
    main = prop_downsample)    
}
dev.off()


    


