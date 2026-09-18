##################
#I/O
##################
.libPaths(c("/n/groups/walsh/indData/Maya/RLibs", .libPaths()))
library(tidyverse)
library(readxl)    
library(ggplot2)
library(fgsea)
library(data.table)
library(Orthology.eg.db)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(RColorBrewer)

#add helper functions 
generate_overlap_p_val <- function(list1, list2, universe)
{
  list1 <- unique(toupper(list1))
  list2 <- unique(toupper(list2))
  universe <- unique(toupper(universe))
	list1 <- intersect(list1, universe)
	list2 <- intersect(list2, universe)


	inList1AndList2 <- length(intersect(list1, list2))
	inList1AndNotList2 <- length(setdiff(list1, list2))
	inList2AndNotList1 <- length(setdiff(list2, list1))
	inNeither <- length(setdiff(universe, c(list1, list2)))


	mat <- matrix(c(inList1AndList2, inList1AndNotList2, inList2AndNotList1, inNeither), nrow = 2)
	pval <- fisher.test(mat, alternative = "greater")$p.value


	#also return the 4 matrix values so they can be reported as source data for reviewers
	return(c(pval, fisher.test(mat, alternative = "greater")$estimate, fisher.test(mat, alternative = "greater")$conf.int[1], inList1AndList2, inList1AndNotList2, inList2AndNotList1, inNeither))
}

read_excel_allsheets <- function(filename, tibble = FALSE) {
    # I prefer straight data.frames
    # but if you like tidyverse tibbles (the default with read_excel)
    # then just pass tibble = TRUE
    sheets <- readxl::excel_sheets(filename)
    x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
    if(!tibble) x <- lapply(x, as.data.frame)
    names(x) <- sheets
    x
}

mapfun <- function(mousegenes)
{
    gns <- mapIds(org.Mm.eg.db, mousegenes, "ENTREZID", "SYMBOL")
    mapped <- select(Orthology.eg.db, gns, "Homo_sapiens","Mus_musculus")
    naind <- is.na(mapped$Homo_sapiens)
    hsymb <- mapIds(org.Hs.eg.db, as.character(mapped$Homo_sapiens[!naind]), "SYMBOL", "ENTREZID")
    out <- data.frame(Mouse_symbol = mousegenes, mapped, Human_symbol = NA)
    out$Human_symbol[!naind] <- hsymb
    out
}

#read in ontologies
##read in SFARI genes 
sfari <- read.csv("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/sfari.csv", header = TRUE) #In Input_Files_Not_Generated_By_Scripts
SFARI_genes <- unique(sfari$gene.symbol)
high_conf_SFARI_genes <- unique((sfari %>% filter(gene.score == 1))$gene.symbol)

##read in ID genes 
id_genes <- read.csv("/n/groups/walsh/indData/becky/SysNDD_figures_250823/Input/SysNDD_ID_genes_definitive.csv", header = TRUE) #In Input_Files_Not_Generated_By_Scripts
just_id_genes <- setdiff(id_genes$symbol, SFARI_genes)

##put ontologies in a single list
pathway_list <- list()
pathway_list[["SFARI"]] <- SFARI_genes
pathway_list[["HC-SFARI"]] <- high_conf_SFARI_genes
pathway_list[["sysID"]] <- just_id_genes

#read in relevant universes
##npcs
basemean_df <-  read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/baseMeanInNTCs.csv", row.names = 1) #In Input_Files_Not_Generated_By_Scripts
npc_universe <- row.names(basemean_df %>% filter(baseMean > 10))

##brainspan
brainspan_data <- as.data.frame(fread("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/brainspan/brainspan.csv")) #In Input_Files_Not_Generated_By_Scripts
row.names(brainspan_data) <- brainspan_data$V1
brainspan_data$V1 <- NULL
rownames <- unname(unlist(read.csv("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/brainspan/rows_metadata.csv")$gene_symbol))
brainspan_data <- brainspan_data[-which(duplicated(rownames)),]
rownames <- rownames[-which(duplicated(rownames))]
row.names(brainspan_data) <-rownames
brainspan_data_vec <- apply(brainspan_data, 1, mean) 
brainspan_universe <- names(brainspan_data_vec)[which(brainspan_data_vec > 1)]

#read in rest data
rest_data <- read_excel_allsheets("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/rockowitz_rest_data/nar-00580-h-2015-File009.xlsx")[c(2,3)] #In Input_Files_Not_Generated_By_Scripts
names(rest_data) <- c("human", "mouse")
rest_data_targets <- lapply(rest_data, function(x) x$"gene name")
rest_data_targets <- lapply(rest_data_targets, function(x) x[-which(x == "NA")])

##convert rest mouse targets to human gene names
rest_data_targets[["mouse"]] <- (mapfun(rest_data_targets[["mouse"]]) %>% filter(!is.na(Human_symbol)))$Human_symbol

##################
#PLOTS
##################
#run enrichment test 
df <- sapply(rest_data_targets, function(x) sapply(pathway_list, function(y) generate_overlap_p_val(x, y, brainspan_universe)[1]))
df <- apply(df, 2, function(x) p.adjust(x, method = "BH"))

#create df for plot 
plot_df <- reshape2::melt(df) %>% arrange(Var2, desc(Var1)) %>% mutate(cleaned_name = paste0(Var1, " (", Var2, ")")) 
colnames(plot_df) <- c("geneSet", "species", "p_adj", "cleaned_name")
plot_df <- plot_df %>% mutate(Pconvert = -log10(p_adj))
plot_df$cleaned_name <- factor(as.character(plot_df$cleaned_name), levels = as.character(plot_df$cleaned_name))
plot_df$sig <- sapply(plot_df$p_adj, function(x) ifelse(x <0.0001, "< 0.0001", ifelse(x < 0.001, "< 0.001", ifelse(x < 0.01, "< 0.01", ifelse(x < 0.05, "< 0.05", "Not Significant")))))
plot_df$sig <- factor(as.character(plot_df$sig), levels = c("< 0.0001", "< 0.001", "< 0.01", "< 0.05", "Not Significant"))

#add OR
plot_df$OR <- apply(plot_df, 1, function(x) generate_overlap_p_val(rest_data_targets[[as.character(x["species"])]], pathway_list[[as.character(x["geneSet"])]], brainspan_universe)[2])
plot_df <- plot_df %>% arrange(desc(OR))
plot_df$lower.ci <- apply(plot_df, 1, function(x) generate_overlap_p_val(rest_data_targets[[as.character(x["species"])]], pathway_list[[as.character(x["geneSet"])]], brainspan_universe)[3])

#add the matrix counts too (list1 = rest targets, list2 = gene set) so reviewers can see what actually went into each fisher test
plot_df$inList1AndList2 <- apply(plot_df, 1, function(x) generate_overlap_p_val(rest_data_targets[[as.character(x["species"])]], pathway_list[[as.character(x["geneSet"])]], brainspan_universe)[4])
plot_df$inList1AndNotList2 <- apply(plot_df, 1, function(x) generate_overlap_p_val(rest_data_targets[[as.character(x["species"])]], pathway_list[[as.character(x["geneSet"])]], brainspan_universe)[5])
plot_df$inList2AndNotList1 <- apply(plot_df, 1, function(x) generate_overlap_p_val(rest_data_targets[[as.character(x["species"])]], pathway_list[[as.character(x["geneSet"])]], brainspan_universe)[6])
plot_df$inNeither <- apply(plot_df, 1, function(x) generate_overlap_p_val(rest_data_targets[[as.character(x["species"])]], pathway_list[[as.character(x["geneSet"])]], brainspan_universe)[7])

#write out the matrix values as source data for the reviewer comment
write.table(plot_df, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Actually_Submitted/For_Third_Submission_And_What_We_Published/5_CLEANED_5_Figure_S8/1_Figure_S8C_REST_Target_Enrichment_Barplot/rest_barplot_source_data.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#plot
cols <- c(brewer.pal(8, "Greens")[c(8, 6, 4, 2)], "lightgray")
pdf("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Actually_Submitted/For_Third_Submission_And_What_We_Published/5_CLEANED_5_Figure_S8/1_Figure_S8C_REST_Target_Enrichment_Barplot/rest_barplot.pdf", height = 3, width = 5)
ggplot(data = plot_df, aes(x = OR, y = forcats::fct_rev(species), fill = sig, linetype = geneSet)) +
  geom_bar(stat = "identity", position=position_dodge2(reverse = TRUE), width = 0.8, color = "black") +
     geom_errorbar(aes(xmin = lower.ci, xmax = OR), position=position_dodge2(reverse = TRUE), width = 0.8, color = "black") +
  scale_fill_manual(values = cols) +
  theme_minimal() +
  ylab("Species") +
  xlab("Odds Ratio") +
  theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8, face = "bold")) +
  guides(fill = guide_legend(title = "Adjusted p-value")) +
  theme(legend.text = element_text(size = 7), legend.title = element_text(size = 8, face = "bold"))
dev.off()