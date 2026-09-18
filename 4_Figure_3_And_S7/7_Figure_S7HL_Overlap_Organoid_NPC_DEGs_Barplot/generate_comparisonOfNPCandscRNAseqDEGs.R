#################################
#I/O
#################################
.libPaths(c("/n/groups/walsh/indData/Maya/RLibs", .libPaths()))
library(tidyverse)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(data.table)
library(pheatmap)

#create function to perform hypergeometric test 
generate_overlap_p_val_and_or <- function(organoid_deg, other_deg, universe_other_deg)
{
  universe <- intersect(row.names(organoid_deg), universe_other_deg)
  organoid_deg <- organoid_deg[universe,]
  other_deg <- intersect(other_deg, universe)

  list1 <- row.names(organoid_deg %>% filter(p_val_adj < 0.05))
  list2 <- other_deg

	inList1AndList2 <- length(intersect(list1, list2))
	inList1AndNotList2 <- length(setdiff(list1, list2))
	inList2AndNotList1 <- length(setdiff(list2, list1))
	inNeither <- length(setdiff(universe, c(list1, list2)))

	mat <- matrix(c(inList1AndList2, inList1AndNotList2, inList2AndNotList1, inNeither), nrow = 2)
	pval <- fisher.test(mat, alternative = "greater")$p.value
  estimate <- unname(fisher.test(mat, alternative = "greater")$estimate)


	#also return the matrix values so they can be reported as source data for reviewers
	return(c(pval, estimate, fisher.test(mat, alternative = "greater")$conf.int[1], inList1AndList2, inList1AndNotList2, inList2AndNotList1, inNeither))
}

#read in npc data
file_path <- "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/DEG_Results/WithManualAdjustment/"
sgRNAs <- unname(unlist(read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/final_targets_to_carry_forward.csv"))) #In Input_Files_Not_Generated_By_Scripts
all_DESeq_res <- lapply(sgRNAs, function(x) read.table(paste0(file_path, x, "_vs_NTC_withDedup_manual_adjustment_Cutoff_Of_10.txt")))
names(all_DESeq_res) <- sgRNAs

all_psDEGs <- lapply(all_DESeq_res, function(x) row.names(x %>% filter(padj < 0.05)))
names(all_psDEGs) <- paste0("NPC_", names(all_psDEGs))

basemean_df <-  read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/baseMeanInNTCs.csv", row.names = 1) #In Input_Files_Not_Generated_By_Scripts
universe_npc <- row.names(basemean_df %>% filter(baseMean > 10))

#################################
#LT DAY 30 RGS AND NPC DATA 
#################################
final_DEG_res_list <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/figures_for_manuscript/generate_fig_6a_6c_s13b_s13c/full_DEG_res_list.RDS")

#get degs of interest 
names_of_interest <- names(final_DEG_res_list)[which(grepl("day30", names(final_DEG_res_list)))]
names_of_interest <- names_of_interest[which(grepl("RG", names_of_interest))]
sgRNAs <- paste0(c("BAZ2B", "CLASP1", "EHMT1", "NR2F1-AS1", "PPP3CA",  "ST7", "WDFY3"), "_")
names_of_interest <- names_of_interest[grepl(paste0(sgRNAs, collapse = "|"), names_of_interest)] #added 7/2/25
curr_DEG_res_list <- final_DEG_res_list[names_of_interest]

#get relevant lists 
guides <- unname(unlist(read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas/guides.txt", header = FALSE))) #In /n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/baseMeanInNTCs.csv

curr_scDEGs <- sapply(guides, function(x) curr_DEG_res_list[which(grepl(x, names(curr_DEG_res_list)))])
names(curr_scDEGs) <- guides
curr_scDEGs <- curr_scDEGs[which(sapply(curr_scDEGs, function(x) nrow(x %>% filter(p_val_adj < 0.05))) >= 20)]
curr_guides <- names(curr_scDEGs)

curr_bulkDEGs <- sapply(curr_guides, function(x) all_psDEGs[which(grepl(paste0(gsub("_KD", "", x), "$"), names(all_psDEGs)))])
names(curr_bulkDEGs) <- curr_guides

#run enrichment test (list1 = organoid scRNA DEGs, list2 = bulk NPC DEGs)
overlap_res <- lapply(curr_guides, function(x) generate_overlap_p_val_and_or(curr_scDEGs[[x]], curr_bulkDEGs[[x]], universe_npc))
names(overlap_res) <- curr_guides
overlap_res_df <- as.data.frame(overlap_res) %>% t() %>% as.data.frame()
colnames(overlap_res_df) <- c("pval", "OR", "lower.ci", "InList1AndList2", "InList1AndNotList2", "InList2AndNotList1", "InNeither")
overlap_res_df <- overlap_res_df %>%
  rownames_to_column(var = "guide") %>%  
  mutate(padj = p.adjust(pval, method = "BH"))  
overlap_res_df$guide <- sapply(overlap_res_df$guide , function(x) gsub(".", "-", x, fixed = TRUE))
overlap_res_df$sig <- sapply(overlap_res_df$padj, function(x) ifelse(x <0.0001, "< 0.0001", ifelse(x < 0.001, "< 0.001", ifelse(x < 0.01, "< 0.01", ifelse(x < 0.05, "< 0.05", "Not Significant")))))
overlap_res_df$sig <- factor(as.character(overlap_res_df$sig), 
                              levels = c("< 0.0001", "< 0.001", "< 0.01", "< 0.05", "Not Significant"))
overlap_res_df <- overlap_res_df %>% arrange(desc(OR))
overlap_res_df$guide <- factor(overlap_res_df$guide, levels = overlap_res_df$guide)

#write out the matrix values as source data for the reviewer comment
write.table(overlap_res_df, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Actually_Submitted/For_Third_Submission_And_What_We_Published/4_CLEANED_Figure_3_And_S7/7_Figure_S7HL_Overlap_Organoid_NPC_DEGs_Barplot/6B_LT_day30_RG_and_NPCs_degOverlapBarplot_source_data.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#create barplot
cols <- c(brewer.pal(8, "Greens")[c(8, 6, 4, 2)], "lightgray")
names(cols) <- c("< 0.0001", "< 0.001", "< 0.01", "< 0.05", "Not Significant")

#clean up names 
overlap_res_df$guide_cleaned <- gsub("_KD", "", overlap_res_df$guide)
overlap_res_df$guide_cleaned <- factor(overlap_res_df$guide_cleaned, levels = unique(as.character(overlap_res_df$guide_cleaned)))

pdf("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Actually_Submitted/For_Third_Submission_And_What_We_Published/4_CLEANED_Figure_3_And_S7/7_Figure_S7HL_Overlap_Organoid_NPC_DEGs_Barplot/6B_LT_day30_RG_and_NPCs_degOverlapBarplot.pdf", height = 3, width = 5)
print(ggplot(data = overlap_res_df, aes(x = OR, y = forcats::fct_rev(guide_cleaned), fill = sig, width = 0.8, color = "black")) + 
        geom_bar(stat = "identity", color = "black") +  geom_errorbar(aes(xmin = lower.ci, xmax = OR),
                position = position_dodge(width = 0.8),  # <-- same dodge here
                width = 0.2,
                color = "black") +
        theme_minimal() +
        ylab("Guide") + 
        xlab("Odds Ratio") + 
        theme(axis.text = element_text(size = 7), 
              axis.title = element_text(size = 8, face = "bold")) + 
        scale_fill_manual(values = cols) + 
        guides(fill = guide_legend(title = "Adjusted p-value")) + 
        theme(legend.text = element_text(size = 7), 
              legend.title = element_text(size = 8, face = "bold")) + 
        guides(color = "none") + ggtitle("LT Day 30 RGs and NPCs"))
dev.off()

#################################
#AD DAY 60 AND NPC DATA 
#################################
#read in deg results
path <- paste0("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas_acuteKD/Output/DEG_Results/")

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
nm <- names(flattened_DEG_res_list)
flattened_DEG_res_list <- lapply(seq_along(flattened_DEG_res_list), function(x) {
  gene_name <- strsplit(names(flattened_DEG_res_list)[x], "_")[[1]][1]
  if (gene_name %in% rownames(flattened_DEG_res_list[[x]])) {
    flattened_DEG_res_list[[x]] <- flattened_DEG_res_list[[x]][-which(rownames(flattened_DEG_res_list[[x]]) == gene_name), ]
  }
  return(flattened_DEG_res_list[[x]])
})
names(flattened_DEG_res_list) <- nm
final_DEG_res_list <- flattened_DEG_res_list
sgRNAs <- paste0(c("BAZ2B", "CLASP1", "EHMT1", "NR2F1-AS1", "PPP3CA",  "ST7", "WDFY3"), "_")
final_DEG_res_list <- final_DEG_res_list[grepl(paste0(sgRNAs, collapse = "|"), names(final_DEG_res_list))] #added 7/2/25
print(names(final_DEG_res_list))

types <- c("ExN", "RG", "IPC")
for (type in types)
{
  print(type)

  #get degs of interest 
  names_of_interest <- names(final_DEG_res_list)[which(grepl(type, names(final_DEG_res_list)))]
  curr_DEG_res_list <- final_DEG_res_list[names_of_interest]

  #get relevant lists 
  guides <- unname(unlist(read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/figures/build_organoid_atlas/guides.txt", header = FALSE)))

  curr_scDEGs <- sapply(guides, function(x) curr_DEG_res_list[which(grepl(x, names(curr_DEG_res_list)))])
  names(curr_scDEGs) <- guides
  curr_scDEGs <- curr_scDEGs[which(sapply(curr_scDEGs, function(x) nrow(x %>% filter(p_val_adj < 0.05))) >= 20)]
  curr_guides <- names(curr_scDEGs)

  curr_bulkDEGs <- sapply(curr_guides, function(x) all_psDEGs[which(grepl(paste0(gsub("_KD", "", x), "$"), names(all_psDEGs)))])
  names(curr_bulkDEGs) <- curr_guides

  #run enrichment test (list1 = organoid scRNA DEGs, list2 = bulk NPC DEGs)
  overlap_res <- lapply(curr_guides, function(x) generate_overlap_p_val_and_or(curr_scDEGs[[x]], curr_bulkDEGs[[x]], universe_npc))
  names(overlap_res) <- curr_guides
  overlap_res_df <- as.data.frame(overlap_res) %>% t() %>% as.data.frame()
  colnames(overlap_res_df) <- c("pval", "OR", "lower.ci", "InList1AndList2", "InList1AndNotList2", "InList2AndNotList1", "InNeither")
  overlap_res_df <- overlap_res_df %>%
    rownames_to_column(var = "guide") %>%  
    mutate(padj = p.adjust(pval, method = "BH"))  
  overlap_res_df$guide <- sapply(overlap_res_df$guide , function(x) gsub(".", "-", x, fixed = TRUE))
  overlap_res_df$sig <- sapply(overlap_res_df$padj, function(x) ifelse(x <0.0001, "< 0.0001", ifelse(x < 0.001, "< 0.001", ifelse(x < 0.01, "< 0.01", ifelse(x < 0.05, "< 0.05", "Not Significant")))))
  overlap_res_df$sig <- factor(as.character(overlap_res_df$sig), 
                                levels = c("< 0.0001", "< 0.001", "< 0.01", "< 0.05", "Not Significant"))
  overlap_res_df <- overlap_res_df %>% arrange(desc(OR))
  overlap_res_df$guide <- factor(overlap_res_df$guide, levels = overlap_res_df$guide)

  #write out the matrix values as source data for the reviewer comment
  write.table(overlap_res_df, paste0("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Actually_Submitted/For_Third_Submission_And_What_We_Published/4_CLEANED_Figure_3_And_S7/7_Figure_S7HL_Overlap_Organoid_NPC_DEGs_Barplot/", type, "_and_NPCs_degOverlapBarplot_source_data.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

  #create barplot
  cols <- c(brewer.pal(8, "Greens")[c(8, 6, 4, 2)], "lightgray")
  names(cols) <- c("< 0.0001", "< 0.001", "< 0.01", "< 0.05", "Not Significant")

  pdf(paste0("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Actually_Submitted/For_Third_Submission_And_What_We_Published/4_CLEANED_Figure_3_And_S7/7_Figure_S7HL_Overlap_Organoid_NPC_DEGs_Barplot/", type, "_and_NPCs_degOverlapBarplot.pdf"), height = 3, width = 5)
  print(ggplot(data = overlap_res_df, aes(x = OR, y = forcats::fct_rev(guide), fill = sig, width = 0.8, color = "black")) + 
          geom_bar(stat = "identity", color = "black") + geom_errorbar(aes(xmin = lower.ci, xmax = OR),
                position = position_dodge(width = 0.8),  # <-- same dodge here
                width = 0.2,
                color = "black") +
          theme_minimal() +
          ylab("Guide") + 
          xlab("Odds Ratio") + 
          theme(axis.text = element_text(size = 7), 
                axis.title = element_text(size = 8, face = "bold")) + 
          scale_fill_manual(values = cols) + 
          guides(fill = guide_legend(title = "Adjusted p-value")) + 
          theme(legend.text = element_text(size = 7), 
                legend.title = element_text(size = 8, face = "bold")) + 
          guides(color = "none") + ggtitle(paste0("AD Day 60 ", type, "s and NPCs")))
  dev.off()
}