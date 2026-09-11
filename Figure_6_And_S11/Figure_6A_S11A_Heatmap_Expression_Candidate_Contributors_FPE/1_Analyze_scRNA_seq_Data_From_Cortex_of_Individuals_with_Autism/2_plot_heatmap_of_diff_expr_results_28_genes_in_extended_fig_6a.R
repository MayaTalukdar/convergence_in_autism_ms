###################
#I/O
###################
.libPaths(c("/n/groups/walsh/indData/Maya/RLibs", .libPaths()))
library(pheatmap)
library(Seurat)
library(tidyverse)
library(RColorBrewer)

#read in comparisons
comparisons <- unname(unlist(read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/analysis_of_mo_dias_15q_data/comparisons.txt", header = F)))
genes_to_analyze <- c(
  "OPHN1", "SKIL", "HNRNPU", "RPL15", "PSIP1", "KPNB1",
  "SRSF11", "PPP1CB", "ST13", "PDIA3", "CENPE", "TOP2B",
  "PHF3", "ZFX", "PAXBP1", "ZNF711", "LRP2", "SMC1A",
  "USP9X", "DDX3X", "TM9SF3", "SUCO", "IPO7", "TCERG1",
  "SENP2", "YTHDC2", "ARHGEF12", "TULP"
)

#first, read in plot df used to create fig s6a so we can appropriately set our color legend
plot_df_curr_fig_s6a <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Revisions/Redo_Figure_6_And_S6_To_Add_Decasien_Data/Figure_6A_S6A_Heatmap_Expression_Candidate_Contributors_FPE/fig_s6a_plot_df_for_heatmap.RDS")

###################
#CREATE HEATMAPS PER CELL TYPE FOR SUPPLEMENT 
###################
#read in data and format it 
path <- "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/analysis_of_mo_dias_15q_data/Output/DEG_Results/MAST/Case_vs_Control_Male/"
case_vs_control_degs <- setNames(
  lapply(comparisons, function(comp) {
    files <- list.files(paste0(path, comp), pattern = ".RDS")
    setNames(
      lapply(files, function(name) {
        readRDS(paste0(path, comp, "/", name))
      }),
      sapply(files, function(name) strsplit(name, "_")[[1]][1])
    )
  }),
  comparisons
)
case_vs_control_degs <- lapply(case_vs_control_degs, function(x) x[c("Neurons", "IN")])

#format it into a compatible version for pheatmap
make_df_for_heatmap <- function(condition_list, genes_to_analyze, mode = c("log2fc", "sig"), sig_cutoff = 0.05) 
{
  mode <- match.arg(mode)
  cell_types <- names(condition_list)
  
  mat <- lapply(cell_types, function(ct) 
  {
    df <- condition_list[[ct]]

    print(paste(
      "These genes are not expressed:",
      paste(setdiff(genes_to_analyze, row.names(df)), collapse = ", ")
    ))    

    updated_genes_to_analyze <- intersect(genes_to_analyze, row.names(df))

    if (mode == "log2fc") 
    {
      vals <- df[updated_genes_to_analyze, "avg_log2FC", drop = FALSE]
      
    } else if (mode == "sig") 
    {
      sig <- df[updated_genes_to_analyze,  "p_val_adj", drop = FALSE] < sig_cutoff
      vals <- ifelse(sig, "*", "")
      vals <- data.frame(vals,
                          row.names = updated_genes_to_analyze,
                          stringsAsFactors = FALSE)
    }
    
    rownames(vals) <- updated_genes_to_analyze
    print(setdiff(genes_to_analyze, row.names(vals)))
    vals
  })
  
  mat <- as.data.frame(t(as.data.frame(do.call(cbind, mat))))
  rownames(mat) <- cell_types
  
  return(mat)

  print('done')
}
log2fc_by_condition <- lapply(case_vs_control_degs, make_df_for_heatmap,
                              genes_to_analyze = genes_to_analyze, mode = "log2fc")
sigmat_by_condition <- lapply(case_vs_control_degs, make_df_for_heatmap,
                              genes_to_analyze = genes_to_analyze, mode = "sig")

#remove any genes not expressed in our data
remove_all_na_cols <- function(df) 
{
  df[, colSums(!is.na(df)) > 0, drop = FALSE]
}
log2fc_by_condition <- lapply(log2fc_by_condition, remove_all_na_cols)
sigmat_by_condition <- lapply(sigmat_by_condition, remove_all_na_cols)
lapply(log2fc_by_condition, function(x) print(ncol(x)))

#impose cutoff of for significance
sigmat_by_condition <- lapply(seq_along(sigmat_by_condition), function(i) 
{
  mat <- sigmat_by_condition[[i]]
  log2fc <- log2fc_by_condition[[i]]
  
  # Ensure matrices are aligned
  mat <- mat[rownames(log2fc), , drop = FALSE]
  
  # Replace values where abs(log2FC) < 0.1
  mat[abs(log2fc) < 0.1] <- ""
  
  return(mat)
})
names(sigmat_by_condition) <- names(log2fc_by_condition)

#plot heatmaps - note,confirmed that global max and min occur in original fig s6a
paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
global_min <- min(plot_df_curr_fig_s6a, na.rm = TRUE)
#[1] -1.001275
global_max <- max(plot_df_curr_fig_s6a, na.rm = TRUE)
#[1] 0.4794934
plot_limit <- max(abs(global_min), abs(global_max))
myBreaks <- c(seq(-plot_limit, 0, length.out = ceiling(paletteLength / 2) + 1),
              seq(plot_limit / paletteLength, plot_limit, length.out = floor(paletteLength / 2)))


for (cond in names(log2fc_by_condition)) 
{
  pheatmap(
    log2fc_by_condition[[cond]],
    border_color = "black",
    legend = TRUE,
    na_col = "lightgray",
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    color = myColor,
    breaks = myBreaks,
    fontsize_col = 5,
    fontsize_row = 5,
    display_numbers = sigmat_by_condition[[cond]],
    width = 4,
    height = 2,
    main = cond,
    filename = paste0(
      "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Revisions/Analysis_of_15q_and_Idiopathic_Autism_Data/Output/Plots/",
      cond,
      "_case_vs_control_heatmap_28_genes_in_extended_fig_6a.pdf"
    )
  )
}
