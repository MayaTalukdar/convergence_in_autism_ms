###################
#I/O
###################
.libPaths(c("/n/groups/walsh/indData/Maya/RLibs", .libPaths()))
library(tidyverse)
library(dplyr)
library(patchwork)

#read in differentially expressed genes across all perturbations
path <- "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Revisions/Analysis_of_shRNA_Data/0_Process_And_Run_Differential_Expression/differential_expression_analysis/DEG_Results/WithManualAdjustment/"
all_DESeq_res_list <- lapply(list.files(path, pattern = ".txt"), function(x) read.table(paste0(path, x), header = TRUE))
names(all_DESeq_res_list) <- sapply(list.files(path, pattern = ".txt"), function(x) gsub("_vs_ctrl_manual_adjustment_Cutoff_Of_10.txt", "", x))
all_DESeq_res_list <- lapply(all_DESeq_res_list, function(x) x %>% filter(!is.na(padj)))

#read in SFARI genes
sfari <- read.csv("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/sfari.csv", header = TRUE) #In Input_Files_Not_Generated_By_Scripts
SFARI_genes <- unique(sfari$gene.symbol)
all_SFARI_genes <- SFARI_genes

###################
#MAKE VOLCANO PLOTS
###################
make_volcano_plot <- function(current_res, current_title, cut_x_outliers = FALSE)
{
     ##convert to data frame
     current_res <- as.data.frame(current_res)
     current_res$FDR <- current_res$padj
     current_res$gene <- rownames(current_res)

     ##calculate -log10(FDR), significance groups, SFARI status, and KD target
     current_res <- current_res %>%
          mutate(
               neg_log10_FDR = -log10(FDR),
               sig_group = case_when(
                    FDR < 0.05 & log2FoldChange > 0 ~ "Upregulated",
                    FDR < 0.05 & log2FoldChange < 0 ~ "Downregulated",
                    TRUE ~ "NS"
               ),
               is_SFARI = gene %in% all_SFARI_genes,
               is_target = gene == current_title
          )

     ##replace infinite values from FDR = 0
     max_y <- max(current_res$neg_log10_FDR[is.finite(current_res$neg_log10_FDR)], na.rm = TRUE)
     current_res$neg_log10_FDR[!is.finite(current_res$neg_log10_FDR)] <- max_y + 1

     ##set point colors
     current_res$fill_group <- current_res$sig_group
     current_res$fill_group[current_res$is_target] <- "KD target"

     ##set plotting data and remove x-axis outliers
     plot_res <- current_res

     if(cut_x_outliers){
          x_limits <- quantile(plot_res$log2FoldChange, c(0.01, 0.99), na.rm = TRUE)
          plot_res <- plot_res %>% filter((log2FoldChange >= x_limits[1] &
                                           log2FoldChange <= x_limits[2]) |
                                          is_target)
     }

     ##count significantly up- and downregulated genes
     n_up <- sum(current_res$FDR < 0.05 & current_res$log2FoldChange > 0, na.rm = TRUE)
     n_down <- sum(current_res$FDR < 0.05 & current_res$log2FoldChange < 0, na.rm = TRUE)

     ##compute downregulation statistics
     sig_genes <- current_res %>% filter(sig_group != "NS")

     overall_down <- sum(sig_genes$log2FoldChange < 0)
     overall_prop <- ifelse(nrow(sig_genes) == 0, NA, overall_down / nrow(sig_genes))

     sig_non_sfari <- sig_genes %>% filter(!is_SFARI)
     non_sfari_down <- sum(sig_non_sfari$log2FoldChange < 0)
     non_sfari_prop <- ifelse(nrow(sig_non_sfari) == 0, NA, non_sfari_down / nrow(sig_non_sfari))

     sig_sfari <- sig_genes %>% filter(is_SFARI)
     sfari_down <- sum(sig_sfari$log2FoldChange < 0)
     sfari_prop <- ifelse(nrow(sig_sfari) == 0, NA, sfari_down / nrow(sig_sfari))

     ##generate volcano plot
     ggplot(plot_res,
            aes(x = log2FoldChange,
                y = neg_log10_FDR)) +
          geom_point(aes(fill = fill_group),
                     shape = 21,
                     color = "gray50",
                     size = 1.8,
                     stroke = 0.2,
                     alpha = 0.8) +
          geom_point(data = subset(plot_res, is_SFARI & sig_group != "NS"),
                     aes(fill = fill_group),
                     shape = 21,
                     color = "black",
                     size = 1.8,
                     stroke = 0.8,
                     alpha = 0.8) +
          geom_point(data = subset(plot_res, is_target),
                     aes(fill = fill_group),
                     shape = 21,
                     color = "black",
                     size = 2.2,
                     stroke = 0.8) +
          geom_hline(yintercept = -log10(0.05),
                     linetype = "dashed") +
          geom_vline(xintercept = 0,
                     linetype = "dashed") +
          scale_fill_manual(values = c(
               "Upregulated" = "red",
               "Downregulated" = "blue",
               "NS" = "gray75",
               "KD target" = "purple"
          ),
          name = NULL) +
          labs(
               title = paste0(
                    current_title,
                    "\nOverall Down: ", overall_down, "/", nrow(sig_genes),
                    " (", scales::percent(overall_prop, accuracy = 0.1), ")",
                    " | Non-SFARI: ", non_sfari_down, "/", nrow(sig_non_sfari),
                    " (", scales::percent(non_sfari_prop, accuracy = 0.1), ")",
                    " | SFARI: ", sfari_down, "/", nrow(sig_sfari),
                    " (", scales::percent(sfari_prop, accuracy = 0.1), ")"
               ),
               x = "Log2 Fold Change",
               y = expression(-log[10]("Adjusted p-value"))
          ) +
          theme_classic() +
          theme(
               plot.title = element_text(hjust = 0.5,
                                         face = "bold",
                                         size = 11),
               legend.position = "right"
          )
}

#generate plot without x outliers
volcano_plot_list_cut <- Map(
     make_volcano_plot,
     all_DESeq_res_list,
     names(all_DESeq_res_list),
     MoreArgs = list(cut_x_outliers = TRUE)
)

combined_volcano_plot_cut <- wrap_plots(
     volcano_plot_list_cut,
     ncol = 3)

pdf("Combined_Volcano_Plots_shRNA_KDs_cut_x_outliers.pdf",
    width = 30,
    height = 18)
print(combined_volcano_plot_cut)
dev.off()

volcano_plot_list_full <- Map(
     make_volcano_plot,
     all_DESeq_res_list,
     names(all_DESeq_res_list),
     MoreArgs = list(cut_x_outliers = FALSE)
)

#generate plot with x outliers
combined_volcano_plot_full <- wrap_plots(
     volcano_plot_list_full,
     ncol = 3)

pdf("Combined_Volcano_Plots_shRNA_KDs_full_x_axis.pdf",
    width = 24,
    height = 6)
print(combined_volcano_plot_full)
dev.off()