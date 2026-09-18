######################
#I/O
######################
.libPaths(c("/n/groups/walsh/indData/Maya/RLibs", .libPaths()))
library(tidyverse)
library(RColorBrewer)
library(data.table)
library(DESeq2)

final_targets <-  c(unname(unlist(read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/final_targets_to_carry_forward.csv", header = FALSE))), "ARHGAP5", "ARID2") #In Input_Files_Not_Generated_By_Scripts

######################
#PARSE DDPCR RESULTS
######################
#read in ddpcr data
ddpcr_df <- read.csv("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/calculate_kd_eff_ddpcr/240420_screen_validations_ddPCR.csv", header = TRUE) #In Input_Files_Not_Generated_By_Scripts
targets <- colnames(ddpcr_df)[-(1:2)]
ddpcr_df <- ddpcr_df %>% mutate(type = ifelse(grepl("NTC", Sample), "NTC", "target"))
ddpcr_df$type <- factor(ddpcr_df$type)
ddpcr_df$Batch <- factor(ddpcr_df$Batch)

#calculate kd efficiency 
calc_kd_eff <- function(target)
{
    dat <- ddpcr_df %>% dplyr::select(Batch, type, target)
    colnames(dat)[3] <- "Expr"

    #get per replicate fold changes to plot
    dat_batch <- dat %>%
    group_by(Batch) %>%
    summarise(
        mean_NTC = mean(Expr[type == "NTC"]),
        target = Expr[type == "target"],
        fold_change = target / mean_NTC)

    mod <- lm(Expr ~ Batch + type, data = dat)
    p_val <- summary(mod)$coefficients["typetarget", "Pr(>|t|)"] #return p value
    
    summary_stats <- dat %>%
    group_by(type) %>%
    summarize(mean_expr = mean(Expr))

    #Calculate fold change
    fold_change <- unique((summary_stats %>%
    mutate(fold_change = mean_expr[type == "target"] / mean_expr[type == "NTC"]))$fold_change)

    #Calculate SE 
    se_prop <- sqrt((fold_change * (1- fold_change))/nrow(dat))

    return (c(fold_change, se_prop, p_val, unlist(dat_batch$fold_change)))


}
kd_eff <- sapply(targets, function(x) calc_kd_eff(x)[1:3]) %>% as.data.frame() %>% t() %>% as.data.frame()
colnames(kd_eff) <- c("foldChange", "se", "pvalue")
row.names(kd_eff) <- sapply(row.names(kd_eff), function(x) gsub(".", "-", x, fixed = TRUE))
ddpcr_kd_eff_final <- kd_eff %>% mutate(expLowCI = foldChange - se) %>% mutate(expHighCI = foldChange + se) %>% select(foldChange, expLowCI, expHighCI, pvalue)
ddpcr_kd_eff_final$sgRNA <- row.names(ddpcr_kd_eff_final)
ddpcr_kd_eff_final <- ddpcr_kd_eff_final %>% arrange(foldChange)
ddpcr_kd_eff_final <- ddpcr_kd_eff_final %>% filter(sgRNA %in% final_targets)
ddpcr_kd_eff_final$padj <- p.adjust(ddpcr_kd_eff_final$pvalue, method = "BH")
ddpcr_kd_raw_data <- sapply(intersect(targets, final_targets), function(x) calc_kd_eff(x)[4:6]) %>% as.data.frame() %>% t() %>% as.data.frame()

######################
#PARSE BULK RNA SEQ RESULTS 
######################
path <- "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/DEG_Results/WithManualAdjustment/"
files <- list.files(path)
files <- files[-which(grepl("^NTC", files))]
all_DESeq2_res <- lapply(files, function(x) fread(paste0(path, x), header = FALSE) %>% as.data.frame())
names(all_DESeq2_res) <- sapply(files, function(x) gsub("_vs_NTC_withDedup_manual_adjustment_Cutoff_Of_10.txt", "", x))
all_DESeq2_res_full <- all_DESeq2_res
all_DESeq2_res <- all_DESeq2_res[-which(names(all_DESeq2_res) %in% row.names(kd_eff))]

#get kd efficiency 
kd_efficiency_df <- do.call(rbind, lapply(names(all_DESeq2_res), function(x) all_DESeq2_res[[x]] %>% filter(V1 == x)))
colnames(kd_efficiency_df) <- c("sgRNA", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj_DEseq", "padj")
kd_efficiency_df <- kd_efficiency_df %>% mutate(lowCI = log2FoldChange - lfcSE) %>% mutate(highCI = log2FoldChange + lfcSE) %>% mutate(expLowCI = 2^lowCI) %>% mutate(expHighCI = 2^highCI) %>% mutate(foldChange = 2^log2FoldChange) %>% arrange(foldChange) 
kd_efficiency_df <- kd_efficiency_df %>% dplyr::select(foldChange, expLowCI, expHighCI, pvalue, padj, sgRNA)
kd_efficiency_df_rna <- kd_efficiency_df

######################
#ADD IN RAW REPLICATE DATA FOR BULK RNA SEQ RESULTS
######################
dds_norm <- counts(readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/dds.RDS"), normalized = TRUE) #In Input_Files_Not_Generated_By_Scripts
metadata_merged <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/metadata_merged.RDS") #In Input_Files_Not_Generated_By_Scripts

kd_efficiency_df_raw <- matrix(NA, nrow = length(final_targets), ncol = 3)
rownames(kd_efficiency_df_raw) <- final_targets
colnames(kd_efficiency_df_raw) <- paste0("V", 1:3)
for (i in seq_along(final_targets))
{
  gene <- final_targets[i]
  batches <- unique(metadata_merged$Batch[metadata_merged$sgRNA == gsub("-", "_", gene)])
  for (j in seq_along(batches))
  {
    batch <- batches[j]
    target_samples <- rownames(metadata_merged)[metadata_merged$Batch == batch & metadata_merged$sgRNA ==  gsub("-", "_", gene)]
    ntc_samples <- rownames(metadata_merged)[metadata_merged$Batch == batch & metadata_merged$sgRNA == "NTC"]
    kd_efficiency_df_raw[i, j] <- mean(unlist(dds_norm[gene, target_samples,drop=TRUE])) / mean(unlist(dds_norm[gene, ntc_samples,drop = TRUE]))
  }
}
kd_efficiency_df_raw <- kd_efficiency_df_raw[-which(row.names(kd_efficiency_df_raw) == "SRCAP"),] #since we assessed this gene with ddpcr

######################
#COMBINE DDPCR AND BULK RNA RESULTS
######################
kd_efficiency_df <- rbind(kd_efficiency_df_rna, ddpcr_kd_eff_final) %>% arrange(foldChange)
row.names(kd_efficiency_df) <- NULL
kd_efficiency_df$sgRNA <- factor(as.character(kd_efficiency_df$sgRNA),levels = as.character(kd_efficiency_df$sgRNA))

#also combine raw data
kd_efficiency_df_raw <- rbind(kd_efficiency_df_raw, ddpcr_kd_raw_data)
cleaned_kd_efficiency_df_raw <- kd_efficiency_df_raw
colnames(cleaned_kd_efficiency_df_raw) <- c("Rep1", "Rep2", "Rep3")
write.table(cleaned_kd_efficiency_df_raw, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Actually_Submitted/For_Third_Submission_And_What_We_Published/0_CLEANED_0_Figure_1_And_S1_And_S2_And_S3/1_Figure_S1A_KD_Eff_Barplot/per_replicate_kd_efficiency_summary_table.csv", sep = ",", quote = FALSE, row.names = TRUE, col.names = TRUE)

######################
#MAKE SUMMARY TABLE
######################
num_degs <- sapply(names(all_DESeq2_res_full), function(x) nrow(all_DESeq2_res_full[[x]] %>% filter(V8 < 0.05)))
kd_efficiency_df$num_degs <- num_degs[as.character(kd_efficiency_df$sgRNA)]

kd_efficiency_df$padj_across_all_tested_targets <- p.adjust(kd_efficiency_df$pvalue, method = "BH")
write.table(kd_efficiency_df, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Actually_Submitted/For_Third_Submission_And_What_We_Published/0_CLEANED_0_Figure_1_And_S1_And_S2_And_S3/1_Figure_S1A_KD_Eff_Barplot/kd_efficiency_summary_table.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)

######################
#MAKE PLOT - PADJ
######################
kd_efficiency_df$sig <- sapply(kd_efficiency_df$padj_across_all_tested_targets, function(x) ifelse(is.na(x), "Not Applicable", ifelse(x <0.0001, "< 0.0001", ifelse(x < 0.001, "< 0.001", ifelse(x < 0.01, "< 0.01", ifelse(x < 0.05, "< 0.05", "Not Significant"))))))
kd_efficiency_df$sig <- factor(as.character(kd_efficiency_df$sig), levels = c("< 0.0001", "< 0.001", "< 0.01", "< 0.05", "Not Significant", "Not Applicable"))

#add number of degs
num_degs <- sapply(names(all_DESeq2_res_full), function(x) nrow(all_DESeq2_res_full[[x]] %>% filter(V8 < 0.05)))
kd_efficiency_df$num_degs <- num_degs[as.character(kd_efficiency_df$sgRNA)]

#create plot
cols <- c(brewer.pal(8, "Greens")[c(8, 6, 4, 2)], "lightgray", "white")

kd_efficiency_df_init <- kd_efficiency_df

kd_efficiency_df <- kd_efficiency_df %>% filter(sgRNA %in% final_targets)

kd_efficiency_df_raw_plot <- as.data.frame(kd_efficiency_df_raw) %>%
  tibble::rownames_to_column("sgRNA") %>%
  tidyr::pivot_longer(cols = -sgRNA, names_to = "Rep", values_to = "foldChange") %>%
  filter(sgRNA %in% final_targets)

pdf("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Actually_Submitted/For_Third_Submission_And_What_We_Published/0_CLEANED_0_Figure_1_And_S1_And_S2_And_S3/1_Figure_S1A_KD_Eff_Barplot/knockdown_efficiency_barplot_using_adj_pvalue_just_genes_used_for_downstream_analysis.pdf", height = 3, width = 6)
ggplot(data = kd_efficiency_df, aes(x = sgRNA, y = foldChange, fill = sig, width = 0.8)) + geom_bar(stat = "identity", color = "black") +
theme_minimal() +
xlab("Gene") +
ylab("Relative Expression of Target") +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
geom_errorbar(aes(ymin=expLowCI, ymax=expHighCI), width=.2, position=position_dodge(.9)) +
geom_point(data = kd_efficiency_df_raw_plot, aes(x = sgRNA, y = foldChange), inherit.aes = FALSE, position = position_jitter(width = 0.12, height = 0), size = 1.5, color = "black") +
geom_hline(yintercept = 1, linetype = "dashed") +
theme(axis.text=element_text(size=7), axis.title=element_text(size=8,face="bold")) +
scale_fill_manual(values = cols) +
guides(fill=guide_legend(title="Adj. P-Value")) +
theme(legend.text=element_text(size=7), legend.title=element_text(size=8,face="bold")) +
guides(color = "none")
dev.off()
