######################
#I/O
######################
.libPaths(c("/n/groups/walsh/indData/Maya/RLibs", .libPaths()))
library(tidyverse)
library(RColorBrewer)
library(data.table)

all_targets <-  unname(unlist(read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/differential_expression_analysis/all_perturbations.csv", header = FALSE)))
final_targets <-  unname(unlist(read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/final_targets_to_carry_forward.csv", header = FALSE)))

######################
#PARSE DDPCR RESULTS
######################
#read in ddpcr data
ddpcr_df <- read.csv("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/calculate_kd_eff_ddpcr/240420_screen_validations_ddPCR.csv", header = TRUE)
targets <- colnames(ddpcr_df)[-(1:2)]
ddpcr_df <- ddpcr_df %>% mutate(type = ifelse(grepl("NTC", Sample), "NTC", "target"))
ddpcr_df$type <- factor(ddpcr_df$type)
ddpcr_df$Batch <- factor(ddpcr_df$Batch)

#calculate kd efficiency 
calc_kd_eff <- function(target)
{
    dat <- ddpcr_df %>% dplyr::select(Batch, type, target)
    colnames(dat)[3] <- "Expr"
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

    return (c(fold_change, se_prop, p_val))

}
kd_eff <- sapply(targets, function(x) calc_kd_eff(x)) %>% as.data.frame() %>% t() %>% as.data.frame()
colnames(kd_eff) <- c("foldChange", "se", "pvalue")
row.names(kd_eff) <- sapply(row.names(kd_eff), function(x) gsub(".", "-", x, fixed = TRUE))
ddpcr_kd_eff_final <- kd_eff %>% mutate(expLowCI = foldChange - se) %>% mutate(expHighCI = foldChange + se) %>% select(foldChange, expLowCI, expHighCI, pvalue)
ddpcr_kd_eff_final$sgRNA <- row.names(ddpcr_kd_eff_final)
ddpcr_kd_eff_final <- ddpcr_kd_eff_final %>% arrange(foldChange)
ddpcr_kd_eff_final <- ddpcr_kd_eff_final %>% filter(sgRNA %in% all_targets)
ddpcr_kd_eff_final$padj <- p.adjust(ddpcr_kd_eff_final$pvalue, method = "BH")

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
#COMBINE DDPCR AND BULK RNA RESULTS
######################
kd_efficiency_df <- rbind(kd_efficiency_df_rna, ddpcr_kd_eff_final) %>% arrange(foldChange)
row.names(kd_efficiency_df) <- NULL
kd_efficiency_df$sgRNA <- factor(as.character(kd_efficiency_df$sgRNA),levels = as.character(kd_efficiency_df$sgRNA))

######################
#MAKE SUMMARY TABLE
######################
num_degs <- sapply(names(all_DESeq2_res_full), function(x) nrow(all_DESeq2_res_full[[x]] %>% filter(V8 < 0.05)))
kd_efficiency_df$num_degs <- num_degs[as.character(kd_efficiency_df$sgRNA)]

kd_efficiency_df$padj_across_38_targets <- p.adjust(kd_efficiency_df$pvalue, method = "BH")
write.table(kd_efficiency_df, "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Initial_Submission/Figure_1_And_S1/Figure_S1A_KD_Eff_Barplot/kd_efficiency_summary_table.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)

# ######################
# #MAKE PLOT - PVALUE
# ######################
# kd_efficiency_df$sig <- sapply(kd_efficiency_df$pvalue, function(x) ifelse(is.na(x), "Not Applicable", ifelse(x <0.0001, "< 0.0001", ifelse(x < 0.001, "< 0.001", ifelse(x < 0.01, "< 0.01", ifelse(x < 0.05, "< 0.05", "Not Significant"))))))
# kd_efficiency_df$sig <- factor(as.character(kd_efficiency_df$sig), levels = c("< 0.0001", "< 0.001", "< 0.01", "< 0.05", "Not Significant", "Not Applicable"))

# #add number of degs
# num_degs <- sapply(names(all_DESeq2_res_full), function(x) nrow(all_DESeq2_res_full[[x]] %>% filter(V8 < 0.05)))
# kd_efficiency_df$num_degs <- num_degs[as.character(kd_efficiency_df$sgRNA)]
# kd_efficiency_df$sig[which(kd_efficiency_df$sgRNA == "THUMPD3-AS1")] <- "Not Significant"

# #create plot (for supplement)
# cols <- c(brewer.pal(8, "Greens")[c(8, 6, 4, 2)], "lightgray", "white")
# #border_cols <- c(rep("white", length(unique(kd_efficiency_df$sig)) - 1), "black")
# pdf("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Initial_Submission/Figure_1_And_S1/Figure_S1A_KD_Eff_Barplot/knockdown_efficiency_barplot_using_raw_pvalue.pdf", height = 3, width = 10)
# ggplot(data = kd_efficiency_df, aes(x = sgRNA, y = foldChange, fill = sig, width = 0.8)) + geom_bar(stat = "identity", color = "black") + 
# theme_minimal() +
# xlab("Gene") + 
# ylab("Relative Expression of Target") + 
# #theme(legend.position = "none") + 
# theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
# geom_errorbar(aes(ymin=expLowCI, ymax=expHighCI), width=.2, position=position_dodge(.9)) + 
# geom_hline(yintercept = 1, linetype = "dashed") + 
# #geom_text(aes(y=expHighCI + 0.02, label=sig)) + 
# theme(axis.text=element_text(size=7), axis.title=element_text(size=8,face="bold")) + 
# scale_fill_manual(values = cols) + 
# #scale_color_manual(values = border_cols) + 
# guides(fill=guide_legend(title="P-value")) + 
# theme(legend.text=element_text(size=7), legend.title=element_text(size=8,face="bold")) + 
# guides(color = "none")
# dev.off()

# #create plot 
# cols <- c(brewer.pal(8, "Greens")[c(8, 6, 4, 2)], "lightgray", "white")
# all_targets <-  unname(unlist(read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/final_targets_to_carry_forward.csv", header = FALSE)))
# kd_efficiency_df_init <- kd_efficiency_df
# kd_efficiency_df <- kd_efficiency_df %>% filter(sgRNA %in% final_targets)
# pdf("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Initial_Submission/Figure_1_And_S1/Figure_S1A_KD_Eff_Barplot/knockdown_efficiency_barplot_using_raw_pvalue_just_genes_used_for_downstream_analysis.pdf", height = 3, width = 6)
# ggplot(data = kd_efficiency_df, aes(x = sgRNA, y = foldChange, fill = sig, width = 0.8)) + geom_bar(stat = "identity", color = "black") + 
# theme_minimal() +
# xlab("Gene") + 
# ylab("Relative Expression of Target") + 
# #theme(legend.position = "none") + 
# theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
# geom_errorbar(aes(ymin=expLowCI, ymax=expHighCI), width=.2, position=position_dodge(.9)) + 
# geom_hline(yintercept = 1, linetype = "dashed") + 
# #geom_text(aes(y=expHighCI + 0.02, label=sig)) + 
# theme(axis.text=element_text(size=7), axis.title=element_text(size=8,face="bold")) + 
# scale_fill_manual(values = cols) + 
# #scale_color_manual(values = border_cols) + 
# guides(fill=guide_legend(title="P-value")) + 
# theme(legend.text=element_text(size=7), legend.title=element_text(size=8,face="bold")) + 
# guides(color = "none")
# dev.off()

######################
#MAKE PLOT - PADJ
######################
kd_efficiency_df$sig <- sapply(kd_efficiency_df$padj_across_38_targets, function(x) ifelse(is.na(x), "Not Applicable", ifelse(x <0.0001, "< 0.0001", ifelse(x < 0.001, "< 0.001", ifelse(x < 0.01, "< 0.01", ifelse(x < 0.05, "< 0.05", "Not Significant"))))))
kd_efficiency_df$sig <- factor(as.character(kd_efficiency_df$sig), levels = c("< 0.0001", "< 0.001", "< 0.01", "< 0.05", "Not Significant", "Not Applicable"))

#add number of degs
num_degs <- sapply(names(all_DESeq2_res_full), function(x) nrow(all_DESeq2_res_full[[x]] %>% filter(V8 < 0.05)))
kd_efficiency_df$num_degs <- num_degs[as.character(kd_efficiency_df$sgRNA)]
kd_efficiency_df$sig[which(kd_efficiency_df$sgRNA == "THUMPD3-AS1")] <- "Not Significant"

#create plot (for supplement)
cols <- c(brewer.pal(8, "Greens")[c(8, 6, 4, 2)], "lightgray", "white")
#border_cols <- c(rep("white", length(unique(kd_efficiency_df$sig)) - 1), "black")
pdf("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Initial_Submission/Figure_1_And_S1/Figure_S1A_KD_Eff_Barplot/knockdown_efficiency_barplot_using_adj_pvalue.pdf", height = 3, width = 10)
ggplot(data = kd_efficiency_df, aes(x = sgRNA, y = foldChange, fill = sig, width = 0.8)) + geom_bar(stat = "identity", color = "black") + 
theme_minimal() +
xlab("Gene") + 
ylab("Relative Expression of Target") + 
#theme(legend.position = "none") + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
geom_errorbar(aes(ymin=expLowCI, ymax=expHighCI), width=.2, position=position_dodge(.9)) + 
geom_hline(yintercept = 1, linetype = "dashed") + 
#geom_text(aes(y=expHighCI + 0.02, label=sig)) + 
theme(axis.text=element_text(size=7), axis.title=element_text(size=8,face="bold")) + 
scale_fill_manual(values = cols) + 
#scale_color_manual(values = border_cols) + 
guides(fill=guide_legend(title="P-value")) + 
theme(legend.text=element_text(size=7), legend.title=element_text(size=8,face="bold")) + 
guides(color = "none")
dev.off()

#create plot 
cols <- c(brewer.pal(8, "Greens")[c(8, 6, 4, 2)], "lightgray", "white")
all_targets <-  unname(unlist(read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/process_perturbation_data_all_perturbations/final_targets_to_carry_forward.csv", header = FALSE)))
kd_efficiency_df_init <- kd_efficiency_df
kd_efficiency_df <- kd_efficiency_df %>% filter(sgRNA %in% final_targets)
pdf("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/For_Initial_Submission/Figure_1_And_S1/Figure_S1A_KD_Eff_Barplot/knockdown_efficiency_barplot_using_adj_pvalue_just_genes_used_for_downstream_analysis.pdf", height = 3, width = 6)
ggplot(data = kd_efficiency_df, aes(x = sgRNA, y = foldChange, fill = sig, width = 0.8)) + geom_bar(stat = "identity", color = "black") + 
theme_minimal() +
xlab("Gene") + 
ylab("Relative Expression of Target") + 
#theme(legend.position = "none") + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
geom_errorbar(aes(ymin=expLowCI, ymax=expHighCI), width=.2, position=position_dodge(.9)) + 
geom_hline(yintercept = 1, linetype = "dashed") + 
#geom_text(aes(y=expHighCI + 0.02, label=sig)) + 
theme(axis.text=element_text(size=7), axis.title=element_text(size=8,face="bold")) + 
scale_fill_manual(values = cols) + 
#scale_color_manual(values = border_cols) + 
guides(fill=guide_legend(title="P-value")) + 
theme(legend.text=element_text(size=7), legend.title=element_text(size=8,face="bold")) + 
guides(color = "none")
dev.off()