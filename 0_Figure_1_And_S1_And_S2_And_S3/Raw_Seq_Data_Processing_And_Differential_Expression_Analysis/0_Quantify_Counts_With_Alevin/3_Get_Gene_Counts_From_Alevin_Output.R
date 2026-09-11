#PURPOSE: used to convert the output of Alevin on dedup mode into a human-readable format for each screen sample (ex: SCR 1A)

paths <- unlist(read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/run_alignment/pathsToAlevinOutput_withDedup.txt", header = FALSE))

for (path in paths)
{
print(path) 

#read in data
library(Matrix)
mat <- readMM(paste0(path, "/quants_mat.mtx.gz"))
rownames <- unlist(read.table(paste0(path, "/quants_mat_rows.txt")))
colnames <- unlist(read.table(paste0(path, "/quants_mat_cols.txt")))

#format counts matrix
row.names(mat) <- rownames
colnames(mat) <- colnames
mat <- as.data.frame(t(as.matrix(mat)))
write.table(mat, paste0(path,"/quants_mat.txt"), col.names = TRUE, row.names = TRUE)

#merge counts from transcripts originating from the same gene
mat$gene <- sapply(colnames, function(x) strsplit(x, "|", fixed = TRUE)[[1]][[6]])
library(tidyverse)
mat_merged_by_gene <- as.data.frame(mat %>% group_by(gene) %>% summarize_each(list(sum)))
row.names(mat_merged_by_gene) <- mat_merged_by_gene$gene
mat_merged_by_gene$gene <- NULL
write.table(mat_merged_by_gene, paste0(path,"/quants_mat_merged_by_gene.txt"), col.names = TRUE, row.names = TRUE)
}
