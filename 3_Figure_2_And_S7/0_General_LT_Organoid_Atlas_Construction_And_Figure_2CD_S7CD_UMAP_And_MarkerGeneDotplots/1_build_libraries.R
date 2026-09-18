##################
#I/O
##################
library(tidyverse)
library(ggplot2)

template_df <- read.csv("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/raw_data/day_60_organoids_scrna/Align_Data/1_BuildLibraries/template_libraries.csv", header = TRUE) #In Input_Files_Not_Generated_By_Scripts
output_dir <- "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/raw_data/day_60_organoids_scrna/Align_Data/1_BuildLibraries/Output"
samples <- unname(unlist(read.table("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/raw_data/day_60_organoids_scrna/Align_Data/1_BuildLibraries/sample_names.txt", header = FALSE))) #In Input_Files_Not_Generated_By_Scripts

for (sample in samples)
{
	current_df <- template_df
	current_df$sample <- sapply(current_df$sample, function(x) gsub("AD1-1", sample, x))
	write.table(current_df, paste0(output_dir, "/libraries_", sample, ".csv"), row.names = FALSE, col.names = TRUE, sep = ",", quote = FALSE)
}

