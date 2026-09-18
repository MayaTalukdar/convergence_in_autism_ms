#!/bin/bash
#SBATCH -c 8                    
#SBATCH -t 12:00:00                    
#SBATCH -p short                
#SBATCH --mem-per-cpu=8G        
#SBATCH -o %j.out               
#SBATCH -e %j.err               
#SBATCH --mail-type=ALL         
#SBATCH --mail-user=talukdar@mit.edu

module load  cellranger/7.0.0

sample=$1

echo $sample

library_file="/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/raw_data/day_60_organoids_scrna/Align_Data/1_BuildLibraries/Output/libraries_"$sample".csv"
cellranger count --id=$sample --libraries=$library_file --feature-ref=/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/raw_data/day_60_organoids_scrna/Align_Data/2_PerformAlignment/Feature_References.csv --transcriptome=/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/raw_data/day_30_organoids_scrna/Align_Data/0_BuildReference/GRCh38 --include-introns=false #From cellranger website

seff %j
