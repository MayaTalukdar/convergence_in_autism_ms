#!/bin/bash
#SBATCH -c 4 
#SBATCH --mem=150GB                  
#SBATCH -t 12:00:00                    
#SBATCH -p short                
#SBATCH -o %j.out               
#SBATCH -e %j.err               
#SBATCH --mail-type=ALL         
#SBATCH --mail-user=talukdar@mit.edu


#conda activate organoid_screen_processing_env
module load bowtie

/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/raw_data/day_15_day_21_organoids_121524/Scripts/F5_UMI_Collapse_TopScript.sh \
/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/raw_data/day_15_day_21_organoids_121524/Scripts/ \
/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/raw_data/day_15_day_21_organoids_121524/FASTQs \
/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/raw_data/day_15_day_21_organoids_121524/Output \
/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/raw_data/day_15_day_21_organoids_121524/F5_Whitelist.txt \
4

seff %j
