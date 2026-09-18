#!/bin/bash
#SBATCH -c 1                    
#SBATCH -t 1:00:00                    
#SBATCH -p short                
#SBATCH --mem=2GB     
#SBATCH -o %j.out               
#SBATCH -e %j.err               
#SBATCH --mail-type=ALL         
#SBATCH --mail-user=talukdar@mit.edu

#Downsamples CRISPR barcode libraries 
module load  cellranger/7.0.0
sample=$1
input_dir="/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/raw_data/day_60_organoids_scrna/FASTQ"
output_dir="/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/raw_data/day_60_organoids_scrna/Align_Data/0_Downsample_FASTQs/FASTQ"
seqtk="/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/software/seqtk/seqtk"

for file in "$input_dir"/*EGFP*"$sample"*; do 
    filename=$(basename "$file")  
    echo $filename
    output_file="$output_dir/downsampled_$filename"  
    $seqtk sample -s100 "$file" 0.0125 | gzip > $output_file
done 

