#this is a sample shell script to submit to a Slurm cluster. 
#Doing this is optional, you can also run this script directly using Rscript, 
#but a Slurm cluster allows you to schedule jobs more efficiently. 

#!/bin/bash
#SBATCH -c 1
#SBATCH -N 1   
#SBATCH -t 0:10:00
#SBATCH -p short
#SBATCH --mem=2GB
#SBATCH -o hostname_%j.out            
#SBATCH -e hostname_%j.err         
#SBATCH --mail-type=ALL
#SBATCH --mail-user=####YOUR EMAIL 

module load gcc #note: may differ based on system 
module load R #note may differ based on system 

Rscript 0_sample_script_to_execute.R
