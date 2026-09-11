#!/bin/bash
#SBATCH -c 1
#SBATCH -N 1   
#SBATCH -t 4:00:00
#SBATCH -p short
#SBATCH --mem=20GB
#SBATCH -o hostname_%j.out            
#SBATCH -e hostname_%j.err         
#SBATCH --mail-type=ALL
#SBATCH --mail-user=talukdar@mit.edu

{
echo input_for_MAGIC.csv
echo 1
} | python3 /n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/magic_analysis/MAGIC_1_1_updated.py

