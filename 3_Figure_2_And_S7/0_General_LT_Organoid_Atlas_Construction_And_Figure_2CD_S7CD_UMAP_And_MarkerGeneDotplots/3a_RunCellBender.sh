#!/bin/bash
#SBATCH -c 1
#SBATCH -N 1   
#SBATCH -t 0-3:0:0 
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=20GB
#SBATCH -o hostname_%j.out            
#SBATCH -e hostname_%j.err         
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=talukdar@mit.edu

#eval "$(conda shell.bash hook)"
#conda activate /home/mt384/.conda/envs/CellBenderEnvWithRC

module purge
module load gcc/9.2.0 cuda/11.7 #prev version was 11.7

# Create output directory
input_path=$(echo "$1" | tr -d '[:space:]')
bname=$(basename $1 | tr -d '[:space:]')
output_dir="/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Organoids/raw_data/day_60_organoids_scrna/Align_Data/3_RunCellBender/Output/"$bname
mkdir -p $output_dir

cellbender remove-background \
                 --input $input_path"/outs/raw_feature_bc_matrix.h5" \
                 --output $output_dir"/cellbender.h5" \
                 --expected-cells 5000 \
                 --cuda \
                 --total-droplets-included 20000 \
                 --fpr 0.01 \
                 --epochs 150

seff %j
