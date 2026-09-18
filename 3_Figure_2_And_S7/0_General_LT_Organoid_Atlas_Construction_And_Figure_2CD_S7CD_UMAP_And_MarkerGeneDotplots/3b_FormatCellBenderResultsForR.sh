#!/bin/bash
#SBATCH -c 1
#SBATCH -N 1   
#SBATCH -t 2:00:00
#SBATCH -p short
#SBATCH --mem=20GB
#SBATCH -o hostname_%j.out            
#SBATCH -e hostname_%j.err         
#SBATCH --mail-type=ALL
#SBATCH --mail-user=talukdar@mit.edu

#conda activate scanpy_env

module purge
wd=$1
cd $wd

ptrepack --complevel 5 cellbender_filtered.h5:/matrix cellbender_filtered_seurat.h5:/matrix

seff %j
