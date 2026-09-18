#!/bin/bash
#SBATCH -c 1
#SBATCH -N 1   
#SBATCH -t 5-00:00:00
#SBATCH -p medium
#SBATCH --mem=5GB
#SBATCH -o hostname_%j.out            
#SBATCH -e hostname_%j.err         
#SBATCH --mail-type=ALL
#SBATCH --mail-user=talukdar@mit.edu

module load gcc/14.2.0
module load java/jdk-23.0.1

#get threshold
java -Xmx5G -jar /n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/full_aracne_network_analysis/ARACNe-AP-master/dist/aracne.jar -e /n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/full_aracne_network_analysis/Input/matrix.txt  -o /n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/full_aracne_network_analysis/Output/aracne_output/ --tfs /n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/full_aracne_network_analysis/Input/RegulatorsList/ProcessedLists/all_NPC_expressed_TFs_and_epigenetic_regulators.csv --pvalue 1E-8 --seed 1 --calculateThreshold #regulator list in Input_Files_Not_Generated_By_Scripts

#run aracne
for i in {1..100}
do
java -Xmx5G -jar /n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/full_aracne_network_analysis/ARACNe-AP-master/dist/aracne.jar -e /n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/full_aracne_network_analysis/Input/matrix.txt  -o /n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/full_aracne_network_analysis/Output/aracne_output/  --tfs /n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/full_aracne_network_analysis/Input/RegulatorsList/ProcessedLists/all_NPC_expressed_TFs_and_epigenetic_regulators.csv --pvalue 1E-8 --seed $i
done

seff %j
