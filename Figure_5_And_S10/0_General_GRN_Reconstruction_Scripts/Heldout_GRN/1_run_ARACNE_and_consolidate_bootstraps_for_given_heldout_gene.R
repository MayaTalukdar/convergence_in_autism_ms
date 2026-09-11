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

heldout_gene=$1
echo $heldout_gene
input_dir="/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/half_heldout_aracne_network_analysis/Input/"$heldout_gene"/matrix.txt"
output_dir="/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/half_heldout_aracne_network_analysis/Output/"$heldout_gene"/All_TFs_And_Epigenetic_Regulators/aracne_output/"
mkdir -p $output_dir

#get threshold
java -Xmx5G -jar /n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/half_heldout_aracne_network_analysis/ARACNe-AP-master/dist/aracne.jar -e $input_dir  -o $output_dir --tfs /n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/half_heldout_aracne_network_analysis/Input/RegulatorsList/ProcessedLists/all_NPC_expressed_TFs_and_epigenetic_regulators.csv --pvalue 1E-8 --seed 1 --calculateThreshold

#run aracne
for i in {1..100}
do
java -Xmx5G -jar /n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/half_heldout_aracne_network_analysis/ARACNe-AP-master/dist/aracne.jar -e $input_dir  -o $output_dir  --tfs /n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/half_heldout_aracne_network_analysis/Input/RegulatorsList/ProcessedLists/all_NPC_expressed_TFs_and_epigenetic_regulators.csv --pvalue 1E-8 --seed $i
done

java -Xmx5G -jar /n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/half_heldout_aracne_network_analysis/ARACNe-AP-master/dist/aracne.jar -o $output_dir --consolidate

seff %j

