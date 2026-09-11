#!/bin/bash
#SBATCH -p medium
#SBATCH --mem=10GB
#SBATCH -t 24:00:00
#SBATCH -c 20
#SBATCH -N 1   
#SBATCH -o withDedup.out        
#SBATCH -e withDedup.err 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=talukdar@mit.edu

#note: script here shown for file 1B - works with any CheapSeq data input (see GEO)
/n/groups/walsh/indData/Maya/neural_lncRNA_CRISPR_screen/CheapSeqPilot/CheapSeq/Software/salmon-v1.9/bin/salmon_v1.9 alevin -l ISR \
-1 /n/groups/walsh/indData/Maya/neural_lncRNA_CRISPR_screen/CheapSeqPerturbations/PreprocessingForAlevin/PreprocessedFiles/1B/Processed_Read_Files/CW_REA_SCR_1B_UMI_IND.fastq.gz \
-2 /n/groups/walsh/indData/Maya/neural_lncRNA_CRISPR_screen/CheapSeqPerturbations/PreprocessingForAlevin/PreprocessedFiles/1B/Processed_Read_Files/CW_REA_SCR_1B_Transcript.fastq.gz \
-i /n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/run_alignment/gencode_transcripts_salmon_index \
-p 20 \
-o /n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/run_alignment/1B/Output \
--bc-geometry 1[1-16] \
--umi-geometry 1[17-24] \
--read-geometry 2[1-end] \
--dumpMtx \
--whitelist /n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/run_alignment/Combinatorial_Index.txt  \
--tgMap /n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/raw_data/run_alignment/txp2gene.tsv
