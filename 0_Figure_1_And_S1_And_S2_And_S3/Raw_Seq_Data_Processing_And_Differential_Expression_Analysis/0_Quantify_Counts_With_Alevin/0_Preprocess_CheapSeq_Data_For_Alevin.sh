#example of slurm parameters used
#!/bin/bash
#SBATCH -c 20
#SBATCH -N 1   
#SBATCH -t 4:00:00
#SBATCH -p short
#SBATCH --mem=250GB
#SBATCH -o %j.out        
#SBATCH -e %j.err 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=talukdar@mit.edu

#define command line args
WORKINGDIR=$1 #note: all fastqs can be found on GEO (GSE305954)
OUTDIR=$2 #example: /n/groups/walsh/indData/Maya/neural_lncRNA_CRISPR_screen/CheapSeqPerturbations/PreprocessingForAlevin/PreprocessedFiles

#create output folders
mkdir $OUTDIR
mkdir $OUTDIR/CB_UMI_Extraction
mkdir $OUTDIR/CB_UMI_Concatenation
mkdir $OUTDIR/Transcript_Trim
mkdir $OUTDIR/UMI_Merge
mkdir $OUTDIR/IND_Merge
mkdir $OUTDIR/Processed_Read_Files

#load necessary software
export PATH=/n/groups/walsh/indData/Maya/neural_lncRNA_CRISPR_screen/CheapSeqPilot/CheapSeq/SOFTWARE/bbmap:$PATH 
module load gcc/6.2.0 
module load java/jdk-11.0.11

#remove _R1.fastq.gz suffix
for i in $(ls $WORKINGDIR | grep "R1.fastq.gz" | rev | cut -c 12- | rev )

do

#using bbduk, extract index 1
bbduk.sh in=$WORKINGDIR/${i}R1.fastq.gz out=$OUTDIR/CB_UMI_Extraction/${i}IND1.fastq ftl=8 ftr=15 minlength=8 &
done 
wait

#remove _R2.fastq.gz suffix
for i in $(ls $WORKINGDIR | grep "R2.fastq.gz" | rev | cut -c 12- | rev )

do

#using bbduk, extract UMI2 (which is the TSO-associated UMI)
bbduk.sh in=$WORKINGDIR/${i}R2.fastq.gz out=$OUTDIR/CB_UMI_Extraction/${i}UMI2.fastq ftl=0 ftr=7 minlength=8 &
#using bbduk, extract index2
bbduk.sh in=$WORKINGDIR/${i}R2.fastq.gz out=$OUTDIR/CB_UMI_Extraction/${i}IND2.fastq ftl=8 ftr=15 minlength=8 &
done
wait

#repair index (IND1 and IND2) read pairs so that they are ordered  (this can take a lot of memory, so I do not fork the processes)
for i in $(ls $OUTDIR/CB_UMI_Extraction | grep "fastq" | rev | cut -c 12- | rev | uniq)
do
repair.sh  in1=$OUTDIR/CB_UMI_Extraction/${i}_IND1.fastq in2=$OUTDIR/CB_UMI_Extraction/${i}_IND2.fastq out1=$OUTDIR/IND_Merge/${i}_IND1.fastq out2=$OUTDIR/IND_Merge/${i}_IND2.fastq repair int=f
done

#append individual IND sequences together (IND1 + IND2)
for i in $(ls $OUTDIR/CB_UMI_Extraction | grep "fastq" | rev | cut -c 12- | rev | uniq)
do
paste -d '\n' $OUTDIR/IND_Merge/${i}_IND1.fastq $OUTDIR/IND_Merge/${i}_IND2.fastq | sed -n 'p;n;n;N;s/\n//p' > $OUTDIR/CB_UMI_Concatenation/${i}_IND_merge.fastq

#repair merged indexes with UMI2 (output called UMI_merge as a hold-over from earlier versions that combined UMI1+UMI2)
#again, repairing the reads  can take a lot of memory, so I don't fork the processes.
repair.sh  in1=$OUTDIR/CB_UMI_Concatenation/${i}_IND_merge.fastq in2=$OUTDIR/CB_UMI_Extraction/${i}_UMI2.fastq out1=$OUTDIR/CB_UMI_Concatenation/${i}_IND_merge_clean.fastq out2=$OUTDIR/CB_UMI_Concatenation/${i}_UMI_merge_clean.fastq repair int=f
done

#append IND and UMI sequences together (IND1+IND2+UMI2)
for i in $(ls $OUTDIR/CB_UMI_Extraction | grep "fastq" | rev | cut -c 12- | rev | uniq)
do
paste -d '\n' $OUTDIR/CB_UMI_Concatenation/${i}_IND_merge_clean.fastq $OUTDIR/CB_UMI_Concatenation/${i}_UMI_merge_clean.fastq | sed -n 'p;n;n;N;s/\n//p' > $OUTDIR/CB_UMI_Concatenation/${i}_IND_UMI_merge.fastq &
done
wait

#trim and extract transcript read
#first line trims TSO sequence, second line trims polyA
for i in $(ls $WORKINGDIR | grep "fastq" | rev | cut -c 13- | rev | uniq)
do
bbduk.sh in=$WORKINGDIR/${i}_R2.fastq.gz out=stdout.fq literal=GAAGTGCCTCATGGG,TGGGACACTCATGGG,ATAACAGGTCATGGG,CTGCCCTTTCATGGG,GGTCTTATTCATGGG,TACCAACTTCATGGG,ATTGGTTCTCATGGG,CTCTGGGATCATGGG ktrim=l k=14 hdist=1 int=f| \
bbduk.sh in=stdin.fq out=$OUTDIR/Transcript_Trim/${i}_Transcript.fastq.gz literal=AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA k=10 hdist=1 ktrim=r int=f &
done
wait

#repair transcript and IND_UMI reads (this can take a lot of memory, so I do not fork the processes)
for i in $(ls $OUTDIR/CB_UMI_Extraction | grep "fastq" | rev | cut -c 12- | rev | uniq)
do
repair.sh  in1=$OUTDIR/Transcript_Trim/${i}_Transcript.fastq.gz in2=$OUTDIR/CB_UMI_Concatenation/${i}_IND_UMI_merge.fastq out1=$OUTDIR/Processed_Read_Files/${i}_UMI_IND.fastq.gz out2=$OUTDIR/Processed_Read_Files/${i}_Transcript.fastq.gz repair int=f
done
wait

