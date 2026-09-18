#!/bin/bash

#Variables
READ1=$1
READ2=$2
OUTDIR=$3
WHITELIST=$4


# Ensure the output directory exists
mkdir -p "$OUTDIR/WorkingDir"
mkdir -p "$OUTDIR/Output"
mkdir -p "$OUTDIR/Log"

# Step 1: Extract UMIs with umi_tools
umi_tools extract --bc-pattern=NNNNNNNNNNNNNNNN \
                  --bc-pattern2=XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXCCCCCCCCCC \
                  --stdin "$READ1" \
                  --read2-in "$READ2" \
                  --stdout "$OUTDIR/WorkingDir/R1_Extract.fastq.gz" \
                  --read2-out "$OUTDIR/WorkingDir/R2_Extract.fastq.gz" \
                  --whitelist "$WHITELIST" \
                  --error-correct-cell

# Step 2: Process Extracted FASTQ
zcat < "$OUTDIR/WorkingDir/R1_Extract.fastq.gz" | \
awk 'BEGIN {OFS="\t"} 
/^@/ {
    # Split the line on spaces
    split($1, fields, "_")
    # Extract read name, CBC, and UMI
    split(fields[1], name, "@")
    read_name = name[2]
    CBC = fields[2]
    UMI = fields[3]
    # Print reformatted output
    print read_name "_" UMI, CBC
}' | sort -k2,2 | gzip -c > "$OUTDIR/WorkingDir/output.txt.gz"

# Step 3: Count Genes with umi_tools
umi_tools count_tab \
    --stdin="$OUTDIR/WorkingDir/output.txt.gz" \
    --stdout="$OUTDIR/Output/gene_counts.tsv" \
    --log="$OUTDIR/Log/umi_tools_count.log"

