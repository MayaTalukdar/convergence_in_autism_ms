#!/bin/bash

# Variables
SCRIPT=$1  # Main script path
SUBSCRIPT="$SCRIPT/F5_UMI_Collapse_MainScript.sh"  # Path to your subscript
INPUTFOLDER=$2  # Input folder containing FASTQ files
OUTDIR=$3  # Output directory
WHITELIST=$4  # Path to whitelist file
CORES=$5  # Number of parallel cores
FILE_SET="$OUTDIR/File_Set.txt"  # Path to File_Set.txt

# Ensure the output directory exists
mkdir -p "$OUTDIR"

# Initialize the File_Set.txt
echo -e "Input_R1\tInput_R2\tOutput_Base" > "$FILE_SET"

# Loop through all R1 files in the input directory
for R1_FILE in "$INPUTFOLDER"/*_R1_*.fastq.gz; do
    # Ensure the file exists
    if [[ ! -e "$R1_FILE" ]]; then
        echo "No R1 files found in $INPUTFOLDER."
        exit 1
    fi
    
    # Derive the corresponding R2 file
    R2_FILE=${R1_FILE/_R1_/_R2_}
    
    # Check if the corresponding R2 file exists
    if [[ ! -e "$R2_FILE" ]]; then
        echo "Matching R2 file for $R1_FILE not found."
        exit 1
    fi
    
    # Extract the base name
    BASENAME=$(basename "$R1_FILE" | sed -E 's/_R1_.*//')
    
    # Append the input and output information to the file
    echo -e "${R1_FILE}\t${R2_FILE}\t${OUTDIR}/${BASENAME}" >> "$FILE_SET"
done

echo "File_Set.txt generated at $FILE_SET"

# Run in parallel
tail -n +2 "$OUTDIR/File_Set.txt" | while IFS=$'\t' read -r READ1 READ2 OUTPUT_BASE; do
  printf "%q " "$SUBSCRIPT" "$READ1" "$READ2" "$OUTPUT_BASE" "$WHITELIST"
  echo
done | parallel -j "$CORES"

