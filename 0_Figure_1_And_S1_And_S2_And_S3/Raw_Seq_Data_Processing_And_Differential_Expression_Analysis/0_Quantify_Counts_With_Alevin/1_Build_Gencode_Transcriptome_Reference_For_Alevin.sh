#download reference files from gencode
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.transcripts.fa.gz

#create salmon index
/n/groups/walsh/indData/Maya/neural_lncRNA_CRISPR_screen/CheapSeqPilot/CheapSeq/Software/salmon-v1.9/bin/salmon_v1.9 index -t gencode.v45.transcripts.fa.gz -i gencode_transcripts_salmon_index -p 12 #or however many cores you can allocate

#create txp2gene file
zcat gencode.v45.transcripts.fa.gz | grep ">" > headers.txt
cut -c 3- headers.txt > headers_cleaned.txt
paste <(cut -f1 headers_cleaned.txt) headers_cleaned.txt > txp2gene.tsv
