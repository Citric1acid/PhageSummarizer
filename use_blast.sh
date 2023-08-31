#!/bin/bash
# usage:
# ./use_blast.sh input output

name=$1
output=$2
outfmt="6 qseqid sseqid stitle pident evalue length sstart send"

blastn -query ${name} -db phage_db -out ${output} -outfmt ${outfmt} -max_target_seqs 10 -num_threads 8

echo "BLAST: ${name} done!"
