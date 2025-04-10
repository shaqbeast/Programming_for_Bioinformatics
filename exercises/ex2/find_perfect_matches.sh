#!/bin/zsh

query_file=$1
subject_file=$2
output_file=$3

blastn -query "$query_file" -subject "$subject_file" -outfmt '6 std sseq' -out output.txt

awk '$3 == 100 && $4 == (($10 - $9) + 1)' output.txt > "$output_file"

wc -l "$output_file" 


