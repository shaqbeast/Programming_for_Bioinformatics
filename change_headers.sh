#!/bin/zsh

input_file=$1
output_file=$2

accession_number=$(basename "$input_file" .fna)

sed "s/>/>${accession_number} /g" "$input_file" > "$output_file"
