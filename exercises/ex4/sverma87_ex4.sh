#!/bin/zsh

query_file=$1
subject_file=$2
bed_file=$3
output_file=$4
basename="${bed_file%.*}"

# added start and send to show where the homologs start and end 
tblastn -query "$query_file" -subject "$subject_file" -outfmt "6 qseqid sseqid pident length qlen sstart send" -out output.txt

awk '($3 > 30) && ($4 / $5 > 0.9)' output.txt > putative_homologs.txt

echo "" > "$output_file"

# new script for ex4 
while read -r s_line; do
	sseqid=$(echo "$s_line" | cut -f2)
	sstart=$(echo "$s_line" | cut -f6)
	send=$(echo "$s_line" | cut -f7)

	while read -r bed_line; do
		bed_seqid=$(echo "$bed_line" | cut -f1)
		bed_start=$(echo "$bed_line" | cut -f2)
		bed_end=$(echo "$bed_line" | cut -f3)
		
		if [ "$sstart" -ge "$bed_start" ] && [ "$send" -le "$bed_end" ]; then
            		echo "$bed_line" | cut -f4 >> "$output_file"
        	fi
	done < "$bed_file"
done < putative_homologs.txt

sort -u "$output_file" -o "$output_file"
