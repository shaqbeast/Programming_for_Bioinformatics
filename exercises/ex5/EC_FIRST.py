#!/usr/bin/env python3
import sys

input_blast_file = sys.argv[1]
input_bed_file = sys.argv[2]
output_file = sys.argv[3]

try:
    blast_file = open(input_blast_file)
    bed_file = open(input_bed_file)
    # keep BLAST hits that are greater than 30% identity and >= 90% length
    with open("processed_blast.txt", 'w') as proc_file:
        for blast_line in blast_file:
            blast_list = blast_line.split("\t")
            pident = blast_list[2] # percentage identity
            length = blast_list[3] # alignment length
            qlength = blast_list[12] # length of query sequence
            if float(pident) > 30 and int(length) / int(qlength) > 0.9:
                    proc_file.write(blast_line)
                
    processed_blast_file = open("processed_blast.txt")
    # sstart and send are present in the 9th and 10th column respectively
    with open(output_file, 'w') as out_file:
        for proc_blast_line in processed_blast_file:
            # get info from blast file (split function will put file line into arr of strings)
            proc_blast_line_list = proc_blast_line.split("\t")
            sstart = proc_blast_line_list[8] # 9 - 1 = 8
            send = proc_blast_line_list[9] # 10 - 1 = 9
            for bed_line in bed_file:
                # get info from bed file
                bed_line_arr = bed_line.split()
                bed_start = bed_line_arr[1]
                bed_end = bed_line_arr[2]
                gene = bed_line_arr[3]
            
                # condition
                if int(sstart) >= int(bed_start) and int(send) <= int(bed_end):
                    print("hello!")
                    out_file.write(gene)
except FileNotFoundError:
    print("File not found.")
    sys.exit(1)