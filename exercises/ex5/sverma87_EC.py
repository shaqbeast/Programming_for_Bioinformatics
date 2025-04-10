#!/usr/bin/env python3
import sys

input_blast_file = sys.argv[1]
input_bed_file = sys.argv[2]
output_file = sys.argv[3]

try:
    # Open BLAST and BED files 
    with open(input_blast_file) as blast_file, open(input_bed_file) as bed_file:
        bed_lines = bed_file.readlines()  # Read all BED file lines into memory

        # Open a new file called processed_blast to put all of the BLAST hits that match the below condition
        with open("processed_blast.txt", 'w') as proc_file:
            for blast_line in blast_file:
                blast_list = blast_line.split("\t")
                pident = float(blast_list[2])  
                length = int(blast_list[3])   
                qlength = int(blast_list[12])  

                # Processing BLAST hits with this condition
                if pident > 30 and length / qlength >= 0.9:
                    proc_file.write(blast_line)

    # Now process the filtered BLAST file and compare with the BED file
    with open("processed_blast.txt") as processed_blast_file, open(output_file, 'w') as out_file:
        
        # Loop through the processed blast file
        for proc_blast_line in processed_blast_file:
            proc_blast_line_list = proc_blast_line.split("\t")
            blast_seq = proc_blast_line_list[1]
            sstart = int(proc_blast_line_list[8])  
            send = int(proc_blast_line_list[9])   

            # Loop through the BED file 
            for bed_line in bed_lines:
                bed_line_arr = bed_line.split()
                bed_seq = bed_line_arr[0]
                bed_start = int(bed_line_arr[1])  
                bed_end = int(bed_line_arr[2])    
                gene = bed_line_arr[3]            
                
                # check if the BLAST and bed seq ids match
                if blast_seq != bed_seq:
                    continue

                # MAIN condition to see if the BLASTt hit lies within the bed file
                if sstart >= bed_start and send <= bed_end:
                    out_file.write(gene + "\n")  # Write gene name to output file

except FileNotFoundError:
    print("File not found.")
    sys.exit(1)
    
# must sort Vc_out.txt with this command: sort -u Vc_out.txt -o Vc_out.txt
# this will get rid of any duplicates
# Can display line counts with: wc -l Vc_out.txt

