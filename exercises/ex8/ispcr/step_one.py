#!/usr/bin/env python3

import subprocess

def run_blast(primer_file: str, assembly_file:str) -> str:
    # formats blast command into readable command
    blast_cmd = [
        "blastn",
        "-task", "blastn-short", 
        "-query", primer_file,
        "-subject", assembly_file,
        "-outfmt", "6 std qlen"  # Corrected format here
    ]

    # subprocess.run runs our command
    # capture_output = True allows us to capture the output
    # text = True turn into a string rather than having it as class byte
    result = subprocess.run(blast_cmd, capture_output=True, text=True)
    
    # must capture the standard out of the blast command
    blast_result = result.stdout
    
    return blast_result

def process_hits(blast_result: str, percent_hits: float) -> list[list[str]]:
    # get the lines of the BLAST output into an array
    # each index is a line of the blast result
    blast_lines = blast_result.strip().splitlines()
    processed_hits = list()
    
    for line in blast_lines:
        # go through each line
        # split the line into columns separated by \t
        columns = line.split("\t")
        percent_alignment = float(columns[2])
        alignment_length = int(columns[3])
        
        if (percent_alignment >= percent_hits) and (alignment_length == 20 or alignment_length == 19):
            processed_hits.append(columns)
            
    processed_hits.sort(key=lambda x: int(x[8]))
    
    return processed_hits
    

def step_one(primer_file: str, assembly_file: str) -> list[list[str]]:
    result = run_blast(primer_file, assembly_file)
    hits = process_hits(result, 80)
    
    return hits


