#!/usr/bin/env python3

import argparse
from collections import defaultdict
import re 
import os
from magnumopus.ispcr import ispcr
from magnumopus.mapping import map_reads_to_ref
import subprocess
from magnumopus.nw import needleman_wunsch
import tempfile

'''SECTION 1'''

def find_reverse_compliment(sequence: str) -> str:
    """ Find the reverse compliment for a particular sequence """
    reverse_compliment = ""
    reverse_seq = sequence[::-1]
    for nucleotide in reverse_seq:
        if nucleotide == "A":
            reverse_compliment += "T"
        elif nucleotide == "T":
            reverse_compliment += "A"
        elif nucleotide == "G":
            reverse_compliment += "C"
        elif nucleotide == "C":
            reverse_compliment += "G"
    
    return reverse_compliment

# Initialize parser
parser = argparse.ArgumentParser(
    "Given zero or more assembly files and zero or more pairs of reads files, utilize" +
    "isPCR, read mapping, and Needleman-Wunsch to extract and align the V4 region of the 16S" +
    "gene using user-provided primers and reference sequences."
)

# Add the arguments
parser.add_argument(
    "-a",
    "--assemblies",
    nargs="+",
    type=str,
    required=False,
    help="Path to assembly files."
)

parser.add_argument(
    "-p",
    "--primers",
    type=str,
    required=True,
    help="Path to primer files."
)

parser.add_argument(
    "-r",
    "--reads",
    nargs="+",
    type=str,
    required=False,
    help="Path to Illumina reads files. Must be given in pairs"
)

parser.add_argument(
    "-s",
    "--ref_seqs",
    type=str, 
    required=False, 
    help="Path of a sequence in the reference FASTA file"
)

# Parse the arguments 
args = parser.parse_args()

assemblies = args.assemblies
primers = args.primers
reads = args.reads
reference = args.ref_seqs

# Do all the regex and glob checks so that any type of input specified in the assignment works
fastq_pattern = re.compile(r"^(.*)_[12]\.fastq$") # pattern that matches the way the Illumina reads will be passed
paired_reads = defaultdict(lambda: {"1": None, "2": None}) # create a default dict that contains the pairs of the seq names
                                                           # default dict automatically creates a new entry when a seq name doesn't exist

for filepath in reads:
    match = fastq_pattern.match(os.path.basename(filepath)) # sees if the pattern matches to the way the file should be passed
    if match:
        sample_name = match.group(1) # the sameple name (e.g. Ecoli)
        if "_1.fastq" in filepath:
            paired_reads[sample_name]["1"] = filepath
        elif "_2.fastq" in filepath:
            paired_reads[sample_name]["2"] = filepath
    else:
        print("Not a valid file path for the Illumina reads. Not using reads for this run.")

# gets the valid pairs that exist within the command line input
valid_pairs = []
for sample_name, reads in paired_reads.items():
    if reads["1"] and reads["2"]:
        valid_pairs.append((reads["1"], reads["2"]))

'''SECTION 2'''
# Perform calculations with the magnumopus methods
# Need 29 results
sequences = {}

'''Assemblies'''
if assemblies:
    for assembly in assemblies:
        amplicons = ispcr(primers, assembly, max_amplicon_size=2000)

        # could also split based ">"
        headers_and_sequences = amplicons.strip().split("\n")
        if headers_and_sequences:
            # Keep only the first amplicon
            header = headers_and_sequences[0]  # The first header
            sequence = headers_and_sequences[1]  # The first sequence
            sequences[header] = sequence
            
'''Reads'''
if valid_pairs:
    for read1, read2 in valid_pairs:
        # get the sam file
        sam = map_reads_to_ref(reference, read1, read2)
        # grabbing consensus sequence and setting fasta=True so that we can get a proper header
        best = sam.best_consensus(fasta=True)
        with tempfile.NamedTemporaryFile(mode='w+', suffix=".fna", delete=False) as best_temp_file:
            best_temp_file.write(best)
            best_temp_path = best_temp_file.name  # Store the file path for later use
        amplicons = ispcr(primers, best_temp_path, max_amplicon_size=2000)
        
        headers_and_sequences = amplicons.strip().split("\n")
        if headers_and_sequences:
                header = headers_and_sequences[0]
                sequence = headers_and_sequences[1]
                sequences[header] = sequence

'''References'''
if reference:
    # need to loop through the entire reference file
    # get each sequence
    # call ispcr on each sequence 
    with open(reference, 'r') as ref_file:
        current_header = None
        current_sequence = []
        
        # loop through the lines of the reference file
        for line in ref_file:
            line = line.strip()
            if line.startswith('>'):
                if current_header and current_sequence:
                    sequence_str = ''.join(current_sequence)
                    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_file:
                        temp_file.write(f"{current_header}\n{sequence_str}")
                        temp_file_path = temp_file.name
                    
                    # call isPCR on the current sequencei n the reference file
                    amplicons = ispcr(primers, temp_file_path, max_amplicon_size=2000)
                    headers_and_sequences = amplicons.strip().split("\n")
                    if headers_and_sequences:
                        header = headers_and_sequences[0]
                        sequence = headers_and_sequences[1]
                        sequences[header] = sequence
                    
                    # need to reset the current sequence
                    current_sequence = []

                # 1st iteration will set automatically set the header to the line
                current_header = line
            else:
                # get the current sequence if we're not at a header 
                current_sequence.append(line)
        
        # final sequence
        if current_header and current_sequence:
            sequence_str = ''.join(current_sequence)
            with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_file:
                temp_file.write(f"{current_header}\n{sequence_str}")
                temp_file_path = temp_file.name
            
            amplicons = ispcr(primers, temp_file_path, max_amplicon_size=2000)
            headers_and_sequences = amplicons.strip().split("\n")
            if headers_and_sequences:
                header = headers_and_sequences[0]
                sequence = headers_and_sequences[1]
                sequences[header] = sequence
    
'''NW'''
amplicons_oriented = {}
first_amplicon = next(iter(sequences.values())) # first value in dictionary (this one is automatically assumed to be oriented correctly)
with open("magop_output.fasta", 'w') as fasta_file:
    for header, amplicon in sequences.items():
        forward_amplicon = amplicon # renaming amplicon to have it make more sense with this code block
        reverse_amplicon = find_reverse_compliment(forward_amplicon)
        
        # Forward Orientation
        _, forward_score = needleman_wunsch(first_amplicon, forward_amplicon, 1, -1, -1)
        
        # Reverse Orientation 
        _, reverse_score = needleman_wunsch(first_amplicon, reverse_amplicon, 1, -1, -1)
        
        # check forward vs reverse scores (whichever is higher is chosen)
        if forward_score > reverse_score:
            fasta_file.write(f"{header}\n")
            fasta_file.write(f"{forward_amplicon}\n")
            print(header)
            print(forward_amplicon)
            print()
        else:
            fasta_file.write(f"{header}\n")
            fasta_file.write(f"{reverse_amplicon}\n")
            print(reverse_amplicon)
            print()

'''SECTION 3'''

'''# Muscle
input_file = "magop_output.fasta"
output_file = "aligned_sequences.aln"
muscle_command = ["muscle", 
                  "-align", input_file, 
                  "-output", output_file]
subprocess.run(muscle_command, check=True)

# IQTree
alignment_file = "aligned_sequences.aln"
output_file = "tree/ml_tree"
iqtree_command = [
            "iqtree2",
            "-s", alignment_file,
            "--prefix", output_file,
            "--quiet",
            "-redo"
            ]
subprocess.run(iqtree_command, check=True)'''

'''# GoTree
input_tree = "tree/ml_tree.treefile"
output_tree = "tree/rooted_tree.treefile"

# Construct the gotree command
gotree_command = [
    "gotree",
    "reroot",
    "midpoint",
    "-i", input_tree,
    "-o", output_tree
]
subprocess.run(gotree_command, check=True)'''
    
    