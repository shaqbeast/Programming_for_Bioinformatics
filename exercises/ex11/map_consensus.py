#!/usr/bin/env python3

import argparse
import subprocess
from magnumopus.sam import SAM
import sys

def run_minimap2(reference, read1, read2, output_sam):
    command = [
        "minimap2",
        "-ax", "sr",
        "-B", "0",
        "-k", "10",
        reference,
        read1,
        read2,
    ]
    
    try:
        with open(output_sam, "w") as sam_file:
            subprocess.run(
                command,
                stdout=sam_file,  # Write minimap2 output directly to aligned.sam
                stderr=subprocess.PIPE,
                text=True,
                check=True,
            )
    except subprocess.CalledProcessError as e:
        print(f"Error running minimap2: {e.stderr}", file=sys.stderr)
        sys.exit(1)



# Initialize parser 
parser = argparse.ArgumentParser(
    description="Print a consensus sequence to the terminal in FASTA format"
)

# Add arguments
parser.add_argument(
    "-1", 
    "--read1",
    type=str, 
    required=True, 
    help="Path to a read file in FASTQ format"
)
parser.add_argument(
    "-2", 
    "--read2",
    type=str, 
    required=True, 
    help="Path to a second read file in FASTQ format"
)
parser.add_argument(
    "-r", 
    "--reference",
    type=str, 
    required=True, 
    help="Path to the reference file"
)
parser.add_argument(
    "-s",
    "--seqname",
    type=str, 
    required=False, 
    help="Name of a sequence in the reference FASTA file"
)

args = parser.parse_args()

# Paths
ref_path = args.reference
read1_path = args.read1
read2_path = args.read2
seq_name = args.seqname
output_sam = "aligned.sam"

run_minimap2(ref_path, read1_path, read2_path, output_sam)

# Parse same file
sam = SAM.from_sam(output_sam)

# boths if and else statements are used when user gives a seq_name
if seq_name:
    consensus_seq = sam.consensus(seq_name)
    if not consensus_seq:
        print(f"No reads mapped to the sequence {seq_name}.", file=sys.stderr)
        sys.exit(1)
    print(f">{seq_name}_consensus")
    print(consensus_seq)
else:
    consensus_seq = sam.best_consensus()
    if not consensus_seq:
        print("No reads mapped to the reference. Exiting.", file=sys.stderr)
        sys.exit(1)
    print(consensus_seq)
