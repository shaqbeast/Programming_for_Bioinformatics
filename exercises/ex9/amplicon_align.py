#!/usr/bin/env python3
import magnumopus
import argparse

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
    description="Perform in-silico PCR on two assemblies and align the amplicons"
)

# Add arguments
parser.add_argument(
    "-1", 
    dest="assembly1",
    type=str, 
    required=True, 
    help="Path to the first assembly file"
)
parser.add_argument(
    "-2", 
    dest="assembly2",
    type=str, 
    required=True, 
    help="Path to the second assembly file"
)
parser.add_argument(
    "-p", 
    dest="primers",
    type=str, 
    required=True, 
    help="Path to the primer file"
)
parser.add_argument(
    "-m",
    dest="max_amplicon_size",
    type=int, 
    required=True, 
    help="Maximum amplicon size for in-silico PCR"
)
parser.add_argument(
    "--match", 
    type=int, 
    required=True, 
    help="Match score to use in alignment"
)
parser.add_argument(
    "--mismatch", 
    type=int, 
    required=True, 
    help="Mismatch penalty to use in alignment"
)
parser.add_argument(
    "--gap", 
    type=int, 
    required=True, 
    help="Gap penalty to use in alignment"
)

args = parser.parse_args()

amplicons1 = magnumopus.ispcr(args.primers, args.assembly1, args.max_amplicon_size)
print(amplicons1)
amplicons2 = magnumopus.ispcr(args.primers, args.assembly2, args.max_amplicon_size)
print(amplicons2)
print()
seq1 = amplicons1.strip().split("\n")[1] # get rid of header 
seq2 = amplicons2.strip().split("\n")[1] # get rid of header
seq2 = find_reverse_compliment(seq2)

aligned, score = magnumopus.needleman_wunsch(seq1, seq2, args.match, args.mismatch, args.gap)

aligned1, aligned2 = aligned

print(aligned1)
print(aligned2)
print(score)


