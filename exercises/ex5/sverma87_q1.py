#!/usr/bin/env python3
import sys

# if a line starts with a > char don't save it
# read the file and save it into two variables
# loop through the sequences and see which ones align
# if they align, then add a | symbol

input_file = sys.argv[1]
seq1 = None
seq2 = None
counter = 1 # keeps track of which line we're on in the file
try:
    file = open(input_file)
    for line in file:
        if line[0] != '<' and counter == 2:
            seq1 = line
        if line[0] != '<' and counter == 4:
            seq2 = line
        counter += 1
except FileNotFoundError:
    print("File not found.")
    sys.exit(1)

print(seq1, end="")
index = 0
while index < len(seq1) and index < len(seq2):
    if seq1[index] == seq2[index]:
        print('|', end="")
    else:
        print(" ", end="")
    index += 1
print("")
print(seq2)
    
    