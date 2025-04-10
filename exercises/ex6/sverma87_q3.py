#!/usr/bin/env python3
import sys

blastfile = sys.argv[1]
bedfile = sys.argv[2]
assemblyfile = sys.argv[3]
outputfile = sys.argv[4]

# read blast file
hits = []
with open(blastfile) as fin:
    for line in fin: # .readlines() is the default iter method for the open file class
        
        # unpack and convert types of desired columns. This is ugly. We'll revisit later...
        _, sid, pcnt, matchlen, _, _, _, _, sstart, send, _, _, qlen = line.split()
        pcnt = float(pcnt)
        matchlen = int(matchlen)
        sstart = int(sstart)
        send = int(send)
        qlen = int(qlen)
    
        # Keep hits that could be homologs
        if pcnt > 30 and matchlen > 0.9*qlen:
            # We could store matches as a list or tuple.
            # We won't want to modify the elements so a tuple is "safer" in that we then can't modify it by mistake
            hits.append((sid, sstart, send))

# Now read the bed file
feats = []
with open(bedfile) as fin:
    for line in fin:
        bed_sid, bed_start, bed_end, gene, score, directionality = line.split() # an asterisk before a variable name when unpacking makes that variable store remaining elements
        bed_start = int(bed_start)
        bed_end = int(bed_end)
        
        feats.append((bed_sid, bed_start, bed_end, gene, directionality))

assembly_dict = {}
with open(assemblyfile) as assembly:
    current_header = ""
    current_sequence = []
    for line in assembly:
        line = line.strip()
        if line.startswith(">"):  # If it's a new header
            if current_header:
                assembly_dict[current_header] = ''.join(current_sequence)
            current_header = line[1:]  # Store the header without '>'
            current_sequence = []
        else:
            current_sequence.append(line)
    # Add the last sequence
    if current_header:
        assembly_dict[current_header] = ''.join(current_sequence)

# Now we have our two datasets read in, we can loop over them to find matches
homologs = []
for blast_sid, blast_sstart, blast_send in hits: # unpack our blast data
    for bed_sid, bed_start, bed_end, gene, directionality in feats:
        # Don't bother checking the rest if the sid doens't match
        if blast_sid != bed_sid:
            continue
        
        # Once we are dealing with features at higher index locations than our hit, go to the next hit (break loop over feats)
        if blast_sstart <= bed_start or blast_send <= bed_start:
            break
        
        # Otherwise, check if the hit is inside the feature
        if (blast_sstart > bed_start
            and blast_sstart <= bed_end
            and blast_send > bed_start
            and blast_send <= bed_end
        ):
            homologs.append(gene)
            break # Each BLAST hit will only be in one feature so move to next hit once you've found it

# Get the unique homologs using a set()
unique_homologs = set(homologs)

# open the bed file
# extract all info from bed file
# use the start and end variables to search within the .fna file to get the sequence 

# finds the reverse compliment of a particular seu
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

# put assembly file into a dictionary
assembly_dict = {}
with open(assemblyfile) as assembly:
    current_header = ""
    current_sequence = []
    for line in assembly:
        line = line.strip()
        if line.startswith(">"):  # If it's a new header
            if current_header:
                assembly_dict[current_header] = ''.join(current_sequence)
            current_header = line[1:].split()[0]  # Store the header without '>' and without everything afterwards
            current_sequence = []
        else:
            current_sequence.append(line)
    # Add the last sequence
    if current_header:
        assembly_dict[current_header] = ''.join(current_sequence)

sequence = ""
with open(outputfile, "w") as outfile:
    # go through all unique homologs
    for homolog_gene in unique_homologs:
        outfile.write(">" + str(homolog_gene) + "\n")
        # go through bed file
        for bed_sid, bed_start, bed_end, gene, directionality in feats:
            if homolog_gene == gene:
                # gets the sequence from the assembly file that's within the valid range of the gene
                sequence = assembly_dict.get(bed_sid, '')[bed_start:bed_end]
                if directionality == "-":
                    reverse_sequence = find_reverse_compliment(sequence)
                    outfile.write(reverse_sequence + "\n")
                else:
                    outfile.write(sequence + "\n")
        


        
    
    
    
    
    
