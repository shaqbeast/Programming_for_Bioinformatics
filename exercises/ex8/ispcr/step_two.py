#!/usr/bin/env python3

import subprocess

def step_two(sorted_hits: list[str], max_amplicon_size: int) -> list[tuple[list[str]]]:
    amplicon_pairs = []
    
    for index1 in range(len(sorted_hits)):
        hit1 = sorted_hits[index1]
        for index2 in range(index1 + 1, len(sorted_hits)):
            hit2 = sorted_hits[index2]
            
            # get positions where primer matches
            hit1_5prime = int(hit1[8])
            hit1_3prime = int(hit1[9])
            hit2_5prime = int(hit2[8])
            hit2_3prime = int(hit2[9])
            
            # condition checks whether the 5' of the first hit < 3' of the second hit 
            # checks if subtracting the 3 prime ends is less than the amplicon size
            if (hit1_5prime < hit2_3prime) and (hit2_3prime - hit1_3prime <= max_amplicon_size):
                pair = (hit1, hit2)
                amplicon_pairs.append(pair)
            
            
        
    return amplicon_pairs