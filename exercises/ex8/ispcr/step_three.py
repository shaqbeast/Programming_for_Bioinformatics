import subprocess

def create_bed_file(hit_pairs: list[tuple[list[str]]]) -> str:
    bed_array = []
    for pair in hit_pairs:
        hit1, hit2 = pair
        
        contig = hit1[1]
        amplicon_start = hit1[9]
        amplicon_end = int(hit2[9]) - 1
        
        bed_line = contig + "\t" + amplicon_start + "\t" + str(amplicon_end) + "\n"
        bed_array.append(bed_line)
    
    bed_file = "data/Vibrio_paired_hits.bed"
    
    with open(bed_file, 'w') as file:
        file.writelines(bed_array)
    
    return bed_file
    
    

def step_three(hit_pairs: list[tuple[list[str]]], assembly_file: str) -> str:
    bed = create_bed_file(hit_pairs)
    
    seqtk_cmd = [
        "seqtk", "subseq",
        f"{assembly_file}",
        f"{bed}"
    ]
    
    result = subprocess.run(seqtk_cmd, capture_output=True, text=True)
    
    result_extraction = result.stdout
        
    return result_extraction