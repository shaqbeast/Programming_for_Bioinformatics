class Read:
    def __init__(self, sam_line):
        '''Creates instance of Read Class'''
        fields = sam_line.strip().split("\t")
        
        self.qname = fields[0]
        self.flag = int(fields[1])
        self.rname = fields[2]
        self.pos = int(fields[3])
        self.mapq = int(fields[4])
        self.cigar = fields[5]
        self.rnext = fields[6]
        self.pnext = int(fields[7])
        self.tlen = int(fields[8])
        self.seq = fields[9]
        self.qual = fields[10]


    @property
    def is_mapped(self) -> bool:
        '''Returns whether the read mapped to the reference'''
        binary_flag = f"{self.flag:012b}"
        mapped = False
        # position 9 indicates the flag attribute "Read unmapped"
        if binary_flag[9] == '0':
            mapped = True
            
        return mapped
    
    @property
    def is_forward(self) -> bool:
        '''Returns whether the read mapped in the forward direction'''
        binary_flag = f"{self.flag:012b}"
        print(binary_flag)
        forward = False
        # position 6 indicates flag attribute "Reverse"
        if binary_flag[7] == '0':
            forward = True
            
        return forward
        
    @property
    def is_reverse(self) -> bool:
        '''Returns whether the read mapped in the reverse direction'''
        binary_flag = f"{self.flag:012b}"
        reverse = False
        # position 6 indicates flag attribute "Reverse"
        if binary_flag[7] == '1':
            reverse = True
            
        return reverse
    
    @property
    def is_primary(self) -> bool:
        '''Returns whether the SAM entry is the primary mapping of read'''
        binary_flag = f"{self.flag:012b}"
        if binary_flag[0] == '1':
            return False
        
        primary = False
        # position 3 indicates flag attribute "Secondary"
        if binary_flag[3] == '0':
            primary = True
        
        return primary
    
    def base_at_pos(self, pos: int) -> str:
        if self.cigar == "*":
            return ""
        
        pos -= 1 
        read_pos = 0
        ref_pos = self.pos - 1
        base = ""
        
        i = 0
        # each iteration will give us a number-letter pair from the cigar string
        while i < len(self.cigar):
            number = "" # number in the num-letter pair
            while i < len(self.cigar) and self.cigar[i].isdigit():
                number += self.cigar[i]
                i += 1
            number = int(number)
            letter = self.cigar[i] # letter in num-letter pair 
            i += 1
            
            if letter == "M":
                for _ in range(number):
                    if ref_pos == pos:
                        base += self.seq[read_pos]
                    read_pos += 1
                    ref_pos += 1
            elif letter == "I":
                if ref_pos == pos + 1:
                    base += self.seq[read_pos]
                read_pos += number
            elif letter == "D":
                ref_pos += number
            elif letter == "S":
                read_pos += number
        
        return base
    
    def qual_at_pos(self, pos: int) -> str:
        pos -= 1 
        read_pos = 0
        ref_pos = self.pos - 1
        score = ""
        
        i = 0
        # each iteration will give us a number-letter pair from the cigar string
        while i < len(self.cigar):
            number = "" # number in the num-letter pair
            while i < len(self.cigar) and self.cigar[i].isdigit():
                number += self.cigar[i]
                i += 1
            number = int(number)
            letter = self.cigar[i] # letter in num-letter pair 
            i += 1
            
            if letter == "M":
                for _ in range(number):
                    if ref_pos == pos:
                        score = self.qual[read_pos]
                    read_pos += 1
                    ref_pos += 1
            elif letter == "I":
                read_pos += number
            elif letter == "D":
                ref_pos += number
            elif letter == "S":
                read_pos += number
        
        return score
        
    
    def mapped_seq(self) -> str:
        mapped = ""
        read_pos = 0
        
        i = 0
        while i < len(self.cigar):
            number = ""
            while i < len(self.cigar) and self.cigar[i].isdigit():
                number += self.cigar[i]
                i += 1
            number = int(number)
            letter = self.cigar[i]
            i += 1
            
            if letter == "M":
                for _ in range(number):
                    mapped += self.seq[read_pos]
                    read_pos += 1
            elif letter == "I":
                for _ in range (number):
                    mapped += self.seq[read_pos]
                read_pos += number
            elif letter == "D":
                for _ in range(number):
                    mapped += "-"
            elif letter == "S":
                read_pos += number
        
        return mapped
                
                
                
            
            
            
            
        
        