import re

class Read:
    def __init__(self, sam_line: str):
        (qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, *tags) = sam_line.strip().split("\t")

        # store basic properties of the read
        self.qname: str = qname
        self.flag: int = int(flag)
        self.rname: str = rname
        self.pos: str = int(pos)
        self.mapq: str = int(mapq)
        self.cigar: str = cigar
        self.rnext: str = rnext
        self.pnext: str = int(pnext)
        self.tlen: str = int(tlen)
        self.seq: str = seq
        self.qual: str = qual
        self.tags: list[str] = tags

        # score mapping properties based on flag
        self.is_mapped: bool = not bool(self.flag & 4)  # 4 bit not in flag
        self.is_forward: bool = not bool(self.flag & 16)
        self.is_reverse: bool = bool(self.flag & 16)
        self.is_primary: bool = not (bool(self.flag & 256) or bool(self.flag & 2048))  # not secondary or supplemental

        # add data for mapped reads only
        self.cigar_bits: tuple[tuple[int, str]] = None
        self.mapped_len: int = None

        if self.is_mapped:
            self.cigar_bits = tuple([(int(n), cig) for n, cig in re.findall(r"(\d+)([A-Z])", self.cigar)])
            self.mapped_len = sum([n for n, cig in self.cigar_bits if cig in {"M", "D"}])

    def read_idx_at_pos(self, pos: int) -> list[None | int]:
        if not self.is_mapped:
            return []

        # adjust by read start
        pos -= self.pos
        if pos < 0:  # If read mapped to the right of requested location
            return []

        # Check if the requested position is right of our read
        if pos >= self.mapped_len:
            return []

        cigar_bits = re.findall(r"(\d+)([A-Z])", self.cigar)
        # pad position based on cigar
        pad = 0
        mapped_count = 0
        for n, (size, cig_type) in enumerate(cigar_bits):
            size = int(size)
            if cig_type in {"S", "H", "I"}:
                pad += size
                continue
            if mapped_count + size >= pos + 1:
                if cig_type == "M":
                    # pad remaining
                    pad += (pos - mapped_count)
                    break
                if cig_type == "D":
                    return []
            else:
                if cig_type == "M":
                    mapped_count += size
                    pad += size
                if cig_type == "D":
                    mapped_count += size

        # Check if next bases are insertion
        if mapped_count + size == pos + 1:
            if n + 1 != len(cigar_bits):
                size, cig_type = cigar_bits[n + 1]
                size = int(size)
                if cig_type == "I":
                    return [i for i in range(pad, pad + size + 1)]

        return [pad]

    def mapped_seq(self) -> str:
        if not self.is_mapped:
            return ""

        idx = 0  # track where we are in read seq
        bases = []  # list to build up over time without costly string concatenation
        for n, cig in self.cigar_bits:
            if cig == "S":
                idx += n
            elif cig == "D":
                bases += ["-"] * n
            elif cig in {"M", "I"}:
                bases += [self.seq[i] for i in range(idx, idx + n)]
                idx += n

        return "".join(bases)

    def base_at_pos(self, pos: int) -> str:
        idx = self.read_idx_at_pos(pos)
        return "".join([self.seq[i] for i in idx])

    def qual_at_pos(self, pos: int) -> str:
        idx = self.read_idx_at_pos(pos)
        return "".join([self.qual[i] for i in idx])

    def cigar_length(self) -> int:
        """Calculates the length of the read's span on the reference sequence"""
        if self.cigar == "*":  # No alignment
            return 0

        ref_length = 0
        i = 0
        while i < len(self.cigar):
            number = ""
            while i < len(self.cigar) and self.cigar[i].isdigit():
                number += self.cigar[i]
                i += 1
            number = int(number)
            letter = self.cigar[i]
            i += 1

            # Only count operations that consume reference positions
            if letter in "MD":
                ref_length += number

        return ref_length


class SAM:
    def __init__(self):
        self.reads = []

    @classmethod
    def from_sam(cls, sam_file_path):
        sam = SAM()

        sam_array = []
        with open(sam_file_path, 'r') as sam_file:
            for line in sam_file:
                line = line.strip()
                if line[0] == '@':
                    continue
                sam_array.append(Read(line))

        for read in sam_array:
            if read.is_mapped and read.is_primary:
                sam.reads.append(read)

        return sam

    def reads_at_pos(self, seq_name, pos) -> list:
        '''Returns a list of all Read instances that map to that position in the reference'''
        read_instances = []

		# go through all the reads
        for read in self.reads:
            if seq_name != read.rname:
                continue

            base = read.base_at_pos(pos)
            if base:
                read_instances.append(read)

        return read_instances

    def pileup_at_pos(self, seq_name, pos) -> tuple[list[str], list[str]]:
        '''Returns tuple of a list of the base calls and a list of quality scores'''
        base_calls = self.reads_at_pos(seq_name, pos)
        quality_scores = []

        for read in self.reads:
            if seq_name != read.rname:
                continue

            base = read.base_at_pos(pos)
            if base:
                base_calls.append(base)
                quality_scores.append(read.qual_at_pos(pos))

        return (base_calls, quality_scores)

    def consensus_at_pos(self, seq_name, pos) -> str:
        '''Given a sequence name, return the majority base call at a specific position'''
        base_calls, _ = self.pileup_at_pos(seq_name, pos)

        # count the occurrences of each base through a dictionary
        base_counts = {"A": 0, "T": 0, "C": 0, "G": 0, "DEL": 0}
        for base in base_calls:
            if base in base_counts:
                base_counts[base] += 1

        # determine the majority base
        total_bases = sum(base_counts.values())
        # base = key, count = value
        for base, count in base_counts.items():
            if count > 0.5 * total_bases:
                return base

        # no majority base
        return "N"

    def consensus(self, seq_name) -> str:
        '''Given a sequence name, return majority bases for ALL positions'''
        # get all reads that map to the sequence
        reads_for_sequence = []
        for read in self.reads:
            if read.rname == seq_name:
                reads_for_sequence.append(read)

        # if no reads map to sequence, return an empty string
        if not reads_for_sequence:
            return ""

        # find the min and max positions of reads mapped to the reference
        min_pos = min(read.pos for read in reads_for_sequence)
        max_pos = max(read.pos + read.cigar_length() - 1 for read in reads_for_sequence)

        # get the majority base for each position in the reference range
        consensus = ""
        for position in range(min_pos, max_pos + 1):
            consensus += self.consensus_at_pos(seq_name, position)

        return consensus

    def best_consensus(self) -> str:
        # gets t
        coverage = {}
        for read in self.reads:
            if read.rname not in coverage:
                coverage[read.rname] = 0
            coverage[read.rname] += read.cigar_length()

        # find the sequence with the highest coverage
        best_mapped = max(coverage, key=coverage.get)

        # generate the consensus for the best-mapped sequence
        best_consensus = self.consensus(best_mapped)

        return best_consensus
