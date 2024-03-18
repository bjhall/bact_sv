class SV:
    def __init__(self, sv_type, start_contig, end_contig, start_pos, end_pos, insertion_sequence = None):
        self.sv_type = sv_type
        self.start_contig = start_contig
        self.end_contig = end_contig
        self.start_pos = start_pos
        self.end_pos = end_pos
        self.insertion_sequence = insertion_sequence

    def __repr__(self):
        out = f"{self.sv_type}\t{self.start_contig}:{self.start_pos}"
        if self.start_contig == self.end_contig:
            out += f"-{self.end_pos}"
        else:
            out += f" - {self.start_contig}:{self.end_pos}"

        if self.sv_type in ["INS", "DELINS"]:
            out += f"\tINS:{self.insertion_sequence} LEN:{len(self.insertion_sequence)}"

        return out
