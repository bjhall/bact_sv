import io
import subprocess

class BlastHit:
    def __init__(self, blast_hit_dict):
        self.q_contig = blast_hit_dict['qaccver']
        self.s_contig = blast_hit_dict['saccver']
        self.s_start = int(blast_hit_dict['sstart'])
        self.s_end = int(blast_hit_dict['send'])
        self.q_start = int(blast_hit_dict['qstart'])
        self.q_end = int(blast_hit_dict['qend'])
        self.mismatches = int(blast_hit_dict['mismatch'])
        self.bitscore = float(blast_hit_dict['bitscore'])
        self.cov = float(blast_hit_dict['length']) / float(blast_hit_dict['qlen'])
        self.hit_count = None

    def fwd(self):
        return self.s_start < self.s_end

    def rev(self):
        return self.s_start > self.s_end

    def direction(self):
        return "fwd" if self.fwd() else "rev"

    def __repr__(self):
        string = f"Q: {self.q_contig}:{self.q_start}-{self.q_end} S: {self.s_contig}:{self.s_start}-{self.s_end} COV:{self.cov} MISMATCHES:{self.mismatches}"
        if self.hit_count and self.hit_count > 1:
            string += f" MULTI:{self.hit_count}"
        return string



def run_blast(query_fasta, subject_fasta):

    blast_header_str = "qaccver saccver evalue bitscore mismatch positive pident length qstart qend sstart send slen qlen";
    proc = subprocess.Popen(
        [
            "blastn",
            "-task", "blastn-short",
            "-subject", subject_fasta,
            "-query", query_fasta,
            "-outfmt", f"6 {blast_header_str}",
        ], stdout=subprocess.PIPE
    )

    hits = []
    blast_header_list = blast_header_str.split(" ");
    for line in io.TextIOWrapper(proc.stdout, encoding="utf-8"):
        #print(line)
        hit_parts = line.strip().split("\t")
        hits.append(BlastHit(dict(zip(blast_header_list, hit_parts))))

    return hits
