class SoftClip:
    def __init__(self, contig, direction, position, contig_length, n_clipped_reads, surrounding_depth, n_multi_mapped, n_double_clipped, clip_consensus):
        self.contig = contig
        self.direction = direction
        self.position = int(position)
        self.contig_length = int(contig_length)
        self.n_clipped_reads = int(n_clipped_reads)
        self.surrounding_depth = int(surrounding_depth)
        self.n_multi_mapped = int(n_multi_mapped)
        self.n_double_clipped = int(n_double_clipped)
        self.clip_consensus_sequence = clip_consensus
        self.blast_results = {}
        self.partner = None
        self.sv_type = None

    def __repr__(self):
        string = f"{self.sv_type if self.sv_type else 'remains'} {self.contig}:{self.position} ({self.direction}) CONTIG_LEN: {self.contig_length}, #CLIP_READS:{self.n_clipped_reads} ({self.surrounding_depth}), SEQ:{self.clip_consensus_sequence}, CLIPSEQLEN:{len(self.clip_consensus_sequence)}, N_\
MULTI:{self.n_multi_mapped} N_DOUBLE:{self.n_double_clipped}"
        if self.partner:
            string += f" PARTNER:{self.partner.contig}:{self.partner.position}"
        string += f"\nQUERY_BLAST: {self.blast('query')}"
        string += f"\n  REF_BLAST: {self.blast('reference')}"
        string += f"\n   TN_BLAST: {self.blast('transposons')}"
        return string

    def num_blast_hits(self, label):
        if label not in self.blast_results:
            return 0
        return len(self.blast_results[label])

    def blast(self, label):
        if label not in self.blast_results:
            return None

        hit_to_return = self.blast_results[label][0]
        hit_to_return.hit_count = len(self.blast_results[label])
        return hit_to_return

    def assign_partner(self, clip):
        self.partner = clip

    def set_sv_type(self, sv_type):
        self.sv_type = sv_type

    def is_assigned(self):
        return True if self.sv_type else False
