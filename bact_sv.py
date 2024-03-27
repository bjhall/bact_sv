#!/usr/bin/env python

# TODOs/FIXMEs:
#  * If both sequences just next to deletions are identical, the "distance" between breakpoints will be high and might fail. Check for these repeated sequences
#  * Repeat deletions (ie nice deletion clips, but the deleted sequence is still elsewhere in the assembly (example /storage/bjorn/hdf3/seq3651_ub651.bam 52832-54027)
#  * Segment shuffling... (/storage/bjorn/b5699/seq3634_ub597.bam, CP060061.1:4745763-4747872)
#  * More

#  MISSED SVs:
#    * Insertion at 4378689 in  /storage/bjorn/analyses/emre_b52_insertions_231109/seq3721_b5220.to.CP024090.bam

import os
import sys
import subprocess
import tempfile
import argparse

from collections import defaultdict
from pprint import pprint as pp
from sv import SV
from blast import BlastHit, run_blast
from softclip import SoftClip

script_dir = os.path.dirname(os.path.realpath(__file__))

MAX_DELETION_SIZE = 500_000



def main(arg):
    clips = get_softclip_consensus_sequences(arg.bam)

    print("BLASTing consensus clips to query assembly genome...", file=sys.stderr)
    blast_all_clips_to_assembly(clips, arg.assembly, "query")

    print("BLASTing consensus clips to reference assembly genome...", file=sys.stderr)
    blast_all_clips_to_assembly(clips, arg.reference, "reference")

    print("BLASTing consensus clips to transposon references...", file=sys.stderr)
    blast_all_clips_to_assembly(clips, f"{script_dir}/transposons.fna", "transposons", max_mismatches = 5)

    query_assembly_sequences = read_fasta(arg.assembly)
    calls = call_svs(clips, query_assembly_sequences)

    for clip in sorted(clips, key=lambda d: d.position):
        if not clip.is_assigned():
            print(clip)

    print()
    for call in calls:
        print("CALL\t",call)
    #pp(calls)


def call_svs(clips, query_assembly_sequences):
    all_calls = []

    # FIXME: This excludes all clips anchored in a repeat. These could probably be used to
    #        strengthen the evidence for Tn inserts.
    flag_clips_anchored_in_repeats(clips)


    # Flag clips with reads that were clipped in both ends. These are typically supposed to map to
    # something that's not present in the reference assembly (insertion, new/unknown plasmid etc)
    flag_double_clipped(clips)

    flag_low_freq_clips(clips)

    # Call different types of SVs
    all_calls.extend(call_circ(clips))
    all_calls.extend(call_inv(clips))
    all_calls.extend(call_del(clips))
    all_calls.extend(call_ins(clips, query_assembly_sequences))
    all_calls.extend(call_delins(clips, query_assembly_sequences))
    all_calls.extend(call_tn_ins(clips))


    # FIXME: This excludes weird varible clips that couldn't be assigned.
    #        I currently don't understand what these are! Technical artifacts?
    flag_variable_short_clips(clips)

    # TODO: Flag various orphans?

    return all_calls


################################
def flag_clips_anchored_in_repeats(clips):
    for clip in clips:
        if not clip.is_assigned() and clip.n_multi_mapped / clip.n_clipped_reads > 0.9:
            clip.set_sv_type("REPEAT_ANCHORED")


def flag_double_clipped(clips):
    for clip in clips:
        if not clip.is_assigned() and clip.n_double_clipped / clip.n_clipped_reads > 0.5:
            clip.set_sv_type("DOUBLE_CLIPPED")


def flag_variable_short_clips(clips):
    for clip in clips:
        n_Ns = clip.clip_consensus_sequence.count("N")
        clip_len = len(clip.clip_consensus_sequence)
        if n_Ns >= 1 and not clip.blast("reference") and not clip.blast("query"):
            clip.set_sv_type("WEIRD")

def flag_low_freq_clips(clips, cutoff = 0.05):
    for clip in clips:
        rough_clip_frequency_estimate = clip.n_clipped_reads / clip.surrounding_depth
        if rough_clip_frequency_estimate < cutoff:
            clip.set_sv_type("LOW_FREQ")


def call_delins(clips, query_assembly_sequences):
    calls = []
    for idx1 in range(0,len(clips)):
        selected_R, selected_L = None, None
        smallest_distance = None
        for idx2 in range(idx1+1,len(clips)):
            clip1 = clips[idx1]
            clip2 = clips[idx2]

            # Skip clips that are already assigned
            if clip1.is_assigned() or clip2.is_assigned():
                continue

            # Skip if clips are on different contigs
            if clip1.contig != clip2.contig:
                continue

            # Skip if clips are in same direction
            if clip1.direction == clip2.direction:
                continue

            # Skip if either clip is mapping to the reference genome
            if clip1.blast("reference") or clip2.blast("reference"):
                continue

            # Skip if either clip is mapping to the reference genome
            if not clip1.blast("query") or not clip2.blast("query"):
                continue

            R, L = order_RL_clips(clip1, clip2)

            # Skip if left clip comes before right clip
            if R.position > L.position:
                continue

            # Find closest deletion clip partner
            distance = L.position - R.position
            if not smallest_distance or distance < smallest_distance:
                selected_R = R
                selected_L = L
                smallest_distance = distance

        if selected_R:
            Rqry = selected_R.blast("query")
            Lqry = selected_L.blast("query")

            # Can we find the insertion sequence?
            insertion_sequence = "Unknown insertion sequence"
            if Rqry.s_contig == Lqry.s_contig and Rqry.direction() == Lqry.direction():
                if Rqry.fwd():
                    insert_sequence_contig = Rqry.s_contig
                    insert_sequence_start = Rqry.s_start - 1
                    insert_sequence_end = Lqry.s_end
                elif Rqry.rev():
                    insert_sequence_contig = Rqry.s_contig
                    insert_sequence_start = Lqry.s_end - 1
                    insert_sequence_end = Rqry.s_start

                insertion_sequence = query_assembly_sequences[insert_sequence_contig][insert_sequence_start:insert_sequence_end]
                if Rqry.direction() == "rev":
                    insertion_sequence = revcomp(insertion_sequence)

            assign_sv_to_clips(R=selected_R, L=selected_L, sv_type="DELINS")
            calls.append(
                SV("DELINS", selected_R.contig, selected_L.contig, selected_R.position+1, selected_L.position, insertion_sequence=insertion_sequence)
            )

    return calls



def call_inv(clips):
    calls = []

    inversion_pairs = []
    for idx1 in range(0,len(clips)):
        for idx2 in range(idx1+1,len(clips)):

            clip1 = clips[idx1]
            clip2 = clips[idx2]

            # Skip clips that are already assigned
            if clip1.is_assigned() or clip2.is_assigned():
                continue

            # Skip if clips are in same direction
            if clip1.direction == clip2.direction:
                continue

            # Skip unless both map to reference
            if not clip1.blast("reference") or not clip2.blast("reference"):
                continue

            R, L = order_RL_clips(clip1, clip2)

            Rref = R.blast("reference")
            Lref = L.blast("reference")

            # Skip if clips mapped to different contigs of the reference
            if Rref.s_contig != Lref.s_contig:
                continue

            if Rref.rev() and Lref.rev():
                if dist(Rref.s_start, L.position) < 20 and dist(Lref.s_end, R.position) < 20:
                    inversion_pairs.append({'R':R, 'L':L})


    for idx1 in range(0, len(inversion_pairs)):
        for idx2 in range(idx1+1, len(inversion_pairs)):
            pair1 = inversion_pairs[idx1]
            pair2 = inversion_pairs[idx2]

            if pair1['R'].contig != pair2['R'].contig:
                continue

            if dist(pair1['L'].position, pair2['R'].position) < 20 and dist(pair1['R'].position, pair2['L'].position) < 20:
                assign_sv_to_clips(R=pair1['R'], L=pair1['L'], sv_type="INV")
                assign_sv_to_clips(R=pair2['R'], L=pair2['L'], sv_type="INV")
                sv_start = min(pair1['L'].position, pair2['R'].position, pair2['L'].position, pair1['R'].position)
                sv_end = max(pair1['L'].position, pair2['R'].position, pair2['L'].position, pair1['R'].position)
                contig = pair1['R'].contig
                calls.append(
                    SV("INV", contig, contig, sv_start, sv_end)
                )

    return calls



def call_tn_ins(clips):
    calls = []

    for clip1, clip2 in triang_iter(clips):

        # Skip clips that are already assigned
        if clip1.is_assigned() or clip2.is_assigned():
            continue

        # Skip if clips are in same direction
        if clip1.direction == clip2.direction:
            continue

        # Requite that at least one clip matches a known transposon or has multiple identical matches to the reference
        transposon_match = clip1.blast("transposons") or clip2.blast("transposons")
        multi_match = clip1.num_blast_hits("reference") > 1 or clip2.num_blast_hits("reference") > 1
        if not transposon_match and not multi_match:
            continue

        R, L = order_RL_clips(clip1, clip2)

        if dist(R.position, L.position) > 20:
            continue

        assign_sv_to_clips(R=R, L=L, sv_type="TN_INS")
        calls.append(
            SV("TN_INS",  R.contig, L.contig, R.position, L.position)
        )

    return calls



def call_ins(clips, query_assembly_sequences):
    calls = []

    for clip1, clip2 in triang_iter(clips):

        # Skip clips that are already assigned
        if clip1.is_assigned() or clip2.is_assigned():
            continue

        # Skip if clips are in same direction
        if clip1.direction == clip2.direction:
            continue

        # Skip if either clip is mapping to the reference genome
        if clip1.blast("reference") or clip2.blast("reference"):
            continue

        # Skip unless both clips map to the query assembly
        if not clip1.blast("query") or not clip2.blast("query"):
            continue

        R, L = order_RL_clips(clip1, clip2)

        if dist(R.position, L.position) > 10:
            continue

        Rqry = R.blast("query")
        Lqry = L.blast("query")

        # Can we find the insertion sequence?
        insertion_sequence = "Unknown insertion sequence"
        if Rqry.s_contig == Lqry.s_contig and Rqry.direction() == Lqry.direction():
            if Rqry.direction() == "rev":
                insert_sequence_contig = Rqry.s_contig
                insert_sequence_start = Lqry.s_end - 1
                insert_sequence_end = Rqry.s_start
            else:
                insert_sequence_contig = Rqry.s_contig
                insert_sequence_start = Rqry.s_start - 1
                insert_sequence_end = Lqry.s_end

            insertion_sequence = query_assembly_sequences[insert_sequence_contig][insert_sequence_start:insert_sequence_end]
            if Rqry.direction() == "rev":
                insertion_sequence = revcomp(insertion_sequence)

        if insertion_sequence == "Unknown insertion sequence":
            pp(Rqry)
            pp(Lqry)

        assign_sv_to_clips(R=R, L=L, sv_type="INS")
        calls.append(
            SV("INS", R.contig, L.contig, R.position, L.position, insertion_sequence=insertion_sequence)
        )

    return calls

def call_circ(clips):
    calls = []

    for clip1, clip2 in triang_iter(clips):

        # Skip clips that are already assigned
        if clip1.is_assigned() or clip2.is_assigned():
            continue

        # Skip clips that are not on the same contig
        if clip1.contig != clip2.contig:
            continue

        # Skip if either of the clips didn't match the reference assembly
        if not clip1.blast("reference") or not clip2.blast("reference"):
            continue

        # Skip unless one is a left clip and one is a right clip
        if clip1.direction == clip2.direction:
            continue

        R, L = order_RL_clips(clip1, clip2)

        Rref = R.blast("reference")
        Lref = L.blast("reference")

        # Skip if clips map to different contigs in the reference
        if Rref.s_contig != Lref.s_contig:
            continue

        # Perfect circularization
        if L.position == 0 and R.position == R.contig_length:
            if Rref.s_start == 1 and Lref.s_end == L.contig_length:
                assign_sv_to_clips(R=R, L=L, sv_type="CIRC")
                calls.append(
                    SV("CIRC", R.contig, L.contig, R.position, L.position)
                )
            else:
                print("WE ShOUDLNT BE HERE")

    return calls

################################
def call_del(clips):
    calls = []
    for clip1, clip2 in triang_iter(clips):

        # Skip clips that are already assigned
        if clip1.is_assigned() or clip2.is_assigned():
            continue

        # Skip if either of the clips didn't match the reference assembly
        if not clip1.blast("reference") or not clip2.blast("reference"):
            continue

        # Skip unless one is a left clip and one is a right clip
        if clip1.direction == clip2.direction:
            continue

        R, L = order_RL_clips(clip1, clip2)

        Rref = R.blast("reference")
        Lref = L.blast("reference")

        if Rref.s_contig == Lref.s_contig:

            # Deletion
            if Rref.fwd() and Lref.fwd() and R.position < L.position and dist(Rref.s_start, L.position) < 30 and dist(Lref.s_end, R.position) < 30:
                assign_sv_to_clips(R=R, L=L, sv_type="DEL")
                calls.append(
                    SV("DEL", R.contig, L.contig, R.position+1, L.position)
                )

            elif Rref.fwd() and Lref.fwd() and R.position > L.position and dist(Rref.s_start, L.position) < 30 and dist(Lref.s_end, R.position) < 30:
                # Deletion spanning circularization point
                # FIXME: This requires a higher distance because one such deletion had identical repeats next to the deletions causing the distance to be
                #        inflated. Fix this by checking for these repeats and take them into account when calculating distances.
                #        Example b7396_a7bcac38_1:4873250-29734 in /storage/bjorn/hdf3/seq3651_ub651.bam
                distance = (R.contig_length - R.position) + L.position
                if distance <= MAX_DELETION_SIZE:
                    assign_sv_to_clips(R=R, L=L, sv_type="DEL")
                    calls.append(
                        SV("CIRC_DEL", R.contig, L.contig, R.position+1, L.position)
                    )
    return calls


def dist(a, b):
    return abs(int(a)-int(b))

##################################
def blast_all_clips_to_assembly(clips, subject_fasta, label, max_mismatches = 1, min_coverage = 0.9):
    with tempfile.NamedTemporaryFile(delete=False) as temp_file:
        for clip_idx, clip in enumerate(clips):
            temp_file.write(bytes(f">{clip_idx}\n{clip.clip_consensus_sequence}\n", 'UTF-8'))

    hits = run_blast(temp_file.name, subject_fasta)
    for hit in hits:

        # Filter out hits that don't fulfill the identity/coverage criteria
        if hit.mismatches > max_mismatches or hit.cov < min_coverage:
            continue

        query_idx = int(hit.q_contig)
        if not clips[query_idx].blast(label) or hit.bitscore > clips[query_idx].blast_results[label][0].bitscore:
            clips[query_idx].blast_results[label] = [hit]
        elif hit.bitscore == clips[query_idx].blast_results[label][0].bitscore:
            clips[query_idx].blast_results[label].append(hit)



##################################
def get_softclip_consensus_sequences(bam_path):
    ''' Collection consensus softclipped sequences from bam file'''
    clips = []
    script_dir = os.path.dirname(os.path.realpath(__file__))
    softclips_out = subprocess.check_output([f"{script_dir}/target/release/bact_sv", bam_path], text=True).split("\n")
    for line in softclips_out:
        if len(line) == 0:
            continue
        direction, contig, position, contig_length, n_clipped_reads, surrounding_depth, n_multi_mapped, n_double_clipped, clip_consensus = line.strip().split("\t")
        clips.append(SoftClip(contig, direction, position, contig_length, n_clipped_reads, surrounding_depth, n_multi_mapped, n_double_clipped, clip_consensus))

    return clips


############### UTILS ####################
def assign_sv_to_clips(R, L, sv_type):
    R.assign_partner(L)
    L.assign_partner(R)
    R.set_sv_type(sv_type)
    L.set_sv_type(sv_type)


def read_fasta(fasta_path):
    seq = defaultdict(str)
    seq_id = None
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                seq_id = line[1:].strip().split(" ")[0]
            else:
                seq[seq_id] += line.strip()
    return seq


def triang_iter(clips):
    for idx1 in range(0,len(clips)):
        for idx2 in range(idx1+1,len(clips)):
            yield clips[idx1], clips[idx2]


def order_RL_clips(clip1, clip2):
    if clip1.direction == "RIGHT":
        return clip1, clip2
    else:
        return clip2, clip1


def revcomp(dna_seq):
    return dna_seq[::-1].upper().replace("A", "t").replace("T", "a").replace("C", "g").replace("G", "c").upper()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='bact_sv', description='Find structural variants in bacterial genomes')
    parser.add_argument('bam')
    parser.add_argument('assembly', help="Assembly of query genome", type=str)
    parser.add_argument('reference', help="Assembly of reference genome", type=str)
    args = parser.parse_args()
    main(args)
