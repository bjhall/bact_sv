//use std::str;
use std::rc::Rc;
use rust_htslib::{bam, bam::Read, bam::record::Aux};
use std::collections::HashMap;
use std::env;

#[derive(Debug)]
struct SoftClip {
    reference_contig: usize,
    reference_start_position: usize,
    clipped_sequence: Vec<u8>,
    is_left_clip: bool,
    multi_mapper: bool,
    double_clip: bool,
}

#[derive(Eq, Hash, PartialEq)]
struct Location {
    contig_name: String,
    position: usize,
    contig_size: i32,
}

fn main() {
    let args: Vec<String> = env::args().collect();

    let mut bam = bam::Reader::from_path(&args[1]).unwrap();

    let contig_sizes = get_contig_sizes_from_bam_header(&bam);
    let contig_names = get_contig_names_from_bam_header(&bam);

    eprintln!("Collecting soft clips...");
    let (softclips, clip_count_by_reference_position) = collect_all_softclips(&mut bam);

    let mut left_clips_by_position: HashMap<Rc<Location>, Vec<Vec<u8>>> = HashMap::new();
    let mut right_clips_by_position: HashMap<Rc<Location>, Vec<Vec<u8>>> = HashMap::new();
    let mut left_clips_multi_count: HashMap<Rc<Location>, usize> = HashMap::new();
    let mut right_clips_multi_count: HashMap<Rc<Location>, usize> = HashMap::new();
    let mut left_clips_double_clip_count: HashMap<Rc<Location>, usize> = HashMap::new();
    let mut right_clips_double_clip_count: HashMap<Rc<Location>, usize> = HashMap::new();


    // Group clips by reference position
    for (_i, softclip) in softclips.iter().enumerate() {

        let location = Rc::new(Location{
            contig_name: contig_names.get(&softclip.reference_contig).unwrap().to_string(),
            position: softclip.reference_start_position,
            contig_size: *contig_sizes.get(&softclip.reference_contig).unwrap(),
        });

        let n_clips_in_position = *clip_count_by_reference_position.get(&softclip.reference_start_position).unwrap();

        // Skip position where we don't have enough soft clipped reads
        if n_clips_in_position < 5 { continue }

        if softclip.is_left_clip {
            let entry = left_clips_by_position.entry(Rc::clone(&location)).or_insert(Vec::new());
            entry.push(softclip.clipped_sequence.clone());
            if softclip.multi_mapper {
                *left_clips_multi_count.entry(Rc::clone(&location)).or_insert(0) += 1;
            }
            if softclip.double_clip {
                *left_clips_double_clip_count.entry(Rc::clone(&location)).or_insert(0) += 1;
            }
        }
        else {
            let entry = right_clips_by_position.entry(Rc::clone(&location)).or_insert(Vec::new());
            entry.push(softclip.clipped_sequence.clone());
            if softclip.multi_mapper {
                *right_clips_multi_count.entry(Rc::clone(&location)).or_insert(0) += 1;
            }
            if softclip.double_clip {
                *right_clips_double_clip_count.entry(Rc::clone(&location)).or_insert(0) += 1;
            }
        }
    }

    for (location, clip_sequences) in &right_clips_by_position {
        let n_multi_mapped = right_clips_multi_count.get(location).unwrap_or(&0);
        let n_double_clipped = right_clips_double_clip_count.get(location).unwrap_or(&0);
        let consensus = generate_consensus_clip(clip_sequences, false);
        if consensus.is_none() { continue };
        println!("RIGHT\t{}\t{}\t{}\t{}\t{}\t{}\t{}", location.contig_name, location.position, location.contig_size, clip_sequences.len(), n_multi_mapped, n_double_clipped, String::from_utf8_lossy(&consensus.unwrap()));
    }

    for (location, clip_sequences) in &left_clips_by_position {
        let n_multi_mapped = left_clips_multi_count.get(location).unwrap_or(&0);
        let n_double_clipped = left_clips_double_clip_count.get(location).unwrap_or(&0);
        let consensus = generate_consensus_clip(clip_sequences, true);
        if consensus.is_none() { continue };
        println!("LEFT\t{}\t{}\t{}\t{}\t{}\t{}\t{}", location.contig_name, location.position, location.contig_size, clip_sequences.len(), n_multi_mapped, n_double_clipped, String::from_utf8_lossy(&consensus.unwrap()));
    }


}
//35M1D98M18S


fn get_contig_sizes_from_bam_header(bam: &bam::Reader) -> HashMap<usize, i32> {
    let header = bam::Header::from_template(bam.header());

    let mut contig_sizes = HashMap::new();
    for (key, records) in header.to_hashmap() {
        if key == "SQ" {
            for (idx, record) in records.iter().enumerate() {
                contig_sizes.insert(idx, record["LN"].parse().unwrap());
            }
        }
    }
    return contig_sizes

}

fn get_contig_names_from_bam_header(bam: &bam::Reader) -> HashMap<usize, String> {
    let header = bam::Header::from_template(bam.header());

    let mut contig_names = HashMap::new();
    for (key, records) in header.to_hashmap() {
        if key == "SQ" {
            for (idx, record) in records.iter().enumerate() {
                contig_names.insert(idx, record["SN"].to_string());
            }
        }
    }
    return contig_names
}


fn collect_all_softclips(bam: &mut bam::Reader) -> (Vec<SoftClip>, HashMap<usize, usize>) {
    let mut softclips: Vec<SoftClip> = Vec::new();
    let mut clip_count_by_reference_position: HashMap<usize, usize> = HashMap::new();

    let mut skipped_due_to_missing_nm = 0;
    let mut skipped_due_to_too_many_mismatches = 0;

    for r in bam.records() {
        let record = r.unwrap();

        // Skip unmapped reads
        if record.pos() < 0 {
            continue
        }

        // Skip reads with too many mismatches
        let mut mismatches = 0;
        match record.aux(b"NM") {
            Ok(value) => {
                if let Aux::U8(v) = value {
                    mismatches = v;
                }
            }
            Err(..) => {
                skipped_due_to_missing_nm += 1;
                continue;
            }
        }
        if mismatches > 3 {
            skipped_due_to_too_many_mismatches += 1;
            continue
        }


        let cigar = record.cigar();
        let mut reference_position: usize = record.pos().try_into().unwrap();
        let mut read_position: usize = 0;

        let mut n_softclips = 0;
        for cigar_op in &cigar{
            let cigar_op_length: usize = cigar_op.len().try_into().unwrap();

            if cigar_op.char() == 'S' {
                if cigar_op_length < 7 {
                    continue
                }
                let seq = record.seq().as_bytes();
                let clipped_sequence = &seq[read_position..read_position + cigar_op_length];
                let is_left_clip = read_position == 0;

                *clip_count_by_reference_position.entry(reference_position).or_insert(0) += 1;
                let softclip = SoftClip {
                    reference_contig: record.tid() as usize,
                    reference_start_position: reference_position,
                    clipped_sequence: clipped_sequence.to_vec(),
                    is_left_clip: is_left_clip,
                    multi_mapper: record.mapq() < 5,
                    double_clip: false,
                };
                softclips.push(softclip);
                n_softclips += 1;
            }

            match cigar_op.char() {
                'M' | '=' | 'X' => {
                    read_position += cigar_op_length;
                    reference_position += cigar_op_length;
                },
                'D' | 'N' => reference_position += cigar_op_length,
                'I' | 'S' => read_position += cigar_op_length,
                _ => (),
            }
        }

        // Flag clips from reads that were clipped in both ends
        if n_softclips > 1 {
            let len = softclips.len();
            softclips[len-1].double_clip = true;
            softclips[len-2].double_clip = true;
        }
    }
    return (softclips, clip_count_by_reference_position)
}




fn generate_consensus_clip(softclip_sequences: &[Vec<u8>], reverse: bool) -> Option<Vec<u8>> {

    let max_len = softclip_sequences.iter().map(|seq| seq.len()).max().unwrap();
    let mut consensus = vec![b'N'; max_len];
    let mut n_bad_sites = 0;

    for i in 0..max_len {
        let mut counts: HashMap<u8, usize> = HashMap::new();
        let mut n_reads_in_position = 0;

        for seq in softclip_sequences {
            if i < seq.len() {
                let pos = if reverse { seq.len() - 1 - i } else { i };
                let base = seq[pos];
                *counts.entry(base).or_insert(0) += 1;
                n_reads_in_position += 1;
            }
        }

        // Stop traversal if only one read remaining ...
        if n_reads_in_position < 2 { break }

        // ... or if two reads remain and they disagree on the current base (this is done to avoid lots of Ns when there are indels in one of the two clips)
        // FIXME: Maybe the clips should be aligned to avoid these problems and allow further extension of the clip consensus
        if n_reads_in_position == 2 && counts.keys().len() > 1 { break}


        if let Some((most_common_base, max_count)) = counts.into_iter().max_by_key(|&(_, count)| count) {
            if max_count as f32 >= 0.65 * n_reads_in_position as f32 {
                consensus[i] = most_common_base;
            } else {
                n_bad_sites += 1;
                consensus[i] = b'N'; // Insert 'N' if less than 75% of reads agree
            }
        }
    }

    // Remove trailing Ns
    while let Some(&b'N') = consensus.last() {
        consensus.pop();
    }

    if n_bad_sites > 4 || consensus.len() < 10 {
        return None
    }

    if reverse {
        return Some(consensus.into_iter().rev().collect())
    }
    return Some(consensus)
}
