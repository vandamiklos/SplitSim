import pysam
import argparse
from intervaltree import IntervalTree
from collections import defaultdict
import re
import os


def load_breakpoints_from_fasta(fasta_file):
    """
    Build interval trees from breakpoint junctions encoded in FASTA headers.

    Only internal junctions are included:
    - Excludes the first junction (edge of the read)
    - Excludes the last junction (edge of the read)
    """
    tree_start = defaultdict(IntervalTree)
    tree_end = defaultdict(IntervalTree)
    tree_start2 = defaultdict(IntervalTree)
    tree_end2 = defaultdict(IntervalTree)
    truth_sites = []
    truth_sites2 = []

    pattern = re.compile(r"(chr[^:]+):(\d+)-(\d+)")

    with open(fasta_file) as f:
        for line in f:
            if not line.startswith(">"):
                continue

            header = line.strip()[1:]

            if "__" not in header:
                continue

            _, coords_part = header.split("__", 1)

            # Extract segments: [(chrom, start, end), ...]
            segments = pattern.findall(coords_part)

            # Need at least 2 segments to have an internal junction
            if len(segments) < 2:
                continue

            # Iterate over junctions
            for i in range(0, len(segments) - 1):
                chrom1, s1, e1 = segments[i]
                chrom2, s2, e2 = segments[i + 1]

                pos1 = int(e1)  # end of segment i
                pos2 = int(s2)  # start of segment i+1

                id = str(chrom1) + ":" + str(pos1) + "-" + str(chrom2) + ":" + str(pos2)

                tree_start[chrom1].addi(pos1, pos1 + 1, id)
                tree_end[chrom2].addi(pos2, pos2 + 1, id)

                truth_sites.append(((chrom1, pos1), (chrom2, pos2)))

    return tree_start, tree_end, truth_sites


def find_truth_overlap(tree, chrom, pos, window):
    """
    Return truth breakpoint if overlap exists
    """
    if chrom not in tree:
        return None

    hits = tree[chrom].overlap(pos - window, pos + window)

    if not hits:
        return None

    return [h.data for h in hits]


def detect_cigar_events(read, min_event_size):
    events = []

    if read.is_unmapped or not read.cigartuples:
        return events

    # --- CIGAR-based events ---
    ref_pos = read.reference_start
    chrom = read.reference_name
    mapq = read.mapping_quality

    for op, length in read.cigartuples:
        if op in (0, 7, 8):  # match, mismatch: consumes query and reference
            ref_pos += length
        elif op in (2, 3):  # deletion / skip: only consumes reference
            if length >= min_event_size:
                events.append(("DELETION_SKIP_START", chrom, ref_pos, mapq))
                events.append(("DELETION_SKIP_END", chrom, ref_pos + length, mapq))
            ref_pos += length
        elif op == 1:  # insertion: only consumes query
            if length >= min_event_size:
                events.append(("INSERTION_START", chrom, ref_pos, mapq))
                events.append(("INSERTION_END", chrom, ref_pos + length, mapq))
        elif op == 4:  # soft clip: consumes query
            if length >= min_event_size:
                events.append(("SOFT_CLIP_START", chrom, ref_pos, mapq))
                events.append(("SOFT_CLIP_END", chrom, ref_pos + length, mapq))
    return events


def detect_split_alignment_events(all_alignments, min_event_size):
    """
    Detect split-read events from multiple alignments (LAST --split).

    Parameters
    ----------
    all_alignments : list of pysam.AlignedSegment
        All alignments for the same read.
    min_event_size : int
        Minimum distance to call a breakpoint.

    Returns
    -------
    events : list of tuples
        Each event is (event_type, chrom, pos, mapq)
    """
    events = []

    # Sort alignments by query start
    aln_sorted = sorted(all_alignments, key=lambda x: x.query_alignment_start)

    for i in range(len(aln_sorted) - 1):
        a1 = aln_sorted[i]

        bp1 = a1.reference_end if a1.is_reverse else a1.reference_start
        bp2 = a1.reference_start if a1.is_reverse else a1.reference_end

        events.append(("SPLIT_READ_START",
                        a1.reference_name,
                        bp1,
                        a1.mapping_quality))
        events.append(("SPLIT_READ_END",
                       a1.reference_name,
                       bp2,
                       a1.mapping_quality))

    return events



def tree_search(predicted, tree, window):
    s = set()
    for chrom, pos in predicted:
        hits = find_truth_overlap(tree, chrom, pos, window)
        if hits is not None:
            for x in hits:
                s.add(x)
    return s


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-b", "--bam", required=True)
    parser.add_argument("-t", "--truth", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-n", "--name", required=True)
    parser.add_argument("-w", "--window", default=50, type=int)
    parser.add_argument("-m", "--min-event-size", default=50, type=int)

    args = parser.parse_args()

    # Load truth breakpoints
    tree_start, tree_end, truth_sites = load_breakpoints_from_fasta(args.truth)

    reads_by_qname = defaultdict(list)
    bam = pysam.AlignmentFile(args.bam)
    for read in bam:
        if read.is_unmapped or read.is_secondary:
            continue
        reads_by_qname[read.query_name].append(read)
    bam.close()


    predicted = []

    for read_id, alignments in reads_by_qname.items():

        events = []

        for aln in alignments:
            events.extend(detect_cigar_events(aln, args.min_event_size))

        events.extend(detect_split_alignment_events(alignments, args.min_event_size))

        for (event_type, chrom, pos, mapq) in events:
            predicted.append((chrom, pos))


    # breakpoints predicted from split-reads
    left = tree_search(predicted, tree_start, args.window)
    right = tree_search(predicted, tree_end, args.window)

    # both sides of the break match
    matched = left & right
    # one side of the break matches only
    partly_matched = left ^ right

#    precision = TP / (TP + FP) if TP + FP > 0 else 0
#    recall = TP / (TP + FN) if TP + FN > 0 else 0
#    fscore = 2 * precision * recall / (precision + recall) if precision + recall > 0 else 0

    # --- Write summary ---
    with open(args.output + ".summary.txt", "w") as out:
        out.write(f"tp\ttp_partial\tn_junctions\tn_events\n{len(matched)}\t{len(partly_matched)}\t{len(truth_sites)}\t{len(predicted)}")

    print("Finished")
    print(f"TP: {len(matched)}, TP partial: {len(partly_matched)}")
    print(f"Number of junctions: {len(truth_sites)}, Number of predicted events: {len(predicted)}")
#    print(f"Precision: {precision:.4f}")
#    print(f"Recall: {recall:.4f}")
#    print(f"F-score: {fscore:.4f}")

    path = "/home/vanda/Documents/simulated_split_reads/50bp_gaps/summary_all.txt"
    write_header = not os.path.exists(path)

    with open(path, "a") as out:
        if write_header:
            out.write(f"name\ttp\ttp_partial\tn_junctions\tn_events\n")

        out.write(f"{args.name}\t{len(matched)}\t{len(partly_matched)}\t{len(truth_sites)}\t{len(predicted)}\n")


if __name__ == "__main__":
    main()
