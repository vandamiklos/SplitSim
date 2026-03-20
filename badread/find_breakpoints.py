import pysam
import argparse
from intervaltree import IntervalTree
from collections import defaultdict
import re


def load_breakpoints_from_fasta(fasta_file, window):
    """
    Build interval trees from breakpoint junctions encoded in FASTA headers,
    excluding the start of the first segment and the end of the last segment.
    """
    trees = defaultdict(IntervalTree)
    truth_sites = defaultdict(set)

    pattern = re.compile(r"(chr[^:]+):(\d+)-(\d+)")

    with open(fasta_file) as f:
        for line in f:
            if not line.startswith(">"):
                continue

            header = line.strip()[1:]

            if "__" not in header:
                continue

            _, coords_part = header.split("__", 1)

            # Extract all segments
            segments = pattern.findall(coords_part)

            # Skip reads with only one segment (no junctions)
            if len(segments) < 2:
                continue

            # Build junctions between consecutive segments
            for i in range(len(segments) - 1):
                chrom1, s1, e1 = segments[i]
                chrom2, s2, e2 = segments[i + 1]

                # Convert to integers
                pos1 = int(e1)  # end of segment 1
                pos2 = int(s2)  # start of segment 2

                # Skip the very first segment's start and last segment's end
                # i > 0 ensures we skip the first segment's end (which would be the start)
                # i < len(segments) - 2 ensures we skip last segment's start
                for j, (chrom, pos) in enumerate([(chrom1, pos1), (chrom2, pos2)]):
                    is_first_segment = i == 0 and j == 0
                    is_last_segment = i == len(segments) - 2 and j == 1
                    if is_first_segment or is_last_segment:
                        continue

                    if pos not in truth_sites[chrom]:
                        trees[chrom].addi(pos - window, pos + window, pos)
                        truth_sites[chrom].add(pos)

    return trees, truth_sites


def find_truth_overlap(trees, chrom, pos, window):
    """
    Return truth breakpoint if overlap exists
    """
    if chrom not in trees:
        return None

    hits = trees[chrom].overlap(pos - window, pos + window)

    if not hits:
        return None

    # return closest truth breakpoint
    return min(hits, key=lambda x: abs(x.data - pos)).data


def detect_events(read, all_alignments=None, min_event_size=50):
    """
    Detect structural variant events from a read.

    Parameters
    ----------
    read : pysam.AlignedSegment
        The read to analyze.
    all_alignments : list of pysam.AlignedSegment, optional
        All alignments for the same read (for split-read detection in LAST).
    min_event_size : int
        Minimum size of indels/soft clips to count as events.

    Returns
    -------
    events : list of tuples
        Each event is (event_type, chrom, pos, mapq)
    """
    events = []

    if read.is_unmapped or not read.cigartuples:
        return events

    # --- CIGAR-based events ---
    ref_pos = read.reference_start
    chrom = read.reference_name
    mapq = read.mapping_quality

    for op, length in read.cigartuples:
        if op == 0:  # match
            ref_pos += length
        elif op in (2, 3):  # deletion / skip
            if length >= min_event_size:
                events.append(("DELETION_SKIP", chrom, ref_pos, mapq))
            ref_pos += length
        elif op == 1:  # insertion
            if length >= min_event_size:
                events.append(("INSERTION", chrom, ref_pos, mapq))
        elif op == 4:  # soft clip
            if length >= min_event_size:
                events.append(("SOFT_CLIP", chrom, ref_pos, mapq))

    # --- Split-read events via SA tag (BWA/Minimap2) ---
    if read.has_tag("SA"):
        sa_entries = read.get_tag("SA").split(";")
        for entry in sa_entries:
            if not entry:
                continue
            chrom2, pos, strand, cigar, mq, nm = entry.split(",")
            events.append(("SPLIT_READ", chrom2, int(pos), int(mq)))

    # --- Split-read events via multiple alignments (LAST --split) ---
    if all_alignments is not None:
        # Sort alignments by reference start
        aln_sorted = sorted(all_alignments, key=lambda x: x.reference_start)
        for i in range(len(aln_sorted) - 1):
            a1 = aln_sorted[i]
            a2 = aln_sorted[i + 1]

            # Only consider alignments on different chromosomes or far apart
            if (a1.reference_name != a2.reference_name or
                    abs(a2.reference_start - (a1.reference_start + a1.query_alignment_length)) > min_event_size):
                # Use the junction at the start of the second alignment
                events.append(("SPLIT_READ", a2.reference_name, a2.reference_start, a2.mapping_quality))

    return events

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-b", "--bam", required=True)
    parser.add_argument("-t", "--truth", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-w", "--window", default=50, type=int)

    args = parser.parse_args()

    # Load truth breakpoints
    trees, truth_sites = load_breakpoints_from_fasta(args.truth, args.window)

    # Store results
    results = []
    matched_truth = set()
    all_read_events = []

    bam = pysam.AlignmentFile(args.bam)

    # Collect all read events
    reads_by_qname = {}
    for read in bam:
        events = detect_events(read)
        all_read_events.extend([(read.query_name, *e) for e in events])
        # For LAST-style split reads, you could pass all alignments here if needed
    bam.close()

    fp_sites = set()

    # --- Compare each read event to truth ---
    for read_id, event_type, chrom, pos, mapq in all_read_events:

        truth_hit = find_truth_overlap(trees, chrom, pos, args.window)

        if truth_hit is not None:
            results.append((read_id, event_type, chrom, pos, mapq, truth_hit, "TP"))
            matched_truth.add((chrom, truth_hit))
        else:
            results.append((read_id, event_type, chrom, pos, mapq, "NA", "FP"))
            fp_sites.add((chrom, pos))

    # --- Compute FN strictly from truth sites ---
    all_truth = {(chrom, pos) for chrom in truth_sites for pos in truth_sites[chrom]}
    fn_sites = all_truth - matched_truth

    TP = len(matched_truth)
    FP = len(fp_sites)
    FN = len(fn_sites)

    precision = TP / (TP + FP) if TP + FP > 0 else 0
    recall = TP / (TP + FN) if TP + FN > 0 else 0
    fscore = 2 * precision * recall / (precision + recall) if precision + recall > 0 else 0

    # --- Write event results ---
    with open(args.output + ".events.csv", "w") as out:
        out.write("read_id\tevent_type\tchrom\tpos\tmapq\ttruth_pos\tclass\n")
        for r in results:
            out.write("\t".join(map(str, r)) + "\n")

    # --- Write FN list ---
    with open(args.output + ".fn.csv", "w") as out:
        out.write("chrom\ttruth_pos\n")
        for chrom, pos in sorted(fn_sites):
            out.write(f"{chrom}\t{pos}\n")

    # --- Write summary ---
    with open(args.output + ".summary.txt", "w") as out:
        out.write(f"TP\tFP\tFN\tPrecision\tRecall\tF-score\n{TP}\t{FP}\t{FN}\t{precision:.4f}\t{recall:.4f}\t{fscore:.4f}")

    print("Finished")
    print(f"TP: {TP}")
    print(f"FP: {FP}")
    print(f"FN: {FN}")
    print(f"Precision: {precision:.4f}")
    print(f"Recall: {recall:.4f}")
    print(f"F-score: {fscore:.4f}")


if __name__ == "__main__":
    main()
