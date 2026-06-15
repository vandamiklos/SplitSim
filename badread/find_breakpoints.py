import pysam
import argparse
from intervaltree import IntervalTree
from collections import defaultdict
import re

from pip._internal.resolution.resolvelib import candidates


def load_breakpoints_from_fasta(fasta_file):
    """
    Build interval trees from breakpoint junctions encoded in FASTA headers.

    Only internal junctions are included:
    - Excludes the first junction (edge of the read)
    - Excludes the last junction (edge of the read)
    """
#    tree_start = defaultdict(IntervalTree)
#    tree_end = defaultdict(IntervalTree)
#    truth_sites = []
#    truth_ids = set()
    truth = []

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

                truth.append((chrom1, int(e1), chrom2, int(s2)))

    return truth


#def find_truth_overlap(tree, chrom, pos, window):
#    """
#    Return truth breakpoint if overlap exists
#    """
#    if chrom not in tree:
#        return []
#
#    hits = tree[chrom].overlap(pos - window, pos + window)
#
#    if not hits:
#        return []
#
#    return [h.data for h in hits]


#def detect_cigar_events(read, min_event_size):
#    events = []
#
#    if read.is_unmapped or not read.cigartuples:
#        return events
#
#    # --- CIGAR-based events ---
#    ref_pos = read.reference_start
#    chrom = read.reference_name
#    mapq = read.mapping_quality
#
#    for op, length in read.cigartuples:
#        if op == 0:  # match
#            ref_pos += length
#        elif op in (2, 3):  # deletion / skip
#            if length >= min_event_size:
#                events.append(("DELETION_SKIP", chrom, ref_pos, mapq))
#            ref_pos += length
#        elif op == 1:  # insertion
#            if length >= min_event_size:
#                events.append(("INSERTION", chrom, ref_pos, mapq))
#        elif op == 4:  # soft clip
#            if length >= min_event_size:
#                events.append(("SOFT_CLIP", chrom, ref_pos, mapq))
#    return events


def detect_split_alignment_events(all_alignments, min_event_size):
    """
    Detect split-read events from multiple alignments (LAST --split, NGMLR, VACmap, ...).

    Parameters
    ----------
    all_alignments : list of pysam.AlignedSegment
        All alignments for the same read.
    min_event_size : int
        Minimum distance to call a breakpoint.

    Returns
    -------
    events
    """
    events = []

    # Sort alignments by query start
    aln_sorted = sorted(all_alignments, key=lambda x: x.query_alignment_start)

    for i in range(len(aln_sorted) - 1):
        a1 = aln_sorted[i]
        a2 = aln_sorted[i + 1]

        if a1.is_supplementary and a2.is_unmapped:
            continue

        # Different chromosome → always a breakpoint
        if a1.reference_name != a2.reference_name:
            events.append((a1.reference_name, a1.reference_end,
                           a2.reference_name, a2.reference_start))
        else:
            # Same chromosome → check gap between end of first and start of second
            gap = a2.reference_start - a1.reference_end
            if gap >= min_event_size:
                events.append((a1.reference_name, a1.reference_end,
                               a2.reference_name, a2.reference_start))

    return events


def detect_sa_events(read):
    """
    Detect structural variant events from the SA tag.

    Parameters
    ----------
    read : pysam.AlignedSegment
        The read to analyze.

    Returns
    -------
    events : list of tuples
        Each event is (event_type, chrom, pos, mapq)
    """
    events = []

    if not read.has_tag("SA"):
        return events

    # --- Split-read events via SA tag (BWA/Minimap2) ---
    sa_entries = read.get_tag("SA").split(";")
    for entry in sa_entries:
        if not entry:
            continue

        chrom, pos, strand, cigar, mq, nm = entry.split(",")

        events.append((read.reference_name, read.reference_end,
                        chrom, int(pos)))

    return events


#def get_primary_alignment(alignments):
#    primaries = [aln for aln in alignments if not aln.is_secondary and not aln.is_supplementary]
#
#    if primaries:
#        return primaries[0]
#
#    return max(alignments, key=lambda x: x.mapping_quality)


def match_bedpe(pred, truth, window):
    """
    Bidirectional tolerant matching
    """
    p1c, p1, p2c, p2 = pred
    t1c, t1, t2c, t2 = truth

    # orientation forward
    fwd = (p1c == t1c) and (abs(p1 - t1) <= window) and \
         (p2c == t2c) and (abs(p2 - t2) <= window)

    # orientation reverse
    rev = (p1c == t2c) and (abs(p1 - t2) <= window) and \
         (p2c == t1c) and (abs(p2 - t1) <= window)

    return fwd or rev

# faster matching is needed
def build_truth_index(truth, bin_size):
    index = defaultdict(list)

    for i, (c1, p1, c2, p2) in enumerate(truth):
        b1 = (c1, p1 // bin_size)
        b2 = (c2, p2 // bin_size)

        index[b1].append(i)
        index[b2].append(i)

    return index


def get_candidates(truth_index, chrom1, pos1, chrom2, pos2, bin_size):
    return set(truth_index.get((chrom1, pos1 // bin_size), [])) | set(truth_index.get((chrom2, pos2 // bin_size), []))


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-b", "--bam", required=True)
    parser.add_argument("-t", "--truth", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-w", "--window", default=50, type=int)
    parser.add_argument("-m", "--min-event-size", default=15, type=int)
    parser.add_argument("-b", "--bin-size", type=int, default=5000)

    args = parser.parse_args()

    # Load truth breakpoints
    truth = load_breakpoints_from_fasta(args.truth)
    truth_index = build_truth_index(truth, args.bin_size)
#    tree_start, tree_end, truth_sites, truth_ids = load_breakpoints_from_fasta(args.truth)

    reads_by_qname = defaultdict(list)
    bam = pysam.AlignmentFile(args.bam)
    for read in bam:
        if read.is_unmapped:
            continue
        reads_by_qname[read.query_name].append(read)
    bam.close()

    predicted_breakpoints = []

    for read_id, alignments in reads_by_qname.items():

         # CIGAR events - will always be FP if both sides of the breakpoint is evaluated

        predicted_breakpoints.extend(detect_split_alignment_events(alignments, args.min_event_size))

        primary = None
        for a in alignments:
            if not a.is_secondary and not a.is_supplementary:
                primary = a
                break
        if primary:
            predicted_breakpoints.extend(detect_sa_events(primary))

    matched_truth = set()
    TP = 0
    FP = 0

    for pred in predicted_breakpoints:
        c1, p1, c2, p2 = pred

        cand_ids = get_candidates(truth_index, c1, p1, c2, p2, args.bin_size)

        found = False

        for i in cand_ids:
            if match_bedpe(pred, truth[i], args.window):
                matched_truth.add(i)
                found = True
                break
        if found:
            TP += 1
        else:
            FP += 1
    FN = len(truth) - len(matched_truth)

    precision = TP / (TP + FP) if TP + FP > 0 else 0
    recall = TP / (TP + FN) if TP + FN > 0 else 0
    fscore = 2 * precision * recall / (precision + recall) if precision + recall > 0 else 0

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
