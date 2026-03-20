import pysam
import argparse
from intervaltree import IntervalTree
from collections import defaultdict


def load_breakpoints(breakpoint_file, window):
    """
    Build interval trees from truth breakpoints
    """
    trees = defaultdict(IntervalTree)
    truth_sites = defaultdict(set)

    with open(breakpoint_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            chrom, pos = line.strip().split()[:2]
            pos = int(pos)

            trees[chrom].add(pos - window, pos + window, pos)
            truth_sites[chrom].add(pos)

    return trees, truth_sites


def find_truth_overlap(trees, chrom, pos):
    """
    Return truth breakpoint if overlap exists
    """
    if chrom not in trees:
        return None

    hits = trees[chrom].overlap(pos)

    if not hits:
        return None

    # return closest truth breakpoint
    return min(hits, key=lambda x: abs(x.data - pos)).data


def detect_events(read):

    events = []

    if read.is_unmapped or not read.cigartuples:
        return events

    ref_pos = read.reference_start
    chrom = read.reference_name
    mapq = read.mapping_quality

    for op, length in read.cigartuples:

        if op == 0:  # match
            ref_pos += length

        elif op in (2, 3):  # deletion / skip
            events.append(("DELETION_SKIP", chrom, ref_pos, mapq))
            ref_pos += length

        elif op == 1:  # insertion
            events.append(("INSERTION", chrom, ref_pos, mapq))

        elif op == 4:  # soft clip
            events.append(("SOFT_CLIP", chrom, ref_pos, mapq))

    # split reads
    if read.has_tag("SA"):

        sa_entries = read.get_tag("SA").split(";")

        for entry in sa_entries:
            if not entry:
                continue

            chrom2, pos, strand, cigar, mq, nm = entry.split(",")

            events.append(("SPLIT_READ", chrom2, int(pos), int(mq)))

    return events


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("-b", "--bam", required=True)
    parser.add_argument("-t", "--truth", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-w", "--window", default=50, type=int)

    args = parser.parse_args()

    trees, truth_sites = load_breakpoints(args.truth, args.window)

    matched_truth = set()
    results = []

    bam = pysam.AlignmentFile(args.bam)

    for read in bam:

        events = detect_events(read)

        for event_type, chrom, pos, mapq in events:

            truth_hit = find_truth_overlap(trees, chrom, pos)

            if truth_hit is not None:

                results.append((
                    read.query_name,
                    event_type,
                    chrom,
                    pos,
                    mapq,
                    truth_hit,
                    "TP"
                ))

                matched_truth.add((chrom, truth_hit))

            else:

                results.append((
                    read.query_name,
                    event_type,
                    chrom,
                    pos,
                    mapq,
                    "NA",
                    "FP"
                ))

    bam.close()

    # compute FN
    all_truth = {(chrom, pos) for chrom in truth_sites for pos in truth_sites[chrom]}
    fn_sites = all_truth - matched_truth

    TP = len(matched_truth)
    FP = sum(1 for r in results if r[-1] == "FP")
    FN = len(fn_sites)

    precision = TP / (TP + FP) if TP + FP > 0 else 0
    recall = TP / (TP + FN) if TP + FN > 0 else 0
    fscore = 2 * precision * recall / (precision + recall) if precision + recall > 0 else 0

    # write event results
    with open(args.output + ".events.csv", "w") as out:

        out.write("read_id\tevent_type\tchrom\tpos\tmapq\ttruth_pos\tclass\n")

        for r in results:
            out.write("\t".join(map(str, r)) + "\n")

    # write FN list
    with open(args.output + ".fn.csv", "w") as out:

        out.write("chrom\ttruth_pos\n")

        for chrom, pos in sorted(fn_sites):
            out.write(f"{chrom}\t{pos}\n")

    # write summary
    with open(args.output + ".summary.txt", "w") as out:

        out.write(f"TP\t{TP}\n")
        out.write(f"FP\t{FP}\n")
        out.write(f"FN\t{FN}\n")
        out.write(f"Precision\t{precision:.4f}\n")
        out.write(f"Recall\t{recall:.4f}\n")
        out.write(f"F-score\t{fscore:.4f}\n")

    print("Finished")
    print(f"TP: {TP}")
    print(f"FP: {FP}")
    print(f"FN: {FN}")
    print(f"Precision: {precision:.4f}")
    print(f"Recall: {recall:.4f}")
    print(f"F-score: {fscore:.4f}")


if __name__ == "__main__":
    main()
