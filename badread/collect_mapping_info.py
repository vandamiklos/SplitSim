import pysam
import pandas as pd
from collections import defaultdict
import click

def get_query_pos_from_cigartuples(r):
    # Infer the position on the query sequence of the alignment using cigar string
    start = 0
    query_length = r.infer_read_length()  # Note, this also counts hard-clips
    end = query_length
    if r.cigartuples[0][0] == 4 or r.cigartuples[0][0] == 5:
        start += r.cigartuples[0][1]
    if r.cigartuples[-1][0] == 4 or r.cigartuples[-1][0] == 5:
        end -= r.cigartuples[-1][1]
    return start, end, query_length


def split_on_gaps(primary, min_gap=15):
    # if an alignment has a gap longer than 15 split that alignment
    fragments = []
    ref_pos = primary.reference_start
    read_pos = 0

    frag_start_ref = ref_pos
    frag_start_read = read_pos

    for cigar, length in primary.cigartuples:

        consumes_ref = cigar in (0, 2, 3, 7, 8)
        consumes_read = cigar in (0, 1, 4, 7, 8)

        # Large reference gap → breakpoint
        if cigar in (2, 3) and length > min_gap:
            # close current fragment (before gap)
            fragments.append({
                "ref_start": frag_start_ref,
                "ref_end": ref_pos,
                "read_start": frag_start_read,
                "read_end": read_pos
            })

            # move reference across gap
            ref_pos += length

            # start new fragment after gap
            frag_start_ref = ref_pos
            frag_start_read = read_pos
            continue

        # advance coordinates normally
        if consumes_ref:
            ref_pos += length
        if consumes_read:
            read_pos += length

    # close final fragment
    fragments.append({
        "ref_start": frag_start_ref,
        "ref_end": ref_pos,
        "read_start": frag_start_read,
        "read_end": read_pos
    })

    return fragments


def mapping_info(f, outf):

    af = pysam.AlignmentFile(f, 'r')
    d = defaultdict(list)
    for a in af.fetch(until_eof=True):
        if not a.flag & 4:
            d[a.qname].append(a)

    res = []
    for qname, v in d.items():
        flag = [(index, i) for index, i in enumerate(v) if not i.flag & 2304]
        # no primary flag set
        if len(flag) == 0:
            print('Flag problem: no primary alignment')
            flag = [(index, i) for index, i in enumerate(v) if not i.flag]
            if len(flag) == 0:
                continue
            flag = flag[0]
        if len(flag) > 1:  # todo check bug in dodi, not currently setting primary alignment flag properly
            flag = [flag[flag.index(max(flag, key=lambda x: x[1].get_tag('AS')))]]
    # try to catch the errors
        #if len(flag) != 1:
        #    print('Error in ', f, 'flag problem', len(flag), [i.flag for i in v])
        #    quit()
        pri_index, pri_read = flag[0]
        primary_reverse = bool(pri_read.flag & 16)
        seq = pri_read.get_forward_sequence()
        fragments = split_on_gaps(pri_read, min_gap=15)
        n_aligns = len(v) + len(fragments) - 1
        any_seq = False

        temp = []

        for index, a in enumerate(v):

            qstart, qend, qlen = get_query_pos_from_cigartuples(a)
            align_reverse = bool(a.flag & 16)

            # strand correction relative to primary
            if primary_reverse != align_reverse:
                start_temp = qlen - qend
                qend = start_temp + qend - qstart
                qstart = start_temp

            chrom = af.get_reference_name(a.rname)
            mapq = a.mapq
            aln_score = a.get_tag('AS') if a.has_tag('AS') else 0

            pri = index == pri_index

            # -------------------------------
            # PRIMARY → replace with fragments if large gap in the alignment
            # -------------------------------
            if pri:

                for frag in fragments:

                    frag_qstart = frag["read_start"]
                    frag_qend = frag["read_end"]

                    # strand correction
                    if primary_reverse:
                        start_temp = qlen - frag_qend
                        frag_qend = start_temp + frag_qend - frag_qstart
                        frag_qstart = start_temp

                    rd = {
                        'qname': a.qname,
                        'chrom': chrom,
                        'rstart': frag["ref_start"] + 1,
                        'rend': frag["ref_end"],
                        'strand': '-' if align_reverse else '+',
                        'qstart': frag_qstart,
                        'qend': frag_qend,
                        'qlen': qlen,
                        'aln_size': frag_qend - frag_qstart,
                        'mapq': mapq,
                        'alignment_score': aln_score,
                        'seq': seq,  # keep full sequence only once
                    }

                    temp.append(rd)

            # ---------------------------------
            # secondary / supplementary
            # ---------------------------------
            else:

                rd = {
                    'qname': a.qname,
                    'chrom': chrom,
                    'rstart': a.reference_start + 1,
                    'rend': a.reference_end,
                    'strand': '-' if align_reverse else '+',
                    'qstart': qstart,
                    'qend': qend,
                    'qlen': qlen,
                    'aln_size': qend - qstart,
                    'mapq': mapq,
                    'alignment_score': aln_score,
                    'seq': '',
                }

                temp.append(rd)

        if not any_seq:
            print('missing', qname, [(len(vv.seq), vv.infer_query_length()) if vv.seq else vv.infer_query_length() for vv in v])
            quit()

        if len(temp) > 1:
            res += temp
            continue

        res += temp

    df = pd.DataFrame.from_records(res).sort_values(['qname', 'qstart'])

    bad_anchors = []
    # flag reads with small anchoring alignments
    for grp, d in df.groupby('qname'):
        aln_s = list(d['aln_size'])
        if aln_s[0] < 50 or aln_s[-1] < 50:
            bad_anchors += [1] * len(d)
        else:
            bad_anchors += [0] * len(d)
    df['short_anchor<50bp'] = bad_anchors

    df = df.sort_values(['n_alignments', 'qname', 'qstart'], ascending=[False, True, True])

    cols = ['chrom', 'rstart', 'rend', 'qname', 'n_alignments', 'aln_size', 'qstart', 'qend', 'strand', 'mapq', 'qlen',
             'alignment_score', 'short_anchor<50bp', 'seq']

    df = df[cols]
    df.to_csv(f"{outf}.bed", index=False, sep="\t")


def collect_mapping_info(args):
    mapping_info(args.bam, args.out)
    print('Done')
