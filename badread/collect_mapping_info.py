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


def mapping_info(f, outf):

    af = pysam.AlignmentFile(f, 'r')
    d = defaultdict(list)
    for a in af.fetch(until_eof=True):
        if not a.flag & 4:
            d[a.qname].append(a)

    res = []
    no = 0
    yes = 0
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
        n_aligns = len(v)
        any_seq = False

        temp = []
        sec = 0
        sup = 0
        for index, a in enumerate(v):
            if a.is_secondary:
                sec = 1
            if a.is_supplementary:
                sup = 1
            qstart, qend, qlen = get_query_pos_from_cigartuples(a)
            align_reverse = bool(a.flag & 16)
            if primary_reverse != align_reverse:
                start_temp = qlen - qend
                qend = start_temp + qend - qstart
                qstart = start_temp
            pri = index == pri_index
            if not pri:
                no += 1
            else:
                yes += 1
                any_seq = len(seq) if seq else 0

            chrom = af.get_reference_name(a.rname)
            start = a.reference_start + 1
            end = a.reference_end
            t = pd.Interval(start, end)

            rd = {'qname': a.qname,
                 'n_alignments': n_aligns,
                 'chrom': chrom,
                 'rstart': start,
                 'rend': end,
                 'strand': '-' if align_reverse else '+',
                 'qstart': qstart,
                 'qend': qend,
                 'qlen': qlen,
                 'aln_size': qend - qstart,
                 'mapq': a.mapq,
                 'alignment_score': (0 if 'AS' not in a.get_tags() else a.get_tag('AS')),
                 'seq': seq if pri else '',
                 'is_secondary': sec,
                 'is_supplementary': sup,
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
             'alignment_score', 'short_anchor<50bp', 'seq', 'is_secondary',
            'is_supplementary']

    df = df[cols]
    df.to_csv(f"{outf}.bed", index=False, sep="\t")


def collect_mapping_info(args):
    mapping_info(args.bam, args.out)
    print('Done')
