import pandas as pd
import matplotlib.pyplot as plt
import argparse
import numpy as np
import seaborn as sns

parse = argparse.ArgumentParser()
parse.add_argument('--input_path', help='')
parse.add_argument('--output_path', help='')
parse.add_argument('--bwa', help='', action='store_true')
parse.add_argument('--bwa_x', help='', action='store_true')
parse.add_argument('--bwa_dodi', help='', action='store_true')
parse.add_argument('--bwa_flags', help='', action='store_true')
parse.add_argument('--bwa_flags_dodi', help='', action='store_true')
parse.add_argument('--minimap2', help='', action='store_true')
parse.add_argument('--minimap2_dodi', help='', action='store_true')
parse.add_argument('--lastal_dodi', help='', action='store_true')
parse.add_argument('--lastal_lastsplit', help='', action='store_true')
parse.add_argument('--lastalsplit', help='', action='store_true')
parse.add_argument('--ngmlr', help='', action='store_true')
parse.add_argument('--ngmlr_x', help='', action='store_true')
parse.add_argument('--vacmap_r', help='', action='store_true')
parse.add_argument('--vacmap_s', help='', action='store_true')
args = parse.parse_args()

data = {}
benchmark_res={}
if args.bwa:
    file = args.input_path + '/bwa.mappings_labelled.csv'
    bwa=pd.read_csv(file, sep='\t')
    data['bwa'] = bwa
    benchmark=pd.read_csv(args.input_path + '/bwa.benchmark_res3.csv', sep='\t')
    benchmark_res['bwa'] = benchmark
if args.bwa_x:
    file = args.input_path + '/bwa_x.mappings_labelled.csv'
    bwa_x=pd.read_csv(file, sep='\t')
    data['bwa_x'] = bwa_x
    benchmark=pd.read_csv(args.input_path + '/bwa_x.benchmark_res3.csv', sep='\t')
    benchmark_res['bwa_x'] = benchmark
if args.bwa_dodi:
    file = args.input_path + '/bwa_dodi.mappings_labelled.csv'
    bwa_dodi=pd.read_csv(file, sep='\t')
    data['bwa_dodi'] = bwa_dodi
    benchmark=pd.read_csv(args.input_path + '/bwa_dodi.benchmark_res3.csv', sep='\t')
    benchmark_res['bwa_dodi'] = benchmark
if args.bwa_flags_dodi:
    file = args.input_path + '/bwa_flags_dodi.mappings_labelled.csv'
    bwa_flags_dodi = pd.read_csv(file, sep='\t')
    data['bwa_flags_dodi'] = bwa_flags_dodi
    benchmark=pd.read_csv(args.input_path + '/bwa_flags_dodi.benchmark_res3.csv', sep='\t')
    benchmark_res['bwa_flags_dodi'] = benchmark
if args.bwa_flags:
    file = args.input_path + '/bwa_flags.mappings_labelled.csv'
    bwa_flags = pd.read_csv(file, sep='\t')
    data['bwa_flags'] = bwa_flags
    benchmark=pd.read_csv(args.input_path + '/bwa_flags.benchmark_res3.csv', sep='\t')
    benchmark_res['bwa_flags'] = benchmark
if args.minimap2:
    file = args.input_path + '/minimap2.mappings_labelled.csv'
    minimap2=pd.read_csv(file, sep='\t')
    data['minimap2'] = minimap2
    benchmark=pd.read_csv(args.input_path + '/minimap2.benchmark_res3.csv', sep='\t')
    benchmark_res['minimap2'] = benchmark
if args.minimap2_dodi:
    file = args.input_path + '/minimap2_dodi.mappings_labelled.csv'
    minimap2_dodi=pd.read_csv(file, sep='\t')
    data['minimap2_dodi'] = minimap2_dodi
    benchmark=pd.read_csv(args.input_path + '/minimap2_dodi.benchmark_res3.csv', sep='\t')
    benchmark_res['minimap2_dodi'] = benchmark
if args.lastalsplit:
    file = args.input_path + '/lastalsplit.mappings_labelled.csv'
    lastalsplit=pd.read_csv(file, sep='\t')
    data['lastalsplit'] = lastalsplit
    benchmark=pd.read_csv(args.input_path + '/lastalsplit.benchmark_res3.csv', sep='\t')
    benchmark_res['lastalsplit'] = benchmark
if args.lastal_lastsplit:
    file = args.input_path + '/lastal_lastsplit.mappings_labelled.csv'
    lastal_lastsplit=pd.read_csv(file, sep='\t')
    data['lastal_lastsplit'] = lastal_lastsplit
    benchmark=pd.read_csv(args.input_path + '/lastal_split.benchmark_res3.csv', sep='\t')
    benchmark_res['lastal_split'] = benchmark
if args.lastal_dodi:
    file = args.input_path + '/lastal_dodi.mappings_labelled.csv'
    lastal_dodi=pd.read_csv(file, sep='\t')
    data['lastal_dodi'] = lastal_dodi
    benchmark=pd.read_csv(args.input_path + '/lastal_dodi.benchmark_res3.csv', sep='\t')
    benchmark_res['lastal_dodi'] = benchmark
if args.ngmlr:
    file = args.input_path + '/ngmlr.mappings_labelled.csv'
    ngmlr=pd.read_csv(file, sep='\t')
    data['ngmlr'] = ngmlr
    benchmark=pd.read_csv(args.input_path + '/ngmlr.benchmark_res3.csv', sep='\t')
    benchmark_res['ngmlr'] = benchmark
if args.ngmlr_x:
    file = args.input_path + '/ngmlr_x.mappings_labelled.csv'
    ngmlr_x = pd.read_csv(file, sep='\t')
    data['ngmlr_x'] = ngmlr_x
    benchmark=pd.read_csv(args.input_path + '/ngmlr_x.benchmark_res3.csv', sep='\t')
    benchmark_res['ngmlr_x'] = benchmark
if args.vacmap_r:
    file = args.input_path + '/vacmap_r.mappings_labelled.csv'
    vacmap_r=pd.read_csv(file, sep='\t')
    data['vacmap_r'] = vacmap_r
    benchmark=pd.read_csv(args.input_path + '/vacmap_r.benchmark_res3.csv', sep='\t')
    benchmark_res['vacmap_r'] = benchmark
if args.vacmap_s:
    file = args.input_path + '/vacmap_s.mappings_labelled.csv'
    vacmap_s=pd.read_csv(file, sep='\t')
    data['vacmap_s'] = vacmap_s
    benchmark=pd.read_csv(args.input_path + '/vacmap_s.benchmark_res3.csv', sep='\t')
    benchmark_res['vacmap_s'] = benchmark

colors = {'bwa': '#e69f00', 'bwa_dodi': 'firebrick',
          'bwa_flags': 'gold', 'bwa_flags_dodi': 'lightcoral',
          'bwa_x': 'red',
          'minimap2': '#009e73', 'minimap2_dodi': 'darkgreen',
          'lastalsplit': '#56b4e9', 'lastal_dodi':'lightblue',
          'ngmlr': '#cc79a7', 'ngmlr_x': 'chocolate',
          'vacmap_r': '#d55e00', 'vacmap_s': 'magenta'}

scale = 0.01
def create_bins(dict, base):
    for name, df in dict.items():
        bins = []
        for i in df['aln_size']:
            bins.append(base * round(i/base))
        dict[name] = df.assign(bins=bins)
    return dict


if 'short' in args.input_path:
    base = 25
    cutoff = 600
    title = 'short'
    total = 99127
if 'medium' in args.input_path:
    base = 50
    cutoff = 1200
    title = 'medium'
    total=99844
if 'long' in args.input_path:
    cutoff = 1500
    base = 100
    title = 'long'
    total=97333

data = create_bins(data, base)
data = create_bins(benchmark_res, base)


def precision_aln_size(data):
    for name, df in data.items():
        bin_precision = []
        bin_id = []
        s = []
        for bid, b in df.groupby('bins'):
            if bid > cutoff:
                break
            if (b['tp'].sum() + b['fp'].sum()) == 0:
                continue
            s.append(len(b)*scale)
            bin_precision.append(b['tp'].sum() / (b['tp'].sum() + b['fp'].sum()))
            bin_id.append(bid)

        plt.plot(bin_id, bin_precision, label=name, c=colors[name], alpha=0.8)
        plt.scatter(bin_id, bin_precision, s=s, alpha=0.25, c=colors[name], linewidths=0)

    plt.legend(loc='best', fontsize='xx-small')
    plt.xscale("log")
    plt.xlabel('Alignment size')
    plt.ylabel('Precision')
    plt.ylim(0, 1.1)
    plt.tight_layout()
    plt.savefig(args.output_path + '/size_vs_precision.png', dpi=600)
    #plt.show()
    plt.close()


precision_aln_size(data)


def wrong_plot_bins(data):
    for name, df in data.items():
        bin_wrong = []
        bin_w = []
        wrong = 0
        s=[]
        for bid, b in df.groupby('bins'):
            bin_wrong.append(wrong / total * 100)
            bin_w.append(bid)
            wrong += b['fp'].sum()
            s.append(len(b) * scale)
        plt.plot(bin_w, bin_wrong, label=name, c=colors[name], alpha=0.8)
        plt.scatter(bin_w, bin_wrong, s=s, alpha=0.25, c=colors[name], linewidths=0)
    plt.legend(loc='best', fontsize='xx-small')
    #plt.xscale("log")
    plt.xlabel('Alignment size')
    plt.ylabel('False positive %')
    plt.tight_layout()
    plt.savefig(args.output_path + '/aln_size_vs_fp_bin.png', dpi=600)
    #plt.show()
    plt.close()


wrong_plot_bins(benchmark_res)



def wrong_plot_bins2(data):
    for name, df in data.items():
        bin_wrong = []
        bin_w = []
        wrong = 0
        s=[]
        for bid, b in df.groupby('bins'):
            bin_wrong.append(wrong / total * 100)
            bin_w.append(bid)
            wrong += b['fn'].sum()
            s.append(len(b)*scale)
        plt.plot(bin_w, bin_wrong, label=name, c=colors[name], alpha=0.8)
        plt.scatter(bin_w, bin_wrong, s=s, alpha=0.25, c=colors[name], linewidths=0)

    plt.legend(loc='best', fontsize='xx-small')
    #plt.xscale("log")
    plt.xlabel('Alignment size')
    plt.ylabel('False negative %')
    plt.tight_layout()
    plt.savefig(args.output_path + '/aln_size_vs_fn_bin.png', dpi=600)
    #plt.show()
    plt.close()


wrong_plot_bins2(benchmark_res)


def precision_mapq(data):
    data2 = {}
    if args.lastal_dodi:
        data2['lastal_dodi'] = lastal_dodi
    if args.lastal_lastsplit:
        data2['lastal_lastsplit'] = lastal_lastsplit
    if args.lastalsplit:
        data2['lastalsplit'] = lastalsplit

    for name, df in data2.items():
        bin_precision = []
        bin_id = []
        s = []
        for bid, b in df.groupby('mapq'):
            if len(b) == 0:
                continue
            prec = b['tp'].sum() / (b['tp'].sum() + b['fp'].sum())
            s.append(len(b)*scale)
            bin_precision.append(prec)
            bin_id.append(bid)

        plt.plot(bin_id, bin_precision, label=name, c=colors[name], alpha=0.8)
        plt.scatter(bin_id, bin_precision, s=s, alpha=0.25, c=colors[name], linewidths=0)

    plt.legend(loc='best', fontsize='xx-small')
    plt.xlabel('MapQ')
    plt.ylabel('Precision')
    plt.ylim(0, 1.1)
    plt.tight_layout()
    plt.savefig(args.output_path + '/mapq_vs_precision1.png', dpi=600)
    #plt.show()
    plt.close()

    data3 = {}
    if args.bwa:
        data3['bwa'] = bwa
    if args.bwa_dodi:
        data3['bwa_dodi'] = bwa_dodi
    if args.bwa_x:
        data3['bwa_x'] = bwa_x
    if args.bwa_flags_dodi:
        data3['bwa_flags_dodi'] = bwa_flags_dodi
    if args.bwa_flags:
        data3['bwa_flags'] = bwa_flags
    if args.minimap2:
        data3['minimap2'] = minimap2
    if args.minimap2_dodi:
        data3['minimap2_dodi'] = minimap2_dodi
    if args.ngmlr:
        data3['ngmlr'] = ngmlr
    if args.ngmlr_x:
        data3['ngmlr_x'] = ngmlr_x
    if args.vacmap_r:
        data3['vacmap_r'] = vacmap_r
    if args.vacmap_s:
        data3['vacmap_s'] = vacmap_s

    for name, df in data3.items():
        bin_precision = []
        bin_id = []
        s = []
        for bid, b in df.groupby('mapq'):
            if len(b) == 0:
                continue
            s.append(len(b)*scale)
            bin_precision.append(b['tp'].sum() / (b['tp'].sum() + b['fp'].sum()))
            bin_id.append(bid)

        plt.plot(bin_id, bin_precision, label=name, c=colors[name], alpha=0.8)
        plt.scatter(bin_id, bin_precision, s=s, alpha=0.25, c=colors[name], linewidths=0)

    plt.legend(loc='best', fontsize='xx-small')
    plt.xlabel('MapQ')
    plt.ylabel('Precision')
    plt.ylim(0, 1.1)
    plt.tight_layout()
    plt.savefig(args.output_path + '/mapq_vs_precision2.png', dpi=600)
    #plt.show()
    plt.close()


precision_mapq(data)


def addlabels(x,y):
    for i in range(len(x)):
        plt.text(i, y[i], y[i], ha='center')


def mapped_alignments(data):
    counts = []
    c = []
    n = []
    for name, df in data.items():
        counts.append(len(df))
        c.append(colors[name])
        n.append(name)

    fig, ax = plt.subplots(figsize=(15,15))
    ax.bar(data.keys(), counts, label=data.keys(), color=c)
    addlabels(data.keys(), counts)
    ax.set_ylabel('Mapped fragments')
    ax.set_title(title)
    ax.legend(title='aligners', bbox_to_anchor=(0.9, 0.75))
    plt.setp(ax.get_xticklabels(), rotation=45, horizontalalignment='right')
    plt.subplots_adjust(bottom=0.2)
    plt.savefig(args.output_path + '/n_fragments.png', dpi=600)
    plt.close()


mapped_alignments(data)


def fragment_length_dist(data):
    plt.figure()
    for name, df in data.items():
        sns.kdeplot(df['aln_size'], alpha=0.7, label = name, color=colors[name])
    plt.xlabel('fragment length')
    plt.ylabel('density function')
    plt.legend(fontsize='xx-small')
    plt.savefig(args.output_path + '/fragment_len_dist.png', dpi=600)
    plt.close()


fragment_length_dist(data)


markers={'bwa': 'p', 'bwa_dodi': 'P',
          'bwa_flags': '.', 'bwa_flags_dodi': 'o',
          'bwa_x': '*',
          'minimap2': 'x', 'minimap2_dodi': 'X',
          'lastalsplit': 'P', 'lastal_dodi':'+',
          'ngmlr': 'd', 'ngmlr_x': 'D',
          'vacmap_r': 'h', 'vacmap_s': 'H'}


def BWA_curve(data):
    # BWA-MEM plot
    for name, df in data.items():
        if name == 'lastalsplit':
            continue
        x = []
        y = []
        s=[]
        tp = df['tp'].sum()
        fp = df['fp'].sum()
        for i, b in df.groupby('mapq'):
            tp -= b['tp'].sum()
            fp -= b['fp'].sum()
            if tp+fp == 0:
                continue
            y.append((fp+tp)/total)
            x.append(fp/(tp+fp))
            s.append(len(b)*scale)
        plt.plot(x, y, alpha=0.8, c=colors[name], label=name)
        plt.scatter(x, y, s=s, alpha=0.25, linewidths=0, c=colors[name])

    plt.ylabel('mapped/total')
    plt.xlabel('wrong/mapped')
    plt.legend(loc='best', fontsize='xx-small')
    plt.savefig(args.output_path + '/bwamempaper_mapq.png', dpi=600)
    plt.close()


BWA_curve(benchmark_res)


def BWA_curve2(data):
    data2={}
    for name,df in data.items():
        if name=='lastalsplit':
            continue
        mapq = []
        tp = []
        fp = []
        fn = []
        for i, b in df.groupby('qname'):
            mapq.append(b['mapq'].mean())
            tp.append(b['tp'].sum())
            fp.append(b['fp'].sum())
            fn.append(b['fn'].sum())
        data2[name] = pd.DataFrame({'mapq':mapq, 'tp':tp, 'fp':fp, 'fn':fn})
    # BWA-MEM plot
    for name, df in data2.items():
        if name == 'lastalsplit':
            continue
        x = []
        y = []
        s=[]
        tp = df['tp'].sum()
        fp = df['fp'].sum()
        for i, b in df.groupby('mapq'):
            if tp+fp == 0:
                continue
            y.append((fp+tp)/total * 100)
            x.append(fp/(tp+fp))
            s.append(len(b)*scale)
            tp -= b['tp'].sum()
            fp -= b['fp'].sum()
        plt.plot(x, y, alpha=0.8, c=colors[name], label=name)
        plt.scatter(x, y, s=s, alpha=0.25, linewidths=0, c=colors[name])

    plt.ylabel('mapped/total')
    plt.xlabel('wrong/mapped')
    plt.legend(loc='best', fontsize='xx-small')
    plt.savefig(args.output_path + '/bwamempaper_mapq2.png', dpi=600)
    plt.close()


BWA_curve2(benchmark_res)


def recall_aln_size(data):
    for name, df in data.items():
        bin_recall = []
        bin_id = []
        s = []
        for bid, b in df.groupby('bins'):
            if bid > cutoff:
                break
            if (b['tp'].sum() + b['fn'].sum()) == 0:
                continue
            s.append(len(b)*scale)
            bin_recall.append(b['tp'].sum() / (b['tp'].sum() + b['fn'].sum()))
            bin_id.append(bid)

        plt.plot(bin_id, bin_recall, label=name, c=colors[name], alpha=0.8)
        plt.scatter(bin_id, bin_recall, s=s, alpha=0.25, c=colors[name], linewidths=0)

    plt.legend(loc='best', fontsize='xx-small')
    plt.xscale("log")
    plt.xlabel('Alignment size')
    plt.ylabel('Recall')
    plt.ylim(0, 1.1)
    plt.tight_layout()
    plt.savefig(args.output_path + '/size_vs_recall.png', dpi=600)
    #plt.show()
    plt.close()


recall_aln_size(benchmark_res)



def precision_recall(data):
    data2={}
    for name,df in data.items():
        if name=='lastalsplit':
            continue
        mapq = []
        tp = []
        fp = []
        fn = []
        for i, b in df.groupby('qname'):
            mapq.append(b['mapq'].mean())
            tp.append(b['tp'].sum())
            fp.append(b['fp'].sum())
            fn.append(b['fn'].sum())
        data2[name] = pd.DataFrame({'mapq':mapq, 'tp':tp, 'fp':fp, 'fn':fn})
    for name, df in data2.items():
        recall = []
        precision = []
        tp = 0
        fp = 0
        fn = 0
        s = []
        for i, b in df.groupby('mapq'):
            if len(b) < 0:
                continue
            tp += b['tp'].sum()
            fp += b['fp'].sum()
            fn += b['fn'].sum()
            if tp + fp == 0 or tp + fn == 0:
                continue
            precision.append(tp / (tp+fp))
            recall.append(tp / (tp+fn))
            s.append(len(b) * scale)
        plt.plot(recall, precision, alpha=0.8, label=name, c=colors[name])
        plt.scatter(recall, precision, s=s, alpha=0.25, linewidths=0, c=colors[name])
        # plt.gca().invert_xaxis()
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.legend(loc='best', fontsize='xx-small')
    plt.savefig(args.output_path + '/Precision-Recall.png', dpi=600)
    plt.close()


precision_recall(benchmark_res)


def tp_fp(data):
    for name, df in data.items():
        if name=='lastalsplit':
            continue
        x = []
        y = []
        s=[]
        tp = df['tp'].sum()
        fp = df['fp'].sum()
        for i, b in df.groupby('mapq'):
            tp -= b['tp'].sum()
            fp -= b['fp'].sum()
            if tp+fp == 0:
                continue
            y.append(tp/total)
            x.append(fp/total)
            s.append(len(b)*scale)
        plt.plot(x, y, alpha=0.8, c=colors[name], label=name)
        plt.scatter(x, y, s=s, alpha=0.25, linewidths=0, c=colors[name])

    plt.ylabel('tp/total')
    plt.xlabel('fp/total')
    plt.legend(loc='best', fontsize='xx-small')
    plt.savefig(args.output_path + '/tp_fp.png', dpi=600)
    plt.close()


tp_fp(benchmark_res)


def ROC(data):
    data2={}
    for name,df in data.items():
        if name=='lastalsplit':
            continue
        mapq = []
        tp = []
        fp = []
        fn = []
        for i, b in df.groupby('qname'):
            mapq.append(b['mapq'].mean())
            tp.append(b['tp'].sum())
            fp.append(b['fp'].sum())
            fn.append(b['fn'].sum())
        data2[name] = pd.DataFrame({'mapq':mapq, 'tp':tp, 'fp':fp, 'fn':fn})
    for name, df in data2.items():
        x = []
        y = []
        tp = 0
        fp = 0
        fn = 0
        s = []
        for i, b in df.groupby('mapq'):
            if len(b) < 0:
                continue
            tp += b['tp'].sum()
            fp += b['fp'].sum()
            fn += b['fn'].sum()
            if tp + fp == 0 or tp + fn == 0:
                continue
            x.append(fp / (fp+len(df)-len(b)))
            y.append(tp / (tp+fn))
            s.append(len(b) * scale)
        plt.plot(x, y, alpha=0.8, label=name, c=colors[name])
        plt.scatter(x, y, s=s, alpha=0.25, linewidths=0, c=colors[name])
        # plt.gca().invert_xaxis()
    plt.xlabel('False positive rate')
    plt.ylabel('True positive rate')
    plt.legend(loc='best', fontsize='xx-small')
    plt.savefig(args.output_path + '/ROC.png', dpi=600)
    plt.close()


ROC(benchmark_res)