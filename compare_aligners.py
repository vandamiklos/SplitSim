import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import argparse
import numpy as np
import seaborn as sns

parse = argparse.ArgumentParser()
parse.add_argument('--input_path', help='')
parse.add_argument('--output_path', help='')
parse.add_argument('--aligner_name', help='List of aligner names',
    nargs='+',  # Accept one or more arguments
    type=str)
parse.add_argument('--simple', action="store_true")
args = parse.parse_args()

data = {}
benchmark_res={}

if args.simple:
    type = args.input_path.split('/')[-1]
    stats = pd.read_csv(args.input_path + '/' + args.aligner_name[0] + '.' + type + '_stats.txt', sep='\t')
    total = int(stats['target_n'].iloc[0])
    for aligner in args.aligner_name:
        data[aligner] = pd.read_csv(args.input_path + '/' + aligner + '.' + type + '_mappings_labelled.csv', sep='\t')
        benchmark_res[aligner] = pd.read_csv(args.input_path + '/' + aligner + '.' + type + '_benchmark_res_fn.csv', sep='\t')

if not args.simple:
    for aligner in args.aligner_name:
        data[aligner] = pd.read_csv(args.input_path + '/' + aligner + '.mappings_labelled.csv', sep='\t')
        benchmark_res[aligner] = pd.read_csv(args.input_path + '/' + aligner + '.benchmark_res_fn.csv', sep='\t')
        stats = pd.read_csv(args.input_path + '/' + args.aligner_name[0] + '.stats.txt', sep='\t')
        total = int(stats['target_n'].iloc[0])

colors = {'bwa': '#F0A430FF', 'bwa_dodi': 'firebrick',
          'bwa_flags': 'gold', 'bwa_flags_dodi': 'lightcoral',
          'bwa_x': 'red',
          'minimap2': '#800000FF', 'minimap2_dodi': 'darkgreen',
          'lastalsplit': '#1B3A54FF', 'lastal_dodi':'lightblue',
          'ngmlr': '#205F4BFF', 'ngmlr_x': 'chocolate',
          'vacmap_r': '#774762FF', 'vacmap_s': '#774762FF'}


markers={'bwa': 'o', 'bwa_dodi': 'P',
          'bwa_flags': '.', 'bwa_flags_dodi': 'o',
          'bwa_x': '*',
          'minimap2': 'x', 'minimap2_dodi': 'X',
          'lastalsplit': '+', 'lastal_dodi':'P',
          'ngmlr': 'p', 'ngmlr_x': 'p',
          'vacmap_r': 's', 'vacmap_s': 's'}

scale = 0.01


def create_bins(dict, base, column):
    for name, df in dict.items():
        bins = []
        for i in df[column]:
            bins.append(base * round(i/base))
        dict[name] = df.assign(bins=bins)
    return dict


if 'short' in args.input_path:
    base = 25
    cutoff = 600
    title = 'short'
if 'medium' in args.input_path:
    base = 25
    cutoff = 1200
    title = 'medium'
if 'long' in args.input_path:
    cutoff = 1500
    base = 25
    title = 'long'

data = create_bins(data, base, 'aln_size')
benchmark_res = create_bins(benchmark_res, base, 'aln_size')


def precision_aln_size(data):
    plt.figure(figsize=(5, 4))
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
        plt.scatter(bin_id, bin_precision, s=s, alpha=0.4, c=colors[name], linewidths=0)

    #plt.legend(loc='best', fontsize='xx-small')
    plt.xscale("log")
    plt.xlabel('Alignment size')
    plt.ylabel('Precision')
    plt.ylim(0, 1.1)
    plt.tight_layout()
    plt.savefig(args.output_path + '/size_vs_precision.png', dpi=600)
    #plt.show()
    plt.close()

    handles = [mpl.patches.Patch(color=colors[x], label=x) for x in colors.keys()]
    plt.legend(handles=handles)
    plt.gca().set_axis_off()
    plt.savefig(args.output_path + '/legend.png', dpi=600)


precision_aln_size(data)


def precision_mapq(data):
    plt.figure(figsize=(5, 4))
    for name, df in data.items():
        bin_precision = []
        bin_id = []
        s = []
        for bid, b in df.groupby('mapq'):
            if len(b) < 5:
                continue
            s.append(len(b)*scale)
            bin_precision.append(b['tp'].sum() / (b['tp'].sum() + b['fp'].sum()))
            bin_id.append(bid)

        plt.plot(bin_id, bin_precision, label=name, c=colors[name], alpha=0.8)
        plt.scatter(bin_id, bin_precision, s=s, alpha=0.4, c=colors[name], linewidths=0)

    #plt.legend(loc='best', fontsize='xx-small')
    plt.xlabel('MapQ')
    plt.ylabel('Precision')
    plt.ylim(0, 1.1)
    plt.tight_layout()
    plt.savefig(args.output_path + '/mapq_vs_precision.png', dpi=600)
    #plt.show()
    plt.close()


precision_mapq(data)


def wrong_plot_bins(data):
    plt.figure(figsize=(5, 4))
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
        plt.scatter(bin_w, bin_wrong, s=s, alpha=0.4, c=colors[name], linewidths=0)
    #plt.legend(loc='best', fontsize='xx-small')
    #plt.xscale("log")
    plt.xlabel('Alignment size')
    plt.ylabel('False positive %')
    plt.tight_layout()
    plt.savefig(args.output_path + '/aln_size_vs_fp_bin.png', dpi=600)
    #plt.show()
plt.close()


#wrong_plot_bins(benchmark_res)


def wrong_plot_mapq(data):
    plt.figure(figsize=(6, 4))
    for name, df in data.items():
        bin_wrong = []
        bin_w = []
        wrong = 0
        s=[]
        for bid, b in df.groupby('mapq'):
            bin_wrong.append(wrong / total * 100)
            bin_w.append(bid)
            wrong += b['fp'].sum()
            s.append(len(b) * scale)
        plt.plot(bin_w, bin_wrong, label=name, c=colors[name], alpha=0.8)
        plt.scatter(bin_w, bin_wrong, s=s, alpha=0.4, c=colors[name], linewidths=0)
    #plt.legend(loc='best', fontsize='xx-small')
    #plt.xscale("log")
    plt.xlabel('Mapq')
    plt.ylabel('False positive %')
    plt.tight_layout()
    plt.savefig(args.output_path + '/fp_mapq.png', dpi=600)
    #plt.show()
    plt.close()


#wrong_plot_mapq(benchmark_res)


def wrong_plot_bins2(data):
    plt.figure(figsize=(5, 4))
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
        plt.scatter(bin_w, bin_wrong, s=s, alpha=0.4, c=colors[name], linewidths=0)

    #plt.legend(loc='best', fontsize='xx-small')
    #plt.xscale("log")
    plt.xlabel('Alignment size')
    plt.ylabel('False negative %')
    plt.tight_layout()
    plt.savefig(args.output_path + '/aln_size_vs_fn_bin.png', dpi=600)
    #plt.show()
    plt.close()


#wrong_plot_bins2(benchmark_res)



def precision_recall_read(data):
    plt.figure(figsize=(4, 4))
    data2={}
    for name,df in data.items():
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
        data2 = create_bins(data2, 1, 'mapq')
    for name, df2 in data2.items():
        recall = []
        precision = []
        tp = 0
        fp = 0
        fn = 0
        for i, b in df2.groupby('bins'):
            tp += b['tp'].sum()
            fp += b['fp'].sum()
            fn += b['fn'].sum()
            if tp + fp == 0 or tp + fn == 0:
                continue
            precision.append(tp / (tp + fp))
            recall.append(tp / (tp + fn))
        plt.plot(recall, precision, alpha=1, label=name, c=colors[name], linewidth=1,
                 markeredgecolor=colors[name], marker=markers[name], markerfacecolor="None", markersize=2, markeredgewidth=0.5)
        # plt.gca().invert_xaxis()
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    #plt.legend(loc='best', fontsize='xx-small')
    plt.savefig(args.output_path + '/Precision_Recall_read.png', dpi=600)
    plt.close()


#precision_recall_read(benchmark_res)



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

    fig, ax = plt.subplots(figsize=(4,4))
    ax.bar(data.keys(), counts, label=data.keys(), color=c)
    addlabels(data.keys(), counts)
    ax.set_ylabel('Mapped fragments')
    ax.set_title(title)
    #ax.legend(title='aligners', bbox_to_anchor=(0.9, 0.75))
    plt.setp(ax.get_xticklabels(), rotation=45, horizontalalignment='right')
    plt.subplots_adjust(bottom=0.2)
    plt.tight_layout()
    plt.savefig(args.output_path + '/n_fragments.png', dpi=600)
    plt.close()


mapped_alignments(data)


def fragment_length_dist(data):
    plt.figure(figsize=(4, 4))
    plt.figure()
    for name, df in data.items():
        sns.kdeplot(df['aln_size'], alpha=0.7, label=name, color=colors[name])
    plt.xlabel('fragment length')
    plt.ylabel('density function')
    #plt.legend(fontsize='xx-small')
    plt.savefig(args.output_path + '/fragment_len_dist.png', dpi=600)
    plt.close()


fragment_length_dist(data)


def BWA_curve(data):
    plt.figure(figsize=(5, 4))
    for name, df in data.items():
        x = []
        y = []
        tp = 0
        fp = 0
        for i, b in sorted(df.groupby('mapq'), key=lambda x: x[0], reverse=True):
            tp += b['tp'].sum()
            fp += b['fp'].sum()
            if tp + fp == 0:
                continue
            y.append((fp+tp)/total)
            x.append(fp/(tp+fp))

        plt.plot(x, y, alpha=1, c=colors[name], label=name, linewidth=1.5, markeredgecolor=colors[name],
                 marker=markers[name], markerfacecolor="None", markersize=4, markeredgewidth=0.5)
    plt.ylabel('mapped/total')
    plt.xlabel('wrong/mapped')
    #plt.legend(loc='best', fontsize='xx-small')
    plt.savefig(args.output_path + '/bwamempaper_mapq.png', dpi=600)
    plt.close()


BWA_curve(benchmark_res)


def BWA_curve2(data):
    plt.figure(figsize=(5, 4))
    data2={}
    for name, df in data.items():
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
        data2 = create_bins(data2, 1, 'mapq')
    # BWA-MEM plot
    for name, df in data2.items():
        x = []
        y = []
        s=[]
        tp = df['tp'].sum()
        fp = df['fp'].sum()
        size = tp + fp
        for i, b in df.groupby('bins'):
            if tp+fp == 0:
                continue
            y.append((fp+tp)/total)
            x.append(fp/(tp+fp))
            # s.append(size*scale)
            tp -= b['tp'].sum()
            fp -= b['fp'].sum()
            # size -= tp + fp
        plt.plot(x, y, alpha=1, c=colors[name], label=name, linewidth=1, markeredgecolor=colors[name],
                 marker=markers[name], markerfacecolor="None", markersize=2, markeredgewidth=0.5)

    plt.ylabel('mapped/total')
    plt.xlabel('wrong/mapped')
    #plt.legend(loc='best', fontsize='xx-small')
    plt.savefig(args.output_path + '/bwamempaper_mapq_read.png', dpi=600)
    plt.close()


#BWA_curve2(benchmark_res)


def recall_aln_size(data):
    plt.figure(figsize=(5, 4))
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

        plt.plot(bin_id, bin_recall, label=name, c=colors[name], alpha=1)
        plt.scatter(bin_id, bin_recall, s=s, alpha=0.4, c=colors[name], linewidths=0)

    #plt.legend(loc='best', fontsize='xx-small')
    plt.xscale("log")
    plt.xlabel('Alignment size')
    plt.ylabel('Recall')
    plt.ylim(0, 1.1)
    plt.tight_layout()
    plt.savefig(args.output_path + '/size_vs_recall.png', dpi=600)
    #plt.show()
    plt.close()


recall_aln_size(benchmark_res)


def ROC(data):
    plt.figure(figsize=(4, 4))
    data2={}
    for name,df in data.items():
        mapq = []
        tp = []
        fp = []
        fn = []
        tn = []
        for i, b in df.groupby('qname'):
            mapq.append(b['mapq'].mean())
            tp.append(b['tp'].sum())
            fp.append(b['fp'].sum())
            fn.append(b['fn'].sum())
            tn.append(total-b['tp'].sum()-b['fn'].sum()-b['fp'].sum())
        data2[name] = pd.DataFrame({'mapq':mapq, 'tp':tp, 'fp':fp, 'fn':fn, 'tn':tn})
        data2 = create_bins(data2, 1, 'mapq')
    for name, df2 in data2.items():
        x = []
        y = []
        tp = 0
        fp = 0
        fn = 0
        tn = 0
        for i, b in df2.groupby('bins'):

            tp += b['tp'].sum()
            fp += b['fp'].sum()
            fn += b['fn'].sum()
            tn += b['tn'].sum()
            if tp + fp == 0 or tp + fn == 0:
                continue

            x.append(fp / total)
            y.append(tp / (tp+fn))
        plt.plot(x, y, alpha=1, label=name, c=colors[name], linewidth=0.8, markeredgecolor=colors[name],
                 marker=markers[name], markerfacecolor="None", markersize=3, markeredgewidth=0.5)
    plt.xlabel('False positive / Total number of alignments')
    plt.ylabel('True positive rate')
    #plt.legend(loc='best', fontsize='xx-small')
    plt.savefig(args.output_path + '/ROC.png', dpi=600)
    plt.close()


ROC(benchmark_res)


# precision - number of alignments
def alignments_precision(data):
    plt.figure(figsize=(4, 4))
    data2={}
    for name,df in data.items():
        mapq = []
        tp = []
        fp = []
        fn = []
        aln_diff = []
        for i, b in df.groupby('qname'):
            mapq.append(b['mapq'].mean())
            tp.append(b['tp'].sum())
            fp.append(b['fp'].sum())
            fn.append(b['fn'].sum())
            aln_diff.append(b['n_alignments'].max() - b['n_target'].max())
        data2[name] = pd.DataFrame({'mapq':mapq, 'tp':tp, 'fp':fp, 'fn':fn, 'aln_diff':aln_diff})
    for name, df2 in data2.items():
        x = []
        y = []
        s = []
        for i, b in df2.groupby('aln_diff'):
            tp = b['tp'].sum()
            fp = b['fp'].sum()
            if tp + fp == 0:
                continue
            if len(b) < 5:
                continue
            y.append(tp / (tp + fp))
            x.append(b['aln_diff'].iloc[0])
            s.append(len(b)*0.1)

        plt.plot(x, y, alpha=0.8, label=name, c=colors[name])
        plt.scatter(x, y, s=s, alpha=0.25, linewidths=0, c=colors[name])

    plt.ylabel('Precision')
    plt.xlabel('Mapped alignments - Expected alignments')
    #plt.legend(loc='best', fontsize='xx-small')
    plt.savefig(args.output_path + '/expected_alns_precision.png', dpi=600)
    plt.close()


alignments_precision(benchmark_res)



def alignments_recall(data):
    plt.figure(figsize=(5, 4))
    data2 = {}
    for name, df in data.items():
        mapq = []
        tp = []
        fp = []
        fn = []
        aln_diff = []
        for i, b in df.groupby('qname'):
            mapq.append(b['mapq'].mean())
            tp.append(b['tp'].sum())
            fp.append(b['fp'].sum())
            fn.append(b['fn'].sum())
            aln_diff.append(b['n_alignments'].max() - b['n_target'].max())
        data2[name] = pd.DataFrame({'mapq': mapq, 'tp': tp, 'fp': fp, 'fn': fn, 'aln_diff': aln_diff})
    for name, df2 in data2.items():
        x = []
        y = []
        s = []
        for i, b in df2.groupby('aln_diff'):
            tp = b['tp'].sum()
            fn = b['fn'].sum()
            if tp + fn == 0:
                continue
            if len(b) < 5:
                continue
            y.append(tp / (tp + fn))
            x.append(b['aln_diff'].iloc[0])
            s.append(len(b)*0.1)

        plt.plot(x, y, alpha=0.8, label=name, c=colors[name])
        plt.scatter(x, y, s=s, alpha=0.25, linewidths=0, c=colors[name])

    plt.ylabel('Recall')
    plt.xlabel('Mapped alignments - Expected alignments')
    #plt.legend(loc='best', fontsize='xx-small')
    plt.figure(figsize=(10, 10))
    plt.savefig(args.output_path + '/expected_alns_recall.png', dpi=600)
    plt.close()


alignments_recall(benchmark_res)


def alignments_f_score(data):
    plt.figure(figsize=(5, 4))
    data2 = {}
    for name, df in data.items():
        mapq = []
        tp = []
        fp = []
        fn = []
        aln_diff = []
        for i, b in df.groupby('qname'):
            mapq.append(b['mapq'].mean())
            tp.append(b['tp'].sum())
            fp.append(b['fp'].sum())
            fn.append(b['fn'].sum())
            aln_diff.append(b['n_alignments'].max() - b['n_target'].max())
        data2[name] = pd.DataFrame({'mapq': mapq, 'tp': tp, 'fp': fp, 'fn': fn, 'aln_diff': aln_diff})
    for name, df2 in data2.items():
        x = []
        y = []
        s = []
        for i, b in df2.groupby('aln_diff'):
            tp = b['tp'].sum()
            fn = b['fn'].sum()
            fp = b['fp'].sum()
            if tp + fn == 0:
                continue
            if len(b) < 5:
                continue
            y.append(tp / (tp + 0.5 * (fp + fn)))
            x.append(b['aln_diff'].iloc[0])
            s.append(len(b)*0.1)

        plt.plot(x, y, alpha=0.8, label=name, c=colors[name])
        plt.scatter(x, y, s=s, alpha=0.25, linewidths=0, c=colors[name])

    plt.ylabel('F-score')
    plt.xlabel('Mapped alignments - Expected')
    #plt.legend(loc='best', fontsize='xx-small')
    plt.savefig(args.output_path + '/expected_alns_f_score.png', dpi=600)
    plt.close()


alignments_f_score(benchmark_res)