import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import argparse
import numpy as np
import seaborn as sns

parse = argparse.ArgumentParser()
parse.add_argument('--input_path_ont', help='')
parse.add_argument('--input_path_pacbio', help='')
parse.add_argument('--output_path', help='')
parse.add_argument('--aligner_name', help='List of aligner names',
    nargs='+',  # Accept one or more arguments
    type=str)
args = parse.parse_args()

data = {}
benchmark_res={}



for aligner in args.aligner_name:
    data[aligner + '_ont'] = pd.read_csv(args.input_path_ont + '/' + aligner + '.mappings_labelled.csv', sep='\t')
    benchmark_res[aligner + '_ont'] = pd.read_csv(args.input_path_ont + '/' + aligner + '.benchmark_res_fn.csv', sep='\t')
    stats_ont = pd.read_csv(args.input_path_ont + '/' + args.aligner_name[0] + '.stats.txt', sep='\t')
    total_ont = int(stats_ont['target_n'].iloc[0])

    data[aligner + '_pacbio'] = pd.read_csv(args.input_path_pacbio + '/' + aligner + '.mappings_labelled.csv', sep='\t')
    benchmark_res[aligner + '_pacbio'] = pd.read_csv(args.input_path_pacbio + '/' + aligner + '.benchmark_res_fn.csv', sep='\t')
    stats_pacbio = pd.read_csv(args.input_path_pacbio + '/' + args.aligner_name[0] + '.stats.txt', sep='\t')
    total_pacbio = int(stats_pacbio['target_n'].iloc[0])

#colors = {'bwa': '#F0A430FF', 'bwa_dodi': 'firebrick',
#          'bwa_flags': 'gold', 'bwa_flags_dodi': 'lightcoral',
#          'bwa_x': 'red',
#          'minimap2': '#800000FF', 'minimap2_dodi': 'darkgreen',
#          'lastalsplit': '#1B3A54FF', 'lastal_dodi':'lightblue',
#          'ngmlr': '#205F4BFF', 'ngmlr_x': 'chocolate',
#          'vacmap_r': '#774762FF', 'vacmap_s': '#774762FF'}

# #90728FFF, #B9A0B4FF, #9D983DFF, #CECB76FF, #E15759FF, #FF9888FF, #6B6B6BFF, #BAB2AEFF, #AA8780FF, #DAB6AFFF
colors = {'bwa_ont': '#665065FF',
          'minimap2_ont': '#666328FF',
          'lastalsplit_ont': '#c82426FF',
          'ngmlr_ont': '#303848FF',
          'vacmap_s_ont': '#664a45FF',
          'bwa_pacbio': '#B9A0B4FF',
          'minimap2_pacbio': '#CECB76FF',
          'lastalsplit_pacbio': '#FF9888FF',
          'ngmlr_pacbio': '#5a6884FF',
          'vacmap_s_pacbio': '#946b63FF'
          }


markers={'bwa_ont': 'o',
          'minimap2_ont': 'x',
          'lastalsplit_ont': '+',
          'ngmlr_ont': 'p',
          'vacmap_s_ont': 's',
          'bwa_pacbio': 'o',
          'minimap2_pacbio': 'x',
          'lastalsplit_pacbio': '+',
          'ngmlr_pacbio': 'p',
          'vacmap_s_pacbio': 's'
          }

scale = 0.01


def create_bins(dict, base, column):
    for name, df in dict.items():
        bins = []
        for i in df[column]:
            bins.append(base * round(i/base))
        dict[name] = df.assign(bins=bins)
    return dict


data = create_bins(data, 25, 'aln_size')
benchmark_res = create_bins(benchmark_res, 25, 'aln_size')


def precision_aln_size_log(data):
    plt.figure(figsize=(5, 4))
    for name, df in data.items():
        bin_precision = []
        bin_id = []
        s = []
        for bid, b in df.groupby('bins'):
            if (b['tp'].sum() + b['fp'].sum()) == 0:
                continue
            s.append(len(b)*scale)
            bin_precision.append(b['tp'].sum() / (b['tp'].sum() + b['fp'].sum()))
            bin_id.append(bid)

        plt.plot(bin_id, bin_precision, label=name, c=colors[name], alpha=0.8)
        plt.scatter(bin_id, bin_precision, s=s, alpha=0.4, c=colors[name], linewidths=0)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)

    #plt.legend(loc='best', fontsize='xx-small')
    plt.xscale("log")
    plt.grid()
    plt.xticks([25, 50, 100, 150, 300, 500, 1000], [25, 50, 100, 150, 300, 500, 1000])
    plt.xlabel('Alignment size', fontsize=13, weight='bold')
    plt.ylabel('Precision', fontsize=13, weight='bold')
    plt.ylim(0, 1.1)
    plt.xlim(0, 1000)
    plt.tight_layout()
    plt.savefig(args.output_path + '/size_vs_precision_log.png', dpi=600)
    #plt.show()
    plt.close()

    handles = [mpl.patches.Patch(color=colors[x], label=x) for x in colors.keys()]
    plt.legend(handles=handles)
    plt.gca().set_axis_off()
    plt.savefig(args.output_path + '/legend.png', dpi=600)


#precision_aln_size_log(data)


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
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.locator_params(axis='x', nbins=10)

    #plt.legend(loc='best', fontsize='xx-small')
    plt.xlabel('MapQ', fontsize=13, weight='bold')
    plt.ylabel('Precision', fontsize=13, weight='bold')
    plt.grid()
    plt.ylim(0, 1.1)
    plt.tight_layout()
    plt.savefig(args.output_path + '/mapq_vs_precision.png', dpi=600)
    #plt.show()
    plt.close()


#precision_mapq(data)

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
    #ax.legend(title='aligners', bbox_to_anchor=(0.9, 0.75))
    plt.setp(ax.get_xticklabels(), rotation=45, horizontalalignment='right')
    plt.subplots_adjust(bottom=0.2)
    plt.tight_layout()
    plt.savefig(args.output_path + '/n_fragments.png', dpi=600)
    plt.close()


#mapped_alignments(data)


def fragment_length_dist(data):
    plt.figure(figsize=(4, 4))
    plt.figure()
    for name, df in data.items():
        sns.kdeplot(df['aln_size'], alpha=0.7, label=name, color=colors[name])
    plt.xlabel('fragment length', fontsize=13, weight='bold')
    plt.ylabel('density function', fontsize=13, weight='bold')
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize='xx-small')
    plt.savefig(args.output_path + '/fragment_len_dist.png', dpi=600)
    plt.close()


#fragment_length_dist(data)


def BWA_curve(data):
    plt.figure(figsize=(5, 4))
    for name, df in data.items():
        if 'ont' in name:
            total = total_ont
        else:
            total = total_pacbio
        x = []
        y = []
        tp = 0
        fp = 0
        for i, b in sorted(df.groupby('mapq'), key=lambda x: x[0], reverse=True):
            tp += b['tp'].sum()
            fp += b['fp'].sum()
            if tp + fp == 0:
                continue
            y.append((fp+tp)/(total))
            x.append(fp/(tp+fp))

        plt.plot(x, y, alpha=1, c=colors[name], label=name, linewidth=1.5, markeredgecolor=colors[name],
                 marker=markers[name], markerfacecolor="None", markersize=6, markeredgewidth=0.5)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)


    plt.ylabel('mapped/total', fontsize=13, weight='bold')
    plt.xlabel('wrong/mapped', fontsize=13, weight='bold')
    plt.grid()
    plt.locator_params(axis='x', nbins=8)
    plt.tight_layout()
    #plt.legend(loc='best', fontsize='xx-small')
    plt.savefig(args.output_path + '/bwamempaper_mapq.png', dpi=600)
    plt.close()


#BWA_curve(benchmark_res)


def recall_aln_size_log(data):
    plt.figure(figsize=(5, 4))
    for name, df in data.items():
        bin_recall = []
        bin_id = []
        s = []
        for bid, b in df.groupby('bins'):
            if (b['tp'].sum() + b['fn'].sum()) == 0:
                continue
            s.append(len(b)*scale)
            bin_recall.append(b['tp'].sum() / (b['tp'].sum() + b['fn'].sum()))
            bin_id.append(bid)

        plt.plot(bin_id, bin_recall, label=name, c=colors[name], alpha=1)
        plt.scatter(bin_id, bin_recall, s=s, alpha=0.4, c=colors[name], linewidths=0)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)

    #plt.legend(loc='best', fontsize='xx-small')
    plt.xscale("log")
    plt.xticks([25, 50, 100, 150, 300, 500, 1000], [25, 50, 100, 150, 300, 500, 1000])
    plt.xlabel('Alignment size', fontsize=13, weight='bold')
    plt.ylabel('Recall', fontsize=13, weight='bold')
    plt.grid()
    plt.ylim(0, 1.1)
    plt.xlim(0, 1000)
    plt.tight_layout()
    plt.savefig(args.output_path + '/size_vs_recall_log.png', dpi=600)
    #plt.show()
    plt.close()


#recall_aln_size_log(benchmark_res)


# precision - number of alignments
def alignments_precision(data):
    plt.figure(figsize=(5, 4))
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
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)

    plt.ylabel('Precision', fontsize=13, weight='bold')
    plt.xlabel('Mapped alignments - Expected alignments', fontsize=13, weight='bold')
    plt.grid()
    plt.tight_layout()
    #plt.legend(loc='best', fontsize='xx-small')
    plt.savefig(args.output_path + '/expected_alns_precision.png', dpi=600)
    plt.close()


#(benchmark_res)



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
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)

    plt.ylabel('Recall', fontsize=16, weight='bold')
    plt.xlabel('Mapped alignments - Expected alignments', fontsize=16, weight='bold')
    plt.grid()
    #plt.legend(loc='best', fontsize='xx-small')
    plt.tight_layout()
    plt.savefig(args.output_path + '/expected_alns_recall.png', dpi=600)
    plt.close()


#alignments_recall(benchmark_res)


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
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)

    plt.ylabel('F-score', fontsize=16, weight='bold')
    plt.xlabel('Mapped alignments - Expected', fontsize=16, weight='bold')
    plt.grid()
    plt.tight_layout()
    #plt.legend(loc='best', fontsize='xx-small')
    plt.savefig(args.output_path + '/expected_alns_f_score.png', dpi=600)
    plt.close()


#alignments_f_score(benchmark_res)


def mapq_aln_size_log(data):
    plt.figure(figsize=(5, 4))
    for name, df in data.items():
        bin_precision = []
        bin_id = []
        s = []
        for bid, b in df.groupby('bins'):
            if (b['mapq'].sum()) == 0:
                continue
            s.append(len(b)*scale)
            bin_precision.append(b['mapq'].sum()/len(b))
            bin_id.append(bid)

        plt.plot(bin_id, bin_precision, label=name, c=colors[name], alpha=0.8)
        plt.scatter(bin_id, bin_precision, s=s, alpha=0.4, c=colors[name], linewidths=0)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)

    #plt.legend(loc='best', fontsize='xx-small')
    plt.xscale("log")
    plt.xticks([25, 50, 100, 150, 300, 500, 1000], [25, 50, 100, 150, 300, 500, 1000])
    plt.xlabel('Alignment size', fontsize=13, weight='bold')
    plt.ylabel('Mapping quality', fontsize=13, weight='bold')
    plt.xlim(0, 1000)
    plt.grid()
    plt.tight_layout()
    plt.savefig(args.output_path + '/size_vs_mapq_log.png', dpi=600)

    #plt.show()
    plt.close()

    handles = [mpl.patches.Patch(color=colors[x], label=x) for x in colors.keys()]
    plt.legend(handles=handles)
    plt.gca().set_axis_off()
    plt.savefig(args.output_path + '/legend.png', dpi=600)


#mapq_aln_size_log(data)



def precision_recall_frag(data):
    plt.figure(figsize=(4, 4))
    for name, df in data.items():
        bin_precision = []
        bin_recall = []
        bin_id = []
        s = []
        for bid, b in df.groupby('bins'):
            if (b['tp'].sum() + b['fp'].sum()) == 0 or (b['tp'].sum() + b['fn'].sum()) == 0:
                continue
            if len(b) < 200:
                continue
            s.append(len(b)*scale)
            bin_precision.append(b['tp'].sum() / (b['tp'].sum() + b['fp'].sum()))
            bin_recall.append(b['tp'].sum() / (b['tp'].sum() + b['fn'].sum()))
            bin_id.append(bid)
        plt.plot(bin_recall, bin_precision, alpha=1, label=name, c=colors[name], linewidth=1,
                 markeredgecolor=colors[name], marker=markers[name], markerfacecolor="None", markersize=2, markeredgewidth=0.5)
    plt.xlabel('Recall', fontsize=13, weight='bold')
    plt.ylabel('Precision', fontsize=13, weight='bold')
    plt.grid()
    plt.tight_layout()
    plt.savefig(args.output_path + '/Precision_Recall.png', dpi=600)
    plt.close()


benchmark_res2 = create_bins(benchmark_res, 50, 'aln_size')


precision_recall_frag(benchmark_res2)