import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import argparse
import numpy as np
import seaborn as sns


parse = argparse.ArgumentParser()
parse.add_argument('--input', help='')
parse.add_argument('--output_path', help='')
args = parse.parse_args()

summary = pd.read_csv(args.input, sep='\t')


def precision_recall_frag(data):
    plt.figure(figsize=(4, 4))
    plt.plot(recall, precision, alpha=1, linewidth=1, markerfacecolor="None", markersize=2, markeredgewidth=0.5)
    plt.xlabel('Recall', fontsize=13, weight='bold')
    plt.ylabel('Precision', fontsize=13, weight='bold')
    plt.grid()
    plt.tight_layout()
    plt.savefig(args.output_path + '/Precision_Recall.png', dpi=600)
    plt.close()


precision_recall_frag(benchmark_res2)