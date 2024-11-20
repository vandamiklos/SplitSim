# Find simpler SVs in the generated events
import argparse
import pysam
import pandas as pd

parse = argparse.ArgumentParser()
parse.add_argument('--input_path', help='')
args = parse.parse_args()

file = args.input_path + '/simulated.fq'

pattern = 'alignments_2'
alignments_2 = []
with pysam.FastxFile(file) as fastq:
    for entry in fastq:
        # Check if the pattern is in the name, sequence, or quality
        if pattern in entry.name:
            name = entry.name.replace('__', ' ')
            parts=name.split(' ')
            qname = parts[0]
            info=parts[2].split('_')
            part_1=info[0].split(':')
            chr_1=part_1[0]
            start_1, end_1=part_1[1].split('-')
            part_2=info[1].split(':')
            chr_2=part_2[0]
            start_2, end_2=part_2[1].split('-')

            alignments_2.append({"qname": qname, "chr_1": chr_1, "start_1":start_1, "end_1":end_1, "chr_2": chr_2, "start_2":start_2, "end_2":end_2, "seq": entry.sequence})

alignments_2_df = pd.DataFrame(alignments_2)

duplications = alignments_2_df[alignments_2_df["chr_1"] == alignments_2_df["chr_2"]]
duplications = duplications[duplications["end_1"] > duplications["start_1"]]