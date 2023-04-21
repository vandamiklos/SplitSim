"""
Generate fusion-seq amplicons

Inputs:
-------
fasta file of the telomere ends (include primer seq)
N fusion seqs

Outputs:
--------
fasta file of telomere fusions with multiple insertions written to stdout

Notes
-----
fasta sequence must have telomere end on right-hand-side. Make sure to reverse complement reverse strand sequence

"""

import click
import random
from badread import misc, fragment_lengths
import pysam
from scipy.stats import poisson
from sys import stderr


__version__ = '0.1'


def read_fasta(args):
    fasta = {}
    fin = pysam.FastxFile(args['fasta'])
    for line in fin:
        fasta[line.name] = line.sequence
    if len(fasta) == 0:
        raise ValueError("Empty fasta file")
    return fasta


def generate_insertions(args, ref, fasta, n_seqs, frag_lengths):
    rev_comps = {k: misc.reverse_complement(a) for k, a in fasta.items()}
    chroms = list(ref.references)
    # make head-to-head fusion with some more fragments in between
    keys = list(fasta.keys())
    for n in range(n_seqs):
        t = random.choice(keys)
        a = fasta[t].upper()
        b = rev_comps[t]
        if args["mu_ins"] == 0:
            raise ValueError("mu-ins must be > 0")
        blocks = 0
        while not blocks:
            blocks = poisson.rvs(mu=3, size=1)[0]
            if not blocks:
                continue

        # print("Insertion number", blocks, file=stderr)
        ins_seqs = []
        tname = f"{t}:0-{len(a)}"
        names = [f">insertion_{blocks}__" + tname]
        blk = 0
        while blk < blocks:
            flen = frag_lengths.get_fragment_length()
            if flen < 15:
                continue
            c = random.choice(chroms)
            pos = random.randint(1, ref.get_reference_length(c) - flen)
            if pos + flen > ref.get_reference_length(c):
                continue  # happens rarely
            blk += 1
            # print("    len", flen, file=stderr)
            ins_seqs.append(ref.fetch(c, pos, pos + flen).upper())
            names.append(f"{c}:{pos}-{pos+flen}")
        final_seq = a + "".join(ins_seqs) + b
        names.append(tname)
        final_name = "_".join(names)
        print(final_name)
        print(final_seq)


def generate_deletions(fasta, n_seqs, frag_lengths):
    rev_comps = {k: misc.reverse_complement(a) for k, a in fasta.items()}
    # make head-to-head fusion with one fragment with a deletion
    keys = list(fasta.keys())
    for n in range(n_seqs):
        t = random.choice(keys)
        a = fasta[t].upper()
        b = rev_comps[t]
        left_hand_side = random.random() > 0.5
        l = frag_lengths.get_fragment_length()
        if left_hand_side:
            flen = max(15, len(a) - l)
            print(f">deletion_{t}:0-{flen}_{t}:0-{len(b)}")
            print(a[0:flen] + b)
        else:
            flen = max(15, len(b) - frag_lengths.get_fragment_length())
            print(f">deletion_{t}:0-{len(a)}_{t}:0-{flen}")
            print(a + b[flen:])
        print(f"deletion {'left' if left_hand_side else 'right'} {l}", file=stderr)


@click.command()
@click.argument('reference')
@click.argument('fasta')
@click.argument('number', type=int)
@click.option("--mu-ins", help="insertion number mean (poisson distribution)", default=3, show_default=True, type=int)
@click.option("--mean-block-len", help="insertion length mean (gamma distribution)", default=150, show_default=True, type=int)
@click.option("--std-block-len", help="insertion length stdev (gamma distribution)", default=150, show_default=True, type=int)
@click.version_option(__version__)
def generate_fusions(**args):
    ref = pysam.FastaFile(args['reference'])
    fasta = read_fasta(args)
    print(f"Generating {args['number']} fusions with insertions", file=stderr)
    frag_lengths = fragment_lengths.FragmentLengths(args['mean_block_len'], args['std_block_len'])
    generate_insertions(args, ref, fasta, args['number'], frag_lengths)
    # generate_deletions(fasta, args['number'], frag_lengths)
    print(f"Done", file=stderr)


if __name__ == "__main__":
    generate_fusions()
