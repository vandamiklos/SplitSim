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


def generate_split_reads(args, ref, n_seqs, mean, frag_lengths):
    chroms = list(ref.references)
    strand = ['forward', 'reverse']
    for n in range(n_seqs):
        if mean == 0:
            raise ValueError("mean must be > 0")
        blocks = 0
        while not blocks:
            blocks = poisson.rvs(mu=mean, size=1)[0]
            if not blocks:
                continue
        ins_seqs = []
        names = [f">split_{blocks}__"]
        blk = 0
        while blk < blocks:
            flen = frag_lengths.get_fragment_length()
            if flen < 15:
                continue
            c = random.choice(chroms)
            s = random.choice(strand)
            pos = random.randint(1, ref.get_reference_length(c) - flen)
            if pos + flen > ref.get_reference_length(c):
                continue  # happens rarely
            blk += 1
            seq = ref.fetch(c, pos, pos + flen).upper()
            if s == 'reverse':
                seq = pysam.reverse_complement(seq)
            ins_seqs.append(seq)
            names.append(f"{c}:{pos}-{pos+flen}")
        final_seq = "".join(ins_seqs)
        final_name = "_".join(names)
        print(final_name)
        print(final_seq)

@click.command()
@click.argument('reference')
@click.argument('number', type=int)
@click.option("--mean", help="alignment number mean (poisson distribution)", default=3, show_default=True, type=int)
@click.option("--mean-block-len", help="alignment length mean (gamma distribution)", default=150, show_default=True, type=int)
@click.option("--std-block-len", help="alignment length stdev (gamma distribution)", default=150, show_default=True, type=int)
@click.version_option(__version__)
def generate_reads(**args):
    ref = pysam.FastaFile(args['reference'])
    print(f"Generating {args['number']} split-reads", file=stderr)
    frag_lengths = fragment_lengths.FragmentLengths(args['mean_block_len'], args['std_block_len'])
    generate_split_reads(args, ref, args['number'], args['mean'], frag_lengths)
    print(f"Done", file=stderr)


if __name__ == "__main__":
    generate_reads()