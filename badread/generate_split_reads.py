import sys
import random
from badread import misc, fragment_lengths
import pysam
from scipy.stats import poisson


def read_fasta(args):
    fasta = {}
    fin = pysam.FastxFile(args.fasta)
    for line in fin:
        fasta[line.name] = line.sequence
    if len(fasta) == 0:
        raise ValueError("Empty fasta file")
    return fasta


def generate_split_reads(args, ref, n_seqs, mean, frag_lengths):
    # names of the reference sequences
    chroms = list(ref.references)
    strand = ['forward', 'reverse']
    #generate split-read sequences and corresponding names
    for n in range(n_seqs):
        if mean == 0:
            raise ValueError("mean must be > 0")
        blocks = 0
        # choose the number of fragments that a read is split into
        # make sure it is split into more than 0 fragments
        while not blocks:
            blocks = poisson.rvs(mu=mean, size=1)[0]
            if not blocks:
                continue
        ins_seqs = []
        # initialize the name of the read with the number of blocks
        names = [f">alignments_{blocks}_"]
        blk = 0
        # generate fragments until the number of fragments reaches the value stored in blocks
        while blk < blocks:
            # get the length of the block
            # check if the length is longer than 15 bases
            flen = frag_lengths.get_fragment_length()
            if flen < 15:
                continue
            # randomly choose a chromosome
            c = random.choice(chroms)
            # randomly choose a strand
            s = random.choice(strand)
            #randomly choose the start coordinate of the fragment
            # make sure the start position is not at the end of the chromosome
            pos = random.randint(1, ref.get_reference_length(c) - flen)
            if pos + flen > ref.get_reference_length(c):
                continue  # happens rarely
            blk += 1
            # fetch the sequence based on the chosen chromosome, start and end coordinates
            seq = ref.fetch(c, pos, pos + flen).upper()
            # reverse complement the sequence if 'reverse' strand was chosen
            if s == 'reverse':
                seq = misc.reverse_complement(seq)
            ins_seqs.append(seq)
            # add the fragment info the the name
            names.append(f"{c}:{pos}-{pos+flen}")
        final_seq = "".join(ins_seqs)
        final_name = "_".join(names)
        print(final_name)
        print(final_seq)


def generate_reads(args):
    ref = pysam.FastaFile(args.reference)

    print(f"Generating {args.number} split-reads", file=sys.stderr)

    frag_lengths = fragment_lengths.FragmentLengths(args.mean_block_len, args.std_block_len)
    generate_split_reads(args, ref, args.number, args.mean, frag_lengths)

    print(f"Done", file=sys.stderr)