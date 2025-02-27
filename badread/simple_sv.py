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


def generate_duplication(args, ref, n_seqs, frag_lengths, valid_chroms, fix_overlap):
    chroms = valid_chroms
    strand = ['forward', 'reverse']
    for n in range(n_seqs):
        seqs = []
        names = [f">duplication_"]
        c = random.choice(chroms)
        s = random.choice(strand)

        flen = []
        blk = 0
        while blk < 2:
            f = frag_lengths.get_fragment_length()
            if f < 15:
                continue
            if ref.get_reference_length(c) - 3 * f < 0:
                c = random.choice(chroms)
                continue
            if blk == 0:
                pos = random.randint(f, ref.get_reference_length(c) - 3*f)
            blk += 1
            flen.append(f)

        seq1 = ref.fetch(c, pos, pos + flen[0]).upper()
        names.append(f"{c}:{pos}-{pos + flen[0]}")
        overlap = random.uniform(fix_overlap, 1)

        if flen[0] > flen[1]:
            seq2 = ref.fetch(c, pos + flen[0] - int(overlap * flen[1]), pos + flen[0] - int(overlap * flen[1]) + flen[1]).upper()
            names.append(f"{c}:{pos + flen[0] - int(overlap * flen[1])}-{pos + flen[0] - int(overlap * flen[1]) + flen[1]}")
        if flen[0] <= flen[1]:
            seq2 = ref.fetch(c, pos + flen[0] - int(overlap * flen[0]), pos + flen[0] - int(overlap * flen[0]) + flen[1]).upper()
            names.append(f"{c}:{pos + flen[0] - int(overlap * flen[0])}-{pos + flen[0] - int(overlap * flen[0]) + flen[1]}")

        if s == 'reverse':
            seq1 = misc.reverse_complement(seq1)
            seq2 = misc.reverse_complement(seq2)

        seqs.append(seq1)
        seqs.append(seq2)

        final_seq = "".join(seqs)
        final_name = "_".join(names)
        print(final_name)
        print(final_seq)

def generate_translocation(args, ref, n_seqs, frag_lengths):
    chroms = list(ref.references)
    strand = ['forward', 'reverse']
    for n in range(n_seqs):
        blocks = 2
        ins_seqs = []
        names = [f">translocation_"]
        blk = 0
        while blk < blocks:
            c = random.choice(chroms)
            s = random.choice(strand)
            flen = frag_lengths.get_fragment_length()
            if flen < 15:
                continue
            pos = random.randint(1, ref.get_reference_length(c) - flen)
            if pos + flen > ref.get_reference_length(c):
                continue  # happens rarely
            blk += 1
            seq = ref.fetch(c, pos, pos + flen).upper()
            if s == 'reverse':
                seq = misc.reverse_complement(seq)
            ins_seqs.append(seq)
            names.append(f"{c}:{pos}-{pos+flen}")
        final_seq = "".join(ins_seqs)
        final_name = "_".join(names)
        print(final_name)
        print(final_seq)


def generate_inversion3(args, ref, n_seqs, frag_lengths, valid_chroms):
    chroms = valid_chroms
    strand = ['forward', 'reverse']
    for n in range(n_seqs):
        c = random.choice(chroms)
        s = random.choice(strand)
        flen = []
        blk = 0
        while blk < 3:
            f = frag_lengths.get_fragment_length()
            if f < 15:
                continue
            if blk == 0:
                if ref.get_reference_length(c) - 3*f < 0:
                    c = random.choice(chroms)
                    continue
                pos = random.randint(1, ref.get_reference_length(c) - 3*f)
            blk += 1
            flen.append(f)

        seq3 = ref.fetch(c, pos + flen[0] + flen[1], pos + flen[0] + flen[1] +flen[2]).upper()
        seq2 = ref.fetch(c, pos + flen[0], pos + flen[0] + flen[1]).upper()
        seq1 = ref.fetch(c, pos, pos + flen[0]).upper()

        if s == 'reverse':
            seq1 = misc.reverse_complement(seq1)
            seq3 = misc.reverse_complement(seq3)
        else:
            seq2 = misc.reverse_complement(seq2)

        ins_seqs = [seq1, seq2, seq3]
        names = [f">inversion3_", f"{c}:{pos}-{pos+flen[0]}",
                 f"{c}:{pos + flen[0]}-{pos + flen[0] + flen[1]}",
                 f"{c}:{pos + flen[0] + flen[1]}-{pos + flen[0] + flen[1] +flen[2]}"]
        final_seq = "".join(ins_seqs)
        final_name = "_".join(names)
        print(final_name)
        print(final_seq)


def generate_inversion2(args, ref, n_seqs, frag_lengths, valid_chroms):
    chroms = valid_chroms
    strand = ['forward', 'reverse']
    for n in range(n_seqs):
        c = random.choice(chroms)
        s = random.choice(strand)
        flen = []
        blk = 0
        while blk < 2:
            f = frag_lengths.get_fragment_length()
            if f < 15:
                continue
            if blk == 0:
                if ref.get_reference_length(c) - 3*f < 0:
                    c = random.choice(chroms)
                    continue
                pos = random.randint(1, ref.get_reference_length(c) - 3*f)
            blk += 1
            flen.append(f)

        seq2 = ref.fetch(c, pos + flen[0], pos + flen[0] + flen[1]).upper()
        seq1 = ref.fetch(c, pos, pos + flen[0]).upper()

        if s == 'reverse':
            seq1 = misc.reverse_complement(seq1)
        else:
            seq2 = misc.reverse_complement(seq2)

        ins_seqs = [seq1, seq2]
        names = [f">inversion2_", f"{c}:{pos}-{pos+flen[0]}",
                 f"{c}:{pos + flen[0]}-{pos + flen[0] + flen[1]}"]
        final_seq = "".join(ins_seqs)
        final_name = "_".join(names)
        print(final_name)
        print(final_seq)


def generate_deletion(args, ref, n_seqs, frag_lengths, valid_chroms):
    chroms = valid_chroms
    strand = ['forward', 'reverse']
    for n in range(n_seqs):
        c = random.choice(chroms)
        s = random.choice(strand)
        flen = []
        blk = 0
        while blk < 3:
            f = frag_lengths.get_fragment_length()
            if f < 15:
                continue
            if blk == 0:
                if ref.get_reference_length(c) - 3*f < 0:
                    c = random.choice(chroms)
                    continue
                pos = random.randint(1, ref.get_reference_length(c) - 3*f)
            blk += 1
            flen.append(f)

        seq1 = ref.fetch(c, pos, pos + flen[0]).upper()
        seq3 = ref.fetch(c, pos + flen[0] + flen[1],  pos + flen[0] + flen[1] + flen[2]).upper()

        if s == 'reverse':
            seq1 = misc.reverse_complement(seq1)
            seq3 = misc.reverse_complement(seq3)

        seqs = [seq1, seq3]
        names = [f">deletion_", f"{c}:{pos}-{pos + flen[0]}",
                 f"{c}:{pos + flen[0] + flen[1]}-{pos + flen[0] + flen[1] + flen[2]}"]
        final_seq = "".join(seqs)
        final_name = "_".join(names)
        print(final_name)
        print(final_seq)


def generate_insertion(args, ref, n_seqs, frag_lengths, valid_chroms):
    chroms = valid_chroms
    strand = ['forward', 'reverse']
    for n in range(n_seqs):
        # randomly choose a chromosome/contig
        c = random.choice(chroms)
        # randomly choose a strand
        s = random.choice(strand)
        flen = []
        pos = []
        blk = 0
        while blk < 3:
            # choose a fragment length
            f = frag_lengths.get_fragment_length()
            # if it's too short choose again
            if f < 15:
                continue
            # use the same chromosome/contig for each fragment in the read
            # choose the position of the first fragment in the read
            # make sure the first fragment is not at the end of the chromosome/contig
            if blk == 0:
                if ref.get_reference_length(c) - 3*f < 0:
                    c = random.choice(chroms)
                    continue
                p = random.randint(1, ref.get_reference_length(c) - 3*f)
            if blk ==1:
                p = random.randint(1, ref.get_reference_length(c) - 3 * f)
            blk += 1
            flen.append(f)
            pos.append(p)
        # fetch the sequences from the reference based on the chosen chromosome/contig and position
        seq1 = ref.fetch(c, pos[0], pos[0] + flen[0]).upper()
        seq2 = ref.fetch(c, pos[1], pos[1] + flen[0]).upper()
        seq3 = ref.fetch(c, pos[0] + flen[0], pos[0] + flen[0] + flen[2]).upper()

        if s == 'reverse':
            seq1 = misc.reverse_complement(seq1)
            seq2 = misc.reverse_complement(seq2)
            seq3 = misc.reverse_complement(seq3)


        seqs = [seq1, seq2, seq3]
        names = [f">insertion_", f"{c}:{pos[0]}-{pos[0]+flen[0]}",
                 f"{c}:{pos[1]}-{pos[1] + flen[1]}",
                 f"{c}:{pos[0] + flen[0]}-{pos[0] + flen[0] + flen[2]}"]
        final_seq = "".join(seqs)
        final_name = "_".join(names)
        print(final_name)
        print(final_seq)


def DNA(length):
    DNA=''
    for l in range(length):
        DNA += random.choice('CGTA')
    return DNA


def generate_random_insertion(args, ref, n_seqs, frag_lengths, valid_chroms):
    chroms = valid_chroms
    strand = ['forward', 'reverse']
    for n in range(n_seqs):
        c = random.choice(chroms)
        s = random.choice(strand)
        flen = []
        blk = 0
        while blk < 3:
            f = frag_lengths.get_fragment_length()
            if f < 15:
                continue
            if blk == 0:
                if ref.get_reference_length(c) - 3*f < 0:
                    c = random.choice(chroms)
                    continue
                pos = random.randint(1, ref.get_reference_length(c) - 3*f)
            blk += 1
            flen.append(f)

        seq1 = ref.fetch(c, pos, pos + flen[0]).upper()
        seq2 = DNA(flen[1])
        seq3 = ref.fetch(c, pos + flen[0], pos + flen[0] + flen[2]).upper()

        if s == 'reverse':
            seq1 = misc.reverse_complement(seq1)
            seq3 = misc.reverse_complement(seq3)

        seqs = [seq1, seq2, seq3]
        names = [f">randominsertion_", f"{c}:{pos}-{pos+flen[0]}",
                 f"randomchr:0-0",
                 f"{c}:{pos + flen[0] + flen[1]}-{pos + flen[0] + flen[1] + flen[2]}"]
        final_seq = "".join(seqs)
        final_name = "_".join(names)
        print(final_name)
        print(final_seq)


def generate_N_insertion(args, ref, n_seqs, frag_lengths, valid_chroms):
    chroms = valid_chroms
    strand = ['forward', 'reverse']
    for n in range(n_seqs):
        c = random.choice(chroms)
        s = random.choice(strand)
        flen = []
        blk = 0
        while blk < 3:
            f = frag_lengths.get_fragment_length()
            if f < 15:
                continue
            if blk == 0:
                if ref.get_reference_length(c) - 3*f < 0:
                    c = random.choice(chroms)
                    continue
                pos = random.randint(1, ref.get_reference_length(c) - 3*f)
            blk += 1
            flen.append(f)

        seq1 = ref.fetch(c, pos, pos + flen[0]).upper()
        seq2 = flen[1] * 'N'
        seq3 = ref.fetch(c, pos + flen[0], pos + flen[0] + flen[2]).upper()

        if s == 'reverse':
            seq1 = misc.reverse_complement(seq1)
            seq3 = misc.reverse_complement(seq3)

        seqs = [seq1, seq2, seq3]
        names = [f">ninsertion_", f"{c}:{pos}-{pos+flen[0]}",
                 f"N:0-0",
                 f"{c}:{pos + flen[0] + flen[1]}-{pos + flen[0] + flen[1] + flen[2]}"]
        final_seq = "".join(seqs)
        final_name = "_".join(names)
        print(final_name)
        print(final_seq)



def generate_svs(args):
    ref = pysam.FastaFile(args.reference)

    print(f"Generating {args.number} SVs", file=sys.stderr)

    frag_lengths = fragment_lengths.FragmentLengths(args.mean_block_len, args.std_block_len)

    sample = [frag_lengths.get_fragment_length() for _ in range(10000)]
    max_len = max(sample)
    chroms = list(ref.references)
    valid_chroms = []
    for c in chroms:
        if 3 * max_len < ref.get_reference_length(c):
            valid_chroms.append(c)

    generate_duplication(args, ref, args.number, frag_lengths, valid_chroms, args.fix_overlap)
    generate_deletion(args, ref, args.number, frag_lengths, valid_chroms)
    generate_random_insertion(args, ref, args.number, frag_lengths, valid_chroms)
    generate_N_insertion(args, ref, args.number, frag_lengths, valid_chroms)
    generate_insertion(args, ref, args.number, frag_lengths, valid_chroms)
    generate_inversion2(args, ref, args.number, frag_lengths, valid_chroms)
    generate_inversion3(args, ref, args.number, frag_lengths, valid_chroms)
    generate_translocation(args, ref, args.number, frag_lengths)

    print(f"Done", file=sys.stderr)