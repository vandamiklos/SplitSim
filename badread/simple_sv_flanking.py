import sys
import random
from badread import misc, fragment_lengths
import pysam


def read_fasta(args):
    fasta = {}
    fin = pysam.FastxFile(args.fasta)
    for line in fin:
        fasta[line.name] = line.sequence
    if len(fasta) == 0:
        raise ValueError("Empty fasta file")
    return fasta


def generate_duplication(args, ref, n_seqs, frag_lengths, valid_chroms, read_lengths):
    # tandem duplications
    chroms = valid_chroms
    strand = ['forward', 'reverse']
    for n in range(n_seqs):
        names = [f">duplication_"]
        c = random.choice(chroms)
        s = random.choice(strand)
        blk = 0

        while blk < 1:
            flen = frag_lengths.get_fragment_length()
            read_len = read_lengths.get_fragment_length()
            if flen < 15 or read_len < 15 or flen > read_len:
                continue
            pos_read = random.randint(read_len, ref.get_reference_length(c) - read_len)
            pos_sv = random.randint(1, read_len-flen)
            blk += 1

        seq_start = ref.fetch(c, pos_read, pos_read + pos_sv - 1).upper()
        seq_dup = ref.fetch(c, pos_read + pos_sv, pos_read + pos_sv + flen - 1).upper()
        seq_end = ref.fetch(c, pos_read + pos_sv + flen, pos_read + read_len).upper()

        names.append(f"{c}:{pos_read}-{pos_read+pos_sv+flen-1}")
        names.append(f"{c}:{pos_read+pos_sv}-{pos_read+pos_sv+flen-1}")
        names.append(f"{c}:{pos_read+pos_sv+flen}-{pos_read+read_len}")

        if s == 'reverse':
            seq_dup_rev = misc.reverse_complement(seq_dup)
            seqs = [seq_start, seq_dup, seq_dup_rev, seq_end]
        else:
            seqs = [seq_start, seq_dup, seq_dup, seq_end]

        final_seq = "".join(seqs)
        final_name = "_".join(names)
        print(final_name)
        print(final_seq)


def generate_translocation(args, ref, n_seqs, frag_lengths, read_lengths):
    chroms = list(ref.references)
    strand = ['forward', 'reverse']
    order = ['short_first', 'long_first']
    for n in range(n_seqs):
        names = [f">translocation_"]
        blk = 0
        while blk < 1:
            c1 = random.choice(chroms)
            c2 = random.choice(chroms)
            s1 = random.choice(strand)
            s2 = random.choice(strand)
            o = random.choice(order)
            read_len = read_lengths.get_fragment_length()
            flen = frag_lengths.get_fragment_length()
            if flen < 15 or read_len < 15 or flen > read_len:
                continue
            if ref.get_reference_length(c2) < read_len:
                continue
            pos = random.randint(1, ref.get_reference_length(c1) - flen)
            pos2 = random.randint(1, ref.get_reference_length(c2) - read_len)
            blk += 1

            seq1 = ref.fetch(c1, pos, pos + flen).upper()
            seq2 = ref.fetch(c2, pos2, pos2 + read_len - flen).upper()

            if s1 == 'reverse':
                seq1 = misc.reverse_complement(seq1)
            if s2 == 'reverse':
                seq2 = misc.reverse_complement(seq2)

        if o == 'short_first':
            ins_seqs = [seq1, seq2]
            names.append(f"{c1}:{pos}-{pos + flen}")
            names.append(f"{c2}:{pos2}-{pos2 + read_len - flen}")
        else:
            ins_seqs = [seq2, seq1]
            names.append(f"{c2}:{pos2}-{pos2 + read_len - flen}")
            names.append(f"{c1}:{pos}-{pos + flen}")

        final_seq = "".join(ins_seqs)
        final_name = "_".join(names)
        print(final_name)
        print(final_seq)


def generate_inversion(args, ref, n_seqs, frag_lengths, valid_chroms, read_lengths):
    chroms = valid_chroms
    strand = ['forward', 'reverse']
    for n in range(n_seqs):
        c = random.choice(chroms)
        s = random.choice(strand)
        blk = 0
        while blk < 1:
            flen = frag_lengths.get_fragment_length()
            read_len = read_lengths.get_fragment_length()
            if flen < 15 or read_len < 15 or flen > read_len:
                continue
            pos_read = random.randint(1, ref.get_reference_length(c) - read_len)
            pos_inv = random.randint(1, read_len-flen)
            blk += 1

        seq1 = ref.fetch(c, pos_read, pos_read + pos_inv - 1).upper()
        seq2 = ref.fetch(c, pos_read + pos_inv, pos_read + pos_inv + flen - 1).upper()
        seq3 = ref.fetch(c, pos_read + pos_inv + flen, pos_read + read_len).upper()

        if s == 'reverse':
            seq1 = misc.reverse_complement(seq1)
            seq3 = misc.reverse_complement(seq3)
        else:
            seq2 = misc.reverse_complement(seq2)

        ins_seqs = [seq1, seq2, seq3]
        names = [f">inversion3_", f"{c}:{pos_read}-{pos_read+pos_inv-1}",
                 f"{c}:{pos_read+pos_inv}-{pos_read+pos_inv+flen-1}",
                 f"{c}:{pos_read+pos_inv+flen}-{pos_read+read_len}"]
        final_seq = "".join(ins_seqs)
        final_name = "_".join(names)
        print(final_name)
        print(final_seq)


def generate_deletion(args, ref, n_seqs, frag_lengths, valid_chroms, read_lengths):
    chroms = valid_chroms
    strand = ['forward', 'reverse']
    for n in range(n_seqs):
        c = random.choice(chroms)
        s = random.choice(strand)
        blk = 0
        while blk < 1:
            flen = frag_lengths.get_fragment_length()
            read_len = read_lengths.get_fragment_length()
            if flen < 15 or read_len < 15 or flen > read_len:
                continue
            pos = random.randint(1, ref.get_reference_length(c) - read_len)
            pos_del = random.randint(1, read_len)
            blk += 1

        seq1 = ref.fetch(c, pos, pos + pos_del).upper()
        seq3 = ref.fetch(c, pos + pos_del + flen,  pos + read_len + flen).upper()

        if s == 'reverse':
            seq1 = misc.reverse_complement(seq1)
            seq3 = misc.reverse_complement(seq3)

        seqs = [seq1, seq3]
        names = [f">deletion_", f"{c}:{pos}-{pos+pos_del}",
                 f"{c}:{pos+pos_del+flen}-{pos+read_len+flen}"]
        final_seq = "".join(seqs)
        final_name = "_".join(names)
        print(final_name)
        print(final_seq)


def generate_insertion(args, ref, n_seqs, frag_lengths, valid_chroms, read_lengths):
    chroms = valid_chroms
    strand = ['forward', 'reverse']
    for n in range(n_seqs):
        c_ins = random.choice(chroms)
        c = random.choice(chroms)
        s = random.choice(strand)
        blk = 0
        while blk < 1:
            flen = frag_lengths.get_fragment_length()
            read_len = read_lengths.get_fragment_length()
            if flen < 15 or read_len < 15 or flen > read_len:
                continue
            p1 = random.randint(1, ref.get_reference_length(c) - read_len)
            p2 = random.randint(1, ref.get_reference_length(c_ins) - flen)
            pos_ins = random.randint(1, read_len-flen)
            blk += 1

        seq1 = ref.fetch(c, p1, p1 + pos_ins - 1).upper()
        seq2 = ref.fetch(c, p2, p2 + flen).upper()
        seq3 = ref.fetch(c, p1 + pos_ins, p1 + read_len - flen).upper()

        if s == 'reverse':
            seq1 = misc.reverse_complement(seq1)
            seq2 = misc.reverse_complement(seq2)
            seq3 = misc.reverse_complement(seq3)

        seqs = [seq1, seq2, seq3]
        names = [f">insertion_", f"{c}:{p1}-{p1+pos_ins-1}",
                 f"{c}:{p2}-{p2+flen}",
                 f"{c}:{p1+pos_ins}-{p1+read_len-flen}"]
        final_seq = "".join(seqs)
        final_name = "_".join(names)
        print(final_name)
        print(final_seq)


def dna(length):
    dna_seq = ''
    for le in range(length):
        dna_seq += random.choice('CGTA')
    return dna_seq


def generate_random_insertion(args, ref, n_seqs, frag_lengths, valid_chroms, read_lengths):
    chroms = valid_chroms
    strand = ['forward', 'reverse']
    for n in range(n_seqs):
        c = random.choice(chroms)
        s = random.choice(strand)
        blk = 0
        while blk < 1:
            flen = frag_lengths.get_fragment_length()
            read_len = read_lengths.get_fragment_length()
            if flen < 15 or read_len < 15 or flen > read_len:
                continue
            pos = random.randint(1, ref.get_reference_length(c) - read_len)
            pos_ins = random.randint(1, read_len - flen)
            blk += 1

        seq1 = ref.fetch(c, pos, pos + pos_ins - 1).upper()
        seq2 = dna(flen)
        seq3 = ref.fetch(c, pos+pos_ins, pos + read_len - flen).upper()

        if s == 'reverse':
            seq1 = misc.reverse_complement(seq1)
            seq3 = misc.reverse_complement(seq3)

        seqs = [seq1, seq2, seq3]
        names = [f">randominsertion_", f"{c}:{pos}-{pos+pos_ins-1}",
                 f"randomchr:0-0",
                 f"{c}:{pos+pos_ins}-{pos+read_len-flen}"]
        final_seq = "".join(seqs)
        final_name = "_".join(names)
        print(final_name)
        print(final_seq)


def generate_n_insertion(args, ref, n_seqs, frag_lengths, valid_chroms, read_lengths):
    chroms = valid_chroms
    strand = ['forward', 'reverse']
    for n in range(n_seqs):
        c = random.choice(chroms)
        s = random.choice(strand)
        blk = 0
        while blk < 1:
            flen = frag_lengths.get_fragment_length()
            read_len = read_lengths.get_fragment_length()
            if flen < 15 and read_len < 15 or flen > read_len:
                continue
            pos = random.randint(1, ref.get_reference_length(c) - read_len)
            pos_ins = random.randint(1, read_len - flen)
            blk += 1

        seq1 = ref.fetch(c, pos, pos + pos_ins - 1).upper()
        seq2 = flen * 'N'
        seq3 = ref.fetch(c, pos + pos_ins, pos + read_len - flen).upper()

        if s == 'reverse':
            seq1 = misc.reverse_complement(seq1)
            seq3 = misc.reverse_complement(seq3)

        seqs = [seq1, seq2, seq3]
        names = [f">ninsertion_", f"{c}:{pos}-{pos+pos_ins-1}",
                 f"N:0-0",
                 f"{c}:{pos+pos_ins}-{pos+read_len-flen}"]
        final_seq = "".join(seqs)
        final_name = "_".join(names)
        print(final_name)
        print(final_seq)


def generate_flanking(args):
    ref = pysam.FastaFile(args.reference)

    print(f"Generating {args.number} SVs", file=sys.stderr)

    frag_lengths = fragment_lengths.FragmentLengths(args.mean_block_len, args.std_block_len)
    read_lengths = fragment_lengths.FragmentLengths(args.mean_read_len, args.std_read_len)

    sample = [read_lengths.get_fragment_length() for _ in range(10000)]
    max_len = max(sample)
    chroms = list(ref.references)
    valid_chroms = []
    for c in chroms:
        if 3 * max_len < ref.get_reference_length(c):
            valid_chroms.append(c)

    generate_duplication(args, ref, args.number, frag_lengths, valid_chroms, read_lengths)
    print(f"Generated {args.number} duplications")
    generate_deletion(args, ref, args.number, frag_lengths, valid_chroms, read_lengths)
    print(f"Generated {args.number} deletions")
    generate_random_insertion(args, ref, args.number, frag_lengths, valid_chroms, read_lengths)
    print(f"Generated {args.number} rand ins")
    generate_n_insertion(args, ref, args.number, frag_lengths, valid_chroms, read_lengths)
    print(f"Generated {args.number} n ins")
    generate_insertion(args, ref, args.number, frag_lengths, valid_chroms, read_lengths)
    print(f"Generated {args.number} insertions")
    generate_inversion(args, ref, args.number, frag_lengths, valid_chroms, read_lengths)
    print(f"Generated {args.number} inversions")
    generate_translocation(args, ref, args.number, frag_lengths, read_lengths)
    print(f"Generated {args.number} translocations")

    print(f"Done", file=sys.stderr)
