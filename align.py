import subprocess
import glob
import argparse
import time

parse = argparse.ArgumentParser()
parse.add_argument('--threads', help='')
parse.add_argument('--ref', help='')
parse.add_argument('--input', help='')
parse.add_argument('--output_path', help='')
parse.add_argument('--bwa', help='', action='store_true')
parse.add_argument('--bwa_dodi', help='', action='store_true')
parse.add_argument('--minimap2', help='', action='store_true')
parse.add_argument('--minimap2_dodi', help='', action='store_true')
parse.add_argument('--last_index', help='', action='store_true')
parse.add_argument('--samtools_index', help='', action='store_true')
parse.add_argument('--lastal', help='', action='store_true')
parse.add_argument('--lastal_dodi', help='', action='store_true')
parse.add_argument('--lastal_lastsplit', help='', action='store_true')
parse.add_argument('--lastalsplit', help='', action='store_true')
parse.add_argument('--ngmlr', help='', action='store_true')
parse.add_argument('--vacmap', help='', action='store_true')
args = parse.parse_args()

###################### BWA
if args.bwa:
    bwa = "bwa mem -t {threads} {ref} {input} | "\
          "samtools view -bh - | "\
          "samtools sort -o {output_path}/bwa.bam; "\
          "samtools index {output_path}/bwa.bam".format(threads=args.threads,
                                                            ref=args.ref,
                                                            input=args.input,
                                                            output_path=args.output_path)
    start_time = time.time()
    subprocess.run(bwa, shell=True)
    elapsed_time = time.time() - start_time
    print(f"bwa took {elapsed_time:.2f} seconds")

###################### BWA | DODI
if args.bwa_dodi:
    # -a, dodi
    bwa_dodi = "bwa mem -a -t {threads} {ref} {input} | "\
               "dodi --paired False - | "\
               "samtools view -bh - | "\
               "samtools sort -o {output_path}/bwa_dodi.bam; "\
               "samtools index {output_path}/bwa_dodi.bam".format(threads=args.threads,
                                                                      ref=args.ref,
                                                                      input=args.input,
                                                                      output_path=args.output_path)
    start_time = time.time()
    subprocess.run(bwa_dodi, shell=True)
    elapsed_time = time.time() - start_time
    print(f"bwa | dodi took {elapsed_time:.2f} seconds")

###################### MINIMAP2
if args.minimap2:
    minimap2 = "minimap2 -ax map-ont -t {threads} {ref} {input} | "\
               "samtools view -bh - | "\
               "samtools sort -o {output_path}/minimap2.bam; "\
               "samtools index {output_path}/minimap2.bam".format(threads=args.threads,
                                                                           ref=args.ref,
                                                                           input=args.input,
                                                                           output_path=args.output_path)
    start_time = time.time()
    subprocess.run(minimap2, shell=True)
    elapsed_time = time.time() - start_time
    print(f"minimap2 took {elapsed_time:.2f} seconds")

###################### MINIMAP2 | DODI
if args.minimap2_dodi:
    #-N 100 --secondary=yes -p 0, dodi
    minimap2_dodi = "minimap2 -N 100 --secondary=yes -p 0 -ax map-ont -t {threads} {ref} {input} | "\
                    "dodi --paired False - | "\
                    "samtools view -bh - | "\
                    "samtools sort -o {output_path}/minimap2_dodi.bam; "\
                    "samtools index {output_path}/minimap2_dodi.bam".format(threads=args.threads,
                                                                            ref=args.ref,
                                                                            input=args.input,
                                                                            output_path=args.output_path)
    start_time = time.time()
    subprocess.run(minimap2_dodi, shell=True)
    elapsed_time = time.time() - start_time
    print(f"minimap2 | dodi took {elapsed_time:.2f} seconds")

###################### LAST
if args.samtools_index:
    awk_command = 'awk \'{print "@SQ\\tSN:"$1"\\tLN:"$2}\' {ref}.fai > '\
                  '{output_path}/sq_lines.txt'.format(ref=args.ref,output_path=args.output_path)
    subprocess.run(awk_command, shell=True)

if args.last_index:
    create_index = "lastdb -P {threads} ref ref".format(threads=args.threads, ref=args.ref)
    subprocess.run(create_index, shell=True)

if args.lastal or args.lastal_dodi or args.lastalsplit or args.lastal_lastsplit:
    if not len(glob.glob(f"args.output_path/trained_parameters.txt")) == 1:
        last_train = "last-train -Q0 {ref} {input} > {output_path}/trained_parameters.txt".format(ref=args.ref,
                                                                                                  input=args.input,
                                                                                                  output_path=args.output_path)
        start_time = time.time()
        subprocess.run(last_train, shell=True)
        elapsed_time = time.time() - start_time
        print(f"last-train took {elapsed_time:.2f} seconds")

###################### LASTAL | LAST-SPLIT
if args.lastal_lastsplit:
    lastal_lastsplit = "lastal -p {output_path}/trained_parameters.txt -P {threads} -Q0 {ref} {input} "\
                       "| last-split | maf-convert sam > {output_path}/lastal_lastsplit.sam; "\
                       "samtools view -H {output_path}/lastal_lastsplit.sam > {output_path}/old_header.sam; "\
                       "cat {output_path}/old_header.sam {output_path}/sq_lines.txt > {output_path}/new_header.sam; "\
                       "cat {output_path}/new_header.sam {output_path}/lastal_lastsplit.sam | "\
                       "samtools view -bh - | "\
                       "samtools sort -o {output_path}/lastal_lastsplit.bam; "\
                       "samtools index {output_path}/lastal_lastsplit.bam".format(ref=args.ref,
                                                                                  input=args.input,
                                                                                  output_path=args.output_path,
                                                                                  threads=args.threads)
    start_time = time.time()
    subprocess.run(lastal_lastsplit, shell=True)
    elapsed_time = time.time() - start_time
    print(f"lastal | last-split took {elapsed_time:.2f} seconds")

###################### LASTAL --SPLIT
if args.lastalsplit:
    # lastal --split
    lastalsplit = "lastal -p {output_path}/trained_parameters.txt -P{threads} -Q0 {ref} {input} "\
                  "| last-split | maf-convert sam > {output_path}/lastalsplit.sam; "\
                  "samtools view -H {output_path}/lastalsplit.sam > {output_path}/old_header.sam; "\
                  "cat {output_path}/old_header.sam {output_path}/sq_lines.txt > {output_path}/new_header.sam; "\
                  "cat {output_path}/new_header.sam {output_path}/lastalsplit.sam |" \
                  "samtools view -bh - | "\
                  "samtools sort -o {output_path}/lastalsplit.bam; "\
                  "samtools index {output_path}/lastalsplit.bam".format(ref=args.ref,
                                                                        input=args.input,
                                                                        output_path=args.output_path,
                                                                        threads=args.threads)
    start_time = time.time()
    subprocess.run(lastalsplit, shell=True)
    elapsed_time = time.time() - start_time
    print(f"lastal --split took {elapsed_time:.2f} seconds")

###################### LASTAL
if args.lastal:
    lastal = "lastal -p {output_path}/trained_parameters.txt -P{threads} -Q0 {ref} {input} | " \
             "maf-convert sam > {output_path}/lastal.sam; " \
             "samtools view -H {output_path}/lastal.sam > {output_path}/old_header.sam; " \
             "cat {output_path}/old_header.sam {output_path}/sq_lines.txt > {output_path}/new_header.sam; " \
             "cat {output_path}/new_header.sam {output_path}/lastal.sam |" \
             "samtools view -bh - | " \
             "samtools sort -o {output_path}/lastal.bam; " \
             "samtools index {output_path}/lastal.bam".format(ref=args.ref,
                                                              input=args.input,
                                                              output_path=args.output_path,
                                                              threads=args.threads)
    start_time = time.time()
    subprocess.run(lastal, shell=True)
    elapsed_time = time.time() - start_time
    print(f"lastal took {elapsed_time:.2f} seconds")

###################### LASTAL | DODI
if args.lastal:
    lastal = "lastal -p {output_path}/trained_parameters.txt -P {threads} -Q0 {ref} {input} | "\
              "maf-convert sam > {output_path}/lastal.sam; "\
              "samtools view -H {output_path}/lastal.sam > {output_path}/old_header.sam; "\
              "cat {output_path}/old_header.sam {output_path}/sq_lines.txt > {output_path}/new_header.sam; "\
              "cat {output_path}/new_header.sam {output_path}/lastal.sam | "\
              "samtools view -bh - | "\
              "samtools sort -o {output_path}/lastal.bam; "\
              "samtools index {output_path}/lastal.bam".format(ref=args.ref,
                                                               input=args.input,
                                                               output_path=args.output_path,
                                                               threads=args.threads)

if args.lastal_dodi and args.lastal:
    lastal_dodi = "dodi --paired False {output_path}/lastal.sam | "\
                  "samtools view -bh - | "\
                  "samtools sort -o {output_path}/lastal_dodi.bam; "\
                  "samtools index {output_path}/lastal_dodi.bam".format(ref=args.ref,
                                                                        input=args.input,
                                                                        output_path=args.output_path,
                                                                        threads=args.threads)
    start_time = time.time()
    subprocess.run(lastal_dodi, shell=True)
    elapsed_time = time.time() - start_time
    print(f"dodi took {elapsed_time:.2f} seconds")

if args.lastal_dodi and not args.lastal:
    lastal_dodi = "lastal -p {output_path}/trained_parameters.txt -P {threads} -Q0 {ref} {input} | "\
                  "dodi --paired False - | maf-convert sam > {output_path}/lastal_dodi.sam; "\
                  "samtools view -H {output_path}/lastal_ldodi.sam > {output_path}/old_header.sam; "\
                  "cat {output_path}/old_header.sam {output_path}/sq_lines.txt > {output_path}/new_header.sam; "\
                  "cat {output_path}/new_header.sam {output_path}/lastal_dodi.sam | "\
                  "samtools view -bh - | "\
                  "samtools sort -o {output_path}/lastal_dodi.bam; "\
                  "samtools index {output_path}/lastal_dodi.bam".format(ref=args.ref,
                                                                        input=args.input,
                                                                        output_path=args.output_path,
                                                                        threads=args.threads)
    start_time = time.time()
    subprocess.run(lastal_dodi, shell=True)
    elapsed_time = time.time() - start_time
    print(f"lastal | dodi took {elapsed_time:.2f} seconds")

###################### NGMLR
if args.ngmlr:
    ngmlr = "ngmlr -t {threads} -r {ref} -q {input} -x ont --verbose | "\
            "samtools view -bh - | "\
            "samtools sort -o {output_path}/ngmlr.bam; "\
            "samtools index {output_path}/ngmlr.bam".format(ref=args.ref,
                                                            input=args.input,
                                                            output_path=args.output_path,
                                                            threads=args.threads)
    start_time = time.time()
    subprocess.run(ngmlr, shell=True)
    elapsed_time = time.time() - start_time
    print(f"ngmlr took {elapsed_time:.2f} seconds")

###################### VACMAP
if args.vacmap:
    vacmap = "vacmap -ref {ref} -read {input} -mode S -t {threads} | samtools sort -@4 > vacmap.bam; "\
             "samtools index -@4 vacmap.bam".format(ref=args.ref,
                                                    input=args.input,
                                                    output_path=args.output_path,
                                                    threads=args.threads)
    start_time = time.time()
    subprocess.run(vacmap, shell=True)
    elapsed_time = time.time() - start_time
    print(f"vacmap took {elapsed_time:.2f} seconds")