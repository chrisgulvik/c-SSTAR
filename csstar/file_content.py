import os
import sys
from collections import Counter

from Bio import SeqIO


def evaluate_nucleotide_composition(infile, min_ACGT):
    min_ACGT /= 100.
    seqlen, c = 0, Counter()
    for rec in SeqIO.parse(infile, 'fasta'):
        seqlen += len(rec.seq)
        c += Counter(list(rec.seq.upper()))
    ACGT_fraction = float(c['A'] + c['C'] + c['G'] + c['T']) / seqlen
    if ACGT_fraction < min_ACGT:
        sys.stderr.write('ERROR: low composition ({:.1f}%) of unambiguous '
                         'nucleotides in {}\n'.
                         format(ACGT_fraction * 100, os.path.basename(infile)))
        sys.exit(1)
    if ACGT_fraction != 1.:
        sys.stderr.write(
            'WARNING: {:.6f}% ambiguous nucleotides detected in {}\n'.
            format(100 - ACGT_fraction * 100, os.path.basename(infile)))
