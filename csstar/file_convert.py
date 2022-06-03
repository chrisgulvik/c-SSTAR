import gzip
import os
import shutil
import sys

from Bio import SeqIO
from Bio.SeqUtils import GC


def decompress_file(infile, outdir):
    uncompressed_file = os.path.basename(infile).rstrip('.gz')
    outfile = os.path.join(outdir, uncompressed_file)
    with gzip.open(infile, 'rb') as ifh, open(outfile, 'wb') as ofh:
        shutil.copyfileobj(ifh, ofh)
    return outfile


def genbank_to_fasta(infile, outdir):
    fasta_file = os.path.basename(infile).rsplit('.', 1)[0] + '.fa'
    outfile = os.path.join(outdir, fasta_file)
    records = []
    for rec in SeqIO.parse(infile, 'genbank'):
        if float(GC(rec.seq)) == 0:
            sys.stderr.write('ERROR: {} appears to lack nucleotides'.
                             format(rec.name))
            sys.exit(1)
        records.append(rec)
    with open(outfile, 'w') as ofh:
        SeqIO.write(records, ofh, 'fasta')
    return outfile
