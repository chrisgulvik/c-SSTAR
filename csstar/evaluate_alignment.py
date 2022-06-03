import re
import sys

from Bio.Alphabet import generic_dna
from Bio.Seq import Seq


def translate_seq(nucl_seq):
    '''use biopython to take in a nucleotide sequence (string) and
    return a protein sequence (string)'''
    protein_seq = str(Seq(nucl_seq, generic_dna).translate(
        cds=False, to_stop=False, stop_symbol='*'))
    return protein_seq


def count_internal_stop_codons(nucl_seq, qstart):
    '''counts number of internal stop codons; requires sequence from database
    to be in frame'''
    if '-' in nucl_seq:
        nucl_seq = nucl_seq.replace('-', '')
    else:
        nucl_seq = nucl_seq
    if int(int(qstart) % 3) == 1:
        frame_start = int(0)
    elif int(int(qstart) % 3) == 0:
        frame_start = int(1)
    elif int(int(qstart) % 3) == 2:
        frame_start = int(2)
    if frame_start > 0:
        nucl_seq = nucl_seq[frame_start:]
    if len(nucl_seq) % 3 == 0:
        protein = translate_seq(nucl_seq)
    elif len(nucl_seq) % 3 > 0:
        protein = translate_seq(nucl_seq[:-(len(nucl_seq) % 3)])
    else:
        sys.stderr.write('ERROR: incorrect nucleotide length ({}) after'
                         ' trimming\n'.format(len(nucl_seq)))
        sys.exit(1)
    num_internal_stop_codons = len(
        re.findall(r'\*[ABCDEFGHIKLMNPQRSTVWYZ]', protein))
    return num_internal_stop_codons


def tag_hit(d, edge):
    '''evaluates the alignment d and returns the same structure containing
    shorthand suffix labels denoting if an allele match was not full length or
    <100%, if the alignment was at the edge of a contig, and if any internal
    stop codons were predicted in the alignment'''
    suffix = ''
    if d['aln_len'] != d['qlen']:
        # incomplete length; SRST2='? indicates that there was
        #   uncertainty in at least one of the alleles'
        suffix += '?'
    if d['aln_len'] == d['qlen'] and d['perc_id'] < 100:
        # full length match but pident!=100%
        # SRST2='* [...] indicates that there were mismatches against at least
        #   one of the alleles. This suggests that you have a novel variant
        #   [...] rather than a precise match'
        suffix += '*'
    if (d['qstart'] - edge) < 0:
        # edge hits on left edge
        if (d['qlen'] - d['qend']) > 0:
            # incomplete alignlen at edge for '$' designation
            #   distinguishing ^ from $ is unnecessary
            suffix += '$'
    elif (d['qend'] + edge) > d['slen']:
        # edge hits on right edge
        if (d['qlen'] + d['qstart']) > d['slen']:
            # require incomplete alignlen
            suffix += '$'
    num_stop_codons = count_internal_stop_codons(d['sseq'], d['qstart'])
    if num_stop_codons > 0:
        # TR indicates 'truncated protein translation'
        suffix += 'TR'
    if suffix:
        d['AR_family'] = d['AR_family'] + suffix
        d['AR_variant'] = d['AR_variant'] + suffix
    return d
