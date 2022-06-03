import argparse
import os
import pathlib
import sys
from multiprocessing import cpu_count

from . import __version__


try:
    import importlib.resources as pkg_resources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as pkg_resources
data_path = pathlib.Path(pkg_resources.resource_filename(__name__, 'data'))


def parse_args():
    parser = argparse.ArgumentParser(
        description='c-SSTAR is a CLI utility'
        ' for rapidly identifying antibiotic resistance gene determinants in'
        ' bacterial genomes', add_help=False,
        epilog='NOTE: report output places special characters on alignments'
        ' to denote special interest. * is >=1 mismatch, ? is incomplete'
        ' align length, TR is truncation predicted with internal stop codon,'
        ' and $ is aligned at contig edge')
    req = parser.add_argument_group('Input')
    req.add_argument(
        '-d', '--database', metavar='FILE',
        default=os.path.join(data_path, 'ResGANNOT_srst2.fasta.gz'),
        help='a SSTAR-formatted FastA database of AR gene sequences,'
        ' optionally gunzip compressed [%(default)s]')
    req.add_argument(
        '-g', '--genome', required=True, metavar='FILE',
        help='query genome in FastA or GenBank format, optionally gunzip'
        ' compressed')
    aln = parser.add_argument_group('Alignment Options')
    aln.add_argument(
        '--cpus', type=require_int_nonnegative,
        metavar='INT', default=cpu_count(),
        help='threads for BLASTn to use [%(default)s]')
    aln.add_argument(
        '--min-aln-frac', type=require_float_1_to_100,
        default=40, metavar='FLOAT',
        help='minimum database alignment length percentage [%(default)s]')
    aln.add_argument(
        '--min-identity', type=require_float_1_to_100,
        metavar='FLOAT', default=95,
        help='minimum alignment nucleotide identity [%(default)s]')
    aln.add_argument(
        '--ungapped', action='store_true', default=False,
        help='store ungapped-only alignments for BLASTn [off]')
    opt = parser.add_argument_group('I/O Options')
    opt.add_argument(
        '--basename', type=str, metavar='STR',
        help='output file base identifier [basename <genome>]')
    opt.add_argument(
        '--edge', type=require_int_nonnegative,
        default=50, metavar='INT',
        help='label alignments at edge of contig (bp) [%(default)s]')
    opt.add_argument(
        '--min-ACGT', type=require_float_1_to_100, metavar='FLOAT',
        default=97,
        help='minimum ACGT percent in each input file to align [%(default)s]')
    opt.add_argument(
        '--outdir', type=str, metavar='PATH',
        default=os.getcwd(),
        help='output path [cwd]')
    opt.add_argument(
        '--report', default=None, metavar='FILE',
        help='summary tab-delimited report [stdout]')
    msc = parser.add_argument_group('Misc Options')
    msc.add_argument(
        '--debug', default=False, action='store_true',
        help='turn on debugging')
    msc.add_argument(
        '-h', '--help', action='help',
        help='show this help message and exit')
    msc.add_argument(
        '-v', '--version', action='version',
        version='%(prog)s v{}'.format(__version__))
    return parser.parse_args()


def require_int_nonnegative(x):
    try:
        if int(x) < 0 or '.' in str(x):
            sys.stderr.write(
                'ERROR: {} must be a non-negative integer\n'.
                format(x))
            sys.exit(1)
    except ValueError:
        sys.stderr.write('ERROR: {} must be an integer\n'.format(x))
        sys.exit(1)
    return int(x)


def require_float_1_to_100(x):
    try:
        x = float(x)
        if 1 < x > 100:
            sys.stderr.write('ERROR: {} must be in 1, 100 range\n'.format(x))
            sys.exit(1)
    except ValueError:
        sys.stderr.write('ERROR: {} must be a float\n'.format(x))
        sys.exit(1)
    return x
