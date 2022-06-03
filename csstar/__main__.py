#!/usr/bin/env python3


__title__ = 'csstar'
__description__ = (
    'Antibiotic resistance gene detection from contigs in FastA or GenBank'
    ' format')
__author__ = 'Christopher A. Gulvik'
__license__ = 'Apache 2.0'


import logging
import ntpath
import os
import pwd
import shutil
import subprocess
import sys
from tempfile import mkdtemp

from .evaluate_alignment import tag_hit
from .file_content import evaluate_nucleotide_composition
from .file_convert import decompress_file, genbank_to_fasta
from .parse_arguments import parse_args
from .utilities_file import require_file_exists_and_test_empty
from .utilities_system import (
    check_python_version, require_dependency, sys_call)
from .version import __version__


def main():
    check_python_version()
    args = parse_args(sys.argv[1:])
    require_dependency('makeblastdb')
    require_dependency('blastn')
    tmp = mkdtemp()
    # NOTE: convert infiles to dict but simple universal way in py2.7 and py3+
    # to iterate items to pass to functions and reassign? iteritems() items()
    infile_query = os.path.realpath(os.path.expanduser(args.genome))
    if infile_query.endswith('.gz'):
        infile_query = decompress_file(infile_query, tmp)
    if infile_query.endswith(('.gbff', '.gbf', '.gbk', '.gb')):
        infile_query = genbank_to_fasta(infile_query, tmp)
    infile_db = os.path.realpath(os.path.expanduser(args.database))
    if infile_db.endswith('.gz'):
        infile_db = decompress_file(infile_db, tmp)
    for file in [infile_db, infile_query]:
        evaluate_nucleotide_composition(file, args.min_ACGT)
    outdir = os.path.realpath(os.path.expanduser(args.outdir))
    min_identity = args.min_identity
    if args.basename is not None:
        base_genome = args.basename
    else:
        base_genome = os.path.splitext(os.path.basename(infile_query))[0]
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    logging.basicConfig(
        filename=os.path.join(outdir, 'c-SSTAR_{}.log'.format(base_genome)),
        format='%(asctime)s: %(levelname)s: %(message)s',
        datefmt='%d-%m-%Y %I:%M:%S %p', level=logging.INFO)
    logging.info('c-SSTAR version: {}'.format(__version__))
    logging.info('user: {}'.format(pwd.getpwuid(os.getuid()).pw_name))
    logging.info('release: {}'.format(os.uname()[3]))
    logging.info('shell env: {}'.format(pwd.getpwuid(os.getuid()).pw_shell))
    logging.info('cwd: {}'.format(pwd.getpwuid(os.getuid()).pw_dir))
    logging.info('python version: {}'.format(sys.version))
    logging.info(subprocess.check_output(
        'command -v blastn', shell=True).rstrip())
    logging.info(subprocess.check_output(
        'blastn -version | tail -n 1', shell=True).rstrip())
    tmp_query_idx = os.path.join(outdir, base_genome)
    sys_call('makeblastdb -in {} -out {} -dbtype nucl'.format(
        infile_query, tmp_query_idx))
    for x in ['.nin', '.nsq', '.nhr']:
        empty = require_file_exists_and_test_empty(tmp_query_idx + x)
        if empty:
            sys.stderr.write(
                'ERROR: failed making BLAST database to'
                ' search\n{} file empty\n'.
                format(os.path.join(outdir, base_genome + x)))
            sys.exit(1)
    cmd = ('blastn -task blastn -query {} -db {} -out {} -evalue 1e-5'
           ' -max_target_seqs 1 -perc_identity {} -culling_limit 1'
           ' -outfmt "6 qseqid sseqid pident length mismatch gaps qstart qend'
           ' sstart send evalue bitscore qlen slen sseq" -num_threads {}'.
           format(infile_db, tmp_query_idx,
                  os.path.join(outdir, base_genome + '.blastn.tsv'),
                  min_identity, args.cpus))
    if args.ungapped:
        cmd += ' -ungapped'
    sys_call(cmd)
    shutil.rmtree(tmp)
    empty = require_file_exists_and_test_empty(
        os.path.join(outdir, base_genome + '.blastn.tsv'))
    if empty:
        sys.stderr.write('INFO: empty alignment output. No antimicrobial'
                         ' resistance genes detected\n')
        sys.exit(0)
    with open(os.path.join(outdir, base_genome + '.blastn.tsv')) as infile:
        # best = [-1, 'a', 'a', 'a', 0., 0, 0., 0, 0, 0, 0, 'a']
        best = {
            'clust_nr': -1,
            'AR_family': 'a',
            'AR_variant': 'a',
            'qry_defln': 'a',
            'perc_id': 0.,
            'aln_len': 0,
            'bits': 0.,
            'qlen': 0,
            'qstart': 0,
            'qend': 0,
            'slen': 0,
            'sseq': 'a'}
        current_cluster_number = -1
        top_hits = []
        for l in infile:
            blast_out = [x for x in l.rstrip().split('\t')]
            threshold = float(blast_out[12]) * (args.min_aln_frac / 100.)
            if float(blast_out[3]) > threshold:
                db_defline_parts = [s for s in blast_out[0].split('__')]
                cluster_number = int(db_defline_parts[0])
                pident = float(blast_out[2])
                bitscore = float(blast_out[11])
                candidate = {
                    'clust_nr': cluster_number,
                    'AR_family': str(db_defline_parts[1]),
                    'AR_variant': str(db_defline_parts[2]),
                    'qry_defln': str(blast_out[1]),
                    'perc_id': pident,
                    'aln_len': int(blast_out[3]),
                    'bits': bitscore,
                    'qlen': int(blast_out[12]),
                    'qstart': int(blast_out[6]),
                    'qend': int(blast_out[7]),
                    'slen': int(blast_out[13]),
                    'sseq': str(blast_out[14])}
                if cluster_number == current_cluster_number:
                    if best['bits'] < candidate['bits']:
                        best = candidate
                else:
                    if best['perc_id'] >= min_identity:
                        top_hits.append(best)
                    current_cluster_number = cluster_number
                    best = candidate
        if best['perc_id'] >= min_identity:
            top_hits.append(best)
    if current_cluster_number == -1:
        sys.stderr.write('INFO: all alignments in output do not meet'
                         ' thresholds. No antimicrobial resistance genes'
                         ' detected\n')
        sys.exit(0)
    output = []
    for i in top_hits:
        hit = tag_hit(i, args.edge)
        output.append(hit['AR_family'] + '\t' + hit['AR_variant'] + '\t' +
                      hit['qry_defln'] + '\t' + str(hit['perc_id']) + '%\t' +
                      str(hit['aln_len']) + '\t' + str(hit['qlen']))
    if args.report is not None:
        report = os.path.realpath(os.path.expanduser(args.report))
        if not os.path.exists(os.path.dirname(outdir)):
            os.mkdir(os.path.dirname(outdir))
        with open(report, 'w') as ofh:
            ofh.write('\n'.join(output))
    else:
        for ln in output:
            print(ln)


if __name__ == '__main__':
    main()
