#!/usr/bin/env python

import argparse, logging, os, pwd, re, subprocess, sys
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

def parseArgs(args=None):
	parser = argparse.ArgumentParser(description='c-SSTAR is a CLI utility for rapidly identifying antibiotic resistance gene determinants in bacterial genomes')
	parser.add_argument('-d', '--database', required=True, help='a SSTAR-formatted FastA database of AR gene sequences')
	parser.add_argument('-g', '--genome', required=True, help='a FastA genome')
	parser.add_argument('-o', '--outdir', default=os.getcwd(), help='output directory')
	parser.add_argument('-s', '--similarity', type=int, default=95, help='minimum percent nucleotide similarity')
	return parser.parse_args()

def syscall(syscmd):
	with open(os.devnull) as dump: 
		returncode = subprocess.call(syscmd, stdout=dump, stderr=dump, shell=True)
		if returncode != 0:
			logging.error('failed syscall ' + syscmd)
			sys.exit('ERROR: failed syscall ' + syscmd)

def translateSeq(nuclSeq):
	proteinSeq = Seq(nuclSeq, generic_dna).translate(cds=False, to_stop=False, stop_symbol='*')
	return proteinSeq

def internalSTOPcodon(candidate):
	'''counts number of internal stop codons; requires sequence from database to be in frame'''
	nucSeq = candidate[10]
	if int(int(candidate[8])%3) == 1:
		frameStart = int(0)
	elif int(int(candidate[8])%3) == 0:
		frameStart = int(1)
	elif int(int(candidate[8])%3) == 2:
		frameStart = int(2)
	if frameStart > 0:
		nucSeq = nucSeq[frameStart:]
	frameEnd = int(int(candidate[9])%3)
	if frameEnd > 0:
		nucSeq = nucSeq[:-frameEnd]
	if len(nucSeq)%3 == 0:
		protein = translateSeq(nucSeq)
	else:
		sys.exit('ERROR: incorrect nucleotide length (%s) after trimming' % len(nucSeq))
	numInternalSTOPcodons = len(re.findall(r'\*[ABCDEFGHIKLMNPQRSTVWYZ]', str(protein)))
	return (protein, str(numInternalSTOPcodons))

def tagHit(l):
	if int(l[5]) != int(l[7]):  #incomplete len; SRST2='? indicates that there was uncertainty in at least one of the alleles'
		l = [ l[0], l[1]+'?', l[2]+'?', l[3], l[4], l[5], l[6], l[7], l[8], l[9], l[10] ]
	if int(l[5]) == int(l[7]) and int(l[4]) < 100:  #full length match but pident!=100%; SRST2='* [...] indicates that there were mismatches against at least one of the alleles. This suggests that you have a novel variant [...] rather than a precise match'
		l = [ l[0], l[1]+'*', l[2]+'*', l[3], l[4], l[5], l[6], l[7], l[8], l[9], l[10] ]
	(prot, numSTOP) = internalSTOPcodon(l)
	if int(numSTOP) > 0:  #TR indicates 'truncated protein translation'
		l = [ l[0], l[1]+'TR', l[2]+'TR', l[3], l[4], l[5], l[6], l[7], l[8], l[9], l[10] ]
	return l

def main(args=None):
	args = parseArgs()
	outdir = args.outdir
	genome = args.genome
	baseGenome = os.path.splitext(os.path.basename(genome))[0]
	if not os.path.exists(outdir): 
		os.mkdir(outdir)
	logging.basicConfig(filename='%s/c-SSTAR_%s.log' % (outdir, baseGenome), format='%(asctime)s: %(levelname)s: %(message)s', datefmt='%d-%m-%Y %I:%M:%S %p', level=logging.INFO)
	logging.info('user: %s' % pwd.getpwuid(os.getuid()).pw_name)
	logging.info('release: %s' % os.uname()[3])
	logging.info('shell env: %s' % pwd.getpwuid(os.getuid()).pw_shell)
	logging.info('cwd: %s' % pwd.getpwuid(os.getuid()).pw_dir)
	logging.info('python version: %s' % sys.version)
	logging.info(subprocess.check_output('command -v blastn', shell=True).rstrip())
	logging.info(subprocess.check_output('blastn -version | tail -n1', shell=True).rstrip())
	database = args.database
	similarity = args.similarity
	syscall('makeblastdb -in ' + genome + ' -out ' + os.path.join(outdir, baseGenome) + ' -dbtype nucl')
	syscall('blastn -task blastn -query {0} -db {out} -out {out}.blastn.tsv -evalue 1e-5 -max_target_seqs 1 -ungapped -perc_identity {1} -culling_limit 1 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen sseq"'.format(database, similarity, out=os.path.join(outdir, baseGenome)))
	syscall('rm ' + os.path.join(outdir, baseGenome) + '.nin')
	syscall('rm ' + os.path.join(outdir, baseGenome) + '.nsq')
	syscall('rm ' + os.path.join(outdir, baseGenome) + '.nhr')
	with open(os.path.join(outdir, baseGenome) + '.blastn.tsv') as infile:
		best = ['-1','a','a','a',0,'-1','a',0,'a','a','a']
		currentClusterNr = '-1'
		topHits = []
		for l in infile:
			blastOut = [x for x in l.split('\t')]
			thresHold = (int(blastOut[12])/5)*2
			if int(blastOut[3]) > thresHold:  #>40% overlap required
				enzymeParts = [s for s in blastOut[0].split('__')]
				clusterNr = enzymeParts[0]
				ident, dec = blastOut[2].split('.')
				pident = int(ident)
				bitscore = blastOut[11]
				candidate = [clusterNr, enzymeParts[1], enzymeParts[2], blastOut[1], pident, blastOut[3], bitscore, blastOut[12], blastOut[6], blastOut[7], blastOut[13].rstrip()]
				if clusterNr == currentClusterNr:
					if int(best[6]) < int(candidate[6]):
						best = candidate
				else:
					if best[4] >= similarity:
						topHits.append(best)
					currentClusterNr = clusterNr
					best = candidate
		if best[4] >= similarity:
			topHits.append(best)
	for i in topHits:
		hit = tagHit(i)
		print hit[1] + '\t' + hit[2] + '\t' + hit[3] + '\t' + str(hit[4]) + '%\t' + str(hit[5]) + '\t' + str(hit[7])

if __name__ == '__main__':
	main()