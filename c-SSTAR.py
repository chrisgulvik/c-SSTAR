#!/usr/bin/env python

import argparse, logging, os, pwd, re, subprocess, sys

def parseArgs(args=None):
	parser = argparse.ArgumentParser(description='c-SSTAR is a CLI for rapidly identifying antimicrobial resistance gene determinants in bacterial genomes')
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
	syscall('makeblastdb -in ' + genome + ' -out ' + os.path.join(outdir, baseGenome) + ' -dbtype nucl')
	syscall('blastn -query {0} -db {out} -out {out}.blastn.tsv -evalue 1e-5 -max_target_seqs 1 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen sseq"'.format(database, out=os.path.join(outdir, baseGenome)))
	syscall('rm -f ' + os.path.join(outdir, baseGenome) + '.{nin,nih,nsq,nhr}')
	similarity = args.similarity
	with open(os.path.join(outdir, baseGenome) + '.blastn.tsv') as infile:
		best = ['-1','a','a','a',1,'-1','a','a']
		currentClusterNr = '-1'
		topHits = []
		for l in infile:
			parts = [x for x in l.split('\t')]
			alignLen = int(parts[3])
			qlen = parts[12]
			thresHold = (int(parts[12])/5)*2
			if alignLen > thresHold:  #>40% overlap required
				enzymeParts = [s for s in parts[0].split('__')]
				clusterNr = enzymeParts[0]
				ARgeneFam = enzymeParts[1]
				ARgeneVar = enzymeParts[2]
				blastTargetLabel = parts[1]
				ident, dec = parts[2].split('.')
				pident = int(ident)
				bitscore = parts[11]
				candidate = [clusterNr, ARgeneFam, ARgeneVar, blastTargetLabel, pident, alignLen, bitscore, qlen, parts[13]]
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
		if i[4] == 100:  #100% nucl similarity
			print i[2] + '\t' + i[3] + '\t' + str(i[4]) + '%\t' + str(i[5]) + '\t' + str(i[7])
		else:
			print 'A potentially new ' + i[1] + ' discovered\t' + i[3] + '\t' + str(i[4]) + '%\t' + str(i[5]) + '\t' + str(i[7])

if __name__ == '__main__':
	main()