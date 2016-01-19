# c-SSTAR
##### CLI-Sequence Search Tool for Antimicrobial Resistance
*for those who lack a Java environment or are simply scared to dive into the deep abyss (GUI)*

![alt tag](https://github.com/chrisgulvik/images/raw/master/c-SSTAR.jpeg)


## System Requirements
- Linux or Mac OS X platform
- BLAST+ (blastn and makeblastdb)
- Python


## Input
1. FastA formatted genome
2. FastA database from SSTAR based on [ARG-ANNOT](https://github.com/tomdeman-bio/Sequence-Search-Tool-for-Antimicrobial-Resistance-SSTAR-/raw/master/ARG-ANNOT.srst2.fasta) or [ResFinder](https://github.com/tomdeman-bio/Sequence-Search-Tool-for-Antimicrobial-Resistance-SSTAR-/raw/master/ResFinder_12-14-2015.srst2.fasta), which are formatted according to Kat Holt's clustering approach for [SRST2](https://github.com/katholt/srst2/tree/master/database_clustering)


## Output
Tab-delimited data are printed to standard out with the following fields:
1. antimicrobial resistance gene name (from database)
2. sequence defline/header where the antimicrobial resistance gene is located (from genome)
3. % nucleotide similarity
4. bp length of alignment
5. bp length of antimicrobial resistance gene (from database)


## Usage
`python c-SSTAR.py -g <genome_file> -d <database_file>`


### Citing SSTAR
Please cite Tom's paper in mSphere: de Man TJB, Limbago BM. 2016. SSTAR, a stand-alone easy-to-use antimicrobial resistance gene predictor. mSphere 1(1): e00050-15. doi: 10.1128/mSphere.00050-15

When using the ARG-ANNOT database please also cite: Gupta SK, Padmanabhan BR, Diene SM, Lopez-Rojas R, Kempf M, Landraud L, Rolain J-M. 2014. ARG-ANNOT (Antibiotic Resistance Gene-ANNOTation), a new bioinformatic tool to discover antibiotic resistance genes in bacterial genomes. Antimicrobial Agents and Chemotherapy 58:212–220. doi: 10.1128/AAC.01310-13

When using the ResFinder database please also cite: Zankari E, Hasman H, Cosentino S, Vestergaard M, Rasmussen S, Lund O, Aarestrup F, Larsen MV. 2012. Identification of acquired antimicrobial resistance genes. Journal of Antimicrobial Chemotherapy 67:2640–2644. doi: 10.1093/jac/dks261
