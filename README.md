# c-SSTAR
##### CLI-Sequence Search Tool for Antimicrobial Resistance
*for those who lack a Java environment or are simply scared to dive into the deep abyss (GUI)*

![alt tag](https://github.com/chrisgulvik/images/raw/master/c-SSTAR.jpeg)


## System Requirements
- Linux or Mac OS X platform
- BLAST+ (blastn and makeblastdb)
- Python and Biopython

## Usage
    python c-SSTAR.py -g <genome_file> -d <database_file>

## Input
1. FastA formatted genome
2. FastA database of antimicrobial resistance (AR) gene sequences from SSTAR. Two databases are available ([ARG-ANNOT](https://github.com/tomdeman-bio/Sequence-Search-Tool-for-Antimicrobial-Resistance-SSTAR-/raw/master/ARG-ANNOT.srst2.fasta) or [ResFinder](https://github.com/tomdeman-bio/Sequence-Search-Tool-for-Antimicrobial-Resistance-SSTAR-/raw/master/ResFinder_12-14-2015.srst2.fasta)), which are formatted according to Kat Holt's clustering approach for [SRST2](https://github.com/katholt/srst2/tree/master/database_clustering)

## Output
#### I) Summary output (to stdout)
A tab-delimited summary is printed to standard out with the following fields:
1. AR gene family (from database)
2. AR gene variant (from database)
3. sequence defline/header where the AR gene is located (from genome)
4. % nucleotide similarity (from blastn output)
5. bp length of alignment (from blastn output)
6. bp length of AR gene (from blastn output)

Columns 1 and 2 will have suffixes appended to denote special interest:
- __*__  indicates the best scoring allele is full-length but has >=1 mismatch (SNP). This often means you have a novel allele.
- __?__  indicates uncertainty in the result due to incomplete length alignment
- __TR__ indicates truncation due to an internal stop codon being present

#### II) Raw blastn output (OUTDIR/BASENAME.blastn.tsv)
The tab-delimited outfmt 6 of blast is saved with two columns added to the right. Column 13 is the query length, and column 14 is the subject sequence.

#### III) Log output (OUTDIR/c-SSTAR_BASENAME.log)
A text file is generated to log the date and time of execution, user ID, shell environment, python version, blastn binary location, and blastn version.

## Example Install
    cd $HOME
    pip install biopython
    git clone https://github.com/chrisgulvik/c-SSTAR.git
    echo 'export PATH="$PATH:$HOME/c-SSTAR"' >> $HOME/.bash_profile    

## Example Usage
###### Run c-SSTAR on several genomes with both databases
`cd ~/genomes && for F in *.fna; do B=$(basename $F .fna); B1="$B"_ARG-ANNOT; B2="$B"_ResFinder; python ~/c-SSTAR/c-SSTAR.py -g $F -b $B1 -d ~/AR/ARG-ANNOT.srst2.fasta -o AR/"$B" > AR/"$B"_"$D1".tab; python ~/c-SSTAR/c-SSTAR.py -g $F -b $B2 -d ~/AR/ResFinder_12-14-2015.srst2.fasta -o AR/"$B" > AR/"$B"_"$D2".tab; done`
###### Filter for full-length hits without protein truncations
`for F in *_ARG-ANNOT.tab; do B=$(basename $F _ARG-ANNOT.tab); echo -ne "$B\t" >> Summary_FullLength_AR_hits.tab; grep -v -e $'TR\t' -e '[Oo]mp' $F | awk '$5 == $6 {print $1}' | sed 's/[\?\*]//1' | tr '\n' ',' >> Summary_FullLength_AR_hits.tab; echo '' >> Summary_FullLength_AR_hits.tab; done; sed -i 's/,$//g' Summary_FullLength_AR_hits.tab`
###### Filter for truncated porins
`for F in *_ARG-ANNOT.tab; do B=$(basename $F _ARG-ANNOT.tab); echo -ne "$B\t" >> Summary_TruncatedPorins_AR_hits.tab; grep -P 'TR\t' $F | awk '/[Oo]mp/ && /TR\t/ {print $1}' | sed 's/[\?\*]//1' | tr '\n' ',' >> Summary_TruncatedPorins_AR_hits.tab; echo '' >> Summary_TruncatedPorins_AR_hits.tab; done;
sed -i 's/TR$//g' Summary_TruncatedPorins_AR_hits.tab; sed -i 's/,$//g' Summary_TruncatedPorins_AR_hits.tab`

## Example Summary Output
|AR_Family | AR_Variant | Query_Defline/Header | Similarity | Align_Len | DB_Gene_Len|
|--------------------|---------|-----------------|-----------|---------|--------------|
|AmpH_Bla__*__ | AmpH__*__ | SAMN04014950_3 | 98% | 1161 | 1161|
|Aac3-Ib_AGly__?__ | Aac3-Ib-Aac6-Ib__?__ | SAMN04014950_94 | 99% | 554 | 1005|
|MphA_MLS | MphA | SAMN04014950_85 | 100% | 906 | 906|
|SHV-OKP-LEN_Bla | SHV-11 | SAMN04014950_2 | 100% | 861 | 861|
|CTX-M-1_Bla | CTX-M-15 | SAMN04014950_84 | 100% | 876 | 876|
|Aac3-IIa_AGly__*__ | Aac3-IIa__*__ | SAMN04014950_122 | 98% | 861 | 861|
|TEM-1D_Bla__*__ | TEM-105__*__ | SAMN04014950_112 | 99% | 861 | 861|
|SulI_Sul | SulI | SAMN04014950_91 | 100% | 840 | 840|
|StrB_AGly__*__ | StrB__*__ | SAMN04014950_69 | 99% | 837 | 837|
|AadA_AGly__*__ | AadA2__*__ | SAMN04014950_91 | 99% | 780 | 780|
|OXA-1_Bla | OXA-1 | SAMN04014950_94 | 100% | 831 | 831|
|SulII_Sul | SulII | SAMN04014950_69 | 100% | 816 | 816|
|StrA_AGly__*__ | StrA__*__ | SAMN04014950_69 | 99% | 804 | 804|
|QnrB_Flq__?__ | QnrB1__?__ | SAMN04014950_79 | 100% | 662 | 681|
|CatA2_Phe__*__ | CatA2__*__ | SAMN04014950_101 | 96% | 642 | 642|
|CatBx_Phe__?__ | CatB4__?__ | SAMN04014950_94 | 100% | 529 | 549|
|DfrA_Tmt | DfrA12 | SAMN04014950_118 | 100% | 498 | 498|
|DfrA5_Tmt | DfrA14 | SAMN04014950_116 | 100% | 474 | 474|
|OmpK35-kpneumoniae__*__ | ompK35__*__ | SAMN04014950_14 | 99% | 1080 | 1080|
|OmpK36-kpneumoniae__*TR__ | ompK36__*TR__ | SAMN04014950_1 | 99% | 1099 | 1099|

### Literature References
_Tom's SSTAR paper:_ de Man TJB, Limbago BM. 2016. SSTAR, a stand-alone easy-to-use antimicrobial resistance gene predictor. mSphere 1(1): e00050-15. doi: 10.1128/mSphere.00050-15

_the ARG-ANNOT database:_ Gupta SK, Padmanabhan BR, Diene SM, Lopez-Rojas R, Kempf M, Landraud L, Rolain J-M. 2014. ARG-ANNOT (Antibiotic Resistance Gene-ANNOTation), a new bioinformatic tool to discover antibiotic resistance genes in bacterial genomes. Antimicrobial Agents and Chemotherapy 58:212–220. doi: 10.1128/AAC.01310-13

_the ResFinder database:_ Zankari E, Hasman H, Cosentino S, Vestergaard M, Rasmussen S, Lund O, Aarestrup F, Larsen MV. 2012. Identification of acquired antimicrobial resistance genes. Journal of Antimicrobial Chemotherapy 67:2640–2644. doi: 10.1093/jac/dks261
