# c-SSTAR
##### CLI-Sequence Search Tool for Antimicrobial Resistance
![alt tag](https://github.com/chrisgulvik/images/raw/master/c-SSTAR.jpeg)


## System Requirements
- Linux or Mac OS X platform
- BLAST+ (blastn and makeblastdb)
- Python 2.7 or 3 with BioPython

## Usage
    c-SSTAR -g <genome_file> -d <database_file>

## Input

1. FastA formatted genome
2. FastA database of antimicrobial resistance (AR) gene sequences from SSTAR. Two databases are available ([ARG-ANNOT](https://github.com/tomdeman-bio/Sequence-Search-Tool-for-Antimicrobial-Resistance-SSTAR-/blob/master/Manuscript_databases/ARG-ANNOT.srst2.fasta) or [ResFinder](https://github.com/tomdeman-bio/Sequence-Search-Tool-for-Antimicrobial-Resistance-SSTAR-/blob/master/Manuscript_databases/Resfinder.srst2.fasta)), which are formatted according to Kat Holt's clustering approach for [SRST2](https://github.com/katholt/srst2/tree/master/database_clustering). A combination of these two databases also exists ["ResGANNOT"](https://github.com/tomdeman-bio/Sequence-Search-Tool-for-Antimicrobial-Resistance-SSTAR-/blob/master/Latest_AR_database/ResGANNOT_srst2.fasta)

## Output
#### I) Summary output (to stdout)
A tab-delimited summary is printed to standard out with the following fields:
1. AR gene family (from database)
2. AR gene variant (from database)
3. sequence defline/header where the AR gene is located (from genome)
4. % nucleotide identity (from blastn output)
5. bp length of alignment (from blastn output)
6. bp length of AR gene (from blastn output)

Columns 1 and 2 will have suffixes appended to denote special interest:
- __*__  indicates the best scoring allele is full-length but has >=1 mismatch (SNP). This often means you have a novel allele.
- __?__  indicates uncertainty in the result due to incomplete length alignment
- __TR__ indicates truncation due to an internal stop codon being present
- __$__ indicates gene detected at edge of contig

#### II) Raw alignment output (OUTDIR/BASENAME.blastn.tsv)
The tab-delimited outfmt 6 of blast is saved with three columns added to the right. Column 13 is the query (AR gene database) length, column 14 is the subject (contig) length, and column 15 is the subject (contig's AR gene) sequence.

#### III) Log output (OUTDIR/c-SSTAR_BASENAME.log)
A text file is generated to log the date and time of execution, user ID, shell environment, python version, blastn binary location, and blastn version.

## Example Install
    pip install biopython
    git clone https://github.com/chrisgulvik/c-SSTAR.git $HOME
    echo 'export PATH="$PATH:$HOME/c-SSTAR"' >> $HOME/.bash_profile

## Example Usage
###### Run c-SSTAR on several genomes with the combo database
```bash
for F in *.fna; do
  B=$(basename $F .fna)
  c-SSTAR -g $F -d ~/c-SSTAR/db/ResGANNOT_srst2.fasta.gz -o $B > "$B"_ResGANNOT.tab
done
```

## Example Summary Output
```bash
c-SSTAR -g ~/c-SSTAR/tests/data/SRR3112344.fa.gz \
 -d ~/c-SSTAR/db/ResGANNOT_srst2.fasta.gz
```
| AR_Family | AR_Variant | Query_Defline | Identity | Aln_Len | DB_Gene_Len |
| --------- | ---------------------- | ------------- | -------- | ------- | ----------- |
| aac(3)* | aac(3)-IId* | tig093 | 99.884% | 861 | 861 |
| aac(3)? | aac(3)-Ib-aac(6')-Ib'? | tig123 | 99.097% | 554 | 1005 |
| aac(6') | aac(6')-Ib-cr | tig123 | 100.0% | 600 | 600 |
| aac(6')? | aac(6')-30-aac(6')-Ib'? | tig123 | 99.309% | 579 | 987 |
| aadA2? | aadA2? | tig104 | 99.875% | 802 | 819 |
| ampH* | ampH* | tig003 | 98.88% | 1161 | 1161 |
| aph(3'')* | aph(3'')-Ib* | tig096 | 99.876% | 804 | 804 |
| aph(6)* | aph(6)-Id* | tig096 | 99.881% | 837 | 837 |
| blaCTX | blaCTX-M-15 | tig089 | 100.0% | 876 | 876 |
| blaOXA | blaOXA-1 | tig123 | 100.0% | 831 | 831 |
| blaSHV | blaSHV-11 | tig016 | 100.0% | 861 | 861 |
| blaSHV* | blaSHV-100* | tig016 | 95.111% | 900 | 900 |
| blaTEM | blaTEM-1B | tig108 | 100.0% | 861 | 861 |
| catA2* | catA2* | tig127 | 96.106% | 642 | 642 |
| catB3?$ | catB3?$ | tig123 | 100.0% | 440 | 633 |
| catB4?$ | catB4?$ | tig123 | 100.0% | 440 | 549 |
| dfrA12 | dfrA12 | tig104 | 100.0% | 498 | 498 |
| dfrA14* | dfrA14* | tig134 | 99.586% | 483 | 483 |
| fosA6?$ | fosA6?$ | tig026 | 98.81% | 420 | 433 |
| mph(A)?TR | mph(A)?TR | tig110 | 99.675% | 922 | 921 |
| oqxA | oqxA | tig024 | 100.0% | 1176 | 1176 |
| oqxB | oqxB | tig024 | 100.0% | 3153 | 3153 |
| qnrB1*TR | qnrB1*TR | tig080 | 99.853% | 681 | 681 |
| sul1* | sul1* | tig104 | 99.885% | 867 | 867 |
| sul2 | sul2 | tig096 | 100.0% | 816 | 816 |
| tet(D) | tet(D) | tig108 | 100.0% | 1185 | 1185 |

### Literature References
[_c-SSTAR_](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5442553/): Cunningham SA, Limbago B, Traczewski M, Anderson K, Hackel M, Hindler J, Sahm D, Alyanak E, Lawsin A, Gulvik CA, de Man TJB, Mandrekar JN, Schuetz AN, Jenkins S, Humphries R, Palavecino E, Vasoo S, Patel R. 2017. Multicenter Performance Assessment of Carba NP Test. J Clin Microbiol 55(6):1954-1960. doi: 10.1128/JCM.00244-17

[_SSTAR_](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4863618/): de Man TJB, Limbago BM. 2016. SSTAR, a stand-alone easy-to-use antimicrobial resistance gene predictor. mSphere 1(1): e00050-15. doi: 10.1128/mSphere.00050-15

[_ARG-ANNOT database_](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3910750/): Gupta SK, Padmanabhan BR, Diene SM, Lopez-Rojas R, Kempf M, Landraud L, Rolain J-M. 2014. ARG-ANNOT (Antibiotic Resistance Gene-ANNOTation), a new bioinformatic tool to discover antibiotic resistance genes in bacterial genomes. Antimicrobial Agents and Chemotherapy 58:212–220. doi: 10.1128/AAC.01310-13

[_ResFinder database_](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3468078/): Zankari E, Hasman H, Cosentino S, Vestergaard M, Rasmussen S, Lund O, Aarestrup F, Larsen MV. 2012. Identification of acquired antimicrobial resistance genes. Journal of Antimicrobial Chemotherapy 67:2640–2644. doi: 10.1093/jac/dks261
