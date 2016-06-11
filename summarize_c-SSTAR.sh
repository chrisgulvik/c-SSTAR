#!/bin/sh


# Filter for full-length hits without protein truncations
for F in *_ARG-ANNOT.tab; do
	B=$(basename $F _ARG-ANNOT.tab);
	echo -ne "$B\t" >> Summary_FullLength_AR_hits.tab;
	grep -v -e $'TR\t' -e '[Oo]mp' $F | awk '$5 == $6 {print $2}' | sed 's/[\?\*]//1' | tr '\n' ',' >> Summary_FullLength_AR_hits.tab;
	echo '' >> Summary_FullLength_AR_hits.tab;
done
sed -i 's/,$//g' Summary_FullLength_AR_hits.tab

# Filter for truncated porins
for F in *_ARG-ANNOT.tab; do
	B=$(basename $F _ARG-ANNOT.tab);
	echo -ne "$B\t" >> Summary_TruncatedPorins_AR_hits.tab;
	grep -P 'TR\t' $F | awk '/[Oo]mp/ && /TR\t/ {print $2}' | sed 's/[\?\*]//1' | tr '\n' ',' >> Summary_TruncatedPorins_AR_hits.tab;
	echo '' >> Summary_TruncatedPorins_AR_hits.tab;
done
sed -i 's/TR$//g' Summary_TruncatedPorins_AR_hits.tab; sed -i 's/,$//g' Summary_TruncatedPorins_AR_hits.tab
