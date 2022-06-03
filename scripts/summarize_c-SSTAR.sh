#!/bin/sh


usage() {
  echo "
Usage: ${0##*/} [FileExtension] [-h|--help]

Summarizes tab-delimited output files from c-SSTAR, which creates
two output files:
  Summary_FullLength_AR_hits.tab       :  full alignments lacking a predicted internal stop codon
  Summary_TruncatedPorins_AR_hits.tab  :  alignments with predicted with internal stop codon

Optional:
  -h | --help      Show this help message and exit.
  <FileExtension>  File extension of c-SSTAR files to summarize
                   Default: _ResGANNOT.tab

Must be ran in the working directory containing files to summarize
and output files will be placed in the same path.
"
}

[[ "${1}" == '--help' || "${1}" == '-h' ]] && { usage; exit 0; }
if [[ -z "${1}" ]]; then
  EXT='_ResGANNOT.tab'
else
  EXT="${2}"
fi

shopt -s nullglob
FILES=( *"${EXT}" )
shopt -u nullglob
if [ "${#FILES[@]}" -lt 1 ]; then
  echo 'ERROR: c-SSTAR output files found to summarize' >&2
  exit 1
fi

# Filter for full-length hits without protein truncations
echo -ne '' > Summary_FullLength_AR_hits.tab
for F in "${FILES[@]}"; do
  B=$(basename "${F}" "${EXT}")
  echo -ne "${B}\t" >> Summary_FullLength_AR_hits.tab
  grep -v -e $'TR\t' -e '[Oo]mp' "${F}" |\
   awk '$5 == $6 {print $2}' |\
   sed 's/[\?\*]//1' |\
   tr '\n' ',' >> Summary_FullLength_AR_hits.tab
  echo '' >> Summary_FullLength_AR_hits.tab
done
sed -i 's/,$//g' Summary_FullLength_AR_hits.tab

# Filter for truncated porins
echo -ne '' > Summary_TruncatedPorins_AR_hits.tab
for F in "${FILES[@]}"; do
  B=$(basename "${F}" "${EXT}")
  echo -ne "${B}\t" >> Summary_TruncatedPorins_AR_hits.tab
  grep -P 'TR\t' "${F}" |\
   awk '/[Oo]mp/ && /TR\t/ {print $2}' |\
   sed 's/[\?\*]//1' |\
   tr '\n' ',' >> Summary_TruncatedPorins_AR_hits.tab
  echo '' >> Summary_TruncatedPorins_AR_hits.tab
done
sed -i 's/TR$//g' Summary_TruncatedPorins_AR_hits.tab
sed -i 's/,$//g' Summary_TruncatedPorins_AR_hits.tab
