#!/bin/bash

set -e

snp_file=""
script_path='/home/hautj/TF_binding'
outname=""
bedfile_prefix="/home/hautj/TF_binding/JASPAR2020_hg19_bedfiles_bychrm/JASPAR2020_hg19.converted"
genomes_loc="/home/hautj/TF_binding/hg19_refgenomes"

print_usage() {
  printf "Usage: ..."
}

while getopts 'o:s:p:b:g:' flag; do
  case "${flag}" in
    o) outname="${OPTARG}" ;;
    s) snp_file="${OPTARG}" ;;
    p) script_path="${OPTARG}" ;;
    b) bedfile_prefix="${OPTARG}" ;;
    g) genomes_loc="${OPTARG}" ;;
    *) print_usage
       exit 1 ;;
  esac
done

lcount=0
while read line; do
  #Skip the first line, since it's the header 
  if [[ "${lcount}" -eq 0 ]]; then
    lcount=1
  else
    curr_pos=(`echo "${line}" | cut -d$'\t' -f 1 | cut -d$':' -f 2`)
    curr_chrm=(`echo "${line}" | cut -d$'\t' -f 1 | cut -d$':' -f 1`)
    curr_ref=(`echo "${line}" | cut -d$'\t' -f 2`)
    curr_alt=(`echo "${line}" | cut -d$'\t' -f 3`)
    curr_outname="${outname}.chrm${curr_chrm}.pos${curr_pos}"
    orig_bedfile="${bedfile_prefix}.chrm${curr_chrm}.bed"
    #Looping through the bedfile (containing the TFs found on that chromosome), and add them to a file.
    awk -v pos="$curr_pos" 'BEGIN{OFS="\t"}($2 <= pos && $3 >= pos) {print $2,$3,$4,$6} $2 > pos {exit}'  "${orig_bedfile}" > "${curr_outname}.TF_genes"
    echo "found TF genes, starting score calculation for snp ${curr_chrm}:${curr_pos}"
    python ${script_path}/get_PWMscores.py -o "${curr_outname}.PWM_scores" -i "${curr_outname}.TF_genes" -m "${script_path}/JASPAR2020_CORE_vertebrates_non-redundant_pfms_transfac" -p "${curr_chrm}:${curr_pos}" -r "${curr_ref}" -a "${curr_alt}" -c "${genomes_loc}/chr${curr_chrm}.fa" -b "${script_path}/ACTG_count.all_chrms.fractions.txt" -z "${script_path}/backgroundZ_forTFs.100reps.txt"
  fi
done <"${snp_file}"

