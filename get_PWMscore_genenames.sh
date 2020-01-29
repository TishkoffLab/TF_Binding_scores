#!/bin/bash

set -e

tf_name=""
script_path='/home/hautj/TF_binding'
outname=""
chrm=""
pos_tocheck=""
seq=""
ref_al=""
alt_al=""

print_usage() {
  printf "Usage: ..."
}

while getopts 'o:c:t:p:q:r:a:' flag; do
  case "${flag}" in
    o) outname="${OPTARG}" ;;
    c) chrm="${OPTARG}" ;;
    t) tf_name="${OPTARG}" ;;
  	p) pos_tocheck="${OPTARG}" ;;
  	q) seq="${OPTARG}" ;;
  	r) ref_al="${OPTARG}" ;;
  	a) alt_al="${OPTARG}" ;;
    *) print_usage
       exit 1 ;;
  esac
done

orig_bedfile="JASPAR2020_hg19.converted.chrm${chrm}.bed"

awk -v pos="$pos_tocheck" 'BEGIN{OFS="\t"}($2 <= pos && $3 >= pos) {print $2,$3,$4} $2 > pos {exit}'  "${script_path}/JASPAR2020_hg19_bedfiles_bychrm/${orig_bedfile}" > "${script_path}/${outname}.TF_genes"
echo "found TF genes, starting score calculation"
python ${script_path}/get_PWMscores.py -o "${script_path}/${outname}.PWM_scores" -i "${script_path}/${outname}.TF_genes" -m "${script_path}/JASPAR2020_CORE_vertebrates_non-redundant_pfms_transfac" -p "${pos_tocheck}" -r "${ref_al}" -a "${alt_al}" -c "${script_path}/hg19_refgenomes/chr${chrm}.fa"

