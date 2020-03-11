#!/bin/bash

set -e

snp_file=""
script_path='/home/hautj/TF_binding'
outname=""
genomes_loc="/home/hautj/TF_binding/hg19_refgenomes"
tfs_touse=""

print_usage() {
  printf "Usage: ..."
}

while getopts 'o:s:p:g:t:' flag; do
  case "${flag}" in
    o) outname="${OPTARG}" ;;
    s) snp_file="${OPTARG}" ;;
    p) script_path="${OPTARG}" ;;
    g) genomes_loc="${OPTARG}" ;;
    t) tfs_touse="${OPTARG}" ;;
    *) print_usage
       exit 1 ;;
  esac
done

argflags_string=''

if [[ $tfs_touse != "" ]]; then
    argflags_string+='-t "${tfs_touse}"'
fi

echo "argflags used: ${argflags_string}"

while read line; do
  curr_pos=(`echo "${line}" | cut -d$'\t' -f 3`)
  curr_chrm=(`echo "${line}" | cut -d$'\t' -f 2`)
  curr_ref=(`echo "${line}" | cut -d$'\t' -f 4`)
  curr_alt=(`echo "${line}" | cut -d$'\t' -f 5`)
  curr_outname="${outname}.chrm${curr_chrm}.pos${curr_pos}"
  python ${script_path}/get_PWMscores.py -o "${curr_outname}.above_threshold.PWM_scores" -m "${script_path}/JASPAR2020_CORE_vertebrates_non-redundant_pfms_transfac" -p "${curr_chrm}:${curr_pos}" -r "${curr_ref}" -a "${curr_alt}" -c "${genomes_loc}/chr${curr_chrm}.fa" -b "${script_path}/ACTG_count.all_chrms.fractions.txt" -z "${script_path}/backgroundZ_forTFs.1000reps.txt" -f "0.01" #"${argflags_string}"
done <"${snp_file}"

