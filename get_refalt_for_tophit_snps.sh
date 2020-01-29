#!/bin/bash

set -e

snp_file=""

assoc_file=""

print_usage() {
  printf "Usage: ..."
}

while getopts 's:a:' flag; do
  case "${flag}" in
    s) snp_file="${OPTARG}" ;;
  	a) assoc_file="${OPTARG}" ;;
    *) print_usage
       exit 1 ;;
  esac
done

# assoc_loc="/local3/mhansen/Projects/CTA/Analysis/GWAS/All_Inds/Results/5M_Filtered/${assoc_file}"
assoc_loc="/local3/mhansen/Projects/CTA/Analysis/GWAS/All_Inds/Results/Meta/${assoc_file}"
head -2  "${assoc_loc}" | tail -1 > "${snp_file}.full_assoc"
while read line; do

  curr_snp=(`echo "${line}" | cut -d$':' -f 2`)
  curr_chrm=(`echo "${line}" | cut -d$':' -f 1`)
  awk -v pos="$curr_snp" -v chrm="$curr_chrm" 'BEGIN{OFS="\t"}($2 == pos && $1 == chrm) {print $0;exit}'  "${assoc_loc}" >> "${snp_file}.full_assoc"
  # awk -v pos="$curr_snp" -v chrm="$curr_chrm" 'BEGIN{OFS="\t"}($2 == pos && $1 == chrm) {print $3,$13,$14;exit}'  "${assoc_loc}" >> "${snp_file}.plus_alleles"
done <"${snp_file}"