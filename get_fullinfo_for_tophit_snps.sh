#!/bin/bash

set -e

arr=(/local3/jake/TF_binding/tophit_SNPs/*.TopVars.LE)


for f in "${arr[@]}"; do
   curr_traitcode=(`echo ${f} | cut -d'/' -f 6 |cut -d'.' -f 1`)
   ./get_refalt_for_tophit_snps.sh -s "$f" -a "${curr_traitcode}.5M.assoc"
done