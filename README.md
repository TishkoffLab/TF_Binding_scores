## Transcription Factor Binding Scores

This repository contains the scripts that are used to calculate the binding affinity of a given sequence, or of a given Reference/Alternative allele found with GWAS.

If given a snp, the scripts will use the set of BED files containing TF start/stop binding sites from the JASPAR database to find the TF genes that the SNP falls within. This script is get_PWMscore_genenames.sh

Then, using the downloaded JASPAR database of Position Weight Matrix files (transfac formatted), it will calculate the score for each TF, taking the sequence from the hg19 reference chromosome and replacing the base at the target SNP position with the REF and ALT alleles. The output is a tab-seperated file with a row for each TF that the SNP falls within, and the following columns: 
TF_Name, PWM Fraction Score (REF allele), PWM Fraction Score (ALT allele), TF Length, TF Counts per position, H (REF), Hprime (ALT)