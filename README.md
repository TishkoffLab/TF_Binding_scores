## Transcription Factor Binding Scores

This repository contains the scripts that are used to calculate the binding affinity of a given sequence, or of a given Reference/Alternative allele found with GWAS.

If given a snp, the scripts will use the set of BED files containing TF start/stop binding sites from the JASPAR database to find the TF genes that the SNP falls within. This script is get_PWMscore_genenames.sh

Then, using the downloaded JASPAR database of Position Weight Matrix files (transfac formatted), it will calculate the score for each TF, taking the sequence from the hg19 reference chromosome and replacing the base at the target SNP position with the REF and ALT alleles. The output is a tab-seperated file with a row for each TF that the SNP falls within, and the following columns: 
TF_Name, PWM Fraction Score (REF allele), PWM Fraction Score (ALT allele), TF Length, TF Counts per position, H (REF), Hprime (ALT)

### Scripts:


get_PWMscore.py

This script takes a particular SNP of interest, the list of TF names that corrisponds to that SNP, and calculates the scores from the PWMs of each TF.
```
Input flags:
-i --input_genes
				input file containing the list of TF gene names, one per row
-m --matrix_loc
				full path of the folder that contains the PWM matrix files
-o --outname
				the name of the file to save the sequence scores to
-r --refallele
				reference allele for the snp of interest, A/C/T/G
-a --altallele
				alternate allele for the snp of interest, A/C/T/G
-p --position
				position, in bp, for the snp of interest
-c --refchrmfasta
				reference fasta file (should just be for a single chromosome) to use for getting the reference sequence
-b --bg_frac_file
				file containing the background frequency of A/C/T/G, for each autosomal chromosome
-z --bg_zscore_file
				file containing the background Z scores for each TF, precalculated using a significant number of replicates
-f --tf_cutoff
				the cutoff for significant pbinding score. If this is provided, the script will check the snp against all tfs and orientations, then save only the results that are above the threshold.

Note: bg_frac_file was made using a seperate script, with frequencies calculated from reference genome hg_19 and is in this repository as "ACTG_count.all_chrms.fractions.txt".

```
get_PWMscore_genenames_multisnpfile.sh

This script runs the python script on each of the SNPs that we have found through GWAS. It finds the TFs from the premade JASPAR bedfile that overlap with each particular SNP and makes a file containing that info, and then runs get_PWMscore.py with the requisite input flags.

```
Input flags:
-o outname
			Name to use for the output files (<outname>.TF_genes and <outname>.PWM_scores)

-s snp_file
			File containing the SNPs that will be analyzed. First line is the column names, followed by one line per SNP. Each line must have the SNP name (chr:pos format), the reference allele, and the alternative allele
-p script_path
			The path where the script is located. This folder must contain get_PWMscore.py, the folder containing 

```
