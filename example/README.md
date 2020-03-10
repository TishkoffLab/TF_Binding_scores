## Running the example

For the standard operation:

```
./get_PWMscore_genenames_multisnpfile.sh -s example/example_SNPs.input\
	 -o example/example_SNPs.test.output\
	 -b example/JASPAR2020_hg19.converted.subset.example \
	 -p ${PWD} \
	 -g ${PWD}/example 

```

Using the functionality of just a single sequence:

```
python get_PWMscores.py -s "AAGACATTTGAAAATTATCTA"\
	 -o example/example_seq.test.output\
	 -m ${PWD}/JASPAR2020_CORE_vertebrates_non-redundant_pfms_transfac \
	 -z ${PWD}/backgroundZ_forTFs.1000reps.txt \
	 -b ${PWD}/ACTG_count.all_chrms.fractions.txt
```