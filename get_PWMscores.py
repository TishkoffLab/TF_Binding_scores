import sys
from pandas import *
import numpy as np
import matplotlib
from matplotlib import pyplot
import random
from scipy.stats import norm
import os
from argparse import ArgumentParser
import pybedtools
# import pdb

parser = ArgumentParser()
parser.add_argument("-i", "--input_genes", dest="input_gene_file",
					help="input file containing the list of TF gene names, one per row")
parser.add_argument("-s", "--sequence", dest="sequence",
					help="sequence to compute score for, A/C/T/G")
parser.add_argument("-t", "--tfname", dest="tf_tocheck",
                    help="name of a specfic transcription factor, or a file containing any number of TFs (one per line). If this argument is supplied, then the script only calculates the score for that TF. Must supply a sequence as well.")
parser.add_argument("-m", "--matrix_loc", dest="matrix_loc",
					help="full path of the folder that contains the PWM matrix files")
parser.add_argument("-o", "--outname", dest="outname",
					help="the name of the file to save the sequence scores to")
parser.add_argument("-r", "--refallele", dest="ref_al",
					help="reference allele for the snp of interest, A/C/T/G")
parser.add_argument("-a", "--altallele", dest="alt_al",
					help="alternate allele for the snp of interest, A/C/T/G")
parser.add_argument("-p", "--position", dest="position",
					help="position, in bp, for the snp of interest")
parser.add_argument("-c", "--refchrmfasta", dest="ref_fasta_file",
					help="reference fasta file (should just be a single chromosome) to use for getting the reference sequence")
parser.add_argument("-b", "--bg_frac_file", dest="bgfrac_file",
                    help="file containing the background frequency of A/C/T/G, for each autosomal chromosome.")

#Reads in the JASPAR PWM file (transfac formatted)
#   infile (str): the PWM file to read in
#Returns:
#   pfms_info (dict): dictionary containig the information about the PWM, from the file (plus an entry with the length of the PWM)
def read_JASPAR_transfac_pfms(infile):
    pfms_file = open(infile,'r')
    pfms_info = {}
    seq_len = 0
    for line in pfms_file:
        line = line.split('\n')[0]
        if(len(line.split(' ')) > 1):
            line = line.split(' ')
            if(line[0] == 'DE'):
                pfms_info['Matrix_ID'] = line[1]
                pfms_info['Matrix_Name'] = line[2]
                seq_len = 0
            elif(line[0] == 'CC'):
                temp = line[1].split(':')
                pfms_info[temp[0]] = temp[1]
                if(seq_len > 0):
                    pfms_info['TF_len'] = seq_len
        elif(len(line.split('\t')) > 1):
            line = line.split('\t')
            if(line[0] == 'PO'):
                curr_matorder = line[1:]
            else:
                curr_vals = {}
                for n,v in enumerate(line[1:]):
                    curr_vals[curr_matorder[n]] = float(v)+1
                pfms_info[int(line[0])] = curr_vals
                seq_len = int(line[0])
        else:
            pass
    pfms_file.close()
    return pfms_info

#Loops through the info_dicts list (of PWM matrix file info), and returns the PWM matrix dict for the given TF
#Inputs:
#   tfname (str): name of the transcription factor
#   info_dicts (list): made by looping over all the JASPAR matrix files; this is a list of all of those matrices as dicts
#Returns:
#   matrix_dict_touse (dict): the dictionary containing the PWM for the given TF
def get_matrix_byTF(tfname,info_dicts):
    matrix_dict_touse = None
    for i in info_dicts:
        if(i['Matrix_Name'] == tfname):
            matrix_dict_touse = i
            break
    if(matrix_dict_touse == None):
        print('Could not find a PWM for Transcription Factor {0}'.format(tfname))
        return None
    return matrix_dict_touse

#Given a matrix dict for a TF, containing the PWM of the counts for each base in the sequence, returns just the PWM with each position entry being recalculated as a fraction of the count at that position
#Inputs:
#   matrix_dict (dict): the dictionary containing the PWM for a given TF, in addition to the other data about that TF
#Returns:
#   PWM_dict (dict): a dicitonary where each key is a position relative to the TF (1-indexed) and each value is a dictionary with keys A/C/G/T and values equal to the raw count divided by the total counts for all four bases at that position.
def get_fracPWM_from_matrixdict(matrix_dict):
    PWM_dict = {}
    for en in range(1,matrix_dict['TF_len']+1):
        temp_matrix = {}
        curr_totcount = sum([float(x) for x in matrix_dict[en].values()])
        for b in 'ACTG':
            if(f == 0.0):
                temp_matrix[b] = 1/curr_totcount
            else:
                temp_matrix[b] = float(matrix_dict[en][b])/curr_totcount
        PWM_dict[en] = temp_matrix
    return PWM_dict

#Given a fractional PWM (made by get_fracPWM_from_matrixdict) and a set of background frequencies of the four bases, calculate a -log PWM
#Inputs:
#   fracPWM (dict): PWM where each entry is a fraction of the counts for each base
#   bgfreqs (DataFrame): a dataframe with a single row, containing columns frac_A/frac_C/frac_G/frac_T which has the fraction of the chromosome/genome corrisponding to that base
#Returns:
#   lnPWM_dict (dict): PWM where each entry is the fracPWM entry, divided by the background base fraction, then taken the negative natural log of it
def get_lnPWM_from_fracPWM(fracPWM,bgfreqs):
    # pdb.set_trace()
    lnPWM_dict = {}
    for en in range(1,len(fracPWM)+1):
        temp_matrix = {}
        for b in 'ACTG':
            f = float(fracPWM[en][b])/bgfreqs['frac_{0}'.format(b)].values[0]
            temp_matrix[b] = -np.log(f)
        lnPWM_dict[en] = temp_matrix
    return lnPWM_dict

#For a given sequence, returns the complementary sequence
#Inputs:
#   seq (str): sequence of A/C/T/G
#Returns:
#   new_seq (str): sequence of A/C/T/G complementary to the original sequence
def get_complseq(seq):
    new_seq = []
    for b in seq:
        if(b == 'A'):
            new_seq.append('T')
        elif(b == 'T'):
            new_seq.append('A')
        elif(b == 'G'):
            new_seq.append('C')
        elif(b == 'C'):
            new_seq.append('G')
        else:
            print('Base pair not A/C/T/G! {0}'.format(b))
    return ''.join(new_seq)

#Given a sequence, replaces the allele at a specified position with the given allele
#Inputs:
#   fullseq (str): sequence of A/C/T/G
#   position (int): position of the allele to be replaced, relative to the length of the input sequence (so it must be <= len(fullseq))
#   allele (str): A/C/T/G, to replace the one at the position in fullseq
#Returns:
#   new_seq (str): fullseq, with the new allele at the given position
def make_seq(fullseq,position,allele):
    new_seq = ''
    for n,b in enumerate(fullseq):
        try:
            if(n==position):
                new_seq = ''.join([new_seq,allele])
            else:
                new_seq = ''.join([new_seq,b])
        except:
            if(b not in 'ACTG'):
                print('Sequence contains a letter, {0}, that is not A/C/G/T at position {1}'.format(b,n))
                return None
            else:
                continue
    return new_seq.upper()

#For a given PWM and sequence length, return the sum of the counts for all four bases at each position
#Inputs:
#   pwm (dict): position weight matrix, with or without the additional info from the matrix file
#   seqlen (int): length of the sequence, so that we can loop over the PWM
#Returns:
#   pos_counts (list): list of counts at each position
def get_matrix_counts(pwm,seqlen):
    pos_counts = []
    for n in range(1,seqlen+1):
        temp = [float(x) for x in pwm[n].values()]
        pos_counts.append(sum(temp))
    return pos_counts

#For a given PWM and sequence, calculates the score at each position
#Inputs:
#   pwm (dict): position weight matrix, with or without the additional info from the matrix file
#   seq (str): A/C/G/T sequence
#Returns:
#   seqval_list (list): list of the values at each position given the specific base in seq
def get_matrix_scores(pwm,seq):
    seqval_list = []
    for n,b in enumerate(seq):
        try:
            seqval_list.append(float(pwm[n+1][b]))
        except:
            if(b not in 'ACTG'):
                print('Sequence contains a letter, {0}, that is not A/C/G/T at position {1}'.format(b,n))
                return None
            else:
                continue
    return seqval_list

if __name__ == "__main__":
    args = parser.parse_args()
    bgfrac_df = read_csv(args.bgfrac_file,delimiter='\t')

    if(args.sequence is not None):
        sequence = args.sequence

        #Read in the matrix files and make dict entries for each one
        transfac_matrix_list = os.listdir(args.matrix_loc)
        infodicts_list = []
        for f in transfac_matrix_list:
            curr_infodict = read_JASPAR_transfac_pfms('{0}/{1}'.format(args.matrix_loc,f))
            infodicts_list.append(curr_infodict)

        score_dict_bytf ={}
        tfs_to_check = []
        try:
            tf_list_file = open(args.tf_tocheck,'r')
            for line in tf_list_file:
                tfs_to_check.append(line.split('\n')[0])
            tf_list_file.close()
        except FileNotFoundError:
            tfs_to_check.append(args.tf_tocheck)

        bgfreqs = bgfrac_df.loc[bgfrac_df['Chrm'] == 'Total'][['frac_A','frac_C','frac_G','frac_T']]
        for tf in tfs_to_check:
            curr_matrix = get_matrix_byTF(tf,infodicts_list)
            curr_fracPWM = get_fracPWM_from_matrixdict(curr_matrix)
            curr_lnfracPWM = get_lnPWM_from_fracPWM(curr_fracPWM,bgfreqs)
            rawscore_list = get_matrix_scores(curr_matrix,sequence)
            pos_counts = get_matrix_counts(curr_matrix,curr_matrix['TF_len'])
            tot_count = sum(pos_counts)
            # pdb.set_trace()
            curr_scoredict = {'raw_score':sum(rawscore_list),'tf_len':curr_matrix['TF_len'],'counts_perpos':min(pos_counts)}
            curr_scoredict['fraction_score'] = sum([rawscore_list[x]/pos_counts[x] for x in range(len(pos_counts))])
            score_ln = get_matrix_scores(curr_lnfracPWM,sequence)
            curr_scoredict['H'] = np.sum(score_ln)
            score_dict_bytf[tf] = curr_scoredict
        #Writing the PWM scores to the output file
        outfile = open(args.outname,'w')
        outfile.write('Scores for Given Transcription Factors for sequence {0}, as a fraction of the total count \nTF_Name\tPWM Fraction Score\tTF Length\tTF Counts per position\tH\n'.format(sequence))
        for tf,scores in sorted(score_dict_bytf.items(), key=lambda k_v: k_v[1]['H'],reverse=True):
            outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(tf,scores['fraction_score'],scores['tf_len'],scores['counts_perpos'],scores['H']))
        outfile.close()
    else:
        position = int(args.position.split(':')[1])
        chromosome = int(args.position.split(':')[0])

        #Reading in the TF_genes file, made by the bash script, which has the TF names, positions, and strandedness
        gene_df = read_csv('{0}'.format(args.input_gene_file),header=None,delimiter='\t')
        gene_df.columns = ['pos_start','pos_end','tf_name','strand']

        #Read in the matrix files and make dict entries for each one
        transfac_matrix_list = os.listdir(args.matrix_loc)
        infodicts_list = []
        for f in transfac_matrix_list:
            curr_infodict = read_JASPAR_transfac_pfms('{0}/{1}'.format(args.matrix_loc,f))
            infodicts_list.append(curr_infodict)

        #Getting the reference sequence that contains all of the TF genes within it, then add the new start/stop coordinates relative to full ref seq to the dataframe
        ref_pos_end = max(gene_df['pos_end'])
        ref_pos_start = min(gene_df['pos_start'])
        ref_full_seq = pybedtools.BedTool.seq('chr{0}:{1}-{2}'.format(chromosome,ref_pos_start,ref_pos_end),args.ref_fasta_file)
        updated_pos = (position-ref_pos_start)
        gene_df['relative_start'] =  gene_df['pos_start']-ref_pos_start
        gene_df['relative_end'] = gene_df['pos_end']-ref_pos_start

        #Creating the final dictionary containing the values for each TF, with the Ref/Alt alleles
        score_dict_bytf ={}
        for i,g in gene_df.iterrows():
            curr_relative_pos = abs(updated_pos-g['relative_start'])
            curr_refseq = make_seq(ref_full_seq[g['relative_start']:(g['relative_end']+1)],curr_relative_pos,args.ref_al)
            curr_altseq = make_seq(ref_full_seq[g['relative_start']:(g['relative_end']+1)],curr_relative_pos,args.alt_al)
            if(g['strand'] == '-'):
                curr_refseq = get_complseq(curr_refseq)
                curr_altseq = get_complseq(curr_altseq)
            curr_matrix = get_matrix_byTF(g['tf_name'],infodicts_list)
            try:
                bgfreqs = bgfrac_df.loc[bgfrac_df['Chrm'] == str(chromosome)][['frac_A','frac_C','frac_G','frac_T']]
            except:
                bgfreqs = bgfrac_df.loc[bgfrac_df['Chrm'] == 'Total'][['frac_A','frac_C','frac_G','frac_T']]
            # pdb.set_trace()
            curr_fracPWM = get_fracPWM_from_matrixdict(curr_matrix)
            curr_lnfracPWM = get_lnPWM_from_fracPWM(curr_fracPWM,bgfreqs)
            refscore_list = get_matrix_scores(curr_matrix,curr_refseq)
            altscore_list = get_matrix_scores(curr_matrix,curr_altseq)
            pos_counts = get_matrix_counts(curr_matrix,curr_matrix['TF_len'])
            tot_count = sum(pos_counts)
            curr_scoredict = {'ref_score':sum(refscore_list),'alt_score':sum(altscore_list),'tf_len':curr_matrix['TF_len'],'counts_perpos':min(pos_counts)}
            curr_scoredict['ref_fraction_score'] = sum([refscore_list[x]/pos_counts[x] for x in range(len(pos_counts))])
            curr_scoredict['alt_fraction_score'] = sum([altscore_list[x]/pos_counts[x] for x in range(len(pos_counts))])

            refscore_ln = get_matrix_scores(curr_lnfracPWM,curr_refseq)
            altscore_ln = get_matrix_scores(curr_lnfracPWM,curr_altseq)
            curr_scoredict['H'] = np.sum(refscore_ln)
            curr_scoredict['Hprime'] = np.sum(altscore_ln)

            score_dict_bytf[g['tf_name']] = curr_scoredict

        #Writing the PWM scores to the output file
        outfile = open(args.outname,'w')
        outfile.write('Scores for Transcription Factors Containing SNP at {0} on chromosome {1}, as a fraction of the total count \nTF_Name\tPWM Fraction Score (REF allele)\tPWM Fraction Score (ALT allele)\tTF Length\tTF Counts per position\tH (REF)\tHprime (ALT)\n'.format(position,chromosome))
        for tf,scores in sorted(score_dict_bytf.items(), key=lambda k_v: k_v[1]['alt_fraction_score'],reverse=True):
            outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(tf,scores['ref_fraction_score'],scores['alt_fraction_score'],scores['tf_len'],
                scores['counts_perpos'],scores['H'],scores['Hprime']))
        outfile.close()
    










