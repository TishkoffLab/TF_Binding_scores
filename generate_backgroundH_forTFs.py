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
import pdb
import math

parser = ArgumentParser()
# parser.add_argument("-i", "--input_genes", dest="input_gene_file",
# 					help="input file containing the list of TF names, one per row")
parser.add_argument("-m", "--matrix_loc", dest="matrix_loc",
					help="full path of the folder that contains the PWM matrix files")
parser.add_argument("-o", "--outname", dest="outname",
					help="the name of the file to save the sequence scores to")
parser.add_argument("-f", "--refchrmfastaloc", dest="ref_fasta_loc",
					help="location of the reference fasta files (should just be a fasta per chromosome) to use for getting the reference sequences")
parser.add_argument("-b", "--bg_frac_file", dest="bgfrac_file",
                    help="file containing the background frequency of A/C/T/G, for each autosomal chromosome.")
parser.add_argument("-r", "--reps", dest="reps",
                    help="number of replicate background binding scores to generate")
parser.add_argument("-g", "--genomesize", dest="genome_size_file",
                    help="file containing the chromosme names and sizes, to pull random sequences from")

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


def get_fraclnPWM_from_matrixdict(matrix_dict):
    lnPWM_dict = {}
    for en in range(1,matrix_dict['TF_len']+1):
        temp_matrix = {}
        curr_totcount = sum([float(x) for x in matrix_dict[en].values()])
        for b in 'ACTG':
            f = float(matrix_dict[en][b])/curr_totcount
            if(f == 0.0):
                temp_matrix[b] = np.log(1)
            else:
                temp_matrix[b] = np.log(f)
        lnPWM_dict[en] = temp_matrix
    return lnPWM_dict

def get_lnPWM_from_fracPWM(fracPWM,bgfreqs):
    lnPWM_dict = {}
    for en in range(1,len(fracPWM)+1):
        temp_matrix = {}
        for b in 'ACTG':
            f = float(fracPWM[en][b])/bgfreqs['frac_{0}'.format(b)].values[0]
            if(f == 0.0):
                print('fraction is 0')
                temp_matrix[b] = np.log(1)
            else:
                temp_matrix[b] = -np.log(f)
        lnPWM_dict[en] = temp_matrix
    return lnPWM_dict

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

def get_matrix_counts(pwm,seqlen):
    pos_counts = []
    for n in range(1,seqlen+1):
        temp = [float(x) for x in pwm[n].values()]
        pos_counts.append(sum(temp))
    return pos_counts

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


def calculate_bgH(seq,ln_pwm,bgfreqs):
    currscore_ln = get_matrix_scores(ln_pwm,seq)
    Y = compute_Y(seq,bgfreqs)
    H = currscore_ln - Y
    return np.sum(H)


def get_random_bgseqs(slen,reps,fastaloc,chrmsizes):
    # x = pybedtools.BedTool()
	chrms_touse = list(chrmsizes['chrom'])
	bgseq_list = []
	while(len(bgseq_list) < int(reps)):
		is_valid_seq = True
		curr_chrm = random.randint(1,22)
		curr_start = random.randint(1,chrmsizes.loc[chrmsizes['chrom'] == 'chr{0}'.format(curr_chrm)]['size'].values[0])
		curr_end = curr_start + slen
		chrmfile = '{0}/chr{1}.fa'.format(fastaloc,curr_chrm)
		curr_seq = pybedtools.BedTool.seq('chr{0}:{1}-{2}'.format(curr_chrm,curr_start,curr_end),chrmfile)
		for b in curr_seq:
			if(b not in 'ACTG'):
				is_valid_seq = False
				continue
		if(is_valid_seq):
			bgseq_list.append({'chrm':curr_chrm,'start':curr_start,'end':curr_end,'seq':curr_seq.upper()})
	return bgseq_list

if __name__ == "__main__":
    args = parser.parse_args()
    bgfrac_df = read_csv(args.bgfrac_file,delimiter='\t')
    chrmsizes_df = read_csv(args.genome_size_file,delimiter='\t')
    transfac_matrix_list = os.listdir(args.matrix_loc)
    outfile = open(args.outname,'w')
    outfile.write('Average Background H score for each TF. Number of replicates: {0}\nTF_name\tBG Z score\n'.format(args.reps))
    bg_H_by_TF = {}
    for f in transfac_matrix_list:
        curr_matrix = read_JASPAR_transfac_pfms('{0}/{1}'.format(args.matrix_loc,f))
        bgseqs = get_random_bgseqs(curr_matrix['TF_len'],args.reps,args.ref_fasta_loc,chrmsizes_df)
        bg_H_list = []
        for s in bgseqs:
            curr_seq = s['seq']
            bgfreqs = bgfrac_df.loc[bgfrac_df['Chrm'] == str(s['chrm'])][['frac_A','frac_C','frac_G','frac_T']]
            curr_fracPWM = get_fracPWM_from_matrixdict(curr_matrix)
            curr_lnfracPWM = get_lnPWM_from_fracPWM(curr_fracPWM,bgfreqs)
            curr_H = np.sum(get_matrix_scores(curr_lnfracPWM,curr_seq))
            bg_H_list.append(curr_H)
        curr_z = [pow(math.e,-x) for x in bg_H_list]
        # bg_H_by_TF[curr_matrix['Matrix_Name']] = sum(curr_z)
        outfile.write('{0}\t{1}\n'.format(curr_matrix['Matrix_Name'],sum(curr_z)))
        # bg_H_by_TF[curr_matrix['Matrix_Name']] = bg_H_list
    
    
    # for tf,h in bg_H_by_TF.items():
        
    outfile.close()



