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

parser = ArgumentParser()
parser.add_argument("-i", "--input_genes", dest="input_gene_file",
					help="input file containing the list of TF gene names, one per row")
parser.add_argument("-s", "--sequence", dest="sequence",
					help="sequence to compute score for, A/C/T/G")
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

def get_bgfracs(bgfrac_file):
    # pdb.set_trace()
    bgfrac_df = read_csv(bgfrac_file,delimiter='\t')
    for row,b in bgfrac_df.iterrows():
        curr_tot = sum([b['count_A'],b['count_C'],b['count_G'],b['count_T']])
        bgfrac_df.at[row,'total'] = curr_tot
        bgfrac_df.at[row,'frac_A'] = int(b['count_A'])/curr_tot
        bgfrac_df.at[row,'frac_C'] = int(b['count_C'])/curr_tot
        bgfrac_df.at[row,'frac_G'] = int(b['count_G'])/curr_tot
        bgfrac_df.at[row,'frac_T'] = int(b['count_T'])/curr_tot
    return bgfrac_df

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

def get_hg19reference_sequence(pos_start,pos_end,filename):
    ref_fasta_file = open(filename,'r')
    seq = ''
    for line in ref_fasta_file:
        if(line[0] == '>'):
            pass
        else:
            seq = ''.join([seq,line.split('\n')[0]])
        if(len(seq) > pos_end):
            break
    ref_fasta_file.close()
    return seq[pos_start:pos_end]


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

def compute_Y(seq,bg_props):
    # pdb.set_trace()
    bpcounts = {'A':0,'C':0,'G':0,'T':0}
    for b in seq:
        try:
            bpcounts[b] += 1
        except:
            if(b not in 'ACTG'):
                print('Sequence contains a letter, {0}, that is not A/C/G/T!'.format(b))
    bpc_list = []
    for o in 'ACGT':
        bpc_list.append(bpcounts[o])
    Y = np.dot(np.log(bg_props),bpc_list)
    return Y

# def get_frac_matrix_scores(tfname,pwm,ref_al,alt_al,position,reference_seq):
#     ref_al = ref_al.upper()
#     alt_al = alt_al.upper()
#     reference_seq = reference_seq.upper()
#     ref_seqval_list,alt_seqval_list = [],[]
#     for n,b in enumerate(reference_seq):
#         try:
#             if(n == position):
#                 ref_seqval_list.append(float(pwm[n+1][ref_al]))
#                 alt_seqval_list.append(float(pwm[n+1][alt_al]))
#             else:
#                 ref_seqval_list.append(float(pwm[n+1][b]))
#                 alt_seqval_list.append(float(pwm[n+1][b]))
#         except:
#             if(b not in 'ACTG'):
#                 print('Sequence contains a letter, {0}, that is not A/C/G/T at position {1}'.format(b,n))
#                 return None
#             else:
#                 continue

    return ref_seqval_list,alt_seqval_list

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

if __name__ == "__main__":
    args = parser.parse_args()

    position = int(args.position.split(':')[1])
    chromosome = int(args.position.split(':')[0])
    # bgfrac_df = get_bgfracs(args.bgfrac_file)
    bgfrac_df = read_csv(args.bgfrac_file,delimiter='\t')

    gene_df = read_csv('{0}'.format(args.input_gene_file),header=None,delimiter='\t')
    gene_df.columns = ['pos_start','pos_end','tf_name','strand']
    transfac_matrix_list = os.listdir(args.matrix_loc)
    infodicts_list = []
    for f in transfac_matrix_list:
        curr_infodict = read_JASPAR_transfac_pfms('{0}/{1}'.format(args.matrix_loc,f))
        infodicts_list.append(curr_infodict)

    ref_pos_end = max(gene_df['pos_end'])
    ref_pos_start = min(gene_df['pos_start'])
    # ref_full_seq = get_hg19reference_sequence(ref_pos_start,ref_pos_end,args.ref_fasta_file)
    ref_full_seq = pybedtools.BedTool.seq('chr{0}:{1}-{2}'.format(chromosome,ref_pos_start,ref_pos_end),args.ref_fasta_file)
    updated_pos = (position-ref_pos_start)
    gene_df['relative_start'] =  gene_df['pos_start']-ref_pos_start
    gene_df['relative_end'] = gene_df['pos_end']-ref_pos_start

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
            bgfreqs = bgfrac_df.loc[bgfrac_df['Chrm'] == chromosome][['frac_A','frac_C','frac_G','frac_T']]
        except:
            bgfreqs = bgfrac_df.loc[bgfrac_df['Chrm'] == 'Total'][['frac_A','frac_C','frac_G','frac_T']]

        curr_fracPWM = get_fracPWM_from_matrixdict(curr_matrix)
        curr_lnfracPWM = get_lnPWM_from_fracPWM(curr_fracPWM,bgfreqs)
        # curr_refscore,curr_altscore,curr_counts = get_matrix_scores(tfname=g['tf_name'],matrix_dict=curr_matrix,ref_al=args.ref_al,alt_al=args.alt_al,position=curr_relative_pos,reference_seq=curr_seq)
        # curr_totcount = sum(curr_counts)
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
        # curr_scoredict['ref_fraction_ln_score'] = sum(refscore_ln)
        # curr_scoredict['alt_fraction_ln_score'] = sum(altscore_ln)
        
        # ref_Y = compute_Y(curr_refseq,bgfreqs)
        # alt_Y = compute_Y(curr_altseq,bgfreqs)
        # ref_H = refscore_ln - ref_Y
        # alt_H = altscore_ln - alt_Y
        # curr_scoredict['H'] = np.sum(ref_H)
        # curr_scoredict['Hprime'] = np.sum(alt_H)

        score_dict_bytf[g['tf_name']] = curr_scoredict

    # score_dict_bytf ={}
    # for i,g in gene_df.iterrows():
    #     curr_relative_pos = abs(updated_pos-g['relative_start'])
    #     curr_seq = ref_full_seq[g['relative_start']:(g['relative_end']+1)]
    #     curr_matrix = get_matrix_byTF(g['tf_name'],infodicts_list)
    #     curr_fracPWM = get_fracPWM_from_matrixdict(curr_matrix)
    #     curr_lnfracPWM = get_fraclnPWM_from_matrixdict(curr_matrix)
    #     curr_refscore,curr_altscore,curr_counts = get_matrix_scores(tfname=g['tf_name'],matrix_dict=curr_matrix,ref_al=args.ref_al,alt_al=args.alt_al,position=curr_relative_pos,reference_seq=curr_seq)
    #     curr_totcount = sum(curr_counts)
    #     curr_scoredict = {'ref_score':curr_refscore,'alt_score':curr_altscore,'ref_fraction_score':(curr_refscore/curr_totcount),'alt_fraction_score':(curr_altscore/curr_totcount),'tf_len':len(curr_seq),'counts_perpos':curr_counts[0]}
    #     frac_refvals,frac_altvals = get_frac_matrix_scores(tfname=g['tf_name'],pwm=curr_fracPWM,ref_al=args.ref_al,alt_al=args.alt_al,position=curr_relative_pos,reference_seq=curr_seq)
    #     fracln_refvals,fracln_altvals = get_frac_matrix_scores(tfname=g['tf_name'],pwm=curr_lnfracPWM,ref_al=args.ref_al,alt_al=args.alt_al,position=curr_relative_pos,reference_seq=curr_seq)
    #     ref_Y = compute_Y(curr_seq,)
    #     score_dict_bytf[g['tf_name']] = curr_scoredict


    outfile = open(args.outname,'w')
    outfile.write('Scores for Transcription Factors Containing SNP at {0} on chromosome {1}, as a fraction of the total count \nTF_Name\tPWM Fraction Score (REF allele)\tPWM Fraction Score (ALT allele)\tTF Length\tTF Counts per position\tH (REF)\tHprime (ALT)\n'.format(position,chromosome))

    # outfile.write('Scores for Transcription Factors Containing SNP at {0} on chromosome {1}, as a fraction of the total count \nTF_Name\tPWM Fraction Score (REF allele)\tPWM Fraction Score (ALT allele)\tTF Length\tTF Counts per position\tREF Log(PWM Fraction Score)\tALT Log(PWM Fraction Score)\tH (REF)\tHprime (ALT)\n'.format(position,chromosome))
    # for tf,scores in score_dict_bytf.items():
    for tf,scores in sorted(score_dict_bytf.items(), key=lambda k_v: k_v[1]['alt_fraction_score'],reverse=True):
        outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(tf,scores['ref_fraction_score'],scores['alt_fraction_score'],scores['tf_len'],
            scores['counts_perpos'],scores['H'],scores['Hprime']))
    outfile.close()









