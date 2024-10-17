# -*- coding: utf-8 -*-
"""
Created on Wed Feb 7 16:53:05 2024

@author: Dan_salv
"""
import os
import numpy as np
import math
import configparser
config = configparser.ConfigParser()

module_dir = os.path.dirname(os.path.realpath(__file__))
config.read(f'{module_dir}/config.ini')

MATRIX_DIR = config['DEFAULT']['MATRIX_DIR']
BLOSUM = config['FILEPATH']['BLOSUM']

matrix_filepath = os.path.join(MATRIX_DIR,BLOSUM)

PSEUDOCOUNT = .0000001

amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '-']
iupac_alphabet = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "U", "V", "W", "Y", "Z", "X", "*", "-"] 

# dictionary to map from amino acid to its row/column in a similarity matrix
aa_to_index = {}
for i, aa in enumerate(amino_acids):	
    aa_to_index[aa] = i

blosum_background_distr = [0.078, 0.051, 0.041, 0.052, 0.024, 0.034, 0.059, 0.083, 0.025, 0.062, 0.092, 0.056, 0.024, 0.044, 0.043, 0.059, 0.055, 0.014, 0.034, 0.072]


# BLOSUM62 background distribution
blosum_background_distr = [0.078, 0.051, 0.041, 0.052, 0.024, 0.034, 0.059, 0.083, 0.025, 0.062, 0.092, 0.056, 0.024, 0.044, 0.043, 0.059, 0.055, 0.014, 0.034, 0.072]

class CONSERVATION_SCORE:    
    
    def __init__(self,
        msa,
        s_matrix_file = matrix_filepath,
        scoring_function: str = 'js_divergence', 
        window_size: int = 3,
        use_gap_penalty: int = 1,
        win_lam: float = .5 ,
        gap_cutoff: float = .3 ,
        bg_distribution: list = blosum_background_distr[:]):
        
        self.msa = msa
        self.s_matrix_file = s_matrix_file
        self.scoring_function = scoring_function
        self.window_size = window_size
        self.use_gap_penalty = use_gap_penalty
        self.win_lam = win_lam
        self.gap_cutoff = gap_cutoff
        self.bg_distribution = bg_distribution
        
          
    def calculate_score(self):
            s_matrix_file = self.s_matrix_file
            msa = self.msa
            gap_cutoff = self.gap_cutoff
            bg_distribution = self.bg_distribution
            use_gap_penalty = self.use_gap_penalty
            names = msa.descriptions
            alignment = modify_alingment(msa.sequences)
            scoring_function = self.scoring_function
            window_size = self.window_size
            win_lam = self.win_lam
            
            if scoring_function == 'js_divergence':
                scoring = js_divergence
            elif scoring_function == 'relative_entropy':
                scoring = relative_entropy
            else:
                raise ValueError('Wrong scoring function')
            
            scores = []
            
            seq_weights = calculate_sequence_weights(alignment)
            s_matrix = read_scoring_matrix(s_matrix_file)
            
            for i in range(len(alignment[0])):
                col = get_column(i, alignment)
                if len(col) == len(alignment):
                    scores.append(scoring(col, s_matrix, bg_distribution, seq_weights, use_gap_penalty))
                else:
                    raise ValueError('alginment does not have proper length')
            
            if window_size > 0:
                scores = window_score(scores, window_size, win_lam)
            
            return scores   
        

def modify_alingment(alignment):
    new_alignment = []
    for seq in alignment:
        for i, aa in enumerate(seq):
            if aa not in iupac_alphabet:
                seq = seq.replace(aa, '-')
            elif aa in ['B' ,'Z','X', '*']:
                seq = seq.replace('B', 'D').replace('Z', 'Q').replace('X', '-').replace('*', '-')
        new_alignment.append(seq)
    return new_alignment
        

def read_scoring_matrix(sm_file):
    """ Read in a scoring matrix from a file, e.g., blosum80.bla, and return it
    as an array. """
    aa_index = 0
    first_line = 1
    row = []
    list_sm = [] # hold the matrix in list form

    try:
        matrix_file = open(sm_file, 'r')

        for line in matrix_file:

            if line[0] != '#' and first_line:
                first_line = 0
                if len(amino_acids) == 0:
                    for c in line.split():
                        aa_to_index[c.lower] = aa_index
                        amino_acids.append((c.lower))
                        aa_index += 1

            elif line[0] != '#' and first_line == 0:
                if len(line) > 1:
                    row = line.split()
                    list_sm.append(row)
      
    except:
        
        return np.identity(20)
        
    # if matrix is stored in lower tri form, copy to upper
    if len(list_sm[0]) < 20:
        for i in range(0,19):
            for j in range(i+1, 20):
                list_sm[i].append(list_sm[j][i])

    for i in range(len(list_sm)):
        for j in range(len(list_sm[i])):
            list_sm[i][j] = float(list_sm[i][j])

    return list_sm

def calculate_sequence_weights(msa):
    
    seq_weights = [0.] * len(msa)
    for i in range(len(msa[0])):
            freq_counts = [0] * len(amino_acids)
            col = []
            for j in range(len(msa)):
                if msa[j][i] != '-' : # ignore gaps
                    freq_counts[aa_to_index[msa[j][i]]] += 1
            num_observed_types = 0
            for j in range(len(freq_counts)):
                if freq_counts[j] > 0: num_observed_types +=1

            for j in range(len(msa)):
                d = freq_counts[aa_to_index[msa[j][i]]] * num_observed_types
                if d > 0:
                    seq_weights[j] += 1. / d

    for w in range(len(seq_weights)):
                seq_weights[w] /= len(msa[0])

    return seq_weights

def get_column(col_num, alignment):
    """Return the col_num column of alignment as a list."""
    col = []
    for seq in alignment:
        if col_num < len(seq):
            col.append(seq[col_num])
    
    return col
def gap_percentage(col):
    
    num_gaps = 0.

    for aa in col:
        if aa == '-': num_gaps += 1

    return num_gaps / len(col)

def weighted_freq_count_pseudocount(col, seq_weights, pc_amount):
    # if the weights do not match, use equal weight
    if len(seq_weights) != len(col):
        seq_weights = [1.] * len(col)
     
    aa_num = 0
    freq_counts = len(amino_acids)*[pc_amount] # in order defined by amino_acids

    for aa in amino_acids:
        for j in range(len(col)):
            if col[j] == aa:
                freq_counts[aa_num] += 1 * seq_weights[j]

        aa_num += 1

    for j in range(len(freq_counts)):
        freq_counts[j] = freq_counts[j] / (sum(seq_weights) + len(amino_acids) * pc_amount)

    return freq_counts

def weighted_gap_penalty(col, seq_weights):
    """ Calculate the simple gap penalty multiplier for the column. If the 
    sequences are weighted, the gaps, when penalized, are weighted 
    accordingly. """

    # if the weights do not match, use equal weight
    if len(seq_weights) != len(col):
        seq_weights = [1.] * len(col)
    
    gap_sum = 0.
    for i in range(len(col)):
        if col[i] == '-':
            gap_sum += seq_weights[i]

    return 1 - (gap_sum / sum(seq_weights))

     
def window_score(scores, window_len, lam=.5):
    
    w_scores = scores[:]
    for i in range(window_len, len(scores) - window_len):
        if scores[i] < 0: 
            continue
        sum = 0.
        num_terms = 0.
        for j in range(i - window_len, i + window_len + 1):
            if i != j and scores[j] >= 0:
                num_terms += 1
                sum += scores[j]
            if num_terms > 0:
                w_scores[i] = (1 - lam) * (sum / num_terms) + lam * scores[i]
    
    return w_scores

def relative_entropy(col, sim_matix, bg_distr, seq_weights, gap_penalty=1):
        distr = bg_distr[:]

        fc = weighted_freq_count_pseudocount(col, seq_weights, PSEUDOCOUNT)
    
        # remove gap count
        if len(distr) == 20:
            new_fc = fc[:-1]
            s = sum(new_fc)
            for i in range(len(new_fc)):
                new_fc[i] = new_fc[i] / s
            fc = new_fc
        
        if len(fc) != len(distr):
            return -1
        
        d = 0.
        for i in range(len(fc)):
            if distr[i] != 0.0:
                d += fc[i] * math.log(fc[i]/distr[i])
        
        d /= math.log(len(fc))

        if gap_penalty == 1: 
            return d * weighted_gap_penalty(col, seq_weights)
        else: 
            return d
    
def js_divergence(col, sim_matix, bg_distr, seq_weights, gap_penalty=1):
    
    distr = bg_distr[:]

    fc = weighted_freq_count_pseudocount(col, seq_weights, PSEUDOCOUNT)
    
    if len(distr) == 20:
        new_fc = fc[:-1]
        s = sum(new_fc)
        for i in range(len(new_fc)):
            new_fc[i] = new_fc[i] / s
        fc = new_fc
    if len(fc) != len(distr):
        d = -1
        return d

   
    r = [.5 * fc[i] + .5 * distr[i]for i in range(len(fc))]
    
    d = 0.
    for i in range(len(fc)):
        if r[i] != 0.0:
            if fc[i] == 0.0:
                d += distr[i] * math.log(distr[i]/r[i], 2)
            elif distr[i] == 0.0:
                d += fc[i] * math.log(fc[i]/r[i], 2) 
            else:
                d += fc[i] * math.log(fc[i]/r[i], 2) + distr[i] * math.log(distr[i]/r[i], 2)
                
                
    d /= 2        
            
        
    if gap_penalty == 1: 
        return d * weighted_gap_penalty(col, seq_weights)
    else: 
        return d