from Bio import SeqIO
import glob
import os
import numpy as np

def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))

def get_sequences(dir):
    files = glob.glob(dir + '/*.fa')
    sequences = []
    for filepath in glob.iglob(os.path.join(dir, "*.fa")):
        with open(filepath) as fp:
            for name, seq in read_fasta(fp):
                sequences.append(seq)
    # print("Read in %d fasta files"%len(files))
    return sequences

def get_negpair_seq(dir):
    '''
    Get sequences for negative pairs.
    Input: Directory for all the sequences
    Output: List of negative pairs sequences
    '''
    files = glob.glob(dir + '/*.fa')
    negpair_list = open('/Users/Snow/Desktop/BMI203_Algorithm/HW3_SNaing_2018/Negpairs.txt').read().splitlines()
    neg_files, neg_files2 = [],[]
    negpair_sequences, negpair2_sequences  = [], []
    for seq_file in files:
        for negpair in negpair_list:
            if negpair[10:22] in seq_file[-12:]:
                neg_files.append(seq_file)
            if negpair[33:] in seq_file[-12:]:
                neg_files2.append(seq_file)
    for neg_file in neg_files:
        with open(neg_file) as neg:
            for name, neg_seq in read_fasta(neg):
                negpair_sequences.append(neg_seq)
    for neg_file in neg_files2:
        with open(neg_file) as neg:
            for name, neg2_seq in read_fasta(neg):
                negpair2_sequences.append(neg2_seq)
    #Some of the sequences has "x" as a gap. I'm replacing those 'x's with '-'
    #so it works with skbio package. Not really a good way but it works.
    without_x = [s for s in negpair2_sequences if 'x' not in s]
    with_x= [s for s in negpair2_sequences if 'x' in s]
    with_x =' '.join(with_x)
    with_x = with_x.replace('x', '-')
    without_x.append(with_x)
    negpair2_sequences = without_x
    # print("Read in %d negative pair fasta files" %len(negpair_sequences))
    # print("Read in %d negative pair fasta files" %len(negpair2_sequences))
    return negpair_sequences, negpair2_sequences

def get_pospair_seq(dir):
    '''
    Get sequences for positive pairs.
    Input: Directory for all the sequences
    Output: List of positive pairs sequences
    '''
    files = glob.glob(dir + '/*.fa')
    pospair_list = open('/Users/Snow/Desktop/BMI203_Algorithm/HW3_SNaing_2018/Pospairs.txt').read().splitlines()
    pos_files,pos2_files = [],[]
    pospair_sequences,pospair2_sequences = [],[]
    for seq_file in files:
        for pospair in pospair_list:
            if pospair[10:22] in seq_file[-12:]:
                pos_files.append(seq_file)
            if pospair[33:] in seq_file[-12:]:
                pos2_files.append(seq_file)
    for pos_file in pos_files:
        with open(pos_file) as pos:
            for name, pos_seq in read_fasta(pos):
                pospair_sequences.append(pos_seq)
    for pos_file in pos2_files:
        with open(pos_file) as pos:
            for name, pos2_seq in read_fasta(pos):
                pospair2_sequences.append(pos2_seq)
    # print("Read in %d positive pair fasta files" %len(pospair_sequences))
    # print("Read in %d positive pair fasta files" %len(pospair2_sequences))
    return pospair_sequences, pospair2_sequences

def read_score_matrix(filename):
    '''
    Read score matrices as dataframe.
    Input: score matrix
    Output: dataframe
    '''
    with open(filename) as mat:
        matrix = mat.read()
        lines = matrix.strip().split('\n')
        score_matrix = []
        for line in lines:
            if not line.startswith("#"):
                score_matrix.append(line)
    aas_column = pd.DataFrame(score_matrix[0].strip().split("  "))
    scores = []
    for i in range(len(score_matrix)-1):
        ls = score_matrix[i+1].strip()
        new = ls.split(" ")
        sub = []
        for n in new:
            if n != "":
                sub.append(int(n))
        scores.append(sub)
        score_df = pd.DataFrame(scores)
        score_df.columns = score_matrix[0].strip().split("  ")
    df_score = aas_column.join(score_df)
    df_score = df_score.set_index([0])
    return df_score

def align(seq1, seq2, gap_open_penalty, gap_extend_penalty, score_matrix):
    #using scikit learn package for local alignment 
    sequence1 = Protein(seq1)
    sequence2 = Protein(seq2)
    return local_pairwise_align(sequence1, sequence2, gap_open_penalty, gap_extend_penalty, score_matrix)
