from Bio import SeqIO
import glob
import os
import numpy as np

sequence_dir = "/Users/Snow/Desktop/BMI203_Algorithm/HW3_due_02_23/sequences"

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
    print("Read in %d fasta files"%len(files))
    return sequences

def get_negpair_seq(dir):
    '''
    Get sequences for negative pairs.
    Input: Directory for all the sequences
    Output: List of negative pairs sequences
    '''
    files = glob.glob(dir + '/*.fa')
    negpair_list = open('/Users/Snow/Desktop/BMI203_Algorithm/HW3_due_02_23/Negpairs.txt').read().splitlines()
    neg_files = []
    negpair_sequences = []
    for seq_file in files:
        for negpair in negpair_list:
            if negpair[10:] in seq_file[-12:]:
                neg_files.append(seq_file)
    for neg_file in neg_files:
        with open(neg_file) as neg:
            for name, neg_seq in read_fasta(neg):
                negpair_sequences.append(neg_seq)
    print("Read in %d negative pair fasta files" %len(negpair_list))
    return negpair_sequences

def get_pospair_seq(dir):
    '''
    Get sequences for positive pairs.
    Input: Directory for all the sequences
    Output: List of positive pairs sequences
    '''
    files = glob.glob(dir + '/*.fa')
    pospair_list = open('/Users/Snow/Desktop/BMI203_Algorithm/HW3_due_02_23/Pospairs.txt').read().splitlines()
    pos_files = []
    pospair_sequences = []
    for seq_file in files:
        for pospair in pospair_list:
            if pospair[10:] in seq_file[-12:]:
                pos_files.append(seq_file)
    for pos_file in pos_files:
        with open(pos_file) as pos:
            for name, pos_seq in read_fasta(pos):
                pospair_sequences.append(pos_seq)
    print("Read in %d positive pair fasta files" %len(pospair_list))
    return pospair_sequences

get_sequences(sequence_dir)
get_negpair_seq(sequence_dir)
get_pospair_seq(sequence_dir)
