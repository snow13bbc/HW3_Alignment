#Write assert statments
from alignment import io, optimization, smithwaterman
import numpy as np
import pandas as pd

def test_alignment():
    seq1 = 'MAAAAAAGAGPEMVRGQVFDVGPRYTNLSYIGEGAYGMVCSAYDNVNKVR'
    seq2 = 'MAAAAAAGAGPEMVRGQVFDVGPRYTNLSYIGEGAYGMVCSAYDNVNKVR'
    blosum62 = io.read_score_matrix('matrices/BLOSUM50')
    _, _1, score, _2 = smithwaterman.align(seq1, seq2, -5, -1, blosum62)

    # update this assertion
    assert score == 333.0



def test_optimization():
    df=pd.read_csv('matrices_mut/pam100-mut.txt', sep=' ', header=0)
    df = df.set_index('Unnamed: 0')
    seq1 = 'MAAAAAAGAGPEMVRGQVFDVGPRYTNLSYIGEGAYGMVCSAYDNVNKVR'
    seq2 = 'MAAAAAAGAGPEMVRGQVFDVGPRYTNLSYIGEGAYGMVCSAYDNVNKVR'
    _, _1, score, _2 = smithwaterman.align(seq1, seq2, -5, -1, df)

    # update this assertion
    assert score == 301.0
