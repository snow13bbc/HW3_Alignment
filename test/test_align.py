#Write assert statments
from alignment import io
import numpy as np
import pandas as pd

def test_alignment():
    seq1 = 'MAAAAAAGAGPEMVRGQVFDVGPRYTNLSYIGEGAYGMVCSAYDNVNKVR'
    seq2 = 'MAAAAAAGAGPEMVRGQVFDVGPRYTNLSYIGEGAYGMVCSAYDNVNKVR'
    blosum62 = io.read_score_matrix('matrices/BLOSUM62')
    aln, score, _ = io.align(seq1, seq2, -5, -1, blosum62)

    # update this assertion
    assert score == 260.0
