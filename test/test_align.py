#Write assert statments
from alignment import io


def test_alignment():
    seq1 = 'MAAAAAAGAGPEMVRGQVFDVGPRYTNLSYIGEGAYGMVCSAYDNVNKVR'
    seq2 = 'MAAAAAAGAGPEMVRGQVFDVGPRYTNLSYIGEGAYGMVCSAYDNVNKVR'
    blosum62 = io.read_score_matrix('/Users/Snow/Desktop/BMI203_Algorithm/HW3_SNaing_2018/ScoreMatrices/BLOSUM62')


    # update this assertion
    assert io.align(seq1, seq2, 5, 1, blosum62) == 260.0
