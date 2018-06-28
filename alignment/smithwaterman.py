import numpy as np
import pandas as pd

def align(seq1, seq2, open_penalty, extend_penalty, substitution_matrix):
    '''
    Given two sequences to algin substitution matrix, gap open penalty and
    extend penalty, output the alignment score matrix.

    '''

    # Initialization of score matrix and state state matrix by making an empty matrix.
    # Score Matrix: Matrix with alignment score calulated in each cell
    # State Matrix: Matrix with info on where maximum score is calculated

    mat_score = [[None]*len(seq2) for _ in range(0,len(seq1))]
    mat_state = [[None]*(len(seq2)+1) for _ in range(0,len(seq1)+1)]

    # Set the first row and column of score matrix to zero
#     mat_score[0,:] = 0
#     mat_score[:,0] = 0
    mat_score.insert(0,[0] * len(seq2))
    # Set the first column equal to zero
    for row in mat_score:
        row.insert(0,0)


#Start scoring matrix by first initializing max score as zero and
#replacing the max score by calculating deletion, insertion or match/mismatch
#score by getting scores from left, right and diagonal respectively.
#And also, update the state as we calculate the score.
    for i in range(1,len(seq1)+1):
        for j in range(1, len(seq2)+1):
            max_score = 0
            state = ""

# Calculate left well deletion score by determining whether it is
#gap open or extension and update the new score if the score calculated
#is greater than zero

            if mat_state[i][j-1] == "del":
                score = mat_score[i][j-1] + extend_penalty
            else:
                score = mat_score[i][j-1] + open_penalty

            if score > max_score:
                max_score = score
                state = "del"

# Calculate up well insertion score by determining whether it is
#gap open or extension and update the new score if the score calculated
#is greater than zero
            if mat_state[i-1][j] == "ins":
                score = mat_score[i-1][j] + extend_penalty
            else:
                score = mat_score[i-1][j] + open_penalty

            if score > max_score:
                max_score = score
                state = "ins"

# Calculate diagonal well alignment score by reffering to substitution
#matrix
            aa_i = seq1[i-1]
            aa_j = seq2[j-1]
            align_score = substitution_matrix.loc[str(aa_i),str(aa_j)]
            score = mat_score[i-1][j-1] + align_score
            if score > max_score:
                max_score = score
                state = "align"

            mat_score[i][j] = max_score
            mat_state[i][j] = state

    return np.asarray(mat_score), np.asanyarray(mat_state), max_score, state
