def obj_fun(pos, neg):
    '''
    return sum of TP rates for FP rates of 0.0, 0.1, 0.2, and 0.3.
    '''

    fp = []
    tp = []

    for t in [0.0,0.1,0.2,0.3]:
        thresh = np.percentile(neg, 100*(1-t))
        false_positives = [i for i in neg if i > thresh]
        fpr = len(false_positives)/float(len(neg))
        true_positives = [i for i in pos if i > thresh]
        tpr = len(true_positives)/float(len(pos))
        tp.append(tpr)
        fp.append(fpr)

    return(sum(tp))

def new_matrix(matrices):
    '''
    Randomoly switches indices of original scoring matrix. (i.e. giving amino acids new scores.)
    Return 5 new matrices.
    '''
    new_matrices = []
    for m in matrices:
        new_mat = m.copy()
        a,b,c,d = np.random.choice(range(len(m.columns)), 4)
        index = list(m.columns)
        index[a], index[b], index[c], index[d]= index[d], index[a], index[b], index[c]
        new_mat.columns = index
        new_mat.index = index
        new_matrices.append(new_mat)
    return new_matrices 
