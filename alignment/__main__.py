import sys
from .io import *
from .smithwaterman import *
from .optimization import *

#Load all the necessary files
directory = '/Users/Snow/Documents/BMI203_Algorithm/HW3_SNaing_2018/sequences/'
negpairlist_filename = '/Users/Snow/Documents/BMI203_Algorithm/HW3_SNaing_2018/Negpairs.txt'
pospairlist_filename = '/Users/Snow/Documents/BMI203_Algorithm/HW3_SNaing_2018/Pospairs.txt'
negpairs1, negpairs2 = get_negpair_seq(directory, negpairlist_filename)
pospairs1, pospairs2 = get_pospair_seq(directory, pospairlist_filename)
blosum50 = read_score_matrix('/Users/Snow/Documents/BMI203_Algorithm/HW3_SNaing_2018/matrices/BLOSUM50')
blosum62 = read_score_matrix('/Users/Snow/Documents/BMI203_Algorithm/HW3_SNaing_2018/matrices/BLOSUM62')
matio = read_score_matrix('/Users/Snow/Documents/BMI203_Algorithm/HW3_SNaing_2018/matrices/MATIO')
pam100 = read_score_matrix('/Users/Snow/Documents/BMI203_Algorithm/HW3_SNaing_2018/matrices/PAM100')
pam250 = read_score_matrix('/Users/Snow/Documents/BMI203_Algorithm/HW3_SNaing_2018/matrices/PAM250')

# Finding best gap panalties scores at TPR = 0.7 for BLOSUM50

possible_gap_extend_penalties = list(range(1,6))
possible_gap_open_penalties = list(range(1,21))
best_fpr = float('inf')
best_gap = None
best_extension = None

for open_penalty in possible_gap_open_penalties:
    for extend_penalty in possible_gap_extend_penalties:
        print ("At open penalty score " + str(open_penalty) + " and extend penalty "+ str(extend_penalty))
        pos_align_score = []
        neg_align_score = []
        for i in range(0,len(pospairs1)):
            mat_score_p,mat_state_p, score_p, state_p = align(pospairs1[i],pospairs2[i], -open_penalty,-extend_penalty, blosum50)
            print ("In " + str(i) + " Pospair: score_p is " + str(score_p))
            mat_score_n,mat_state_n, score_n, state_n = align(negpairs1[i],negpairs2[i], -open_penalty,-extend_penalty, blosum50)
            print ("In " + str(i) + " Negpair: score_n is " + str(score_n))
            pos_align_score.append(score_p)
            neg_align_score.append(score_n)
        pos_align_score.sort()
        neg_align_score.sort()
        tpr = 0.7 #tpr = #>thresh / total scores
        cutoff_index = int((1-tpr)*len(pos_align_score))# find the index that makes it so tpr*lenscores is > threshold
        threshold = pos_align_score[cutoff_index]
        neg_score_np = np.array(neg_align_score)
        fpr = len(neg_score_np[neg_score_np > threshold])/len(neg_score_np) #False Positive Rate
        if fpr < best_fpr:
            best_fpr = fpr
            best_gap = open_penalty
            best_extension = extend_penalty

            print ("Best gap now is " + str(best_gap))
            print ("Best extend now is "+ str(best_extension))

print('Best Open and extension gap penalty scores are %d and %d.' % (best_gap,best_extension))
print('Best False positive rate is %d.' %(best_fpr))

# Finding best scoring matrix; Plotting ROC graph

matrices = [blosum50,blosum62,matio,pam100,pam250]
matrix_names = ['blosum50','blosum62','matio','pam100','pam250']
score_list = [] #Combined negative and positive scores as numpy array
true_pos = [1]*50 # ture_pos = [1]*len(pos)
true_neg = [0]*50 # ture_neg = [0]*len(neg)
combined_trues = np.concatenate((true_pos,true_neg)) #Telling sckit-roc curve which are pos, which are neg

for m in matrices:
    pos_align_score =[]
    neg_align_score = []
    for i in range(0,len(pospairs1)):
        mat_score_p,mat_state_p, score_p, state_p = align(pospairs1[i],pospairs2[i], -best_gap ,-best_extension, m)
        mat_score_n,mat_state_n, score_n, state_n = align(negpairs1[i],negpairs2[i], -best_gap ,-best_extension, m)
        pos_align_score.append(score_p)
        neg_align_score.append(score_n)
    score_list.append(np.concatenate((np.asarray(pos_align_score), np.asarray(neg_align_score))))

plt.figure()
plt.plot([0,1],[0,1])
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("ROC")
plt.axhline(0.7, color = 'm', linestyle = '--')
for score, label in zip(score_list, matrix_names):
    fpr,tpr, treshold = roc_curve(combined_trues,score, pos_label=1)
    plt.plot(fpr, tpr, label=label)
plt.legend()
plt.show()
plt.close()

# Normalize SW score with shortest sequence length

norm_pos_score, norm_neg_score = [], []
for i in range(0,len(pospairs1)):
    mat_score_p,mat_state_p, score_p, state_p = align(pospairs1[i],pospairs2[i], -best_gap ,-best_extension, blosum62)
    mat_score_n,mat_state_n, score_n, state_n = align(negpairs1[i],negpairs2[i], -best_gap ,-best_extension, blosum62)
    norm_pos = score_p/min(len(pospairs1[i]),len(pospairs2[i]))
    norm_neg = score_n/min(len(negpairs1[i]),len(negpairs2[i]))
    norm_pos_score.append(norm_pos)
    norm_neg_score.append(norm_neg)

norm_pos_np = np.asarray(norm_pos_score)
norm_neg_np = np.asarray(norm_neg_score)
norm_combined = np.concatenate((norm_pos_np,norm_neg_np))
norm_true_pos = [1]*len(norm_pos_score)
norm_true_neg = [0]*len(norm_neg_score)
norm_combined_trues = np.concatenate((norm_true_pos,norm_true_neg))

norm_fpr,norm_tpr, treshold = roc_curve(norm_combined_trues,norm_combined, pos_label=1) #This is for normalized ROC
fpr,tpr, treshold = roc_curve(combined_trues,score_list[1], pos_label=1) # This is for regular ROC
plt.figure()
plt.plot(norm_fpr,norm_tpr, 'r', label = 'BLOSUM62-Normalized')
plt.plot(fpr, tpr, label = 'BLOSUM62')
plt.plot([0,1],[0,1])
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.legend()
plt.title("ROC: Best Matrix - Blosum 62")
plt.axhline(0.7, color = 'm', linestyle = '--')
plt.show()
plt.close()

# Optimization
new_matrix_names= ['blosum50-mut','blosum62-mut','matio-mut','pam100-mut','pam250-mut']
new_score_list = [] #Combined negative and positive scores as numpy array
new_true_pos= [1]*50 # ture_pos = [1]*len(pos)
new_true_neg= [0]*50 # ture_neg = [0]*len(neg)
new_combined_trues= np.concatenate((new_true_pos,new_true_neg)) #Telling sckit-roc curve which are pos, which are neg
obj_sum_tps= [] #For calculating sum of tps using objective function
max_fitness = 0 # To calculate max fitness
for m in new_matrix(matrices):
    new_pos_align_score =[]
    new_neg_align_score = []
    print ('At matrix' + str(m))
    for i in range(0,len(pospairs1)):
        mat_score_p,mat_state_p, new_score_p, state_p = align(pospairs1[i],pospairs2[i], -best_gap ,-best_extension, m)
        mat_score_n,mat_state_n, new_score_n, state_n = align(negpairs1[i],negpairs2[i], -best_gap ,-best_extension, m)

        new_pos_align_score.append(new_score_p)
        new_neg_align_score.append(new_score_n)
        print('At score!')
        sum_tps = obj_fun(new_pos_align_score, new_neg_align_score)
        obj_sum_tps.append(sum_tps)
        if sum_tps > max_fitness:
            max_fitness = sum_tps
            print ('calculated Max Fitness!')
    new_score_list.append(np.concatenate((np.asarray(new_pos_align_score), np.asarray(new_neg_align_score))))

plt.figure()
plt.plot([0,1],[0,1])
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("ROC using new mutated score matrices")
plt.axhline(0.7, color = 'y', linestyle = '--')
for new_score, new_label in zip(new_score_list, new_matrix_names):
    new_fpr,new_tpr, treshold = roc_curve(new_combined_trues,new_score, pos_label=1)
    plt.plot(new_fpr, new_tpr, label=new_label)
plt.legend()
plt.show()
plt.close()

#save new matrices
new_mat = new_matrix(matrices)
blosum50_mut = pd.DataFrame(new_mat[0])
blosum62_mut = pd.DataFrame(new_mat[1])
matio_mut = pd.DataFrame(new_mat[2])
pam100_mut = pd.DataFrame(new_mat[3])
pam250_mut = pd.DataFrame(new_mat[4])
blosum50_mut.to_csv('blosum50-mut.txt', index=True, sep=' ', header=True)
blosum62_mut.to_csv('blosum62-mut.txt', index=True, sep=' ', header=True)
matio_mut.to_csv('matio_mut.txt', index=True, sep=' ', header=True)
pam100_mut.to_csv('pam100-mut.txt', index=True, sep=' ', header=True)
pam250_mut.to_csv('pam250-mut.txt', index=True, sep=' ', header=True)
