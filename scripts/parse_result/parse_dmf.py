import numpy as np
import pandas as pd
import os
import re
from sklearn.metrics import f1_score, matthews_corrcoef, confusion_matrix, average_precision_score, balanced_accuracy_score
from scipy.special import softmax

RESULT_DIR = 'DeepMicrobeFinder_results'
# RESULT_DIR = 'results_mu0_delta0.005'
# RESULT_DIR = 'results_2000'

results_fn = os.listdir(RESULT_DIR)
results_fn.sort()

def construct_target(nums):
    order = [[3], [4], [2], [0], [1]]
    return np.array(sum([order[i]*nums[i] for i in range(5)], []))

def construct_euk_array(array):
    euk_array = np.zeros((array.shape[0], 2))
    euk_array[np.logical_or.reduce((array==1, array==2, array==3, array==4)), 0] = 1
    euk_array[array==0, 1] = 1
    return euk_array

def construct_vir_array(array):
    euk_array = np.zeros((array.shape[0], 2))
    euk_array[np.logical_or.reduce((array==0, array==2, array==3)), 0] = 1
    euk_array[np.logical_or(array==1, array==4), 1] = 1

    return euk_array

def construct_array(array, target_idx):
    target_array = np.zeros((array.shape[0], 2))
    target_array[array==target_idx, 1] = 1
    target_array[array!=target_idx, 0] = 1
    return target_array

f1_scores = []

summary_df = pd.DataFrame(columns=['filename', 'category', 'f1_score', 'accuracy'])
categories = ['eukaryote', 'eukaryote_virus', 'plasmid', 'prokaryote', 'prokaryote_virus', 'virus', 'multiclass']

accs = []
f1s = []
for f in results_fn:
    # print(f[:-27])
    # if not f.endswith('_single_2000.txt'):
    if not f.endswith('_hybrid.txt'):
        continue
    result = pd.read_table(os.path.join(RESULT_DIR, f))
    result = result.iloc[:, 1:].to_numpy()

    idx = result.sum(axis=1)!=0
    result = result.argmax(axis=1)

    nums_str = re.findall(r'\d+', f)[1:]
    nums = [int(n) for n in nums_str]
    target = construct_target(nums)

    # result[result==3] = 2
    # result[result==4] = 3

    # result = construct_euk_array(result)
    # result = construct_vir_array(result)
    for i in range(5):
        result_class = construct_array(result, i)
        result_class = result_class.argmax(axis=1)
        target_class = construct_array(target, i)
        target_class = target_class.argmax(axis=1)
        acc = (result_class==target_class).sum()/len(result_class)
        f1 = f1_score(target_class, result_class)
        summary_df = pd.concat([summary_df, pd.DataFrame([['_'.join(nums_str), categories[i], f1, acc]], columns=['filename', 'category', 'f1_score', 'accuracy'])])

    summary_df = pd.concat([summary_df, pd.DataFrame([['_'.join(nums_str), 'multiclass', f1_score(target, result, average='weighted'), (target==result).sum()/len(result)]], columns=['filename', 'category', 'f1_score', 'accuracy'])])

    result_prok = result[np.logical_or.reduce((target==2, target==3, target==4))]
    target_prok = target[np.logical_or.reduce((target==2, target==3, target==4))]
    summary_df = pd.concat([summary_df, pd.DataFrame([['_'.join(nums_str), 'prok_groups', f1_score(target_prok, result_prok, average='weighted'), (target_prok==result_prok).sum()/len(result_prok)]], columns=['filename', 'category', 'f1_score', 'accuracy'])])

    result = construct_vir_array(result)
    result = result.argmax(axis=1)

    target = construct_vir_array(target)
    target = target.argmax(axis=1)

    summary_df = pd.concat([summary_df, pd.DataFrame([['_'.join(nums_str), 'virus', f1_score(target, result), (result==target).sum()/len(result)]], columns=['filename', 'category', 'f1_score', 'accuracy'])])


    # print(f1[0], f1[1], f1[2], f1[3], f1[4])
    # print(f'{f1}\t{1-idx.sum()/len(idx)}')
    # print(f'{1-idx.sum()/len(idx)}')
    # print(f1_score(target, result), balanced_accuracy_score(target, result), average_precision_score(target, result), matthews_corrcoef(target, result), acc)
    # print(confusion_matrix(target, result, normalize='true'))
    # print('')

# df = pd.DataFrame(f1_scores, columns=['Euk', 'EukVir', 'Plasmid', 'Prok', 'ProkVir'])
# df.to_csv(f'results/{RESULT_DIR}.csv', index=False)
# df.to_csv(f'results/results_mu0_delta0.csv', index=False)

# print(', '.join([str(i) for i in accs]))
# print(', '.join([str(i) for i in f1s]))

summary_df.to_csv(f'perf_summary/dmf.csv', index=False)
