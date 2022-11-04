import numpy as np
import pandas as pd
import os
import re
from sklearn.metrics import f1_score, matthews_corrcoef, confusion_matrix, average_precision_score, balanced_accuracy_score
from scipy.special import softmax

# RESULT_DIR = 'DeepMicrobeFinder_results'
# RESULT_DIR = 'results_mu0_delta0.005'
RESULT_DIR = 'results_2000'

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

    # euk_array[np.logical_or.reduce((array==0, array==2, array==3, array==1)), 0] = 1
    # euk_array[array==4, 1] = 1

    # euk_array[np.logical_or.reduce((array==0, array==2, array==4, array==1)), 0] = 1
    # euk_array[array==3, 1] = 1

    # euk_array[np.logical_or.reduce((array==0, array==1, array==3, array==4)), 0] = 1
    # euk_array[array==2, 1] = 1

    return euk_array

def construct_array(array, target_idx):
    target_array = np.zeros((array.shape[0], 2))
    target_array[array==target_idx, 1] = 1
    target_array[array!=target_idx, 0] = 1
    return target_array

f1_scores = []

accs = []
f1s = []
for f in results_fn:
    # print(f[:-27])
    # if not f.endswith('_single_2000.txt'):
    if not f.endswith('_hybrid.txt'):
        continue
    result = pd.read_table(os.path.join(RESULT_DIR, f))
    result = result.iloc[:, 1:].to_numpy()
    # result[:, :2] = 0
    idx = result.sum(axis=1)!=0
    #result = softmax(result, axis=1)
    # idx = np.logical_and(result.max(axis=1) > t, idx)
    result = result.argmax(axis=1)

    # result[result==3] = 2
    # result[result==4] = 3

    # result = construct_euk_array(result)
    # result = construct_vir_array(result)
    result = construct_array(result, 4)
    result = result.argmax(axis=1)
    #result = result[idx]

    nums = re.findall(r'\d+', f)[1:]
    nums = [int(n) for n in nums]
    target = construct_target(nums)

    # target[target==3] = 2
    # target[target==4] = 3

    # selection = np.logical_or(target==2, target==3, target==4)

    # result = result[selection]
    # target = target[selection]

    # target = construct_euk_array(target)
    # target = construct_vir_array(target)
    target = construct_array(target, 4)
    target = target.argmax(axis=1)
    # target = target[idx]
    try:
        acc = (result==target).sum() / len(target)
        # f1 = f1_score(target, result, average=None)
        f1 = f1_score(target, result)
        # f1 = matthews_corrcoef(target, result)

    except:
        print(f)
#     scores.append(f1)
    # print(f)
    # print(f'Acc: {acc}\tF1: {f1}')
    # print(acc, end=', ')
    accs.append(acc)
    f1s.append(f1)

    # print(acc)
    # print(f1)
    # f1_scores.append(f1)

    # print(f1[0], f1[1], f1[2], f1[3], f1[4])
    # print(f'{f1}\t{1-idx.sum()/len(idx)}')
    # print(f'{1-idx.sum()/len(idx)}')
    # print(f1_score(target, result), balanced_accuracy_score(target, result), average_precision_score(target, result), matthews_corrcoef(target, result), acc)
    # print(confusion_matrix(target, result, normalize='true'))
    # print('')

# df = pd.DataFrame(f1_scores, columns=['Euk', 'EukVir', 'Plasmid', 'Prok', 'ProkVir'])
# df.to_csv(f'results/{RESULT_DIR}.csv', index=False)
# df.to_csv(f'results/results_mu0_delta0.csv', index=False)

print(', '.join([str(i) for i in accs]))
print(', '.join([str(i) for i in f1s]))
