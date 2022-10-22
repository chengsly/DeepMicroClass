import numpy as np
import pandas as pd
import os
import re
from sklearn.metrics import f1_score

RESULT_DIR = 'data/result_vibrant'

results_fn = os.listdir(RESULT_DIR)
# results_fn = [i for i in results_fn if i.startswith('VIBRANT_phages_')]
results_fn.sort()

TARGET_DIR = 'data/sampled_sequence_target'

target_fn = os.listdir(TARGET_DIR)
target_fn.sort()

def construct_result(df):
    result = [1 for _ in range(len(set(df.index)))]
    # for i, row in df.iterrows():
    #     result.append(1)
    return np.array(result)

summary_df = pd.DataFrame(columns=['filename', 'f1_score', 'accuracy'])

accs = []
f1s = []
for f in results_fn:
    # result = pd.read_table(os.path.join(RESULT_DIR, f, f'{f[15:]}.phages_combined.txt'), index_col=0, header=None)
    result = pd.read_table(os.path.join(RESULT_DIR, f, f'VIBRANT_phages_{f[8:]}', f'{f[8:]}.phages_combined.txt'), index_col=0, header=None)
    
    index = list(result.index)
    index = [i.split()[0] for i in index]
    result.index = index
    
    # target_pkl = pd.read_pickle(os.path.join(TARGET_DIR, f'{f[15:]}.fa.pkl'))
    target_pkl = pd.read_pickle(os.path.join(TARGET_DIR, f'{f[8:]}.fa.pkl'))
    # target = np.array([target_pkl.loc[i] for i in index]).flatten()
    target = np.array([target_pkl.loc[i] for i in set(index)]).flatten()

    missed = len(target_pkl) - len(target)
    missed_label = []
    for i in target_pkl.index:
        if i not in result.index:
            missed_label.append(target_pkl.loc[i, 0])
    missed_label = np.array(missed_label).flatten()
    target = np.concatenate((target, missed_label))

    target_binary = np.zeros(len(target))
    target_binary[target==4] = 1

    result = construct_result(result)

    result = np.concatenate((result, np.zeros(missed)))
    f1 = f1_score(target_binary, result)
    try:
        acc = (result==target_binary).sum() / len(target_binary)
        f1 = f1_score(target_binary[result!=-1], result[result!=-1])
    except:
        print(f)
    accs.append(acc)
    f1s.append(f1)

    nums = re.findall(r'\d+', f)[1:]
    summary_df = pd.concat([summary_df, pd.DataFrame([['_'.join(nums), f1, acc]], columns=['filename', 'f1_score', 'accuracy'])])

    # print(f)
    # print(f'Acc: {acc}\tF1: {f1}\tUnknown: {(result==-1).sum()/len(result)}')
    # print(acc, end=', ')
    # print(f1, end=', ')
# print(', '.join([str(i) for i in accs]))
# print(', '.join([str(i) for i in f1s]))
summary_df.to_csv('perf_summary/vibrant.csv', index=False)