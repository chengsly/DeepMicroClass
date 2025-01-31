import numpy as np
import pandas as pd
import os
import re
from sklearn.metrics import f1_score

RESULT_DIR = 'data/result_others/result_whokaryote'

results_fn = os.listdir(RESULT_DIR)
results_fn = [f'{f}/whokaryote_predictions_T.tsv' for f in results_fn if not f.startswith('log')]
results_fn.sort()

TARGET_DIR = 'data/sampled_sequence_target'

target_fn = os.listdir(TARGET_DIR)
target_fn.sort()

# a = pd.read_table(os.path.join(RESULT_DIR, results_fn[0]), index_col=0)

def construct_result(df):
    result = []
    prok = set(('archaea', 'bacteria', 'prokaryote'))
    euk = set(('eukaryote', 'organelle'))
    for i, row in df.iterrows():
        if row[0] in euk:
            result.append(1)
        elif row[0] in prok:
            result.append(0)
        else:
            result.append(0)
    return np.array(result)

mistakes = []

summary_df = pd.DataFrame(columns=['filename', 'f1_score', 'accuracy'])

for f in results_fn:
    nums = re.findall(r'\d+', f)[1:]
    result = pd.read_table(os.path.join(RESULT_DIR, f), index_col=0)
    
    index = list(result.index)
    index = [i.split()[0] for i in index]
    result.index = index
    
    target_pkl = pd.read_pickle(os.path.join(TARGET_DIR, f'{f[:-29]}.fa.pkl'))
    target = np.array([target_pkl.loc[i] for i in index]).flatten()

    missed = len(target_pkl) - len(target)
    # print(missed)
    missed_label = []
    for i in target_pkl.index:
        if i not in result.index:
            missed_label.append(target_pkl.loc[i, 0])
    missed_label = np.array(missed_label).flatten()
    target = np.concatenate((target, missed_label))

    target_binary = np.zeros(len(target))
    target_binary[target==0] = 1

    result = construct_result(result)

    result = np.concatenate((result, np.zeros(missed)))

    try:
        acc = (result==target_binary).sum() / len(target_binary)
        # f1 = f1_score(target_binary[result!=-1], result[result!=-1])
        f1 = f1_score(target_binary, result)
    except:
        print(f)

    mistake = [(target[result==1]==3).sum(), (target[result==1]==4).sum(), (target[result==0]==0).sum(), (target[result==1]==1).sum(), (target[result==1]==2).sum()]
    mistakes.append(mistake)

    # print(acc, end=', ')
    # print(f1, end=', ')
    summary_df = pd.concat([summary_df, pd.DataFrame([['_'.join(nums), f1, acc]], columns=['filename', 'f1_score', 'accuracy'])], axis=0)
summary_df.to_csv('perf_summary/whokaryote.csv', index=False)

misclassified = pd.DataFrame(mistakes, columns=['Prok->Euk', 'ProkVir->Euk', 'Euk->Non-Euk', 'EukVir->Euk', 'Plas->Euk'])
misclassified.to_csv('perf_summary/misclassified_whokaryote.csv', index=False)