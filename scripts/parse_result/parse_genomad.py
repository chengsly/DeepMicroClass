import numpy as np
import pandas as pd
import os
import re
from sklearn.metrics import f1_score

RESULT_DIR = 'data/result_others_10000/result_genomad'

results_fn = os.listdir(RESULT_DIR)
results_fn = [f for f in results_fn if f.endswith('aggregated_classification')]
results_fn.sort()
results_fn = [f'{fn}/{fn}.tsv' for fn in results_fn]


TARGET_DIR = 'data/sampled_sequence_target'

target_fn = os.listdir(TARGET_DIR)
target_fn.sort()


def construct_result(df):
    result = []
    for i, row in df.iterrows():
        if row['label'].startswith('plasmid'):
            result.append(1)
        else:
            result.append(0)
    return np.array(result)

misclassified = pd.DataFrame(columns=['Prok->Plas', 'ProkVir->Plas', 'Euk->Plas', 'EukVir->Plas', 'Plas->NonPlas'])
summary_df = pd.DataFrame(columns=['filename', 'f1_score', 'accuracy'])

mistakes = []

plasmid_summary_df = pd.DataFrame(columns=['filename', 'f1_score', 'accuracy'])
prokvirus_summary_df = pd.DataFrame(columns=['filename', 'f1_score', 'accuracy'])
multiclass_summary_df = pd.DataFrame(columns=['filename', 'f1_score', 'accuracy'])

accs = []
f1s = []
for f in results_fn:
    result = pd.read_table(os.path.join(RESULT_DIR, f), sep='\t', index_col=0)
    
    index = list(result.index)

    result = result.to_numpy().argmax(axis=1)
    
    target_filename = [fn for fn in target_fn if fn.startswith(f[:2])][0]
    target_pkl = pd.read_pickle(os.path.join(TARGET_DIR, f'{target_filename}'))
    target = np.array([target_pkl.loc[i] for i in index]).flatten()

    nums = re.findall(r'\d+', f)[1:6]

    target_binary = np.zeros(len(target))
    target_binary[target==2] = 1
    result_binary = np.zeros(len(result))
    result_binary[result==1] = 1
    plasmid_summary_df = pd.concat([plasmid_summary_df, pd.DataFrame([['_'.join(nums), f1_score(target_binary, result_binary), (target_binary==result_binary).sum()/len(target_binary)]], columns=['filename', 'f1_score', 'accuracy'])])

    target_binary = np.zeros(len(target))
    # target_binary[np.logical_or(target==4, target==2)] = 1
    target_binary[target==4] = 1
    result_binary = np.zeros(len(result))
    result_binary[result==2] = 1
    prokvirus_summary_df = pd.concat([prokvirus_summary_df, pd.DataFrame([['_'.join(nums), f1_score(target_binary, result_binary), (target_binary==result_binary).sum()/len(target_binary)]], columns=['filename', 'f1_score', 'accuracy'])])

    # genomad_target = np.zeros(len(target)) + 3
    # genomad_target[target==2] = 1
    # genomad_target[target==4] = 2
    # genomad_target[target==3] = 0



    # summary_df = pd.concat([summary_df, pd.DataFrame([['_'.join(nums), f1, acc]], columns=['filename', 'f1_score', 'accuracy'])])

    mistake = [(target[result==1]==3).sum(), (target[result==1]==4).sum(), (target[result==1]==0).sum(), (target[result==1]==1).sum(), (target[result!=1]==2).sum()] # For plasmid
    # mistake = [(target[result==2]==3).sum(), (target[result==2]==0).sum(), (target[result==2]==1).sum(), (target[result==2]==2).sum(), (target[result!=2]==4).sum()] # For prokaryotic virus
    mistakes.append(mistake)

    # print(f)
    # print(f'Acc: {acc}\tF1: {f1}\tUnknown: {(result==-1).sum()/len(result)}')
    # print(acc, end=', ')
    # print(f1, end=', ')
    # accs.append(acc)
    # f1s.append(f1)
# print(', '.join([str(i) for i in accs]))
# print(', '.join([str(i) for i in f1s]))
# summary_df.to_csv('perf_summary/plasflow.csv', index=False)

# plasmid_summary_df.to_csv('perf_summary/genomad_plasmid.csv', index=False)
# prokvirus_summary_df.to_csv('perf_summary/genomad_vir.csv', index=False)

misclassified = pd.DataFrame(mistakes, columns=['Prok->Plas', 'ProkVir->Plas', 'Euk->Plas', 'EukVir->Plas', 'Plas->NonPlas']) # For plasmid
misclassified.to_csv('perf_summary/misclassified_genomad_plasmid.csv', index=False)
# misclassified = pd.DataFrame(mistakes, columns=['Prok->ProkVir', 'Euk->ProkVir', 'EukVir->ProkVir', 'Plas->ProkVir', 'ProkVir->NonProkVir']) # For prokaryotic virus
# misclassified.to_csv('perf_summary/misclassified_genomad_prokvir.csv', index=False)