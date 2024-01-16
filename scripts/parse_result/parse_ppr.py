import numpy as np
import pandas as pd
import os
import re
from sklearn.metrics import f1_score, confusion_matrix

RESULT_DIR = 'data/result_others/result_ppr'
# RESULT_DIR = 'result_other_2000/result_ppr'

results_fn = os.listdir(RESULT_DIR)
results_fn = [f for f in results_fn if not f.endswith('fasta')]
results_fn.sort()

TARGET_DIR = 'data/sampled_sequence_target'

target_fn = os.listdir(TARGET_DIR)
target_fn.sort()


def construct_result(df, target):
    result = []
    for i, row in df.iterrows():
        if row['Possible_source'] == target:
            result.append(1)
        else:
            result.append(0)
    return np.array(result)

result_mapping = {
    'chromosome': 0,
    'plasmid': 1,
    'phage': 2,
}

misclassified = pd.DataFrame(columns=['Prok->Plas', 'ProkVir->Plas', 'Euk->Plas', 'EukVir->Plas', 'Plas->NonPlas'])
mistakes = []

plasmid_summary_df = pd.DataFrame(columns=['filename', 'f1_score', 'accuracy'])
prokvirus_summary_df = pd.DataFrame(columns=['filename', 'f1_score', 'accuracy'])
multiclass_summary_df = pd.DataFrame(columns=['filename', 'f1_score', 'accuracy'])
accs = []
f1s = []
for f in results_fn:
    org_result = pd.read_table(os.path.join(RESULT_DIR, f), sep=',')
    
    index = list(org_result['Header'])
    index = [i.split()[0] for i in index]
    org_result.index = index
    
    target_pkl = pd.read_pickle(os.path.join(TARGET_DIR, f'{f[:-4]}.fa.pkl'))
    target = np.array([target_pkl.loc[i] for i in index]).flatten()

    missed = len(target_pkl) - len(target)
    missed_label = []
    for i in target_pkl.index:
        if i not in org_result.index:
            missed_label.append(target_pkl.loc[i, 0])
    missed_label = np.array(missed_label).flatten()
    target = np.concatenate((target, missed_label))


    target_binary = np.zeros(len(target))
    # target_binary[target==2] = 1
    target_binary[target==4] = 1

    # result = construct_result(result, 'plasmid')
    result = construct_result(org_result, 'phage')
    result = np.concatenate((result, np.zeros(missed)))

    nums = re.findall('\d+', f)[1:]
    prokvirus_summary_df = pd.concat([prokvirus_summary_df, pd.DataFrame([['_'.join(nums), f1_score(target_binary, result), (result==target_binary).sum()/len(result)]], columns=['filename', 'f1_score', 'accuracy'])])

    target_binary = np.zeros(len(target))
    target_binary[target==2] = 1
    result = construct_result(org_result, 'plasmid')
    result = np.concatenate((result, np.zeros(missed)))
    plasmid_summary_df = pd.concat([plasmid_summary_df, pd.DataFrame([['_'.join(nums), f1_score(target_binary, result), (result==target_binary).sum()/len(result)]], columns=['filename', 'f1_score', 'accuracy'])])

    prok_idx = np.logical_and(target!=0, target!=1)
    ppr_target = np.zeros(len(target))
    ppr_target[target==2] = 1
    ppr_target[target==4] = 2
    ppr_target = ppr_target[prok_idx]
    result = list(org_result['Possible_source'])
    result = [result_mapping[i] for i in result]
    result = np.array(result)
    result = np.concatenate((result, np.zeros(missed)))
    result = result[prok_idx]
    multiclass_summary_df = pd.concat([multiclass_summary_df, pd.DataFrame([['_'.join(nums), f1_score(ppr_target, result, average='weighted'), (result==ppr_target).sum()/len(result)]], columns=['filename', 'f1_score', 'accuracy'])])

    # try:
    #     acc = (result==target_binary).sum() / len(target_binary)
    #     f1 = f1_score(target_binary[result!=-1], result[result!=-1])
    # except:
    #     print(f)



    # accs.append(acc)
    # f1s.append(f1)
    
    # mistake = [(target[result==1]==3).sum(), (target[result==1]==4).sum(), (target[result==1]==0).sum(), (target[result==1]==1).sum(), (target_binary[result==0]==1).sum()] # For plasmid
    result = construct_result(org_result, 'phage')
    result = np.concatenate((result, np.zeros(missed)))
    mistake = [(target[result==1]==3).sum(), (target[result==1]==0).sum(), (target[result==1]==1).sum(), (target[result==1]==2).sum(), (target[result==0]==4).sum()] # For prokaryotic virus
    mistakes.append(mistake)

    # misclassified = pd.concat([misclassified, pd.DataFrame([mistake], columns=['Prok->Plas', 'ProkVir->Plas', 'Euk->Plas', 'EukVir->Plas', 'Plas->NonPlas'])], ignore_index=True)
# misclassified = pd.DataFrame(mistakes, columns=['Prok->Plas', 'ProkVir->Plas', 'Euk->Plas', 'EukVir->Plas', 'Plas->NonPlas']) # For plasmid
misclassified = pd.DataFrame(mistakes, columns=['Prok->ProkVir', 'Euk->ProkVir', 'EukVir->ProkVir', 'Plas->ProkVir', 'ProkVir->NonProkVir']) # For prokaryotic virus
misclassified.to_csv('perf_summary/misclassified_ppr_prokvir.csv', index=False)

# print(', '.join([str(i) for i in accs]))
# print(', '.join([str(i) for i in f1s]))
plasmid_summary_df.to_csv('perf_summary/ppr_plasmid.csv', index=False)
prokvirus_summary_df.to_csv('perf_summary/ppr_vir.csv', index=False)
multiclass_summary_df.to_csv('perf_summary/ppr_multiclass.csv', index=False)