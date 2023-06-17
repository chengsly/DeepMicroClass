import numpy as np
import pandas as pd
import os
import re
from sklearn.metrics import f1_score, matthews_corrcoef, confusion_matrix, average_precision_score, balanced_accuracy_score

RESULT_DIR = 'data/result_others/result_tiara'

results_fn = os.listdir(RESULT_DIR)
results_fn = [f for f in results_fn if not f.startswith('log')]
results_fn.sort()

TARGET_DIR = 'data/sampled_sequence_target'

target_fn = os.listdir(TARGET_DIR)
target_fn.sort()

def construct_target(nums):
    order = [[3], [4], [2], [0], [1]]
    return np.array(sum([order[i]*nums[i] for i in range(5)], []))

def construct_result(df):
    result = []
    prok = set(('archaea', 'bacteria', 'prokarya'))
    euk = set(('eukarya', 'organelle'))
    for i, row in df.iterrows():
        if row[0] in euk:
            result.append(1)
        # elif row[0] in prok:
            # result.append(0)
        else:
            result.append(0)
    return np.array(result)

def construct_euk_array(array):
    euk_array = np.zeros((array.shape[0], 2))
    euk_array[np.logical_or.reduce((array==1, array==2, array==3, array==4)), 0] = 1
    euk_array[array==0, 1] = 1
    return euk_array

mistakes = []

summary_df = pd.DataFrame(columns=['filename', 'f1_score', 'accuracy'])

for f in results_fn:
    nums = re.findall(r'\d+', f)[1:] # get the numbers in the filename
    result = pd.read_table(os.path.join(RESULT_DIR, f), index_col=0) # read the result
    
    index = list(result.index)
    index = [i.split()[0] for i in index] # get the sequence name
    result.index = index # set the sequence name as index
    
    target_pkl = pd.read_pickle(os.path.join(TARGET_DIR, f'{f[:-4]}.fa.pkl')) # read the target
    target = np.array([target_pkl.loc[i, 0] for i in index]).flatten() # get the target sequence name
    missed = len(target) - len(target_pkl) # get the number of missed sequences

    # add the missed sequences to the target
    missed_label = []
    for i in target_pkl.index:
        if i not in result.index:
            missed_label.append(target_pkl.loc[i, 0])
    missed_label = np.array(missed_label).flatten()
    target = np.concatenate((target, missed_label))

    # convert the target to binary
    target_binary = np.zeros(len(target))
    target_binary[target==0] = 1
    # target_binary[np.logical_or(target==0, target==1)] = 1
    # target_binary = np.concatenate((target_binary, np.ones(missed)))

    result = construct_result(result)
    result = np.concatenate((result, np.zeros(missed))) # add the missed sequences as 0 to the result

    try:
        acc = (result==target_binary).sum() / len(target_binary)
        # f1 = f1_score(target_binary, result, average='weighted')
        f1 = f1_score(target_binary, result)
        # f1 = matthews_corrcoef(target_binary, result)
    except:
        print(f)

    summary_df = pd.concat([summary_df, pd.DataFrame([['_'.join(nums), f1, acc]], columns=['filename', 'f1_score', 'accuracy'])])

    mistake = [(target[result==1]==3).sum(), (target[result==1]==4).sum(), (target[result==0]==0).sum(), (target[result==1]==1).sum(), (target[result==1]==2).sum()]
    mistakes.append(mistake)

    # print(f)
    # print(f'Acc: {acc}\tF1: {f1}\tUnknown: {(result==-1).sum()/len(result)}')
    # print(f'{f1:.2f}, Unknown: {(result==-1).sum()/len(result):.2f}')
    # dmf_result = pd.read_table(os.path.join('DeepMicrobeFinder_results', f'{f[:-4]}.fa_pred_one-hot_hybrid.txt'), index_col=0)
    # dmf_index = [i.split()[0] for i in dmf_result.index]
    # dmf_result.index = dmf_index
    # dmf_result = np.array([dmf_result.loc[i].values for i in index]).argmax(axis=1)

    # dmf_euk = construct_euk_array(dmf_result)
    # dmf_euk = dmf_euk.argmax(axis=1)

    # dmf_f1 = f1_score(target_binary, dmf_euk)
    # dmf_f1 = matthews_corrcoef(target_binary, dmf_euk)
    # dmf_target = construct_euk_array(target_pkl.iloc[:, 0]).argmax(axis=1)
    # dmf_f1 = f1_score(dmf_target, dmf_euk)
    # print(f'Tiara F1: {f1:.2f}\t DMF F1: {dmf_f1:.2f}')
    # print(f'{f1:.2f}\t{dmf_f1:.2f}')
    # print(acc, end=', ')
    # print(f1, end=', ')
    # print(f1_score(target_binary, result), balanced_accuracy_score(target_binary, result), average_precision_score(target_binary, result), matthews_corrcoef(target_binary, result), acc)
    # print(confusion_matrix(target_binary, result, normalize='true'))
    # print('')
# a = [0.9818181818181818,0.9815950920245399,0.9868421052631579,0.9852941176470589,0.9940119760479043,0.997920997920998,0.9977827050997783,0.9925187032418954,0.996389891696751,0.9962546816479401,0.9919786096256684,0.990990990990991,0.9974248927038627,0.9946332737030412,0.9961977186311786,0.9946409431939979,0.9939879759519038,0.9958275382475661,0.9940740740740741,0.9925062447960034,]
# b = [0.9299363057324841,0.9271523178807948,0.9219858156028369,0.9448818897637796,0.9561586638830899,0.9521739130434782,0.9634703196347033,0.9690721649484536,0.9588014981273408,0.9542483660130718,0.9527777777777778,0.9674418604651163,0.9580731489741303,0.96875,0.9624505928853755,0.945945945945946,0.9510489510489509,0.9602888086642599,0.9575289575289575,0.9483814523184603,]
# from scipy.stats import mannwhitneyu
# mannwhitneyu(a, b)

# summary_df.to_csv('perf_summary/tiara.csv', index=False)

misclassified = pd.DataFrame(mistakes, columns=['Prok->Euk', 'ProkVir->Euk', 'Euk->Non-Euk', 'EukVir->Euk', 'Plas->Euk'])
misclassified.to_csv('perf_summary/misclassified_tiara.csv', index=False)
