import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

TOTAL_SEQ_NUM = 1000

PROK_EUK_RATIO = np.array([
    [9, 1],
    [7, 3],
    [5, 5],
    [3, 7],
    [1, 9],
], dtype=float)
PROK_EUK_RATIO /= PROK_EUK_RATIO.sum(axis=1, keepdims=True)

PROK_PROKVIR_PLASMID_RATIO = np.array([
    [5, 1, 1],
    [4, 1, 1],
    [3, 1, 1],
    [2, 1, 1],
], dtype=float)
PROK_PROKVIR_PLASMID_RATIO /= PROK_PROKVIR_PLASMID_RATIO.sum(axis=1, keepdims=True)

EUK_EUKVIR_RATIO = np.array([
    [5, 1],
    [4, 1],
    [3, 1],
    [2, 1],
], dtype=float)
EUK_EUKVIR_RATIO /= EUK_EUKVIR_RATIO.sum(axis=1, keepdims=True)

ratio = np.concatenate([np.concatenate([PROK_EUK_RATIO[i, 0] * PROK_PROKVIR_PLASMID_RATIO, PROK_EUK_RATIO[i, 1] * EUK_EUKVIR_RATIO], axis=1) for i in range(PROK_EUK_RATIO.shape[0])], axis=0) * 100

tiara_acc = [0.989010989010989, 0.989, 0.989, 0.993, 0.979, 0.978021978021978, 0.984, 0.988, 0.9669669669669669, 0.964964964964965, 0.966, 0.979, 0.953, 0.966, 0.962, 0.952, 0.92992992992993, 0.945054945054945, 0.945, 0.941]
tiara_f1 = [0.9299363057324841, 0.9271523178807948, 0.9219858156028369, 0.9448818897637796, 0.9561586638830899, 0.9521739130434782, 0.9634703196347033, 0.9690721649484536, 0.9588014981273408, 0.9542483660130718, 0.9527777777777778, 0.9674418604651163, 0.9580731489741303, 0.96875, 0.9624505928853755, 0.945945945945946, 0.9510489510489509, 0.9602888086642599, 0.9575289575289575, 0.9483814523184603,]

whok_acc = [0.9585062240663901, 0.9611344537815126, 0.9512448132780082, 0.9473140495867769, 0.9332591768631813, 0.9394618834080718, 0.9144444444444444, 0.9058956916099773, 0.9076354679802956, 0.882640586797066, 0.8824969400244798, 0.8544776119402985, 0.8688741721854305, 0.8481182795698925, 0.8297872340425532, 0.789261744966443, 0.8418740849194729, 0.8179148311306902, 0.7800586510263929, 0.7246376811594203,]
whok_f1 = [0.6875, 0.6725663716814159, 0.656934306569343, 0.6277372262773723, 0.8333333333333333, 0.8411764705882353, 0.7843137254901961, 0.7067137809187279, 0.8648648648648649, 0.826086956521739, 0.8160919540229886, 0.7557411273486431, 0.8745247148288973, 0.8515111695137977, 0.8293333333333334, 0.7617602427921093, 0.8943248532289629, 0.8729508196721312, 0.8414376321353066, 0.7850678733031675,]

dmf_acc = [0.998001998001998, 0.998, 1.0, 0.998, 0.996, 0.999000999000999, 0.998, 0.998, 0.995995995995996, 0.995995995995996, 0.995, 0.995, 0.994, 0.991, 0.996, 0.991, 0.992992992992993, 0.993006993006993, 0.991, 0.987,]
dmf_f1 = [0.9878048780487805, 0.9876543209876543, 1.0, 0.9850746268656716, 0.992, 0.997920997920998, 0.9955357142857144, 0.995, 0.9951923076923077, 0.9949874686716792, 0.9933065595716198, 0.992503748125937, 0.9948453608247423, 0.9919282511210762, 0.9961977186311786, 0.9903743315508022, 0.9953364423717521, 0.9951287404314543, 0.9933184855233853, 0.9891213389121339,]

tiara_acc = np.array(tiara_acc)
tiara_f1 = np.array(tiara_f1)
whok_acc = np.array(whok_acc)
whok_f1 = np.array(whok_f1)
dmf_acc = np.array(dmf_acc)
dmf_f1 = np.array(dmf_f1)

xticks = [f'DS_{i}' for i in range(1, 21)]

ratio_df = pd.DataFrame(ratio, columns=['Prok', 'ProkVirus', 'Plasmid', 'Euk', 'EukVirus'], index=xticks)
ax = ratio_df.plot(kind='bar',
    stacked=True,
    width=0.8,
    color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd'],
    # xlabel="Dataset",
    ylabel="Percentage (%)"
    )
plt.xticks(rotation=45)
ax.grid(axis='y')
ax.set_axisbelow(True)
ax.margins(0.01)
plt.savefig('figures/ds_composition.pdf', bbox_inches='tight')

acc_df = pd.DataFrame(columns=['Tiara', 'Whok', 'DMF'], index=xticks)
acc_df['Tiara'] = tiara_acc
acc_df['Whok'] = whok_acc
acc_df['DMF'] = dmf_acc

ax = acc_df.plot(kind='bar',
    width=0.8,
)
# ax.set_xlabel("Dataset")
ax.set_ylabel("Accuracy")
plt.xticks(rotation=45)
ax.legend(loc='lower left')
ax.grid(axis='y')
ax.set_axisbelow(True)
ax.margins(0.01)
plt.savefig('figures/euk_acc.pdf', bbox_inches='tight')

f1_df = pd.DataFrame(columns=['Tiara', 'Whok', 'DMF'], index=xticks)
f1_df['Tiara'] = tiara_f1
f1_df['Whok'] = whok_f1
f1_df['DMF'] = dmf_f1

ax = f1_df.plot(kind='bar',
    width=0.8, 
)
# ax.set_xlabel("Dataset")
ax.set_ylabel("F1 Score")
plt.xticks(rotation=45)
ax.set_ylim(top=1)
ax.legend(loc='lower left')
ax.grid(axis='y')
ax.set_axisbelow(True)
ax.margins(0.01)
plt.savefig('figures/euk_f1.pdf', bbox_inches='tight')

# import seaborn as sns
# sns.catplot(data=acc_df, kind='bar',
#     height=4,
#     aspect=2
#     )