import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

from scipy.stats import wilcoxon

results = [pd.read_csv('results/results_mu0_delta0.csv'), pd.read_csv('results/results_mu0_delta0.005.csv'), pd.read_csv('results/results_mu0.005_delta0.csv'), pd.read_csv('results/results_mu0.005_delta0.005.csv')]

xticks = [f'DS_{i}' for i in range(1, 21)]

name = ['Euk', 'EukVir', 'Plasmid', 'Prok', 'ProkVir']

for i in range(5):
    out = np.concatenate([result.iloc[:, i].to_numpy()[:, None] for result in results], axis=1)
    print(name[i])
    print(wilcoxon(out[:, 0], out[:, 1], alternative='greater'))
    print(wilcoxon(out[:, 0], out[:, 2], alternative='greater'))
    print(wilcoxon(out[:, 0], out[:, 3], alternative='greater'))
    print()
    # df = pd.DataFrame(out, columns=['mu=0, delta=0', 'mu=0, delta=0.05', 'mu=0.05, delta=0', 'mu=0.05, delta=0.05'], index=xticks)
    # ax = df.plot(kind='bar', width=0.8,)
    # ax.set_ylabel("F1 score")
    # plt.xticks(rotation=45)
    # ax.legend(loc='lower left')
    # ax.grid(axis='y')
    # ax.set_axisbelow(True)
    # ax.margins(0.01)