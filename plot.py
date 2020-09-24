try:
    import matplotlib.pyplot as plt
except:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

from argparse import ArgumentParser
import pandas as pd
import numpy as np


def load_csv(filename, dtype=float):
    df = pd.read_csv(filename, header=None)
    try:
        df = df.astype(dtype)
        data = df.values
    except Exception:
        print('[Warning ({})] "*+" exists, and replace to 0'.format(filename))
        for i in range(10, 20):
            df = df.replace('*' * i, 0)
        df = df.astype(dtype)
        data = df.values
    return data


parser = ArgumentParser()
parser.add_argument('skips', type=int, default=5)
args = parser.parse_args()


rhoes = load_csv('rhoe.csv')
# rhois = load_csv('rhoi.csv')
rhois = [np.min(rhoes), np.max(rhoes)]
phis = load_csv('phi.csv')
exs = load_csv('ex.csv')

min_val = np.min([np.min(rhoes), np.min(rhois)])
max_val = np.max([np.max(rhoes), np.max(rhois)])

fig = plt.figure()
ax1 = plt.subplot(3, 1, 1)
ax2 = plt.subplot(3, 1, 2)
ax3 = plt.subplot(3, 1, 3)
for i in range(0, len(rhoes), args.skips):
    print('\r{}'.format(i), end='')
    ax1.cla()
    ax1.set_ylim(min_val, max_val)
    ax1.plot(rhoes[i, :], label='rhoe')
    # ax1.plot(rhois[i, :], label='rhoi')
    ax1.legend()
    ax2.cla()
    ax2.set_ylim(np.min(phis), np.max(phis))
    ax2.plot(phis[i, :], label='phi')
    ax2.legend()
    ax3.cla()
    ax3.set_ylim(np.min(exs), np.max(exs))
    ax3.plot(exs[i, :], label='ex')
    ax3.legend()
    plt.pause(0.1)
