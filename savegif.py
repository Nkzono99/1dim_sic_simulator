try:
    import matplotlib.pyplot as plt
except:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

import matplotlib.animation as anm
import pandas as pd
import numpy as np
from argparse import ArgumentParser
import f90nml
import os


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
parser.add_argument('--output', '-o', default='output.gif')
parser.add_argument('--input', '-i', default='plasma.in')
parser.add_argument('--datadir', '-d', default='.')
args = parser.parse_args()


params = f90nml.read(args.input)

# シミュレーションパラメータ ###
nspec = params['plasma']['nspec']
dt = params['simulation']['dt']
nsteps = params['simulation']['nsteps']
output_steps = params['output']['output_steps']
output_skips = nsteps // output_steps
dt = dt * output_skips
dx = params['simulation']['dx']
#######################


rhos = np.array([load_csv(os.path.join(
    args.datadir, 'rho{}.csv'.format(i + 1))) for i in range(nspec)])
phis = load_csv(os.path.join(args.datadir, 'phi.csv'))
exs = load_csv(os.path.join(args.datadir, 'ex.csv'))

fig = plt.figure()
ax1 = plt.subplot(3, 1, 1)
ax2 = plt.subplot(3, 1, 2)
ax3 = plt.subplot(3, 1, 3)


def update(i):
    i = i * args.skips
    istep = i * output_skips
    t = i * dt
    # print('\r{:16.8f} s ({} steps)'.format(t, istep), end='')
    fig.suptitle('{:16.8f} s ({} steps)'.format(t, istep))
    ax1.cla()
    ax1.set_ylim(np.min(rhos), np.max(rhos))
    for ispec in range(nspec):
        ax1.plot(rhos[ispec, i, :], label='rho{}'.format(ispec+1))
    ax1.legend()
    ax2.cla()
    ax2.set_ylim(np.min(phis), np.max(phis))
    ax2.plot(phis[i, :], label='phi')
    ax2.legend()
    ax3.cla()
    ax3.set_ylim(np.min(exs), np.max(exs))
    ax3.plot(exs[i, :], label='ex')
    ax3.legend()


frames = rhos.shape[1] // args.skips
ani = anm.FuncAnimation(
    fig,
    update,
    interval=100,
    frames=frames
)
ani.save(args.output)
