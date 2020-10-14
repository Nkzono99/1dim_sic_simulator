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
parser.add_argument('--datadir', '-d', default='.')
parser.add_argument('--show', action='store_true', default=False)
args = parser.parse_args()


params = f90nml.read(os.path.join(args.datadir, 'plasma.in'))

# シミュレーションパラメータ ###
nspec = params['plasma']['nspec']
dt = params['simulation']['dt']
nsteps = params['simulation']['nsteps']
output_steps = params['output']['output_steps']
output_skips = nsteps // output_steps
dt = dt * output_skips
dx = params['simulation']['dx']
ngrid = params['simulation']['ngrid']
sim_size = dx * ngrid
#######################

datadir = os.path.join(args.datadir, 'simplex')
if not os.path.exists(datadir):
    os.mkdir(datadir)

for ispec in range(nspec):
    i = 0
    path = os.path.join(args.datadir, 'simplex{}.csv'.format(ispec+1))
    lines = open(path, 'r', encoding='utf-8').readlines()
    fig = None
    for line in lines:
        if line.startswith('Time'):
            i = 0
            time = int(line.replace('Time', '').strip())
            if fig is not None:
                save_path = os.path.join(
                    datadir, 'simplex{}_{:04}.png'.format(ispec+1, time))
                fig.savefig(save_path)
                plt.close(fig)
            fig = plt.figure()
            plt.title('simplex{} {}steps'.format(ispec+1, time))
            continue
        i += 1
        y = np.ones((2, ))*i
        px1, px2 = np.array(list(map(float, line.split(','))))
        ncycle1, lpx1 = divmod(px1, sim_size)
        ncycle2, lpx2 = divmod(px2, sim_size)
        dcycles = abs(px1 - px2) // sim_size

        if dcycles == 0:
            color = 'blue'
        elif dcycles == 1:
            color = 'green'
        else:
            color = 'red'

        plt.plot([lpx1, lpx2], y, 'o', color=color)
        if ncycle1 == ncycle2:
            plt.plot([lpx1, lpx2], y, '-', color=color)
        else:
            if abs(px1 - px2) > sim_size:
                plt.plot([0, sim_size], y, '-', color=color)
            elif px1 < px2:
                plt.plot([lpx1, sim_size], y, '-', color=color)
                plt.plot([0, lpx2], y, '-', color=color)
            else:
                plt.plot([lpx2, sim_size], y, '-', color=color)
                plt.plot([0, lpx1], y, '-', color=color)

