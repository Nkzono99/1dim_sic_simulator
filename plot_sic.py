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
#######################


distance = load_csv(os.path.join(args.datadir, 'distance_between_tracers.csv'))
distance /= dx
nt = distance.shape[0]
t = np.linspace(0, dt * nt, nt)

fig = plt.figure()

for ispec in range(nspec):
    plt.subplot(nspec, 1, ispec+1)
    plt.plot(t, distance[:, ispec], label='p{}'.format(ispec+1))
    plt.legend()
    plt.xlabel('t [s]')
    plt.ylabel('mean distance / grid width')

if args.show:
    plt.show()
else:
    plt.savefig(os.path.join(args.datadir, 'distance.png'))
