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
parser.add_argument('--output', '-o', default=None)
parser.add_argument('--input', '-i', default='plasma.in')
parser.add_argument('--datadir', '-d', default='.')
args = parser.parse_args()


params = f90nml.read(args.input)

# シミュレーションパラメータ ###
dt = params['simulation']['dt']
nsteps = params['simulation']['nsteps']
output_steps = params['output']['output_steps']
output_skips = nsteps // output_steps
dt = dt * output_skips
dx = params['simulation']['dx']
#######################


kenergy = load_csv(os.path.join(args.datadir, 'kinetic_energy.csv'))
esenergy = load_csv(os.path.join(args.datadir, 'es_energy.csv'))
energy = kenergy + esenergy

# kenergy -= kenergy[0]
# esenergy -= esenergy[0]
# energy -= energy[0]

nt = len(kenergy)

t = np.linspace(0, dt * nt, nt)

plt.figure(figsize=(10, 9))

plt.subplot(3, 1, 1)
plt.plot(t, energy, label='sum')
plt.xlabel('t [s]')
plt.ylabel('energy [J]')
plt.legend()

plt.subplot(3, 1, 2)
plt.plot(t, kenergy, label='kinetic')
plt.xlabel('t [s]')
plt.ylabel('energy [J]')
plt.legend()

plt.subplot(3, 1, 3)
plt.plot(t, esenergy, label='electrostatic')
plt.xlabel('t [s]')
plt.ylabel('energy [J]')
plt.legend()

if args.output is None:
    plt.show()
else:
    plt.savefig(args.output)
