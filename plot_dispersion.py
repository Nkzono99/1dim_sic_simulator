try:
    import matplotlib.pyplot as plt
except:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import f90nml
from argparse import ArgumentParser
import os


c = 2.99792458e8


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


def plot_spectram(data, d, label=None):
    data_fft = np.fft.fft(data)
    N = len(data_fft)
    freq = np.linspace(0, 1.0/d, N)
    amp = np.abs(data_fft/(N/2))
    plt.plot(freq, amp, label=label)


def plot_2dmap(data, mesh=None, vmin=None, vmax=None, nsigma=3, label=None):
    if vmin == 'fit':
        vmin = np.min(data)
    if vmax == 'fit':
        vmax = np.max(data)
    if isinstance(vmin, str) and vmin.startswith('tile'):
        percent = int(vmin.replace('tile', ''))
        vmin = np.percentile(data.reshape(-1), percent)
    if isinstance(vmax, str) and vmax.startswith('tile'):
        percent = int(vmax.replace('tile', ''))
        vmax = np.percentile(data.reshape(-1), percent)
    if vmin is None:
        mean = np.mean(data)
        var = np.var(data)
        vmin = mean - var * nsigma
    if vmax is None:
        mean = np.mean(data)
        var = np.var(data)
        vmax = mean + var * nsigma
    if mesh is None:
        X, Y = np.meshgrid(np.arange(data.shape[1]), np.arange(data.shape[0]))
    else:
        X, Y = mesh
    print(np.min(data), np.max(data), np.mean(data))
    plt.pcolormesh(X, Y, data, vmin=vmin, vmax=vmax, label=label, cmap='jet')
    plt.colorbar(orientation='vertical')


def myfft2(data, dx, dt):
    """
    二次元方向のフーリエ変換を行う.

    データの中心が周波数領域の原点となる.
              w ^
                |
                |
                |
                | 
    -----------------------> k
         (0, 0) |
                |
                |
                |
    Parameters
    ----------
    data: 2d-array(shape: (nt, nx))
        データ
    dx: x方向の幅
    dt: t方向の幅
    """
    nt, nx = data.shape

    x = np.linspace(0, (nx - 1) * dx, nx)
    t = np.linspace(0, (nt - 1) * dt, nt)

    nyq_k = 1 / (2 * dx)
    nyq_f = 1 / (2 * dt)
    dk = 1 / (nx * dx)
    df = 1 / (nt * dt)
    k = np.linspace(-nyq_k, nyq_k-dk, nx)
    f = np.linspace(-nyq_f, nyq_f-df, nt)

    dispersion = np.zeros_like(data, dtype=complex)

    for i1 in range(nx):
        for j1 in range(nt):
            arr_k = np.repeat((k[i1] * x).reshape(1, -1), nt, axis=0)
            arr_f = np.repeat((f[j1] * t).reshape(-1, 1), nx, axis=1)
            dispersion[j1, i1] = np.sum(
                data * np.exp(-1j * (2 * np.pi) * (arr_k + arr_f))) * dx * dt
            # 上記のコードは以下のコードの最適化ver
            # for i2 in range(nx):
            #     for j2 in range(nt):
            #         dispersion[j1, i1] += data[j2, i2] * np.exp(-1j * (2 * np.pi) * (k[i1] * x[i2] + f[j1] * t[j2])) * dx * dt


def plot_dispersion(data,
                    dx=1,
                    dt=1,
                    label=None,
                    vmin=None,
                    vmax=None,
                    nsigma=3,
                    xstart=0.5,
                    xend=1.0,
                    tstart=0.5,
                    tend=1.0):
    dispersion = np.fft.fftshift(np.fft.fft2(data))
    # dispersion = myfft2(data, dx=dx, dt=dt)

    amp = np.abs(dispersion) / (data.shape[0] / 2) / (data.shape[1] / 2)
    # amp = np.log(amp)

    nt, nx = data.shape
    t_start, t_end = int(tstart * nt), int(tend * nt)
    x_start, x_end = int(xstart * nx), int(xend * nx)
    amp_focused = amp[t_start:t_end, x_start:x_end]

    # nyq_k = 1 / (2 * dx) * 2 * np.pi
    # nyq_w = 1 / (2 * dt) * 2 * np.pi
    # dk = 2 * nyq_k / nx
    # dw = 2 * nyq_w / nt
    # k = np.linspace(0, nyq_k-dk, nx // 2)[:x_end-x_start]
    # w = np.linspace(0, nyq_w-dw, nt // 2)[:t_end-t_start]
    k = np.fft.fftfreq(nx, d=dx)[:x_end-x_start] * 2 * np.pi
    w = np.fft.fftfreq(nt, d=dt)[:t_end-t_start] * 2 * np.pi
    w /= wpe
    mesh = np.meshgrid(k, w)

    plot_2dmap(amp_focused,
               mesh=mesh,
               vmin=vmin,
               vmax=vmax,
               nsigma=nsigma,
               label=label)
    plt.xlabel('k')
    plt.ylabel('w / wpe')

    gamma = 3
    electron_disp = wpe * np.sqrt(1 + gamma * k ** 2 * debye ** 2)
    electron_disp /= wpe
    ion_disp = np.sqrt(gamma * k ** 2 * Ti * kB / mi + gamma *
                       k ** 2 * Te * kB / mi / (1 + gamma * k ** 2 * debye ** 2))
    ion_disp /= wpe
    if not args.notheory:
        plt.plot(k, electron_disp, color='red', label='electron plasma wave')
        plt.plot(k, ion_disp, color='yellow', label='ion sound wave')
        plt.axhline(y=wpe/wpe, xmin=0, xmax=1, color='white', label='wpe')
        # plt.axhline(y=wpi/wpe, xmin=0, xmax=1, color='green', label='wpi')

    plt.ylim(np.min(w), np.max(w))


parser = ArgumentParser()
parser.add_argument('--xstart', '-xs', type=float, default=0.5)
parser.add_argument('--xend', '-xe', type=float, default=1)
parser.add_argument('--tstart', '-ts', type=float, default=0.5)
parser.add_argument('--tend', '-te', type=float, default=1)
parser.add_argument('--nsamples_t', '-nt', type=int, default=1000)
parser.add_argument('--nsamples_x', '-nx', type=int, default=4096)
parser.add_argument('--input', '-i', default='plasma.in')
parser.add_argument('--output', '-o', default=None)
parser.add_argument('--datadir', '-d', default='.')
parser.add_argument('--notheory', action='store_true', default=False)
args = parser.parse_args()

params = f90nml.read(args.input)

# シミュレーションパラメータ ###
dt = params['simulation']['dt']
nsteps = params['simulation']['nsteps']
output_steps = params['output']['output_steps']
dt = dt * nsteps / output_steps
dx = params['simulation']['dx']
#######################

t_start = 2000
t_end = output_steps
nsamples = args.nsamples_t

x_start = 0
x_end = params['simulation']['ngrid']
nsamples_x = min(args.nsamples_x, x_end - x_start)

skips = (t_end - t_start) // nsamples
skips_x = (x_end - x_start) // nsamples_x

exs = load_csv(os.path.join(args.datadir, 'ex.csv'))

debye = params['plasma']['lambda']
Te = params['plasma']['Ts'][0]
Ti = params['plasma']['Ts'][1]
kB = 1.380649e-23
me = 9.10938356e-31
mi = me * params['plasma']['m_ratio'][1]

wpe = 1 / debye * np.sqrt(kB * Te / me)
wpi = 1 / debye * np.sqrt(kB * Ti / mi)

plt.figure()

plot_dispersion(exs[t_start:t_end:skips, x_start:x_end:skips_x],
                dt=dt*skips,
                dx=dx*skips_x,
                vmin='tile50',
                vmax='tile100',
                nsigma=1e2,
                label='ex',
                xstart=args.xstart,
                xend=args.xend,
                tstart=args.tstart,
                tend=args.tend)
# plt.axhline(y=wpi, xmin=0, xmax=1, color='yellow', label='wpi')
plt.legend()

plt.plot()
plt.title('ex')

if args.output is None:
    plt.show()
else:
    plt.savefig(args.output)
