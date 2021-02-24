import argparse

import matplotlib.pyplot
import netCDF4

def compare(offline_path, gotm_path, figure_path=None):
    with netCDF4.Dataset(offline_path) as ncoff, netCDF4.Dataset(gotm_path) as ncgotm:
        time1 = netCDF4.num2date(ncoff['time'][:], ncoff['time'].units, only_use_cftime_datetimes=False)
        time2 = netCDF4.num2date(ncgotm['time'][:], ncgotm['time'].units, only_use_cftime_datetimes=False)
        bm_tot1 = ncoff['biomass'][:]
        h1 = ncoff['h'][:]
        bm_tot2 = ncgotm['fish_c_tot'][:, 0, 0]
        h2 = ncgotm['fish_w_integrator_result'][:, 0, 0]**2 / ncgotm['fish_w2_integrator_result'][:, 0, 0]

    fig = matplotlib.pyplot.figure(figsize=(8, 6))
    ax = fig.add_subplot(211)
    ax.plot(time1, h1, label='offline')
    ax.plot(time2, h2, label='online-uncoupled')
    ax.set_ylabel('effective depth range (m)')
    ax.autoscale(enable=True, axis='x', tight=True)
    ax.grid()
    ax.legend()

    ax = fig.add_subplot(212)
    ax.plot(time1, bm_tot1, label='offline')
    ax.plot(time2, bm_tot2, label='online-uncoupled')
    ax.set_ylabel('total fish biomass (g WM/m2)')
    ax.autoscale(enable=True, axis='x', tight=True)
    ax.grid()
    ax.legend()
    if figure_path is None:
        matplotlib.pyplot.show()
    else:
        fig.savefig(figure_path)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('offline_path')
    parser.add_argument('gotm_path')
    parser.add_argument('--out')
    args = parser.parse_args()
    compare(args.offline_path, args.gotm_path, args.out)