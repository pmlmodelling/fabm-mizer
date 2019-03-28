import os
import argparse

import numpy
import netCDF4
from matplotlib import pyplot

parser = argparse.ArgumentParser()
parser.add_argument('source')
parser.add_argument('--minw', type=float, default=1e-16)
parser.add_argument('--maxw', type=float, default=1e6)
args = parser.parse_args()

size_grid_bnds = numpy.linspace(numpy.log10(args.minw), numpy.log10(args.maxw), 1000)
size_grid_ctr = (size_grid_bnds[1:] + size_grid_bnds[:-1]) / 2
dw = 10.**size_grid_bnds[1:] - 10.**size_grid_bnds[:-1]
w = 10.**size_grid_ctr
biomass_ctr = None

def add_group(data, minw, maxw):
   global biomass_ctr
   log10_minw = numpy.log10(minw)
   log10_maxw = numpy.log10(maxw)
   fractions = numpy.zeros_like(size_grid_ctr)
   for i in xrange(size_grid_ctr.size):
      fractions[i] = max(0., min(size_grid_bnds[i + 1], log10_maxw) - max(size_grid_bnds[i], log10_minw))/(log10_maxw - log10_minw)
   curbiomass = (h * data/12).sum(axis=1)[:, numpy.newaxis] * fractions[numpy.newaxis, :]
   if biomass_ctr is None:
      biomass_ctr = curbiomass
   else:
      biomass_ctr += curbiomass

with netCDF4.Dataset(args.source) as nc:
   nctime = nc.variables['time']
   dt = netCDF4.num2date(nctime[:], nctime.units)
   h = nc.variables['h'][:, :, 0, 0]
   add_group(nc.variables['P1_c'][:, :, 0, 0], 4.188e-9, 4.188e-6)
   add_group(nc.variables['P2_c'][:, :, 0, 0], 4.188e-12, 4.188e-9)
   add_group(nc.variables['P3_c'][:, :, 0, 0], 4.188e-15, 4.188e-12)
   add_group(nc.variables['P4_c'][:, :, 0, 0], 4.188e-9, 4.188e-6)
   add_group(nc.variables['Z4_c'][:, :, 0, 0], 4.188e-6, 1e-3)
   add_group(nc.variables['Z5_c'][:, :, 0, 0], 4.188e-9, 4.188e-6)
   add_group(nc.variables['Z6_c'][:, :, 0, 0], 4.188e-12, 4.188e-9)
   offset = nc.variables['fish_pelagic_size_spectrum_offset'][:, 0, 0]
   slope = nc.variables['fish_pelagic_size_spectrum_slope'][:, 0, 0]
   fish = []
   i = 1
   while True:
      name = 'fish_c%i' % i
      if name not in nc.variables:
         break
      print 'Reading %s...' % name
      fish.append(nc.variables[name][:, 0, 0])
      i = i + 1
   fish_log10_w = numpy.linspace(-3, 6, len(fish))
   fish_log10_w_bnds = numpy.linspace(fish_log10_w[0] - (fish_log10_w[1] - fish_log10_w[0])*0.5, fish_log10_w[-1] + (fish_log10_w[1] - fish_log10_w[0])*0.5, len(fish)+1)
   fish_dw = 10.**fish_log10_w_bnds[1:] - 10.**fish_log10_w_bnds[:-1]
   fish_w = 10.**fish_log10_w

fish = numpy.vstack(fish)

biomass_ctr /= dw[numpy.newaxis, :]
fish /= fish_dw[:, numpy.newaxis]

fig = pyplot.figure()
ax = fig.gca()
plankton_spectrum, = ax.loglog(w, biomass_ctr[0, :]/dw, '-')
fitted_plankton_spectrum, = ax.loglog((1e-7, 1e-3), (numpy.exp(offset[0] + slope[0] * numpy.log(1e-7)), numpy.exp(offset[0] + slope[0] * numpy.log(1e-3))), '--r')
fish_spectrum, = ax.loglog(fish_w, fish[:, 0], '-g')
ax.set_xlim(10.**size_grid_bnds[0], 10.**size_grid_bnds[-1])
startspectrum = numpy.exp(offset + slope * numpy.log(1e-3))
ax.set_ylim(startspectrum.min()*1e-9, max(biomass_ctr.max(), fish.max()))
ax.set_xlabel('wet mass (g)')
ax.set_ylabel('mass density (g/g)')
ax.grid()

def plot(i):
   plankton_spectrum.set_ydata(biomass_ctr[i, :])
   fitted_plankton_spectrum.set_ydata((numpy.exp(offset[i] + slope[i] * numpy.log(1e-7)), numpy.exp(offset[i] + slope[i] * numpy.log(1e-3))))
   fish_spectrum.set_ydata(fish[:, i])
   ax.set_title(dt[i].strftime('%Y-%m-%d'))

if not os.path.isdir('ani_spectrum'):
   os.makedirs('ani_spectrum')
for itime in xrange(biomass_ctr.shape[0]):
   if itime % 100 == 0:
      print 'Processing time %i...' % itime
   plot(itime)
   fig.savefig('ani_spectrum/%05i.png' % itime, dpi=96)

#pyplot.show()
