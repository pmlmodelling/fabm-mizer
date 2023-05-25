#!/usr/bin/env python
from __future__ import print_function
import argparse
import numpy
import netCDF4


#Command line: python3 prepare_restart.py <offline_model.nc> --bmsource=<NEMO_restart_file>
parser = argparse.ArgumentParser()
parser.add_argument('target') #Restart file in which to add fish data too
parser.add_argument('--sourcevar', default='Y3_c')
parser.add_argument('--bmsource')  # OFfline model nc file
parser.add_argument('--bmscale', type=float, default=1./0.12)  # Data comes in g WM. Convert to mmol C. 0.012 g C/mmol C = 0.12 g WM/mmol C
parser.add_argument('-N', type=int, default=100)
args = parser.parse_args()

assignments = dict([('fish_c%i' % (i+1), 0.) for i in range(args.N)])
assignments['fish_landings'] = 0
prefixes = 'TRB', 'TRN'
prefixes = 'fabm_st2Db', 'fabm_st2Dn'

if args.bmsource is not None:
   with netCDF4.Dataset(args.bmsource) as nc:
      data = nc.variables['Nw_final'][...] * args.bmscale
      assert data.ndim == 3
     # assert data.shape[0] == 1
     # assert data.shape[3] == args.N
      print('Source biomass array has dimensions %i x %i' % data.shape[0:2])
      for i in range(args.N):
         assignments['fish_c%i' % (i+1)] = data[ :, :, i]

with netCDF4.Dataset(args.target, 'r+') as nc:
   print('%s contains the following variables:' % args.target)
   #for name in sorted(nc.variables.keys()):
   #   print('  %s' % name)
   for prefix in prefixes:
      ncsrc = nc.variables[prefix + args.sourcevar]
      for name, value in assignments.items():
         if prefix + name in nc.variables:
            nctgt = nc.variables[prefix + name]
         else:
            nctgt = nc.createVariable(prefix + name, ncsrc.dtype, ncsrc.dimensions)
         print('Writing %s (min = %.4g, max = %.4g, mean = %.4g)...' % (prefix + name, numpy.min(value), numpy.max(value), numpy.mean(value)))
         nctgt[0,:] = value
