from __future__ import print_function

import argparse
import glob
import os.path
import numpy
import netCDF4

#Example run command: python3 <source_data_file_path_(i.e AMM7 run output> <output_filename>.nc  --start=2000 --stop=2003

parser = argparse.ArgumentParser()
parser.add_argument('source') #Directory of baseline nemo-ersem data
parser.add_argument('target') #outputfile
parser.add_argument('--start', type=int, default=None)
parser.add_argument('--stop', type=int, default=None)
arguments = parser.parse_args()

#root = '/nerc/n01/n01/momme/AMM7-HINDCAST-v0'

assert os.path.isdir(arguments.source), 'First argument must be the path to a directory (%s is not)' % arguments.source

start_year = None
stop_year = None
for path in glob.glob(os.path.join(arguments.source, '????')):
    if os.path.isdir(path):
        try:
            year = int(os.path.basename(path))
            if start_year is None:
                start_year = year
                stop_year = year
            start_year = min(start_year, year)
            stop_year = max(stop_year, year)
        except:
            pass

if arguments.start is None:
    arguments.start = start_year
    print('Start year: %04i' % arguments.start)
if arguments.stop is None:
    arguments.stop = stop_year
    print('Stop year: %04i' % arguments.stop)

filename = 'amm7_1d_*_ptrc_T.nc'
filename_temp = 'amm7_1d_*_grid_T.nc'

variables = 'P1_c', 'P2_c', 'P3_c', 'P4_c', 'Z4_c', 'Z5_c', 'Z6_c','Y3_c','Y2_c','Y4_c'
ben_vars = 'Y3_c','Y2_c','Y4_c'
pel_vars='P1_c', 'P2_c', 'P3_c', 'P4_c', 'Z4_c', 'Z5_c', 'Z6_c'



compress = False
contiguous = False
chunk = True

def copyVariable(ncout, ncvar, dimensions=None, copy_data=True):
   copy_data = copy_data and dimensions is None
   if dimensions is None:
      dimensions = ncvar.dimensions
   kwargs = {}
   if chunk and 'time_counter' in dimensions:
      chunksizes = [1]*len(dimensions)
      chunksizes[dimensions.index('time_counter')] = ntime
      kwargs['chunksizes'] = chunksizes
   ncvar_out = ncout.createVariable(ncvar.name, ncvar.dtype, dimensions, zlib=compress, contiguous=contiguous, **kwargs)
   ncvar_out.units = ncvar.units
   ncvar_out.long_name = ncvar.long_name
   if copy_data:
      ncvar_out[...] = ncvar[...]
   return ncvar_out

ncvariables_out = None
iout = 0

month2length = (31, (28, 29), 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
def getInputs():
   for year in range(arguments.start, arguments.stop+1):
       for month in range(1, 13):
           paths = glob.glob(os.path.join(arguments.source, '%04i' % year, '%02i' % month, filename))
           assert len(paths) == 1, 'Multiple paths found: %s' % paths
           path = paths[0]
           paths_temp = glob.glob(os.path.join(arguments.source, '%04i' % year, '%02i' % month, filename_temp))
           assert len(paths_temp) == 1, 'Multiple paths found: %s' % paths_temp
           path_temp = paths_temp[0]
           print('Processing %s...' % path)
           nc, nc_temp = netCDF4.Dataset(path), netCDF4.Dataset(path_temp)
           nctime = nc.variables['time_counter']
           expected_length = month2length[month-1]
           if isinstance(expected_length, int):
              expected_length = (expected_length,)
           if nctime.size not in expected_length:
              print('Month %i has invalid length %i (expected %s)' % (month, nctime.size, expected_length))
              return
           yield nc, nc_temp, year, month

ntime = 0
for nc, nc_temp, year, month in getInputs():
   with nc, nc_temp:
      nctime = nc.variables['time_counter']
      ntime += nctime.size
print('Time count: %i' % ntime)
mode = 'r+' if os.path.isfile(arguments.target) else 'w'
mode = 'w'
with netCDF4.Dataset(arguments.target, mode, clobber=False, diskless=True, persist=True) as ncout:
   istart = 0
   if mode != 'w':
      nctime_out = ncout.variables['time_counter']
      nctemp_out = ncout.variables['votemper']
      ncw_int_out = ncout.variables['bm_int']
      ncw2_int_out = ncout.variables['bm2_int']
      istart = nctime_out.size-2
      ncvariables_out = [ncout.variables[variable] for variable in variables]
   for nc, nc_temp, year, month in getInputs():
           with nc, nc_temp:
               nctime = nc.variables['time_counter']
               nctemp = nc_temp['votemper']
               first = ncvariables_out is None
               if first:
                  ncout.createDimension('time_counter', ntime)
                  ncout.createDimension('x', nc.dimensions['x'].size)
                  ncout.createDimension('y', nc.dimensions['y'].size)
                  copyVariable(ncout, nc.variables['nav_lat'])
                  copyVariable(ncout, nc.variables['nav_lon'])
                  nctemp_out = copyVariable(ncout, nctemp, dimensions=('time_counter', 'y', 'x'), copy_data = False)

                  nctime_out = copyVariable(ncout, nctime, copy_data=False)
                  kwargs = {}
                  if chunk:
                     kwargs['chunksizes'] = (ntime, 1, 1)
                  ncw_int_out = ncout.createVariable('bm_int', 'f', ('time_counter', 'y', 'x'), zlib=compress, contiguous=contiguous, **kwargs)
                  ncw2_int_out = ncout.createVariable('bm2_int', 'f', ('time_counter', 'y', 'x'), zlib=compress, contiguous=contiguous, **kwargs)
                  ncomega_out  = ncout.createVariable('omega', 'f', ('time_counter', 'y', 'x'), zlib=compress, contiguous=contiguous, **kwargs)
                  nctemp_bot_out = ncout.createVariable('votemper_bot', 'f', ('time_counter', 'y', 'x'), zlib=compress, contiguous=contiguous, **kwargs)
                  nctemp_out.coordinates = 'nav_lon nav_lat'
                  ncw_int_out.coordinates = 'nav_lon nav_lat'
                  ncw2_int_out.coordinates = 'nav_lon nav_lat'
                  nctemp_bot_out.coordinates = 'nav_lon nav_lat'
                  ncomega_out.coordinates = 'nav_lon nav_lat'
                  
                  ncvariables_out = []
               ncvariables_pel = []
               ncvariables_ben = []
               for variable in variables:
                  assert variable in nc.variables, 'Variable %s not found in %s. Available: %s' % (variable, path, ', '.join(nc.variables.keys()))
                  ncvar = nc.variables[variable]
                  if variable in pel_vars:
                      ncvariables_pel.append(ncvar)
                  elif variable in ben_vars:    
                      ncvariables_ben.append(ncvar)
                  else:
                      print ('{} is no defined as benthic or pelagic variable'.format (variable))
                      exit
                  if first:
                     ncvar_out = copyVariable(ncout, ncvar, dimensions=('time_counter', 'y', 'x'), copy_data = False)
                     ncvar_out.coordinates = 'nav_lon nav_lat'
                     ncvariables_out.append(ncvar_out)
               ncthickness = nc.variables['e3t']
               for day in range(nctime.size):
                  print('  %04i-%02i-%02i' % (year, month, day+1))
                  n=0
                  if iout >= istart:
                     #Pelagic data
                     alldata = [ncvariable[day, ...] for ncvariable in ncvariables_pel]
                     biomass = sum(alldata)
                 #    print (biomass.shape)
                     thickness = ncthickness[day, ...]
                     weights = biomass * thickness
                     weights_int = weights.sum(axis=0)
                     weights2_int = (biomass**2 * thickness).sum(axis=0)
                     ncw_int_out[iout, :, :] = weights_int
                     ncw2_int_out[iout, :, :] = weights2_int
                     weights /= weights_int[numpy.newaxis, :, :]
 #                    for ncvariable_out, data in zip(ncvariables_out, alldata):
                     for data in alldata:
                        ncvariables_out[n][iout, :, :] = (weights*data).sum(axis=0)
                        n+=1
                     nctemp_out[iout, :, :] = (weights*nctemp[day, ...]).sum(axis=0)
                     nctemp_bot_out[iout, :, :] = (nctemp[day,-1, ...]).sum(axis=0)
                     #benthic data
                     all_bendata = [ncvariable[day, ...] for ncvariable in ncvariables_ben]
                     ben_biomass = sum(all_bendata)
                     for data in all_bendata:
                        ncvariables_out[n][iout, :, :] = data
                        n+=1
                     omega = weights_int/(weights_int+ben_biomass)
                     ncomega_out  [iout, :, :] = omega      
                     nctime_out[iout] = nctime[day]
                  iout += 1
               if iout >= istart:
                  ncout.sync()

