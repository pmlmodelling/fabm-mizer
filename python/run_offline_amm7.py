from __future__ import print_function

import sys
import os
import glob
import datetime
import argparse
import re
import shutil
import gc

gc.set_debug(gc.DEBUG_UNCOLLECTABLE)

import yaml

import numpy
from matplotlib import pyplot
from matplotlib.dates import datestr2num, date2num, num2date
import netCDF4

fabm_root = '../../fabm'
fabm_root = '../../fabm-git'
build_dir = '../../build/pyfabm'
build_dir = '../build'

sys.path.insert(0, os.path.realpath(os.path.join(os.path.dirname(__file__), fabm_root, 'src/drivers/python')))
sys.path.insert(0, os.path.realpath(os.path.join(os.path.dirname(__file__), build_dir, 'Release')))  # Visual Studio/Windows
sys.path.insert(0, os.path.realpath(os.path.join(os.path.dirname(__file__), build_dir)))             # Linux

start_time = datetime.datetime(2009, 1, 1)
stop_time = datetime.datetime(2012, 1, 1)

start_time = None
stop_time = None

import mizer

# Function for converting from Equivalent Spherical Diameter (micrometer) to wet mass in g
def esd2mass(*d): # d: equivalent spherical diameter in micrometer
    V = 4./3.*numpy.pi*(numpy.array(d)/2e6)**3  # V: volume in m3
    return V*1e6  # mass in g approximately equals volume in m3 multiplied by 1e6 (assumes density of 1000 kg/m3)

additional_outputs = []

preylist = []
preylist.append(('diatoms', 'P1_c', esd2mass(20., 200.)))
preylist.append(('nanophytoplankton', 'P2_c',   esd2mass(2., 20.)))
preylist.append(('picophytoplankton', 'P3_c',   esd2mass(.2, 2.)))
preylist.append(('microphytoplankton', 'P4_c', esd2mass(20., 200.)))
preylist.append(('mesozooplankton', 'Z4_c', (4.188e-6, 1e-3)))
preylist.append(('microzooplankton', 'Z5_c', esd2mass(20., 200.)))
preylist.append(('heterotrophic nanoflagellates', 'Z6_c', esd2mass(2., 20.)))
temp_name = 'votemper'
time_name = 'time_counter'

# mizer parameters
parameters = {
    'w_min': 1e-3,
    'w_inf': 1e6,
    'nclass': 100,
    'T_dependence': 1,
    'T_ref': 13.,
    'E_a': 0.63,
    'beta': 100,
    'sigma': float(numpy.log(10.)),   # paper has log10 units, we use ln units
    'gamma': 156, # clearance in m3/yr for single individual of mass 1 g. Blanchard et al 2009: 640 m3/yr; Blanchard et al 2012: 64 ?UNITS? [times kappa=0.5 for time spent in pelagic]; Faking giants paper gives 10^14.076 * W^0.926 * exp(-Ea/(kT) L d-1, which is 428 L d-1 = 156 m3 yr-1
    'q': 0.82,
    'alpha': 0.2,
    'z0_type': 1,
    'z0pre': 0.1,
    'z0exp': -0.25,
    'w_s': 1000.,
    'z_s': 0.3,
    'ks': 0.,
    'SRR': 0,
    'recruitment': 0.,
    'h': 1e9,
    'fishing_type': 1,
    'w_minF': 1.25, # Blanchard et al 2012
    'F': 0.4
}

def addVariable(nc, name, long_name, units, data=None, dimensions=None, zlib=False, contiguous=True, dtype='f4'):
    if dimensions is None:
        dimensions = (time_name,)
    chunksizes = None
    if not contiguous:
        chunksizes = []
        for dim in dimensions:
           chunksizes.append(1 if dim in ('x', 'y') else len(nc.dimensions[dim]))
    ncvar = nc.createVariable(name, dtype, dimensions, zlib=zlib, fill_value=-2e20, contiguous=contiguous, chunksizes=chunksizes)
    ncvar.set_var_chunk_cache(0, 0, 0)
    if data is not None:
        ncvar[...] = data
    ncvar.long_name = long_name
    ncvar.units = units
    if 'x' in dimensions and 'y' in dimensions and 'nav_lon' in nc.variables and 'nav_lat' in nc.variables:
       ncvar.coordinates = 'nav_lon nav_lat'
    return ncvar

def copyVariable(nc, ncvar, **kwargs):
   ncvar_out = nc.createVariable(ncvar.name, ncvar.dtype, ncvar.dimensions, fill_value=getattr(ncvar, '_FillValue', None), **kwargs)
   for key in ncvar.ncattrs():
      if key != '_FillValue':
         setattr(ncvar_out, key, getattr(ncvar, key))
   if 'x' in ncvar.dimensions and 'y' in ncvar.dimensions and 'nav_lon' in nc.variables and 'nav_lat' in nc.variables:
      ncvar_out.coordinates = 'nav_lon nav_lat'
   ncvar_out[...] = ncvar[...]
   return ncvar_out

def processLocation(args):
    path, i, j = args
    print('Processing %s for i=%i, j=%i...' % (path, i, j))

    # prey (currently from GOTM-ERSEM simulation) - scale to g WM/m3
    scale_factor = 10*0.001 # 10 g wet mass/g carbon * 0.001 g C/mg C
    prey = []
    for name, ncname, size_range in preylist:
        timeseries = mizer.datasources.TimeSeries(path, ncname, scale_factor=scale_factor, time_name=time_name, x=i, y=j, stop=stop_time)
        times = timeseries.times
        prey.append(mizer.Prey(name, size_range, timeseries))
    prey_collection = mizer.PreyCollection(*prey)
    prey_collection = mizer.GriddedPreyCollection(prey_collection)

    # environment
    temp = mizer.datasources.TimeSeries(path, temp_name, time_name=time_name, x=i, y=j, stop=stop_time)
    depth = mizer.datasources.TimeSeries(path, 'bm_int**2/bm2_int', time_name=time_name, x=i, y=j, stop=stop_time)
    #temp = 12.

    # create mizer model
    m = mizer.Mizer(prey=prey_collection, parameters=parameters, temperature=temp, recruitment_from_prey=1, depth=depth)

    # Time-integrate
    spinup = 50
    istart, istop = 0, times.size
    if start_time is not None:
        istart = times.searchsorted(date2num(start_time))
    if stop_time is not None:
        istop = times.searchsorted(date2num(stop_time))
    times = times[istart:istop]

    result = m.run(times, spinup=spinup, verbose=True, save_spinup=False,dt=1/(24*4))

    if result is None:
        return
    #result.plot_spectrum()
    #result.plot_lfi_timeseries(500., 1.25)
    #result.plot_biomass_timeseries(0., 500.)
    #result.plot_timeseries('landings')
    #result.plot_annual_mean('landings', plot_change=True)
    #pyplot.show()

    biomass = result.get_biomass_timeseries()
    landings_var, landings = result.get_timeseries('landings')
    lfi10 = result.get_lfi_timeseries(10.)
    lfi80 = result.get_lfi_timeseries(80.)
    lfi500 = result.get_lfi_timeseries(500.)
    lfi10000 = result.get_lfi_timeseries(10000.)
    landings[1:] = landings[1:] - landings[:-1]
    landings[0] = 0
    return path, i, j, times, biomass, landings, lfi10, lfi80, lfi500, lfi10000, m.bin_masses, result.spectrum

def ppProcessLocation(args, p):
    import run_offline_amm7
    run_offline_amm7.parameters = p
    return run_offline_amm7.processLocation(args)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('source_path')
    parser.add_argument('output_path')
    parser.add_argument('--method', choices=('serial', 'multiprocessing', 'pp'), default='pp')
    parser.add_argument('--ncpus', type=int, default=None)
    parser.add_argument('--ppservers', default=None)
    parser.add_argument('--secret', default=None)
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--parameters', default=None)
    parser.add_argument('--ifirst', type=int, default=None)
    parser.add_argument('--shm', action='store_true')
    parser.add_argument('--profile', action='store_true')
    args = parser.parse_args()

    if args.parameters is not None:
        with open(args.parameters, 'rU') as f:
            args.parameters = yaml.load(f)
        parameters = args.parameters

    if isinstance(args.ppservers, (str, u''.__class__)):
        match = re.match(r'(.*)\[(.*)\](.*)', args.ppservers)
        if match is not None:
            # Hostnames in PBS/SLURM notation, e.g., node[01-06]
            ppservers = []
            left, middle, right = match.groups()
            for item in middle.split(','):
                if '-' in item:
                    start, stop = item.split('-')
                    for i in range(int(start), int(stop)+1):
                        ppservers.append('%s%s%s' % (left, str(i).zfill(len(start)), right))
                else:
                    ppservers.append('%s%s%s' % (left, item, right))
        else:
            # Comma-separated hostnames
            ppservers = args.ppservers.split(',')
        ppservers = tuple(ppservers)
    else:
        assert args.ppservers is None
        ppservers = ()
    if args.ncpus is None:
        args.ncpus = 'autodetect'

    tasks = []
    if not os.path.isdir(args.output_path):
       os.mkdir(args.output_path)
    paths = glob.glob(args.source_path)
    assert len(paths) > 0, 'no files found at %s' % args.source_path
    for path in paths:
        print('Opening %s...' % path)
        with netCDF4.Dataset(path) as nc:
            if 'mask' in nc.variables:
                mask = nc.variables['mask'][...] > 0
            else:
                mask = (nc.variables['bm_int'][...] > 0).any(axis=0)
            for i in range(len(nc.dimensions['x'])):
                for j in range(len(nc.dimensions['y'])):
                    if mask[j, i]:
                        tasks.append((path, i, j))
    if args.ifirst is not None:
        tasks = tasks[args.ifirst:]

    source2output = {}
    source2vars = {}
    def getOutput(source, times, w, compress=False, add_biomass_per_bin=False, contiguous=False):
        if source not in source2output:
           output_path = os.path.join(args.output_path, os.path.basename(source))
           if args.ifirst is not None:
              assert os.path.isfile(output_path)
              ncout = netCDF4.Dataset(output_path, 'r+')
           else:
              with netCDF4.Dataset(path) as nc:
                print('Creating output file %s...' % output_path, end='')
                ncout = netCDF4.Dataset(output_path, 'w') #, persist=True, diskless=True)
                nctime_in = nc.variables[time_name]
                ncout.createDimension(time_name, len(times))
                ncout.createDimension('x', len(nc.dimensions['x']))
                ncout.createDimension('y', len(nc.dimensions['y']))
                nctime_out = ncout.createVariable(time_name, nctime_in.datatype, nctime_in.dimensions, zlib=compress, contiguous=contiguous)
                nctime_out.units = nctime_in.units
                dates = [dt.replace(tzinfo=None) for dt in num2date(times)]
                nctime_out[...] = netCDF4.date2num(dates, nctime_out.units)
                if 'nav_lon' in nc.variables:
                    copyVariable(ncout, nc.variables['nav_lon'], zlib=compress)
                    copyVariable(ncout, nc.variables['nav_lat'], zlib=compress)
                vardict = {}
                vardict['mask'] = ncout.createVariable('mask', 'i1', ('y', 'x'), zlib=compress, contiguous=contiguous)
                vardict['mask'][...] = 0
                vardict['biomass'] = addVariable(ncout, 'biomass', 'biomass', 'g WM/m2', dimensions=(time_name, 'y', 'x'), zlib=compress, contiguous=contiguous)
                vardict['landings'] = addVariable(ncout, 'landings', 'landings', 'g WM', dimensions=(time_name, 'y', 'x'), zlib=compress, contiguous=contiguous)
                vardict['lfi10'] = addVariable(ncout, 'lfi10', 'fraction of fish > 10 g', '-', dimensions=(time_name, 'y', 'x'), zlib=compress, contiguous=contiguous)
                vardict['lfi80'] = addVariable(ncout, 'lfi80', 'fraction of fish > 80 g', '-', dimensions=(time_name, 'y', 'x'), zlib=compress, contiguous=contiguous)
                vardict['lfi500'] = addVariable(ncout, 'lfi500', 'fraction of fish > 500 g', '-', dimensions=(time_name, 'y', 'x'), zlib=compress, contiguous=contiguous)
                vardict['lfi10000'] = addVariable(ncout, 'lfi10000', 'fraction of fish > 10000 g', '-', dimensions=(time_name, 'y', 'x'), zlib=compress, contiguous=contiguous)
                if add_biomass_per_bin:
                    ncout.createDimension('bin', w.size)
                    addVariable(ncout, 'w', 'individual mass', 'g WM', dimensions=('bin',), zlib=compress)[:] = w
                    vardict['Nw'] = addVariable(ncout, 'Nw', 'biomass per bin', 'g WM/m2', dimensions=(time_name, 'y', 'x', 'bin'), zlib=compress, contiguous=contiguous)
                    vardict['Nw_final'] = addVariable(ncout, 'Nw_final', 'final biomass per bin', 'g WM/m2', dimensions=('y', 'x', 'bin'), zlib=compress, contiguous=contiguous)
                print('done')
           source2output[source] = ncout
           source2vars[source] = vardict
        return source2output[source], source2vars[source]

    def saveResult(result, sync=True, add_biomass_per_bin=False):
        source, i, j, times, biomass, landings, lfi10,lfi80, lfi500, lfi10000, w, spectrum = result
        ncout, vardict = getOutput(source, times, w, add_biomass_per_bin=add_biomass_per_bin)
        print('saving results from %s, i=%i, j=%i (mean biomass = %.3g)' % (source, i, j, biomass.mean()))
        vardict['biomass'][:, j, i] = biomass
        vardict['landings'][:, j, i] = landings
        vardict['lfi10'][:,j,i]= lfi10
        vardict['lfi80'][:, j, i] = lfi80
        vardict['lfi500'][:, j, i] = lfi500
        vardict['lfi10000'][:, j, i] = lfi10000
        vardict['mask'][j, i] = 1
        if add_biomass_per_bin:
           vardict['Nw'][:, j, i, :] = spectrum
           vardict['Nw_final'][j, i, :] = spectrum[-1, :]
        if sync:
           print('Synchronizing NetCDF output to disk...', end='')
           ncout.sync()
           print('done')

    job_server = None
    final_output_path = None
    if args.method == 'serial':
        def runSerial(n):
            for i in range(n):
                saveResult(processLocation(tasks[i]))
        if args.profile:
            import cProfile
            import pstats
            cProfile.run('runSerial(%s)' % min(len(tasks), 3), 'mizerprof')
            p = pstats.Stats('mizerprof')
            p.strip_dirs().sort_stats('cumulative').print_stats()
        else:
            runSerial(min(len(tasks), 3))
    elif args.method == 'multiprocessing':
        # Process all EEZs using all available cores
        # Kill child process after processing a single EEZ (maxtasksperchild=1) to prevent ever increasing memory consumption.
        import multiprocessing
        pool = multiprocessing.Pool(processes=None, maxtasksperchild=1)

        #results = pool.map(processLocation, tasks)
        #for result in  results:
        #    saveResult(result)

        #result = pool.map_async(processLocation, tasks, callback=saveResult)
        #result.wait()

        for result in pool.imap(processLocation, tasks):
            saveResult(result, sync=False, add_biomass_per_bin=True)
    else:
        if args.debug:
            import logging
            logging.basicConfig( level=logging.DEBUG)
        import pp
        if args.shm:
            final_output_path = args.output_path
            args.output_path = '/dev/shm'
        job_server = pp.Server(ncpus=args.ncpus, ppservers=ppservers, restart=True, secret=args.secret)
        jobs = []
        for task in tasks:
            jobs.append(job_server.submit(ppProcessLocation, (task, args.parameters)))
        ijob = 0
        nfailed = 0
        while jobs:
            job = jobs.pop(0)
            result = job()
            sync = ijob % 1000 == 0
            if result is not None:
               print('job %i: saving result...' % ijob)
               saveResult(result, sync=sync, add_biomass_per_bin=True)
               if sync:
                  gc.collect()
                  print(gc.garbage)
            else:
               print('job %i: FAILED!' % ijob)
               nfailed += 1
            ijob += 1
        print('%i tasks out of %i FAILED.' % (nfailed, len(tasks)))
        job_server.print_stats()
 
    for source, nc in source2output.items():
        name = os.path.basename(source)
        print('Closing %s...' % os.path.join(args.output_path, name))
        nc.close()
        if final_output_path is not None:
           target = os.path.join(final_output_path, name)
           if os.path.isfile(target):
              os.remove(target)
           shutil.move(os.path.join(args.output_path, name), target)

    if job_server is not None:
       job_server.destroy()

