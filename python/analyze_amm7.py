import sys
import os
import glob
import datetime
import argparse
import re
import shutil

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
preylist.append(('mesozooplankton', 'Z4_c', (1e-5, 1e-3)))
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
    'gamma': 365., # QuestFish paper=64. times kappa=0.5; Faking giants paper gives approx 1 m3/d at 13 degrees, i.e., 365 m3/yr for a fish of 1 g. For 1 fish of 1 g per m3 the specific clearance rate is therefore 365 yr-1
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
    'w_minF': 1.25, # Blanchard et al 2012
    'F': 0.8  # note: need to put double the intended value due to fisheries equation!
}

def addVariable(nc, name, long_name, units, data=None, dimensions=None, zlib=False, contiguous=True):
    if dimensions is None:
        dimensions = (time_name,)
    chunksizes = [1]*len(dimensions)
    if time_name in dimensions:
        chunksizes[dimensions.index(time_name)] = len(nc.dimensions[time_name])
    ncvar = nc.createVariable(name, float, dimensions, zlib=zlib, fill_value=-2e20, contiguous=contiguous, chunksizes=chunksizes)
    if data is not None:
        ncvar[:] = data
    ncvar.long_name = long_name
    ncvar.units = units
    return ncvar

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
    m = mizer.Mizer(prey=prey_collection, parameters=parameters, temperature=temp, recruitment_from_prey=True, depth=depth)

    # Time-integrate
    spinup = 50
    #istart = times.searchsorted(date2num(start_time))
    istart = 0
    istop = times.searchsorted(date2num(stop_time))
    times = times[istart:istop]

    result = m.run(times, spinup=spinup, verbose=True, save_spinup=False)

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
    lfi80 = result.get_lfi_timeseries(80.)
    lfi500 = result.get_lfi_timeseries(500.)
    lfi10000 = result.get_lfi_timeseries(10000.)
    landings[1:] = landings[1:] - landings[:-1]
    landings[0] = 0
    return path, i, j, times, biomass, landings, lfi80, lfi500, lfi10000, result.spectrum

def ppProcessLocation(args, p):
    import analyze_amm7
    analyze_amm7.parameters = p
    return analyze_amm7.processLocation(args)

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
    args = parser.parse_args()

    if args.parameters is not None:
        with open(args.parameters, 'rU') as f:
            args.parameters = yaml.load(f)

    if isinstance(args.ppservers, basestring):
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
    for path in glob.glob(args.source_path):
        with netCDF4.Dataset(path) as nc:
            if 'mask' in nc.variables:
                mask = nc.variables['mask'][...] > 0
            else:
                mask = (nc.variables['bm_int'][...] > 0).any(axis=0)
            for i in xrange(len(nc.dimensions['x'])):
                for j in xrange(len(nc.dimensions['y'])):
                    if mask[j, i]:
                        tasks.append((path, i, j))

    source2output = {}
    def getOutput(source, times, nbins, compress=False, add_biomass_per_bin=False, contiguous=False):
        if source not in source2output:
            with netCDF4.Dataset(path) as nc:
                output_path = os.path.join(args.output_path, os.path.basename(source))
                ncout = netCDF4.Dataset(output_path, 'w')
                nctime_in = nc.variables[time_name]
                ncout.createDimension(time_name, len(times))
                ncout.createDimension('x', len(nc.dimensions['x']))
                ncout.createDimension('y', len(nc.dimensions['y']))
                nctime_out = ncout.createVariable(time_name, nctime_in.datatype, nctime_in.dimensions, zlib=compress, contiguous=contiguous)
                nctime_out.units = nctime_in.units
                dates = [dt.replace(tzinfo=None) for dt in num2date(times)]
                nctime_out[...] = netCDF4.date2num(dates, nctime_out.units)
                ncbiomass = addVariable(ncout, 'biomass', 'biomass', 'g WM/m2', dimensions=(time_name, 'y', 'x'), zlib=compress, contiguous=contiguous)
                nclandings = addVariable(ncout, 'landings', 'landings', 'g WM', dimensions=(time_name, 'y', 'x'), zlib=compress, contiguous=contiguous)
                nclfi80 = addVariable(ncout, 'lfi80', 'fraction of fish > 80 g', '-', dimensions=(time_name, 'y', 'x'), zlib=compress, contiguous=contiguous)
                nclfi500 = addVariable(ncout, 'lfi500', 'fraction of fish > 500 g', '-', dimensions=(time_name, 'y', 'x'), zlib=compress, contiguous=contiguous)
                nclfi10000 = addVariable(ncout, 'lfi10000', 'fraction of fish > 10000 g', '-', dimensions=(time_name, 'y', 'x'), zlib=compress, contiguous=contiguous)
                if add_biomass_per_bin:
                    for i in range(nbins):
                        ncbm = addVariable(ncout, 'Nw%i' % (i+1), 'biomass in bin %i' % (i + 1), 'g WM/m2', dimensions=(time_name, 'y', 'x'), zlib=compress, contiguous=contiguous)
            source2output[source] = ncout
        return source2output[source]

    def saveResult(result, sync=True, add_biomass_per_bin=False):
        source, i, j, times, biomass, landings, lfi80, lfi500, lfi10000, spectrum = result
        print('saving results from %s, i=%i, j=%i' % (source, i, j))
        ncout = getOutput(source, times, spectrum.shape[1], add_biomass_per_bin=add_biomass_per_bin)
        ncout.variables['biomass'][:, j, i] = biomass
        ncout.variables['landings'][:, j, i] = landings
        ncout.variables['lfi80'][:, j, i] = lfi80
        ncout.variables['lfi500'][:, j, i] = lfi500
        ncout.variables['lfi10000'][:, j, i] = lfi10000
        if add_biomass_per_bin:
            for i in range(spectrum.shape[1]):
                ncout.variables['Nw%i' % (i+1)][:, j, i] = spectrum[:, i]
        if sync:
           ncout.sync()

    job_server = None
    final_output_path = None
    if args.method == 'serial':
        import cProfile
        import pstats
        def runSerial(n):
            for i in range(n):
                saveResult(processLocation(tasks[i]))
        cProfile.run('runSerial(3)', 'mizerprof')
        p = pstats.Stats('mizerprof')
        p.strip_dirs().sort_stats('cumulative').print_stats()
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
        final_output_path = args.output_path
        args.output_path = '/dev/shm'
        job_server = pp.Server(ncpus=args.ncpus, ppservers=ppservers, restart=True, secret=args.secret)
        jobs = []
        for task in tasks:
            jobs.append(job_server.submit(ppProcessLocation, (task, args.parameters)))
        for ijob, job in enumerate(jobs):
            result = job()
            if result is not None:
               print('job %i: saving result...' % ijob)
               saveResult(result, sync=False, add_biomass_per_bin=True)
            else:
               print('job %i: FAILED!' % ijob)
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

