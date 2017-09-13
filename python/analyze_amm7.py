import sys
import os
import glob
import datetime
import argparse
import re

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
temp_name = 'tos-273.15'
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

def addVariable(nc, name, long_name, units, data=None, dimensions=None, zlib=False):
    if dimensions is None:
        dimensions = (time_name,)
    ncvar = nc.createVariable(name, float, dimensions, zlib=zlib, fill_value=-2e20)
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
    #temp = mizer.datasources.TimeSeries(forcing_file, temp_name)
    temp = 12.

    # create mizer model
    m = mizer.Mizer(prey=prey_collection, parameters=parameters, temperature=temp, recruitment_from_prey=True)

    # Time-integrate
    spinup = 50
    istart = times.searchsorted(date2num(start_time))
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
    return path, i, j, times, biomass, landings, lfi80, lfi500, lfi10000

def ppProcessLocation(args):
    import analyze_amm7
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
    args = parser.parse_args()

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
            mask = nc.variables['mask']
            for i in xrange(len(nc.dimensions['x'])):
                for j in xrange(len(nc.dimensions['y'])):
                    if mask[j, i] > 0:
                        tasks.append((path, i, j))

    source2output = {}
    def getOutput(source, times, compress=False):
        if source not in source2output:
            with netCDF4.Dataset(path) as nc:
                output_path = os.path.join(args.output_path, os.path.basename(source))
                ncout = netCDF4.Dataset(output_path, 'w')
                nctime_in = nc.variables[time_name]
                ncout.createDimension(time_name)
                ncout.createDimension('x', len(nc.dimensions['x']))
                ncout.createDimension('y', len(nc.dimensions['y']))
                nctime_out = ncout.createVariable(time_name, nctime_in.datatype, nctime_in.dimensions, zlib=compress)
                nctime_out.units = nctime_in.units
                dates = [dt.replace(tzinfo=None) for dt in num2date(times)]
                nctime_out[...] = netCDF4.date2num(dates, nctime_out.units)
                ncbiomass = addVariable(ncout, 'biomass', 'biomass', 'g WM/m3', dimensions=(time_name, 'y', 'x'), zlib=compress)
                nclandings = addVariable(ncout, 'landings', 'landings', 'g WM', dimensions=(time_name, 'y', 'x'), zlib=compress)
                nclfi80 = addVariable(ncout, 'lfi80', 'fraction of fish > 80 g', '-', dimensions=(time_name, 'y', 'x'), zlib=compress)
                nclfi500 = addVariable(ncout, 'lfi500', 'fraction of fish > 500 g', '-', dimensions=(time_name, 'y', 'x'), zlib=compress)
                nclfi10000 = addVariable(ncout, 'lfi10000', 'fraction of fish > 10000 g', '-', dimensions=(time_name, 'y', 'x'), zlib=compress)
            source2output[source] = ncout
        return source2output[source]

    def saveResult(result, sync=True):
        source, i, j, times, biomass, landings, lfi80, lfi500, lfi10000 = result
        print('saving results from %s, i=%i, j=%i' % (source, i, j))
        ncout = getOutput(source, times)
        ncout.variables['biomass'][:, j, i] = biomass
        ncout.variables['landings'][:, j, i] = landings
        ncout.variables['lfi80'][:, j, i] = lfi80
        ncout.variables['lfi500'][:, j, i] = lfi500
        ncout.variables['lfi10000'][:, j, i] = lfi10000
        if sync:
           ncout.sync()

    if args.method == 'serial':
        import cProfile
        import pstats
        cProfile.run('saveResult(processLocation(tasks[0]))', 'mizerprof')
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
            saveResult(resulti, sync=False)
    else:
        if args.debug:
            import logging
            logging.basicConfig( level=logging.DEBUG)
        import pp
        job_server = pp.Server(ncpus=args.ncpus, ppservers=ppservers, restart=True, secret=args.secret)
        jobs = []
        for task in tasks:
            jobs.append(job_server.submit(ppProcessLocation, (task,)))
        for ijob, job in enumerate(jobs):
            result = job()
            if result is not None:
               print('job %i: saving result...' % ijob)
               saveResult(result, sync=False)
            else:
               print('job %i: FAILED!' % ijob)
        job_server.print_stats()
        job_server.destroy()
 
    for nc in source2output.values():
        nc.close()
