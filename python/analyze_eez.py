import sys
import os
import glob
import datetime

import numpy
from matplotlib import pyplot
from matplotlib.dates import datestr2num, date2num, num2date
import netCDF4

import slurm

platform = 'ceto'
source = 'IPSL-CM5A-LR'

if platform == 'ceto':
    root = '/work/jbr/FAO/EEZ-data'
    build_dir = '../build'
    fabm_root = '../../fabm-git'
else:
    root = 'C:/Users/jbr/OneDrive/PML/FAO/EEZ-data'
    root = 'C:/Users/Jorn/OneDrive/PML/FAO/EEZ-data'
    root = 'D:/temp/EEZ-data'
    fabm_root = '../../fabm'
    build_dir = '../../build/pyfabm'
    build_dir = '../build-python'

sys.path.insert(0, os.path.join(os.path.dirname(__file__), fabm_root, 'src/drivers/python'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), build_dir, 'Release')) # windows-only
sys.path.insert(0, os.path.join(os.path.dirname(__file__), build_dir))

import mizer

additional_outputs = []

preylist = []

if source.startswith('IPSL-'):
    # PISCES size classes, Olivier Maury pers comm 20/22 Jun 2017
    preylist.append(('diatoms', 'phydiat', (10., 100.)))
    #preylist.append(('miscellaneous phytoplankton', 'phymisc', (2., 20.)))
    preylist.append(('microzooplankton', 'zmicro', (20., 200.)))
    preylist.append(('mesozooplankton', 'zmeso', (200., 2000.)))
    temp_name = 'tos-273.15' # note IPSL models express SST in Kelvin!
else:
    assert False, 'don\'t know prey and temperature variable names for %s' % source

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

def addVariable(nc, name, long_name, units, data):
    ncvar = nc.createVariable(name, float, ('time',), zlib=True)
    ncvar[:] = data
    ncvar.long_name = long_name
    ncvar.units = units

# Function for converting from Equivalent Spherical Diameter (micrometer) to wet mass in g
def esd2mass(d): # d: equivalent spherical diameter in micrometer
    V = 4./3.*numpy.pi*(numpy.array(d)/2e6)**3  # V: volume in m3
    return V*1e6  # mass in g approximately equals volume in m3 multiplied by 1e6 (assumes density of 1000 kg/m3)

def processEEZ(eez_name):
    prefix = 'EEZ%s_Omon_%s_' % (eez_name, source)
    postfix = '_r1i1p1.nc'
    files = map(os.path.abspath, glob.glob(os.path.join(root, source, 'merged/%s*%s' % (prefix, postfix))))
    historical = map(os.path.abspath, glob.glob(os.path.join(root, source, 'merged/%shistorical%s' % (prefix, postfix))))[0]
    files.remove(historical)
    files.insert(0, historical)

    initial_state = None
    for forcing_file in files:
        print('Processing %s...' % forcing_file)

        # prey (currently from GOTM-ERSEM simulation)
        scale_factor = 10*12.011 # 10 g wet mass/g carbon * 12.011 g C/mol C
        prey = []
        for name, ncname, size_range in preylist:
            prey.append(mizer.Prey(name, esd2mass(size_range), mizer.datasources.TimeSeries(forcing_file, ncname, scale_factor=scale_factor)))
        prey_collection = mizer.PreyCollection(*prey)
        prey_collection = mizer.GriddedPreyCollection(prey_collection)

        # environment
        temp = mizer.datasources.TimeSeries(forcing_file, temp_name)

        # create mizer model
        m = mizer.Mizer(prey=prey_collection, parameters=parameters, temperature=temp, recruitment_from_prey=True)

        times = temp.times
        if initial_state is not None:
            ilast = times.searchsorted(date2num(datetime.datetime(2100, 1, 1)))
            times = times[:ilast]

        # Time-integrate
        spinup = 50 if initial_state is None else 0
        result = m.run(times, spinup=spinup, verbose=True, save_spinup=False, initial_state=initial_state)

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

        if initial_state is None:
            initial_state = result.y[-1, :]

        output_dir = os.path.join(root, source, 'fish')
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        with netCDF4.Dataset(os.path.join(output_dir, os.path.basename(forcing_file)), 'w') as ncout, netCDF4.Dataset(forcing_file) as nc:
            nctime_in = nc.variables['time']
            ncout.createDimension('time')
            nctime_out = ncout.createVariable('time', nctime_in.datatype, nctime_in.dimensions, zlib=True)
            nctime_out.units = nctime_in.units
            dates = [dt.replace(tzinfo=None) for dt in num2date(times)]
            nctime_out[...] = netCDF4.date2num(dates, nctime_out.units)
            addVariable(ncout, 'biomass', 'biomass', 'g WM/m3', biomass)
            addVariable(ncout, 'landings', 'landings', 'g WM', landings)
            addVariable(ncout, 'lfi80', 'fraction of fish > 80 g', '-', lfi80)
            addVariable(ncout, 'lfi500', 'fraction of fish > 500 g', '-', lfi500)
            addVariable(ncout, 'lfi10000', 'fraction of fish > 10000 g', '-', lfi10000)

def ppProcessEEZ(*args):
    import analyze_eez
    return analyze_eez.processEEZ(*args)

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--method', choices=('serial', 'multiprocessing', 'pp'), default='serial')
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--ncpus', type=int, default=None)
    parser.add_argument('--ppservers', default=None)
    parser.add_argument('--secret', default=None)
    args = parser.parse_args()

    # Build list of EEZ names
    eez_names = []
    for path in glob.glob(os.path.join(root, source, 'merged/EEZ*_Omon_%s_historical_r1i1p1.nc' % source)):
        eez_names.append(os.path.basename(path).split('_', 1)[0][3:])

    #processEEZ('231')
    #processEEZ('228')

    # Process all EEZs using all available cores
    # Kill child process after processing a single EEZ (maxtasksperchild=1) to prevent ever increasing memory consumption.
    if args.method == 'serial':
        for eez_name in eez_names:
            processEEZ(eez_name)
    elif args.method == 'multiprocessing':
        import multiprocessing
        pool = multiprocessing.Pool(processes=None, maxtasksperchild=1)
        pool.map(processEEZ, eez_names)
    else:
        if args.ncpus is None:
            args.ncpus = 'autodetect'
        if args.debug:
            import logging
            logging.basicConfig( level=logging.DEBUG)
        import pp
        job_server = pp.Server(ncpus=args.ncpus, ppservers=slurm.getNodes(args.ppservers), restart=True, secret=args.secret)
        jobs = []
        for eez_name in eez_names:
            jobs.append(job_server.submit(ppProcessEEZ, (eez_name,)))
        for ijob, job in enumerate(jobs):
            result = job()
        job_server.print_stats()
        job_server.destroy()
