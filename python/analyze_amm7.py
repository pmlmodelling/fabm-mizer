import sys
import os
import glob
import datetime
import argparse

import numpy
from matplotlib import pyplot
from matplotlib.dates import datestr2num, date2num, num2date
import netCDF4

fabm_root = '../../fabm'
fabm_root = '../../fabm-git'
build_dir = '../../build/pyfabm'
build_dir = '../build'

sys.path.insert(0, os.path.realpath(os.path.join(os.path.dirname(__file__), fabm_root, 'src/drivers/python')))
sys.path.insert(0, os.path.realpath(os.path.join(os.path.dirname(__file__), build_dir, 'Release')))

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

def addVariable(nc, name, long_name, units, data=None, dimensions=None):
    if dimensions is None:
        dimensions = (time_name,)
    ncvar = nc.createVariable(name, float, dimensions, zlib=True, fill_value=-2e20)
    if data is not None:
        ncvar[:] = data
    ncvar.long_name = long_name
    ncvar.units = units
    return ncvar

def processLocation(args):
    assert len(args) == 4
    path, i, j, output_path = args
    print('Processing %s for i=%i, j=%i...' % (path, i, j))

    # prey (currently from GOTM-ERSEM simulation) - scale to g WM/m3
    scale_factor = 10*0.001 # 10 g wet mass/g carbon * 0.001 g C/mg C
    prey = []
    for name, ncname, size_range in preylist:
        timeseries = mizer.datasources.TimeSeries(path, ncname, scale_factor=scale_factor, time_name=time_name, x=i, y=j)
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
    result = m.run(times, spinup=spinup, verbose=True, save_spinup=False)

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
    return (times, biomass, landings, lfi80, lfi500, lfi10000)

    with netCDF4.Dataset(output_path, 'w') as ncout, netCDF4.Dataset(path) as nc:
        nctime_in = nc.variables[time_name]
        ncout.createDimension(time_name)
        nctime_out = ncout.createVariable(time_name, nctime_in.datatype, nctime_in.dimensions, zlib=True)
        nctime_out.units = nctime_in.units
        dates = [dt.replace(tzinfo=None) for dt in num2date(times)]
        nctime_out[...] = netCDF4.date2num(dates, nctime_out.units)
        addVariable(ncout, 'biomass', 'biomass', 'g WM/m3', biomass)
        addVariable(ncout, 'landings', 'landings', 'g WM', landings)
        addVariable(ncout, 'lfi80', 'fraction of fish > 80 g', '-', lfi80)
        addVariable(ncout, 'lfi500', 'fraction of fish > 500 g', '-', lfi500)
        addVariable(ncout, 'lfi10000', 'fraction of fish > 10000 g', '-', lfi10000)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('source_path')
    parser.add_argument('output_path')
    args = parser.parse_args()

    tasks = []
    if not os.path.isdir(args.output_path):
       os.mkdir(args.output_path)
    for path in glob.glob(args.source_path):
        with netCDF4.Dataset(path) as nc:
            mask = nc.variables['mask']
            for i in xrange(len(nc.dimensions['x'])):
                for j in xrange(len(nc.dimensions['y'])):
                    if mask[j, i] > 0:
                        tasks.append((args.source_path, i, j, os.path.join(args.output_path, 'i=%i,j=%i.nc' % (i,j))))

    #processLocation(tasks[0])
    #sys.exit(0)

    # Process all EEZs using all available cores
    # Kill child process after processing a single EEZ (maxtasksperchild=1) to prevent ever increasing memory consumption.
    import multiprocessing
    pool = multiprocessing.Pool(processes=None, maxtasksperchild=1)
    results = pool.map(processLocation, tasks)
    output_path = os.path.join(args.output_path, os.path.basename(args.source_path))
    with netCDF4.Dataset(path) as nc, netCDF4.Dataset(output_path, 'w') as ncout:
        nctime_in = nc.variables[time_name]
        ncout.createDimension(time_name)
        ncout.createDimension('x', len(nc.dimensions['x']))
        ncout.createDimension('y', len(nc.dimensions['y']))
        nctime_out = ncout.createVariable(time_name, nctime_in.datatype, nctime_in.dimensions, zlib=True)
        nctime_out.units = nctime_in.units
        times = results[0][0]
        dates = [dt.replace(tzinfo=None) for dt in num2date(times)]
        nctime_out[...] = netCDF4.date2num(dates, nctime_out.units)
        ncbiomass = addVariable(ncout, 'biomass', 'biomass', 'g WM/m3', dimensions=(time_name, 'y', 'x'))
        nclandings = addVariable(ncout, 'landings', 'landings', 'g WM', dimensions=(time_name, 'y', 'x'))
        nclfi80 = addVariable(ncout, 'lfi80', 'fraction of fish > 80 g', '-', dimensions=(time_name, 'y', 'x'))
        nclfi500 = addVariable(ncout, 'lfi500', 'fraction of fish > 500 g', '-', dimensions=(time_name, 'y', 'x'))
        nclfi10000 = addVariable(ncout, 'lfi10000', 'fraction of fish > 10000 g', '-', dimensions=(time_name, 'y', 'x'))
        for (source, i, j, output), (times, biomass, landings, lfi80, lfi500, lfi10000) in zip(tasks, results):
            ncbiomass[:, j, i] = biomass
            nclandings[:, j, i] = landings
            nclfi80[:, j, i] = lfi80
            nclfi500[:, j, i] = lfi500
            nclfi10000[:, j, i] = lfi10000

