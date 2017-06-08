import sys
import os

import numpy
from matplotlib import pyplot
from matplotlib.dates import datestr2num, date2num, num2date
import netCDF4

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../fabm/src/drivers/python'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../build/pyfabm/Release'))

import mizer

# ERSEM outputs to save (in tab-separated text output)
ersem_outputs = 'P1_c', 'P2_c', 'P3_c', 'P4_c', 'Z4_c', 'Z5_c', 'Z6_c', 'Y2_c', 'Y3_c', 'total_bioturbation_activity_at_bottom_calculator_result', 'total_bioirrigation_activity_at_bottom_calculator_result', 'Q17_c', 'Q17_n', 'Q17_p', 'sst', 'salt[-1]', 'mld_surf'

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
    'gamma': 0.5*64./50., #paper=64., we incorporate kappa=0.5, and convert from density to concentration for a column of 50 m!
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
    'F': 1.6
}

scenarios = 'reference', 'N3n=half', 'N3n=double', 'N1p=half', 'N1p=double'

for scenario_name in scenarios:
    print('Processing scenario %s...' % scenario_name)

    gotm_result = os.path.join(os.path.dirname(__file__), '../../gotm-fabm-ersem-configurations/L4/%s/L4_time_daily_mean_16.06.nc' % scenario_name)

    # Function for converting from Equivalent Spherical Diameter (micrometer) to wet mass in g
    def esd2mass(d): # d: equivalent spherical diameter in micrometer
        V = 4./3.*numpy.pi*(numpy.array(d)/2e6)**3  # V: volume in m3
        return V*1e6  # mass in g approximately equals volume in m3 multiplied by 1e6 (assumes density of 1000 kg/m3)

    # prey (currently from GOTM-ERSEM simulation)
    scale_factor = 0.01 # 10 g wet mass/g carbon * 0.001 g C/mg C
    prey = (
        mizer.Prey('diatoms', esd2mass((20,200)), mizer.datasources.TimeSeries(gotm_result, 'P1_c', weights='h', scale_factor=scale_factor, z='sum')),
        mizer.Prey('nanophy', esd2mass((2,20)), mizer.datasources.TimeSeries(gotm_result, 'P2_c', weights='h', scale_factor=scale_factor, z='sum')),
        mizer.Prey('picophy', esd2mass((.2,2)), mizer.datasources.TimeSeries(gotm_result, 'P3_c', weights='h', scale_factor=scale_factor, z='sum')),
        #mizer.Prey('microphy', esd2mass((20,200)), mizer.datasources.TimeSeries(gotm_result, 'P4_c', weights='h', scale_factor=scale_factor, z='sum')),
        mizer.Prey('microzoo', esd2mass((20,200)), mizer.datasources.TimeSeries(gotm_result, 'Z5_c', weights='h', scale_factor=scale_factor, z='sum')),
        mizer.Prey('nanoflag', esd2mass((2,20)), mizer.datasources.TimeSeries(gotm_result, 'Z6_c', weights='h', scale_factor=scale_factor, z='sum')),
        mizer.Prey('mesozoo', (1e-5,1e-3), mizer.datasources.TimeSeries(gotm_result, 'Z4_c', weights='h', scale_factor=scale_factor, z='sum')),
    )
    prey_collection = mizer.PreyCollection(*prey)
    prey_collection = mizer.GriddedPreyCollection(prey_collection)

    # environment
    temp = mizer.datasources.TimeSeries(gotm_result, 'temp', weights='h', z='mean')

    # create mizer model
    m = mizer.Mizer(prey=prey_collection, parameters=parameters, temperature=temp, recruitment_from_prey=True)

    # Time-integrate
    result = m.run(temp.times, spinup=50, verbose=True, save_spinup=False)

    #result.plot_spectrum()
    #result.plot_lfi_timeseries(500., 1.25)
    #result.plot_biomass_timeseries(0., 500.)
    #result.plot_timeseries('landings')
    #result.plot_annual_mean('landings', plot_change=True)
    #pyplot.show()

    biomass = result.get_biomass_timeseries()

    # Write tab-separated text file with ERSEM and mizer variables.
    # Depth-explicit ERSEM variables will be depth integrated.
    # The selection of ERSEM variables thta is written is taken from "ersem_outputs", defined near the top of this file.
    with open('%s.dat' % scenario_name, 'w') as f, netCDF4.Dataset(gotm_result) as nc:
        f.write('date')
        ersem_data = []
        h = nc.variables['h'][:, :, 0, 0]
        nctime = nc.variables['time']
        gotmtime = date2num(netCDF4.num2date(nctime[:], nctime.units))
        for name in ersem_outputs:
            prefix, postfix, location = '', '', None
            if name.endswith(']'):
                name, location = name.split('[')
                location = location[:-1]
            ncvar = nc.variables[name]
            dat = ncvar[..., 0, 0]
            units = ncvar.units
            if dat.ndim == 2:
                if location is None:
                    # depth integrate
                    dat = (dat*h).sum(axis=1)
                    prefix = 'depth-integrated '
                    units = '%s*m' % units
                else:
                    dat = dat[:, int(location)]
                    postfix = '@k=%i' % int(location)
            f.write('\t%s%s%s (%s)' % (prefix, ncvar.long_name, postfix, units))
            ersem_data.append(dat)
        f.write('\tfish biomass (g WM/m2)')
        f.write('\n')
        for time, bm in zip(result.t, biomass):
            f.write(num2date(time).strftime('%Y-%m-%d'))
            for dat in ersem_data:
                f.write('\t%s' % numpy.interp(time, gotmtime, dat))
            f.write('\t%s' % bm)
            f.write('\n')

