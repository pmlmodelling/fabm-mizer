from __future__ import print_function

import sys
import os.path
import numpy
import scipy.integrate
from matplotlib import pyplot
from matplotlib import animation
import yaml
from matplotlib.dates import num2date

try:
    import pyfabm
except ImportError:
    print('Unable to load pyfabm. See https://github.com/fabm-model/code/wiki/python.')
    sys.exit(1)

import datasources

class Prey(object):
    def __init__(self, name, mass, value):
        self.name = name
        self.mass = mass
        if not isinstance(value, datasources.ValueProvider):
            value = datasources.Constant(value)
        self.value_provider = value

class Mizer(object):
    def __init__(self, path=os.path.join(os.path.dirname(__file__), 'mizer.yaml'), prey=(), parameters={}, temperature=None):
        with open(path, 'rU') as f:
            configuration = yaml.load(f)

        mizer_params = configuration['parameters']
        mizer_params.update(parameters)
        mizer_coupling = {'waste': 'zero_hz'}
        mizer_yaml = {'model': 'mizer/size_structured_population', 'parameters': mizer_params, 'coupling': mizer_coupling}
        fabm_yaml = {'instances': {'fish': mizer_yaml}}

        self.prey_items = list(prey)
        for name, data in configuration.get('prey', {}):
            self.prey_items.append(Prey(name, data['wet_mass'], datasources.Constant(data['constant_value'])))

        iprey = 1
        for prey_item in self.prey_items:
            fabm_yaml['instances'][prey_item.name] = {'model': 'mizer/prey', 'parameters': {'w': prey_item.mass}}
            mizer_coupling['Nw_prey%i' % iprey] = '%s/Nw' % prey_item.name
            iprey += 1

        mizer_params['nprey'] = len(self.prey_items)
        with open('fabm.yaml', 'w') as f:
            yaml.dump(fabm_yaml, f, default_flow_style=False)

        self.fabm_model = pyfabm.Model('fabm.yaml')

        self.prey_indices = []
        for prey_item in self.prey_items:
            for i, variable in enumerate(self.fabm_model.state_variables):
                if variable.path == '%s/Nw' % prey_item.name:
                    break
            else:
                assert False, 'Prey %s/Nw not found in model.' % prey_item.name
            self.prey_indices.append(i)

        self.bin_indices = []
        self.bin_masses = []
        while 1:
            for i, variable in enumerate(self.fabm_model.state_variables):
                if variable.path == 'fish/Nw%i' % (len(self.bin_indices)+1,):
                    self.bin_masses.append(variable.getRealProperty('particle_mass'))
                    break
            else:
                break
            self.bin_indices.append(i)
        self.bin_masses = numpy.array(self.bin_masses)
        log10masses = numpy.log10(self.bin_masses)
        log10binwidth = log10masses[1] - log10masses[0]
        self.bin_widths = 10.**(log10masses+log10binwidth/2) - 10.**(log10masses-log10binwidth/2)

        # Used to convert between depth-integrated fluxes and depth-explicit fluxes (not used if prey is prescribed)
        self.fabm_model.findDependency('bottom_depth').value = 1.

        if temperature is not None and not isinstance(temperature, datasources.ValueProvider):
            temperature = datasources.Constant(temperature)
        self.temperature_provider = temperature
        if self.temperature_provider is None and 'temperature' in configuration:
            self.temperature_provider = datasources.Constant(configuration['temperature']['constant_value'])
        if self.temperature_provider is not None:
            self.temperature = self.fabm_model.findDependency('temperature')
            self.temperature.value = self.temperature_provider.mean()

        # Verify the model is ready to be used
        assert self.fabm_model.checkReady(), 'One or more model dependencies have not been fulfilled.'

        for parameter in self.fabm_model.parameters:
            if parameter.path.startswith('fish/'):
                print('%s: %s %s' % (parameter.long_name, parameter.value, parameter.units))

    def run(self, t, verbose=False, spinup=0):
        def dy(y, current_time):
            self.fabm_model.state[:] = y
            if not in_spinup:
                for i, prey_item in zip(self.prey_indices, self.prey_items):
                    self.fabm_model.state[i] = prey_item.value_provider.get(current_time)
                if self.temperature is not None:
                    self.temperature.value = self.temperature_provider.get(current_time)
            return self.fabm_model.getRates()

        # Time-integrate (note: FABM's internal time unit is seconds!)
        if spinup > 0:
            in_spinup = True
            t_spinup = numpy.arange(t[0]-365.23*spinup, t[0], 1.)
            for i, prey_item in zip(self.prey_indices, self.prey_items):
                self.fabm_model.state[i] = prey_item.value_provider.mean()
            if self.temperature is not None:
                self.temperature.value = self.temperature_provider.mean()
            if verbose:
                print('Spinning up from %s to %s d' % (num2date(t_spinup[0]), num2date(t_spinup[-1])))
            y = scipy.integrate.odeint(dy, self.fabm_model.state, t_spinup*86400)

        if verbose:
            print('Time integrating from %s to %s d' % (num2date(t[0]), num2date(t[-1])))
        in_spinup = False
        y = scipy.integrate.odeint(dy, self.fabm_model.state, t*86400)
        return MizerResult(self, t, y)

class MizerResult(object):
    def __init__(self, model, t, y):
        self.model = model
        self.t = t
        self.y = y
        self.spectrum = y[:, self.model.bin_indices]
        self.biomass_density = self.spectrum/self.model.bin_widths
        self.abundance_density = self.spectrum/self.model.bin_widths/self.model.bin_masses

    def plot_spectrum(self, itime=-1, fig=None, normalization=0, global_range=False):
        if fig is None:
            fig = pyplot.figure()
        ax = fig.gca()
        style = '.' if normalization == 0 else '-'
        prey_masses, prey_values, lines = None, None, []
        if normalization == 0:
            values = self.spectrum
            ax.set_ylabel('wet mass (g)')
            prey_masses = numpy.array([prey_item.mass for prey_item in self.model.prey_items])
            prey_values = numpy.array([prey_item.value_provider.get(self.t[itime]) for prey_item in self.model.prey_items])
        elif normalization == 1:
            values = self.biomass_density
            ax.set_ylabel('wet mass density (g/g)')
        elif normalization == 2:
            values = self.abundance_density
            ax.set_ylabel('abundance density (#/g)')
        if global_range:
            minval, maxval = values.min(), values.max()
            if prey_masses is not None:
                minval, maxval = min(minval, prey_values.min()), max(maxval, prey_values.max())
            ax.set_ylim(minval/10, maxval*10)
        if prey_masses is not None:
            line, = ax.loglog(prey_masses, prey_values, '.')
            lines.append(line)
        line, = ax.loglog(self.model.bin_masses, values[itime, :], style)
        lines.append(line)
        ax.grid(True)
        ax.set_xlabel('wet mass (g)')
        return lines

    def plot_biomass_timeseries(self, fig=None):
        if fig is None:
            fig = pyplot.figure()
        ax = fig.gca()
        line, = ax.plot_date(self.t, self.spectrum[:, :].sum(axis=1), '-')
        ax.set_xlabel('time (d)')
        ax.set_ylabel('biomass (%s)' % self.model.fabm_model.state_variables[self.model.bin_indices[0]].units)
        ax.grid(True)
        return line,

    def plot_timeseries(self, name, fig=None):
        for i, variable in enumerate(self.model.fabm_model.state_variables):
            if variable.path == 'fish/%s' % name:
                break
        else:
            raise Exception('Variable "%s" not found in mizer output' % name)
        if fig is None:
            fig = pyplot.figure()
        ax = fig.gca()
        line, = ax.plot_date(self.t, self.y[:, i], '-')
        ax.set_xlabel('time (d)')
        ax.set_ylabel('%s (%s)' % (self.model.fabm_model.state_variables[i].long_name, self.model.fabm_model.state_variables[i].units))
        ax.grid(True)
        return line,

    def animate_spectrum(self, dir='.', normalization=0):
        fig = pyplot.figure()
        lines = self.plot_spectrum(0, fig=fig, normalization=normalization, global_range=True)
        prey_values = None
        if normalization == 0:
            values = self.spectrum
            prey_values = numpy.empty((len(self.t), len(self.model.prey_items)))
            for iprey, prey_item in enumerate(self.model.prey_items):
                prey_values[:, iprey] = prey_item.value_provider.get(self.t)
        elif normalization == 1:
            values = self.biomass_density
        elif normalization == 2:
            values = self.abundance_density
        def new_frame(itime):
            if prey_values is not None:
                lines[0].set_ydata(prey_values[itime, :])
            lines[-1].set_ydata(values[itime, :])
            return lines
        return animation.FuncAnimation(fig, new_frame, frames=self.spectrum.shape[0], interval=1000./30)

if __name__ == '__main__':
    # Time-integrate over 200 days (note: FABM's internal time unit is seconds!)
    m = Mizer()
    t = numpy.arange(0, 365*100., 1.)
    result = m.run(t, verbose=True)

    # Plot results
    result.animate_spectrum()
