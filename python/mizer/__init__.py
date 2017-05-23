from __future__ import print_function

import sys
import os.path
import yaml
import numpy
import scipy.integrate
import datetime
from matplotlib import pyplot
from matplotlib import animation
from matplotlib.dates import num2date, date2num

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
        self.value_provider = datasources.asValueProvider(value)

class BasePreyCollection(object):
    def __init__(self):
        self.names = None
        self.masses = None

class PreyCollection(BasePreyCollection):
    def __init__(self, *items):
        BasePreyCollection.__init__(self)
        self.items = items
        self.names = [item.name for item in self.items]
        self.masses = numpy.array([10.**numpy.mean(numpy.log10(item.mass)) for item in self.items])

    def getValues(self, time):
        result = numpy.empty(numpy.shape(time)+self.masses.shape)
        for i, item in enumerate(self.items):
            result[..., i] = item.value_provider.get(time)
        return result

    def getMean(self):
        return numpy.array([item.value_provider.mean() for item in self.items])

class GriddedPreyCollection(BasePreyCollection):
    def __init__(self, source):
        BasePreyCollection.__init__(self)
        self.source = source

        min_prey_mass, max_prey_mass = None, None
        for prey_item in self.source.items:
            current_min = prey_item.mass if isinstance(prey_item.mass, float) else prey_item.mass[0]
            min_prey_mass = current_min if min_prey_mass is None else min(min_prey_mass, current_min)
            current_max = prey_item.mass if isinstance(prey_item.mass, float) else prey_item.mass[1]
            max_prey_mass = current_max if max_prey_mass is None else max(max_prey_mass, current_max)
        self.delta_log10mass = 0.1
        prey_mass_bounds = numpy.arange(numpy.log10(min_prey_mass)-self.delta_log10mass/2, numpy.log10(max_prey_mass)+self.delta_log10mass, self.delta_log10mass)
        prey_mass_centers = (prey_mass_bounds[1:]+prey_mass_bounds[:-1])/2
        self.prey_bin_weights = []
        for prey_item in self.source.items:
            log10mass = numpy.log10(prey_item.mass)
            weights = numpy.zeros_like(prey_mass_centers)
            if isinstance(prey_item.mass, float):
                i = numpy.argmin(numpy.abs(prey_mass_centers-log10mass))
                weights[i] = 1.
            else:
                for ibin in range(len(prey_mass_centers)):
                    left  = max(log10mass[0], prey_mass_bounds[ibin])
                    right = min(log10mass[1], prey_mass_bounds[ibin+1])
                    weights[ibin] = max(0., right-left)/(log10mass[1]-log10mass[0])
            assert abs(weights.sum()-1.) < 1e-12, '%s: weights should add up to 1, but currently add up to %s' % (prey_item.name, weights.sum())
            self.prey_bin_weights.append(weights)

        self.names = ['prey%i' % i for i in range(len(prey_mass_centers))]
        self.masses = 10.**prey_mass_centers

    def getValues(self, time):
        result = numpy.zeros(numpy.shape(time) + self.masses.shape)
        for item, weights in zip(self.source.items, self.prey_bin_weights):
            result += weights*numpy.expand_dims(item.value_provider.get(time), -1)
        return result

    def getMean(self):
        result = numpy.zeros_like(self.masses)
        for item, weights in zip(self.source.items, self.prey_bin_weights):
            result += weights*item.value_provider.mean()
        return result

class Mizer(object):
    def __init__(self, parameters={}, prey=(), temperature=None, recruitment_from_prey=False):
        fabm_yaml_path = 'fabm.yaml'
        mizer_params = dict(parameters)
        mizer_coupling = {'waste': 'zero_hz'}
        mizer_yaml = {'model': 'mizer/size_structured_population', 'parameters': mizer_params, 'coupling': mizer_coupling}
        fabm_yaml = {'instances': {'fish': mizer_yaml}}

        if not isinstance(prey, BasePreyCollection):
            prey = PreyCollection(*prey)
        self.prey = prey

        iprey = 0
        for name, mass in zip(self.prey.names, self.prey.masses):
            iprey += 1
            fabm_yaml['instances'][name] = {'model': 'mizer/prey', 'parameters': {'w': float(mass)}}
            mizer_coupling['Nw_prey%i' % iprey] = '%s/Nw' % name
        mizer_params['nprey'] = iprey
        with open(fabm_yaml_path, 'w') as f:
            yaml.dump(fabm_yaml, f, default_flow_style=False)

        self.fabm_model = pyfabm.Model('fabm.yaml')

        self.prey_indices = numpy.empty((iprey,), dtype=int)
        for iprey, name in enumerate(self.prey.names):
            for self.prey_indices[iprey], variable in enumerate(self.fabm_model.state_variables):
                if variable.path == '%s/Nw' % name:
                    break
            else:
                assert False, 'Prey %s/Nw not found in model.' % name

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
        self.log10bin_width = log10masses[1] - log10masses[0]   # assume equal log10 spacing between mizer size classes
        self.bin_widths = 10.**(log10masses+self.log10bin_width/2) - 10.**(log10masses-self.log10bin_width/2)

        # Used to convert between depth-integrated fluxes and depth-explicit fluxes (not used if prey is prescribed)
        self.fabm_model.findDependency('bottom_depth').value = 1.

        self.temperature_provider = None
        if temperature is not None:
            self.temperature_provider = datasources.asValueProvider(temperature)
            self.temperature = self.fabm_model.findDependency('temperature')
            self.temperature.value = self.temperature_provider.mean()

        # Verify the model is ready to be used
        assert self.fabm_model.checkReady(), 'One or more model dependencies have not been fulfilled.'

        for parameter in self.fabm_model.parameters:
            if parameter.path.startswith('fish/'):
                print('%s: %s %s' % (parameter.long_name, parameter.value, parameter.units))

        self.initial_state = numpy.copy(self.fabm_model.state)
        self.recruitment_from_prey = recruitment_from_prey

    def run(self, t, verbose=False, spinup=0, save_spinup=False):
        # Shortcuts to objects used during time integration
        state = self.fabm_model.state
        getRates = self.fabm_model.getRates

        # Build list of indices of state variables that must be kept constant.
        iconstant_states = set(self.prey_indices)
        if self.recruitment_from_prey:
            iconstant_states.add(self.bin_indices[0])
        iconstant_states = list(sorted(iconstant_states))

        def dy(y, current_time):
            state[:] = y
            if not in_spinup:
                state[self.prey_indices] = self.prey.getValues(current_time)
                if self.temperature is not None:
                    self.temperature.value = self.temperature_provider.get(current_time)
                if self.recruitment_from_prey:
                    state[self.bin_indices[0]] = state[self.prey_indices].mean()*self.log10bin_width/self.prey.delta_log10mass
            rates = getRates()*86400
            rates[iconstant_states] = 0.
            return rates

        state[:] = self.initial_state

        # Spin up with time-averaged prey abundances
        t_spinup, y_spinup = None, None
        if spinup > 0:
            in_spinup = True
            t_spinup = t[0] - numpy.arange(0., 365.23*spinup, 1.)[::-1]
            state[self.prey_indices] = self.prey.getMean()
            if self.recruitment_from_prey:
                state[self.bin_indices[0]] = state[self.prey_indices].mean()*self.log10bin_width/self.prey.delta_log10mass
            if self.temperature is not None:
                self.temperature.value = self.temperature_provider.mean()
            if verbose:
                print('Spinning up from %s to %s' % (num2date(t_spinup[0]), num2date(t_spinup[-1])))
            y_spinup = scipy.integrate.odeint(dy, state, t_spinup)
            state[:] = y_spinup[-1, :]

        if verbose:
            print('Time integrating from %s to %s' % (num2date(t[0]), num2date(t[-1])))
        in_spinup = False
        y = scipy.integrate.odeint(dy, state, t)

        # Overwrite prey masses with imposed values.
        y[:, self.prey_indices] = self.prey.getValues(t)
        if self.recruitment_from_prey:
            y[:, self.bin_indices[0]] = y[:, self.prey_indices].mean(axis=1)*self.log10bin_width/self.prey.delta_log10mass

        if spinup > 0 and save_spinup:
            # Prepend spinup to results
            t = numpy.hstack((t_spinup[:-1], t))
            y = numpy.vstack((y_spinup[:-1, :], y))

        if verbose:
            print('Done.')
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
            prey_masses = self.model.prey.masses
            prey_values = self.y[itime, self.model.prey_indices]
        elif normalization == 1:
            values = self.biomass_density
            ax.set_ylabel('wet mass density (g/g)')
        elif normalization == 2:
            values = self.abundance_density
            ax.set_ylabel('abundance density (#/g)')
        if global_range:
            minval, maxval = values.min(), values.max()
            if prey_masses is not None:
                prey_min = prey_values.min()
                if prey_min > 0:
                    minval = min(minval, prey_min)
                maxval = max(maxval, prey_values.max())
            ax.set_ylim(minval/10, maxval*10)
        if prey_masses is not None:
            line, = ax.loglog(prey_masses, prey_values, '.')
            lines.append(line)
        line, = ax.loglog(self.model.bin_masses, values[itime, :], style)
        lines.append(line)
        ax.grid(True)
        ax.set_xlabel('wet mass (g)')
        title = ax.set_title(num2date(self.t[itime]).strftime('%Y-%m-%d'))
        return tuple(lines) + (title,)

    def get_biomass_timeseries(self, min_weight=None, max_weight=None):
        istart = 0
        if min_weight is not None:
            istart = self.model.bin_masses.searchsorted(min_weight)
        istop = self.spectrum.shape[1]
        if max_weight is not None:
            istop = self.model.bin_masses.searchsorted(max_weight)
        return self.spectrum[:, istart:istop].sum(axis=1)

    def get_lfi_timeseries(self, lfi_weight, min_weight=None):
        assert min_weight is None or min_weight < lfi_weight
        return self.get_biomass_timeseries(lfi_weight)/self.get_biomass_timeseries(min_weight)

    def plot_biomass_timeseries(self, min_weight=None, max_weight=None, fig=None):
        if fig is None:
            fig = pyplot.figure()
        ax = fig.gca()
        line, = ax.plot_date(self.t, self.get_biomass_timeseries(min_weight, max_weight), '-')
        ax.set_xlabel('time (d)')
        ax.set_ylabel('biomass (%s)' % self.model.fabm_model.state_variables[self.model.bin_indices[0]].units)
        ax.grid(True)
        return line,

    def plot_lfi_timeseries(self, lfi_weight, min_weight=None, fig=None):
        if fig is None:
            fig = pyplot.figure()
        ax = fig.gca()
        line, = ax.plot_date(self.t, self.get_lfi_timeseries(lfi_weight, min_weight), '-')
        ax.set_xlabel('time (d)')
        ax.set_ylabel('fraction of popluation > %s g WM' % lfi_weight)
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

    def plot_annual_mean(self, name, fig=None, plot_change=False):
        for i, variable in enumerate(self.model.fabm_model.state_variables):
            if variable.path == 'fish/%s' % name:
                break
        else:
            raise Exception('Variable "%s" not found in mizer output' % name)
        all_data = self.y[:, i]

        start_year, stop_year = num2date(self.t[0]).year, num2date(self.t[-1]).year
        year_bounds = [date2num(datetime.datetime(year, 1, 1, 0, 0, 0)) for year in range(start_year, stop_year+1)]
        if year_bounds[0] < self.t[0]:
            year_bounds.pop(0)
            start_year += 1
        #print('Years completely spanned by simulation: %i - %i' % (start_year, stop_year-1))
        year_data = numpy.empty((stop_year - start_year))
        istart = 0
        for iyear in range(year_data.size):
            istop = self.t.searchsorted(year_bounds[iyear+1])
            if plot_change:
                year_data[iyear] = all_data[min(istop, self.t.size-1)] - all_data[istart]
            else:
                year_data[iyear] = all_data[istart:istop].mean()
            istart = istop

        if fig is None:
            fig = pyplot.figure()
        ax = fig.gca()
        bars = ax.bar(numpy.arange(start_year, stop_year), year_data, align='center')
        ax.set_xlabel('year')
        ax.get_xaxis().get_major_formatter().set_useOffset(False)
        if plot_change:
            ax.set_ylabel('change in %s (%s)' % (self.model.fabm_model.state_variables[i].long_name, self.model.fabm_model.state_variables[i].units))
        else:
            ax.set_ylabel('annual mean %s (%s)' % (self.model.fabm_model.state_variables[i].long_name, self.model.fabm_model.state_variables[i].units))
        ax.grid(True, axis='y')
        return bars

    def animate_spectrum(self, dir='.', normalization=0):
        fig = pyplot.figure()
        objects = self.plot_spectrum(0, fig=fig, normalization=normalization, global_range=True)
        lines, title = objects[:-1], objects[-1]
        prey_values = None
        if normalization == 0:
            values = self.spectrum
            prey_values = self.y[:, self.model.prey_indices]
        elif normalization == 1:
            values = self.biomass_density
        elif normalization == 2:
            values = self.abundance_density
        dates = num2date(self.t)
        def new_frame(itime):
            if prey_values is not None:
                lines[0].set_ydata(prey_values[itime, :])
            lines[-1].set_ydata(values[itime, :])
            title.set_text(dates[itime].strftime('%Y-%m-%d'))
            return objects
        return animation.FuncAnimation(fig, new_frame, frames=self.spectrum.shape[0], interval=1000./30, blit=True)

if __name__ == '__main__':
    # Time-integrate over 200 days (note: FABM's internal time unit is seconds!)
    m = Mizer()
    t = numpy.arange(0, 365*100., 1.)
    result = m.run(t, verbose=True)

    # Plot results
    result.animate_spectrum()
