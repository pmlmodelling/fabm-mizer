from __future__ import print_function

import sys
import os.path
import numpy
import scipy.integrate
from matplotlib import pyplot
import yaml

try:
    import pyfabm
except ImportError:
    print('Unable to load pyfabm. See https://github.com/fabm-model/code/wiki/python.')
    sys.exit(1)

with open('mizer.yaml', 'rU') as f:
    configuration = yaml.load(f)

mizer_params = configuration['parameters']
prey_data = configuration.get('prey', {})
mizer_params['nprey'] = len(prey_data)
mizer_coupling = {'waste': 'zero_hz'}
mizer_yaml = {'model': 'mizer/size_structured_population', 'parameters': mizer_params, 'coupling': mizer_coupling}
fabm_yaml = {'instances': {'fish': mizer_yaml}}
iprey = 1
for name, data in prey_data.items():
    fabm_yaml['instances'][name] = {'model': 'mizer/prey', 'parameters': {'w': data['wet_mass']}}
    mizer_coupling['Nw_prey%i' % iprey] = '%s/Nw' % name
    iprey += 1
with open('fabm.yaml', 'w') as f:
    yaml.dump(fabm_yaml, f, default_flow_style=False)

fabm_model = pyfabm.Model('fabm.yaml')

prey_indices = []
prey_values = []
for name, data in prey_data.items():
    for iprey, prey in enumerate(fabm_model.state_variables):
        if prey.path == '%s/Nw' % name:
            break
    else:
        assert False, 'Prey %s/Nw not found in model.' % name
    prey_indices.append(iprey)
    prey_values.append(data['constant_value'])

# Used to convert between depth-integrated fluxes and depth-explicit fluxes (not used if rpey is prescribed)
fabm_model.findDependency('bottom_depth').value = 1.
if 'temperature' in configuration:
    fabm_model.findDependency('temperature').value = configuration['temperature']['constant_value']

# Verify the model is ready to be used
assert fabm_model.checkReady(), 'One or more model dependencies have not been fulfilled.'

class Mizer(object):
    def __init__(self):
        pass

    def run(self, t, verbose=False):
        def dy(y, t0):
            fabm_model.state[:] = y
            for i, value in zip(prey_indices, prey_values):
                fabm_model.state[i] = value
            return fabm_model.getRates()

        # Time-integrate (note: FABM's internal time unit is seconds!)
        if verbose:
            print('Time integrating from %s to %s d' % (t[0], t[-1]))
        y = scipy.integrate.odeint(dy, fabm_model.state, t*86400)
        return MizerResult(self, t, y)

class MizerResult(object):
    def __init__(self, model, t, y):
        self.model = model
        self.t = t
        self.y = y

        n = 0
        istart = fabm_model.findStateVariable('fish/Nw1').index
        while 1:
            try:
                fabm_model.findStateVariable('fish/Nw%i' % (n+1))
            except KeyError:
                break
            n += 1

        self.spectrum = y[:, istart:istart+n]

    def plot_spectrum(self, itime=-1, fig=None):
        if fig is None:
            fig = pyplot.figure()
        ax = fig.gca()
        line, = ax.semilogy(self.spectrum[itime, :])
        ax.grid(True)
        return line

    def animate_spectrum(self, dir='.', prefix='still'):
        fig = pyplot.figure()
        line = self.plot_spectrum(0, fig=fig)
        fig.gca().set_ylim(self.spectrum.min()/10, self.spectrum.max()*10)
        print('Saving animation to %s/%s?????.png...' % (dir, prefix))
        for itime in range(self.spectrum.shape[0]):
            line.set_ydata(self.spectrum[itime, :])
            path = os.path.join(dir, '%s%05i.png' % (prefix, itime+1))
            print('  saving %s...' % path)
            fig.savefig(path)

if __name__ == '__main__':
    # Time-integrate over 200 days (note: FABM's internal time unit is seconds!)
    m = Mizer()
    t = numpy.arange(0, 365*100., 1.)
    result = m.run(t, verbose=True)

    # Plot results
    result.animate_spectrum()
